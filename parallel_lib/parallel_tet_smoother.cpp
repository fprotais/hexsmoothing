#include "parallel_tet_smoother.h"    
#include <numeric>

#include "ultimaille/syntactic-sugar/assert.h"
#include "lbfgs_wrapper.h"

#include <utils/logTime.h>

using namespace Eigen;

inline double chi(double eps, double det) {
	const double eps2 = eps * eps;
    return det > 0 ? 
            (det + std::sqrt(eps2 + det * det)) * .5 :
            .5 * eps2 / (std::sqrt(eps2 + det * det) - det);
}

inline double chi_deriv(double eps, double det) {
    return .5 + det / (2. * std::sqrt(eps * eps + det * det));
}

Parallel_tet_smoother::Parallel_tet_smoother(utilities::TetrahedralMesh &mesh) 
: _mesh(mesh)
, _refs(mesh._tets.size(), equilateralTetRef)
, _vertDimLock(3*mesh._pts.size(), false)
{}

void Parallel_tet_smoother::setRefs(std::vector<std::array<utilities::vec3, 4>> const &refs) {
    um_assert(refs.size() == _refs.size());
    _refs.assign(refs.begin(), refs.end());
}

void Parallel_tet_smoother::setRef(unsigned tet, std::array<utilities::vec3, 4> const &ref) {
    um_assert(tet < _refs.size());
    _refs[tet] = ref;
}



void Parallel_tet_smoother::setVertsDimLocks(std::vector<bool> const &vertDimLock) {
    um_assert(vertDimLock.size() == _vertDimLock.size());
    _vertDimLock.assign(vertDimLock.begin(), vertDimLock.end());
}

void Parallel_tet_smoother::setVertDimLocks(unsigned vert, unsigned dim, bool lock) {
    um_assert(vert < _mesh._pts.size());
    um_assert(dim < 3);
    _vertDimLock[3*vert+dim] = lock;
}

void Parallel_tet_smoother::setVertsLocks(std::vector<bool> const &vertLock) {
    um_assert(vertLock.size() == _mesh._pts.size());
    for (unsigned v = 0; v < _mesh._pts.size(); ++v) {
        for (unsigned d = 0; d < 3; ++d) {
            _vertDimLock[3*v+d] = vertLock[v];
        }
    }
}

void Parallel_tet_smoother::setVertLock(unsigned vert, bool lock) {
    um_assert(vert < _mesh._pts.size());
    _vertDimLock[3*vert+0] = lock;
    _vertDimLock[3*vert+1] = lock;
    _vertDimLock[3*vert+2] = lock;
}

bool Parallel_tet_smoother::vertIsLocked(unsigned v) const {
    return _vertDimLock[3*v + 0] && _vertDimLock[3*v + 1] && _vertDimLock[3*v + 2];
}



void Parallel_tet_smoother::checkRefs() {
    for (auto ref : _refs) {
        um_assert(cross(ref[1]-ref[0], ref[2]-ref[0]) * (ref[3]-ref[0]) > 0);
    }
    if (verbose) std::cout << "-All refs are valid with positive det.-" << std::endl;
}

void Parallel_tet_smoother::createCompressSortedData() {
    if (verbose) std::cout << "    Pre-processing data from tetrahedra... " << std::endl;
    unsigned nbPoints = _mesh._pts.size();
    std::vector<unsigned> nbTetsOnVert(nbPoints);
    _tetData.reserve(_mesh._tets.size());
    for (unsigned i = 0; i < _mesh._tets.size(); ++i) {
        bool hasUnlocked = false;
        for (int v : _mesh._tets[i]) {
            if (!vertIsLocked(v)) {
                hasUnlocked = true;
                break;
            }
        }
        if (!hasUnlocked) continue;
        TetStorage storage;
        for (unsigned j = 0; j < 4; ++j) {
            storage.verts[j] = _mesh._tets[i][j];
            ++nbTetsOnVert[storage.verts[j]];
        }

        utilities::vec3 e0 = _refs[i][1] - _refs[i][0];
        utilities::vec3 e1 = _refs[i][2] - _refs[i][0];
        utilities::vec3 e2 = _refs[i][3] - _refs[i][0];

        Matrix3d M {{e0.x, e0.y, e0.z}, {e1.x, e1.y, e1.z}, {e2.x, e2.y, e2.z}};
        Matrix3d invM = M.inverse();
        storage.ig = { -invM.col(0) - invM.col(1) - invM.col(2), invM.col(0), invM.col(1), invM.col(2) };

        _tetData.push_back(storage);
    }
    _tetData.shrink_to_fit();
    _tetGrad.resize(_tetData.size());

    if (verbose) std::cout << "    Compressing vertices... " << std::endl;
    _originalToCompressed.resize(nbPoints);
    unsigned nbCompressed = 0;
    std::vector<unsigned> compressedToOriginal;
    compressedToOriginal.reserve(nbPoints);
    for (unsigned v = 0; v < nbPoints; ++v) {
        if (nbTetsOnVert[v] == 0) {
            _originalToCompressed[v] = -1;
            continue;
        }
        compressedToOriginal.push_back(v);
        _originalToCompressed[v] = nbCompressed++;
    }
    compressedToOriginal.shrink_to_fit();
    _compressedCoord.resize(3*nbCompressed);
    _compressedLocks.resize(3*nbCompressed);
    for (unsigned v = 0; v < nbPoints; ++v) {
        if (_originalToCompressed[v] == -1) continue;
        for (unsigned d = 0; d < 3; ++d) {
            _compressedCoord[3*_originalToCompressed[v]+d] = _mesh._pts[v][d];
            _compressedLocks[3*_originalToCompressed[v]+d] = _vertDimLock[3*v+d];
        }
    }

    for (TetStorage &tet : _tetData) {
        for (unsigned i = 0; i < 4; ++i) {
            tet.verts[i] = _originalToCompressed[tet.verts[i]];
            um_assert(tet.verts[i] != -1);
        }
    }

    if (verbose) std::cout << "    Storing vert to tet connectivity for parallelism..." << std::endl;
    std::vector<unsigned> currId(nbCompressed,0);
    _vert2tetCorner.reserve(nbCompressed);
    for (unsigned v = 0; v < nbCompressed; ++v) {
        _vert2tetCorner.push_back(std::vector<unsigned>(nbTetsOnVert[compressedToOriginal[v]]));
    }
    for (unsigned t = 0; t < _tetData.size(); ++t) {
        for (unsigned i = 0; i < 4; ++i) {
            unsigned v = _tetData[t].verts[i];
            _vert2tetCorner[v][currId[v]++] = 4*t+i;
        }
    }

}

bool Parallel_tet_smoother::go() {
    if (verbose) std::cout << "==== Running Parallel Elliptic smoother 3d. ====" << "\n";
    TimeLog logging("Parallel_tet_smoother");
    checkRefs();
    createCompressSortedData();
    if (verbose) logging.logSubStep("Initialization");
    if (_compressedCoord.size() == 0) {
        std::cout << "No variables to optimize." << std::endl;
        return true;
    }
    if (verbose) std::cout << "Nb of optimization variables (vertices x3): " << _compressedCoord.size()  << std::endl;
    if (verbose) std::cout << "Nb of energy terms (tetrahedra): " << _tetData.size() << std::endl;

    _eps0 = 1e-3;
    _eps = _eps0;

    bool res = _detMin > 0;
    if (verbose) std::cout << "Initial energy: " << elliptic_energy(_compressedCoord) << " | eps: " << _eps << " detmin: " << _detMin << " ninv: " << _nbInverted << std::endl;

    updateEps(0.5);
    _eps0 = _eps;

    for (unsigned iter = 0; iter < max_untangling_iter; ++iter) {
        if (verbose) std::cout << "Smoothing iteration #" << iter << "\n";
        if (verbose) std::cout << "    curr eps: " << _eps << std::endl;;
        
        double e_prev = elliptic_energy(_compressedCoord); // both with same eps
        if (!fineTimeLogging && verbose) std::cout << "    optimizing"<<std::flush;
        TimeLog fineLogging("untangling iter " + std::to_string(iter));
        LBFGS_wrapper opt ([&](Eigen::VectorXd const &x, Eigen::VectorXd &g){
            return elliptic_energy(x, &g);
        });
        opt._callBack = [&](VectorXd const &x, VectorXd const &g,double f,unsigned iternb,unsigned nbEval) {
            if (fineTimeLogging)  fineLogging.logSubStep("LBFGS iter " + std::to_string(iternb), "Energy: " + std::to_string(f));
            else if (verbose) std::cout << "." << std::flush;
            return false;
        };
        opt.setMaxIter(max_lbfgs_iter);
        bool optRes = opt.optimize(_compressedCoord);
        if (!fineTimeLogging && verbose) std::cout << " done." << std::endl;
        if (fineTimeLogging) fineLogging.logTotalTime();
        double e = elliptic_energy(_compressedCoord);
        if (verbose) std::cout << "    E: " << e_prev << " -> " << e  << " | " << " detmin: " << _detMin << " ninv: " << _nbInverted << "\n";
        if (verbose) std::cout << "    Status: " << opt.getMessage() << std::endl;

        updateEps(e/e_prev);
        res = _detMin > 0; 
        if (verbose) logging.logSubStep("untangling iter " + std::to_string(iter+1) );

        if (res && optRes) {
            if (verbose) std::cout << "Optimization converged and no tangled elements." << std::endl;
            break;
        };
    }
    if (verbose) logging.logTotalTime();
    
    updateVertCoord();
    return res;
}

void Parallel_tet_smoother::updateEps(double rate) {
    // mix of foldover and 2002 epsilons.
    if (_detMin > 0) {
        _eps = 1e-16;
        return;
    }
    double sigma = std::max(1. - rate, 0.1);
    double mu = (1 - sigma) * chi(_eps, _detMin);
    _eps = _detMin < mu ? 2 * std::sqrt(mu * (mu - _detMin)) : _noUntanglingEps;
    _eps = std::min(_eps, std::sqrt(_eps0 * _eps0 + 0.04 * _detMin * _detMin));
    _eps0 *= 0.85;
}


void Parallel_tet_smoother::updateVertCoord() {
    for (unsigned v = 0; v < _mesh._pts.size(); ++v) {
        if (_originalToCompressed[v] == -1) continue;
        unsigned comp = _originalToCompressed[v];
        for (unsigned d = 0; d < 3; ++d) {
            if (_vertDimLock[3*v+d]) continue;
            _mesh._pts[v][d] = _compressedCoord(3*comp + d);
        }
    }
}

inline Matrix<double,1,3> subVector(Eigen::VectorXd const &x, unsigned i) {
    i *= 3;
    return {x(i), x(i+1), x(i+2)};
}

double Parallel_tet_smoother::elliptic_energy(Eigen::VectorXd const &x, Eigen::VectorXd *G) {
    _detMin = std::numeric_limits<double>::max();
    _nbInverted = 0;
    double F = 0;
#pragma omp parallel for reduction(+:F,_nbInverted) reduction(min: _detMin)
    for (unsigned t = 0; t < _tetData.size(); ++t) {
        TetStorage const &tet = _tetData[t];

        Matrix3d J = tet.ig[0] * subVector(x,tet.verts[0])
                   + tet.ig[1] * subVector(x,tet.verts[1])
                   + tet.ig[2] * subVector(x,tet.verts[2])
                   + tet.ig[3] * subVector(x,tet.verts[3]);
                   

        double det = J.determinant();
        _detMin = std::min(det, _detMin);
        _nbInverted += det <= 0;

        double c1 = chi(_eps, det);
        double c2 = std::pow(c1, 2. / 3.);

        double f = J.squaredNorm() / c2;

        F +=  f;

        if (G == nullptr) continue;

        double c3 = chi_deriv(_eps, det);
        Matrix3d K; // dual basis
        K.col(0) = J.col(1).cross(J.col(2)); 
        K.col(1) = J.col(2).cross(J.col(0)); 
        K.col(2) = J.col(0).cross(J.col(1)); 
        Matrix3d dfdJ = J * (2. / c2) - K * ((2. * f * c3) / (3. * c1));
        dfdJ.transposeInPlace();
        for (unsigned v = 0; v < 4; ++v) {
            _tetGrad[t][v] = dfdJ * tet.ig[v];
        }
    }
    if (G != nullptr) {
#pragma omp parallel for
        for (unsigned v = 0; v < _vert2tetCorner.size();++v) {
            for (unsigned d = 0; d < 3; ++d) {
                (*G)(3*v+d) = 0;
                if (_compressedLocks[3*v+d]) continue;
                for (unsigned tc : _vert2tetCorner[v]) {
                    (*G)(3*v+d) += _tetGrad[tc/4][tc%4](d);
                }
            }
        }
    }

    return F;
}


