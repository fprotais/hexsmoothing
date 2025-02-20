#include <ultimaille/all.h>
#include <lib/mesh_converter.h>
#include <lib/elliptic_smoothing.h>
#include <lib/hilbertMeshSort.h>

#include <parallel_lib/parallel_tet_smoother.h>

#include <Eigen/Eigen>
#include <omp.h>

#include <utils/logTime.h>
#include <random>


#define FOR(i, n) for(int i = 0; i < n; i++)

using namespace UM;


double getAvgEdgeSize(utilities::TetrahedralMesh const &mesh) {
    double avgEdgeSize = 0;
    FOR(i, mesh._tets.size()) FOR(j, 3) {
            utilities::vec3 e = mesh._pts[mesh._tets[i][j+1]] - mesh._pts[mesh._tets[i][0]];
            avgEdgeSize += e.norm();
    }
    return avgEdgeSize / (3*mesh._tets.size());
}

double det(utilities::TetrahedralMesh const &mesh, int t) {
    std::array<utilities::vec3,3> b;
    FOR(j, 3) {
        b[j] = mesh._pts[mesh._tets[t][j+1]] - mesh._pts[mesh._tets[t][0]];
    }
    return cross(b[0], b[1]) * b[2];
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " mixedMesh.ext smoothed_hexmesh.ext nbiter" << std::endl;
        std::cout << "- (input) hexmesh.ext must contain a mesh. It can include Tetrahedra, Hexahedra, Wedges and Pyramids." << std::endl;
        std::cout << "- (output) smoothed_hexmesh.ext, result smoothed mesh'" << std::endl;
        std::cout << "- (input - optional) nbiter: number of smoothing iter. Default is 20." << std::endl;
        std::cout << std::endl;
        std::cout << "ext formats can be '.mesh' -medit- or '.vtk'." << std::endl;
        std::cout << "contact: francoisprotais@gmail.com" << std::endl;
        return 1;
    }
    std::string inputMesh(argv[1]);
    std::string outputMesh(argv[2]);
    int nbIter = 20;
    if (argc > 3) {
        try {
            nbIter = std::stoi(argv[3]);
        }
        catch (std::invalid_argument const& ex) {}
        if (nbIter < 0) nbIter = 20;
    }

    std::vector<UM::vec3> verts;
    std::vector<int> edges;
    std::vector<int> tris; 
    std::vector<int> quads;
    std::vector<int> tets;
    std::vector<int> hexes;
    std::vector<int> wedges; 
    std::vector<int> pyramids;

    std::cout << "PARALLEL MIXED MESH SMOOTHING" << std::endl;
    std::cout << "Nb of available threads: " << omp_get_max_threads() << std::endl;
    std::cout << "Nb threads used by Eigen: " << Eigen::nbThreads() << std::endl;

    TimeLog logging("Parallel smoothing");
    UM::read_mixedMesh_byExtension(inputMesh, verts, edges, tris, quads, tets, hexes, wedges, pyramids);

    // UM::write_medit_format("input.mesh", verts, edges, tris, quads, tets, hexes, wedges, pyramids);

    

    logging.logSubStep("Reading mesh");
	utilities::TetrahedralMesh proxy_mesh;
	std::vector<std::array<utilities::vec3, 4>> refs;
    std::vector<bool> bndVert;

    convertedToTetMesh(verts, tets, hexes, wedges, pyramids, proxy_mesh, bndVert, refs);



    // for re-ordering
    std::vector<int> vert_old2new(proxy_mesh._pts.size());
    std::vector<int> tet_old2new(proxy_mesh._tets.size());
    std::iota(vert_old2new.begin(), vert_old2new.end(), 0.);
    std::iota(tet_old2new.begin(), tet_old2new.end(), 0.);

    constexpr bool doHilbertSort = false; // Should improve cache access. Did not see impact for meshes with 200k hexes. Maybe will have with larger? 

    if (doHilbertSort) { 
        hilbertSort(proxy_mesh, vert_old2new, tet_old2new);
        {
            std::vector<std::array<utilities::vec3, 4>> tmprefs = refs;
            for (unsigned t = 0; t < tet_old2new.size(); ++t) {
                tmprefs[tet_old2new[t]] = refs[t];
            }
            refs = tmprefs;
            std::vector<bool> tmpbndVert = bndVert;
            for (unsigned v = 0; v < vert_old2new.size(); ++v) {
                tmpbndVert[vert_old2new[v]] = bndVert[v];
            }
            bndVert = tmpbndVert;
        }
        logging.logSubStep("Hilbert sort");
    }



    double avgEdgeSize = getAvgEdgeSize(proxy_mesh);

    FOR(i, proxy_mesh._pts.size()) {
        // proxy_mesh._pts[i] /= avgEdgeSize;
    }
    // avgEdgeSize = getAvgEdgeSize(proxy_mesh);


    double pert = 100.;

    // std::uniform_real_distribution<float> distrib(-1., 1.);
    // std::mt19937 gen(42);
    // FOR(v, verts.size()) if (!bndVert[v]) proxy_mesh._pts[v] +=  pert * avgEdgeSize * utilities::vec3{distrib(gen), distrib(gen), distrib(gen)};

    FOR(i, verts.size()) FOR(d, 3) verts[i][d] = proxy_mesh._pts[vert_old2new[i]][d];

    // UM::write_medit_format("input_pert.mesh", verts, edges, tris, quads, tets, hexes, wedges, pyramids);

    logging.logSubStep("Convert to tet elements");

    double mindet = std::numeric_limits<double>::max();
    FOR(t, proxy_mesh._tets.size()) mindet = std::min(det(proxy_mesh, t), mindet);

    bool untangleFirst = (mindet < 0);
    std::cout << "Min det is: " << mindet << std::endl;

    if (untangleFirst) {
        std::cout << "===================================================================================" << std::endl;
        std::cout << "Starting by untangling mesh" << std::endl;
        std::cout << "===================================================================================" << std::endl;
        std::vector<bool> locks = bndVert;
        constexpr unsigned maxUntangleIter = 10;
        unsigned untangleIter = 1;
        unsigned nbVar = locks.size();
        while (mindet < 0) {
            std::cout << "Untangling iter: " << untangleIter << std::endl;
            Parallel_tet_smoother untangler(proxy_mesh);
            untangler.setRefs(refs);
            untangler.setVertsLocks(locks);
            untangler.max_untangling_iter =  untangleIter * untangleIter + 2;
            untangler.max_lbfgs_iter = 5000;
            untangler.verbose = nbVar > 5000;
            untangler.fineTimeLogging = (nbVar > 500000); // enable to see the progress on large meshes
            bool res = untangler.go();

            std::fill(locks.begin(), locks.end(), true); // on second iter untangling only flipped elements
            nbVar = 0;
            mindet = std::numeric_limits<double>::max();
            FOR(t, proxy_mesh._tets.size()) {
                double tetdet = det(proxy_mesh, t);
                mindet = std::min(tetdet, mindet);
                if (tetdet > 0) continue;
                FOR(tv, 4) if (!bndVert[proxy_mesh._tets[t][tv]] && locks[proxy_mesh._tets[t][tv]]) {
                    ++nbVar;
                    locks[proxy_mesh._tets[t][tv]] = false;
                }
            }
            std::cout << "mindet: " << mindet << std::endl;

            if (mindet > 0) break;

            if (maxUntangleIter == untangleIter++) break;
        }

        FOR(i, verts.size()) FOR(d, 3) verts[i][d] = proxy_mesh._pts[vert_old2new[i]][d];
        
        std::cout << "Pre-saving mesh: " << outputMesh << std::endl;
        UM::write_mixedMesh_byExtension(outputMesh, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        std::cout << "===================================================================================" << std::endl;
        std::cout << "Untangling done: " << ((mindet > 0) ? "Success" : "Fail") << std::endl;
        std::cout << "===================================================================================" << std::endl;
        logging.logSubStep("Untangling");
    }

    Parallel_tet_smoother smoother(proxy_mesh);
    smoother.setRefs(refs);
    smoother.setVertsLocks(bndVert);
    smoother.max_untangling_iter = nbIter;
    smoother.max_lbfgs_iter = 50;
    smoother.fineTimeLogging = (proxy_mesh._tets.size() > 3000000); // enable to see the progress on large meshes
    bool res = smoother.go();

    FOR(i, verts.size()) FOR(d, 3) verts[i][d] = proxy_mesh._pts[vert_old2new[i]][d];
    logging.logSubStep("Smoothing");
    

    std::cout << "Done: " << (res ? "success" : "failed untangling") << std::endl;
    UM::write_mixedMesh_byExtension(outputMesh, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    logging.logSubStep("Saving mesh");
    logging.logTotalTime();


   
    return 0;
}

