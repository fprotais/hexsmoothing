#pragma once
#include <utils/meshStructures.h>
#include <Eigen/Eigen>

static const std::array<utilities::vec3, 4> equilateralTetRef = {{{std::sqrt(8./9),0,-1./3}, {-std::sqrt(2./9),std::sqrt(2./3),-1./3}, {0,0,1}, {-std::sqrt(2./9),-std::sqrt(2./3),-1./3}}};

class Parallel_tet_smoother {
public:   
    Parallel_tet_smoother(utilities::TetrahedralMesh &mesh);

    void setRefs(std::vector<std::array<utilities::vec3, 4>> const &refs);
    void setRef(unsigned tet, std::array<utilities::vec3, 4> const &ref = equilateralTetRef);

    void setVertsDimLocks(std::vector<bool> const &vertDimLock);
    void setVertDimLocks(unsigned vert, unsigned dim, bool lock = true);
    void setVertsLocks(std::vector<bool> const &vertLock);
    void setVertLock(unsigned vert, bool lock = true);
    
    unsigned max_untangling_iter = 50;
    unsigned max_lbfgs_iter = 100;

    bool verbose = true;
    bool fineTimeLogging = false;

    bool go();

private:

    bool vertIsLocked(unsigned v) const;

    void checkRefs();

    void updateEps(double rate);
    double _eps = 1e-12;
    double _eps0 = 1e-5; 
    double _noUntanglingEps = 1e-14;
    unsigned _nbInverted = 0;
    double _detMin = std::numeric_limits<double>::max();

    utilities::TetrahedralMesh &_mesh;
    std::vector<std::array<utilities::vec3, 4>> _refs;
    std::vector<bool> _vertDimLock;

    void createCompressSortedData();
    void hilbertSort();
    Eigen::VectorXd _compressedCoord;
    std::vector<bool> _compressedLocks;
    std::vector<std::vector<unsigned>> _vert2tetCorner;
    std::vector<int> _originalToCompressed;

    struct TetStorage {
        std::array<unsigned, 4> verts;
        std::array<Eigen::Vector3d, 4> ig;
    };
    std::vector<TetStorage> _tetData;

    void updateVertCoord();

    std::vector<std::array<Eigen::Vector3d, 4>> _tetGrad; // used inside elliptic_energy for parallel
    double elliptic_energy(Eigen::VectorXd const &x, Eigen::VectorXd *G = nullptr);

};


