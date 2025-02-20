#include "ultimaille/all.h"
#include "lib/standard_smoother.h"


#define FOR(i, n) for(int i = 0; i < n; i++)

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " hexmesh.ext smoothed_hexmesh.ext nbiter" << std::endl;
        std::cout << "Input:" << std::endl;
        std::cout << "- hexmesh.ext must contain a hexahedral mesh" << std::endl;
        std::cout << "Output: " << std::endl;
        std::cout << "- smoothed_hexmesh.ext, result smoothed mesh'" << std::endl;
        std::cout << "Optional input: " << std::endl;
        std::cout << "- nbiter: number of smoothing iter. Default is 50." << std::endl;
        std::cout << std::endl;
        std::cout << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cout << "contact: francoisprotais@gmail.com" << std::endl;
        return 1;
    }
    std::string inputMesh(argv[1]);
    std::string outputMesh(argv[2]);
    int nbIter = 50;
    if (argc > 3) {
        try {
            nbIter = std::stoi(argv[3]);
        }
        catch (std::invalid_argument const& ex) {}
        if (nbIter < 0) nbIter = 50;
    }


    UM::Hexahedra hexMesh;
    UM::read_by_extension(inputMesh, hexMesh);


	std::vector<bool> locks(hexMesh.nverts() * 3, false);
    UM::VolumeConnectivity vec(hexMesh);
    for (int c : UM::range(hexMesh.ncells())) {
        for (int cf : UM::range(6)) {
            if (vec.adjacent[6 * c + cf] != -1) continue;
            for (int cfv : UM::range(4)) {
                FOR(i, 3) locks[3*hexMesh.facet_vert(c,cf,cfv) + i] = true;
            }
        }
    }

    utilities::HexahedralMesh m;
    m._pts.resize(hexMesh.points.size());
    FOR(i, m._pts.size()) FOR(d, 3) m._pts[i][d] = hexMesh.points[i][d];
    m._hexes.resize(hexMesh.ncells());
    FOR(i, m._hexes.size()) FOR(j, 8) m._hexes[i][j] = hexMesh.vert(i, j);
    std::cout << "Nb verts: " << m._pts.size() << std::endl;
    std::cout << "Nb hexes: " << m._hexes.size() << std::endl;
    std::cout << "Nb iter planned: " << nbIter << std::endl;
    std::cout << "Now smoothing:" << std::endl;
    smooth_hex_mesh(m, locks, nbIter, true);
    std::cout << "Done." << std::endl;

    FOR(i, m._pts.size()) FOR(d, 3) hexMesh.points[i][d] = m._pts[i][d];
    UM::write_by_extension(outputMesh, hexMesh);
    return 0;
}

