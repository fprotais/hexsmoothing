#include <ultimaille/all.h>
#include <lib/mesh_converter.h>
#include <lib/elliptic_smoothing.h>
#include <lib/hilbertMeshSort.h>

#include <parallel_lib/parallel_tet_smoother.h>

#include <Eigen/Eigen>
#include <omp.h>

#include <utils/logTime.h>


#define FOR(i, n) for(int i = 0; i < n; i++)

using namespace UM;

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

    UM::write_medit_format("inputPar.mesh", verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    
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



    double avgEdgeSize = 0;
    FOR(i, proxy_mesh._tets.size()) FOR(j, 3) {
            utilities::vec3 dir = proxy_mesh._pts[proxy_mesh._tets[i][j+1]] - proxy_mesh._pts[proxy_mesh._tets[i][0]];
            avgEdgeSize += dir.norm();
    }
    avgEdgeSize /= 3*proxy_mesh._tets.size();
    FOR(i, proxy_mesh._tets.size()) FOR(j, 4) {
        refs[i][j] *= avgEdgeSize;
    }
    logging.logSubStep("Convert to tet elements");

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

    // UM::write_medit_format("outputPar.mesh", verts, edges, tris, quads, tets, hexes, wedges, pyramids);

   
    return 0;
}

