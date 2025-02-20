#include "ultimaille/all.h"
#include <lib/mesh_converter.h>
#include "lib/elliptic_smoothing.h"
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

    TimeLog logging("Mixed smoothing");

    UM::read_mixedMesh_byExtension(inputMesh, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    std::cout << "MIXED MESH SMOOTHING" << std::endl;
    // UM::write_medit_format("input.mesh", verts, edges, tris, quads, tets, hexes, wedges, pyramids);


	utilities::TetrahedralMesh proxy_mesh;
	std::vector<std::array<utilities::vec3, 4>> refs;
    std::vector<bool> bndVert;

    convertedToTetMesh(verts, tets, hexes, wedges, pyramids, proxy_mesh, bndVert, refs);

    std::vector<double> serializedVerts(verts.size() * 3);
    std::vector<bool> vertLocks(verts.size() * 3);

	FOR(v, verts.size()) FOR(d, 3) serializedVerts[3 * v + d] = verts[v][d];
	FOR(v, verts.size()) FOR(d, 3) vertLocks[3 * v + d] = bndVert[v];

    double avgEdgeSize = 0;
    FOR(i, proxy_mesh._tets.size()) FOR(j, 3) {
            utilities::vec3 dir = proxy_mesh._pts[proxy_mesh._tets[i][j+1]] - proxy_mesh._pts[proxy_mesh._tets[i][0]];
            avgEdgeSize += dir.norm();
    }
    avgEdgeSize /= 3*proxy_mesh._tets.size();
    FOR(i, proxy_mesh._tets.size()) FOR(j, 4) {
        refs[i][j] *= avgEdgeSize;
    }

	Tets_id_with_lock var(serializedVerts, proxy_mesh._tets, vertLocks);


	smoother_options options = _3D_default;
	options.static_threshold = 1e-7;
	options.bfgs_threshold = 1e-14;
	options.theta = 0.;
	options.bfgs_maxiter = 500;
	options.eps_from_theorem = true;
	options.maxiter = nbIter;
	options.debug = true;
    std::cout << "Smoothing..." << std::endl;
    std::cout << "Proxy tet mesh: " << proxy_mesh._tets.size() << " elements" << std::endl;
    std::cout << "                " << proxy_mesh._pts.size() << " vertices"<< std::endl;


	Elliptic_smoother_3D opt(var, proxy_mesh._tets.size(), options);
	FOR(i, refs.size()) opt.set_tet_ref(i, refs[i]);
    
	bool res = opt.go();
	var.get_verts(serializedVerts);
	FOR(v, verts.size()) FOR(d, 3) verts[v][d] = serializedVerts[3 * v + d];


    std::cout << "Done: " << (res ? "success" : "failed untangling") << std::endl;
    UM::write_mixedMesh_byExtension(outputMesh, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    logging.logTotalTime();

    // UM::write_medit_format("output.mesh", verts, edges, tris, quads, tets, hexes, wedges, pyramids);

   
    return 0;
}

