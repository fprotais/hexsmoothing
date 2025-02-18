#include "ultimaille/all.h"
#include "lib/elliptic_smoothing.h"
#include <utils/meshStructures.h>
#include <map>
#include <set>

#define INVERSE_ORIENTATION 0
#define SMOOTH_TETS 1
#define SMOOTH_HEXES 1
#define SMOOTH_PYRAMIDS 1
#define SMOOTH_WEDGES 1

#define FOR(i, n) for(int i = 0; i < n; i++)

using namespace UM;

bool smoothMixedMesh(
    std::vector<UM::vec3> &mixedVerts,
    std::vector<int> const &tets,
    std::vector<int> const &hexes,
    std::vector<int> const &wedges, 
    std::vector<int> const &pyramids,
    int maxIter
) {
    std::cout << "Creating proxy mesh..." << std::endl;

	utilities::TetrahedralMesh proxy_mesh;

    proxy_mesh._pts.resize(mixedVerts.size());
    FOR(i, mixedVerts.size()) FOR(d, 3) proxy_mesh._pts[i][d] = mixedVerts[i][d];


	proxy_mesh._tets.reserve(4 * tets.size() + 8 * hexes.size() + 6 * wedges.size() + 4 * pyramids.size());
	std::vector<std::array<utilities::vec3, 4>> refs;
	refs.reserve(4 * tets.size() + 8 * hexes.size() + 6 * wedges.size() + 4 * pyramids.size());

    std::map<std::set<int>, int> faceCnt;

#if SMOOTH_TETS
    { // tets
        std::cout << "copying tets. (" << tets.size() / 4 << ")" << std::endl;
        #if INVERSE_ORIENTATION
        const std::array<utilities::vec3,4> tet_ref = {{
            {std::sqrt(8./9),0,-1./3}, {-std::sqrt(2./9),std::sqrt(2./3),-1./3}, {0,0,1}, {-std::sqrt(2./9),-std::sqrt(2./3),-1./3}
        }};
        #else 
        const std::array<utilities::vec3,4> tet_ref = {{
            {std::sqrt(8./9),0,-1./3}, {-std::sqrt(2./9),std::sqrt(2./3),-1./3}, {-std::sqrt(2./9),-std::sqrt(2./3),-1./3}, {0,0,1}
        }};
        #endif
        
        FOR(c, tets.size() / 4) {
            std::array<int, 4> tet;
            FOR(j, 4) tet[j] =  tets[4* c + j];
            proxy_mesh._tets.push_back(tet);
            refs.push_back(tet_ref);
            static constexpr int facet_vertex[4][3] = {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}};
            FOR(i, 4) {
                std::set<int> tri;
                FOR(j, 3) tri.insert(tets[4* c + facet_vertex[i][j]]);
                if (faceCnt.find(tri) == faceCnt.end()) {
                    faceCnt.insert({tri, 0});
                }
                ++faceCnt[tri];
            }
        }

    }
#endif

#if SMOOTH_HEXES
    { // hexes
        std::cout << "copying hexes. (" <<  hexes.size() / 8 << ")" << std::endl;
        constexpr int HEX_CORNER_SPLITTING[8][4] = {
            {0,1,2,4}, {1,3,0,5}, {2,0,3,6}, {3,2,1,7},
            {4,6,5,0}, {5,4,7,1}, {6,7,4,2}, {7,5,6,3},
        };
        #if INVERSE_ORIENTATION
        const utilities::vec3 hex_ref[8] = {
            {0,0,0}, {0,1,0}, {1,0,0}, {1,1,0},
            {0,0,1}, {0,1,1}, {1,0,1}, {1,1,1},
        };
        #else 
        const utilities::vec3 hex_ref[8] = {
            {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
            {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1},
        };
        #endif 

        
        FOR(c, hexes.size() / 8) {
            FOR(i, 8) {
                std::array<int, 4> tet;
                std::array<utilities::vec3, 4> ref;
                FOR(j, 4) tet[j] = hexes[8* c + HEX_CORNER_SPLITTING[i][j]];
                FOR(j, 4) ref[j] = hex_ref[HEX_CORNER_SPLITTING[i][j]];
                proxy_mesh._tets.push_back(tet);
                refs.push_back(ref);
            }
            static constexpr int facet_vertex[6][4] = {{0,4,6,2}, {1,3,7,5}, {0,1,5,4}, {2,6,7,3}, {0,2,3,1}, {4,5,7,6}};
            FOR(i, 6) {
                std::set<int> quad;
                FOR(j, 4) quad.insert(hexes[8* c + facet_vertex[i][j]]);
                if (faceCnt.find(quad) == faceCnt.end()) {
                    faceCnt.insert({quad, 0});
                }
                ++faceCnt[quad];
            }
        }
    }
#endif

#if SMOOTH_PYRAMIDS
    { // pyramid
    // for definition, see volume.h in ultimaille
        std::cout << "copying pyramids. (" <<  pyramids.size() / 5 << ")" << std::endl;
        constexpr int pyramid_CORNER_SPLITTING[4][4] = {
            {0,1,3,4}, {1,2,3,4}, {2,3,1,4}, {3,0,2,4},
        };

        #if INVERSE_ORIENTATION
        const utilities::vec3 pyramid_ref[5] = {
            {0,0,0}, {0,1,0}, {1,1,0}, {1,0,0}, {0.5,0.5,0.5}
        };
        #else 
        const utilities::vec3 pyramid_ref[5] = {
            {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0.5,0.5,0.5}
        };        
        #endif
        FOR(c, pyramids.size() / 5) {
            FOR(i, 4) {
                std::array<int, 4> tet;
                std::array<utilities::vec3, 4> ref;
                FOR(j, 4) tet[j] = pyramids[5* c + pyramid_CORNER_SPLITTING[i][j]];
                FOR(j, 4) ref[j] = pyramid_ref[pyramid_CORNER_SPLITTING[i][j]];
                proxy_mesh._tets.push_back(tet);
                refs.push_back(ref);
            }
            static constexpr int facet_vertex[5][4] = {{0,3,2,1}, {0,1,4,-1}, {0,4,3,-1}, {2,3,4,-1}, {1,2,4,-1}};
            FOR(i, 5) {
                std::set<int> face;
                FOR(j, 4) if (facet_vertex[i][j] != -1) face.insert(pyramids[5* c + facet_vertex[i][j]]);
                if (faceCnt.find(face) == faceCnt.end()) {
                    faceCnt.insert({face, 0});
                }
                ++faceCnt[face];
            }
        }
    }
#endif 


#if SMOOTH_WEDGES
    { // wedge
        std::cout << "copying wedges. (" <<  wedges.size() / 6 << ")" << std::endl;
        constexpr int WEDGE_CORNER_SPLITTING[6][4] = {
            {0,1,2,3}, {1,2,0,4}, {2,0,1,5}, 
            {3,5,4,0}, {4,3,5,1}, {5,4,3,2}, 
        };

        #if INVERSE_ORIENTATION
        const utilities::vec3 wedge_ref[6] = {
            {0,0,0}, {1,0,0}, {0.5,std::sqrt(3)/2,0}, 
            {0,0,-1}, {1,0,-1}, {0.5,std::sqrt(3)/2,-1}, 
        };
        #else 
        const utilities::vec3 wedge_ref[6] = {
            {0,0,0}, {1,0,0}, {0.5,std::sqrt(3)/2,0}, 
            {0,0,1}, {1,0,1}, {0.5,std::sqrt(3)/2,1}, 
        };
        #endif
        FOR(c, wedges.size() / 6) {
            FOR(i, 6) {
                std::array<int, 4> tet;
                std::array<utilities::vec3, 4> ref;
                FOR(j, 4) tet[j] = wedges[6* c + WEDGE_CORNER_SPLITTING[i][j]];
                FOR(j, 4) ref[j] = wedge_ref[WEDGE_CORNER_SPLITTING[i][j]];
                proxy_mesh._tets.push_back(tet);
                refs.push_back(ref);
            }
            static constexpr int facet_vertex[5][4] = {{0,2,1,-1}, {3,4,5,-1}, {0,1,4,3}, {0,3,5,2}, {1,2,5,4}};
            FOR(i, 5) {
                std::set<int> face;
                FOR(j, 4) if (facet_vertex[i][j] != -1) face.insert(wedges[6* c + facet_vertex[i][j]]);
                if (faceCnt.find(face) == faceCnt.end()) {
                    faceCnt.insert({face, 0});
                }
                ++faceCnt[face];
            }
        }
    }
#endif 

	std::vector<double> verts(mixedVerts.size() * 3);
	FOR(v, mixedVerts.size()) FOR(d, 3) verts[3 * v + d] = mixedVerts[v][d];

    std::cout << "Initializing boundaries..." << std::endl;
    std::vector<bool> locks (mixedVerts.size() * 3, false);
    for (auto [s, cnt] : faceCnt) {
        if (cnt == 1) {
            for (int v : s) {
                FOR(d, 3) locks[3*v+d] = true;
            }
        }
    }

    //// debug quadratures
    // FOR(i, proxy_mesh._tets.size()) {
    //     std::cout << "===================================" << std::endl;
    //     std::cout << "REF:" << std::endl;
    //     FOR(j, 4) {
    //         std::cout << refs[i][j] << std::endl;
    //     }
    //     std::cout << cross(refs[i][1] - refs[i][0], refs[i][2] - refs[i][0]) * (refs[i][3] - refs[i][0]) << std::endl; 
    //     std::cout << "TET:" << std::endl;
    //     FOR(j, 4) {
    //         std::cout << proxy_mesh._tets[i][j] << ",";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "coords:" << std::endl;
    //     FOR(j, 4) {
    //         std::cout << proxy_mesh._pts[proxy_mesh._tets[i][j]] << std::endl;
    //     }
    //     utilities::vec3 v0 = proxy_mesh._pts[proxy_mesh._tets[i][1]] - proxy_mesh._pts[proxy_mesh._tets[i][0]];
    //     utilities::vec3 v1 = proxy_mesh._pts[proxy_mesh._tets[i][2]] - proxy_mesh._pts[proxy_mesh._tets[i][0]];
    //     utilities::vec3 v2 = proxy_mesh._pts[proxy_mesh._tets[i][3]] - proxy_mesh._pts[proxy_mesh._tets[i][0]];
    //     std::cout << cross(v0, v1) * v2 << std::endl; 
    // }

    // scaling refs to avoid epsilon issues
    double avgEdgeSize = 0;
    FOR(i, proxy_mesh._tets.size()) FOR(j, 4) {
            utilities::vec3 dir = proxy_mesh._pts[proxy_mesh._tets[i][j]] - proxy_mesh._pts[proxy_mesh._tets[i][0]];
            avgEdgeSize += dir.norm();
    }
    avgEdgeSize /= proxy_mesh._tets.size();
    FOR(i, proxy_mesh._tets.size()) FOR(j, 4) {
        refs[i][j] *= avgEdgeSize;
    }
    // FOR(v, mixedVerts.size()) {
    //     if (locks[3*v]) continue;
    //     FOR(d, 3) {
    //         double pert = (6. * std::rand()) / RAND_MAX;
    //         verts[3*v+d] += pert;
    //         mixedVerts[v][d] += pert;
    //     }
    // }
    // std::vector<int> empty;
    // UM::write_medit_format("pertub.mesh", mixedVerts, empty, empty, empty, tets, hexes, wedges, pyramids);


	Tets_id_with_lock var(verts, proxy_mesh._tets, locks);


	smoother_options options = _3D_default;
	options.static_threshold = 1e-7;
	options.bfgs_threshold = 1e-14;
	options.theta = 0.;
	options.bfgs_maxiter = 50;
	options.eps_from_theorem = true;
	options.maxiter = maxIter;
	options.debug = true;
    std::cout << "Smoothing..." << std::endl;
    std::cout << "Proxy tet mesh: " << proxy_mesh._tets.size() << " elements" << std::endl;
    std::cout << "                " << proxy_mesh._pts.size() << " vertices"<< std::endl;



	Elliptic_smoother_3D opt(var, proxy_mesh._tets.size(), options);
	FOR(i, refs.size()) opt.set_tet_ref(i, refs[i]);
    
	bool res = opt.go();
	var.get_verts(verts);
	FOR(v, mixedVerts.size()) FOR(d, 3) mixedVerts[v][d] = verts[3 * v + d];

    return res;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " hexmesh.ext smoothed_hexmesh.ext nbiter" << std::endl;
        std::cerr << "Input:" << std::endl;
        std::cerr << "- hexmesh.ext must contain a hexahedral mesh" << std::endl;
        std::cerr << "Output: " << std::endl;
        std::cerr << "- smoothed_hexmesh.ext, result smoothed mesh'" << std::endl;
        std::cerr << "Optional input: " << std::endl;
        std::cerr << "- nbiter: number of smoothing iter. Default is 20." << std::endl;
        std::cerr << std::endl;
        std::cerr << "ext formats are ultimaille supported volumic formats (geogram, medit -.mesh- and vtk)." << std::endl;
        std::cerr << "contact: francoisprotais@gmail.com" << std::endl;
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

    UM::read_vtk_format(inputMesh, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    std::cout << "MIXED MESH SMOOTHING" << std::endl;
    // UM::write_medit_format("input.mesh", verts, edges, tris, quads, tets, hexes, wedges, pyramids);

    bool res = smoothMixedMesh(verts, tets, hexes, wedges, pyramids, nbIter);

    std::cout << "Done: " << (res ? "success" : "failed untangling") << std::endl;
    UM::write_vtk_format(outputMesh, verts, edges, tris, quads, tets, hexes, wedges, pyramids);

    // UM::write_medit_format("output.mesh", verts, edges, tris, quads, tets, hexes, wedges, pyramids);

   
    return 0;
}

