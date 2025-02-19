#include "mesh_converter.h"

#include <map>
#include <set>

#define INVERSE_ORIENTATION 0
#define CONVERT_TETS 1
#define CONVERT_HEXES 1
#define CONVERT_PYRAMIDS 1
#define CONVERT_WEDGES 1

#define FOR(i, n) for(int i = 0; i < n; i++)

void convertedToTetMesh(
    std::vector<UM::vec3> const &verts,
    std::vector<int> const &tets,
    std::vector<int> const &hexes,
    std::vector<int> const &wedges, 
    std::vector<int> const &pyramids,
    utilities::TetrahedralMesh &proxy_mesh, 
    std::vector<bool> &bndVerts,
    std::vector<std::array<utilities::vec3, 4>> &idealUnitShape
) {
    std::cout << "Converting Mixed mesh to tets " << std::endl;
    proxy_mesh._pts.clear();
    proxy_mesh._tets.clear();
    idealUnitShape.clear();
    
    proxy_mesh._pts.resize(verts.size());
    FOR(i, verts.size()) FOR(d, 3) proxy_mesh._pts[i][d] = verts[i][d];

    unsigned nbTets =  CONVERT_TETS * (tets.size() / 4) 
                     + CONVERT_HEXES * hexes.size()
                     + CONVERT_PYRAMIDS * 4 * (pyramids.size() / 5)
                     + CONVERT_WEDGES * wedges.size();
    proxy_mesh._tets.reserve(nbTets);
	idealUnitShape.reserve(nbTets);

    std::map<std::set<int>, int> faceCnt;


    #if CONVERT_TETS
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
            idealUnitShape.push_back(tet_ref);
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

#if CONVERT_HEXES
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
                idealUnitShape.push_back(ref);
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

#if CONVERT_PYRAMIDS
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
                idealUnitShape.push_back(ref);
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


#if CONVERT_WEDGES
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
                idealUnitShape.push_back(ref);
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
    std::cout << "Initializing boundaries..." << std::endl;
    bndVerts.clear();
    bndVerts.resize(verts.size(), false);
    for (auto [s, cnt] : faceCnt) {
        if (cnt == 1) {
            for (int v : s) {
                bndVerts[v] = true;
            }
        }
    }
}
