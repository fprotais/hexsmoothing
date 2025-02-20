#include "hilbertMeshSort.h"
#include <numeric>
#include "ultimaille/helpers/hilbert_sort.h"

void hilbertSort(utilities::TetrahedralMesh &mesh, std::vector<int> &vert_old2new, std::vector<int> &tet_old2new) {
    vert_old2new.resize(mesh._pts.size());
    tet_old2new.resize(mesh._tets.size());
    
    // vertices
    {
        std::vector<UM::vec3> umPoints(mesh._pts.size());
        for (unsigned v = 0; v < mesh._pts.size(); ++v) {
            for (unsigned d = 0; d < 3; ++d) {
                umPoints[v][d] = mesh._pts[v][d];
            }
        }
        std::vector<int> new2old(mesh._pts.size());
        std::iota(new2old.begin(), new2old.end(), 0.);
        UM::HilbertSort sort(umPoints);
        sort.apply(new2old);
        for (unsigned v = 0; v < mesh._pts.size(); ++v) {            
            vert_old2new[new2old[v]] = v;
        }    
        std::vector<utilities::vec3> tmp_pts = mesh._pts; 
        for (unsigned v = 0; v < mesh._pts.size(); ++v) {            
            tmp_pts[vert_old2new[v]] = mesh._pts[v];
        }    
        mesh._pts = tmp_pts;
        for (unsigned t = 0; t < mesh._tets.size(); ++t) {
            for (unsigned tv = 0; tv < 4; ++tv) {
                mesh._tets[t][tv] = vert_old2new[mesh._tets[t][tv]];
            }
        }
    }

    // tets
    {
        std::vector<UM::vec3> umPoints(mesh._tets.size());
        for (unsigned t = 0; t < mesh._tets.size(); ++t) {
            for (unsigned tv = 0; tv < 4; ++tv) {
                for (unsigned d = 0; d < 3; ++d) {
                    umPoints[t][d] = mesh._pts[mesh._tets[t][tv]][d];
                }
            }
            umPoints[t] /= 4;
        }
        std::vector<int> new2old(mesh._tets.size());
        std::iota(new2old.begin(), new2old.end(), 0.);
        UM::HilbertSort sort(umPoints);
        sort.apply(new2old);
        for (unsigned t = 0; t < mesh._tets.size(); ++t) {            
            tet_old2new[new2old[t]] = t;
        }  
        std::vector<std::array<int, 4>> tmp_tets = mesh._tets; 
        for (unsigned t = 0; t < mesh._tets.size(); ++t) {
            tmp_tets[tet_old2new[t]] = mesh._tets[t];
        }    
        mesh._tets = tmp_tets;
    }
}
