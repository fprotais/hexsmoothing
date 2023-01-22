#include <utils/meshStructures.h>



bool smooth_tet_mesh(utilities::TetrahedralMesh& m, std::vector<bool>& locks);
bool smooth_tet_mesh(utilities::TetrahedralMesh& m, const utilities::TetrahedralMesh& ref, std::vector<bool>& locks);

bool smooth_hex_mesh(utilities::HexahedralMesh& m, const std::vector<bool>& locks, unsigned max_iter = 50, bool only_SJ_tets = true);
