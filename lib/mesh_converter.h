#pragma once

#include <utils/meshStructures.h>
#include "ultimaille/algebra/vec.h"

void convertedToTetMesh(
    std::vector<UM::vec3> const &verts,
    std::vector<int> const &tets,
    std::vector<int> const &hexes,
    std::vector<int> const &wedges, 
    std::vector<int> const &pyramids,
    utilities::TetrahedralMesh &proxy_mesh, 
    std::vector<bool> &bndVerts,
    std::vector<std::array<utilities::vec3, 4>> &idealUnitShape
);