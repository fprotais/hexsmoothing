#pragma once

#include <vector>
#include <array>
#include "vec.h"

// numbering convention is the same as https://github.com/ssloy/ultimaille
namespace utilities {

    struct CurveMesh {
        std::vector<vec3> _pts;
        std::vector<std::array<int, 2>> _edges;
    };

    // the 3 points form a 2D direct base with tri[0] the orign
    struct TriangleMesh {
        std::vector<vec3> _pts;
        std::vector<std::array<int, 3>> _tris;
    };

    // the 4 points form a direct base with tet[0] the origin
    struct TetrahedralMesh {
        std::vector<vec3> _pts;
        std::vector<std::array<int, 4>> _tets;
    };

       /* 
        *           6-----------7                                               
        *          /|          /|                                           
        *         / |         / |                                            
        *        /  |        /  |                                               
        *       4-----------5   |                                             
        *       |   |       |   |                                               
        *       |   2-------|---3                                               
        *       |  /        |  /                                                 
        * Z     | /         | /                                            
        * ^  Y  |/          |/                                             
        * | /   0-----------1                                              
        * |/                                                               
        * o----> X   
        */   
    struct HexahedralMesh {
        std::vector<vec3> _pts;
        std::vector<std::array<int, 8>> _hexes;
    };

    static constexpr std::array<std::array<int, 4>, 8> HEX_TO_SJ_TET_SPLIT = { {
        {0,1,2,4},
        {1,3,0,5},
        {2,0,3,6},
        {3,2,1,7},
        {4,6,5,0},
        {5,4,7,1},
        {6,7,4,2},
        {7,5,6,3},
    } };    

}

