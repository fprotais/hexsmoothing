#ifndef __MEDIT_H__
#define __MEDIT_H__

#include <cstring>
#include "ultimaille/attributes.h"
#include "ultimaille/volume.h"
#include "ultimaille/surface.h"
#include "ultimaille/polyline.h"

// TODO: colors and scalar fields are yet to be implemented

/**
 * N.B.: Medit local numbering convention differs, thus the vertices are reordered on import/export
 *
 *     ultimaille:                Medit:
 *     6--------7               7--------6
 *    /|       /|              /|       /|
 *   / |      / |             / |      / |
 *  4--------5  |            4--------5  |
 *  |  |     |  |            |  |     |  |
 *  |  2-----|--3            |  3-----|--2
 *  | /      | /             | /      | /
 *  |/       |/              |/       |/
 *  0--------1               0--------1
 */

namespace UM {
    void write_medit_format(const std::string& filename, const std::vector<vec3>& verts_, const std::vector<int>& edges_, const std::vector<int>& tris_, const std::vector<int>& quads_, const std::vector<int>& tets_, const  std::vector<int>& hexes_, const std::vector<int>& wedges_, const  std::vector<int>& pyramids_);
    void read_medit_format(const std::string& filename, std::vector<int>& vcolors_, std::vector<vec3>& verts_, std::vector<int>& edges_, std::vector<int>& tris_, std::vector<int>& quads_, std::vector<int>& tets_, std::vector<int>& hexes_, std::vector<int>& wedges_, std::vector<int>& pyramids_);

    void write_medit(const std::string filename, const PolyLine &pl);
    void write_medit(const std::string filename, const Surface &m);
    void write_medit(const std::string filename, const Volume  &m);

    PolyLineAttributes read_medit(const std::string filename, PolyLine   &m);
    SurfaceAttributes  read_medit(const std::string filename, Triangles  &m);
    SurfaceAttributes  read_medit(const std::string filename, Quads      &m);
    SurfaceAttributes  read_medit(const std::string filename, Polygons   &m);
    VolumeAttributes   read_medit(const std::string filename, Tetrahedra &m);
    VolumeAttributes   read_medit(const std::string filename, Hexahedra  &m); // cf. N.B.
    VolumeAttributes   read_medit(const std::string filename, Wedges     &m);
    VolumeAttributes   read_medit(const std::string filename, Pyramids   &m);
}

#endif // __MEDIT_H__

