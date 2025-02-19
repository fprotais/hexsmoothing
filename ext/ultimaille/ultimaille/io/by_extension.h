#ifndef __BY_EXTENSION_H__
#define __BY_EXTENSION_H__

#include <filesystem>
#include <type_traits>

#include "ultimaille/attributes.h"
#include "ultimaille/io/geogram.h"
#include "ultimaille/io/medit.h"
#include "ultimaille/io/vtk.h"
#include "ultimaille/io/obj.h"

namespace UM {
    inline PolyLineAttributes empty_attr(const PolyLine &m) {
        (void)m;
        return {};
    }

    inline SurfaceAttributes empty_attr(const Surface &m) {
        (void)m;
        return {};
    }

    inline VolumeAttributes empty_attr(const Volume &m) {
        (void)m;
        return {};
    }

    template <class M> void write_by_extension(const std::string path, const M &m, const decltype(empty_attr(m)) a = {}) {
        std::string ext = std::filesystem::path(path).extension().string();
        if (ext == ".geogram")
            write_geogram(path, m, a);
        if (ext == ".mesh")
            write_medit(path, m);
        if (ext == ".vtk")
            write_vtk(path, m);
        if constexpr (std::is_same_v<decltype(empty_attr(m)), SurfaceAttributes>) {
            if (ext == ".obj")
                write_wavefront_obj(path, m, a);
        }
    }

    template <class M> auto read_by_extension(const std::string path, M& m) -> decltype(empty_attr(m)) {
        std::string ext = std::filesystem::path(path).extension().string();
        if (ext == ".geogram")
            return read_geogram(path, m);
        if (ext == ".mesh")
            return read_medit(path, m);
        if (ext == ".vtk")
            return read_vtk(path, m);
        if constexpr (std::is_same_v<decltype(empty_attr(m)), SurfaceAttributes>) {
            if (ext == ".obj")
                return read_wavefront_obj(path, m);
        }
        return {};
    }


    inline void write_mixedMesh_byExtension(
        const std::string& filename, 
        const std::vector<vec3>& verts, 
        const std::vector<int>& edges, 
        const std::vector<int>& tris, 
        const std::vector<int>& quads, 
        const std::vector<int>& tets, 
        const std::vector<int>& hexes, 
        const std::vector<int>& wedges, 
        const  std::vector<int>& pyramids
    ) {
        std::string ext = std::filesystem::path(filename).extension().string();
        if (ext == ".vtk") {
            UM::write_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        }
        else if (ext == ".mesh") {
            UM::write_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        }
        else {
            throw std::runtime_error("Error: unsupported extension for writing mixed mesh: '" + ext+ "'. Use only '.mesh' (medit) or '.vtk'");
        }
    }

    inline void read_mixedMesh_byExtension(
        const std::string& filename, 
        std::vector<vec3>& verts, 
        std::vector<int>& edges, 
        std::vector<int>& tris, 
        std::vector<int>& quads, 
        std::vector<int>& tets, 
        std::vector<int>& hexes, 
        std::vector<int>& wedges, 
        std::vector<int>& pyramids
    ) {
        std::string ext = std::filesystem::path(filename).extension().string();
        std::vector<int> color;
        if (ext == ".vtk") {
            UM::read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        }
        else if (ext == ".mesh") {
            UM::read_medit_format(filename, color, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        }
        else {
            throw std::runtime_error("Error: unsupported extension for reading mixed mesh: '" + ext+ "'. Use only '.mesh' (medit) or '.vtk'");
        }
    }
}

#endif // __BY_EXTENSION_H__

