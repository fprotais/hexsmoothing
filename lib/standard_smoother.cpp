#include "standard_smoother.h"
#include "elliptic_smoothing.h"

#include <assert.h>

using utilities::vec3;

using utilities::TetrahedralMesh;
using utilities::HexahedralMesh;

#define FOR(i, n) for(int i = 0; i < n; i++)

bool smooth_tet_mesh(TetrahedralMesh& m, std::vector<bool>& locks) {
	smoother_options options = _3D_default;
	options.static_threshold = 1e-7;
	options.bfgs_threshold = 1e-14;
	options.theta = 1e-3;
	options.bfgs_maxiter = 1000;
	options.eps_from_theorem = false;
	options.maxiter = 100;
	options.debug = false;

    assert(m._pts.size() * 3 == locks.size());

    std::vector<double> verts(m._pts.size() * 3);
    FOR(v, m._pts.size()) FOR(d, 3) verts[3 * v + d] = m._pts[v][d];
    std::vector<std::array<int, 4>> tets(m._tets.size());
    FOR(f, m._tets.size()) FOR(fv, 4) tets[f][fv] = m._tets[f][fv];
    Tets_id_with_lock var(verts, tets, locks);

    Elliptic_smoother_3D opt(var, tets.size(), options);
    bool res = opt.go();

    var.get_verts(verts);
    FOR(v, m._pts.size()) FOR(d, 3) m._pts[v][d] = verts[3 * v + d];
    return res;
}
bool smooth_tet_mesh(TetrahedralMesh& m, const TetrahedralMesh& ref, std::vector<bool>& locks){
	smoother_options options = _3D_default;
	options.static_threshold = 1e-7;
	options.bfgs_threshold = 1e-14;
	options.theta = 1e-3;
	options.bfgs_maxiter = 1000;
	options.eps_from_theorem = false;
	options.maxiter = 100;
	options.debug = false;
	
	assert(m._pts.size() == ref._pts.size());
    assert(m._pts.size() * 3 == locks.size());

    std::vector<double> verts(m._pts.size() * 3);
    FOR(v, m._pts.size()) FOR(d, 3) verts[3 * v + d] = m._pts[v][d];
    std::vector<std::array<int, 4>> tets(m._tets.size());
    FOR(f, m._tets.size()) FOR(fv, 4) tets[f][fv] = m._tets[f][fv];
    Tets_id_with_lock var(verts, tets, locks);

    Elliptic_smoother_3D opt(var, tets.size(), options);
    FOR(c, tets.size()) {
        std::array<vec3, 4> tet_ref;
        FOR(cv, 4) tet_ref[cv] = ref._pts[tets[c][cv]];
        opt.set_tet_ref(c, tet_ref);
    }
    bool res = opt.go();

    var.get_verts(verts);
    FOR(v, m._pts.size()) FOR(d, 3) m._pts[v][d] = verts[3 * v + d];
    return res;
}

constexpr int HEX_CORNER_SPLITTING[8][4] = {
	{0,1,2,4}, {1,3,0,5}, {2,0,3,6}, {3,2,1,7},
	{4,6,5,0}, {5,4,7,1}, {6,7,4,2}, {7,5,6,3},
};

constexpr int TET_MATRIX_HEX_LEXICOGRAPHIC_SPLIT[24][4] = {
    {0,1,3,4},{0,1,2,5},{0,2,6,1},{0,2,4,3},
    {0,4,1,6},{0,4,5,2},{1,3,0,7},{1,3,2,5},
    {1,5,3,4},{1,5,7,0},{2,3,7,0},{2,3,6,1},
    {2,6,0,7},{2,6,4,3},{3,7,6,1},{3,7,2,5},
    {4,5,0,7},{4,5,1,6},{4,6,5,2},{4,6,7,0},
    {5,7,1,6},{5,7,3,4},{6,7,5,2},{6,7,4,3}
};

static const vec3 hex_ref[8] = {
	{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
	{0,0,1}, {1,0,1}, {0,1,1}, {1,1,1},
};

bool smooth_hex_mesh(HexahedralMesh& m, const std::vector<bool>& locks, unsigned max_iter, bool only_SJ_tets) {
	assert(m._pts.size() * 3 == locks.size());
	smoother_options options = _3D_default;
	options.static_threshold = 1e-7;
	options.bfgs_threshold = 1e-14;
	options.theta = 1e-3;
	options.bfgs_maxiter = 100;
	options.eps_from_theorem = false;
	options.maxiter = max_iter;
	options.debug = true;

	std::vector<double> verts(m._pts.size() * 3);
	FOR(v, m._pts.size()) FOR(d, 3) verts[3 * v + d] = m._pts[v][d];

	TetrahedralMesh proxy_mesh;
	proxy_mesh._pts = m._pts;
	proxy_mesh._tets.reserve(only_SJ_tets ? 8 : 32);
	std::vector<std::array<vec3, 4>> refs;
	refs.reserve(only_SJ_tets ? 8 : 32);


	std::vector<int> map(m._pts.size(), -1);
	const std::array<vec3, 8> default_cell = { { {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {1.,1.,0.}, {0.,0.,1.}, {1.,0.,1.}, {0.,1.,1.}, {1.,1.,1.} } };
	FOR(c, m._hexes.size()) {
		std::array<bool, 8> has_on_in = { 0 };
		int start_c = proxy_mesh._tets.size();
		FOR(i, 8) {
			bool has_in = false;
			FOR(j, 4) FOR(d, 3) if (!locks[3 * m._hexes[c][HEX_CORNER_SPLITTING[i][j]] + d]) has_in = true;
			if (!has_in) continue;
			std::array<int, 4> tet;
			std::array<vec3, 4> ref;
			FOR(j, 4) tet[j] = m._hexes[c][HEX_CORNER_SPLITTING[i][j]];
			FOR(j, 4) ref[j] = hex_ref[HEX_CORNER_SPLITTING[i][j]];
			proxy_mesh._tets.push_back(tet);
			refs.push_back(ref);
		}
		if (only_SJ_tets) continue;

		FOR(i, 8) {
			bool has_in = false;
			FOR(j, 4) FOR(d, 3) if (!locks[3 * m._hexes[c][TET_MATRIX_HEX_LEXICOGRAPHIC_SPLIT[i][j]] + d]) has_in = true;
			if (!has_in) continue;
			std::array<int, 4> tet;
			std::array<vec3, 4> ref;
			FOR(j, 4) tet[j] = m._hexes[c][TET_MATRIX_HEX_LEXICOGRAPHIC_SPLIT[i][j]];
			FOR(j, 4) ref[j] = hex_ref[TET_MATRIX_HEX_LEXICOGRAPHIC_SPLIT[i][j]];
			proxy_mesh._tets.push_back(tet);
			refs.push_back(ref);
		}

	}

	Tets_id_with_lock var(verts, proxy_mesh._tets, locks);

	Elliptic_smoother_3D opt(var, proxy_mesh._tets.size(), options);
	FOR(i, refs.size()) opt.set_tet_ref(i, refs[i]);

	bool res = opt.go();
	var.get_verts(verts);
	FOR(v, m._pts.size()) FOR(d, 3) m._pts[v][d] = verts[3 * v + d];
	std::cout << "Done\n";
	return res;
}
