#include "vec.h"
#include <array>

namespace utilities {

    double point_segment_squared_distance(const vec3& P, const std::array<vec3, 2>& AB, vec3& closest_point, std::array<double, 2>& l);

    double point_triangle_squared_distance(const vec3& P, const std::array<vec3, 3>& ABC, vec3& closest_point, std::array<double, 3>& l);

}
