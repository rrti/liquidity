#ifndef LIQUIDITY_MATERIAL_HDR
#define LIQUIDITY_MATERIAL_HDR

#include "coor_tuple.hpp"
#include "color.hpp"

namespace sph {
	struct t_material {
	public:
		t_material() {}
		t_material(const t_vec3f& _params, const t_color& _color): params(_params), color(_color) {}

		float  mass() const { return (params.x()); }
		float scale() const { return (params.y()); }
		float  bias() const { return (params.z()); }

	public:
		t_vec3f params; // .x := mass, .y := scale, .z := bias
		t_color color;
	};
};

#endif

