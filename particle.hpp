#ifndef LIQUIDITY_PARTICLE_HDR
#define LIQUIDITY_PARTICLE_HDR

#include "coor_tuple.hpp"

namespace sph {
	struct t_particle {
		t_vec3f cur_pos;
		t_vec3f cur_vel;

		t_vec3f prv_pos;
		t_vec3f rlx_pos;

		// pressures
		float rest_press;
		float near_press;

		// densities
		float rest_dens;
		float near_dens;

		// mass, reciprocal mass
		float raw_mass;
		float inv_mass;

		uint32_t mat_index;
	};
};

#endif

