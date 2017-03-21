#ifndef LIQUIDITY_PLANE_HDR
#define LIQUIDITY_PLANE_HDR

#include "coor_tuple.hpp"
#include "global_consts.hpp"

namespace math {
	struct t_plane {
	public:
		t_plane(float nx = 0.0f, float ny = 0.0f, float nz = 0.0f, float d = 0.0f): eq({nx, ny, nz, d}) {}

		void collision(sph::t_particle& pi) {
			const t_vec3f& n = t_vec3f({eq.x(), eq.y(), eq.z()});

			const float center_dist = (n * pi.cur_pos) - eq.w();
			const float radius_dist = center_dist - consts::PARTICLE_RADIUS;
			const float collis_mask = (radius_dist < 0.0f);

			pi.cur_pos += (n * radius_dist * -1.0f * collis_mask); // clamp
			pi.cur_vel -= (n * ((n * pi.cur_vel) * 2.0f * collis_mask)); // reflect
		}

	public:
		t_vec4f eq;
	};

	#if 0
	struct t_edge {
		t_vec3f p0;
		t_vec3f p1;
		t_vec3f n;
	};
	#endif
};

#endif

