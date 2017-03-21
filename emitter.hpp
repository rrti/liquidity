#ifndef LIQUIDITY_EMITTER_HDR
#define LIQUIDITY_EMITTER_HDR

#include "coor_tuple.hpp"
#include "global_consts.hpp"
#include "global_types.hpp"

namespace sph {
	struct t_particle;
	struct t_emitter {
	public:
		t_emitter(
			uint32_t mat_idx = 0,
			const t_vec3f& _pos = t_vec3f({0.0f, 0.0f, 0.0f}),
			const t_vec3f& _dir = t_vec3f({1.0f, 0.0f, 0.0f}),
			const t_vec3f& _rsd = t_vec3f({1.0f, 1.0f, 1.0f})
		): pos(_pos), dir(_dir), rsd(_rsd), mat_index(mat_idx), emit_count(0) {
			dir *= (1.0f / (rsd.y() = dir.len()));
		}

		float radius() const { return (rsd.x()); }
		float speed() const { return (rsd.y()); }
		float spawn_count() const { return (radius() / (2.0f * consts::PARTICLE_RADIUS)); }

		bool enabled(uint32_t max_emit_count) const { return (emit_count < max_emit_count); }
		bool delayed(float emit_dt) { return ((rsd.z() -= emit_dt) > 0.0f); }

		void update();
		void emit(t_particle& pi, float mass, float offs);

	public:
		t_vec3f pos;
		t_vec3f dir;
		t_vec3f rsd; // .x := radius, .y := speed, .z := delay-time

		uint32_t mat_index;
		uint32_t emit_count;
	};
};

#endif

