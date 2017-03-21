#ifndef LIQUIDITY_SPHSIM_STATE_HDR
#define LIQUIDITY_SPHSIM_STATE_HDR

#include <array>

#include "coor_tuple.hpp"
#include "global_consts.hpp"
#include "global_types.hpp"

#include "color.hpp"
#include "emitter.hpp"
#include "material.hpp"
#include "particle.hpp"
#include "plane.hpp"

namespace ldty {
	struct t_sphsim_state {
	public:
		void init();
		void update();

		void emit_particles(uint32_t update_step);
		void update_particle_positions();
		void update_particle_cell_coors();
		void calc_particle_pressures();
		void calc_relaxed_particle_positions();
		void relax_particle_positions();
		void handle_plane_collisions();
		void move_planes(float dist);
		void add_explosion();

		bool updated() {
			if (num_sim_updates == num_vis_updates)
				return false;

			num_vis_updates = num_sim_updates;
			return true;
		}

	public:
		// for every particle, stores its <x,y> cell-indices
		std::array<sph::t_cell_coor, consts::MAX_PARTICLES> cell_coors;
		// for all k, index_grid[k][0] is the number of particles in cell <k>
		std::array<std::array<uint32_t, 1 + consts::MAX_NEIGHBORS>, consts::NUM_CELLS> index_grid;

		std::array<sph::t_particle, consts::MAX_PARTICLES> particles;
		std::array<sph::t_ngb_cache, consts::MAX_PARTICLES> ngb_caches;

		std::array<sph::t_material, consts::MAX_MATERIALS> materials = {{
			sph::t_material(t_vec3f({1.0f, 0.080f, 0.9f}), sph::t_color(0.2f, 0.90f, 0.3f, 1.0f)),
			sph::t_material(t_vec3f({1.4f, 0.075f, 1.5f}), sph::t_color(0.9f, 0.10f, 0.1f, 1.0f)),
			sph::t_material(t_vec3f({1.2f, 0.085f, 1.2f}), sph::t_color(0.8f, 0.70f, 0.2f, 1.0f)),
			sph::t_material(t_vec3f({0.8f, 0.090f, 0.6f}), sph::t_color(0.4f, 0.10f, 0.5f, 1.0f)),
		}};

		std::array<math::t_plane, consts::MAX_EMITTERS> planes;
		std::array<sph::t_emitter, consts::MAX_EMITTERS> emitters;

	public:
		// total number of emitted particles
		uint32_t num_particles = 0;

		// total number of simulation frames executed
		uint32_t num_sim_updates = 0;
		uint32_t num_vis_updates = 0;
	};
};

#endif

