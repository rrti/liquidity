#include "sphsim_state.hpp"
#include "uniform_rng.hpp"

extern util::t_uniform_f32_rng g_uniform_rng;


template<typename type> static inline type clamp(const type v, const type a, const type b) {
	return (std::max(a, std::min(v, b)));
}


void ldty::t_sphsim_state::init() {
	num_particles = 0;

	num_sim_updates = 0;
	num_vis_updates = 0;

	cell_coors.fill({0, 0});
	index_grid.fill({0});

	particles.fill(sph::t_particle{{}, {}, {}, {}, 0.0f, 0.0f,  0.0f, 0.0f,  0.0f, 0.0f,  0});
	ngb_caches.fill(sph::t_ngb_cache{std::pair<uint32_t, float>(0, 0.0f)});

	constexpr float r = consts::WORLD_RADIUS * 0.625f;

	for (size_t i = 0; i < planes.size(); i++) {
		const float angle = ((consts::MATH_PI * 2.0f) / planes.size()) * i;

		const t_vec3f p = t_vec3f({std::cos(angle) * r, std::sin(angle) * r, 0.0f});
		const t_vec3f n = -p * (1.0f / r);

		// place an emitter in front of every plane
		planes[i] = {n.x(), n.y(), n.z(), -r};
		emitters[i] = sph::t_emitter(i % consts::MAX_MATERIALS, p + n * 0.1f, n * 5.0f, t_vec3f({0.2f, 0.0f, 6.0f * i}));
	}
}


void ldty::t_sphsim_state::update() {
	emit_particles(num_sim_updates++);
	update_particle_positions();
	update_particle_cell_coors();
	calc_particle_pressures();
	calc_relaxed_particle_positions();
	relax_particle_positions();
	update_particle_cell_coors();
	handle_plane_collisions();
}


void ldty::t_sphsim_state::update_particle_cell_coors() {
	// clear grid particle-counters
	for (uint32_t n = 0; n < index_grid.size(); n++) {
		index_grid[n][0] = 0;
	}

	// add particle(-indices) into cells
	for (uint32_t i = 0; i < num_particles; ++i) {
		sph::t_particle& pi = particles[i];

		const float px = clamp(pi.cur_pos.x(), consts::WORLD_XORIG, consts::WORLD_XORIG + consts::WORLD_XSIZE);
		const float py = clamp(pi.cur_pos.y(), consts::WORLD_YORIG, consts::WORLD_YORIG + consts::WORLD_YSIZE);

		// never map a particle into a border-cell; simplifies indexing in calc_particle_pressures
		// const uint32_t cx = clamp(uint32_t((px - WORLD_XORIG) / CELL_SIZE), 1u, GRID_XSIZE - 2);
		// const uint32_t cy = clamp(uint32_t((py - WORLD_YORIG) / CELL_SIZE), 1u, GRID_YSIZE - 2);
		const uint32_t cx = clamp(uint32_t((px - consts::WORLD_XORIG) / consts::CELL_SIZE), 0u, consts::GRID_XSIZE - 1);
		const uint32_t cy = clamp(uint32_t((py - consts::WORLD_YORIG) / consts::CELL_SIZE), 0u, consts::GRID_YSIZE - 1);

		const uint32_t  j = cy * consts::GRID_XSIZE + cx;
		      uint32_t& n = index_grid[j][0];

		assert(j < index_grid.size());
		assert(n < index_grid[j].size());

		index_grid[j][ ++n ] = i;
		cell_coors[i] = {cx, cy};
	}
}


void ldty::t_sphsim_state::update_particle_positions() {
	for (uint32_t i = 0; i < num_particles; ++i) {
		sph::t_particle& pi = particles[i];

		// remember current position
		pi.prv_pos      = pi.cur_pos;
		pi.cur_vel.y() -= (consts::GRAVITATIONAL_ACC * consts::SIM_STEP_SIZE);
		pi.cur_pos     += (pi.cur_vel * consts::SIM_STEP_SIZE);
	}
}


void ldty::t_sphsim_state::calc_particle_pressures() {
	for (uint32_t i = 0; i < num_particles; ++i) {
		sph::t_particle& pi = particles[i];
		sph::t_ngb_cache& nc = ngb_caches[i];

		uint32_t& nc_count = nc[0].first;

		const uint32_t grid_col_idx = cell_coors[i].first;
		const uint32_t grid_row_idx = cell_coors[i].second;

		const uint32_t min_cell_row_idx = grid_row_idx - (1 * (grid_row_idx > (                     0)));
		const uint32_t min_cell_col_idx = grid_col_idx - (1 * (grid_col_idx > (                     0)));
		const uint32_t max_cell_row_idx = grid_row_idx + (1 * (grid_row_idx < (consts::GRID_YSIZE - 1)));
		const uint32_t max_cell_col_idx = grid_col_idx + (1 * (grid_col_idx < (consts::GRID_XSIZE - 1)));

		float rest_dens = 0.0f;
		float near_dens = 0.0f;

		nc_count = 0;

		// examine all immediately-adjacent neighbor cells
		for (uint32_t cell_row_idx = min_cell_row_idx; cell_row_idx <= max_cell_row_idx; cell_row_idx++) {
			for (uint32_t cell_col_idx = min_cell_col_idx; cell_col_idx <= max_cell_col_idx; cell_col_idx++) {
				const auto& cell = index_grid[cell_row_idx * consts::GRID_XSIZE + cell_col_idx];

				for (uint32_t k = 1; k <= cell[0]; k++) {
					const sph::t_particle& pj = particles[ cell[k] ];
					const t_vec3f dp = pj.cur_pos - pi.cur_pos;

					const float rr = dp * dp;

					if (rr < consts::EPSILON_SQ || rr > consts::NEIGHBOR_RADIUS_SQ)
						continue;

					const float r = std::sqrt(rr);
					const float a = 1.0f - (r / consts::NEIGHBOR_RADIUS);
					const float v = pj.raw_mass * a*a*a;

					rest_dens += (v *     consts::REST_DENSITY_NORM);
					near_dens += (v * a * consts::NEAR_DENSITY_NORM);

					if (nc_count >= consts::MAX_NEIGHBORS)
						continue;

					// consider all particles in this cell within the
					// interval [EPSILON, NGB_RADIUS] neighbours of <i>
					nc[++nc_count] = {cell[k], r};
				}
			}
		}

		pi.rest_dens = rest_dens;
		pi.near_dens = near_dens;

		pi.rest_press = consts::REST_PRES_COEFF * (rest_dens - pi.raw_mass * consts::REST_MASS_COEFF);
		pi.near_press = consts::NEAR_PRES_COEFF * (near_dens                                        );
	}
}


void ldty::t_sphsim_state::calc_relaxed_particle_positions() {
	for (uint32_t i = 0; i < num_particles; ++i) {
		sph::t_particle& pi = particles[i];
		sph::t_ngb_cache& nc = ngb_caches[i];

		t_vec3f& rp = pi.rlx_pos;
		t_vec3f  cp = pi.cur_pos;

		const uint32_t nc_count = nc[0].first;

		for (uint32_t j = 1; j <= nc_count; ++j) {
			const sph::t_particle& pj = particles[ nc[j].first ];

			const t_vec3f dp = pj.cur_pos - pi.cur_pos;
			const t_vec3f dv = pi.cur_vel - pj.cur_vel;

			const float nps = pi.near_press + pj.near_press;
			const float rps = pi.rest_press + pj.rest_press;

			const float r  = nc[j].second;
			const float a  = 1.0f - (r / consts::NEIGHBOR_RADIUS);
			const float aa = a * a;
			const float k1 = aa * a * consts::NEAR_DENSITY_NORM;
			const float k2 = aa *     consts::REST_DENSITY_NORM;
			const float d  = consts::SIM_STEP_SIZE_SQ * (nps * k1 + rps * k2) * 0.5f;

			const float st = (std::fabs(pi.raw_mass - pj.raw_mass) < consts::EPSILON);

			const float ir = 1.0f / r;
			const float u  = std::max(0.0f, (dv * dp) * ir);
			const float I  = 0.5f * consts::SIM_STEP_SIZE * a * (consts::L_VISCOSITY_COEFF * u + consts::Q_VISCOSITY_COEFF * u*u) * pi.inv_mass;

			// relaxation; x / (r*m) == x * (inv_m * inv_r)
			cp -= (dp * (d * pi.inv_mass * ir));
			// surface-tension
			cp += (dp * (consts::SURF_TENS_COEFF * pi.inv_mass) * pj.raw_mass * k2 * st);
			// viscosity; no contribution when u=dot(dv,dp) <= 0
			cp -= (dp * I * consts::SIM_STEP_SIZE);
        }

		rp = cp;
    }
}


void ldty::t_sphsim_state::relax_particle_positions() {
	for (uint32_t i = 0; i < num_particles; ++i) {
		sph::t_particle& pi = particles[i];

		pi.cur_pos = pi.rlx_pos;
		pi.cur_vel = (pi.rlx_pos - pi.prv_pos) * consts::SIM_STEP_RATE;
	}
}


void ldty::t_sphsim_state::handle_plane_collisions() {
	for (uint32_t i = 0; i < num_particles; ++i) {
		sph::t_particle& pi = particles[i];

		for (uint32_t j = 0; j < planes.size(); ++j) {
			planes[j].collision(pi);
		}
    }
}


void ldty::t_sphsim_state::move_planes(float dist) {
	for (size_t i = 0; i < planes.size(); i++) {
		planes[i].eq.w() += dist;
	}
}


void ldty::t_sphsim_state::add_explosion() {
	if (num_particles == 0)
		return;

	// map [0.9, 1.1] back to [0.0, 1.0]; pick a random particle
	const float rscale = 1.0f / g_uniform_rng.val_range();
	const float rindex = ((g_uniform_rng() - g_uniform_rng.min_range()) * rscale) * num_particles;
	const float radius = ((g_uniform_rng() - g_uniform_rng.min_range()) * rscale) * consts::WORLD_RADIUS * 0.5f;

	const uint32_t i = clamp(uint32_t(rindex), 0u, num_particles - 1);
	const uint32_t k = std::max(1.0f, radius / consts::CELL_SIZE);

	const sph::t_particle& pi = particles[i];

	const uint32_t col_idx = cell_coors[i].first;
	const uint32_t row_idx = cell_coors[i].second;

	const uint32_t min_row_idx = row_idx - std::min(                           row_idx, k);
	const uint32_t min_col_idx = col_idx - std::min(                           col_idx, k);
	const uint32_t max_row_idx = row_idx + std::min((consts::GRID_YSIZE - 1) - row_idx, k);
	const uint32_t max_col_idx = col_idx + std::min((consts::GRID_XSIZE - 1) - col_idx, k);

	// visit all particles within the explosion radius
	for (uint32_t row_idx = min_row_idx; row_idx <= max_row_idx; row_idx++) {
		for (uint32_t col_idx = min_col_idx; col_idx <= max_col_idx; col_idx++) {
			const auto& cell = index_grid[row_idx * consts::GRID_XSIZE + col_idx];

			for (uint32_t k = 1; k <= cell[0]; k++) {
				sph::t_particle& pj = particles[ cell[k] ];

				const t_vec3f dp = pj.cur_pos - pi.cur_pos;
				const t_vec3f dv = dp / std::max(dp.len(), consts::EPSILON);

				assert(!std::isinf(dv.x()));
				assert(!std::isnan(dv.x()));

				// should modify pressure, but that is recalculated next frame
				const float dist_sqr = dp.sq_len();
				const float dist_att = 1.0f - std::min(1.0f, (dist_sqr / (radius * radius)));

				pj.cur_vel += (dv * dist_att * pi.near_press);
			}
		}
	}
}


void ldty::t_sphsim_state::emit_particles(uint32_t update_step) {
	if (num_particles == particles.size())
		return;

	// emitted particles can interfere with one another; space them out
	if ((update_step % consts::EMITTER_DELAY) != 0)
		return;

	for (sph::t_emitter& emitter: emitters) {
		if (!emitter.enabled(particles.size() / emitters.size()))
			continue;

		emitter.update();
	}
}

