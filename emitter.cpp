#include "emitter.hpp"
#include "particle.hpp"
#include "sphsim_state.hpp"
#include "uniform_rng.hpp"

extern ldty::t_sphsim_state g_sphsim_state;
extern util::t_uniform_f32_rng g_uniform_rng;

void sph::t_emitter::update() {
	if (delayed(consts::SIM_STEP_SIZE * consts::EMITTER_DELAY))
		return;

	auto& particles = g_sphsim_state.particles;

	const uint32_t max_spawnees = spawn_count();
	const uint32_t max_particles = particles.size();

	const float mass = g_sphsim_state.materials[mat_index].mass();
	const float norm = 1.0f / max_spawnees;

	uint32_t  num_spawnees = 0;
	uint32_t& num_particles = g_sphsim_state.num_particles;

	// emit each spawnee symmetrically around pos, i.e. at offset ((n / N) - 0.5) * radius
	for (; num_spawnees <= max_spawnees && num_particles < max_particles; emit_count++) {
		emit(particles[num_particles++], mass, ((num_spawnees++) * norm) - 0.5f);
	}
}

void sph::t_emitter::emit(sph::t_particle& pi, float mass, float offs) {
	const float emit_offset = radius() * offs;
	const float emit_speed = speed();

	pi.cur_pos.x() = pos.x() - (dir.y() * emit_offset);
	pi.cur_pos.y() = pos.y() + (dir.x() * emit_offset);
	pi.cur_vel.x() = dir.x() * (emit_speed * g_uniform_rng());
	pi.cur_vel.y() = dir.y() * (emit_speed * g_uniform_rng());

	pi.raw_mass = mass;
	pi.inv_mass = 1.0f / pi.raw_mass;

	pi.mat_index = mat_index;
}

