#include <GL/glew.h>
#include <GL/freeglut.h>

#include "render_state.hpp"
#include "sphsim_state.hpp"

extern ldty::t_sphsim_state g_sphsim_state;
// extern t_render_state g_render_state;

bool ldty::t_render_state::init() {
	glClearColor(0.02f, 0.01f, 0.01f, 1.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(min_fcps.x(), max_fcps.x(),  min_fcps.y(), max_fcps.y(),  min_fcps.z(), max_fcps.z());

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPointSize(2.5f * consts::PARTICLE_RADIUS * consts::WINDOW_XSIZE / consts::WORLD_XSIZE);


	glGenBuffers(1, &vbo_id);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_id);

	// glBufferData(GL_ARRAY_BUFFER, g_sphsim_state.particles.size() * vbo_elem_size, nullptr, GL_DYNAMIC_DRAW);
	glBufferStorage(GL_ARRAY_BUFFER, g_sphsim_state.particles.size() * vbo_elem_size, nullptr, vbo_flag_bits);
	// NOTE:
	//   "it is an error to specify GL_MAP_PERSISTENT_BIT if the buffer's data store
	//   was not allocated through a call to the glBufferStorage command in which the
	//   GL_MAP_PERSISTENT_BIT was also set"
	vbo_ptr = reinterpret_cast<float*>(glMapBufferRange(GL_ARRAY_BUFFER, 0, g_sphsim_state.particles.size() * vbo_elem_size, vbo_flag_bits));

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	return ((vbo_ptr != nullptr) && (glGetError() == GL_NO_ERROR));
}

void ldty::t_render_state::kill() {
	glBindBuffer(GL_ARRAY_BUFFER, vbo_id);
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDeleteBuffers(1, &vbo_id);
	assert(glGetError() == GL_NO_ERROR);
}


void ldty::t_render_state::swap_buffers() const {
	glutSwapBuffers();
	glutPostRedisplay();
}


void ldty::t_render_state::update() {
	const auto& particles = g_sphsim_state.particles;
	const auto& materials = g_sphsim_state.materials;

	for (uint32_t i = 0, j = 0; i < g_sphsim_state.num_particles; ++i) {
		const sph::t_particle& pi = particles[i];
		const sph::t_material& mat = materials[pi.mat_index];
		const sph::t_color col = mat.color * (mat.bias() + mat.scale() * pi.rest_press * pi.rest_press);

		vbo_ptr[j++] = pi.cur_pos.x();
		vbo_ptr[j++] = pi.cur_pos.y();
		vbo_ptr[j++] = pi.cur_pos.z();
		vbo_ptr[j++] = col.r;
		vbo_ptr[j++] = col.g;
		vbo_ptr[j++] = col.b;
		vbo_ptr[j++] = mat.color.a;

		assert(j <= (g_sphsim_state.num_particles * vbo_elem_size));
	}
}

void ldty::t_render_state::render() {
	if (g_sphsim_state.updated())
		update();

	submit(g_sphsim_state.num_particles);
	swap_buffers();
}

void ldty::t_render_state::submit(uint32_t num_particles) const {
	glClear(GL_COLOR_BUFFER_BIT);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_id);

	glVertexPointer(3, GL_FLOAT, vbo_elem_size, vbo_vertex_offset);
	glColorPointer(4, GL_FLOAT, vbo_elem_size, vbo_colour_offset);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDrawArrays(GL_POINTS, 0, num_particles);

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
}

