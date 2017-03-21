#include <GL/glew.h>
#include <GL/freeglut.h>

#include "global_consts.hpp"
#include "render_state.hpp"
#include "sphsim_state.hpp"
#include "uniform_rng.hpp"
#include "system_timer.hpp"
#include "wall_clock.hpp"


ldty::t_sphsim_state g_sphsim_state;
ldty::t_render_state g_render_state;
util::t_uniform_f32_rng g_uniform_rng;

static util::t_wall_clock g_wall_clock;
static util::t_system_timer g_update_timer;
static util::t_system_timer g_render_timer;


static void update_frame() {
	while (g_wall_clock.get_render_time_ns() >= consts::SIM_STEP_TIME_NS) {
		g_update_timer.tick_time();
		g_sphsim_state.update();
		g_wall_clock.add_render_time_ns(-consts::SIM_STEP_TIME_NS);
		g_wall_clock.add_update_time_ns(g_update_timer.tock_time());
		g_wall_clock.add_system_time_ns(g_update_timer.tock_time());
	}

	g_wall_clock.add_update_call();
}

static void render_frame() {
	g_render_timer.tick_time();
	g_render_state.render();

	if (g_wall_clock.update(g_render_timer.tock_time())) {
		g_wall_clock.add_system_time_ns(-consts::WALL_SEC_TIME_NS);
		g_wall_clock.output_timings(stdout);
	}
}


void key_pressed(uint8_t key, int32_t, int32_t) {
	switch (key) {
		case 'p': { g_wall_clock.inv_render_dt_mult(                        ); } break;
		case 's': { g_wall_clock.add_render_time_ns(consts::SIM_STEP_TIME_NS); } break;

		case 'w': { g_sphsim_state.move_planes( consts::WORLD_RADIUS * 0.01f); } break;
		case 'q': { g_sphsim_state.move_planes(-consts::WORLD_RADIUS * 0.01f); } break;

		case 'r': { g_sphsim_state.init         (); } break;
		case 'e': { g_sphsim_state.add_explosion(); } break;

		default : {} break;
	}
}


static void run(int argc, char** argv) {
	g_uniform_rng.set_seed((argc > 1)? std::atoi(argv[1]): 1); // deterministic emitters
	g_uniform_rng.set_range(0.9f, 1.1f); // a + (b - a) * random01

	g_sphsim_state.init();
	g_render_state.init();

	glutMainLoop();

	g_render_state.kill();
}


int main(int argc, char** argv) {
	glutInitWindowSize(consts::WINDOW_XSIZE, consts::WINDOW_YSIZE);
	glutInit(&argc, argv);

	// TODO: GL4 geometry-shader for particle quads
	// glutInitContextVersion(3, 3);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
	glutInitContextProfile(GLUT_CORE_PROFILE);

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutCreateWindow("liquidity");
	glewInit();

	// register callbacks
	glutDisplayFunc(render_frame);
	glutIdleFunc(update_frame);
	glutKeyboardFunc(key_pressed);

	run(argc, argv);
	return 0;
}

