#ifndef LIQUIDITY_GLOBAL_CONSTS_HDR
#define LIQUIDITY_GLOBAL_CONSTS_HDR

#include <cmath>

namespace consts {
	static constexpr uint32_t WINDOW_XSIZE = 800;
	static constexpr uint32_t WINDOW_YSIZE = 800;

	// origin is at bottom-left of projection volume
	static constexpr float WORLD_XSIZE  = 10.0f;
	static constexpr float WORLD_YSIZE  = (WINDOW_YSIZE * WORLD_XSIZE / WINDOW_XSIZE);
	static constexpr float WORLD_XORIG  = -(WORLD_XSIZE * 0.5f);
	static constexpr float WORLD_YORIG  = -(WORLD_YSIZE * 0.5f);
	static constexpr float WORLD_RADIUS = std::sqrt((WORLD_XSIZE * 0.5f * WORLD_XSIZE * 0.5f) + (WORLD_YSIZE * 0.5f * WORLD_YSIZE * 0.5f));

	static constexpr float MATH_PI = 3.1415926535f;

	static constexpr float REST_MASS_COEFF = 80.0f;
	static constexpr float REST_PRES_COEFF = 0.08f;
	static constexpr float NEAR_PRES_COEFF = 0.1f;
	static constexpr float SURF_TENS_COEFF = 0.0004f;

	static constexpr float L_VISCOSITY_COEFF = 0.5f;
	static constexpr float Q_VISCOSITY_COEFF = 1.0f;
	static constexpr float GRAVITATIONAL_ACC = 9.81f;

	static constexpr float PARTICLE_RADIUS    = 0.05f;
	static constexpr float NEIGHBOR_RADIUS    = 5.0f * PARTICLE_RADIUS;
	static constexpr float NEIGHBOR_RADIUS_SQ = NEIGHBOR_RADIUS * NEIGHBOR_RADIUS;


	static constexpr uint32_t SIM_STEP_RATE    = 60 * 4; // Hz
	static constexpr    float SIM_STEP_SIZE    = 1.0f / SIM_STEP_RATE; // dt (ms)
	static constexpr    float SIM_STEP_SIZE_SQ = SIM_STEP_SIZE * SIM_STEP_SIZE;
	static constexpr uint64_t SIM_STEP_TIME_NS = (1000.0f / SIM_STEP_RATE) * 1000 * 1000;
	static constexpr uint64_t WALL_SEC_TIME_NS = 1000 * 1000 * 1000;

	static constexpr float REST_DENSITY_NORM = 20.0f / (2.0f * MATH_PI * NEIGHBOR_RADIUS_SQ);
	static constexpr float NEAR_DENSITY_NORM = 30.0f / (2.0f * MATH_PI * NEIGHBOR_RADIUS_SQ);


	static constexpr float EPSILON    = 0.0000001f;
	static constexpr float EPSILON_SQ = EPSILON * EPSILON;

	static constexpr uint32_t   MAX_EMITTERS = 4;
	static constexpr uint32_t  MAX_PARTICLES = 1500 * MAX_EMITTERS;
	static constexpr uint32_t  MAX_MATERIALS = 4;
	static constexpr uint32_t  MAX_NEIGHBORS = 64 * 16;
	static constexpr uint32_t  EMITTER_DELAY = 5;

	static constexpr float     CELL_SIZE = NEIGHBOR_RADIUS;
	static constexpr uint32_t GRID_XSIZE = WORLD_XSIZE / CELL_SIZE;
	static constexpr uint32_t GRID_YSIZE = WORLD_YSIZE / CELL_SIZE;
	static constexpr uint32_t  NUM_CELLS = GRID_XSIZE * GRID_YSIZE;
};

#endif

