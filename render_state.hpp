#ifndef LIQUIDITY_RENDER_STATE_HDR
#define LIQUIDITY_RENDER_STATE_HDR

#include "coor_tuple.hpp"
#include "color.hpp"
#include "global_consts.hpp"

namespace ldty {
	struct t_render_state {
	public:
		bool init();
		void kill();

		void update();
		void render();
		void submit(uint32_t num_particles) const;
		void swap_buffers() const;

	public:
		uint32_t vbo_id = 0;

		static constexpr void* vbo_vertex_offset = reinterpret_cast<void*>(              0);
		static constexpr void* vbo_colour_offset = reinterpret_cast<void*>(sizeof(t_vec3f));

		static constexpr uint32_t vbo_elem_size = sizeof(t_vec3f) + sizeof(sph::t_color);
		static constexpr uint32_t vbo_flag_bits = GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT;

		float* vbo_ptr = nullptr;

		// frustum clipping-plane offsets
		const t_vec3f min_fcps = t_vec3f({consts::WORLD_XORIG                      ,  consts::WORLD_YORIG                      , 0.0f});
		const t_vec3f max_fcps = t_vec3f({consts::WORLD_XORIG + consts::WORLD_XSIZE,  consts::WORLD_YORIG + consts::WORLD_YSIZE, 1.0f});
	};
};

#endif

