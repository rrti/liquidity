#ifndef LIQUIDITY_GLOBAL_TYPES_HDR
#define LIQUIDITY_GLOBAL_TYPES_HDR

#include "global_consts.hpp"

namespace sph {
	// for all i, t_ngb_cache[i][0].first is the number of neighbors for particle <i>
	typedef std::array<std::pair<uint32_t, float>, 1 + consts::MAX_NEIGHBORS> t_ngb_cache;
	typedef std::pair<uint32_t, uint32_t> t_cell_coor;
};

#endif

