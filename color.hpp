#ifndef LIQUIDITY_COLOR_HDR
#define LIQUIDITY_COLOR_HDR

namespace sph {
	struct t_color {
	public:
		t_color(float _r = 0.0f, float _g = 0.0f, float _b = 0.0f, float _a = 0.0f): r(_r), g(_g), b(_b), a(_a) {}

		t_color  operator *  (const t_color& c) const { return {r * c.r, g * c.g, b * c.b, a * c.a}; }
		t_color& operator *= (const t_color& c)       { r *= c.r; g *= c.g; b *= c.b; a *= c.a; return *this; }

		t_color  operator *  (const float s) const { return {r * s, g * s, b * s, a * s}; }
		t_color& operator *= (const float s)       { r *= s; g *= s; b *= s; a *= s; return *this; }

	public:
		float r;
		float g;
		float b;
		float a;
	};
};

#endif

