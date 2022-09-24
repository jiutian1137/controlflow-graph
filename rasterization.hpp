#pragma once

/// Geometry Rasterization.
///@license Free 
///@review 2022-9-23 
///@author LongJiangnan, Jiang1998Nan@outlook.com 
#define _MATH_GEOMETRY_RASTERIZATION_

#include <cassert>
#include <utility>

#include <algorithm>// std::clamp

///@diagram 
///      1    2    3    4    5    6    7    8    9    10   11 unit length = 1.0[raster] = each lattice width;
/// 0----+----+----+----+----+----+----+----+----+----+----+----+----->
/// |    |    |    |    |    |    |    |    |    |    |    |    |
/// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- 1
/// |    |    |    |    |    |    |    |    |    |    |    |    |
/// + -- + -- _-*- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- 2
/// |    |   *|  \ |    |    |    |    | *  |    |    |    |    |
/// + -- + -- \ --\+ -- + -- + -- + -- +/-\ + -- + -- + -- + -- + -- 3
/// |    |    |\   \    |    |    |   /|   \|    |    |    |    |
/// + -- + -- + \- +\-- + -- + -- + /- + -- \ -- + -- + -- + -- + -- 4
/// |    |    |  \ | \  |    |    /    |    |\   |    |    |    |
/// + -- + -- + --\+ -\ + -- + -* + -- + -- + \- + -- + -- + -- + -- 5
/// |    |    |    \   \|    |    |  \ |    |  \ |    |    |    |
/// + -- + -- + -- +\-- \ -- + -- + -- + -\ + --\+ -- + -- + -- + -- 6
/// |    |    |    | \  |\   |    |    |    |  \ \    |    |    |
/// + -- + -- + -- + -\ + \- + -- + -- + -- + -- +*-- + -- + -- + -- 7
/// |    |    |    |   \|  \ |    |    |    |    |    |    |    |
/// + -- + -- + -- + -- \ _-*+ -- + -- + -- + -- + -- + -- + -- + -- 8
/// |    |    |    |    |* width  |    |    |    |    |    |    |
/// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- 9
/// |    |    |    |    |    |    |    |    |    |    |    |    |
/// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- 10
/// |    |    |    |    |    |    |    |    |    |    |    |    |
///\|
/// 
///@theory 
/// What color displayed at a point ?
///		A point have 1.0[raster] width, we sample color many times in the width,
///		then average these colors will get the point color.
///		In other word, color of a point is statistical phenomenon.
/// 
///@note 
/// Somethings to remember:
///		1. single sampling use midpoint, because midpoint implies average .
/// 
///		2. two form of raster: 
///			one is "each square" is raster, that midpoint is { Point(M.5,N.5), M,N is integer }.
///			the other is "each point" is raster, that midpoint on the each point { Point(M,N), M,N is integer }. 
///			 (this is we used, because we can set boundary value and we can avoid coordinate transform)
/// 
///		3. infinitely thin line cannot display, because cannot sample infinitely thin line.
/// 
///		4. raster in 2D also named "pixel", it in 3D also named "voxel".
///			pixel and voxel requires positive[0,inf), but raster allow negative on (-inf,inf).
/// 
namespace math { namespace geometry {
	///@link "https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm"
	template<typename function_ixy>
	void rfor_each_line(int x0, int y0, int x1, int y1, const function_ixy& raster_proc) {
		int dx = x0 < x1 ? (x1 - x0) : (x0 - x1);
		int sx = x0 < x1 ? 1 : -1;
		int dy = -( y0 < y1 ? (y1 - y0) : (y0 - y1) );
		int sy = y0 < y1 ? 1 : -1;
		int err = dx + dy;  /* error value e_xy */
		while (true) {   /* loop */
			raster_proc(x0,y0);
			if (x0 == x1 && y0 == y1) {
				break;
			}
			int e2 = err * 2;
			if (e2 >= dy) { /* e_xy+e_x > 0 */
				err += dy;
				x0 += sx;
			}
			if (e2 <= dx) { /* e_xy+e_y < 0 */
				err += dx;
				y0 += sy;
			}
		}
	}

	template<typename function_ixy>
	void rfor_each_line(int x0, int y0, int x1, int y1,
		int xlowest, int xmax, int ylowest, int ymax, const function_ixy& raster_proc) {
		rfor_each_line(std::clamp(x0, xlowest, xmax), std::clamp(x1, xlowest, xmax), std::clamp(y0, ylowest, ymax), std::clamp(y1, ylowest, ymax),
			raster_proc);
	}

	///@diagram
	/// 0--------------->
	/// |    p0 *
	/// |      /  \
	/// |     /     * p2
	/// |    /   /
	/// |   / /
	/// |  * p1
	///\|
	/// 
	/// Clockwise or From Small To Large.
	///@note many problems about triangles only have two cases, which can be perfectly solved if you work slightly harder.
	/// 

	template<typename float_, typename function_ixy>
	inline void rfor_each_horizonline(float_ x, float_ xupper, float_ y, const function_ixy& raster_proc) {
		int ix      = int(ceil(x));
		int ixupper = int(floor(xupper));
		int iy      = int(y);
		for ( ; ix <= ixupper; ++ix) 
			raster_proc(ix, iy);
	}

	template<typename float_, typename function_ixy>
	void rfor_each_triangle(float_ x0, float_ y0, float_ x1, float_ y1, float_ x2, float_ y2, bool cullface, const function_ixy& raster_proc) {
		float_ orient_or_area = /* cross(p1-p0,p2-p0)[2] = */
			/* (p1-p0)[0]*(p2-p0)[1] - (p1-p0)[1]*(p2-p0)[0] = */
			(x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);

		/// ignore zero-area.
		if (orient_or_area == 0) 
			return;
		
		/// process cullface.(when we look forward, 
		if (orient_or_area > 0)//then triangle faced me requires that backward.)
			if (cullface) {
				return;
			} else {
				float_ tmp = x1; x1 = x2; x2 = tmp;
				       tmp = y1; y1 = y2; y2 = tmp;
			}

		/// rotate order, because cannot modify orient.
		if (y1 < y0 && y1 <= y2) {
			float_ tmp = x2; x2 = x0; x0 = x1; x1 = tmp;
			       tmp = y2; y2 = y0; y0 = y1; y1 = tmp;
		} else if (y2 < y0 && y2 <= y1) {
			float_ tmp = x1; x1 = x0; x0 = x2; x2 = tmp;
			       tmp = y1; y1 = y0; y0 = y2; y2 = tmp;
		}

		/// divide into binary cases is to optimization. (it can combined into one) 
		if (y1 < y2) {
			float_ y = ceil(y0);
			float_ ymid = floor(y1);
			float_ yend = floor(y2);
			assert( y <= ymid && y <= yend && ymid <= yend );

			float_ slope1 = (x1 - x0)/(y1 - y0);
			float_ slope2 = (x2 - x0)/(y2 - y0);
			///@note 
			/// positive triangle forloop[0,mid][mid+1,..] and negative triangle forloop[0,mid)[mid,..] are mutually exclusive,
			/// so we must branch these two cases, we use the else case transform negative triangle (process first line) into positive triangle forloop.
			if (!isinf(slope1)) {
				for ( ; y <= ymid; ++y) {
					///@note 
					/// cannot pre-substract ystart,
					///  because the y0 in "(y - y0)" contains fraction-part, use integer float still need to substract in each step. and
					/// cannot use fraction float "++(y - y0)" to do inc as integer,
					///  because can avoid all errors why allow errors.
					/*int ix      = int(ceil(slope1*(y - y0) + x0));
					int ixupper = int(floor(slope2*(y - y0) + x0));
					int iy      = int(y);
					for ( ; ix <= ixupper; ++ix) 
						raster_proc(ix, iy);*/
					rfor_each_horizonline(slope1*(y - y0) + x0, slope2*(y - y0) + x0, y,
						raster_proc);
				}
			} else {
				///    1 -- 0     0 -- 1
				/// (1) \      (2)      \
				///      2               2
				/// these are only two cases, the other(2) is backfaced.
				assert( x0 >= x1 );
					rfor_each_horizonline(x1, x0, y,
						raster_proc);
				++y;
			}

			if (ymid != yend) {
				slope1 = (x2 - x1)/(y2 - y1);
				for ( ; y <= yend; ++y) {
					rfor_each_horizonline(slope1*(y - y1) + x1, slope2*(y - y0) + x0, y,
						raster_proc);
				}
			}
		} else {
			float_ y = ceil(y0);
			float_ ymid = floor(y2);
			float_ yend = floor(y1);
			assert( y <= ymid && y <= yend && ymid <= yend );

			float_ slope1 = (x1 - x0)/(y1 - y0);
			float_ slope2 = (x2 - x0)/(y2 - y0);
			if (!isinf(slope2)) {
				for ( ; y <= ymid; ++y) {
					rfor_each_horizonline(slope1*(y - y0) + x0, slope2*(y - y0) + x0, y,
						raster_proc);
				}
			} else {
				///    2 -- 0     0 -- 2
				/// (1) \      (2)      \
				///      1               1
				/// these are only two cases, the other(1) is backfaced.
				assert( x2 >= x0 );
					rfor_each_horizonline(x0, x2, y,
						raster_proc);
				++y;
			}

			if (ymid != yend) {
				slope2 = (x1 - x2)/(y1 - y2);
				for ( ; y <= yend; ++y) {
					rfor_each_horizonline(slope1*(y - y0) + x0, slope2*(y - y2) + x2, y,
						raster_proc);
				}
			}
		}
	};

	template<typename float_, typename function_ixy>
	inline void rfor_each_horizonline(float_ x, float_ xupper, float_ y, float_ xlowest, float_ xmax, const function_ixy& raster_proc) {
		using std::min, std::max;
		if (x > xmax || xupper < xlowest) 
			return;
		int ix      = int(ceil(max(x,xlowest)));
		int ixupper = int(floor(min(xupper,xmax)));
		int iy      = int(y);
		for ( ; ix <= ixupper; ++ix) 
			raster_proc(ix, iy);
	}

	template<typename float_, typename function_ixy>
	void rfor_each_triangle(float_ x0, float_ y0, float_ x1, float_ y1, float_ x2, float_ y2, 
		float_ xlowest, float_ xmax, float_ ylowest, float_ ymax, bool cullface, const function_ixy& raster_proc) {
		using std::min, std::max;

		float_ orient_or_area = /* cross(p1-p0,p2-p0)[2] = */
			/* (p1-p0)[0]*(p2-p0)[1] - (p1-p0)[1]*(p2-p0)[0] = */
			(x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);

		/// ignore zero-area.
		if (orient_or_area == 0) 
			return;
		
		/// process cullface.(when we look forward, 
		if (orient_or_area > 0)//then triangle faced me requires that backward.)
			if (cullface) {
				return;
			} else {
				float_ tmp = x1; x1 = x2; x2 = tmp;
				       tmp = y1; y1 = y2; y2 = tmp;
			}

		/// rotate order, because cannot modify orient.
		if (y1 < y0 && y1 <= y2) {
			float_ tmp = x2; x2 = x0; x0 = x1; x1 = tmp;
			       tmp = y2; y2 = y0; y0 = y1; y1 = tmp;
		} else if (y2 < y0 && y2 <= y1) {
			float_ tmp = x1; x1 = x0; x0 = x2; x2 = tmp;
			       tmp = y1; y1 = y0; y0 = y2; y2 = tmp;
		}

		/// divide into binary cases is to optimization. (it can combined into one) 
		if (y1 < y2) {
			if (y0 > ymax || y2 < ylowest) 
				return;
			float_ y = ceil(max(y0,ylowest));
			float_ ymid = floor(min(y1,ymax));
			float_ yend = floor(min(y2,ymax));

			float_ slope1 = (x1 - x0)/(y1 - y0);
			float_ slope2 = (x2 - x0)/(y2 - y0);
			if (!isinf(slope1)) {
				for ( ; y <= ymid; ++y) {
					rfor_each_horizonline(slope1*(y - y0) + x0, slope2*(y - y0) + x0, y,
						xlowest, xmax, raster_proc);
				}
			} else {
				///    1 -- 0     0 -- 1
				/// (1) \      (2)      \
				///      2               2
				/// these are only two cases, the other(2) is backfaced.
				assert( x0 >= x1 );
					rfor_each_horizonline(x1, x0, y,
						xlowest, xmax, raster_proc);
				++y;
			}

			if (ymid != yend) {
				slope1 = (x2 - x1)/(y2 - y1);
				for ( ; y <= yend; ++y) {
					rfor_each_horizonline(slope1*(y - y1) + x1, slope2*(y - y0) + x0, y,
						xlowest, xmax, raster_proc);
				}
			}
		} else {
			if (y0 > ymax || y1 < ylowest)
				return;
			float_ y = ceil(max(y0,ylowest));
			float_ ymid = floor(min(y2,ymax));
			float_ yend = floor(min(y1,ymax));

			float_ slope1 = (x1 - x0)/(y1 - y0);
			float_ slope2 = (x2 - x0)/(y2 - y0);
			if (!isinf(slope2)) {
				for ( ; y <= ymid; ++y) {
					rfor_each_horizonline(slope1*(y - y0) + x0, slope2*(y - y0) + x0, y,
						xlowest, xmax, raster_proc);
				}
			} else {
				///    2 -- 0     0 -- 2
				/// (1) \      (2)      \
				///      1               1
				/// these are only two cases, the other(1) is backfaced.
				assert( x2 >= x0 );
					rfor_each_horizonline(x0, x2, y,
						xlowest, xmax, raster_proc);
				++y;
			}

			if (ymid != yend) {
				slope2 = (x1 - x2)/(y1 - y2);
				for ( ; y <= yend; ++y) {
					rfor_each_horizonline(slope1*(y - y0) + x0, slope2*(y - y2) + x2, y,
						xlowest, xmax, raster_proc);
				}
			}
		}
	}

#if 0
int main() {
	char bitmap[32][32];

	for (size_t i = 0; i != 32; ++i)
		for (size_t j = 0; j != 32; ++j)
			bitmap[i][j] = '0';

	/// 0.0  1.0  2.0  3.0  4.0
	/// |    |    |    |    |    |    |
	/// + -- 1 -- X -- X -- 0 -- + -- + -- 5.0
	/// |    |    |    |    |    |    |
	/// + -- + -- + -- + -- + -- 2 -- + -- 6.0
	/// |    |    |    |    |    |    |
	math::geometry::rfor_each_triangle<float>(4,5, 1,5, 5,6, true, [&](int x,int y){ bitmap[y][x] = '1'; }); //infinity slope test1
	// print...

	for (size_t i = 0; i != 32; ++i)
		for (size_t j = 0; j != 32; ++j)
			bitmap[i][j] = '0';

	/// 0.0  1.0  2.0  3.0  4.0  5.0  6.0  7.0
	/// |    |    |    |    |    |    |    |    |    |
	/// + -- + -- + -- + -- 0 -- X -- X -- X -- 2 -- + -- 5.0
	/// |    |    |    |    |    |    |    |    |    |
	/// + -- + -- + -- + -- + -- 1 -- + -- + -- + -- + -- 6.0
	/// |    |    |    |    |    |    |    |    |    |
	math::geometry::rfor_each_triangle<float>(4,5, 5,6, 8,5, true, [&](int x, int y) { bitmap[y][x] = '1'; }); //infinity slope test2
	// print...

	math::geometry::rfor_each_triangle<float>(4*2 - 10,5*2-15, 7*2 - 10,7*2, 8*2 - 10,3*2, 
		0, 63, 0, 63, true, [&](int x,int y){ bitmap[y][x] = '1'; });
	// print...
}
#endif
} }//end of namespace math::geometry

#if 0
	auto rastertraignel = [&](int x0, int y0, int x1, int y1, int x2, int y2) {
		/// desiable order. 
		/// 0--------------->
		/// |    p0 *
		/// |      /  \
		/// |     /     * p2
		/// |    /   /
		/// |   / /
		/// |  * p1
		/// |/
		if (y1 < y0 && y1 < y2) {
			std::swap(x0, x1);
			std::swap(y0, y1);
		} else if (y2 < y0 && y2 < y1) {
			std::swap(x0, x2);
			std::swap(y0, y2);
		}
		if (x2 < x1) {
			std::swap(x1, x2);
			std::swap(y1, y2);
		}

		int dx0 = x0 < x1 ? (x1 - x0) : (x0 - x1);
		int sx0 = x0 < x1 ? 1 : -1;
		assert( y0 < y1 );
		int dy0 = -(y1 - y0);
		//int sy0 = 1;
		int err0 = dx0 + dy0;  /* error value e_xy */

		int dx1 = x0 < x2 ? (x2 - x0) : (x0 - x2);
		int sx1 = x0 < x2 ? 1 : -1;
		assert( y0 < y2 );
		int dy1 = -(y2 - y0);
		//int sy1 = 1;
		int err1 = dx1 + dy1;  /* error value e_xy */

		int xi = x0, xiend = x0;
		int y = y0;

		int ymid, yend;
		if (y1 < y2) {
			ymid = y1;
			yend = y2;
		} else {
			ymid = y2;
			yend = y1;
		}

		while (true) {
			bool b = false;
			while (true) {
				if (xiend == x2/* && y == y2*/) {
					break;
				}
				int e2 = err1 * 2;
				if (e2 >= dy1) { /* e_xy+e_x > 0 */
					if (e2 <= dx1) { /* e_xy+e_y < 0 */
						b = true;
						break;
					}
					err1 += dy1;
					xiend += sx1;
				}
				if (e2 <= dx1) { /* e_xy+e_y < 0 */
					/*err1 += dx1;
					y1 += sy1;*/
					break;
				}
			}

			for (int i = xi; i <= xiend; ++i) {
				bitmap[y][i] = '1';
			}

			if (y == ymid)
				break;

			while (true) {
				assert( xi != x1/* && y != y1*/ );
				int e2 = err0 * 2;
				if (e2 >= dy0) { /* e_xy+e_x > 0 */
					err0 += dy0;
					xi += sx0;
				}
				if (e2 <= dx0) { /* e_xy+e_y < 0 */
					/*err0 += dx0;
					y0 += sy0;*/
					break;
				}
			}

			err0 += dx0;
			err1 += dx1;
			if (b) {
				err1 += dy1;
				xiend += sx1;
			}
			y += 1;
		}

#if 1
		if (y1 != y2) {
			if (y1 < y2) {
				dx0 = dx1;
				sx0 = sx1;
				dy0 = dy1;
				err0 = err1;
			}
			assert( x1 < x2 );
			dx1 = y1 < y2 ? (x2 - x1) : -(x2 - x1);
			sx1 = y1 < y2 ? 1 : -1;
			dy1 = -(y1 < y2 ? (y2 - y1) : (y1 - y2));
			err1 = dx1 + dy1;  /* error value e_xy */

			while (true) {
				//bool b = false;
				//while (true) {
				//	if (xiend == (sx1 == -1 ? x1 : x2)/* && y == y2*/) {
				//		break;
				//	}
				//	int e2 = err1 * 2;
				//	if (e2 >= dy1) { /* e_xy+e_x > 0 */
				//		if (e2 <= dx1) { /* e_xy+e_y < 0 */
				//			b = true;
				//			break;
				//		}
				//		err1 += dy1;
				//		xiend += sx1;
				//	}
				//	if (e2 <= dx1) { /* e_xy+e_y < 0 */
				//		/*err1 += dx1;
				//		y1 += sy1;*/
				//		break;
				//	}
				//}

				for (int i = xi; i <= xiend; ++i) {
					bitmap[y][i] = '1';
				}

				if (y == yend)
					break;

				while (true) {
					assert( xi != x1/* && y != y1*/ );
					int e2 = err0 * 2;
					if (e2 >= dy0) { /* e_xy+e_x > 0 */
						err0 += dy0;
						xi += sx0;
					}
					if (e2 <= dx0) { /* e_xy+e_y < 0 */
						/*err0 += dx0;
						y0 += sy0;*/
						break;
					}
				}

				while (true) {
					//if (xiend == (sx1 == -1 ? x1 : x2)/* && y == y2*/) {
					//	break;
					//}
					int e2 = err1 * 2;
					if (e2 >= dy1) { /* e_xy+e_x > 0 */
						err1 += dy1;
						xiend += sx1;
					}
					if (e2 <= dx1) { /* e_xy+e_y < 0 */
						/*err1 += dx1;
						y1 += sy1;*/
						break;
					}
				}

				err0 += dx0;
				err1 += dx1;
				/*if (b) {
					err1 += dy1;
					xiend += sx1;
				}*/
				y += 1;
			}
		}
#endif
	};
#endif