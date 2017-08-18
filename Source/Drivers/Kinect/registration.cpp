/*
* This file is part of the OpenKinect Project. http://www.openkinect.org
*
* Copyright (c) 2014 individual OpenKinect contributors. See the CONTRIB file
* for details.
*
* This code is licensed to you under the terms of the Apache License, version
* 2.0, or, at your option, the terms of the GNU General Public License,
* version 2.0. See the APACHE20 and GPL2 files for the text of the licenses,
* or the following URLs:
* http://www.apache.org/licenses/LICENSE-2.0
* http://www.gnu.org/licenses/gpl-2.0.txt
*
* If you redistribute this file in source form, modified or unmodified, you
* may:
*   1) Leave this header intact and distribute it under the same terms,
*      accompanying it with the APACHE20 and GPL20 files, or
*   2) Delete the Apache 2.0 clause and accompany it with the GPL2 file, or
*   3) Delete the GPL v2 clause and accompany it with the APACHE20 file
* In all cases you must keep the copyright notice intact and include a copy
* of the CONTRIB file.
*
* Binary distributions must follow the binary distribution requirements of
* either License.
*/

/** @file Implementation of merging depth and color images. */

#define _USE_MATH_DEFINES
#include <math.h>
#include "registration.h"
#include <limits>


/*
* most information, including the table layout in command_response.h, was
* provided by @sh0 in https://github.com/OpenKinect/libfreenect2/issues/41
*/


class RegistrationImpl
{
public:
	RegistrationImpl();

	void apply(int dx, int dy, float dz, float& cx, float &cy) const;
	void apply(RGBQUAD* rgb, RGBQUAD* registered);
	void apply(UINT16* depth, UINT16* registered);
	void distort(int mx, int my, float& dx, float& dy) const;
	void depth_to_color(float mx, float my, float& rx, float& ry) const;

private:
	IrCameraParams depth;    ///< Depth camera parameters.
	ColorCameraParams color; ///< Color camera parameters.

	int distort_map[512 * 424];
	float depth_to_color_map_x[512 * 424];
	float depth_to_color_map_y[512 * 424];
	int depth_to_color_map_yi[512 * 424];

	const int filter_width_half;
	const int filter_height_half;
	const float filter_tolerance;
};

void RegistrationImpl::distort(int mx, int my, float& x, float& y) const
{
	// see http://en.wikipedia.org/wiki/Distortion_(optics) for description
	float dx = ((float)mx - depth.cx) / depth.fx;
	float dy = ((float)my - depth.cy) / depth.fy;
	float dx2 = dx * dx;
	float dy2 = dy * dy;
	float r2 = dx2 + dy2;
	float dxdy2 = 2 * dx * dy;
	float kr = 1 + ((depth.k3 * r2 + depth.k2) * r2 + depth.k1) * r2;
	x = depth.fx * (dx * kr + depth.p2 * (r2 + 2 * dx2) + depth.p1 * dxdy2) + depth.cx;
	y = depth.fy * (dy * kr + depth.p1 * (r2 + 2 * dy2) + depth.p2 * dxdy2) + depth.cy;
}

void RegistrationImpl::depth_to_color(float mx, float my, float& rx, float& ry) const
{
	mx = (mx - depth.cx) * depth.mq;
	my = (my - depth.cy) * depth.mq;

	float wx =
		(mx * mx * mx * color.mx_x3y0) + (my * my * my * color.mx_x0y3) +
		(mx * mx * my * color.mx_x2y1) + (my * my * mx * color.mx_x1y2) +
		(mx * mx * color.mx_x2y0) + (my * my * color.mx_x0y2) + (mx * my * color.mx_x1y1) +
		(mx * color.mx_x1y0) + (my * color.mx_x0y1) + (color.mx_x0y0);

	float wy =
		(mx * mx * mx * color.my_x3y0) + (my * my * my * color.my_x0y3) +
		(mx * mx * my * color.my_x2y1) + (my * my * mx * color.my_x1y2) +
		(mx * mx * color.my_x2y0) + (my * my * color.my_x0y2) + (mx * my * color.my_x1y1) +
		(mx * color.my_x1y0) + (my * color.my_x0y1) + (color.my_x0y0);

	rx = (wx / (color.fx * color.mq)) - (color.shift_m / color.shift_d);
	ry = (wy / color.mq) + color.cy;
}

void Registration::apply(int dx, int dy, float dz, float& cx, float &cy) const
{
	impl_->apply(dx, dy, dz, cx, cy);
}

void RegistrationImpl::apply(int dx, int dy, float dz, float& cx, float &cy) const
{
	const int index = dx + dy * 512;
	float rx = depth_to_color_map_x[index];
	cy = depth_to_color_map_y[index];

	//rx += (color.shift_m / dz);
	cx = rx * color.fx + color.cx;
}

void Registration::apply(RGBQUAD *rgb, RGBQUAD *registered)
{
	impl_->apply(rgb, registered);
}
void Registration::apply(UINT16 *depth, UINT16 *registered)
{
	impl_->apply(depth, registered);
}

void RegistrationImpl::apply(RGBQUAD *rgb, RGBQUAD *registered)
{
	// Check if all frames are valid and have the correct size
	int *color_depth_map = 0;
	const int size_depth = 512 * 424;
	const int size_color = 1920 * 1080;

	unsigned int *rgb_data = reinterpret_cast<unsigned int*>(rgb);
	unsigned int *registered_data = reinterpret_cast<unsigned int*>(registered);
	const int *map_dist = distort_map;
	const float *map_x = depth_to_color_map_x;
	const int *map_yi = depth_to_color_map_yi;

	const float color_cx = color.cx + 0.5f; // 0.5f added for later rounding

	const int size_filter_map = size_color + 1920 * filter_height_half * 2;
	const int offset_filter_map = 1920 * filter_height_half;
	float *filter_map = NULL;
	float *p_filter_map = NULL;
	int *depth_to_c_off = color_depth_map ? color_depth_map : new int[size_depth];
	int *map_c_off = depth_to_c_off;

	for (int i = 0; i < size_depth; ++i, ++map_dist, ++map_x, ++map_yi, ++map_c_off){

		const int index = *map_dist;

		if (index < 0){
			*map_c_off = -1;
			continue;
		}
		const float rx = (*map_x) * color.fx + color_cx;
		const int cx = (int)rx; // same as round for positive numbers
		const int cy = *map_yi;
		const int c_off = cx + cy * 1920;
		if (c_off < 0 || c_off >= size_color){
			*map_c_off = -1;
			continue;
		}
		*map_c_off = c_off;
	}

	map_c_off = depth_to_c_off;

	//for (int j = 0, i = 0; i < size_depth; ++i, ++map_c_off, ++j, registered_data++)
	map_c_off = map_c_off + 26 * 512;
	int i = 26 * 512;
	for (; i < 389*512; ++i, ++map_c_off,  registered_data++)
	{
		const int c_off = *map_c_off;
		*registered_data = c_off < 0 ? 0 : *(rgb_data + c_off);
	}
	if (!color_depth_map) delete[] depth_to_c_off;

}

void RegistrationImpl::apply(UINT16 *depth, UINT16 *undistorted)
{
	// Check if all frames are valid and have the correct size

	//float *depth_data = reinterpret_cast<float*>(depth);
	//float *undistorted_data = reinterpret_cast<float*>(undistorted);
	UINT16 *depth_data = depth;
	UINT16 *undistorted_data = undistorted;
	const int *map_dist = distort_map;
	const float *map_x = depth_to_color_map_x;
	const int *map_yi = depth_to_color_map_yi;

	const int size_depth = 512 * 424;
	const int size_color = 1920 * 1080;
	const float color_cx = color.cx + 0.5f; // 0.5f added for later rounding

	const int size_filter_map = size_color + 1920 * filter_height_half * 2;
	const int offset_filter_map = 1920 * filter_height_half;

	int *depth_to_c_off = new int[size_depth];
	int *map_c_off = depth_to_c_off;
	for (int i = 0; i < size_depth; ++i,  ++map_dist, ++map_x, ++map_yi, ++map_c_off)
	{
		const int index = *map_dist;

		if (index < 0)
		{
			*map_c_off = -1;
			continue;
		}
		if (depth_data[index] <= 0)
		{
			*map_c_off = -1;
			continue;
		}
		const float rx = (*map_x) * color.fx + color_cx;
		const int cx = (int)rx; // same as round for positive numbers
		const int cy = *map_yi;
		const int c_off = cx + cy * 1920;
		if (c_off < 0 || c_off >= size_color)
		{
			*map_c_off = -1;
			continue;
		}
		*map_c_off = c_off;

	}
	map_c_off = depth_to_c_off;
	map_c_off = map_c_off + 26 * 512;
	depth_data = depth_data + 26 * 512;
	int i = 26 * 512;
	for (; i < 389 * 512; ++i, ++map_c_off, ++undistorted_data, depth_data++)
	{
		const int c_off = *map_c_off;

		// check for allowed depth noise
		*undistorted_data = c_off < 0 ? 0 : *depth_data;
	}

	delete[] depth_to_c_off;

}

Registration::Registration() :impl_(new RegistrationImpl()) {}

Registration::~Registration()
{
	delete impl_;
}

RegistrationImpl::RegistrationImpl() :
filter_width_half(2), filter_height_half(1), filter_tolerance(0.01f)
{
	float mx, my;
	int ix, iy, index;
	float rx, ry;
	int *map_dist = distort_map;
	float *map_x = depth_to_color_map_x;
	float *map_y = depth_to_color_map_y;
	int *map_yi = depth_to_color_map_yi;

	for (int y = 0; y < 424; y++) {
		for (int x = 0; x < 512; x++) {
			// compute the dirstored coordinate for current pixel
			distort(x, y, mx, my);
			// rounding the values and check if the pixel is inside the image
			ix = (int)(mx + 0.5f);
			iy = (int)(my + 0.5f);
			if (ix < 0 || ix >= 512 || iy < 0 || iy >= 424)
				index = -1;
			else
				// computing the index from the coordianted for faster access to the data
				index = iy * 512 + ix;
			*map_dist++ = index;

			// compute the depth to color mapping entries for the current pixel
			depth_to_color(x, y, rx, ry);
			*map_x++ = rx;
			*map_y++ = ry;
			// compute the y offset to minimize later computations
			*map_yi++ = (int)(ry + 0.5f);
		}
	}
}

