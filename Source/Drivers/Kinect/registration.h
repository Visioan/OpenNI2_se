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

/** @file registration.h Class for merging depth and color frames. */

#ifndef REGISTRATION_H_
#define REGISTRATION_H_

#include <string>
#include "BaseKinect2Stream.h"
#include "XnList.h"

class RegistrationImpl;

	/** @defgroup registration Registration and Geometry
	* Register depth to color, create point clouds. */

	/** Combine frames of depth and color camera. @ingroup registration
	* Right now this class uses a reverse engineered formula that uses factory
	* preset extrinsic parameters.  We do not have a clear understanding of these
	* particular extrinsic parameters and do not know how to calibrate them by
	* hand.
	*
	* If you want to perform registration with standard camera extrinsic matrix,
	* you probably need something else.
	*/
	struct IrCameraParams
	{
		float fx = (float)351.447; ///< Focal length x (pixel)
		float fy = (float)354.899; ///< Focal length y (pixel)
		float cx = (float)256.486694; ///< Principal point x (pixel)
		float cy = (float)207.852905; ///< Principal point y (pixel)
		//float fx = (float)365.053497; ///< Focal length x (pixel)
		//float fy = (float)365.053497; ///< Focal length y (pixel)
		//float cx = (float)256.486694; ///< Principal point x (pixel)
		//float cy = (float)207.852905; ///< Principal point y (pixel)
		float k1 = (float)0.0944473594; ///< Radial distortion coefficient, 1st-order
		float k2 = (float)-0.272574991; ///< Radial distortion coefficient, 2nd-order
		float k3 = (float)0.0929763764; ///< Radial distortion coefficient, 3rd-order
		float p1 = (float)0.000000000; ///< Tangential distortion coefficient
		float p2 = (float)0.000000000; ///< Tangential distortion coefficient

		float mq = (float)0.01;
	};

	struct ColorCameraParams
	{
		/** @name Intrinsic parameters */
		///@{
		//float fx = (float)1081.37207; ///< Focal length x (pixel)
		//float fy = (float)1081.37207; ///< Focal length y (pixel)
		//float cx = (float)959.500000; ///< Principal point x (pixel)
		//float cy = (float)539.500000; ///< Principal point y (pixel)
		float fx = (float)1081.37207; ///< Focal length x (pixel)
		float fy = (float)1081.37207; ///< Focal length y (pixel)
		float cx = (float)957.425; ///< Principal point x (pixel)
		float cy = (float)540.0; ///< Principal point y (pixel)
		///@}

		/** @name Extrinsic parameters
		* These parameters are used in [a formula](https://github.com/OpenKinect/libfreenect2/issues/41#issuecomment-72022111) to map coordinates in the
		* depth camera to the color camera.
		*
		* They cannot be used for matrix transformation.
		*/
		///@{
		float shift_d = (float)863.000000, shift_m = (float)52.0000000;

		float mx_x3y0 = (float)0.000770506682; // xxx
		float mx_x0y3 = (float)1.20775903e-005; // yyy
		float mx_x2y1 = (float)3.76318094e-005; // xxy
		float mx_x1y2 = (float)0.000614197692; // yyx
		float mx_x2y0 = (float)0.000680659898; // xx
		float mx_x0y2 = (float)3.44224718e-005; // yy
		float mx_x1y1 = (float)6.82990285e-005; // xy
		float mx_x1y0 = (float)0.640516996; // x
		float mx_x0y1 = (float)-0.00459693093; // y
		float mx_x0y0 = (float)0.145085305; // 1

		float my_x3y0 = (float)2.47515800e-006; // xxx
		float my_x0y3 = (float)0.000991923735; // yyy
		float my_x2y1 = (float)0.000699647586; // xxy
		float my_x1y2 = (float)3.94420204e-005; // yyx
		float my_x2y0 = (float)-4.04973107e-005; // xx
		float my_x0y2 = (float)0.000106151798; // yy
		float my_x1y1 = (float)0.000555969891; // xy
		float my_x1y0 = (float)0.00499383500; // x
		float my_x0y1 = (float)0.639892220; // y
		float my_x0y0 = (float)0.000404909311; // 1

		float mq = (float)0.002199;
		///@}
	};
	
	class  Registration
	{
	public:
		/**
		* @param depth_p Depth camera parameters. You can use the factory values, or use your own.
		* @param rgb_p Color camera parameters. Probably use the factory values for now.
		*/
		Registration();
		~Registration();

		/** Undistort and register a single depth point to color camera.
		* @param dx Distorted depth coordinate x (pixel)
		* @param dy Distorted depth coordinate y (pixel)
		* @param dz Depth value (millimeter)
		* @param[out] cx Undistorted color coordinate x (normalized)
		* @param[out] cy Undistorted color coordinate y (normalized)
		*/
		void apply(int dx, int dy, float dz, float& cx, float &cy) const;

		/** Map color images onto depth images
		* @param rgb Color image (1920x1080 BGRX)
		* @param depth Depth image (512x424 float)
		* @param[out] undistorted Undistorted depth image
		* @param[out] registered Color image for the depth image (512x424)
		* @param enable_filter Filter out pixels not visible to both cameras.
		* @param[out] bigdepth If not `NULL`, return mapping of depth onto colors (1920x1082 float). **1082** not 1080, with a blank top and bottom row.
		* @param[out] color_depth_map Index of mapped color pixel for each depth pixel (512x424).
		*/
		void apply( RGBQUAD* rgb, RGBQUAD* registered) ;
		void apply(UINT16* depth, UINT16* registered);
		/** Undistort depth
		* @param depth Depth image (512x424 float)
		* @param[out] undistorted Undistorted depth image
		*/
	
	private:
		RegistrationImpl *impl_;

		/* Disable copy and assignment constructors */
		Registration(const Registration&);
		Registration& operator=(const Registration&);
	};
	template<typename T>
	void resize(T* src,T* dst,int src_width,int src_height,int dst_width,int dst_height )
	{
		double inv_scale_x = (double)dst_width / src_width;
		double inv_scale_y = (double)dst_height / src_height;
		double scale_x = 1. / inv_scale_x, scale_y = 1. / inv_scale_y;

		int* x_ofs = new int[dst_width];

		for (int x = 0; x < dst_width; x++)
		{
			int sx = std::floor(x*scale_x);
			x_ofs[x] = min(sx, src_width - 1);
		}

		for (int y = 0; y <dst_height; y++)
		{
			T* D = dst + y*dst_width;
			int sy = min(std::floor(y*scale_y), src_height - 1);
			const T* S = src + sy*src_width;
			int x;
			for (x = 0; x <= dst_width - 2; x += 2)
				{
					T t0 = S[x_ofs[x]];
					T t1 = S[x_ofs[x + 1]];
					D[x] = t0;
					D[x + 1] = t1;
				}

			for (; x <dst_width; x++)
					D[x] = S[x_ofs[x]];

			//for (x = 0; x < dst_width; x++)
			//	*(T*)(D + x) = *(T*)(S + x_ofs[x]);

			//switch (pix_size)
			//{
			//case 1:
			//	for (x = 0; x <= dsize.width - 2; x += 2)
			//	{
			//		uchar t0 = S[x_ofs[x]];
			//		uchar t1 = S[x_ofs[x + 1]];
			//		D[x] = t0;
			//		D[x + 1] = t1;
			//	}

			//	for (; x < dsize.width; x++)
			//		D[x] = S[x_ofs[x]];
			//	break;
			//case 2:
			//	for (x = 0; x < dsize.width; x++)
			//		*(ushort*)(D + x * 2) = *(ushort*)(S + x_ofs[x]);
			//	break;
			//case 3:
			//	for (x = 0; x < dsize.width; x++, D += 3)
			//	{
			//		const uchar* _tS = S + x_ofs[x];
			//		D[0] = _tS[0]; D[1] = _tS[1]; D[2] = _tS[2];
			//	}
			//	break;
			//case 4:
			//	for (x = 0; x < dsize.width; x++)
			//		*(int*)(D + x * 4) = *(int*)(S + x_ofs[x]);
			//	break;
			//case 6:
			//	for (x = 0; x < dsize.width; x++, D += 6)
			//	{
			//		const ushort* _tS = (const ushort*)(S + x_ofs[x]);
			//		ushort* _tD = (ushort*)D;
			//		_tD[0] = _tS[0]; _tD[1] = _tS[1]; _tD[2] = _tS[2];
			//	}
			//	break;
			//case 8:
			//	for (x = 0; x < dsize.width; x++, D += 8)
			//	{
			//		const int* _tS = (const int*)(S + x_ofs[x]);
			//		int* _tD = (int*)D;
			//		_tD[0] = _tS[0]; _tD[1] = _tS[1];
			//	}
			//	break;
			//case 12:
			//	for (x = 0; x < dsize.width; x++, D += 12)
			//	{
			//		const int* _tS = (const int*)(S + x_ofs[x]);
			//		int* _tD = (int*)D;
			//		_tD[0] = _tS[0]; _tD[1] = _tS[1]; _tD[2] = _tS[2];
			//	}
			//	break;
			//default:
			//	for (x = 0; x < dsize.width; x++, D += pix_size)
			//	{
			//		const int* _tS = (const int*)(S + x_ofs[x]);
			//		int* _tD = (int*)D;
			//		for (int k = 0; k < pix_size4; k++)
			//			_tD[k] = _tS[k];
			//	}
			//}
		}

		delete[] x_ofs;

	}
#endif /* REGISTRATION_H_ */
