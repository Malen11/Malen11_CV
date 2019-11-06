#pragma once
#include "stdafx.h"
#include "Image.hpp"
#include "ImageFilters.hpp"

namespace CV_labs {

	class ImageDetectors {

		//static const int kHarrisResponseDirect = 401;		//response for Harris (Direct by calculate own numbers)
		static const int kHarrisResponseBase = 402;			//response for Harris (Base with k)
		static const int kHarrisResponseForstner = 403;		//response for Harris (Förstner and Gülch)
	public:
		ImageDetectors();
		~ImageDetectors();

		//Moravec interesting points detector. wk - k for window, d - move vector value, p - local max radius
		static std::vector<Point> Moravec(const Image &image, int wk, int d, int p, double Threshold);

		//Moravec interesting points detector. wk - k for window, d - move vector value, p - local max radius
		static std::vector<Point> MoravecRaw(int rows, int cols, double * data, int wk, int d, int p, double Threshold);

		//apply Harris to Image
		static std::vector<Point> Harris(Image& image, int wk, int localMinK, double Threshold, int pointsNeeded = -1, int PartDerivativeType = ImageFilters::kPartDerivativeTypeSobelCore);

		//apply Harris to Image data
		static std::vector<Point> HarrisRaw(int rows, int cols, double * data, int wk, int localMinK, double Threshold, int pointsNeeded = -1, int PartDerivativeType = ImageFilters::kPartDerivativeTypeSobelCore);

		//calculate response for Harris
		static double* HarrisResponse(int rows, int cols, double* A, double* B, double* C, int harrisResponseType = kHarrisResponseBase);

		//Use adaptive non-max supression
		static std::vector<Point> ANMS(std::vector<Point> points, int rows, int cols, double* responseMap, int pointsNeeded, double c = 0.9);
	};
}