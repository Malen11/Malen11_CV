#pragma once
#include "stdafx.h"
#include "Image.hpp"
#include "ImageFilters.hpp"
#include "ScaleSpace.hpp"

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

		//apply Harris to Image
		static std::vector<Point> Harris(Image& image, double sigmaD, int wk, double sigmaI, int localMinK, double Threshold, int pointsNeeded = -1, int PartDerivativeType = ImageFilters::kPartDerivativeTypeSobelCore);

		//apply Harris to Image data (doesn't smooth data)
		static std::vector<Point> HarrisRaw(int rows, int cols, double * data, double sigmaD, int wk, double sigmaI, int localMinK, double Threshold, int pointsNeeded = -1, int PartDerivativeType = ImageFilters::kPartDerivativeTypeSobelCore);

		//calculate response for Harris
		static double* HarrisResponseRaw(int rows, int cols, double * data, double sigmaD, int wk, double sigmaI, int PartDerivativeType = ImageFilters::kPartDerivativeTypeSobelCore, int harrisResponseType = kHarrisResponseBase);

		//Use adaptive non-max supression
		static std::vector<Point> ANMS(std::vector<Point> points, int rows, int cols, double* responseMap, int pointsNeeded, double c = 0.9);

		//find points of interest via DoG
		//static std::vector<ScalePoint> DifferenceOfGaussians(const Image &gaussians);

		//find points of interest via DoG
		static double* DifferenceOfGaussiansRaw(const Image &gaussian0, const Image &gaussian1);

		//find points of interest via DoG for array of Images
		static std::vector<Point> DifferenceOfGaussians(const std::vector<Image>  &gaussians);

		//find points of interest via DoG and Harris
		static std::vector<ScalePoint> FindPointsOfInterest(const ScaleSpace &imagePyramid, int wk, int localMinK, double harrisThreshold, double dogThreshold);

		//find points of interest via DoG and Harris
		static std::vector<ScalePoint> HarrisLaplass(const ScaleSpace &imagePyramid, int wk, int localMinK, double harrisThreshold, double dogThreshold);
	private:
		//calculate difference of arrays
		static double* arrayDifference(int size, const double* const data0, const double* const data1);

		//find extremum in 3d (for DoG)
		static std::vector<Point> ExtremumIn3DRaw(int rows, int cols, const double* const dataPrevious, const double* const dataTarget, const double* const dataNext, double threshold);
	};
}