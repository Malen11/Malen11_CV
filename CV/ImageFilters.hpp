#pragma once
#include "stdafx.h"
#include "Image.hpp"

namespace CV_labs {

	//Core (for filter's functions)
	struct Core {

		int xk = 0;
		int yk = 0;
		double* data = NULL;
	};

	//Separable core (for filter's functions)
	struct SeparableCore {

		Core X;
		Core Y;
	};

	//Gradient (for filter's functions)
	struct Gradient {

		double value;
		double phi;
	};

	//class, that contains filters for images
	class ImageFilters {
	public:
		static const int kInterpolateZero = 201;				//Interpolate with 0
		static const int kInterpolateBorder = 202;				//Interpolate with border value
		static const int kInterpolateReflection = 203;			//Interpolate with border reflection value
		static const int kInterpolateWrap = 204;				//Interpolate by wrap border

		static const int kPartDerivativeDirectionXY = 301;		//Partial Derivative direction by X and Y
		static const int kPartDerivativeDirectionX = 302;		//Partial Derivative direction by X
		static const int kPartDerivativeDirectionY = 303;		//Partial Derivative direction by Y

		static const int kPartDerivativeTypeSobelCore = 306;	//Sobel core for Part Derivative 
		static const int kPartDerivativeTypePrewittCore = 307;	//Prewitt core for Part Derivative 
		static const int kPartDerivativeTypeScharrCore = 308;	//Scharr core for Part Derivative 

		static const double PI() { return std::atan(1.0) * 4; }

		//Apply filter from core to image
		static Image ApplyFilter(const Image& image, Core core, int interpolateType = kInterpolateZero, bool normalizeData = true);

		//Apply filter from core to image
		static Image ApplyFilter(const Image& image, SeparableCore core, int interpolateType = kInterpolateZero, bool normalizeData = true);

		//Apply filter from core to image data
		static double* ApplyFilterRaw(int rows, int cols, double* data, Core core, int interpolateType = kInterpolateZero);

		//Apply filter from core to image data
		static double* ApplyFilterRaw(int rows, int cols, double* data, SeparableCore core, int interpolateType = kInterpolateZero);

		//Calculate gradient value to image
		static Image CalculateGradientValue(const Image& image, int partDerivativeType = kPartDerivativeTypeSobelCore, int interpolateType = kInterpolateZero);

		//Calculate gradient value to image data
		static double* CalculateGradientValueRaw(int rows, int cols, double* data, int partDerivativeType = kPartDerivativeTypeSobelCore,int interpolateType = kInterpolateZero);

		//Calculate gradient value and angle to image data
		static Gradient* CalculateGradientRaw(int rows, int cols, double* data, int partDerivativeType = kPartDerivativeTypeSobelCore, int interpolateType = kInterpolateZero);

		static double* NormalizeData(int size, double* data);

		//Create Gauss core
		static Core GenerateGaussCore(double sigma);

		//Create Gauss core
		static Core GenerateGaussCore(double sigma, int k);

		//Create Gauss core
		static SeparableCore GenerateGaussSeparableCore(double sigma);

		//Create Gauss core
		static SeparableCore GenerateGaussSeparableCore(double sigma, int k);

		//Create Sobel core for gradient
		static Core GenerateSobelCore(int PartDerivativeDirection);

		//Create Sobel core for gradient
		static SeparableCore GenerateSobelSeparableCore(int PartDerivativeDirection);

		//Create Prewitt core for gradient
		static Core GeneratePrewittCore(int PartDerivativeDirection);

		//Create Prewitt core for gradient
		static SeparableCore GeneratePrewittSeparableCore(int PartDerivativeDirection);

		//Create Scharr core for gradient
		static Core GenerateScharrCore(int PartDerivativeDirection);

		//Create Scharr core for gradient
		static SeparableCore GenerateScharrSeparableCore(int PartDerivativeDirection);

		//Calculate virtual pixel value for interpolation 
		static double GetVirtualPixel(int row, int col, int rows, int cols, const double* const data, int interpolateType = kInterpolateZero);

		//Calculate virtual pixel value for interpolation 
		static double GetVirtualPixelZero(int row, int col, int rows, int cols, const double* const data);

		//Calculate virtual pixel value for interpolation 
		static double GetVirtualPixelBorder(int row, int col, int rows, int cols, const double* const data);
		
		//Calculate virtual pixel value for interpolation 
		static double GetVirtualPixelReflection(int row, int col, int rows, int cols, const double* const data);

		//Calculate virtual pixel value for interpolation 
		static double GetVirtualPixelWrap(int row, int col, int rows, int cols, const double* const data);

		//Calculate virtual pixel gradient value for interpolation 
		static Gradient GetVirtualPixelGradient(int row, int col, int rows, int cols, const Gradient* const data, int interpolateType = kInterpolateZero);

		//Calculate virtual pixel gradient value for interpolation 
		static Gradient GetVirtualPixelGradientZero(int row, int col, int rows, int cols, const Gradient* const data);

		//Calculate virtual pixel gradient value for interpolation 
		static Gradient GetVirtualPixelGradientBorder(int row, int col, int rows, int cols, const Gradient* const data);

		//Calculate virtual pixel gradient value for interpolation 
		static Gradient GetVirtualPixelGradientReflection(int row, int col, int rows, int cols, const Gradient* const data);

		//Calculate virtual pixel gradient value for interpolation 
		static Gradient GetVirtualPixelGradientWrap(int row, int col, int rows, int cols, const Gradient* const data);
	};

	inline double ImageFilters::GetVirtualPixel(int row, int col, int rows, int cols, const double* const data, int interpolateType) {

		if (interpolateType == kInterpolateZero) {
			return GetVirtualPixelZero(row, col, rows, cols, data);
		}
		else if (interpolateType == kInterpolateBorder) {
			return GetVirtualPixelBorder(row, col, rows, cols, data);
		}
		else if (interpolateType == kInterpolateReflection) {
			return GetVirtualPixelReflection(row, col, rows, cols, data);
		}
		else if (interpolateType == kInterpolateWrap) {
			return GetVirtualPixelWrap(row, col, rows, cols, data);
		}
		else {
			throw new std::invalid_argument("Unexpected partDerivativeType value");
		}
	}

	inline double ImageFilters::GetVirtualPixelZero(int row, int col, int rows, int cols, const double* const data) {

		if (row >= 0 && row < rows - 1 && col >= 0 && col < cols - 1) {

			return data[row * cols + col];
		}

		return 0;
	}

	inline double ImageFilters::GetVirtualPixelBorder(int row, int col, int rows, int cols, const double* const data) {

		if (row >= 0 && row < rows - 1 && col >= 0 && col < cols - 1) {

			return data[row * cols + col];
		}

		return data[std::min(rows - 1, std::max(row, 0)) *cols + std::min(cols - 1, std::max(col, 0))];
	}

	inline double ImageFilters::GetVirtualPixelReflection(int row, int col, int rows, int cols, const double* const data) {

		if (row >= 0 && row < rows - 1 && col >= 0 && col < cols - 1)
			return data[row * cols + col];

		int trow = row;
		int tcol = col;

		if (trow < 0)
			trow = 0 - row;
		else if (trow >= rows)
			trow = rows - 1 - trow % (rows - 1);

		if (tcol < 0)
			tcol = 0 - col;
		else if (tcol >= cols)
			tcol = cols - 1 - tcol % (cols - 1);

		return data[trow * cols + tcol];
	}

	inline double ImageFilters::GetVirtualPixelWrap(int row, int col, int rows, int cols, const double* const data) {

		if (row >= 0 && row < rows - 1 && col >= 0 && col < cols - 1)
			return data[row * cols + col];

		return data[((rows + row) % rows)*cols + (cols + col) % cols];
	}

	inline Gradient ImageFilters::GetVirtualPixelGradient(int row, int col, int rows, int cols, const Gradient* const data, int interpolateType) {

		if (interpolateType == kInterpolateZero) {
			return GetVirtualPixelGradientZero(row, col, rows, cols, data);
		}
		else if (interpolateType == kInterpolateBorder) {
			return GetVirtualPixelGradientBorder(row, col, rows, cols, data);
		}
		else if (interpolateType == kInterpolateReflection) {
			return GetVirtualPixelGradientReflection(row, col, rows, cols, data);
		}
		else if (interpolateType == kInterpolateWrap) {
			return GetVirtualPixelGradientWrap(row, col, rows, cols, data);
		}
		else {
			throw new std::invalid_argument("Unexpected partDerivativeType value");
		}
	}

	inline Gradient ImageFilters::GetVirtualPixelGradientZero(int row, int col, int rows, int cols, const Gradient* const data) {

		if (row >= 0 && row < rows - 1 && col >= 0 && col < cols - 1)
			return data[row * cols + col];

		return { 0, 0 };
	}

	inline Gradient ImageFilters::GetVirtualPixelGradientBorder(int row, int col, int rows, int cols, const Gradient* const data) {

		if (row >= 0 && row < rows - 1 && col >= 0 && col < cols - 1)
			return data[row * cols + col];

		return data[std::min(rows - 1, std::max(row, 0)) *cols + std::min(cols - 1, std::max(col, 0))];
	}

	inline Gradient ImageFilters::GetVirtualPixelGradientReflection(int row, int col, int rows, int cols, const Gradient* const data) {

		if (row >= 0 && row < rows - 1 && col >= 0 && col < cols - 1)
			return data[row * cols + col];

		int trow = row;
		int tcol = col;

		if (trow < 0)
			trow = 0 - row;
		else if (trow >= rows)
			trow = rows - 1 - trow % (rows - 1);

		if (tcol < 0)
			tcol = 0 - col;
		else if (tcol >= cols)
			tcol = cols - 1 - tcol % (cols - 1);

		return data[trow * cols + tcol];
	}

	inline Gradient ImageFilters::GetVirtualPixelGradientWrap(int row, int col, int rows, int cols, const Gradient* const data) {

		if (row >= 0 && row < rows - 1 && col >= 0 && col < cols - 1)
			return data[row * cols + col];

		return data[((rows + row) % rows)*cols + (cols + col) % cols];
	}

	inline Core ImageFilters::GenerateGaussCore(double sigma) {

		return GenerateGaussCore(sigma, 3 * sigma);
	}
	
	inline SeparableCore ImageFilters::GenerateGaussSeparableCore(double sigma) {

		return GenerateGaussSeparableCore(sigma, 3 * sigma);
	}
}