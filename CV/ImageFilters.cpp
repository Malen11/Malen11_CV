#include "stdafx.h"
#include "ImageFilters.hpp"

using namespace std;
using namespace CV_labs;

Image ImageFilters::ApplyFilter(Image & image, Core core, int interpolateType) {
	
	double* data = image.GetDataD();

	//double* normalizedData = NormalizeData(image.GetSize(), data);

	double* resultData = ApplyFilterRaw(image.GetRowsNumber(), image.GetColsNumber(), data, core, interpolateType);

	double* normalizedResultData = NormalizeData(image.GetSize(), resultData);

	Image result(image.GetRowsNumber(), image.GetColsNumber(), normalizedResultData);

	delete[] data;
	//delete[] normalizedData;
	delete[] resultData;
	delete[] normalizedResultData;

	return result;
}

Image ImageFilters::ApplyFilter(Image & image, SeparableCore core, int interpolateType) {

	double* data = image.GetDataD();

	//double* normalizedData = NormalizeData(image.GetSize(), data);

	double* resultData = ApplyFilterRaw(image.GetRowsNumber(), image.GetColsNumber(), data, core, interpolateType);

	double* normalizedResultData = NormalizeData(image.GetSize(), resultData);

	Image result(image.GetRowsNumber(), image.GetColsNumber(), normalizedResultData);

	delete[] data;
	//delete[] normalizedData;
	delete[] resultData;
	delete[] normalizedResultData;

	return result;
}

double * CV_labs::ImageFilters::ApplyFilterRaw(int rows, int cols, double * data, Core core, int interpolateType) {

	double* result = new double[rows * cols];

	int coreRowSize = 2 * core.xk + 1;
	int center = core.yk * coreRowSize + core.xk;

	int pos;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			result[pos] = 0;

			for (int u = -core.yk; u <= core.yk; u++) {
				for (int v = -core.xk; v <= core.xk; v++) {

					if (0 <= (i + u) && (i + u) < rows && 0 <= (j + v) && (j + v) < cols) {
						result[pos] += data[pos + u * cols + v] * core.data[center + u * coreRowSize + v];
					}
					else {
						result[pos] += GetVirtualPixel(i + u, j + v, rows, cols, data, interpolateType) * core.data[center + u * coreRowSize + v];
					}
				}
			}
		}
	}

	return result;
}

double * CV_labs::ImageFilters::ApplyFilterRaw(int rows, int cols, double * data, SeparableCore core, int interpolateType) {

	double* tmpData = ApplyFilterRaw(rows, cols, data, core.X, interpolateType);
	double* result = ApplyFilterRaw(rows, cols, tmpData, core.Y, interpolateType);

	delete[] tmpData;

	return result;
}

Image CV_labs::ImageFilters::CalculateGradientValue(Image& image, int partDerivativeType, int interpolateType) {

	double * data = image.GetDataD();

	double * resultData = CalculateGradientValueRaw(image.GetRowsNumber(), image.GetColsNumber(), data, partDerivativeType, interpolateType);

	double * normalizedData = NormalizeData(image.GetSize(), resultData);

	Image result(image.GetRowsNumber(), image.GetColsNumber(), normalizedData);

	delete[] data, resultData, normalizedData;

	return result;
}

double * CV_labs::ImageFilters::CalculateGradientValueRaw(int rows, int cols, double * data, int partDerivativeType, int interpolateType) {

	SeparableCore coreX, coreY;

	switch (partDerivativeType) {
	case kPartDerivativeTypeSobelCore:
		coreX = GenerateSobelSeparableCore(kPartDerivativeDirectionX);
		coreY = GenerateSobelSeparableCore(kPartDerivativeDirectionY);
		break;

	case kPartDerivativeTypeScharrCore:
		coreX = GenerateScharrSeparableCore(kPartDerivativeDirectionX);
		coreY = GenerateScharrSeparableCore(kPartDerivativeDirectionY);
		break;

	case kPartDerivativeTypePrewittCore:
		coreX = GeneratePrewittSeparableCore(kPartDerivativeDirectionX);
		coreY = GeneratePrewittSeparableCore(kPartDerivativeDirectionY);
		break;

	default:
		throw std::invalid_argument("Unexpected partDerivativeType value");
		break;
	}

	double* partDerX = ApplyFilterRaw(rows, cols, data, coreX, interpolateType);
	double* partDerY = ApplyFilterRaw(rows, cols, data, coreY, interpolateType);

	int size = rows * cols;
	double* result = new double[size];

	for (int i = 0; i < size; i++) {

		result[i] = std::sqrt(partDerX[i] * partDerX[i] + partDerY[i] * partDerY[i]);
	}

	delete[] partDerX, partDerY;

	return result;
}

double * CV_labs::ImageFilters::NormalizeData(int size, double * data) {

	double min, max, k;

	if (data == NULL || size == NULL) {

		throw std::invalid_argument("Data is empty");
	}

	min = max = data[0];

	for (int i = 0; i < size; i++) {

		if (data[i] > max) {
			max = data[i];
		}

		if (data[i] < min) {
			min = data[i];
		}
	}

	k = 1 / (max - min);

	double* result = new double[size];

	for (int i = 0; i < size; i++) {

		result[i] = (double)(data[i] - min) * k;
	}

	return result;
}

Core CV_labs::ImageFilters::GenerateGaussCore(double sigma, int k) {

	int rows = 2 * k + 1;
	int cols = rows;

	Core gaussCore;
	gaussCore.xk = k;
	gaussCore.yk = k;
	gaussCore.data = new double[rows * cols];

	double alpha = 1.0 / (2 * PI() * pow(sigma, 2));
	int pos;

	for (int i = -k; i <= k; i++) {
		for (int j = -k; j <= k; j++) {

			pos = (k + i) * cols + (k + j);
			gaussCore.data[pos] = alpha * std::exp(-(pow(i, 2) + pow(j, 2)) / (2.0 * sigma * sigma));
		}
	}

	return gaussCore;
}

SeparableCore CV_labs::ImageFilters::GenerateGaussSeparableCore(double sigma, int k) {

	SeparableCore separableGauss;

	separableGauss.X.xk = k;
	separableGauss.X.yk = 0;
	separableGauss.X.data = new double[2 * k + 1];

	separableGauss.Y.xk = 0;
	separableGauss.Y.yk = k;
	separableGauss.Y.data = new double[2 * k + 1];
	
	double alpha = 1.0 / (sqrt(2 * PI()) * sigma);
		
	for (int i = -k; i <= k; i++) {
		separableGauss.Y.data[k + i] = separableGauss.X.data[k + i] = alpha * std::exp(-(i * i) / (2.0 * sigma * sigma));
	}
	
	/*for (int i = -k; i <= k; i++) {
		separableGauss.Y.data[k + i] = alpha * std::exp(-(i * i) / (2.0 * sigma * sigma));
	}*/
	
	return separableGauss;
}

Core CV_labs::ImageFilters::GenerateSobelCore(int PartDerivativeDirection) {

	Core sobelCore;
	sobelCore.xk = 1;
	sobelCore.yk = 1;
	sobelCore.data = new double[9];

	double sobelCoreValsX[] = { 1, 0, -1, 2, 0, -2, 1, 0, -1 };
	double sobelCoreValsY[] = { 1, 2, 1, 0, 0, 0, -1, -2, -1 };

	switch (PartDerivativeDirection) {
	case kPartDerivativeDirectionX:
		for (int i = 0; i < 9; i++) {
			sobelCore.data[i] = sobelCoreValsX[i];
		}
		break;

	case kPartDerivativeDirectionY:
		for (int i = 0; i < 9; i++) {
			sobelCore.data[i] = sobelCoreValsY[i];
		}

		break;

	default:
		throw std::invalid_argument("Unexpected PartDerivativeDirection");
		break;
	}

	return sobelCore;
}

SeparableCore CV_labs::ImageFilters::GenerateSobelSeparableCore(int PartDerivativeDirection) {

	SeparableCore sobelSeparableCore;

	sobelSeparableCore.X.xk = 1;
	sobelSeparableCore.X.yk = 0;
	sobelSeparableCore.X.data = new double[3];

	sobelSeparableCore.Y.xk = 0;
	sobelSeparableCore.Y.yk = 1;
	sobelSeparableCore.Y.data = new double[3];

	double sobelCoreValsLine1[] = { 1, 0, -1 };
	double sobelCoreValsLine2[] = { 1, 2, 1 };

	switch (PartDerivativeDirection) {
	case kPartDerivativeDirectionX:
		for (int i = 0; i < 3; i++) {
			sobelSeparableCore.X.data[i] = sobelCoreValsLine1[i];
		}
		for (int i = 0; i < 3; i++) {
			sobelSeparableCore.Y.data[i] = sobelCoreValsLine2[i];
		}
		break;

	case kPartDerivativeDirectionY:
		for (int i = 0; i < 3; i++) {
			sobelSeparableCore.X.data[i] = sobelCoreValsLine2[i];
		}
		for (int i = 0; i < 3; i++) {
			sobelSeparableCore.Y.data[i] = sobelCoreValsLine1[i];
		}
		break;

	default:
		throw std::invalid_argument("Unexpected PartDerivativeDirection");
		break;
	}

	return sobelSeparableCore;
}

Core CV_labs::ImageFilters::GeneratePrewittCore(int PartDerivativeDirection) {

	Core prewittCore;
	prewittCore.xk = 1;
	prewittCore.yk = 1;
	prewittCore.data = new double[9];

	double prewittCoreValsX[] = { 1, 0, -1, 1, 0, -1, 1, 0, -1 };
	double prewittCoreValsY[] = { 1, 1, 1, 0, 0, 0, -1, -1, -1 };

	switch (PartDerivativeDirection) {
	case kPartDerivativeDirectionX:
		for (int i = 0; i < 9; i++) {
			prewittCore.data[i] = prewittCoreValsX[i];
		}
		break;

	case kPartDerivativeDirectionY:
		for (int i = 0; i < 9; i++) {
			prewittCore.data[i] = prewittCoreValsY[i];
		}
		break;

	default:
		throw std::invalid_argument("Unexpected PartDerivativeDirection");
		break;
	}

	return prewittCore;
}

SeparableCore CV_labs::ImageFilters::GeneratePrewittSeparableCore(int PartDerivativeDirection) {

	SeparableCore prewittSeparableCore;

	prewittSeparableCore.X.xk = 1;
	prewittSeparableCore.X.yk = 0;
	prewittSeparableCore.X.data = new double[3];

	prewittSeparableCore.Y.xk = 0;
	prewittSeparableCore.Y.yk = 1;
	prewittSeparableCore.Y.data = new double[3];

	double prewittCoreValsLine1[] = { 1, 0, -1 };
	double prewittCoreValsLine2[] = { 1, 1, 1 };

	switch (PartDerivativeDirection) {
	case kPartDerivativeDirectionX:
		for (int i = 0; i < 3; i++) {
			prewittSeparableCore.X.data[i] = prewittCoreValsLine1[i];
		}
		for (int i = 0; i < 3; i++) {
			prewittSeparableCore.Y.data[i] = prewittCoreValsLine2[i];
		}
		break;

	case kPartDerivativeDirectionY:
		for (int i = 0; i < 3; i++) {
			prewittSeparableCore.X.data[i] = prewittCoreValsLine2[i];
		}
		for (int i = 0; i < 3; i++) {
			prewittSeparableCore.Y.data[i] = prewittCoreValsLine1[i];
		}
		break;

	default:
		throw std::invalid_argument("Unexpected PartDerivativeDirection");
		break;
	}

	return prewittSeparableCore;
}

Core CV_labs::ImageFilters::GenerateScharrCore(int PartDerivativeDirection) {

	Core scharrCore;
	scharrCore.xk = 1;
	scharrCore.yk = 1;
	scharrCore.data = new double[9];

	double scharrCoreValsX[] = { 3, 0, -3, 10, 0, -10, 3, 0, -3 };
	double scharrCoreValsY[] = { 3, 10, 3, 0, 0, 0, -3, -10, -3 };

	switch (PartDerivativeDirection) {
	case kPartDerivativeDirectionX:
		for (int i = 0; i < 9; i++) {
			scharrCore.data[i] = scharrCoreValsX[i];
		}
		break;

	case kPartDerivativeDirectionY:
		for (int i = 0; i < 9; i++) {
			scharrCore.data[i] = scharrCoreValsY[i];
		}
		break;

	default:
		throw std::invalid_argument("Unexpected PartDerivativeDirection");
		break;
	}

	return scharrCore;
}

SeparableCore CV_labs::ImageFilters::GenerateScharrSeparableCore(int PartDerivativeDirection) {

	SeparableCore scharrSeparableCore;

	scharrSeparableCore.X.xk = 1;
	scharrSeparableCore.X.yk = 0;
	scharrSeparableCore.X.data = new double[3];

	scharrSeparableCore.Y.xk = 0;
	scharrSeparableCore.Y.yk = 1;
	scharrSeparableCore.Y.data = new double[3];

	double scharrCoreValsLine1[] = { std::sqrt(3), 0, -std::sqrt(3) };
	double scharrCoreValsLine2[] = { std::sqrt(3), 10. / std::sqrt(3), std::sqrt(3) };

	switch (PartDerivativeDirection) {
	case kPartDerivativeDirectionX:
		for (int i = 0; i < 3; i++) {
			scharrSeparableCore.X.data[i] = scharrCoreValsLine1[i];
		}
		for (int i = 0; i < 3; i++) {
			scharrSeparableCore.Y.data[i] = scharrCoreValsLine2[i];
		}
		break;

	case kPartDerivativeDirectionY:
		for (int i = 0; i < 3; i++) {
			scharrSeparableCore.X.data[i] = scharrCoreValsLine2[i];
		}
		for (int i = 0; i < 3; i++) {
			scharrSeparableCore.Y.data[i] = scharrCoreValsLine1[i];
		}
		break;

	default:
		throw std::invalid_argument("Unexpected PartDerivativeDirection");
		break;
	}

	return scharrSeparableCore;
}
