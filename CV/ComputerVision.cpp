#include "stdafx.h"
#include "ComputerVision.hpp"

using namespace std;

Image ComputerVision::Sobel(Image & img, int PartDerivativeType, int interpolateType) {

	double* data = img.GetNormalizeDataF();
	double* temp = SobelRaw(img.GetRowsNumber(), img.GetColsNumber(), data, PartDerivativeType, interpolateType);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);

	delete[] temp;
	delete[] data;

	return result;
}

Image ComputerVision::ApplyFilter(Image & img, Core core, int type) {

	uchar* data = img.GetData();

	double* temp = ApplyFilterRaw(img.GetRowsNumber(), img.GetColsNumber(), data, core, type);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp);

	delete[] temp;
	delete[] data;

	return result;
}

Image ComputerVision::Canny(Image & img, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType) {

	double* data = img.GetNormalizeDataF();
	double* temp = CannyRaw(img.GetRowsNumber(), img.GetColsNumber(), data, sigma, k, lowThreshhold, highThreshhold, interpolateType);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);

	delete[] temp;
	delete[] data;

	return result;
}

Image ComputerVision::GaussDefault(Image & img, double sigma, int interpolateType) {
	return Gauss(img, sigma, 3*sigma, interpolateType);
}

Image ComputerVision::Gauss(Image & img, double sigma, int k, int interpolateType) {
	
	double* data = img.GetNormalizeDataF();
	double* temp = GaussRaw(img.GetRowsNumber(), img.GetColsNumber(), data, sigma, k, interpolateType);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);

	delete[] temp;
	delete[] data;

	return result;
}
