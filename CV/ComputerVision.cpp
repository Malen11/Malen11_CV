#include "stdafx.h"
#include "ComputerVision.hpp"

using namespace std;

std::vector<Dot> ComputerVision::Harris(Image & img, int wk, int localMinK, double Threshold, int ANMSNeeded, int PartDerivativeType) {

	double* data = img.GetNormalizeDataF();
	vector<Dot> result = HarrisRaw(img.GetRowsNumber(), img.GetColsNumber(), data, wk, localMinK, Threshold, ANMSNeeded, PartDerivativeType);

	delete[] data;

	return result;
}

double * ComputerVision::HarrisResponse(int rows, int cols, double * A, double * B, double * C, int harrisResponseType) {
	
	int size = rows * cols;
	double* result = new double[size];
	int k = 0.05; // for base

	switch (harrisResponseType) {
	//in work!
	case kHarrisResponseDirect:

		delete[] result;
		result = NULL;
		break;

	case kHarrisResponseBase:

		for (int i = 0; i < size; i++) {
			result[i] = A[i] * C[i] - pow(B[i], 2) - k * pow(A[i] + C[i], 2);
		}

		break;

	case kHarrisResponseForstner:

		for (int i = 0; i < size; i++) {
			if (A[i] * C[i] == 0)
				result[i] = 0;
			else
				result[i] = (A[i] * C[i] - pow(B[i], 2)) / (A[i] + C[i]);
		}

		break;
	}

	return result;
}

std::vector<Dot> ComputerVision::Moravec(Image & img, int wk, int d, int p, double Threshold) {
	
	double* data = img.GetNormalizeDataF();
	vector<Dot> result = MoravecRaw(img.GetRowsNumber(), img.GetColsNumber(), data, wk, d, p, Threshold);

	delete[] data;

	return result;
}

Image ComputerVision::Sobel(Image & img, int PartDerivativeType, int interpolateType) {

	double* data = img.GetNormalizeDataF();
	double* temp = SobelRaw(img.GetRowsNumber(), img.GetColsNumber(), data, PartDerivativeType, interpolateType);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);

	delete[] temp;
	delete[] data;

	return result;
}

std::vector<Dot> ComputerVision::ANMS(std::vector<Dot> points, int rows, int cols, double * responseMap, int num, double c) {

	int size = rows * cols;
	Dot point;

	if (points.size() <= num) {

		return points;
	}
	else {

		int* radius = new int[points.size()];
		int row1, col1, row2, col2;

		for (int i = 0; i < points.size(); i++) {

			radius[i] = size;
			row1 = points[i].y;
			col1 = points[i].x;

			for (int j = 0; j < points.size(); j++) {

				if (i != j) {
					row2 = points[j].y;
					col2 = points[j].x;

					if (responseMap[row1*cols + col1] < c*responseMap[row2*cols + col2]
						&& radius[i]>std::sqrt(std::pow((row1 - row2), 2) + std::pow((col1 - col2), 2))) {

						radius[i] = std::ceil(std::sqrt(std::pow((row1 - row2), 2) + std::pow((col1 - col2), 2)));
					}
				}
			}
		}

		int* tempR = new int[points.size()];
		copy(radius, &(radius[points.size()]), stdext::checked_array_iterator<int*>(tempR, points.size()));
		int tempVal;


		//for (int i = 0; i < points.size(); i++)
		//	cout<<tempR[i]<<endl;

		for (int i = 0; i < points.size(); i++) {
			for (int j = 0; j < points.size()-1-i; j++) {

				if (tempR[j] < tempR[j+1]) {
					tempVal = tempR[j];
					tempR[j] = tempR[j+1];
					tempR[j+1] = tempVal;
				}
			}
		}

		//cout  << "########";
		//for (int i = 0; i < points.size(); i++)
		//	cout << tempR[i] << endl;

		int radiusNeeded = tempR[num];
		if (num < points.size())
			if(radiusNeeded == tempR[(num-1) + 1])
				++radiusNeeded;

		delete[] tempR;

		vector<Dot> result;
		for (int i = 0; i < points.size(); i++) {
			if (radius[i] >= radiusNeeded)
				result.push_back(points[i]);
		}

		delete[] radius;

		return result;
	}
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

Core ComputerVision::CreateGaussCore(double sigma, int k) {

	int rows = 2 * k + 1;
	int cols = rows;

	Core gaussCore;
	gaussCore.xk = k;
	gaussCore.yk = k;
	gaussCore.data = new double[rows*cols];

	double alpha = 1.0 / (2 * PI()*pow(sigma, 2));
	int pos;

	for (int i = -k; i <= k; i++) {
		for (int j = -k; j <= k; j++) {

			pos = (k + i)*cols + (k + j);
			gaussCore.data[pos] = alpha * std::exp(-(pow(i, 2) + pow(j, 2)) / (2.0 * sigma * sigma));
		}
	}

	return gaussCore;
}

Image ComputerVision::Gauss(Image & img, double sigma, int k, int interpolateType) {
	
	double* data = img.GetNormalizeDataF();
	double* temp = GaussRaw(img.GetRowsNumber(), img.GetColsNumber(), data, sigma, k, interpolateType);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);

	delete[] temp;
	delete[] data;

	return result;
}
