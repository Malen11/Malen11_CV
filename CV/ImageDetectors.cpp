#include "stdafx.h"
#include "ImageDetectors.hpp"

using namespace std;
using namespace CV_labs;

ImageDetectors::ImageDetectors()
{
}


ImageDetectors::~ImageDetectors()
{
}

//Moravec interesting points detector. wk - k for window, d - move vector value, p - local max radius
std::vector<Point> CV_labs::ImageDetectors::Moravec(const Image & image, int wk, int d, int p, double Threshold) {

	double* data = image.GetNormalizedImageDataD();
	vector<Point> result = MoravecRaw(image.GetRowsNumber(), image.GetColsNumber(), data, wk, d, p, Threshold);
	
	delete[] data;
	
	return result;
}

//Moravec interesting points detector. wk - k for window, d - move vector value, p - local max radius
std::vector<Point> ImageDetectors::MoravecRaw(int rows, int cols, double * data, int wk, int d, int p, double Threshold) {

	int size = rows * cols;
	int pos;

	double** shifts = new double*[8];

	//смещения, начиная с левого и против часовой
	for (int i = 0; i < 8; i++) {

		shifts[i] = new double[size];
	}

	//сначала считаем все смещения для каждого пикселя (их всего 8)
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			shifts[0][pos] = pow(data[pos] - ImageFilters::GetVirtualPixelReflection(i, j + d, rows, cols, data), 2);
			shifts[1][pos] = pow(data[pos] - ImageFilters::GetVirtualPixelReflection(i - d, j + d, rows, cols, data), 2);
			shifts[2][pos] = pow(data[pos] - ImageFilters::GetVirtualPixelReflection(i - d, j, rows, cols, data), 2);
			shifts[3][pos] = pow(data[pos] - ImageFilters::GetVirtualPixelReflection(i - d, j - d, rows, cols, data), 2);
			shifts[4][pos] = pow(data[pos] - ImageFilters::GetVirtualPixelReflection(i, j - d, rows, cols, data), 2);
			shifts[5][pos] = pow(data[pos] - ImageFilters::GetVirtualPixelReflection(i + d, j - d, rows, cols, data), 2);
			shifts[6][pos] = pow(data[pos] - ImageFilters::GetVirtualPixelReflection(i + d, j, rows, cols, data), 2);
			shifts[7][pos] = pow(data[pos] - ImageFilters::GetVirtualPixelReflection(i + d, j + d, rows, cols, data), 2);
		}
	}

	double temp;
	double* s = new double[size];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			for (int shift = 0; shift < 8; shift++) {

				temp = 0;

				for (int u = -wk; u <= wk; u++) {
					for (int v = -wk; v <= wk; v++) {

						if (u != 0 && v != 0)
							temp += ImageFilters::GetVirtualPixelBorder(i + u, j + v, rows, cols, shifts[shift]);
					}
				}

				if (shift == 0)
					s[pos] = temp;
				else
					s[pos] = std::min(s[pos], temp);
			}
		}
	}

	std::vector<Point> result;
	Point point;
	bool isPoint;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			isPoint = true;
			pos = i * cols + j;

			if (s[pos] > Threshold) {

				for (int u = -p; u <= p; u++) {
					for (int v = -p; v <= p; v++) {

						if (!(u == 0 && v == 0) && s[pos] <= ImageFilters::GetVirtualPixelZero(i + u, j + v, rows, cols, s))
							isPoint = false;
					}
				}
			}
			else {
				isPoint = false;
			}

			if (isPoint) {

				point.x = j;
				point.y = i;

				result.push_back(point);
			}
		}
	}

	for (int shift = 0; shift < 8; shift++) {
		delete[] shifts[shift];
	}

	delete[] shifts;
	delete[] s;

	return result;
}

std::vector<Point> ImageDetectors::Harris(Image & image, int wk, int localMinK, double Threshold, int pointsNeeded, int PartDerivativeType) {

	double* data = image.GetNormalizedImageDataD();
	vector<Point> result = HarrisRaw(image.GetRowsNumber(), image.GetColsNumber(), data, wk, localMinK, Threshold, pointsNeeded, PartDerivativeType);

	delete[] data;

	return result;
}

std::vector<Point> ImageDetectors::HarrisRaw(int rows, int cols, double * data, int wk, int localMinK, double Threshold, int pointsNeeded, int PartDerivativeType) {

	int size = rows * cols;
	int pos;

	//part derivative X, Y
	double *partDerX = NULL, *partDerY = NULL;

	if (PartDerivativeType == ImageFilters::kPartDerivativeTypeSobelCore) {

		partDerX = ImageFilters::ApplyFilterRaw(rows, cols, data, ImageFilters::GenerateSobelSeparableCore(ImageFilters::kPartDerivativeDirectionX), ImageFilters::kInterpolateBorder);
		partDerY = ImageFilters::ApplyFilterRaw(rows, cols, data, ImageFilters::GenerateSobelSeparableCore(ImageFilters::kPartDerivativeDirectionY), ImageFilters::kInterpolateBorder);
	}
	else if (PartDerivativeType == ImageFilters::kPartDerivativeTypePrewittCore) {

		partDerX = ImageFilters::ApplyFilterRaw(rows, cols, data, ImageFilters::GeneratePrewittSeparableCore(ImageFilters::kPartDerivativeDirectionX), ImageFilters::kInterpolateBorder);
		partDerY = ImageFilters::ApplyFilterRaw(rows, cols, data, ImageFilters::GeneratePrewittSeparableCore(ImageFilters::kPartDerivativeDirectionY), ImageFilters::kInterpolateBorder);
	}
	else if (PartDerivativeType == ImageFilters::kPartDerivativeTypeScharrCore) {

		partDerX = ImageFilters::ApplyFilterRaw(rows, cols, data, ImageFilters::GenerateScharrSeparableCore(ImageFilters::kPartDerivativeDirectionX), ImageFilters::kInterpolateBorder);
		partDerY = ImageFilters::ApplyFilterRaw(rows, cols, data, ImageFilters::GenerateScharrSeparableCore(ImageFilters::kPartDerivativeDirectionY), ImageFilters::kInterpolateBorder);
	}

	double* A = new double[size];
	double* B = new double[size];
	double* C = new double[size];

	//gauss core for gauss weight
	int gaussSizeD = 2 * wk + 1;//gause size of dimension (gaussSizeD x gaussSizeD)-square matrix
	Core gaussCore = ImageFilters::GenerateGaussCore(wk / 3.0, wk);
	double Ix, Iy;
	double* gaussWeight = gaussCore.data;

	int gpos, iB, jB;	//B-Biased

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			A[pos] = 0;
			B[pos] = 0;
			C[pos] = 0;

			for (int u = -wk; u <= wk; u++) {
				for (int v = -wk; v <= wk; v++) {

					iB = i + u;
					jB = j + v;

					if (0 <= iB && iB < rows - 1 && 0 <= jB && jB < cols - 1) {

						Ix = partDerX[iB * cols + jB];
						Iy = partDerY[iB * cols + jB];
					}
					else {

						Ix = ImageFilters::GetVirtualPixelZero(iB, jB, rows, cols, partDerX);
						Iy = ImageFilters::GetVirtualPixelZero(iB, jB, rows, cols, partDerY);
					}

					gpos = (wk + u) * gaussSizeD + (wk + v);

					A[pos] += Ix * Ix * gaussWeight[gpos];
					B[pos] += Ix * Iy * gaussWeight[gpos];
					C[pos] += Iy * Iy * gaussWeight[gpos];
				}
			}
		}
	}

	//calculate response map
	double* responseMap = HarrisResponse(rows, cols, A, B, C, kHarrisResponseForstner);

	//show response map
	/*double* responseMapNormalized = ImageFilters::NormalizeData(rows * cols, responseMap);
	Image test(rows, cols, responseMapNormalized);
	cv::imshow("test" + std::to_string(rand()), test.GetMat());*/

	/*double* responseMapNormalized = ImageFilters::NormalizeData(gaussSizeD * gaussSizeD, gaussWeight);
	Image test((2 * gaussCore.xk + 1), (2 * gaussCore.yk + 1), responseMapNormalized);
	cv::imwrite("pyramid/gauss.png", test.GetMat());*/

	std::vector<Point> result;
	Point point;
	bool isPoint;
	int r;

	//check is some point greater then threshold AND local max
	double localValue;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			isPoint = true;
			pos = i * cols + j;

			if (responseMap[pos] > Threshold) {
				for (int u = -localMinK; u <= localMinK && isPoint; u++) {
					for (int v = -localMinK; v <= localMinK && isPoint; v++) {
						if ((localMinK >= sqrt(pow(u, 2) + pow(v, 2))) && !(u == 0 && v == 0) && responseMap[pos] <= ImageFilters::GetVirtualPixelZero(i + u, j + v, rows, cols, responseMap)) {
							isPoint = false;
						}
					}
				}
			}
			else {
				isPoint = false;
			}

			if (isPoint) {

				point.x = j;
				point.y = i;

				result.push_back(point);
			}
		}
	}

	//clear arrays
	delete[] A;
	delete[] B;
	delete[] C;

	std::vector<Point> resultANMS;

	//calculate ANMS or not
	if (pointsNeeded > 0) {

		resultANMS = ANMS(result, rows, cols, responseMap, pointsNeeded, 1);
		result.clear();
		result = resultANMS;
		resultANMS.clear();
	}

	delete[] responseMap;

	return result;
}

double * ImageDetectors::HarrisResponse(int rows, int cols, double * A, double * B, double * C, int harrisResponseType) {

	int size = rows * cols;
	double* result = new double[size];
	double k = 0.05; // for base
					 
	//in work!
	/*if (harrisResponseType == kHarrisResponseDirect) {

		delete[] result;
		result = NULL;
	}
	else */
	if (harrisResponseType == kHarrisResponseBase) {
		for (int i = 0; i < size; i++) {
			result[i] = A[i] * C[i] - pow(B[i], 2) - k * pow(A[i] + C[i], 2);
		}
	}
	else if(harrisResponseType == kHarrisResponseForstner) {

		for (int i = 0; i < size; i++) {
			if (A[i] * C[i] == 0)
				result[i] = 0;
			else
				result[i] = (A[i] * C[i] - pow(B[i], 2)) / (A[i] + C[i]);
		}
	}
	else {
		throw new invalid_argument("Unexpected harrisResponseType value");
	}

	return result;
}

//Use adaptive non-max supression
std::vector<Point> ImageDetectors::ANMS(std::vector<Point> points, int rows, int cols, double * responseMap, int pointsNeeded, double c) {

	int size = rows * cols;
	int pointsNum = points.size();

	if (points.size() <= pointsNeeded) {

		return points;
	}

	double* radius = new double[points.size()];
	double radiusTemp;
	int row1, col1, row2, col2;

	for (int i = 0; i < pointsNum; i++) {

		radius[i] = size;
		row1 = points[i].y;
		col1 = points[i].x;

		for (int j = 0; j < pointsNum; j++) {

			if (i != j) {

				row2 = points[j].y;
				col2 = points[j].x;

				radiusTemp = std::sqrt(std::pow((row1 - row2), 2) + std::pow((col1 - col2), 2));

				if (responseMap[row1 * cols + col1] < c * responseMap[row2 * cols + col2] && radius[i] > radiusTemp) {

					radius[i] = radiusTemp;
				}
			}
		}
	}

	double* tempR = new double[pointsNum];
	copy(radius, &(radius[pointsNum]), stdext::checked_array_iterator<double*>(tempR, pointsNum));

	double tempVal;
	for (int i = 0; i < pointsNum; i++) {
		for (int j = 0; j < pointsNum - 1 - i; j++) {

			if (tempR[j] < tempR[j + 1]) {
				tempVal = tempR[j];
				tempR[j] = tempR[j + 1];
				tempR[j + 1] = tempVal;
			}
		}
	}

	double radiusNeeded = tempR[pointsNeeded];
	if (pointsNeeded < pointsNum) {
		if (radiusNeeded == tempR[(pointsNeeded - 1) + 1]) {
			radiusNeeded += 0.1;
		}
	}

	delete[] tempR;

	vector<Point> result;
	for (int i = 0; i < pointsNum; i++) {
		if (radius[i] >= radiusNeeded) {
			result.push_back(points[i]);
		}
	}

	delete[] radius;

	return result;
}