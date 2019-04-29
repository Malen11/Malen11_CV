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

std::vector<Descriptor> ComputerVision::CalculateDescriptors(Image & img, std::vector<Dot> points, int gridSizeK, int histogramSize, int intervalsNum) {
	
	return std::vector<Descriptor>();
}

Histogram ComputerVision::CreateHistogram(int row0, int row1, int col0, int col1, int rows, int cols, double * partDerivX, double * partDerivY, double * weight, int intervalsNum) {
	
	Histogram result;
	result.intervalsNum = intervalsNum;
	result.intervals = new Interval[intervalsNum];
	
	double step = (1 - (-1)) / intervalsNum;	//because val in[-1; 1]
	double startValue = -1 + step / 2;		//because we use mid value (half step)

	for (int i = 0; i < intervalsNum; i++) {
		
		result.intervals[i].data = 0;
		result.intervals[i].midVal = startValue + i * step;
	}

	double Th, G, k;
	int pos, column0, column1;
	for (int i = row0; i <= row1; i++)	{
		for (int j = col0; j <= col1; j++) {

			pos = i * cols + j;
			column0 = -1;			//test
			column1 = -1;			//test
			k = 0;					//test
			G = std::sqrt(std::pow(partDerivX[pos], 2) + std::pow(partDerivY[pos], 2));		//gradient value
			Th = atan2(partDerivY[pos], partDerivX[pos]) / PI();							//gradient vector


			//edged value
			if (Th >= result.intervals[intervalsNum - 1].midVal|| Th < result.intervals[0].midVal) {

				column0 = intervalsNum - 1;
				column1 = 0;

				if (Th < 0) {

					k = (Th - result.intervals[0].midVal)/(step/2);							//!!!
				}
				else {

					k = 1 - (result.intervals[intervalsNum - 1].midVal - Th) / (step / 2);		//!!!
				}
			}
			else {
				//find nearest column (basket)
				for (int b = 0; b < intervalsNum - 1; b++) {

					if (Th >= result.intervals[b].midVal && Th < result.intervals[(b + 1)].midVal) {

						column0 = b;
						column1 = (b + 1);
					}
				}

				k = 1 - (Th - result.intervals[column0].midVal) / (result.intervals[column1].midVal - result.intervals[column0].midVal);	//!!!
			}

			result.intervals[column0].midVal += G * k;
			result.intervals[column1].midVal += G * (k-1);
		}
	}

	return result;
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
