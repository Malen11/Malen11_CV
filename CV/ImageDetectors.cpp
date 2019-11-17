#include "stdafx.h"
#include "ImageDetectors.hpp"

using namespace std;
using namespace CV_labs;

namespace cv {
	/** function MatPlotPoints */
	cv::Mat MatPlotPoints(Image& img, vector<CV_labs::Point> points, int color = 0) {

		cv::Mat temp = img.GetMat();
		cv::Mat colored;
		cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);

		for (vector<CV_labs::Point>::iterator it = points.begin(); it != points.end(); it++) {

			circle(colored, cv::Point((*it).x, (*it).y), 1, CV_RGB(255, 0, 0), 3);
		}

		return colored;
	}

	cv::Mat MatPlotLines(Image& img, vector<CV_labs::Line> lines, int color = 0) {

		cv::Mat temp = img.GetMat();
		cv::Mat colored;
		cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);
		RNG rng(12345);

		for (vector<CV_labs::Line>::iterator it = lines.begin(); it != lines.end(); it++) {
			Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
			line(colored, cv::Point((*it).pointA.x, (*it).pointA.y), cv::Point((*it).pointB.x, (*it).pointB.y), color, 2);
		}

		return colored;
	}

	cv::Mat MatPlotBlobs(Image& img, vector<CV_labs::ScalePoint> points, int color = 0) {

		cv::Mat temp = img.GetMat();
		cv::Mat colored;
		cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);

		for (vector<CV_labs::ScalePoint>::iterator it = points.begin(); it != points.end(); it++) {

			//проверка, что все значения допустимы (костыль, потом удалить)
			img.GetValueAt((*it).point);

			//if ((*it).point.x > 300 && (*it).point.x < 340 && (*it).point.y > 140 && (*it).point.y < 180)
			circle(colored, cv::Point((*it).point.x, (*it).point.y), sqrt(2) * (*it).scale, CV_RGB(255, 0, 0), 1);
		}

		return colored;
	}
}

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

//apply Harris to Image
std::vector<Point> ImageDetectors::Harris(Image & image, int wk, int localMinK, double Threshold, int pointsNeeded, int PartDerivativeType) {

	double* data = image.GetNormalizedImageDataD();

	//double* dataSmoothed = ImageFilters::ApplyFilterRaw(image.GetRowsNumber(), image.GetColsNumber(), data, ImageFilters::GenerateGaussSeparableCore((double)wk / 3, wk), ImageFilters::kInterpolateReflection);

	vector<Point> result = HarrisRaw(image.GetRowsNumber(), image.GetColsNumber(), data, (double)wk / 3, wk, (double)wk / 3, localMinK, Threshold, pointsNeeded, PartDerivativeType);

	delete[] data;
	//delete[] dataSmoothed;

	return result;
}

//apply Harris to Image
std::vector<Point> ImageDetectors::Harris(Image & image, double sigmaD,int wk, double sigmaI, int localMinK, double Threshold, int pointsNeeded,  int PartDerivativeType) {

	double* data = image.GetNormalizedImageDataD();

	//double* dataSmoothed = ImageFilters::ApplyFilterRaw(image.GetRowsNumber(), image.GetColsNumber(), data, ImageFilters::GenerateGaussSeparableCore((double)wk / 3, wk), ImageFilters::kInterpolateReflection);

	vector<Point> result = HarrisRaw(image.GetRowsNumber(), image.GetColsNumber(), data, sigmaD, wk, sigmaI, localMinK, Threshold, pointsNeeded, PartDerivativeType);

	delete[] data;
	//delete[] dataSmoothed;

	return result;
}

//apply Harris to Image data
std::vector<Point> ImageDetectors::HarrisRaw(int rows, int cols, double * data, double sigmaD, int wk, double sigmaI, int localMinK, double Threshold, int pointsNeeded,  int PartDerivativeType) {

	//calculate response map
	double* responseMap = HarrisResponseRaw(rows, cols, data, sigmaD, wk, sigmaI, PartDerivativeType, kHarrisResponseForstner);

	std::vector<Point> result;
	bool isPoint;
	int pos;

	//check is some point greater then threshold AND local max
	//double localValue;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			isPoint = true;
			pos = i * cols + j;

			if (responseMap[pos] > Threshold) {
				for (int u = -localMinK; isPoint && (u <= localMinK); u++) {
					for (int v = -localMinK; isPoint && (v <= localMinK); v++) {

						if (u == 0 && v == 0) {
							continue;
						}

						if (responseMap[pos] <= ImageFilters::GetVirtualPixelZero(i + u, j + v, rows, cols, responseMap)) {

							isPoint = false;
						}
					}
				}
			}
			else {
				isPoint = false;
			}

			if (isPoint) {

				result.push_back({ j, i });
			}
		}
	}

	std::vector<Point> resultANMS;

	//calculate ANMS or not
	if (pointsNeeded > 0) {

		resultANMS = ANMS(result, rows, cols, responseMap, pointsNeeded, 1);
		result.clear();
		result = resultANMS;
		resultANMS.clear();
	}

	delete[] responseMap;

	resultANMS.clear();
	resultANMS.shrink_to_fit();

	return result;
}

//calculate response for Harris
double * ImageDetectors::HarrisResponseRaw(int rows, int cols, double * data, double sigmaD, int wk, double sigmaI, int PartDerivativeType, int harrisResponseType) {

	int size = rows * cols;
	int pos;

	//part derivative X, Y
	double *partDerX = NULL, *partDerY = NULL;

	double * dataSmoothed = ImageFilters::ApplyFilterRaw(rows, cols, data, ImageFilters::GenerateGaussSeparableCore(sigmaD), ImageFilters::kInterpolateBorder);

	if (PartDerivativeType == ImageFilters::kPartDerivativeTypeSobelCore) {

		partDerX = ImageFilters::ApplyFilterRaw(rows, cols, dataSmoothed, ImageFilters::GenerateSobelSeparableCore(ImageFilters::kPartDerivativeDirectionX), ImageFilters::kInterpolateBorder);
		partDerY = ImageFilters::ApplyFilterRaw(rows, cols, dataSmoothed, ImageFilters::GenerateSobelSeparableCore(ImageFilters::kPartDerivativeDirectionY), ImageFilters::kInterpolateBorder);
	}
	else if (PartDerivativeType == ImageFilters::kPartDerivativeTypePrewittCore) {

		partDerX = ImageFilters::ApplyFilterRaw(rows, cols, dataSmoothed, ImageFilters::GeneratePrewittSeparableCore(ImageFilters::kPartDerivativeDirectionX), ImageFilters::kInterpolateBorder);
		partDerY = ImageFilters::ApplyFilterRaw(rows, cols, dataSmoothed, ImageFilters::GeneratePrewittSeparableCore(ImageFilters::kPartDerivativeDirectionY), ImageFilters::kInterpolateBorder);
	}
	else if (PartDerivativeType == ImageFilters::kPartDerivativeTypeScharrCore) {

		partDerX = ImageFilters::ApplyFilterRaw(rows, cols, dataSmoothed, ImageFilters::GenerateScharrSeparableCore(ImageFilters::kPartDerivativeDirectionX), ImageFilters::kInterpolateBorder);
		partDerY = ImageFilters::ApplyFilterRaw(rows, cols, dataSmoothed, ImageFilters::GenerateScharrSeparableCore(ImageFilters::kPartDerivativeDirectionY), ImageFilters::kInterpolateBorder);
	}

	double* A = new double[size];
	double* B = new double[size];
	double* C = new double[size];

	//gauss core for gauss weight
	int gaussSizeD = 2 * wk + 1;//gause size of dimension (gaussSizeD x gaussSizeD)-square matrix
	Core gaussCore = ImageFilters::GenerateGaussCore(sigmaI, wk);
	double* gaussWeight = gaussCore.data;

	double* Ix2 = new double[size];
	double* Iy2 = new double[size];
	double* IxIy = new double[size];

	for (int i = 0; i < size; i++) {

		Ix2[i] = partDerX[i] * partDerX[i];
		Iy2[i] = partDerY[i] * partDerY[i];
		IxIy[i] = partDerX[i] * partDerY[i];
	}

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
					gpos = (wk + u) * gaussSizeD + (wk + v);

					if (0 <= iB && iB < rows - 1 && 0 <= jB && jB < cols - 1) {

						A[pos] += Ix2[iB * cols + jB] * gaussWeight[gpos];
						B[pos] += IxIy[iB * cols + jB] * gaussWeight[gpos];
						C[pos] += Iy2[iB * cols + jB] * gaussWeight[gpos];
					}
					else {

						A[pos] += ImageFilters::GetVirtualPixelZero(iB, jB, rows, cols, Ix2) * gaussWeight[gpos];
						B[pos] += ImageFilters::GetVirtualPixelZero(iB, jB, rows, cols, IxIy) * gaussWeight[gpos];
						C[pos] += ImageFilters::GetVirtualPixelZero(iB, jB, rows, cols, Iy2) * gaussWeight[gpos];
					}
				}
			}
		}
	}

	double* result = new double[size];
					 
	if (harrisResponseType == kHarrisResponseBase) {
		for (int i = 0; i < size; i++) {
			result[i] = A[i] * C[i] - pow(B[i], 2) - 0.05 * pow(A[i] + C[i], 2);
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

	//clear arrays
	delete[] dataSmoothed;
	delete[] gaussWeight;

	delete[] Ix2, Iy2, IxIy;
	delete[] partDerX, partDerY;

	delete[] A;
	delete[] B;
	delete[] C;

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

//vector<ScalePoint> ImageDetectors::DifferenceOfGaussians(const Image & image) {
//
//	ScaleSpace imagePyramid(image, 0.5, 1.6, 5, 5, 3);
//
//	vector<ScalePoint> pointsOfInterest = ImageDetectors::DifferenceOfGaussiansRaw(imagePyramid);
//
//	return pointsOfInterest;
//}

//double* ImageDetectors::DifferenceOfGaussiansRaw(const Image &gaussian0, const Image &gaussian1) {
//	
//	int rowsNum = gaussian0.GetRowsNumber();
//	int colsNum = gaussian0.GetColsNumber();
//	
//	if (rowsNum != gaussian1.GetRowsNumber() || colsNum != gaussian1.GetColsNumber()) {
//		throw new std::invalid_argument("Image dimension not match");
//	}
//	int size = rowsNum * colsNum;
//
//	double* dataDif = new double[size];
//	double* data0 = gaussian0.GetDataD();
//	double* data1 = gaussian1.GetDataD();
//
//	for (int i = 0; i < size; i++) {
//		
//		dataDif[i] = data1[i] - data0[i];
//	}
//	/*
//	if (gaussians.GetLayersNum() + gaussians.GetCrossLayersNum() < 4) {
//		throw new std::invalid_argument("Gaussians must have at least 4 layers and cross layers");
//	}
//
//	double* dataScale0 = nullptr;
//	double* dataScale1 = nullptr;
//
//	double* dogPrevious = nullptr;
//	double* dogTarget = nullptr;
//	double* dogNext = nullptr;
//
//	vector<CV_labs::Point> extremums;
//
//	int rowsNum, colsNum, size;
//	double sigma;
//
//	for (int octave = 0; octave < gaussians.GetOctavesNum(); octave++) {
//
//		rowsNum = gaussians.GetImage(octave, 0).GetRowsNumber();
//		colsNum = gaussians.GetImage(octave, 0).GetColsNumber();
//		size = rowsNum * colsNum;
//
//		for (int layer = 1; layer < gaussians.GetLayersNum() + gaussians.GetCrossLayersNum() - 2; layer++) {
//
//			if (layer == 1) {
//
//				dataScale0 = gaussians.GetImage(octave, layer - 1).GetDataD();
//
//				dataScale1 = gaussians.GetImage(octave, layer).GetDataD();
//				dogPrevious = arrayDifference(size, dataScale1, dataScale0);
//				delete[] dataScale0;
//				dataScale0 = dataScale1;
//
//				dataScale1 = gaussians.GetImage(octave, layer + 1).GetDataD();
//				dogTarget = arrayDifference(size, dataScale0, dataScale1);
//				delete[] dataScale0;
//				dataScale0 = dataScale1;
//			}
//
//			dataScale1 = gaussians.GetImage(octave, layer + 2).GetDataD();
//			dogNext = arrayDifference(size, dataScale0, dataScale1);
//			delete[] dataScale0;
//			dataScale0 = dataScale1;
//
//			/////////////////////////////
//			//if (octave != 0) {
//			////	Image tempImage0(rowsNum, colsNum, dataScale0);
//			////	imshow("dataScale0", tempImage0.GetUpsampledImage(std::pow(2, octave)).GetNormalizedImage().GetMat());
//			////	Image tempImage1(rowsNum, colsNum, dataScale1);
//			////	imshow("dataScale1", tempImage1.GetUpsampledImage(std::pow(2, octave)).GetNormalizedImage().GetMat());
//			//	
//			////	double max = dogPrevious[0], min = dogPrevious[0];
//			////	for (int i = 0; i < size; i++) {
//			////		if (max < dogPrevious[i])
//			////			max = dogPrevious[i];
//			////		if (min > dogPrevious[i])
//			////			min = dogPrevious[i];
//			////	}
//			//	
//			////	Image tempImageP(rowsNum, colsNum, ImageFilters::NormalizeData(size, dogPrevious));
//			////	imshow("dogPrevious", tempImageP.GetMat());
//			////	Image tempImageT(rowsNum, colsNum, ImageFilters::NormalizeData(size, dogTarget));
//			////	imshow("dogTarget", tempImageT.GetMat());
//			////	Image tempImageN(rowsNum, colsNum, ImageFilters::NormalizeData(size, dogNext));
//			////	imshow("dogNext", tempImageN.GetMat());
//			////	cvWaitKey();
//			//}
//			//////////////////////////////
//			extremums = findExtremumIn3D(rowsNum, colsNum, dogTarget, dogPrevious, dogNext);
//			sigma = gaussians.GetImageSigma(octave, layer);
//
//			for (int pos = 0; pos < extremums.size(); pos++) {
//				pointsOfInterest.push_back({ gaussians.RestoreCoordinate(extremums[pos], octave, layer), sigma });
//			}
//
//			extremums.clear();
//			
//			delete[] dogPrevious;
//			dogPrevious = dogTarget;
//			dogTarget = dogNext;
//			dogNext = nullptr;
//		}
//
//		delete[] dogTarget, dogNext, dataScale0, dataScale1;
//
//		////////
//		//break;
//		///////
//		
//	}*/
//
//	return dataDif;
//}

//std::vector<Point> CV_labs::ImageDetectors::DifferenceOfGaussians(const vector<Image> & gaussians) {
//	
//	vector<Point> pointsOfInterest;
//
//	if (gaussians.size() < 4) {
//		throw new std::invalid_argument("Vector must contain at least 4 elements");
//	}
//
//	double* dataScale0 = nullptr;
//	double* dataScale1 = nullptr;
//
//	double* dogPrevious = nullptr;
//	double* dogTarget = nullptr;
//	double* dogNext = nullptr;
//
//	vector<CV_labs::Point> extremums;
//
//	int rowsNum = gaussians[0].GetRowsNumber();
//	int colsNum = gaussians[0].GetRowsNumber();
//	int size = rowsNum * colsNum;
//
//	for (int index = 1; index <= gaussians.size() - 1; index++) {
//
//			if (index == 1) {
//
//				dataScale0 = gaussians[index - 1].GetDataD();
//
//				dataScale1 = gaussians[index].GetDataD();
//				dogPrevious = arrayDifference(size, dataScale1, dataScale0);
//				delete[] dataScale0;
//				dataScale0 = dataScale1;
//
//				dataScale1 = gaussians[index + 1].GetDataD();
//				dogTarget = arrayDifference(size, dataScale0, dataScale1);
//				delete[] dataScale0;
//				dataScale0 = dataScale1;
//			}
//
//			dataScale1 = gaussians[index + 2].GetDataD();
//			dogNext = arrayDifference(size, dataScale0, dataScale1);
//			delete[] dataScale0;
//			dataScale0 = dataScale1;
//
//			extremums = findExtremumIn3D(rowsNum, colsNum, dogTarget, dogPrevious, dogNext);
//
//			for (int pos = 0; pos < extremums.size(); pos++) {
//				pointsOfInterest.push_back(extremums[pos]);
//			}
//
//			extremums.clear();
//
//			delete[] dogPrevious;
//			dogPrevious = dogTarget;
//			dogTarget = dogNext;
//			dogNext = nullptr;
//		}
//
//		delete[] dogTarget, dogNext, dataScale0, dataScale1;
//
//		////////
//		//break;
//		////////
//
//	return pointsOfInterest;
//}

//find points of interest by via DoG and Harris
std::vector<ScalePoint> CV_labs::ImageDetectors::FindPointsOfInterest(const ScaleSpace & imagePyramid, int wk, int localMinK, double harrisThreshold, double dogThreshold) {

	int /*wk = 2, localMinK = 2,*/ pointsNeeded = 1;
	//double threshold = 0.01;

	vector<ScalePoint> pointsOfInterest;

	vector<Point> dogPointsTemp;
	vector<Point> harrisPoints;
	double sigma;

	double* dogPrev, *dogTarget, *dogNext;

	for (int octaveIndex = 0; octaveIndex < imagePyramid.GetOctavesNum(); octaveIndex++) {
		for (int layerIndex = 0; layerIndex < imagePyramid.GetLayersNum(); layerIndex++) {

			Image imageTemp = imagePyramid.GetImage(octaveIndex, layerIndex);
			int cols = imageTemp.GetColsNumber();
			int rows = imageTemp.GetRowsNumber();
			int size = cols * rows;

			vector<Point> harrisPointsTemp = Harris(imageTemp, wk, localMinK, harrisThreshold, pointsNeeded);

			dogPrev		= imagePyramid.DoG(octaveIndex, layerIndex - 1);//ImageFilters::NormalizeData(size, DifferenceOfGaussiansRaw(imagePyramid.GetImage(octaveIndex, layerIndex - 1), imageTemp));
			dogTarget	= imagePyramid.DoG(octaveIndex, layerIndex);//ImageFilters::NormalizeData(size, DifferenceOfGaussiansRaw(imageTemp, imagePyramid.GetImage(octaveIndex, layerIndex + 1)));
			dogNext		= imagePyramid.DoG(octaveIndex, layerIndex + 1);//ImageFilters::NormalizeData(size, DifferenceOfGaussiansRaw(imagePyramid.GetImage(octaveIndex, layerIndex + 1), imagePyramid.GetImage(octaveIndex, layerIndex + 2)));
			
			vector<Point> extremumIn3DPointsTemp = ImageDetectors::ExtremumIn3DRaw(
				imageTemp.GetRowsNumber(),
				imageTemp.GetColsNumber(),
				dogPrev,
				dogTarget,
				dogNext,
				dogThreshold
			);

			for (int harrisPointIndex = 0; harrisPointIndex < harrisPointsTemp.size(); harrisPointIndex++) {
				for (int exptremumPointIndex = 0; exptremumPointIndex < extremumIn3DPointsTemp.size(); exptremumPointIndex++) {

					if (harrisPointsTemp[harrisPointIndex].x == extremumIn3DPointsTemp[exptremumPointIndex].x && 
						harrisPointsTemp[harrisPointIndex].y == extremumIn3DPointsTemp[exptremumPointIndex].y) {

						double sigma = imagePyramid.GetImageSigma(octaveIndex, layerIndex);
						pointsOfInterest.push_back({ imagePyramid.RestoreCoordinate(harrisPointsTemp[harrisPointIndex], octaveIndex, layerIndex), sigma });
						break;
					}
				}
			}

			/////////////////////
			/*string name = "image " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex - 1) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", imageTemp.GetMat());

			Image prev(
				rows,
				cols,
				ImageFilters::NormalizeData(size, DifferenceOfGaussiansRaw(imagePyramid.GetImage(octaveIndex, layerIndex - 1), imageTemp))
			);
			name = "dog " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex - 1) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", prev.GetMat());

			Image target(
				rows,
				cols,
				ImageFilters::NormalizeData(size, DifferenceOfGaussiansRaw(imageTemp, imagePyramid.GetImage(octaveIndex, layerIndex + 1)))
			);
			name = "dog " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", target.GetMat());

			Image next(
				rows,
				cols,
				ImageFilters::NormalizeData(size, DifferenceOfGaussiansRaw(imagePyramid.GetImage(octaveIndex, layerIndex + 1), imagePyramid.GetImage(octaveIndex, layerIndex + 2)))
			);
			name = "dog " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex + 1) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", next.GetMat());

			//cv::imshow("prev ", prev.GetNormalizedImage().GetMat());
			//cv::imshow("target ", target.GetNormalizedImage().GetMat());
			//cv::imshow("next ", next.GetNormalizedImage().GetMat());

			//cv::imshow("harris", cv::MatPlotPoints(imageTemp, harrisPointsTemp));
			name = "harris " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", cv::MatPlotPoints(imageTemp, harrisPointsTemp));

			//cv::imshow("blob's centers", cv::MatPlotPoints(imageTemp, extremumIn3DPointsTemp));
			name = "blob " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", cv::MatPlotPoints(imageTemp, extremumIn3DPointsTemp));

			//cv::imshow("harris-dog", cv::MatPlotBlobs(imageTemp, pointsOfInterest));

			//cv::waitKey();*/
			////////////////////

			harrisPointsTemp.clear();
			harrisPointsTemp.shrink_to_fit();
			extremumIn3DPointsTemp.clear();
			extremumIn3DPointsTemp.shrink_to_fit();
			delete[] dogPrev, dogTarget, dogNext;
		}
	}

	//vector<ScalePoint> pointsOfInterest;
	//
	//double eps = 0.0000001;
	//
	//for (int dogIndex = 0; dogIndex < dogPoints.size(); dogIndex++) {
	//	for (int harrisIndex = 0; harrisIndex < harrisPoints.size(); harrisIndex++) {
	//
	//		if (std::abs(dogPoints[dogIndex].scale - harrisPoints[harrisIndex].scale) < eps) {
	//			if (dogPoints[dogIndex].point.x == harrisPoints[harrisIndex].point.x && dogPoints[dogIndex].point.y == harrisPoints[harrisIndex].point.y) {
	//				pointsOfInterest.push_back(dogPoints[dogIndex]);
	//			}
	//		}
	//	}
	//}

	return pointsOfInterest;
}

//find points of interest by via DoG and Harris
std::vector<ScalePoint> CV_labs::ImageDetectors::HarrisLaplass(const ScaleSpace & imagePyramid, int wk, int localMinK, double harrisThreshold, double dogThreshold, int pointsNeeded) {

	vector<ScalePoint> pointsOfInterest;

	double* dogPrev, *dogTarget, *dogNext;
	double sigma;
	int pos;

	for (int octaveIndex = 0; octaveIndex < imagePyramid.GetOctavesNum(); octaveIndex++) {
		for (int layerIndex = 0; layerIndex < imagePyramid.GetLayersNum(); layerIndex++) {

			Image imageTemp = imagePyramid.GetImage(octaveIndex, layerIndex);
			double sigma = imagePyramid.GetImageSigma(octaveIndex, layerIndex);

			double * data = imageTemp.GetDataD();
			int cols = imageTemp.GetColsNumber();
			int rows = imageTemp.GetRowsNumber();
			int size = cols * rows;
			
			//calculate response map
			double* responseMap = HarrisResponseRaw(rows, cols, data, wk / 3.0, wk, wk / 3.0);

			std::vector<Point> harrisPoints;
			std::vector<Point> harrisPointsTemp;
			bool isPoint;
			int pos;

			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {

					isPoint = true;
					pos = i * cols + j;

					if (responseMap[pos] > harrisThreshold) {
						for (int u = -localMinK; isPoint && (u <= localMinK); u++) {
							for (int v = -localMinK; isPoint && (v <= localMinK); v++) {

								if (u == 0 && v == 0) {
									continue;
								}

								if (responseMap[pos] <= ImageFilters::GetVirtualPixelZero(i + u, j + v, rows, cols, responseMap)) {

									isPoint = false;
								}
							}
						}
					}
					else {
						isPoint = false;
					}

					if (isPoint) {

						harrisPoints.push_back({ j, i });
					}
				}
			}

			//calculate DoG
			dogPrev = imagePyramid.DoG(octaveIndex, layerIndex - 1);
			dogTarget = imagePyramid.DoG(octaveIndex, layerIndex);
			dogNext = imagePyramid.DoG(octaveIndex, layerIndex + 1);

			for (int harrisPointIndex = 0; harrisPointIndex < harrisPoints.size(); harrisPointIndex++) {

				pos = harrisPoints[harrisPointIndex].y * cols + harrisPoints[harrisPointIndex].x;

				if ((std::abs(dogTarget[pos]) > dogThreshold) && (dogTarget[pos] > dogPrev[pos]) && (dogTarget[pos] > dogNext[pos])) {

					harrisPointsTemp.push_back(harrisPoints[harrisPointIndex]);
				}
			}

			harrisPoints.clear();
			harrisPoints.shrink_to_fit();

			if (pointsNeeded > 0) {

				harrisPoints = ANMS(harrisPointsTemp, rows, cols, responseMap, pointsNeeded, 1);
			}
			else {
				harrisPoints = harrisPointsTemp;
			}
			
			for (int harrisPointIndex = 0; harrisPointIndex < harrisPoints.size(); harrisPointIndex++) {

				//pointsOfInterest.push_back({ harrisPoints[harrisPointIndex], sigma });
				pointsOfInterest.push_back({ harrisPoints[harrisPointIndex], sigma });
			}

			/////////////////////
			/*string name = "image " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", imageTemp.GetMat());

			Image prev(
				rows,
				cols,
				ImageFilters::NormalizeData(size, dogPrev)
			);
			name = "dog " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex - 1) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", prev.GetMat());

			Image target(
				rows,
				cols,
				ImageFilters::NormalizeData(size, dogTarget)
			);
			name = "dog " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", target.GetMat());

			Image next(
				rows,
				cols,
				ImageFilters::NormalizeData(size, dogNext)
			);
			name = "dog " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex + 1) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", next.GetMat());

			//cv::imshow("prev ", prev.GetNormalizedImage().GetMat());
			//cv::imshow("target ", target.GetNormalizedImage().GetMat());
			//cv::imshow("next ", next.GetNormalizedImage().GetMat());

			//cv::imshow("harris", cv::MatPlotPoints(imageTemp, harrisPointsTemp));
			name = "harris " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", cv::MatPlotPoints(imageTemp, harrisPointsTemp));

			//cv::imshow("blob's centers", cv::MatPlotPoints(imageTemp, extremumIn3DPointsTemp));
			name = "blob " + std::to_string(octaveIndex) + " layer " + std::to_string(layerIndex) + " rand " + std::to_string(std::rand());
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", cv::MatPlotPoints(imageTemp, harrisPoints));

			//cv::imshow("harris-dog", cv::MatPlotBlobs(imageTemp, pointsOfInterest));

			//cv::waitKey();*/
			////////////////////

			harrisPoints.clear();
			harrisPoints.shrink_to_fit();
			harrisPointsTemp.clear();
			harrisPointsTemp.shrink_to_fit();

			delete[] dogPrev, dogTarget, dogNext, responseMap, data;
		}
	}

	return pointsOfInterest;
}

//find points of interest by via DoG and Harris
//std::vector<ScalePoint> CV_labs::ImageDetectors::FindPointsOfInterest(const Image & image) {
//	
//	ScaleSpace imagePyramid(image, 0.5, 1.6, 5, 5, 3);
//
//	vector<ScalePoint> dogPoints = ImageDetectors::DifferenceOfGaussiansRaw(imagePyramid);
//
//	vector<ScalePoint> harrisPoints;
//	
//	vector<Point> harrisPointsTemp;
//	double sigma;
//
//	for (int octaveIndex = 0; octaveIndex < imagePyramid.GetOctavesNum(); octaveIndex++) {
//		for (int layerIndex = 0; layerIndex < imagePyramid.GetLayersNum() + imagePyramid.GetCrossLayersNum(); layerIndex++) {
//
//			Image imageTemp = imagePyramid.GetImage(octaveIndex, layerIndex);
//			harrisPointsTemp = Harris(imageTemp, 5, 5, 0.001, 100);
//			sigma = imagePyramid.GetImageSigma(octaveIndex, layerIndex);
//
//			for (int i = 0; i < harrisPointsTemp.size(); i++) {
//
//				harrisPoints.push_back({ harrisPointsTemp[i], sigma });
//			}
//		}
//	}
//
//	vector<ScalePoint> pointsOfInterest;
//
//	double eps = 0.0000001;
//
//	for (int dogIndex = 0; dogIndex < dogPoints.size(); dogIndex++) {
//		for (int harrisIndex = 0; harrisIndex < harrisPoints.size(); harrisIndex++) {
//
//			if (std::abs(dogPoints[dogIndex].scale - harrisPoints[harrisIndex].scale) < eps) {
//				if (dogPoints[dogIndex].point.x == harrisPoints[harrisIndex].point.x && dogPoints[dogIndex].point.y == harrisPoints[harrisIndex].point.y) {
//					pointsOfInterest.push_back(dogPoints[dogIndex]);
//				}
//			}
//		}
//	}
//
//	return pointsOfInterest;
//}

//calculate difference of arrays
double * ImageDetectors::arrayDifference(int size, const double * const data0, const double * const data1) {

	double* arrDif = new double[size];

	for (int i = 0; i < size; i++) {

		arrDif[i] = data0[i] - data1[i];
	}

	return arrDif;
}

//find extremum in 3d (for DoG)
vector<Point> ImageDetectors::ExtremumIn3DRaw(int rows, int cols, const double* const dataPrevious, const double* const dataTarget, const double* const dataNext, double threshold) {
	
	vector<Point> extremums;
	double val;
	bool maxFlag;
	bool minFlag;
	int pos;
	const double** data = new const double *[3];
	data[0] = dataPrevious;
	data[1] = dataTarget;
	data[2] = dataNext;

	double eps = dataTarget[0];

	for (int row = 1; row < rows - 1; row++) {
		for (int col = 1; col < cols - 1; col++) {

			val = dataTarget[row * cols + col];
				if (eps < val)
					eps = val;
			if (std::abs(val) < threshold) {
				continue;
			}

			maxFlag = true;
			minFlag = true;

			if (row >= 47 && row <= 47 && col >= 57 && col <= 59)//38-74, 39-74???
				int stop = 1;
			for (int index = 0; (index < 3) && (maxFlag || minFlag); index++) {
				for (int rowOff = row - 1; (rowOff <= row + 1) && (maxFlag || minFlag); rowOff++) {
					for (int colOff = col - 1; (colOff <= col + 1) && (maxFlag || minFlag); colOff++) {

						if (index == 1 && rowOff == row && colOff == col) {
							continue;
						}

						if (maxFlag && ImageFilters::GetVirtualPixelZero(rowOff, colOff, rows, cols, data[index]) - val >= 0) {
							maxFlag = false;
						}

						if (minFlag && ImageFilters::GetVirtualPixelZero(rowOff, colOff, rows, cols, data[index]) - val <= 0) {
							minFlag = false;
						}
					}
				}
			}

			if (maxFlag || minFlag) {
				extremums.push_back({ col, row });
			}
		}
	}

	delete[] data;

	return extremums;
}
