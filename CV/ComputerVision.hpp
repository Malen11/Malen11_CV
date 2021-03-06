﻿#pragma once
#include "stdafx.h"
//
////definition 
//
//
////Point
//struct Dot {
//
//	int x;
//	int y;
//};
//
////Gradients and angles array
//struct GradientsAnglesArray {
//
//	double* gradients;
//	double* angles;
//};
//
////Point
//struct PairDot {
//
//	Dot point0;
//
//	Dot point1;
//};
//
////Column (for Histogram)
//struct Interval {
//
//	double data;
//	double midVal;
//};
//
////Histogram (for descriptor)
//struct Histogram {
//
//	int intervalsNum;
//	Interval* intervals;
//};
//
////Descriptor
//struct Descriptor {
//
//	Dot point;
//	int featuresNum;
//	double* features;
//};
//
//
////class, that contains CV functions
//class ComputerVision {
//private:
//public:
//
//	//static const int kHarrisResponseDirect = 401;		//response for Harris (Direct by calculate own numbers)
//	static const int kHarrisResponseBase = 402;			//response for Harris (Base with k)
//	static const int kHarrisResponseForstner = 403;		//response for Harris (Förstner and Gülch)
//
//	static const int kDescriptorSimple = 501;			//Descriptor type (simple square descriptor)
//
//	static const int kDescriptorNormalization2Times = 601;	//Descriptor normalization type (normalize, thresh, normalize)
//	static const int kDescriptorNormalizationSimple = 602;	//Descriptor normalization type (normalize 1 time only)
//	static const int kDescriptorNormalizationNone = 603;	//Descriptor normalization type (not normalize)
//
//	static const int kDescriptorsComparisonEuclid = 701;	//Descriptor comparison type (Euclid distance)
//	static const int kDescriptorsComparisonManhattan = 702;	//Descriptor comparison type (Manhattan metric)
//	static const int kDescriptorsComparisonSSD = 703;		//Descriptor comparison type (Sum of squared distances)
//
//	static const int kDescriptorsMatchingBase = 701;	//Descriptor matching type (by max matching, have multiselection!)
//	static const int kDescriptorsMatchingNNDR = 702;	//Descriptor matching type (Next Nearest Distance Ratio)
//	static const int kDescriptorsMatchingMutal = 703;	//Descriptor matching type (Mutal choice)
//
//
//	template<typename srcT, typename dstT>
//
//	static dstT * LinearNormalization(int dataSize, srcT* data, dstT newMin, dstT newMax) {
//
//		if (dataSize != 0 && data != NULL) {
//
//			double min = data[0], max = data[0], k;
//
//			for (int i = 0; i < dataSize; i++) {
//				if (data[i] > max) {
//					max = data[i];
//				}
//
//				if (data[i] < min) {
//					min = data[i];
//				}
//			}
//
//			k = (newMax - newMin) / (max - min);
//
//			dstT* result = new dstT[dataSize];
//
//			for (int i = 0; i < dataSize; i++) {
//
//				result[i] = (dstT)((data[i] - min)*k + newMin);
//			}
//
//			return result;
//		}
//		else {
//			return NULL;
//		}
//	}
//
//	static std::vector<Dot> ANMS(std::vector<Dot> points, int rows, int cols, double* responseMap, int pointsNeeded, double c=0.9);
//
//	//apply Canny to Image
//	static Image Canny(Image& img, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType = kInterpolateZero);
//	
//	//apply Canny to Image data (or array) 
//	template <typename T>
//	static double* CannyRaw(int rows, int cols, T* data, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType = kInterpolateZero);
//	
//	//create Descriptor
//	static std::vector<Descriptor> CreateDescriptors(Image& img, std::vector<Dot> points, int DescriptorSizeX, int DescriptorSizeY, int histogramNumX, int histogramNumY, int intervalsNum, int descriptorType = kDescriptorSimple, int descriptorNormalizationType = kDescriptorNormalizationNone, int PartDerivativeType = kPartDerivativeSobel);
//
//	//create Descriptors
//	template <typename T>
//	static std::vector<Descriptor> CreateDescriptorsRaw(int rows, int cols, T * data, std::vector<Dot> points, int DescriptorSizeX, int DescriptorSizeY, int histogramNumX, int histogramNumY, int intervalsNum, int descriptorType = kDescriptorSimple, int descriptorNormalizationType = kDescriptorNormalizationNone, int partDerivativeType = kPartDerivativeSobel);
//	
//	//create Simple Descriptor
//	static Descriptor CreateSimpleDescriptorRaw(Dot point, int rows, int cols, double* gradients, double* angles, int histogramNumX, int histogramNumY, int intervalsNum, double alpha = 0);
//
//	//calculate descriptor orientation's angles
//	static std::vector<double> CalculateDescriptorOrientationAnglesRaw(int rows, int cols, double* gradientsWeighted, double* angles, int intervalsNum, double thresh=0.8, int alphaNumMax =2);
//
//	//calculate points' gradient AND angle on part of image
//	static GradientsAnglesArray CalculateGradientsFullRaw(int rows, int cols, double* partDerX, double* partDerY);
//
//	//calculate points' gradient AND angle on part of image
//	static GradientsAnglesArray CalculateGradientsFullRaw(int rows, int cols, double* partDerX, double* partDerY, Dot PointLT, Dot PointRB);
//	
//	//create Histogram
//	static Histogram CreateHistogramRaw(std::vector<double> values, std::vector<double> positions, double intervalsMinVal, double intervalsMaxVal,  int intervalsNum);
//
//
//	//calculate distanse beatween 2 descriptors 
//	static double DescriptorsDifference(Descriptor desc0, Descriptor desc1, int descriptorsComparisonType = kDescriptorsComparisonEuclid);
//	
//	//matching descriptors from descs1 to descs2
//	static int* DescriptorsMatchingRaw(std::vector<Descriptor> descs1, std::vector<Descriptor> descs2, int descriptorsComparisionType = kDescriptorsComparisonEuclid,int descriptorsMatchingType=kDescriptorsMatchingBase, double thresh=0);
//
//	//matching descriptors from descs1 to descs2
//	static std::vector<PairDot> DescriptorsMatching(std::vector<Descriptor> descs1, std::vector<Descriptor> descs2, int descriptorsComparisionType = kDescriptorsComparisonEuclid, int descriptorsMatchingType = kDescriptorsMatchingBase, double thresh = 0);
//
//
//	
//	//apply Harris to Image
//	static std::vector<Dot> Harris(Image& img, int wk, int localMinK, double Threshold, int pointsNeeded = -1, int PartDerivativeType = kPartDerivativeSobel);
//
//	//apply Harris to Image data (or array) 
//	template <typename T>
//	static std::vector<Dot> HarrisRaw(int rows, int cols, T* data, int wk, int localMinK, double Threshold, int pointsNeeded = -1, int PartDerivativeType= kPartDerivativeSobel);
//
//	//calculate response for Harris
//	static double* HarrisResponse(int rows, int cols, double* A, double* B, double* C, int harrisResponseType = kHarrisResponseBase);
//
//	template<typename srcT, typename dstT>
//	static dstT* LinearNormalization(int dataSize, srcT* data, dstT newMin, dstT newMax);	//функция для линейной нормализации
//	
//	//apply Moravec to Image data (or array) //apply Sobel to Image
//	static std::vector<Dot> Moravec(Image& img, int wk, int d, int p, double Threshold);
//
//	//apply Moravec to Image data (or array) 
//	template <typename T>
//	static std::vector<Dot> MoravecRaw(int rows, int cols, T* data, int wk, int d, int p, double Threshold);
//
//	///////////////////////////////////////////////
//
//	//need modify!
//	template <typename T>
//	static T* ExpandImgByZero(int rows, int cols, T* data, int addX, int addY);
//};
//
/////////////////////////////////////
//// template function realization //
/////////////////////////////////////
//
//
///////////In work!!!
//template<typename T>
//double * ComputerVision::CannyRaw(int rows, int cols, T * data, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType) {
//	
//	double* temp;
//	double* gaussImage = GaussRaw(rows, cols, data, sigma, k, interpolateType);
//
//	//double* G = SobelRaw(rows, cols, data, kPartDerivativeSobel);
//	double* Gx = PartDerivative(rows, cols, gaussImage, kPartDerivativeXSobel);
//	double* Gy = PartDerivative(rows, cols, gaussImage, kPartDerivativeYSobel);
//
//	int size = rows * cols;
//	double* G = new double[size];
//	int* Th = new int[size];
//
//	for (int i = 0; i < size; i++) {
//
//		G[i] = sqrt(Gx[i] * Gx[i] + Gy[i] * Gy[i]);
//		Th[i] = 45 * (int)round(4 * atan2(Gy[i], Gx[i]) / PI());
//	}
//
//	bool* extremum = new bool[size];
//	int pos;
//	double neighbourA, neighbourB;
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//
//
//			pos = i * cols + j;
//
//			if (i == 0 || i == rows - 1 || j == 0 || j == cols - 1) {
//				extremum[pos] = false;
//			}
//			else {
//				pos = i * cols + j;
//				neighbourA = 0;
//				neighbourB = 0;
//
//				switch (Th[pos]) {
//
//				case 0:
//				case 180:
//				case -180:
//
//					if (j != 0)
//						neighbourA = G[pos - 1];
//					if (j != cols - 1)
//						neighbourB = G[pos + 1];
//					break;
//
//				case 90:
//				case -90:
//
//					if (i != 0)
//						neighbourA = G[pos - cols];
//					if (i != rows - 1)
//						neighbourB = G[pos + cols];
//					break;
//
//				case 45:
//				case -135:
//
//					if (i != 0 && j != cols - 1)
//						neighbourA = G[pos - cols + 1];
//					if (i != rows - 1 && j != 0)
//						neighbourA = G[pos + cols - 1];
//					break;
//
//				case -45:
//				case 135:
//
//					if (i != 0 && j != 0)
//						neighbourA = G[pos - cols - 1];
//					if (i != rows - 1 && j != cols - 1)
//						neighbourA = G[pos + cols + 1];
//					break;
//
//				default:
//					break;
//				}
//
//				extremum[pos] = (G[pos] > neighbourA && G[pos] > neighbourB);
//			}
//		}
//	}
//
//	temp = LinearNormalization<double, double>(size, G, 0, 100);
//
//	double* result = new double[size];
//	for (int i = 0; i < size; i++) {
//
//		if (extremum[i] && temp[i] > highThreshhold)
//			result[i] = 1;
//		else
//			result[i] = 0;
//	}
//
//	return result;
//}
//
//template<typename T>
//std::vector<Descriptor> ComputerVision::CreateDescriptorsRaw(int rows, int cols, T * data, std::vector<Dot> points, int windowSizeX, int windowSizeY, int histogramNumX, int histogramNumY, int intervalsNum, int descriptorType, int descriptorNormalizationType, int partDerivativeType) {
//
//	std::vector<Descriptor> descsVec;
//
//	double *partDerX = NULL, *partDerY = NULL;
//
//	if (partDerivativeType == kPartDerivativeSobel) {
//
//		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXSobel, kInterpolateBorder);
//		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYSobel, kInterpolateBorder);
//	}
//	else if (partDerivativeType == kPartDerivativePrewitt) {
//
//		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXPrewitt, kInterpolateBorder);
//		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYPrewitt, kInterpolateBorder);
//	}
//	else if (partDerivativeType == kPartDerivativeScharr) {
//
//		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXScharr, kInterpolateBorder);
//		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYScharr, kInterpolateBorder);
//	}
//
//	Dot pointLT, pointRB;
//
//	//if create simple descriptor
//	if (descriptorType = kDescriptorSimple) {
//
//		//check dimensions
//		if (windowSizeX != windowSizeY)
//			throw std::invalid_argument("windowSizeX must match windowSizeY for descriptorType = kDescriptorSimple");
//
//		//for every point create descriptor
//		for (int i = 0; i < points.size(); i++) {
//
//			//create window
//			pointLT.x = points[i].x - windowSizeX / 2;
//			pointLT.y = points[i].y - windowSizeY / 2;
//
//			pointRB.x = pointLT.x + windowSizeX - 1;
//			pointRB.y = pointLT.y + windowSizeY - 1;
//
//			//calculate gradient values and angles
//			GradientsAnglesArray gaaTemp = CalculateGradientsFullRaw(rows, cols, partDerX, partDerY, pointLT, pointRB);
//			double* gradients = gaaTemp.gradients;
//			double* angles = gaaTemp.angles;
//
//			//calculate gauss weight for this part of image with 2d-size = rowsNew x rowsNew
//			double* gaussWeight = GetGaussWeight(windowSizeX / 6.0, windowSizeX);//!!! is sigma correct?
//
//			//calculate weighted gradient
//			double* gradientsWeighted = new double[windowSizeX*windowSizeX];
//			for (int j = 0; j <= windowSizeX * windowSizeX; j++)
//				gradientsWeighted[j] = gaussWeight[j] * gradients[j];
//
//			//calculate dominating angles
//			std::vector<double> alpha = CalculateDescriptorOrientationAnglesRaw(windowSizeX, windowSizeX, gradientsWeighted, angles, 36);
//			
//			//for each dominating angle calculate point descriptor
//			for (int j = 0; j < alpha.size(); j++) {
//				
//				Descriptor desc = CreateSimpleDescriptorRaw(points[i], pointRB.y - pointLT.y + 1, pointRB.x - pointLT.x + 1, gradientsWeighted, angles, histogramNumX, histogramNumY, intervalsNum, alpha[j]);
//				descsVec.push_back(desc);
//				/*
//				double* gradientsRotated = new double[windowSizeX*windowSizeX];
//				double* anglesRotated = new double[windowSizeX*windowSizeX];
//				for (int rs = 0; rs < windowSizeX*windowSizeX; rs++) {
//					gradientsRotated[rs] = 0;
//					anglesRotated[rs] = 0;
//				}
//
//				double colMid = windowSizeX / 2.0, rowMid = windowSizeX / 2.0;
//
//				//int rowNew, colNew;
//				int rowOld, colOld;
//
//				for (int ri = 0; ri < windowSizeX; ri++) {
//					for (int rj = 0; rj < windowSizeX; rj++) {
//
//						colOld = std::round(std::cos(alpha[j])*(rj - colMid) - std::sin(alpha[j])*(ri - rowMid) + colMid);
//						rowOld = std::round(std::sin(alpha[j])*(rj - colMid) + std::cos(alpha[j])*(ri - rowMid) + rowMid);
//
//						if (colOld >= 0 && colOld < windowSizeX && rowOld >= 0 && rowOld < windowSizeX) {
//
//							gradientsRotated[ri*windowSizeX + rj] = gradients[rowOld*windowSizeX + colOld];
//							anglesRotated[ri*windowSizeX + rj] = angles[rowOld*windowSizeX + colOld] - alpha[j];
//
//							if (anglesRotated[ri*windowSizeX + rj] < 0)
//								anglesRotated[ri*windowSizeX + rj] += 2 * PI();
//						}
//					}
//				}
//
//				Descriptor desc = CreateSimpleDescriptorRaw(points[i], pointRB.y - pointLT.y + 1, pointRB.x - pointLT.x + 1, gradientsRotated, anglesRotated, histogramNumX, histogramNumY, intervalsNum);
//				descsVec.push_back(desc);*/
//			}
//			delete[] gradients;
//			delete[] angles;
//			delete[] gaussWeight;
//		}
//	}
//	else {
//		throw std::invalid_argument("Wrong descriptorType");
//	}
//
//	delete[] partDerX;
//	delete[] partDerY;
//
//	if (descriptorNormalizationType == kDescriptorNormalization2Times) {
//
//		double* temp;
//		for (int i = 0; i < descsVec.size(); i++) {
//
//			temp = LinearNormalization<double, double>(descsVec[i].featuresNum, descsVec[i].features, 0, 1);
//			delete[] descsVec[i].features;
//			descsVec[i].features = temp;
//
//			for (int j = 0; j < descsVec[i].featuresNum; j++) {
//				if (descsVec[i].features[j] < 0.2)
//					descsVec[i].features[j] = 0;
//			}
//
//			temp = LinearNormalization<double, double>(descsVec[i].featuresNum, descsVec[i].features, 0, 1);
//			delete[] descsVec[i].features;
//			descsVec[i].features = temp;
//		}
//	}
//	else if (descriptorNormalizationType == kDescriptorNormalizationSimple) {
//
//		double* temp;
//		for (int i = 0; i < descsVec.size(); i++) {
//
//			temp = LinearNormalization<double, double>(descsVec[i].featuresNum, descsVec[i].features, 0, 1);
//			delete[] descsVec[i].features;
//			descsVec[i].features = temp;
//		}
//	}
//	else {
//		cout << "descriptorNormalizationType not set. Descriptors not normalize";
//	}
//
//
//	return descsVec;
//}
//
//template<typename T>
//double * ComputerVision::GaussRaw(int rows, int cols, T * data, double sigma, int k, int interpolateType) {
//	
//	double *temp = NULL, *result = NULL;
//
//	Core gaussX;
//	gaussX.xk = k;
//	gaussX.yk = 0;
//	gaussX.data = new double[2 * k + 1];
//
//	Core gaussY;
//	gaussY.xk = 0;
//	gaussY.yk = k;
//	gaussY.data = new double[2 * k + 1];
//
//	double alpha = 1.0 / (sqrt(2 * PI())*sigma);
//	//double val;
//	
//	for (int i = -k; i <= k; i++)
//		gaussX.data[k + i] = alpha * std::exp(-(i*i) / (2.0 * sigma * sigma));
//
//	for (int i = -k; i <= k; i++)
//		gaussY.data[k + i] = alpha * std::exp(-(i*i) / (2.0 * sigma * sigma));
//
//	temp = ApplyFilterRaw(rows, cols, data, gaussX, interpolateType);
//	result = ApplyFilterRaw(rows, cols, temp, gaussY, interpolateType);
//
//
//	delete[] temp;
//	return result;
//}
//
//
//template<typename T>
//std::vector<Dot> ComputerVision::HarrisRaw(int rows, int cols, T * data, int wk, int localMinK, double Threshold, int pointsNeeded, int PartDerivativeType) {
//
//	int size = rows * cols;
//	int pos;
//
//	//part derivative X, Y
//	double *partDerX = NULL, *partDerY = NULL;
//
//	if (PartDerivativeType == kPartDerivativeSobel) {
//
//		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXSobel, kInterpolateBorder);
//		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYSobel, kInterpolateBorder);
//	}
//	else if (PartDerivativeType == kPartDerivativePrewitt) {
//
//		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXPrewitt, kInterpolateBorder);
//		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYPrewitt, kInterpolateBorder);
//	}
//	else if (PartDerivativeType == kPartDerivativeScharr) {
//
//		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXScharr, kInterpolateBorder);
//		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYScharr, kInterpolateBorder);
//	}
//
//	double* A = new double[size];
//	double* B = new double[size];
//	double* C = new double[size];
//
//	//Core gaussCore = CreateGaussCore(wk / 3, wk); //gauss core for gauss weight
//	double Ix, Iy;
//	int gaussSizeD = 2 * wk + 1;//gause size of dimension (gaussSizeD x gaussSizeD)-square matrix
//	double* gaussWeight = GetGaussWeight(wk / 3.0, gaussSizeD);
//	
//	int gpos, iB, jB;	//B-Biased
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//
//			pos = i * cols + j;
//
//			A[pos] = 0;
//			B[pos] = 0;
//			C[pos] = 0;
//
//			for (int u = -wk; u <= wk; u++) {
//				for (int v = -wk; v <= wk; v++) {
//
//					iB = i + u;
//					jB = j + v;
//
//					if (iB >= 0 && iB < rows - 1 && jB >= 0 && jB < cols - 1) {
//
//						Ix = partDerX[iB*cols + jB];
//						Iy = partDerY[iB*cols + jB];
//					}
//					else {
//
//						Ix = GetVirtualPixel<double>(iB, jB, rows, cols, partDerX, kInterpolateZero);
//						Iy = GetVirtualPixel<double>(iB, jB, rows, cols, partDerY, kInterpolateZero);
//					}
//
//					gpos = (wk + u)*gaussSizeD + (wk + v);
//
//					A[pos] += Ix * Ix * gaussWeight[gpos];
//					B[pos] += Ix * Iy * gaussWeight[gpos];
//					C[pos] += Iy * Iy * gaussWeight[gpos];
//				}
//			}
//		}
//	}
//
//	//calculate response map
//	double* responseMap = HarrisResponse(rows, cols, A, B, C, ComputerVision::kHarrisResponseForstner);
//
//	//show response map
//	//Image test(rows, cols, responseMap, true);
//	//test.GetMaxValue();
//	//cv::imshow("test"+std::to_string(rand()), test.GetMat());
//
//	std::vector<Dot> result;
//	Dot point;
//	bool isPoint;
//	int r;
//
//	//check is some point greater then threshold AND local max
//	double localValue;
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//
//			isPoint = true;
//			pos = i * cols + j;
//
//			if (responseMap[pos] > Threshold) {
//
//				for (int u = -localMinK; u <= localMinK && isPoint; u++) {
//					for (int v = -localMinK; v <= localMinK && isPoint; v++) {
//
//						if ((localMinK >= sqrt(pow(u, 2) + pow(v, 2))) && !(u == 0 && v == 0) && responseMap[pos] <= GetVirtualPixel<T>(i + u, j + v, rows, cols, responseMap, kInterpolateZero))
//							isPoint = false;
//					}
//				}
//			}
//			else {
//				isPoint = false;
//			}
//
//			if (isPoint) {
//
//				point.x = j;
//				point.y = i;
//
//				result.push_back(point);
//			}
//		}
//	}
//	
//	//clear arrays
//	delete[] A;
//	delete[] B;
//	delete[] C;
//
//	std::vector<Dot> resultANMS;
//
//	//calculate ANMS or not
//	if (pointsNeeded > 0) {
//	
//		resultANMS = ANMS(result, rows, cols, responseMap, pointsNeeded, 1);
//		result.clear();
//		result = resultANMS;
//		resultANMS.clear();
//	}
//
//	delete[] responseMap;
//
//	return result;
//}
//
//template<typename srcT, typename dstT>
//dstT * ComputerVision::LinearNormalization(int dataSize, srcT * data, dstT newMin, dstT newMax) {
//
//	if (dataSize != 0 && data != NULL) {
//
//		double min = data[0], max = data[0], k;
//
//		for (int i = 0; i < dataSize; i++) {
//			if (data[i] > max) {
//				max = data[i];
//			}
//
//			if (data[i] < min) {
//				min = data[i];
//			}
//		}
//
//		k = (newMax - newMin) / (max - min);
//
//		dstT* result = new dstT[dataSize];
//
//		for (int i = 0; i < dataSize; i++) {
//
//			result[i] = (dstT)((data[i] - min)*k + newMin);
//		}
//
//		return result;
//	}
//	else {
//		return NULL;
//	}
//}
//
//template <typename T>
//double * ComputerVision::SobelRaw(int rows, int cols, T* data, int PartDerivativeType, int interpolateType) {
//
//	double *temp1 = NULL, *temp2 = NULL, *result = NULL;
//
//	if (PartDerivativeType == kPartDerivativeSobel) {
//
//		temp1 = PartDerivative(rows, cols, data, kPartDerivativeXSobel, interpolateType);
//		temp2 = PartDerivative(rows, cols, data, kPartDerivativeYSobel, interpolateType);
//	}
//	else if (PartDerivativeType == kPartDerivativePrewitt) {
//
//		temp1 = PartDerivative(rows, cols, data, kPartDerivativeXPrewitt, interpolateType);
//		temp2 = PartDerivative(rows, cols, data, kPartDerivativeYPrewitt, interpolateType);
//	}
//	else if (PartDerivativeType == kPartDerivativeScharr) {
//
//		temp1 = PartDerivative(rows, cols, data, kPartDerivativeXScharr, interpolateType);
//		temp2 = PartDerivative(rows, cols, data, kPartDerivativeYScharr, interpolateType);
//	}
//
//	//Image test1(rows, cols, temp1, true);
//	//imshow("Gx image", test1.GetMat());
//
//	//Image test2(rows, cols, temp2, true);
//	//imshow("Gy image", test2.GetMat());
//
//	int size = rows * cols;
//	result = new double[size];
//
//	for (int i = 0; i < size; i++) {
//
//		result[i] = sqrt(temp1[i] * temp1[i] + temp2[i] * temp2[i]);
//	}
//
//	delete[] temp1;
//	delete[] temp2;
//	return result;
//}
//
//template<typename T>
//std::vector<Dot> ComputerVision::MoravecRaw(int rows, int cols, T * data, int wk, int d, int p, double Threshold) {
//	
//	int size = rows*cols;
//	int pos;
//	
//	double** shifts = new double*[8];
//
//	//смещения, начиная с левого и против часовой
//	for (int i = 0; i < 8; i++) {
//
//		shifts[i] = new double[size];
//	}
//
//	//сначала считаем все смещения для каждого пикселя (их всего 8)
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//
//			pos = i * cols + j;
//
//			shifts[0][pos] = pow(data[pos] - GetVirtualPixel<T>(i, j + d, rows, cols, data, kInterpolateReflection), 2);
//			shifts[1][pos] = pow(data[pos] - GetVirtualPixel<T>(i - d, j + d, rows, cols, data, kInterpolateReflection), 2);
//			shifts[2][pos] = pow(data[pos] - GetVirtualPixel<T>(i - d, j, rows, cols, data, kInterpolateReflection), 2);
//			shifts[3][pos] = pow(data[pos] - GetVirtualPixel<T>(i - d, j - d, rows, cols, data, kInterpolateReflection), 2);
//			shifts[4][pos] = pow(data[pos] - GetVirtualPixel<T>(i, j - d, rows, cols, data, kInterpolateReflection), 2);
//			shifts[5][pos] = pow(data[pos] - GetVirtualPixel<T>(i + d, j - d, rows, cols, data, kInterpolateReflection), 2);
//			shifts[6][pos] = pow(data[pos] - GetVirtualPixel<T>(i + d, j, rows, cols, data, kInterpolateReflection), 2);
//			shifts[7][pos] = pow(data[pos] - GetVirtualPixel<T>(i + d, j + d, rows, cols, data, kInterpolateReflection), 2);
//		}
//	}
//
//	double temp;
//	double* s = new double[size];
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//
//			pos = i * cols + j;
//
//			for (int shift = 0; shift < 8; shift++) {
//
//				temp = 0;
//
//				for (int u = -wk; u <= wk; u++) {
//					for (int v = -wk; v <= wk; v++) {
//
//						if (u != 0 && v != 0)
//							temp += GetVirtualPixel<T>(i + u, j + v, rows, cols, shifts[shift], kInterpolateBorder);
//					}
//				}
//
//				if (shift == 0)
//					s[pos] = temp;
//				else
//					s[pos] = std::min(s[pos], temp);
//			}
//		}
//	}
//
//	//Image test(rows, cols, s, true);
//	//cv::imshow("test", test.GetMat());
//
//	std::vector<Dot> result;
//	Dot point;
//	bool isPoint;
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//
//			isPoint = true;
//			pos = i * cols + j;
//
//			if (s[pos] > Threshold) {
//
//				for (int u = -p; u <= p; u++) {
//					for (int v = -p; v <= p; v++) {
//
//						if (!(u == 0 && v == 0) && s[pos] <= GetVirtualPixel<T>(i + u, j + v, rows, cols, s, kInterpolateZero))
//							isPoint = false;
//					}
//				}
//			}
//			else {
//				isPoint = false;
//			}
//
//			if (isPoint) {
//
//				point.x = j;
//				point.y = i;
//
//				result.push_back(point);
//			}
//		}
//	}
//
//	for (int shift = 0; shift < 8; shift++)
//		delete[] shifts[shift];
//
//	delete[] shifts;
//	delete[] s;
//
//	return result;
//}
//
//template<typename T>
//double * ComputerVision::PartDerivative(int rows, int cols, T * data, int PartDerivativeType, int interpolateType) {
//
//	Core Gp1, Gp2;// p(x or y) gradient's filters (separable)
//
//	Gp1.xk = 1;
//	Gp1.yk = 0;
//	Gp1.data = new double[3];
//
//	Gp2.xk = 0;
//	Gp2.yk = 1;
//	Gp2.data = new double[3];
//
//
//	switch (PartDerivativeType) {
//	case kPartDerivativeXSobel:
//
//		Gp1.data[0] = 1; Gp1.data[1] = 0; Gp1.data[2] = -1;
//		Gp2.data[0] = 1; Gp2.data[1] = 2; Gp2.data[2] = 1;
//		break;
//
//	case kPartDerivativeYSobel:
//
//		Gp1.data[0] = 1; Gp1.data[1] = 2; Gp1.data[2] = 1;
//		Gp2.data[0] = 1; Gp2.data[1] = 0; Gp2.data[2] = -1;
//		break;
//
//	case kPartDerivativeXPrewitt:
//
//		Gp1.data[0] = 1; Gp1.data[1] = 0; Gp1.data[2] = -1;
//		Gp2.data[0] = 1; Gp2.data[1] = 1; Gp2.data[2] = 1;
//		break;
//
//	case kPartDerivativeYPrewitt:
//
//		Gp1.data[0] = 1; Gp1.data[1] = 1; Gp1.data[2] = 1;
//		Gp2.data[0] = 1; Gp2.data[1] = 0; Gp2.data[2] = -1;
//		break;
//
//	case kPartDerivativeXScharr:
//
//		Gp1.data[0] = sqrt(3); Gp1.data[1] = 0; Gp1.data[2] = -sqrt(3);
//		Gp2.data[0] = sqrt(3); Gp2.data[1] = 10.0 / (3 * sqrt(3)); Gp2.data[2] = sqrt(3);
//		break;
//
//	case kPartDerivativeYScharr:
//
//		Gp1.data[0] = sqrt(3); Gp1.data[1] = 10.0 / (3 * sqrt(3)); Gp1.data[2] = sqrt(3);
//		Gp2.data[0] = sqrt(3); Gp2.data[1] = 0; Gp2.data[2] = -sqrt(3);
//		break;
//
//	default:
//		break;
//	}
//	double *temp, *result;
//
//	temp = ApplyFilterRaw(rows, cols, data, Gp1, interpolateType);
//	result = ApplyFilterRaw(rows, cols, temp, Gp2, interpolateType);
//
//	delete[] temp;
//
//	return result;
//}
//
//template<typename T>
//T* ComputerVision::ExpandImgByZero(int rows, int cols, T* data, int addX, int addY) {
//
//	int trows = rows + 2 * addX;
//	int tcols = cols + 2 * addY;
//	T* temp = new T[trows*tcols];
//
//	//center
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//
//			temp[(i + addY)*tcols + (j + addX)] = data[i*cols + j];
//		}
//	}
//
//	//top&bottom + angles
//	for (int i = 0; i < addX; i++) {
//		for (int j = 0; j < tcols; j++) {
//
//			temp[i*tcols + j] = 0;
//			temp[(trows - 1 - i)*tcols + j] = 0;
//		}
//	}
//
//	//left&right
//	for (int i = addX; i < trows - addX; i++) {
//		for (int j = 0; j < addY; j++) {
//
//			temp[i*tcols + j] = 0;
//			temp[i*tcols + (tcols - 1 - j)] = 0;
//		}
//	}
//
//	return temp;
//}
