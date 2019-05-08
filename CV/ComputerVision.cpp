#include "stdafx.h"
#include "ComputerVision.hpp"

using namespace std;

std::vector<Dot> ComputerVision::ANMS(std::vector<Dot> points, int rows, int cols, double * responseMap, int pointsNeeded, double c) {

	int size = rows * cols;
	int pointsNum = points.size();

	if (points.size() <= pointsNeeded) {

		return points;
	}

	int* radius = new int[points.size()], radiusTemp;
	int row1, col1, row2, col2, pos1;

	for (int i = 0; i < pointsNum; i++) {

		radius[i] = size;
		row1 = points[i].y;
		col1 = points[i].x;

		for (int j = 0; j < pointsNum; j++) {

			if (i != j) {

				row2 = points[j].y;
				col2 = points[j].x;

				radiusTemp = std::sqrt(std::pow((row1 - row2), 2) + std::pow((col1 - col2), 2));

				if (responseMap[row1*cols + col1] < c*responseMap[row2*cols + col2] && radius[i] > radiusTemp) {

					radius[i] = radiusTemp;
				}
			}
		}
	}

	int* tempR = new int[pointsNum];
	copy(radius, &(radius[pointsNum]), stdext::checked_array_iterator<int*>(tempR, pointsNum));
	int tempVal;


	//for (int i = 0; i < points.size(); i++)
	//	cout<<tempR[i]<<endl;

	for (int i = 0; i < pointsNum; i++) {
		for (int j = 0; j < pointsNum - 1 - i; j++) {

			if (tempR[j] < tempR[j + 1]) {
				tempVal = tempR[j];
				tempR[j] = tempR[j + 1];
				tempR[j + 1] = tempVal;
			}
		}
	}

	//cout  << "########";
	//for (int i = 0; i < points.size(); i++)
	//	cout << tempR[i] << endl;

	int radiusNeeded = tempR[pointsNeeded];
	if (pointsNeeded < pointsNum)
		if(radiusNeeded == tempR[(pointsNeeded -1) + 1])
			++radiusNeeded;

	delete[] tempR;

	vector<Dot> result;
	for (int i = 0; i < pointsNum; i++) {
		if (radius[i] >= radiusNeeded)
			result.push_back(points[i]);
	}

	delete[] radius;
	
	return result;
	/*
	double* values = new double[points.size()];
	for (int i = 0; i < points.size(); i++)
		values[i] = responseMap[points[i].y*cols + points[i].x];

	double temp;
	Dot point;
	for (int i = 0; i < points.size()-1; i++) {
		for (int j = 0; j < points.size()-i-1; j++) {
			if (values[j] < values[j + 1]) {
				
				temp = values[j];
				values[j] = values[j + 1];
				values[j + 1] = temp;

				point = points[j];
				points[j] = points[j + 1];
				points[j + 1] = point;
			}
		}
	}

	//bool* notRemoved = new bool[points.size()];
	//for (int i = 0; i < points.size(); i++)
	//	notRemoved[i] = true;

	int remain = points.size();
	double radius, radiusSelected = 1.0, radiusMin;
 	while (remain > pointsNeeded) {

		radiusMin = cols * rows;

		for (int i = 0; i < points.size()-1; i++) {

			if (values[i] > 0) {
				for (int j = 1; j < points.size(); j++) {

					if (values[j] > 0) {
						if (values[i] * c > values[j]) {

							radius = std::sqrt(std::pow(points[i].x - points[j].x, 2) + std::pow(points[i].y - points[j].y, 2));
							
							if (radius <= radiusSelected) {
								values[j] = -1;
								--remain;
							}
							else {
								if (radius == 4.0)
									int wtf=0;

								radiusMin = std::min(radiusMin, ceil(radius));
							}
						}

					}
				}
			}
		}

		radiusSelected = radiusMin;
	}

	vector<Dot> pointsResult(remain);
	for (int i = 0; i < points.size(); i++)
		if (values[i] > 0)
			pointsResult.push_back(points[i]);

	delete[] values;
	
	return pointsResult;
	*/
}

Image ComputerVision::ApplyFilter(Image & img, Core core, int type) {

	uchar* data = img.GetData();

	double* temp = ApplyFilterRaw(img.GetRowsNumber(), img.GetColsNumber(), data, core, type);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp);

	delete[] temp;
	delete[] data;

	return result;
}

Descriptor ComputerVision::CreateSimpleDescriptorRaw(Dot point, int rows, int cols, double* partDerX, double* partDerY, Dot pointLT, Dot pointRB, int histogramNumX, int histogramNumY, int intervalsNum) {

	int rowsNew = (pointRB.y - pointLT.y + 1);
	int colsNew = (pointRB.x - pointLT.x + 1);
	int size = rowsNew * colsNew;
	//int histogramsNum;
	int pos, posNew;

	Descriptor desc;
	desc.point = point;
	desc.featuresNum = histogramNumX * histogramNumY*intervalsNum;
	desc.features = new double[desc.featuresNum];

	//calculate gradient value and angle [0;2PI]
	double* gradients = new double[size];
	double* angles = new double[size];

	double valueX, valueY;
	double angle;
	for (int i = pointLT.y; i <= pointRB.y; i++) {
		for (int j = pointLT.x; j <= pointRB.x; j++) {

			//pos = i * cols + j;
			posNew = (i - pointLT.y)*colsNew + (j - pointLT.x);

			valueX = GetVirtualPixel(i, j, rows, cols, partDerX);
			valueY = GetVirtualPixel(i, j, rows, cols, partDerY);

			gradients[posNew] = sqrt(std::pow(valueX, 2) + std::pow(valueY, 2));
			angle =atan2(valueY, valueX);
			angles[posNew] = angle >= 0 ? angle : 2*PI() + angle;
				//cout << angles[posNew] << endl;
		}
	}

	//calculate gauss weight for this part of image with size = sizeXY x SizeXY
	double* gaussWeight = GetGaussWeight(double(rowsNew) / 6.0, rowsNew);//!!! is sigma correct?

	//get all pixels, that belong to specifecal histogram
	vector<double> pixelValues;
	vector<double> pixelAngles;

	//grid division by histograms number
	double histogramStepX = double(colsNew) / histogramNumX;
	double histogramStepY = double(rowsNew) / histogramNumY;
	
	//double pixelValue;
	double pixelAffil;
	int iStart, jStart;
	double iMax, jMax;
	
	//for each hist by OY
	for (int hy = 0; hy < histogramNumY; hy++) {

		iStart = int(hy * histogramStepY);
		iMax = (hy + 1) * histogramStepY;

		//for each hist by OX
		for (int hx = 0; hx < histogramNumX; hx++) {

			jStart = int(hx * histogramStepX);
			jMax = (hx + 1) * histogramStepX;

			pixelValues.clear();
			pixelAngles.clear();

			//check all affiled pixels
			for (int i = iStart; i < iMax; i++) {
				for (int j = jStart; j < jMax; j++) {

					posNew = i * colsNew + j;
					pixelAffil = 1;

					//check edged value and calculate % of affil
					if (i == iStart)
						pixelAffil *= 1 - (hy * histogramStepY - i);
					if (i+1 > iMax)
						pixelAffil *= (hy + 1) * histogramStepY - i;
					if (j == jStart)
						pixelAffil *= 1 - (hx * histogramStepX - j);
					if (j + 1 > jMax)
						pixelAffil *= (hx + 1) * histogramStepX - j;

					//calculate values for pixel as gradient val * gauss weight * pixel affil (mostly for border's values)
					pixelValues.push_back(gradients[posNew] * gaussWeight[posNew] * pixelAffil);
					pixelAngles.push_back(angles[posNew]);
				}
			}

			Histogram hist = CreateHistogramRaw(pixelValues, pixelAngles, intervalsNum);
			
			for (int  i = 0; i < intervalsNum; i++) {
				desc.features[(hy*histogramNumX + hx)*intervalsNum + i] = hist.intervals[i].data;
			}
			//desc.features[(hy*histogramNumX + hx)*intervalsNum] 
		}
	}

	delete[] gradients;
	delete[] angles;
	delete[] gaussWeight;

	return desc;
}

Histogram ComputerVision::CreateHistogramRaw(std::vector<double> pixelValues, std::vector<double> pixelAngles, int intervalsNum) {
	
	Histogram hist;
	hist.intervalsNum = intervalsNum;
	hist.intervals = new Interval[intervalsNum];

	double step = (2*PI()) / intervalsNum;	//because val in[0; 2*PI]
	double startValue = step / 2;		//because we use mid value (half step)

	for (int i = 0; i < intervalsNum; i++) {

		hist.intervals[i].data = 0;
		hist.intervals[i].midVal = startValue + i * step;
	}

	//double Th, G;
	double k;
	int pos, column0, column1;

	for (int i = 0; i < pixelValues.size(); i++) {
		//edged value
		if (pixelAngles[i] >= hist.intervals[intervalsNum - 1].midVal || pixelAngles[i] < hist.intervals[0].midVal) {

			column0 = intervalsNum - 1;
			column1 = 0;

			if (pixelAngles[i] < PI()) {

				k = (hist.intervals[0].midVal-pixelAngles[i]) / step;							//!!!
			}
			else {

				k = 1 - (pixelAngles[i] - hist.intervals[intervalsNum - 1].midVal) / step;		//!!!
			}
		}
		else {
			//find nearest column (basket)
			for (int b = 0; b < intervalsNum - 1; b++) {

				if (pixelAngles[i] >= hist.intervals[b].midVal && pixelAngles[i] < hist.intervals[(b + 1)].midVal) {

					column0 = b;
					column1 = (b + 1);
					break;
				}
			}

			k = 1 - (pixelAngles[i] - hist.intervals[column0].midVal) / (hist.intervals[column1].midVal - hist.intervals[column0].midVal);	//!!!
		}

		hist.intervals[column0].data += pixelValues[i] * k;
		hist.intervals[column1].data += pixelValues[i] * (1 - k);
		//if (hist.intervals[column0].data < 0 || hist.intervals[column1].data < 0)
		//	int wtf = 0;
	}

	return hist;
}

Image ComputerVision::Canny(Image & img, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType) {

	double* data = img.GetNormalizeDataF();
	double* temp = CannyRaw(img.GetRowsNumber(), img.GetColsNumber(), data, sigma, k, lowThreshhold, highThreshhold, interpolateType);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);

	delete[] temp;
	delete[] data;

	return result;
}

std::vector<Descriptor> ComputerVision::CreateDescriptors(Image & img, std::vector<Dot> points, int windowSizeX, int windowSizeY, int histogramNumX, int histogramNumY, int intervalsNum, int descriptorType, int descriptorNormalizationType, int PartDerivativeType) {
	
	double* data = img.GetNormalizeDataF(); 
	
	vector<Descriptor> descriptors = CreateDescriptorsRaw(img.GetRowsNumber(), img.GetColsNumber(), data, points, windowSizeX, windowSizeY, histogramNumX, histogramNumY, intervalsNum, descriptorType, descriptorNormalizationType, PartDerivativeType);
	
	delete[] data;

	return descriptors;
}

double ComputerVision::DescriptorsDifference(Descriptor desc0, Descriptor desc1, int descriptorsComparisonType) {
	
	double difference = 0;

	if (descriptorsComparisonType == kDescriptorsComparisonEuclid) {

		for (int i = 0; i < desc0.featuresNum; i++) {
				difference += std::pow(desc0.features[i] - desc1.features[i], 2);
		}

		difference = std::sqrt(difference);
	}
	else if (descriptorsComparisonType == kDescriptorsComparisonManhattan) {

		for (int i = 0; i < desc0.featuresNum; i++) {
			difference += std::abs(desc0.features[i] - desc1.features[i]);
		}
	}
	else if (descriptorsComparisonType == kDescriptorsComparisonSSD) {

		for (int i = 0; i < desc0.featuresNum; i++) {
			difference += std::pow(desc0.features[i] - desc1.features[i], 2);
		}
	}

	return difference;
}

int * ComputerVision::DescriptorsMatchingRaw(std::vector<Descriptor> descs1, std::vector<Descriptor> descs2, int descriptorsComparisionType, int descriptorsMatchingType, double thresh) {
	
	int descs1Num = descs1.size(), descs2Num = descs2.size();
	double* descsDiff = new double[descs1Num * descs2Num];

	for (int i = 0; i < descs1Num; i++) {
		for (int j = 0; j < descs2Num; j++) {
			descsDiff[i*descs2Num + j] = DescriptorsDifference(descs1[i], descs2[j], descriptorsComparisionType);
		}
	}

	int* desc1PairInDesc2 = new int[descs1Num];
	int* desc2PairInDesc1;
	int* desc1PairInDesc2NND;

	switch (descriptorsMatchingType) {
	case kDescriptorsMatchingNNDR:

		desc1PairInDesc2NND = new int[descs1Num];

		for (int i = 0; i < descs1Num; i++) {

			desc1PairInDesc2[i] = 0;
			desc1PairInDesc2NND[i] = -1;

			//find nearest descriptor and second nearest
			for (int j = 1; j < descs2Num; j++) {
				if (descsDiff[i*descs2Num + desc1PairInDesc2[i]] > descsDiff[i*descs2Num + j]) {
					
					desc1PairInDesc2NND[i] = desc1PairInDesc2[i];
					desc1PairInDesc2[i] = j;
				}
			}

			if (desc1PairInDesc2NND[i] != -1) {
				if (descsDiff[i*descs2Num + desc1PairInDesc2[i]] / descsDiff[i*descs2Num + desc1PairInDesc2NND[i]] > thresh)
					desc1PairInDesc2[i] = -1;
			}
		}

		delete[] desc1PairInDesc2NND;

		break;

	case kDescriptorsMatchingMutal:

		for (int i = 0; i < descs1Num; i++) {

			desc1PairInDesc2[i] = 0;

			for (int j = 1; j < descs2Num; j++) {

				if (descsDiff[i*descs2Num + desc1PairInDesc2[i]] > descsDiff[i*descs2Num + j])
					desc1PairInDesc2[i] = j;
			}
		}

		desc2PairInDesc1 = new int[descs2Num];

		for (int i = 0; i < descs2Num; i++) {

			desc2PairInDesc1[i] = 0;

			for (int j = 1; j < descs1Num; j++) {

				if (descsDiff[desc2PairInDesc1[i] *descs2Num + i] > descsDiff[j*descs2Num + i])
					desc2PairInDesc1[i] = j;
			}
		}

		for (int i = 0; i < descs1Num; i++) {
			if (i != desc2PairInDesc1[desc1PairInDesc2[i]])
				desc1PairInDesc2[i] = -1;
		}

		delete[] desc2PairInDesc1;

		break;

	case kDescriptorsMatchingBase:
	default:

		for (int i = 0; i < descs1Num; i++) {
			
			desc1PairInDesc2[i] = 0;

			for (int j = 1; j < descs2Num; j++) {
				if (descsDiff[i*descs2Num + desc1PairInDesc2[i]] > descsDiff[i*descs2Num + j])
					desc1PairInDesc2[i] = j;
			}
		}

		break;
	}

	return desc1PairInDesc2;
}

std::vector<PairDot> ComputerVision::DescriptorsMatching(std::vector<Descriptor> descs1, std::vector<Descriptor> descs2, int descriptorsComparisionType, int descriptorsMatchingType, double thresh) {

	int* desc1PairInDesc2 = DescriptorsMatchingRaw(descs1, descs2, descriptorsComparisionType, descriptorsMatchingType, thresh);

	std::vector<PairDot> matchedPoints(descs1.size());
	PairDot pairPoints;
	
	for (int i = 0; i < descs1.size(); i++) {
		if (desc1PairInDesc2[i] != -1) {

			pairPoints.point0 = descs1[i].point;
			pairPoints.point1 = descs2[desc1PairInDesc2[i]].point;
			matchedPoints.push_back(pairPoints);
		}
	}

	return matchedPoints;
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

Image ComputerVision::CutOutImage(Image img, int row0, int col0, int row1, int col1) {

	uchar* data = img.GetData();
	uchar* temp = CutOutImageRaw(img.GetRowsNumber(), img.GetColsNumber(), data, row0, col0, row1, col1);

	Image result(row1 - row0 + 1, col1 - col0 + 1, temp, true);

	delete[] temp;
	delete[] data;

	return result;
}

Image ComputerVision::GaussDefault(Image & img, double sigma, int interpolateType) {
	return Gauss(img, sigma, 3*sigma, interpolateType);
}

double * ComputerVision::GetGaussWeight(double sigma, int size) {

	double* data = new double[size*size];

	double alpha = 1.0 / (2 * PI()*pow(sigma, 2));
	double center = double(size - 1) / 2.0;
	//int pos;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {

			data[i* size + j] = alpha * std::exp(-(pow(i - center, 2) + pow(j - center, 2)) / (2.0 * sigma * sigma));
		}
	}

	return data;
}

Image ComputerVision::Gauss(Image & img, double sigma, int k, int interpolateType) {
	
	double* data = img.GetNormalizeDataF();
	double* temp = GaussRaw(img.GetRowsNumber(), img.GetColsNumber(), data, sigma, k, interpolateType);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);

	delete[] temp;
	delete[] data;

	return result;
}

std::vector<Dot> ComputerVision::Harris(Image & img, int wk, int localMinK, double Threshold, int pointsNeeded, int PartDerivativeType) {

	double* data = img.GetNormalizeDataF();
	vector<Dot> result = HarrisRaw(img.GetRowsNumber(), img.GetColsNumber(), data, wk, localMinK, Threshold, pointsNeeded, PartDerivativeType);

	delete[] data;

	return result;
}

double * ComputerVision::HarrisResponse(int rows, int cols, double * A, double * B, double * C, int harrisResponseType) {
	
	int size = rows * cols;
	double* result = new double[size];
	int k = 0.05; // for base

	switch (harrisResponseType) {
	//in work!
	/*
	case kHarrisResponseDirect:

		delete[] result;
		result = NULL;
		break;
	*/

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
