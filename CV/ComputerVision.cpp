#include "stdafx.h"
//
//using namespace std;
//
//std::vector<Dot> ComputerVision::ANMS(std::vector<Dot> points, int rows, int cols, double * responseMap, int pointsNeeded, double c) {
//
//	int size = rows * cols;
//	int pointsNum = points.size();
//
//	if (points.size() <= pointsNeeded) {
//
//		return points;
//	}
//
//	double* radius = new double[points.size()], radiusTemp;
//	int row1, col1, row2, col2;// , pos1;
//
//	for (int i = 0; i < pointsNum; i++) {
//
//		radius[i] = size;
//		row1 = points[i].y;
//		col1 = points[i].x;
//
//		for (int j = 0; j < pointsNum; j++) {
//
//			if (i != j) {
//
//				row2 = points[j].y;
//				col2 = points[j].x;
//
//				radiusTemp = std::sqrt(std::pow((row1 - row2), 2) + std::pow((col1 - col2), 2));
//
//				if (responseMap[row1*cols + col1] < c*responseMap[row2*cols + col2] && radius[i] > radiusTemp) {
//
//					radius[i] = radiusTemp;
//				}
//			}
//		}
//	}
//
//	double* tempR = new double[pointsNum];
//	copy(radius, &(radius[pointsNum]), stdext::checked_array_iterator<double*>(tempR, pointsNum));
//	
//
//
//	//for (int i = 0; i < points.size(); i++)
//	//	cout<<tempR[i]<<endl;
//	double tempVal;
//	for (int i = 0; i < pointsNum; i++) {
//		for (int j = 0; j < pointsNum - 1 - i; j++) {
//
//			if (tempR[j] < tempR[j + 1]) {
//				tempVal = tempR[j];
//				tempR[j] = tempR[j + 1];
//				tempR[j + 1] = tempVal;
//			}
//		}
//	}
//
//	//cout  << "########";
//	//for (int i = 0; i < points.size(); i++)
//	//	cout << tempR[i] << endl;
//
//	double radiusNeeded = tempR[pointsNeeded];
//	if (pointsNeeded < pointsNum)
//		if(radiusNeeded == tempR[(pointsNeeded -1) + 1])
//			radiusNeeded+=0.1;
//
//	delete[] tempR;
//
//	vector<Dot> result;
//	for (int i = 0; i < pointsNum; i++) {
//		if (radius[i] >= radiusNeeded)
//			result.push_back(points[i]);
//	}
//
//	delete[] radius;
//	
//	return result;
//	/*
//	double* values = new double[points.size()];
//	for (int i = 0; i < points.size(); i++)
//		values[i] = responseMap[points[i].y*cols + points[i].x];
//
//	double temp;
//	Dot point;
//	for (int i = 0; i < points.size()-1; i++) {
//		for (int j = 0; j < points.size()-i-1; j++) {
//			if (values[j] < values[j + 1]) {
//				
//				temp = values[j];
//				values[j] = values[j + 1];
//				values[j + 1] = temp;
//
//				point = points[j];
//				points[j] = points[j + 1];
//				points[j + 1] = point;
//			}
//		}
//	}
//
//	//bool* notRemoved = new bool[points.size()];
//	//for (int i = 0; i < points.size(); i++)
//	//	notRemoved[i] = true;
//
//	int remain = points.size();
//	double radius, radiusSelected = 1.0, radiusMin;
// 	while (remain > pointsNeeded) {
//
//		radiusMin = cols * rows;
//
//		for (int i = 0; i < points.size()-1; i++) {
//
//			if (values[i] > 0) {
//				for (int j = 1; j < points.size(); j++) {
//
//					if (values[j] > 0) {
//						if (values[i] * c > values[j]) {
//
//							radius = std::sqrt(std::pow(points[i].x - points[j].x, 2) + std::pow(points[i].y - points[j].y, 2));
//							
//							if (radius <= radiusSelected) {
//								values[j] = -1;
//								--remain;
//							}
//							else {
//								if (radius == 4.0)
//									int wtf=0;
//
//								radiusMin = std::min(radiusMin, ceil(radius));
//							}
//						}
//
//					}
//				}
//			}
//		}
//
//		radiusSelected = radiusMin;
//	}
//
//	vector<Dot> pointsResult(remain);
//	for (int i = 0; i < points.size(); i++)
//		if (values[i] > 0)
//			pointsResult.push_back(points[i]);
//
//	delete[] values;
//	
//	return pointsResult;
//	*/
//}
//
//Descriptor ComputerVision::CreateSimpleDescriptorRaw(Dot point, int rows, int cols, double* gradients, double* angles, int histogramNumX, int histogramNumY, int intervalsNum, double alpha = 0) {
//
//	int size = rows * cols;
//	int histNum = histogramNumX * histogramNumY;
//	int pos;
//
//	Descriptor desc;
//	desc.point = point;
//	desc.featuresNum = histNum * intervalsNum;
//	desc.features = new double[desc.featuresNum];
//
//	//grid division step by histograms
//	double histStepX = double(cols) / histogramNumX;
//	double histStepY = double(rows) / histogramNumY;
//
//	//arrays of pixels, that belong to specifical histogram
//	vector<double>* pixelValues = new vector<double>[histNum];
//	vector<double>* pixelAngles = new vector<double>[histNum];
//
//	double x, y, xRotated, yRotated;	//(x, y) before and after rotate
//	int histXIndex, histYIndex, histPos;//histogram index by x and y, position in array
//
//	//cos(alpha), sin(alpha)
//	double cosAlpha = std::cos(alpha);
//	double sinAlpha = std::sin(alpha);
//
//	double descriptorMinX = 0, descriptorMaxX = cols;	// descriptor border
//	double  descriptorMinY = 0, descriptorMaxY = rows;	// descriptor border
//	double  descriptorCenterX = (descriptorMaxX - descriptorMinX) / 2;	// descriptor center
//	double  descriptorCenterY = (descriptorMaxY - descriptorMinY) / 2;	// descriptor center
//
//	//how much pixel cross border line
//	double crossLeftBorderK, crossTopBorderK, crossBotBorderK, crossRightBorderK;
//
//	//error for comparing
//	double eps = 0.0001;
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//
//			//pixel position on image
//			pos = i * cols + j;
//
//			//0.5 - pixel's center
//			x = j + 0.5;
//			y = i + 0.5;
//
//			//new pixel center
//			xRotated = descriptorCenterX + (x - descriptorCenterX) * cosAlpha - (y - descriptorCenterY) * sinAlpha;
//			yRotated = descriptorCenterY + (y - descriptorCenterY) * cosAlpha + (x - descriptorCenterX) * sinAlpha;
//
//			//check if we went beyound the boundaries 
//			if (xRotated <= -0.5 + descriptorMinX || xRotated >= descriptorMaxX + 0.5)
//				continue;
//
//			if (yRotated <= -0.5 + descriptorMinY || yRotated >= descriptorMaxY + 0.5)
//				continue;
//
//			//calculate hist index
//			histXIndex = std::floor(xRotated / histStepX);
//			histYIndex = std::floor(yRotated / histStepY);
//			histPos = histYIndex * histogramNumX + histXIndex;
//
//			//calculate how much pixel cross border line
//			crossLeftBorderK = std::max(histXIndex * histStepX - xRotated - 0.5, 0.);
//			crossRightBorderK = std::max(xRotated + 0.5 - histXIndex * histStepX, 0.);
//			crossTopBorderK = std::max(histYIndex * histStepY - yRotated - 0.5, 0.);
//			crossBotBorderK = std::max(yRotated + 0.5 - histYIndex * histStepY, 0.);
//
//			//if pixel has crossed left border AND we have histogram on left 
//			if (crossLeftBorderK > 0 + eps && histXIndex - 1 > 0) {
//
//				//get gradient value for this hist
//				pixelValues[histPos - 1].push_back(gradients[pos] * crossLeftBorderK * (1 - crossTopBorderK)* (1 - crossBotBorderK));
//				pixelAngles[histPos - 1].push_back(angles[pos] - alpha);
//			}
//
//			//if pixel has crossed Right border AND we have histogram on Right 
//			if (crossRightBorderK > 0 + eps && histXIndex + 1 < histogramNumX) {
//
//				//get gradient value for this hist
//				pixelValues[histPos + 1].push_back(gradients[pos] * crossRightBorderK * (1 - crossTopBorderK)* (1 - crossBotBorderK));
//				pixelAngles[histPos + 1].push_back(angles[pos] - alpha);
//			}
//
//			//if pixel has crossed Top border AND we have histogram on Top 
//			if (crossTopBorderK > 0 + eps && histYIndex - 1 > 0) {
//
//				//get gradient value for this hist
//				pixelValues[histPos - histogramNumX].push_back(gradients[pos] * crossLeftBorderK * (1 - crossTopBorderK)* (1 - crossBotBorderK));
//				pixelAngles[histPos - histogramNumX].push_back(angles[pos] - alpha);
//			}
//
//			//if pixel has crossed Right border AND we have histogram on Right 
//			if (crossRightBorderK > 0 + eps && histXIndex + 1 < histogramNumX) {
//
//				//get gradient value for this hist
//				pixelValues[histPos + 1].push_back(gradients[pos] * crossRightBorderK * (1 - crossTopBorderK)* (1 - crossBotBorderK));
//				pixelAngles[histPos + 1].push_back(angles[pos] - alpha);
//			}
//
//			///////////////////
//			////////////////////
//			////////////////////
//			//////////////////////////
//			//////////////////////////
//			//////////////////////////
//		}
//	}
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//			for(int hi = 0; hi<hist)
//			
//		}
//	}
//	
//	//double pixelValue;
//	double pixelAffil;
//	int iStart, jStart;
//	double iMax, jMax;
//	
//	//for each hist by OY
//	for (int hy = 0; hy < histogramNumY; hy++) {
//
//		iStart = int(hy * histogramStepY);
//		iMax = (hy + 1) * histogramStepY;
//
//		//for each hist by OX
//		for (int hx = 0; hx < histogramNumX; hx++) {
//
//			jStart = int(hx * histogramStepX);
//			jMax = (hx + 1) * histogramStepX;
//
//			pixelValues.clear();
//			pixelAngles.clear();
//
//			//check all affiled pixels
//			for (int i = iStart; i < iMax; i++) {
//				for (int j = jStart; j < jMax; j++) {
//
//					pos = i * cols + j;
//					pixelAffil = 1;
//
//					//check edged value and calculate % of affil
//					if (i == iStart)
//						pixelAffil *= 1 - (hy * histogramStepY - i);
//					if (i+1 > iMax)
//						pixelAffil *= (hy + 1) * histogramStepY - i;
//					if (j == jStart)
//						pixelAffil *= 1 - (hx * histogramStepX - j);
//					if (j + 1 > jMax)
//						pixelAffil *= (hx + 1) * histogramStepX - j;
//
//					//calculate values for pixel as gradient val * gauss weight * pixel affil (mostly for border's values)
//					pixelValues.push_back(gradients[pos] * pixelAffil);
//					pixelAngles.push_back(angles[pos]);
//				}
//			}
//
//			Histogram hist = CreateHistogramRaw(pixelValues, pixelAngles, 0, 2*PI(), intervalsNum);
//			
//			for (int  i = 0; i < intervalsNum; i++) {
//				desc.features[(hy*histogramNumX + hx)*intervalsNum + i] = hist.intervals[i].data;
//			}
//			//desc.features[(hy*histogramNumX + hx)*intervalsNum] 
//		}
//	}
//
//	return desc;
//}
//
//std::vector<double> ComputerVision::CalculateDescriptorOrientationAnglesRaw(int rows, int cols, double* gradientsWeighted, double* angles, int intervalsNum, double thresh, int alphaNumMax) {
//	
//	int size = rows * cols;
//	
//	//convert arrays of pixel values and angles to vectors
//	std::vector<double> pixelValues(gradientsWeighted, gradientsWeighted + size);
//	std::vector<double> pixelAngles(angles, angles + size);
//	
//	//create histogram
//	Histogram hist = CreateHistogramRaw(pixelValues, pixelAngles, 0, 2 * PI(), intervalsNum);
//
//	//find max value in interval
//	double maxVal = hist.intervals[0].data;
//
//	for (int i = 1; i < intervalsNum; i++) {
//
//		if (maxVal < hist.intervals[i].data)
//			maxVal = hist.intervals[i].data;
//
//	}
//
//	//find all intervals, than > when maxVal*thresh
//	vector<int> pretendetIntervals;
//
//	for (int i = 0; i < intervalsNum; i++) {
//
//		if (thresh*maxVal < hist.intervals[i].data) {
//
//			pretendetIntervals.push_back(i);
//		}
//	}
//
//	//orientation of decriptor
//	vector<double> alphas;
//
//	int tempInterval;
//	for (int i = 0; i < pretendetIntervals.size()-1; i++) {
//		for (int j = 0; j < pretendetIntervals.size()-i-1; j++) {
//			if (hist.intervals[pretendetIntervals[j]].data < hist.intervals[pretendetIntervals[j + 1]].data) {
//
//				tempInterval = pretendetIntervals[j];
//				pretendetIntervals[j] = pretendetIntervals[j + 1];
//				pretendetIntervals[j + 1] = tempInterval;
//			}
//		}
//	}
//
//	for (int i = 0; i < std::min(alphaNumMax, (int)pretendetIntervals.size()); i++)
//		alphas.push_back(hist.intervals[pretendetIntervals[i]].midVal);
//
//	return alphas;
//}
//
//
//GradientsAnglesArray ComputerVision::CalculateGradientsFullRaw(int rows, int cols, double * partDerX, double * partDerY) {
//	
//	Dot pointLT, pointRB;
//
//	pointLT.x = 0;
//	pointLT.y = 0;
//
//	pointRB.x = rows - 1;
//	pointRB.y = cols - 1;
//
//	return CalculateGradientsFullRaw(rows,cols,partDerX, partDerY, pointLT, pointRB);
//}
//
//GradientsAnglesArray ComputerVision::CalculateGradientsFullRaw(int rows, int cols, double * partDerX, double * partDerY, Dot pointLT, Dot pointRB) {
//	
//	int rowsNew = pointRB.y - pointLT.y+1;
//	int colsNew = pointRB.x - pointLT.x+1;
//
//	int size = rowsNew * colsNew;
//	double* gradients = new double[size];
//	double* angles = new double[size];
//
//	double valueX, valueY, angle;
//	int posNew;
//
//	for (int i = pointLT.y; i <= pointRB.y; i++) {
//		for (int j = pointLT.x; j <= pointRB.x; j++) {
//
//			posNew = (i - pointLT.y)*colsNew + (j - pointLT.x);
//
//			if (0 <= i && i < rows && 0 <= j && j < cols) {
//				valueX = partDerX[i * cols + j];
//				valueY = partDerY[i * cols + j];
//			}
//			else {
//				valueX = GetVirtualPixel(i, j, rows, cols, partDerX);
//				valueY = GetVirtualPixel(i, j, rows, cols, partDerY);
//			}
//
//			gradients[posNew] = sqrt(std::pow(valueX, 2) + std::pow(valueY, 2));
//			angle = atan2(valueY, valueX);
//			angles[posNew] = angle >= 0 ? angle : 2 * PI() + angle;
//		}
//	}
//
//	GradientsAnglesArray gaa;
//	gaa.gradients = gradients;
//	gaa.angles = angles;
//
//	return gaa;
//}
//
//Histogram ComputerVision::CreateHistogramRaw(std::vector<double> values, std::vector<double> positions, double intervalsMinVal, double intervalsMaxVal, int intervalsNum) {
//	
//	Histogram hist;
//	hist.intervalsNum = intervalsNum;
//	hist.intervals = new Interval[intervalsNum];
//
//	double step = (intervalsMaxVal - intervalsMinVal) / intervalsNum;
//	double startValue = step / 2;		//because we use mid value (half step)
//
//	for (int i = 0; i < intervalsNum; i++) {
//
//		hist.intervals[i].data = 0;
//		hist.intervals[i].midVal = startValue + i * step;
//	}
//
//	double k;
//	int  column0, column1;// pos;
//
//	for (int i = 0; i < values.size(); i++) {
//		//edged right value
//		if (positions[i] >= hist.intervals[intervalsNum - 1].midVal) {
//
//			column0 = intervalsNum - 1;
//			column1 = 0;
//
//			k = 1 - (positions[i] - hist.intervals[intervalsNum - 1].midVal) / step;
//
//		}//edged left value
//		else if(positions[i] < hist.intervals[0].midVal) {
//
//			column0 = intervalsNum - 1;
//			column1 = 0;
//				
//			k = (hist.intervals[0].midVal- positions[i]) / step;							//!!!
//			
//		}
//		else {
//			//find nearest column (basket)
//			for (int b = 0; b < intervalsNum - 1; b++) {
//
//				if (positions[i] >= hist.intervals[b].midVal && positions[i] < hist.intervals[(b + 1)].midVal) {
//
//					column0 = b;
//					column1 = (b + 1);
//					break;
//				}
//			}
//
//			k = 1 - (positions[i] - hist.intervals[column0].midVal) / (hist.intervals[column1].midVal - hist.intervals[column0].midVal);	//!!!
//		}
//
//		hist.intervals[column0].data += values[i] * k;
//		hist.intervals[column1].data += values[i] * (1 - k);
//		//if (hist.intervals[column0].data < 0 || hist.intervals[column1].data < 0)
//		//	int wtf = 0;
//	}
//
//	return hist;
//}
//
//Image ComputerVision::Canny(Image & img, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType) {
//
//	double* data = img.GetNormalizeDataF();
//	double* temp = CannyRaw(img.GetRowsNumber(), img.GetColsNumber(), data, sigma, k, lowThreshhold, highThreshhold, interpolateType);
//
//	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);
//
//	delete[] temp;
//	delete[] data;
//
//	return result;
//}
//
//std::vector<Descriptor> ComputerVision::CreateDescriptors(Image & img, std::vector<Dot> points, int windowSizeX, int windowSizeY, int histogramNumX, int histogramNumY, int intervalsNum, int descriptorType, int descriptorNormalizationType, int PartDerivativeType) {
//	
//	double* data = img.GetNormalizeDataF(); 
//	
//	vector<Descriptor> descriptors = CreateDescriptorsRaw(img.GetRowsNumber(), img.GetColsNumber(), data, points, windowSizeX, windowSizeY, histogramNumX, histogramNumY, intervalsNum, descriptorType, descriptorNormalizationType, PartDerivativeType);
//	
//	delete[] data;
//
//	return descriptors;
//}
//
//double ComputerVision::DescriptorsDifference(Descriptor desc0, Descriptor desc1, int descriptorsComparisonType) {
//	
//	double difference = 0;
//
//	if (descriptorsComparisonType == kDescriptorsComparisonEuclid) {
//
//		for (int i = 0; i < desc0.featuresNum; i++) {
//				difference += std::pow(desc0.features[i] - desc1.features[i], 2);
//		}
//
//		difference = std::sqrt(difference);
//	}
//	else if (descriptorsComparisonType == kDescriptorsComparisonManhattan) {
//
//		for (int i = 0; i < desc0.featuresNum; i++) {
//			difference += std::abs(desc0.features[i] - desc1.features[i]);
//		}
//	}
//	else if (descriptorsComparisonType == kDescriptorsComparisonSSD) {
//
//		for (int i = 0; i < desc0.featuresNum; i++) {
//			difference += std::pow(desc0.features[i] - desc1.features[i], 2);
//		}
//	}
//
//	return difference;
//}
//
//int * ComputerVision::DescriptorsMatchingRaw(std::vector<Descriptor> descs1, std::vector<Descriptor> descs2, int descriptorsComparisionType, int descriptorsMatchingType, double thresh) {
//	
//	int descs1Num = descs1.size(), descs2Num = descs2.size();
//	double* descsDiff = new double[descs1Num * descs2Num];
//
//	for (int i = 0; i < descs1Num; i++) {
//		for (int j = 0; j < descs2Num; j++) {
//			descsDiff[i*descs2Num + j] = DescriptorsDifference(descs1[i], descs2[j], descriptorsComparisionType);
//		}
//	}
//
//	int* desc1PairInDesc2 = new int[descs1Num];
//	int* desc2PairInDesc1;
//	int* desc1PairInDesc2NND;
//
//	switch (descriptorsMatchingType) {
//	case kDescriptorsMatchingNNDR:
//
//		desc1PairInDesc2NND = new int[descs1Num];
//
//		for (int i = 0; i < descs1Num; i++) {
//
//			desc1PairInDesc2[i] = 0;
//			desc1PairInDesc2NND[i] = -1;
//
//			//find nearest descriptor and second nearest
//			for (int j = 1; j < descs2Num; j++) {
//				if (descsDiff[i*descs2Num + desc1PairInDesc2[i]] > descsDiff[i*descs2Num + j]) {
//					
//					desc1PairInDesc2NND[i] = desc1PairInDesc2[i];
//					desc1PairInDesc2[i] = j;
//				}
//			}
//
//			if (desc1PairInDesc2NND[i] != -1) {
//				if (descsDiff[i*descs2Num + desc1PairInDesc2[i]] / descsDiff[i*descs2Num + desc1PairInDesc2NND[i]] > thresh)
//					desc1PairInDesc2[i] = -1;
//			}
//		}
//
//		delete[] desc1PairInDesc2NND;
//
//		break;
//
//	case kDescriptorsMatchingMutal:
//
//		for (int i = 0; i < descs1Num; i++) {
//
//			desc1PairInDesc2[i] = 0;
//
//			for (int j = 1; j < descs2Num; j++) {
//
//				if (descsDiff[i*descs2Num + desc1PairInDesc2[i]] > descsDiff[i*descs2Num + j])
//					desc1PairInDesc2[i] = j;
//			}
//		}
//
//		desc2PairInDesc1 = new int[descs2Num];
//
//		for (int i = 0; i < descs2Num; i++) {
//
//			desc2PairInDesc1[i] = 0;
//
//			for (int j = 1; j < descs1Num; j++) {
//
//				if (descsDiff[desc2PairInDesc1[i] *descs2Num + i] > descsDiff[j*descs2Num + i])
//					desc2PairInDesc1[i] = j;
//			}
//		}
//
//		for (int i = 0; i < descs1Num; i++) {
//			if (i != desc2PairInDesc1[desc1PairInDesc2[i]])
//				desc1PairInDesc2[i] = -1;
//		}
//
//		delete[] desc2PairInDesc1;
//
//		break;
//
//	case kDescriptorsMatchingBase:
//	default:
//
//		for (int i = 0; i < descs1Num; i++) {
//			
//			desc1PairInDesc2[i] = 0;
//
//			for (int j = 1; j < descs2Num; j++) {
//				if (descsDiff[i*descs2Num + desc1PairInDesc2[i]] > descsDiff[i*descs2Num + j])
//					desc1PairInDesc2[i] = j;
//			}
//		}
//
//		break;
//	}
//
//	return desc1PairInDesc2;
//}
//
//std::vector<PairDot> ComputerVision::DescriptorsMatching(std::vector<Descriptor> descs1, std::vector<Descriptor> descs2, int descriptorsComparisionType, int descriptorsMatchingType, double thresh) {
//
//	int* desc1PairInDesc2 = DescriptorsMatchingRaw(descs1, descs2, descriptorsComparisionType, descriptorsMatchingType, thresh);
//
//	std::vector<PairDot> matchedPoints(descs1.size());
//	PairDot pairPoints;
//	
//	for (int i = 0; i < descs1.size(); i++) {
//		if (desc1PairInDesc2[i] != -1) {
//
//			pairPoints.point0 = descs1[i].point;
//			pairPoints.point1 = descs2[desc1PairInDesc2[i]].point;
//			matchedPoints.push_back(pairPoints);
//		}
//	}
//
//	return matchedPoints;
//}
//
//Image ComputerVision::CutOutImage(Image img, int row0, int col0, int row1, int col1) {
//
//	uchar* data = img.GetData();
//	uchar* temp = CutOutImageRaw(img.GetRowsNumber(), img.GetColsNumber(), data, row0, col0, row1, col1);
//
//	Image result(row1 - row0 + 1, col1 - col0 + 1, temp, true);
//
//	delete[] temp;
//	delete[] data;
//
//	return result;
//}
//
//Image ComputerVision::GaussDefault(Image & img, double sigma, int interpolateType) {
//	return Gauss(img, sigma, 3*sigma, interpolateType);
//}
//
//double * ComputerVision::GetGaussWeight(double sigma, int size) {
//
//	int sizeArr = size * size;
//	double* data = new double[sizeArr];
//
//	double alpha = 1.0 / (2 * PI()*pow(sigma, 2));
//	double center = double(size - 1) / 2.0;
//	//int pos;
//
//	for (int i = 0; i < size; i++) {
//		for (int j = 0; j < size; j++) {
//
//			data[i* size + j] = alpha * std::exp(-(pow(i - center, 2) + pow(j - center, 2)) / (2.0 * sigma * sigma));
//		}
//	}
//
//	return data;
//}
//
//Image ComputerVision::Gauss(Image & img, double sigma, int k, int interpolateType) {
//	
//	double* data = img.GetNormalizeDataF();
//	double* temp = GaussRaw(img.GetRowsNumber(), img.GetColsNumber(), data, sigma, k, interpolateType);
//
//	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);
//
//	delete[] temp;
//	delete[] data;
//
//	return result;
//}
//
//std::vector<Dot> ComputerVision::Harris(Image & img, int wk, int localMinK, double Threshold, int pointsNeeded, int PartDerivativeType) {
//
//	double* data = img.GetNormalizeDataF();
//	vector<Dot> result = HarrisRaw(img.GetRowsNumber(), img.GetColsNumber(), data, wk, localMinK, Threshold, pointsNeeded, PartDerivativeType);
//
//	delete[] data;
//
//	return result;
//}
//
//double * ComputerVision::HarrisResponse(int rows, int cols, double * A, double * B, double * C, int harrisResponseType) {
//	
//	int size = rows * cols;
//	double* result = new double[size];
//	double k = 0.05; // for base
//
//	switch (harrisResponseType) {
//	//in work!
//	/*
//	case kHarrisResponseDirect:
//
//		delete[] result;
//		result = NULL;
//		break;
//	*/
//
//	case kHarrisResponseBase:
//
//		for (int i = 0; i < size; i++) {
//			result[i] = A[i] * C[i] - pow(B[i], 2) - k * pow(A[i] + C[i], 2);
//		}
//
//		break;
//
//	case kHarrisResponseForstner:
//
//		for (int i = 0; i < size; i++) {
//			if (A[i] * C[i] == 0)
//				result[i] = 0;
//			else
//				result[i] = (A[i] * C[i] - pow(B[i], 2)) / (A[i] + C[i]);
//		}
//
//		break;
//	}
//
//	return result;
//}
//
//std::vector<Dot> ComputerVision::Moravec(Image & img, int wk, int d, int p, double Threshold) {
//	
//	double* data = img.GetNormalizeDataF();
//	vector<Dot> result = MoravecRaw(img.GetRowsNumber(), img.GetColsNumber(), data, wk, d, p, Threshold);
//
//	delete[] data;
//
//	return result;
//}
//
//Image ComputerVision::Sobel(Image & img, int PartDerivativeType, int interpolateType) {
//
//	double* data = img.GetNormalizeDataF();
//	double* temp = SobelRaw(img.GetRowsNumber(), img.GetColsNumber(), data, PartDerivativeType, interpolateType);
//
//	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);
//
//	delete[] temp;
//	delete[] data;
//
//	return result;
//}
