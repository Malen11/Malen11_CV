#include "stdafx.h"
#include "ImageDescriptors.hpp"

using namespace std;
using namespace CV_labs;

//extern double testAlpha;

//create Descriptors
std::vector<Descriptor> CV_labs::ImageDescriptorMethods::CreateDescriptors(Image & image, std::vector<Point> points, int DescriptorSizeD, int histogramNumD, int intervalsNum, int descriptorType, int descriptorNormalizationType) {

	double* data = image.GetDataD();

	//double* normalizedData = NormalizeData(image.GetSize(), data);

	std::vector<Descriptor> result = CreateDescriptorsRaw(image.GetRowsNumber(), image.GetColsNumber(), data, points, DescriptorSizeD, histogramNumD, intervalsNum, descriptorType, descriptorNormalizationType);

	delete[] data;

	return result;
}

std::vector<Descriptor> ImageDescriptorMethods::CreateDescriptorsRaw(int rows, int cols, double * data, std::vector<Point> points, int windowSizeD, int histogramNumD, int intervalsNum, int descriptorType, int descriptorNormalizationType) {

	if (windowSizeD / histogramNumD == 0) {
		throw new std::invalid_argument("HistogramNumD must be grater then windowSizeD");
	}

	std::vector<Descriptor> descriptorsVector;

	//if create simple descriptor
	if (descriptorType = kDescriptorTypeSquare) {

		Gradient* gradients = ImageFilters::CalculateGradientRaw(rows, cols, data, ImageFilters::kPartDerivativeTypeSobelCore, ImageFilters::kInterpolateBorder);

		//for every point create descriptor
		for (int i = 0; i < points.size(); i++) {

			//calculate dominating angles
			std::vector<double> alpha = CalculateDominatingAngle(points[i], rows, cols, gradients, windowSizeD);

			for (int j = 0; j < alpha.size(); j++) {

				cout << '(' << points[i].x << ", " << points[i].y << ") => " << alpha[j] << "\n";
				Descriptor descriptor = CreateSquareDescriptorRaw(points[i], rows, cols, gradients, windowSizeD, histogramNumD, intervalsNum, alpha[j]);
				descriptorsVector.push_back(descriptor);
			}
		}

		delete[] gradients;
	}
	else {
		throw std::invalid_argument("Wrong descriptorType");
	}

	if (descriptorNormalizationType == kDescriptorNormalization2Times) {

		double* temp;
		for (int i = 0; i < descriptorsVector.size(); i++) {

			temp = ImageFilters::NormalizeData(descriptorsVector[i].featuresNum, descriptorsVector[i].features);
			delete[] descriptorsVector[i].features;
			descriptorsVector[i].features = temp;

			for (int j = 0; j < descriptorsVector[i].featuresNum; j++) {
				//> 0.2 или <0.2?
				if (descriptorsVector[i].features[j] > 0.2)
					descriptorsVector[i].features[j] = 0.2;
			}

			temp = ImageFilters::NormalizeData(descriptorsVector[i].featuresNum, descriptorsVector[i].features);
			delete[] descriptorsVector[i].features;
			descriptorsVector[i].features = temp;
		}
	}
	else if (descriptorNormalizationType == kDescriptorNormalizationSimple) {

		double* temp;
		for (int i = 0; i < descriptorsVector.size(); i++) {

			temp = ImageFilters::NormalizeData(descriptorsVector[i].featuresNum, descriptorsVector[i].features);
			delete[] descriptorsVector[i].features;
			descriptorsVector[i].features = temp;
		}
	}
	else {
		cout << "descriptorNormalizationType not set. Descriptors not normalize";
	}


	return descriptorsVector;
}

Descriptor CV_labs::ImageDescriptorMethods::CreateSquareDescriptorRaw(Point point, int rows, int cols, Gradient * gradients, int DescriptorSizeD, int histogramNumD, int intervalsNum, double alpha) {

	int size = rows * cols;
	int histogramsNum = histogramNumD * histogramNumD;

	//create window
	int x0, x1, y0, y1;

	x0 = point.x - DescriptorSizeD / 2;
	x1 = x0 + DescriptorSizeD;
	y0 = point.y - DescriptorSizeD / 2;
	y1 = y0 + DescriptorSizeD;

	//calculate gauss weight for this part of image with 2d-size = rowsNew x rowsNew
	double* gaussWeight = ImageFilters::GenerateGaussCore(0.5 * DescriptorSizeD, DescriptorSizeD / 2).data;

	Descriptor descriptor;
	descriptor.point = point;
	descriptor.featuresNum = histogramsNum * intervalsNum;
	descriptor.features = new double[descriptor.featuresNum];

	//grid division step by histograms
	double histStepX = double(DescriptorSizeD) / histogramNumD;
	double histStepY = double(DescriptorSizeD) / histogramNumD;

	//arrays of gradients, that belong to specifical histogram
	vector<Gradient>* histogramsGradients = new vector<Gradient>[histogramsNum];

	//cos(alpha), sin(alpha)
	double cosAlpha = cos(alpha);
	double sinAlpha = sin(alpha);

	//error for comparing
	double eps = 0.0001;

	double  descriptorCenterX = x0 + (double)(x1 - x0) / 2.0;	// descriptor center
	double  descriptorCenterY = y0 + (double)(y1 - y0) / 2.0;	// descriptor center

	int histogramPos, histogramXIndex, histogramYIndex;
	int posAbs, posRel, rowAbs, colAbs;
	double xCenter, yCenter, xCenterRotated, yCenterRotated, xCoef, yCoef;

	for (int rowRel = 0; rowRel < DescriptorSizeD; rowRel++) {
		for (int colRel = 0; colRel < DescriptorSizeD; colRel++) {

			rowAbs = y0 + rowRel;
			colAbs = x0 + colRel;

			//pixel position on image
			posAbs = rowAbs * cols + colAbs;
			posRel = rowRel * DescriptorSizeD + colRel;

			//0.5 - pixel's center
			xCenter = colAbs + 0.5;
			yCenter = rowAbs + 0.5;

			//pixel center after rotate
			xCenterRotated = descriptorCenterX + (xCenter - descriptorCenterX) * cosAlpha - (yCenter - descriptorCenterY) * sinAlpha;
			yCenterRotated = descriptorCenterY + (xCenter - descriptorCenterX) * sinAlpha + (yCenter - descriptorCenterY) * cosAlpha;

			//calculate left-top hist index
			histogramXIndex = (xCenterRotated - 0.5 - x0) / histStepX;
			histogramYIndex = (yCenterRotated - 0.5 - y0) / histStepY;

			if (histogramXIndex < - 1 || histogramXIndex > histogramNumD - 1 || histogramYIndex < - 1 || histogramYIndex > histogramNumD - 1) {
				continue;
			}

			//calc distance between left-top point and right border of gist (max 1)
			xCoef = std::min(x0 + (histogramXIndex + 1) * histStepX - (xCenterRotated - 0.5), 1.0);
			yCoef = std::min(y0 + (histogramYIndex + 1) * histStepY - (yCenterRotated - 0.5), 1.0);

			//gradient by coordinate
			Gradient pixelGradient = ImageFilters::GetVirtualPixelGradientZero(rowAbs, colAbs, rows, cols, gradients);
			pixelGradient.value = pixelGradient.value * gaussWeight[posRel];
			pixelGradient.phi = pixelGradient.phi - alpha;

			if (pixelGradient.phi < 0) {

				pixelGradient.phi += 2 * ImageFilters::PI();
			}
			else if (pixelGradient.phi > 2 * ImageFilters::PI()) {

				pixelGradient.phi -= 2 * ImageFilters::PI();
			}

			//left top
			if (histogramXIndex >= 0 && histogramXIndex < histogramNumD && histogramYIndex >= 0 && histogramYIndex < histogramNumD) {
				
				Gradient pixelPartGradient = { pixelGradient.value * xCoef * yCoef, pixelGradient.phi };

				histogramsGradients[histogramYIndex * histogramNumD + histogramXIndex].push_back(pixelPartGradient);
			}

			//right top
			if (histogramXIndex + 1 >= 0 && histogramXIndex + 1 < histogramNumD && histogramYIndex >= 0 && histogramYIndex < histogramNumD) {

				Gradient pixelPartGradient = { pixelGradient.value * (1 - xCoef) * yCoef, pixelGradient.phi };

				histogramsGradients[histogramYIndex * histogramNumD + histogramXIndex + 1].push_back(pixelPartGradient);
			}

			//left bot
			if (histogramXIndex >= 0 && histogramXIndex < histogramNumD && histogramYIndex + 1 >= 0 && histogramYIndex + 1 < histogramNumD) {

				Gradient pixelPartGradient = { pixelGradient.value * xCoef * (1 - yCoef), pixelGradient.phi };

				histogramsGradients[(histogramYIndex + 1) * histogramNumD + histogramXIndex].push_back(pixelPartGradient);
			}

			//right bot
			if (histogramXIndex + 1 >= 0 && histogramXIndex + 1 < histogramNumD && histogramYIndex + 1 >= 0 && histogramYIndex + 1 < histogramNumD) {

				Gradient pixelPartGradient = { pixelGradient.value * (1 - xCoef) * (1 - yCoef), pixelGradient.phi };
				
				histogramsGradients[(histogramYIndex + 1) * histogramNumD + histogramXIndex + 1].push_back(pixelPartGradient);
			}
		}
	}

	for (int i = 0; i < histogramsNum; i++) {

		Histogram hist = CreateHistogram(histogramsGradients[i], 0, 2 * ImageFilters::PI(), intervalsNum);

		for (int j = 0; j < hist.intervalsNum; j++) {
			descriptor.features[i * intervalsNum + j] = hist.intervals[j].data;
		}
	}

	return descriptor;
}

Histogram ImageDescriptorMethods::CreateHistogram(std::vector<Gradient> gradients, double intervalsMinVal, double intervalsMaxVal, int intervalsNum) {

	Histogram hist;
	hist.intervalsNum = intervalsNum;
	hist.intervals = new HistogramInterval[intervalsNum];

	double step = (intervalsMaxVal - intervalsMinVal) / intervalsNum;
	double startValue = step / 2;		//because we use mid value (half step)

	for (int i = 0; i < intervalsNum; i++) {

		hist.intervals[i].data = 0;
		hist.intervals[i].midVal = startValue + i * step;
	}

	double k;
	int  column0, column1;// pos;

	for (int i = 0; i < gradients.size(); i++) {

		//edged right value
		if (gradients[i].phi >= hist.intervals[intervalsNum - 1].midVal) {

			column0 = intervalsNum - 1;
			column1 = 0;

			k = (gradients[i].phi - hist.intervals[intervalsNum - 1].midVal) / step;

		}//edged left value
		else if (gradients[i].phi < hist.intervals[0].midVal) {

			column0 = intervalsNum - 1;
			column1 = 0;

			k = 1 - (hist.intervals[0].midVal - gradients[i].phi) / step;							//!!!

		}
		else {
			//find nearest column (basket)
			column0 = (gradients[i].phi - startValue) / step;
			column1 = column0 + 1;

			k = (gradients[i].phi - hist.intervals[column0].midVal) / step;	//!!!
		}

		hist.intervals[column0].data += gradients[i].value * (1 - k);
		hist.intervals[column1].data += gradients[i].value * k;

		//if (hist.intervals[column0].data < 0 || hist.intervals[column1].data < 0)
		//	int wtf = 0;
	}

	return hist;
}

//calculate distanse beatween 2 descriptors 
double CV_labs::ImageDescriptorMethods::DescriptorsDifference(Descriptor desc0, Descriptor desc1, int descriptorsComparisonType) {

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
	else {
		 throw new std::invalid_argument("Unexpected descriptorsComparisonType value");
	}

	return difference;
}

//matching descriptors from descs1 to descs2
int * ImageDescriptorMethods::DescriptorsMatchingRaw(std::vector<Descriptor> descriptors1, std::vector<Descriptor> descriptors2, int descriptorsComparisionType, int descriptorsMatchingType, double thresh) {

	int descriptors1Num = (int)descriptors1.size();
	int descriptors2Num = (int)descriptors2.size();
	double* descsDiff = new double[descriptors1Num * descriptors2Num];

	for (int i = 0; i < descriptors1Num; i++) {
		for (int j = 0; j < descriptors2Num; j++) {
			descsDiff[i * descriptors2Num + j] = DescriptorsDifference(descriptors1[i], descriptors2[j], descriptorsComparisionType);
		}
	}

	int* descriptor1PairIn2 = new int[descriptors1Num];
	int* descriptor2PairIn1;
	int* descriptor1PairIn2NND;

	int matchedIn2, matchedIn1;

	switch (descriptorsMatchingType) {
	case kDescriptorsMatchingNNDR:

		descriptor1PairIn2NND = new int[descriptors1Num];

		for (int i = 0; i < descriptors1Num; i++) {

			descriptor1PairIn2[i] = 0;
			descriptor1PairIn2NND[i] = -1;

			//find nearest descriptor and second nearest
			for (int j = 1; j < descriptors2Num; j++) {
				if (descsDiff[i * descriptors2Num + descriptor1PairIn2[i]] > descsDiff[i * descriptors2Num + j]) {

					descriptor1PairIn2NND[i] = descriptor1PairIn2[i];
					descriptor1PairIn2[i] = j;
				}
			}

			if (descriptor1PairIn2NND[i] != -1) {
				if (descsDiff[i*descriptors2Num + descriptor1PairIn2[i]] / descsDiff[i * descriptors2Num + descriptor1PairIn2NND[i]] > thresh)
					descriptor1PairIn2[i] = -1;
			}
		}

		delete[] descriptor1PairIn2NND;

		break;

	case kDescriptorsMatchingMutal:
		
		descriptor2PairIn1 = new int[descriptors2Num];

		for (int i = 0; i < descriptors1Num; i++) {

			matchedIn2 = 0;

			for (int j = 1; j < descriptors2Num; j++) {
				if (descsDiff[i * descriptors2Num + matchedIn2] > descsDiff[i * descriptors2Num + j]) {
					matchedIn2 = j;
				}
			}

			descriptor1PairIn2[i] = matchedIn2;
		}

		for (int i = 0; i < descriptors2Num; i++) {

			matchedIn1 = 0;

			for (int j = 1; j < descriptors1Num; j++) {
				if (descsDiff[matchedIn1 * descriptors2Num + i] > descsDiff[j * descriptors2Num + i]) {
					matchedIn1 = j;
				}
			}
			descriptor2PairIn1[i] = matchedIn1;
		}

		for (int i = 0; i < descriptors1Num; i++) {
			if (i != descriptor2PairIn1[descriptor1PairIn2[i]]) {
				descriptor1PairIn2[i] = -1;
			}
		}

		delete[] descriptor2PairIn1;

		break;

	case kDescriptorsMatchingBase:

		for (int i = 0; i < descriptors1Num; i++) {

			descriptor1PairIn2[i] = 0;

			for (int j = 1; j < descriptors2Num; j++) {
				if (descsDiff[i * descriptors2Num + descriptor1PairIn2[i]] > descsDiff[i * descriptors2Num + j])
					descriptor1PairIn2[i] = j;
			}
		}

		break;

	default:
		throw new std::invalid_argument("Unexpected descriptorsMatchingType value");
	}

	return descriptor1PairIn2;
}

//matching descriptors from descs1 to descs2
std::vector<Line> ImageDescriptorMethods::DescriptorsMatching(std::vector<Descriptor> descs1, std::vector<Descriptor> descs2, int descriptorsComparisionType, int descriptorsMatchingType, double thresh) {

	int* desc1PairInDesc2 = DescriptorsMatchingRaw(descs1, descs2, descriptorsComparisionType, descriptorsMatchingType, thresh);

	std::vector<Line> matchedPoints;
	Line pairPoints;

	for (int i = 0; i < descs1.size(); i++) {

		cout << "desc 1: " << i;

		if (desc1PairInDesc2[i] != -1) {

			cout << " desc2: " << desc1PairInDesc2[i];

			pairPoints.pointA = descs1[i].point;
			pairPoints.pointB = descs2[desc1PairInDesc2[i]].point;

			matchedPoints.push_back(pairPoints);
		}

		cout << "\n";
	}

	return matchedPoints;
}

//calculate dominating angles
std::vector<double> ImageDescriptorMethods::CalculateDominatingAngle(Point point, int rows, int cols, const Gradient * gradients, int windowSizeD, int intervalsNum, double thresh, int alphaNumMax) {

	int size = rows * cols;

	//create window
	int x0, x1, y0, y1;

	x0 = point.x - windowSizeD / 2;
	x1 = x0 + windowSizeD;
	y0 = point.y - windowSizeD / 2;
	y1 = y0 + windowSizeD;

	//calculate gauss weight for this part of image with 2d-size = rowsNew x rowsNew
	double* gaussWeight = ImageFilters::GenerateGaussCore(1.5 * windowSizeD, windowSizeD / 2.0).data;

	//place all pixels, that belong to window into vector
	vector<Gradient> histogramsGradients;

	for (int row = y0; row < y1; row++) {
		for (int col = x0; col < x1; col++) {

			Gradient pixelGradient = ImageFilters::GetVirtualPixelGradientZero(row, col, rows, cols, gradients);
			pixelGradient.value = pixelGradient.value * gaussWeight[(row - y0) * (x1 - x0) + (col - x0)];

			histogramsGradients.push_back(pixelGradient);
			//histogramsGradients.push_back(ImageFilters::GetVirtualPixelGradientZero(row, col, rows, cols, gradients));
		}
	}

	//create histogram
	Histogram hist = CreateHistogram(histogramsGradients, 0, 2 * ImageFilters::PI(), intervalsNum);

	//find max value in interval
	double maxVal = hist.intervals[0].data;

	for (int i = 1; i < intervalsNum; i++) {

		if (maxVal < hist.intervals[i].data) {

			maxVal = hist.intervals[i].data;
		}
	}

	//find all intervals, that > when maxVal * thresh
	vector<int> pretendetIntervals;

	for (int i = 0; i < intervalsNum; i++) {

		if (thresh * maxVal < hist.intervals[i].data) {

			pretendetIntervals.push_back(i);
		}
	}

	//orientation of decriptor
	vector<double> alphas;

	int tempInterval;
	for (int i = 0; i < pretendetIntervals.size() - 1; i++) {
		for (int j = 0; j < pretendetIntervals.size() - i - 1; j++) {

			if (hist.intervals[pretendetIntervals[j]].data < hist.intervals[pretendetIntervals[j + 1]].data) {

				tempInterval				= pretendetIntervals[j];
				pretendetIntervals[j]		= pretendetIntervals[j + 1];
				pretendetIntervals[j + 1]	= tempInterval;
			}
		}
	}

	//lower value for alpha 
	//int		alphaNumMaxMinPos	= std::min(alphaNumMax, (int)pretendetIntervals.size()) - 1;
	//double	alphaNumMaxMinValue	= hist.intervals[pretendetIntervals[alphaNumMaxMinPos]].data;
	//
	////if we can check next element
	//if ((alphaNumMaxMinPos == alphaNumMax - 1) && (alphaNumMaxMinPos + 1 < pretendetIntervals.size())) {

	//	//we check, is it same as alphaNumMaxMinValue or not
	//	if (hist.intervals[pretendetIntervals[alphaNumMaxMinPos + 1]].data == alphaNumMaxMinValue) {

	//		alphaNumMaxMinValue = alphaNumMaxMinValue + 0.0000001;
	//	}

	//	for (int i = 0; i < pretendetIntervals.size(); i++) {
	//		
	//		if (hist.intervals[pretendetIntervals[i]].data > alphaNumMaxMinValue) {

	//			alphaNumMaxMinPos = i;
	//			break;
	//		}
	//	}
	//}

	for (int i = 0; i < std::min(alphaNumMax, (int)pretendetIntervals.size()); i++) {

		alphas.push_back(hist.intervals[pretendetIntervals[i]].midVal);
	}

	return alphas;
}
