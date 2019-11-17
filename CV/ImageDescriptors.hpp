#pragma once
#include "stdafx.h"
#include "Image.hpp"
#include "ImageFilters.hpp"
#include "ScaleSpace.hpp"

namespace CV_labs {

	struct Descriptor {

		Point point;
		int featuresNum;
		double* features;
	};

	//Column (for Histogram)
	struct HistogramInterval {

		double data;
		double midVal;
	};

	//Histogram (for descriptor)
	struct Histogram {

		int intervalsNum;
		HistogramInterval* intervals;
	};

	class ImageDescriptorMethods {
	public:
		static const int kDescriptorTypeSquare = 101;			//Descriptor type (simple square descriptor)

		static const int kDescriptorNormalization2Times = 201;	//Descriptor normalization type (normalize, thresh, normalize)
		static const int kDescriptorNormalizationSimple = 202;	//Descriptor normalization type (normalize 1 time only)
		static const int kDescriptorNormalizationNone = 203;	//Descriptor normalization type (not normalize)

		static const int kDescriptorsComparisonEuclid = 301;	//Descriptor comparison type (Euclid distance)
		static const int kDescriptorsComparisonManhattan = 302;	//Descriptor comparison type (Manhattan metric)
		static const int kDescriptorsComparisonSSD = 303;		//Descriptor comparison type (Sum of squared distances)

		static const int kDescriptorsMatchingBase = 401;	//Descriptor matching type (by max matching, have multiselection!)
		static const int kDescriptorsMatchingNNDR = 402;	//Descriptor matching type (Next Nearest Distance Ratio)
		static const int kDescriptorsMatchingMutal = 403;	//Descriptor matching type (Mutal choice)

		//create Descriptors
		static std::vector<Descriptor> CreateDescriptors(ScaleSpace& scaleSpace, std::vector<ScalePoint> points, int DescriptorSizeD, int histogramNumD, int intervalsNum, int descriptorType = kDescriptorTypeSquare, int descriptorNormalizationType = kDescriptorNormalizationNone);

		//create Descriptors
		static std::vector<Descriptor> CreateDescriptors(Image& image, std::vector<Point> points, int DescriptorSizeD, int histogramNumD, int intervalsNum, int descriptorType = kDescriptorTypeSquare, int descriptorNormalizationType = kDescriptorNormalizationNone);
		
		//create Descriptors
		static std::vector<Descriptor> CreateDescriptorsRaw(int rows, int cols, double * data, std::vector<Point> points, int DescriptorSizeD, int histogramNumD, int intervalsNum, int descriptorType = kDescriptorTypeSquare, int descriptorNormalizationType = kDescriptorNormalizationNone);
		
		//create Square Descriptor
		static Descriptor CreateSquareDescriptorRaw(Point point, int rows, int cols, Gradient * gradients, int DescriptorSizeD, int histogramNumD, int intervalsNum, double alpha = 0);
	
		//create Histigram
		static Histogram CreateHistogram(std::vector<Gradient> gradients, double intervalsMinVal, double intervalsMaxVal, int intervalsNum);

		//calculate distanse beatween 2 descriptors 
		static double DescriptorsDifference(Descriptor desc0, Descriptor desc1, int descriptorsComparisonType = kDescriptorsComparisonEuclid);

		//matching descriptors from descs1 to descs2
		static int* DescriptorsMatchingRaw(std::vector<Descriptor> descriptors1, std::vector<Descriptor> descriptors2, int descriptorsComparisionType = kDescriptorsComparisonEuclid, int descriptorsMatchingType = kDescriptorsMatchingBase, double thresh = 0.8);

		//matching descriptors from descs1 to descs2
		static std::vector<Line> DescriptorsMatching(std::vector<Descriptor> descriptors1, std::vector<Descriptor> descriptors2, int descriptorsComparisionType = kDescriptorsComparisonEuclid, int descriptorsMatchingType = kDescriptorsMatchingBase, double thresh = 0.8);

		//calculate dominating angles
		static std::vector<double> CalculateDominatingAngle(Point point, int rows, int cols, const Gradient * gradients, int windowSizeD, int intervalsNum = 36, double thresh = 0.8, int alphaNumMax = 2);
	};
}