#include "stdafx.h"
#include "ImagePyramid.hpp"

using namespace std;
using namespace CV_labs;

#pragma region Constructors & Destructor

//Default constructor.
ImagePyramid::ImagePyramid(const Image& image, double sigmaA, double sigma0, int octavesNum, int layersNum) {

	if (octavesNum <= 0)
		throw std::invalid_argument("octavesNum must be >= 0");
	if (layersNum <= 0)
		throw std::invalid_argument("layersNum must be >= 0");
	if (sigma0 <= 0)
		throw std::invalid_argument("sigma0 must be != 0");

	this->octavesNum	= octavesNum;
	this->layersNum		= layersNum;
	this->imagesNum		= octavesNum * (layersNum + 1);//lyersNum + base layer image
	this->images		= new Image[imagesNum];
	this->sigmaA		= sigmaA;
	this->sigma0		= sigma0;
	//this->realSigmas	= new double[imagesNum];

	this->sigmaInterval = std::pow(2, 1. / layersNum);
	double sigmaLocal, sigmaReal;
	int pos;

	//Blur image to sigma0
	images[0] = ImageFilters::ApplyFilter(image, ImageFilters::GenerateGaussSeparableCore(sqrt(sigma0 * sigma0 - sigmaA * sigmaA)), ImageFilters::kInterpolateReflection);
	sigmaReal = sigmaLocal = sigma0;

	for (int i = 0; i < octavesNum; i++) {

		//+1 because each level have layersNum images + original layer image
		for (int j = 1; j < layersNum + 1; j++) {

			pos = i * (layersNum + 1) + j;

			images[pos] = ImageFilters::ApplyFilter(
				images[pos - 1],
				ImageFilters::GenerateGaussSeparableCore(sqrt(pow(sigmaLocal * sigmaInterval, 2) - pow(sigmaLocal, 2))),
				ImageFilters::kInterpolateReflection
			); 

			sigmaLocal *= sigmaInterval;
			sigmaReal *= sigmaInterval;
		}

		//if we have next octave, we downsample image
		if (i < octavesNum - 1) {

			//save downsample image in next position
			++pos;

			images[pos]	= images[pos - 1].GetDownsampledImage(2);
			sigmaLocal	= sigma0;
		}
	}
}

ImagePyramid::ImagePyramid(const Image& img, double sigmaA, double sigma0, int layersNum)
	:ImagePyramid(img, sigmaA, sigma0, std::min(std::log2(img.GetRowsNumber()), std::log2(img.GetColsNumber())), layersNum) {
}

ImagePyramid::~ImagePyramid() {

}

//Get image by sigma
Image ImagePyramid::GetImage(double sigma) const {

	int pos = getImageIndex(sigma);

	return images[pos];
}

//Get image at octave and layer
Image ImagePyramid::GetImage(int octave, int layer) const {

	if (octave >= octavesNum) {
		throw new std::invalid_argument("Try to get image at nonexistent octave");
	}
	if (layer >= layersNum) {
		throw new std::invalid_argument("Try to get image at nonexistent layer");
	}

	return images[octave * (layersNum + 1) + layer];
}

#pragma endregion
//
//Image ImagePyramid::GetNearestImage(double sigma) {
//
//	return images[FindNearestSigma(sigma)];
//}
//
//uchar ImagePyramid::L(int x, int y, double sigma) {
//
//	int pos = FindNearestSigma(sigma);
//	int octave = pos / (1 + layersNum + layersAdditionalNum);
//
//	return images[pos].GetValueAt(y / (1+octave), x/(1+octave));
//}
//
//int ImagePyramid::FindNearestSigma(double sigma) {
//	
//	//border values
//	if (sigma < sigmaArray[0])
//		return 0;
//	else if (sigma > sigmaArray[imagesNum - 1 - layersAdditionalNum])
//		return imagesNum - 1 - layersAdditionalNum;
//	
//	//calculate octave and layer
//	int octave = std::log2(sigma / sigma0);
//	int layer = std::log2(sigma / (sigma0 * std::pow(2, octave))) / std::log2(sigmaInterval);
//
//	//if next value more simillary, then take next
//	if (layer < layersNum) {
//
//		double left = std::abs(sigma - sigmaArray[octave*(1 + layersNum + layersAdditionalNum) + layer]);
//		double right = std::abs(sigma - sigmaArray[octave*(1 + layersNum + layersAdditionalNum) + layer + 1]);
//		
//		if (left > right)
//			++layer;
//	}
//
//	//if we move to last layer on octave, take next octave's first layer
//	if (layer == layersNum && octave != octavesNum-1) {
//		
//		++octave;
//		layer = 0;
//	}
//
//	return octave * (1 + layersNum + layersAdditionalNum) + layer;
//}
