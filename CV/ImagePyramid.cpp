#include "stdafx.h"
#include "ImagePyramid.hpp"


ImagePyramid::ImagePyramid() {
	
	octavesNum = 0;
	layersNum = 0;
	sigmaA = 0.;
	sigma0 = 0.;

	images = NULL;

	sigmaInterval = 0.;
	imagesNum = 0;
}

ImagePyramid::ImagePyramid(Image img, int octavesNum, int layersNum, double sigmaA, double sigma0, int layersAdditionalNum) {

	if (octavesNum <= 0)
		throw std::invalid_argument("octavesNum must be >= 0");
	if (layersNum <= 0)
		throw std::invalid_argument("layersNum must be >= 0");
	if (sigma0 <= 0)
		throw std::invalid_argument("sigma0 must be != 0");

	this->octavesNum = octavesNum;
	this->layersNum = layersNum;
	this->layersAdditionalNum = layersAdditionalNum;
	this->imagesNum = octavesNum * (1 + layersNum + layersAdditionalNum);
	this->images = new Image[imagesNum];
	this->sigma0 = sigma0;
	this->sigmaA = sigmaA;
	this->sigmaArray = new double[imagesNum];
	this->sigmaInterval = std::pow(2, 1. / layersNum);
	
	double sigmaLocal, sigmaReal;
	int octavePos;

	for (int i = 0; i < octavesNum; i++) {

		octavePos = i * (1+layersNum + layersAdditionalNum);

		if (i == 0) {

			images[0] = ComputerVision::GaussDefault(img, sqrt(sigma0*sigma0 - sigmaA * sigmaA),ComputerVision::kInterpolateReflection);
			sigmaLocal = sigmaReal = sigmaArray[0] = sigma0;
		}
		else {

			images[octavePos] = images[octavePos - 1- layersAdditionalNum].GetDownsampleImage(2);
			sigmaReal = sigmaArray[octavePos] = 2 * sigmaArray[octavePos - 1 - (layersNum + layersAdditionalNum)];
			sigmaLocal = sigma0;
		}


		for (int j = 1; j <= layersNum + layersAdditionalNum; j++) {

			images[octavePos + j] = ComputerVision::GaussDefault(images[octavePos + j - 1], sqrt(pow(sigmaLocal*sigmaInterval, 2) - pow(sigmaLocal, 2)), ComputerVision::kInterpolateReflection);

			sigmaLocal *= sigmaInterval;
			sigmaReal *= sigmaInterval;
			sigmaArray[octavePos + j] = sigmaReal;
		}
	}
}

ImagePyramid::ImagePyramid(Image img, int layersNum, double sigmaA, double sigma0, int layersAdditionalNum)
	:ImagePyramid(img, std::min(std::log2(img.GetRowsNumber()), std::log2(img.GetColsNumber())), layersNum, sigmaA, sigma0, layersAdditionalNum) {
}


ImagePyramid::~ImagePyramid() {

}

Image ImagePyramid::GetImage(int octave, int layer) const {

	int pos = octave * (layersNum+1+layersAdditionalNum) + layer;

	if (pos >= imagesNum || octave >= octavesNum || layer >= (1+layersNum+layersAdditionalNum)) {
		throw std::out_of_range("Called not existend Image");
	}

	return images[pos];
}

Image ImagePyramid::GetNearestImage(double sigma) {

	return images[FindNearestSigma(sigma)];
}

int ImagePyramid::GetOctavesNum() const {
	
	return octavesNum;
}

int ImagePyramid::GetLayersNum() const {

	return layersNum;
}

int ImagePyramid::GetLayersAdditionalNum() const {

	return layersAdditionalNum;
}

int ImagePyramid::GetImagesNum() const {

	return imagesNum;
}

double ImagePyramid::GetSigma0() const {
	
	return sigma0;
}

double ImagePyramid::GetSigmaInterval() const {

	return sigmaInterval;
}

uchar ImagePyramid::L(int x, int y, double sigma) {

	int pos = FindNearestSigma(sigma);
	int octave = pos / (1 + layersNum + layersAdditionalNum);

	return images[pos].GetValueAt(y / (1+octave), x/(1+octave));
}

int ImagePyramid::FindNearestSigma(double sigma) {
	
	//border values
	if (sigma < sigmaArray[0])
		return 0;
	else if (sigma > sigmaArray[imagesNum - 1 - layersAdditionalNum])
		return imagesNum - 1 - layersAdditionalNum;
	
	//calculate octave and layer
	int octave = std::log2(sigma / sigma0);
	int layer = std::log2(sigma / (sigma0 * std::pow(2, octave))) / std::log2(sigmaInterval);

	//if next value more simillary, then take next
	if (layer < layersNum) {

		double left = std::abs(sigma - sigmaArray[octave*(1 + layersNum + layersAdditionalNum) + layer]);
		double right = std::abs(sigma - sigmaArray[octave*(1 + layersNum + layersAdditionalNum) + layer + 1]);
		
		if (left > right)
			++layer;
	}

	//if we move to last layer on octave, take next octave's first layer
	if (layer == layersNum && octave != octavesNum-1) {
		
		++octave;
		layer = 0;
	}

	return octave * (1 + layersNum + layersAdditionalNum) + layer;
}
