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

ImagePyramid::ImagePyramid(Image img, int octavesNum, int layersNum, double sigmaA, double sigma0) {

	if (octavesNum <= 0)
		throw std::invalid_argument("octavesNum must be >= 0");
	if (layersNum <= 0)
		throw std::invalid_argument("layersNum must be >= 0");
	if (sigma0 <= 0)
		throw std::invalid_argument("sigma0 must be != 0");

	this->octavesNum = octavesNum;
	this->layersNum = layersNum;
	this->imagesNum = octavesNum * layersNum;
	this->images = new Image[imagesNum];
	this->sigmas = new double[imagesNum];
	this->sigmaInterval = std::pow(2, 1. / layersNum);
	
	double sigma;
	int octavePos;

	for (int i = 0; i < octavesNum; i++) {

		octavePos = i * layersNum;
		if (i == 0) {
			images[0] = ComputerVision::GaussDefault(img, sqrt(sigma0*sigma0 - sigmaA * sigmaA),ComputerVision::kInterpolateReflection);
			sigmas[0] = sigma0;
		}
		else {
			images[octavePos] = images[octavePos - 1].GetDownsampleImage(2);
			sigmas[octavePos] = 2 * sigmas[octavePos - layersNum];
		}
		sigma = sigma0;
		for (int j = 1; j < layersNum; j++) {

			images[octavePos + j] = ComputerVision::GaussDefault(images[octavePos + j - 1], sqrt(pow(sigma*sigmaInterval,2) - pow(sigma,2)), ComputerVision::kInterpolateReflection);

			sigma *= sigmaInterval;
			sigmas[octavePos + j] = sigma;
		}
	}
}

ImagePyramid::ImagePyramid(Image img, int layersNum, double sigmaA, double sigma0)
	:ImagePyramid(img, std::min(std::log2(img.GetRowsNumber()), std::log2(img.GetColsNumber())), layersNum, sigmaA, sigma0) {
}


ImagePyramid::~ImagePyramid() {

}

Image ImagePyramid::GetImage(int octave, int layer) const {

	int pos = octave * layersNum + layer;

	if (pos >= imagesNum || octave >= octavesNum || layer >= layersNum) {
		throw std::out_of_range("Called not existend Image");
	}

	return images[pos];
}

int ImagePyramid::GetOctavesNum() const {
	
	return octavesNum;
}

int ImagePyramid::GetLayersNum() const {

	return layersNum;
}

int ImagePyramid::GetImagesNum() const {

	return imagesNum;
}

double ImagePyramid::GetSigmaInterval() const {

	return sigmaInterval;
}

Image ImagePyramid::L(double sigma) {
	/*
	int k=0;
	double tempSigma = sigma;

	while (tempSigma > sigmaInterval) {
		tempSigma /= sigmaInterval;
		++k;
	}

	int octave = k / layersNum;
	int layer = k % layersNum;
	//double sigmaPos = sigma0;
	*/

	for (int i = 0; i < imagesNum; i++) {

		//sig
	}
	

	return GetImage(0,0);
}
