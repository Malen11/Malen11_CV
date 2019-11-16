#include "stdafx.h"
#include "ScaleSpace.hpp"

using namespace std;
using namespace CV_labs;

//Default constructor.
CV_labs::ScaleSpace::ScaleSpace() {

	prependOctavesNum = 0;
	prependLayersNum = 0;
	appendLayersNum = 0;

	octavesNum = 0;
	layersNum = 0;	

	realOctavesNum = 0;
	realLayersNum = 0;	

	sigmaA = 0;			
	sigma0 = 0;		

	sigmaInterval = 0;

	images = nullptr;	
	imagesNum = 0;	
}

ScaleSpace::ScaleSpace(const Image &image, int octavesNum, int layersNum, double sigmaA, double sigma0, int prependOctavesNum, int prependLayersNum, int appendLayersNum) {

	if (octavesNum <= 0)
		throw std::invalid_argument("octavesNum must be > 0");
	if (layersNum <= 0)
		throw std::invalid_argument("layersNum must be > 0");
	if (prependOctavesNum < 0)
		throw std::invalid_argument("prependOctavesNum must be >= 0");
	if (prependLayersNum < 0)
		throw std::invalid_argument("prependLayers must be >= 0");
	if (appendLayersNum < 0)
		throw std::invalid_argument("appendLayers must be >= 0");
	if (sigma0 <= 0)
		throw std::invalid_argument("sigma0 must be > 0");
	if (sigma0 <= sigmaA)
		throw std::invalid_argument("sigma0 must be > sigmaA");

	this->prependOctavesNum	= prependOctavesNum;
	this->prependLayersNum	= prependLayersNum;
	this->appendLayersNum	= appendLayersNum;

	this->octavesNum		= octavesNum;
	this->layersNum			= layersNum;
	this->realOctavesNum	= prependOctavesNum + octavesNum;
	this->realLayersNum		= prependLayersNum + layersNum + appendLayersNum;

	this->imagesNum			= (prependOctavesNum + octavesNum) * (prependLayersNum + layersNum + appendLayersNum);
	this->images			= new Image[imagesNum];
	this->sigmaA			= sigmaA;
	this->sigma0			= sigma0;
	//this->realSigmas		= new double[imagesNum];

	this->sigmaInterval = std::pow(2, 1. / layersNum);
	//calculate started local sigma for each octave
	//double sigmaStartOffset =  std::pow(sigmaInterval, prependLayersNum);
	//double sigmaEndOffset =  std::pow(sigmaInterval, appendLayersNum);

	if (sigma0 / std::pow(sigmaInterval, prependLayersNum) <= sigmaA) {
		throw std::invalid_argument("sigmaStart must be > sigmaA");
	}

	double sigmaLocal, sigmaReal;
	int pos;

	//Create first image
	sigmaLocal = sigma0 / std::pow(sigmaInterval, prependLayersNum);
	sigmaReal = sigmaLocal;

	if (prependOctavesNum > 0) {

		sigmaReal = sigma0 / std::pow(sigmaInterval, prependLayersNum) / std::pow(2, prependOctavesNum);

		images[0] = ImageFilters::ApplyFilter(
			image.GetUpsampledImage(prependOctavesNum), 
			ImageFilters::GenerateGaussSeparableCore(sqrt(sigmaReal * sigmaReal - sigmaA * sigmaA)),
			ImageFilters::kInterpolateZero
		);
	}
	else {

		images[0] = ImageFilters::ApplyFilter(
			image,
			ImageFilters::GenerateGaussSeparableCore(sqrt(sigmaReal * sigmaReal - sigmaA * sigmaA)),
			ImageFilters::kInterpolateBorder
		);
	}


	for (int i = 0; i < realOctavesNum; i++) {

		//sigmaReal = sigmaStart * std::pow(2, i);

		for (int j = 1; j < realLayersNum; j++) {

			pos = i * (realLayersNum) + j;


			images[pos] = ImageFilters::ApplyFilter(
				images[pos - 1],
				ImageFilters::GenerateGaussSeparableCore(sqrt(pow(sigmaLocal * sigmaInterval, 2) - pow(sigmaLocal, 2))),
				ImageFilters::kInterpolateZero
			); 

			sigmaLocal *= sigmaInterval;
			sigmaReal *= sigmaInterval;
		}

		//if we have next octave, we downsample image
		if (i < realOctavesNum - 1) {

			//save downsample image in next position
			++pos;

			sigmaLocal /= std::pow(sigmaInterval, prependLayersNum + appendLayersNum - 1);
			sigmaReal /= std::pow(sigmaInterval, prependLayersNum + appendLayersNum - 1);

			images[pos]	= images[pos - 1 - appendLayersNum - prependLayersNum].GetDownsampledImage(2);

			sigmaLocal = sigma0 / std::pow(sigmaInterval, prependLayersNum);
		}
	}
}

ScaleSpace::ScaleSpace(const Image &image, int layersNum, double sigmaA, double sigma0, int prependOctavesNum, int prependLayersNum, int appendLayersNum)
	:ScaleSpace(image, std::min(std::log2(image.GetRowsNumber()), std::log2(image.GetColsNumber())), sigmaA, sigma0, prependOctavesNum, prependLayersNum, appendLayersNum) {
}

ScaleSpace::~ScaleSpace() {

	if (imagesNum > 0) {
		delete[] images;
	}
}

//Get image by sigma
Image ScaleSpace::GetImage(double sigma) const {

	int pos = getImageIndex(sigma);

	return images[pos];
}

//Get image at octave and layer
Image ScaleSpace::GetImage(int octave, int layer) const {

	if (octave < -prependOctavesNum || octave >= octavesNum) {
		throw new std::invalid_argument("Try to get image at nonexistent octave");
	}
	if (layer < -prependLayersNum || layer >= layersNum + appendLayersNum) {
		throw new std::invalid_argument("Try to get image at nonexistent layer");
	}

	return images[(prependOctavesNum + octave) * realLayersNum + (prependLayersNum + layer)];
}

//calculate DoG
double * ScaleSpace::DoG(int octave, int layer) const {

	if (octave < -prependOctavesNum || octave >= octavesNum) {
		throw new std::invalid_argument("Try to get image at nonexistent octave");
	}
	if (layer < -prependLayersNum || layer >= layersNum + appendLayersNum) {
		throw new std::invalid_argument("Try to get image at nonexistent layer");
	}

	int imageInd0, imageInd1, size;
	double sigma0, sigma1;
	double* data0, *data1, *dataDif;

	imageInd0 = (prependOctavesNum + octave) * realLayersNum + (prependLayersNum + layer);
	data0 = images[imageInd0].GetDataD();
	size = images[imageInd0].GetSize();
	
	if (layer >= layersNum + appendLayersNum - 1) {
		throw new std::invalid_argument("Can't calculate dog for last image in layer");
		//imageInd1 = (prependOctavesNum + octave + 1) * realLayersNum + (prependLayersNum + appendLayersNum + 1);
	}
	else {

		imageInd1 = imageInd0 + 1;
		data1 = images[imageInd1].GetDataD();
	}

	sigma0 = GetImageSigma(octave, layer);
	sigma1 = GetImageSigma(octave, layer + 1);

	dataDif = new double[size];
	double k = sigma0 / sigmaInterval;
	////38-74, 39-74???
	//double stop1 = k *(data1[38 * 400 + 74] - data0[38 * 400 + 74]);
	//double stop2 = k* (data1[39 * 400 + 74] - data0[39 * 400 + 74]);
	for (int i = 0; i < size; i++) {
		dataDif[i] = k * (data1[i] - data0[i]);
	}

	delete[] data0, data1;

	return dataDif;
}

//
////Calculate difference of gaussians
//double ScaleSpace::DoG(Point point, int octave, int layer) const {
//	ScaleSpace DoG;
//
//	DoG.prependOctavesNum = this->prependOctavesNum;
//	DoG.prependLayersNum = this->prependLayersNum;
//
//	if (DoG.appendLayersNum <= 0) {
//		throw new std::invalid_argument("Scale space must contain atleast 1 append layer");
//	}
//
//	DoG.appendLayersNum	= this->appendLayersNum - 1;
//	DoG.layersNum		= this->layersNum;
//	DoG.sigmaInterval	= this->sigmaInterval;
//	DoG.sigma0			= this->sigma0;
//	DoG.sigmaA			= this->sigmaA;
//	DoG.imagesNum		= this->imagesNum;
//	DoG.realLayersNum	= DoG.prependLayersNum + DoG.layersNum + DoG.appendLayersNum;
//	DoG.realOctavesNum	= DoG.prependOctavesNum + DoG.octavesNum;
//	DoG.images			= new Image[imagesNum];
//
//	double * imageData0, * imageData1;
//
//	for (int octaveIndex = 0; octaveIndex < this->realOctavesNum; octaveIndex++) {
//
//		imageData0 = this->GetImage(octaveIndex, 0).GetDataD();
//
//		for (int layerIndex = 0; layerIndex < this->realLayersNum - 1; layerIndex++) {
//
//			imageData1 = this->GetImage(octaveIndex, layerIndex + 1).GetDataD();
//
//				for (int i = 0; i < this->GetImage(octaveIndex, layerIndex).GetSize(); i++) {
//
//			}
//		}
//	}
//
//	return DoG;
//}

//
//Image ScaleSpace::GetNearestImage(double sigma) {
//
//	return images[FindNearestSigma(sigma)];
//}
//
//uchar ScaleSpace::L(int x, int y, double sigma) {
//
//	int pos = FindNearestSigma(sigma);
//	int octave = pos / (1 + layersNum + layersAdditionalNum);
//
//	return images[pos].GetValueAt(y / (1+octave), x/(1+octave));
//}
//
//int ScaleSpace::FindNearestSigma(double sigma) {
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
