#pragma once
#include "Image.hpp"
#include "ComputerVision.hpp"

class ImagePyramid {
private:
	int octavesNum;
	int layersNum;
	int layersAdditionalNum;
	double sigmaA;
	double sigma0;

	Image* images;

	double sigmaInterval;
	int imagesNum;
public:
	double* sigmaArray;
	ImagePyramid();
	ImagePyramid(Image img, int octavesNum, int layersNum, double sigmaA, double sigma0, int layersAdditionalNum = 0);
	ImagePyramid(Image img, int layersNum, double sigmaA, double sigma0, int layersAdditionalNum = 0);
	~ImagePyramid();

	Image GetImage(int octave, int layer) const;
	Image GetNearestImage(double sigma);
	int GetOctavesNum() const;
	int GetLayersNum() const;
	int GetLayersAdditionalNum() const;
	int GetImagesNum() const;
	double GetSigma0() const;
	double GetSigmaInterval() const;

	uchar L(int x, int y, double sigma);

	int FindNearestSigma(double sigma);
	//Image L(double sigmaArray);
};

