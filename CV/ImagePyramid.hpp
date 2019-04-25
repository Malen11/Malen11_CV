#pragma once
#include "Image.hpp"
#include "ComputerVision.hpp"

class ImagePyramid {
private:
	int octavesNum;
	int layersNum;
	double sigmaA;
	double sigma0;

	Image* images;

	double sigmaInterval;
	int imagesNum;
public:
	ImagePyramid();
	ImagePyramid(Image img, int octavesNum, int layersNum, double sigmaA, double sigma0);
	ImagePyramid(Image img, int layersNum, double sigmaA, double sigma0);
	~ImagePyramid();

	Image GetImage(int octave, int layer) const;
	int GetOctavesNum() const;
	int GetLayersNum() const;
	int GetImagesNum() const;
	double GetSigmaInterval() const;
};

