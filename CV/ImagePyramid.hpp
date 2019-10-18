#pragma once
#include "stdafx.h"
#include "Image.hpp"
#include "ImageFilters.hpp"

namespace CV_labs {
	
	//Class for storing collection of images in images pyramid's format
	class ImagePyramid {

	private:
		int octavesNum;				//Number of octaves
		int layersNum;				//Number of layers(not count base layer image)
		//int layersAdditionalNum;	//Number of cross layers (not count in layersNum)
		double sigmaA;				//Original sigma
		double sigma0;				//Starter sigma
		
		double sigmaInterval;		//Sigma interval between layers
		//double* realSigmas;			//Array of real sigmas for each image

		Image* images;				//Images
		int imagesNum;				//Number of Images


	public:

#pragma region Constructors & Destructor

		//Default constructor.
		ImagePyramid(const Image &img, double sigmaA, double sigma0, int octavesNum, int layersNum);

		//Auto-octave constructor.
		ImagePyramid(const Image &img, double sigmaA, double sigma0, int layersNum);

		//Destructor
		~ImagePyramid();

#pragma endregion

		//Get number of octaves
		int GetOctavesNum() const;

		//Get number of layers (base image on each level not count)
		int GetLayersNum() const;

		//Get number of images
		int GetImagesNum() const;

		//Get sigma 0
		double GetSigma0() const;

		//Get image by sigma
		Image GetImage(double sigma) const;

		//Get image at octave and layer
		Image GetImage(int octave, int layer) const;

		//Get sigma at octave and layer
		double GetImageSigma(int octave, int layer) const;

		//Get pixel for sigma at row, col
		uchar L(int row, int col, double sigma) const;

		//Get pixel for sigma at point
		uchar L(Point point, double sigma) const;

		//Image L(double sigmaArray);
	protected:

		//Get converted coordinate to specific sigma
		CV_labs::Point ConvertCoordinate(int row, int col, double sigma) const;

		//Get converted coordinate to specific sigma
		CV_labs::Point ConvertCoordinate(Point point, double sigma) const;
		
		//Get nearest image index to specific sigma
		int getImageIndex(double sigma) const;
	};

	//Get number of octaves
	inline int ImagePyramid::GetOctavesNum() const {

		return this->octavesNum;
	}

	//Get number of layers (base image on each level not count)
	inline int ImagePyramid::GetLayersNum() const {

		return this->layersNum;
	}

	//Get number of images
	inline int ImagePyramid::GetImagesNum() const {

		return this->imagesNum;
	}

	//Get sigma 0
	inline double ImagePyramid::GetSigma0() const {

		return this->sigma0;
	}

	//Get sigma at octave and layer
	inline double ImagePyramid::GetImageSigma(int octave, int layer) const {
		
		if (octave >= octavesNum) {
			throw new std::invalid_argument("Try to get image at nonexistent octave");
		}
		if (layer >= layersNum) {
			throw new std::invalid_argument("Try to get image at nonexistent layer");
		}

		return std::pow(sigmaInterval, layer) * std::pow(2, octave) * sigma0;
	}

	//Get converted coordinate to specific sigma
	inline CV_labs::Point ImagePyramid::ConvertCoordinate(Point point, double sigma) const {

		//each octave contain layersNum + 1 image, so, if we take pos of image with specific sigma, we can calculate octave
		int octave = getImageIndex(sigma) / (layersNum + 1);
		double k = 1 / std::pow(2, octave);

		CV_labs::Point result;
		point.x = std::round(point.x * k);
		point.y = std::round(point.y * k);

		return result;
	}

	//Get converted coordinate to specific sigma
	inline CV_labs::Point ImagePyramid::ConvertCoordinate(int row, int col, double sigma) const {
		
		return ConvertCoordinate({ col, row }, sigma);
	}

	//Get nearest image index to specific sigma
	inline int ImagePyramid::getImageIndex(double sigma) const {

		//if sigma lower then first image's sigma, return first image
		if(sigma < sigma0) {
			return 0;
		}

		//if sigma greater, then sigma of last image, return last image index
		if (sigma > GetImageSigma(octavesNum - 1, layersNum - 1)) {
			return imagesNum - 1;
		}

		//calc coefficient
		double k = sigma / sigma0;
		
		int octave = std::log2(k);
		
		//local sigma's coefficient
		double kl = k / pow(2, octave);

		//round layer to near int
		int layer = std::round(std::log(kl) / std::log(sigmaInterval));

		//if layer rounded to layer > layer - 1, then sigma closer to next octave
		if (layer > layersNum - 1) {
			octave = octave + 1;
			layer = 0;
		}

		return octave * (this->layersNum + 1) + layer;
	}

	//Get pixel for sigma at row, col
	inline uchar ImagePyramid::L(int row, int col, double sigma) const {

		return L({ col, row }, sigma);
	}

	//Get pixel for sigma at point
	inline uchar ImagePyramid::L(Point point, double sigma) const {

		int pos = getImageIndex(sigma);

		//each octave contain layersNum + 1 image, so, if we take pos of image with specific sigma, we can calculate octave
		int octave = pos / (layersNum + 1);
		double k = 1 / std::pow(2, octave);

		CV_labs::Point result;
		point.x = std::round(point.x * k);
		point.y = std::round(point.y * k);

		return images[pos].GetValueAt(point);
	}
}