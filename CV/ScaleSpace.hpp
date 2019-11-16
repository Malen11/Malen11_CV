#pragma once
#include "stdafx.h"
#include "Image.hpp"
#include "ImageFilters.hpp"

namespace CV_labs {
	
	struct ScalePoint {

		Point point;
		double scale;
	};
	//Class for storing collection of images in images pyramid's format
	class ScaleSpace {

	private:
		int prependOctavesNum;	
		int prependLayersNum;			
		int appendLayersNum;

		int octavesNum;				//Number of octaves
		int layersNum;				//Number of layers

		int realOctavesNum;			//Number of octaves
		int realLayersNum;			//Number of layers

		double sigmaA;				//Original sigma
		double sigma0;				//Starter sigma
		
		double sigmaInterval;		//Sigma interval between layers
		//double* realSigmas;			//Array of real sigmas for each image

		Image* images;				//Images
		int imagesNum;				//Number of Images


	public:

#pragma region Constructors & Destructor

		//Default constructor.
		ScaleSpace();

		//Base constructor.
		ScaleSpace(const Image &image, int octavesNum, int layersNum, double sigmaA, double sigma0, int prependOctavesNum = 0, int prependLayersNum = 0, int appendLayersNum = 0);

		//Auto-octave constructor.
		ScaleSpace(const Image &image, int layersNum, double sigmaA, double sigma0, int prependOctavesNum = 0, int prependLayersNum = 0, int appendLayersNum = 0);

		//Destructor
		~ScaleSpace();

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

		//Get pixel for sigma at point
		uchar L(Point point, double sigma) const;

		//Get pixel for sigma at row, col
		uchar L(int row, int col, double sigma) const;

		//Get pixel for octave and layer at point
		uchar L(Point point, int octave, int layer) const;

		//Get pixel for octave and layer at row, col
		uchar L(int row, int col, int octave, int layer) const;

		//Get converted coordinate to specific sigma
		CV_labs::Point ConvertCoordinate(int row, int col, double sigma) const;

		//Get converted coordinate to specific sigma
		CV_labs::Point ConvertCoordinate(Point point, double sigma) const;

		//Get converted coordinate to specific octave
		CV_labs::Point ConvertCoordinate(int row, int col, int octave, int layer) const;

		//Get converted coordinate to specific octave
		CV_labs::Point ConvertCoordinate(Point point, int octave, int layer) const;

		//Get restored coordinate to specific sigma
		CV_labs::Point RestoreCoordinate(int row, int col, double sigma) const;

		//Get restored coordinate to specific sigma
		CV_labs::Point RestoreCoordinate(Point point, double sigma) const;

		//Get restored coordinate to specific octave
		CV_labs::Point RestoreCoordinate(int row, int col, int octave, int layer) const;

		//Get restored coordinate to specific octave
		CV_labs::Point RestoreCoordinate(Point point, int octave, int layer) const;

		//Calculate difference of gaussians fpr point
		//double DoG(Point point, int octave, int layer) const;

		//calculate DoG
		double* DoG(int octave, int layer) const;
	protected:
		
		//Get nearest image index to specific sigma
		int getImageIndex(double sigma) const;
	};

	//Get number of octaves
	inline int ScaleSpace::GetOctavesNum() const {

		return this->octavesNum;
	}

	//Get number of layers
	inline int ScaleSpace::GetLayersNum() const {

		return this->layersNum;
	}

	//Get number of images
	inline int ScaleSpace::GetImagesNum() const {

		return this->imagesNum;
	}

	//Get sigma 0
	inline double ScaleSpace::GetSigma0() const {

		return this->sigma0;
	}

	//Get sigma at octave and layer
	inline double ScaleSpace::GetImageSigma(int octave, int layer) const {

		if (octave < -prependOctavesNum || octave >= octavesNum) {
			throw new std::invalid_argument("Try to get image at nonexistent octave");
		}
		if (layer < -prependLayersNum || layer >= layersNum + appendLayersNum) {
			throw new std::invalid_argument("Try to get image at nonexistent layer");
		}

		return std::pow(sigmaInterval, layer) * std::pow(2, octave) * sigma0;
	}
	
	//Get converted coordinate to specific sigma
	inline CV_labs::Point ScaleSpace::ConvertCoordinate(Point point, double sigma) const {

		int index = getImageIndex(sigma);

		int octave = index / realLayersNum;
		double k = 1 / std::pow(2, octave);

		CV_labs::Point result;
		result.x = std::round(point.x * k);
		result.y = std::round(point.y * k);

		return result;
	}

	//Get restored coordinate to specific sigma
	inline CV_labs::Point ScaleSpace::RestoreCoordinate(int row, int col, double sigma) const {

		return RestoreCoordinate({ col, row }, sigma);
	}

	//Get restored coordinate to specific sigma
	inline CV_labs::Point ScaleSpace::RestoreCoordinate(Point point, double sigma) const {

		int index = getImageIndex(sigma);

		int octave = index / realLayersNum;
		double k = std::pow(2, octave);

		CV_labs::Point result;
		result.x = std::round(point.x * k);
		result.y = std::round(point.y * k);

		return result;
	}

	//Get restored coordinate to specific octave
	inline CV_labs::Point ScaleSpace::RestoreCoordinate(int row, int col, int octave, int layer) const
	{
		return RestoreCoordinate({ col, row }, octave, layer);
	}

	//Get restored coordinate to specific octave
	inline CV_labs::Point ScaleSpace::RestoreCoordinate(Point point, int octave, int layer) const {

		if (octave < -prependOctavesNum || octave >= octavesNum) {
			throw new std::invalid_argument("Try to get image at nonexistent octave");
		}
		if (layer < -prependLayersNum || layer >= layersNum + appendLayersNum) {
			throw new std::invalid_argument("Try to get image at nonexistent layer");
		}

		double k = std::pow(2, octave);

		CV_labs::Point result;
		result.x = std::round(point.x * k);
		result.y = std::round(point.y * k);

		return result;
	}

	//Get converted coordinate to specific sigma
	inline CV_labs::Point ScaleSpace::ConvertCoordinate(int row, int col, double sigma) const {
		
		return ConvertCoordinate({ col, row }, sigma);
	}

	//Get nearest image index to specific sigma
	inline int ScaleSpace::getImageIndex(double sigma) const {

		//if sigma lower then first image's sigma, return first image
		if(sigma < sigma0) {
			return prependOctavesNum * realLayersNum + prependLayersNum;
		}

		//if sigma greater, then sigma of last image, return last image index
		if (sigma > GetImageSigma(octavesNum - 1, layersNum - 1)) {
			return (prependOctavesNum + octavesNum - 1)  * realLayersNum + (prependLayersNum + layersNum - 1);
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

		return (prependOctavesNum + octave) * realLayersNum + (prependLayersNum + layer);
	}

	//Get pixel for sigma at row, col
	inline uchar ScaleSpace::L(int row, int col, double sigma) const {

		return L({ col, row }, sigma);
	}

	//Get pixel for sigma at point
	inline uchar ScaleSpace::L(Point point, double sigma) const {

		int index = getImageIndex(sigma);

		int octave = index / realLayersNum;
		double k = 1 / std::pow(2, octave);

		CV_labs::Point result;
		result.x = std::round(point.x * k);
		result.y = std::round(point.y * k);

		return images[index].GetValueAt(result);
	}
}