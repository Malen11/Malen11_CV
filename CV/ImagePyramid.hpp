#pragma once
#include "stdafx.h"

//namespace CV_labs {
//
//	//Parameters for ImagePyramid class object
//	class ImagePyramidParams {
//
//	public:
//		bool setDefaultOctavesNum = false;
//		int octavesNum = 1;
//		int layersNum = 1;
//		double sigmaA = 0;
//		double sigma0 = 1;
//		int layersAdditionalNum = 0;
//
//	private:
//		std::string error;
//
//	public:
//
//		//Check is values valid
//		bool IsValid() {
//
//			if (octavesNum <= 0 && setDefaultOctavesNum == false) {
//				error = "octavesNum must be >= 0";
//			}
//			else if (layersNum <= 0) {
//				error = "layersNum must be >= 0";
//			}
//			else if (sigma0 <= 0) {
//				error = "sigma0 must be > 0";
//			}
//			else if (sigmaA < 0) {
//				error = "sigmaA must be >= 0";
//			}
//			else if (layersAdditionalNum < 0) {
//				error = "layersAdditionalNum must be >= 0";
//			}
//			else {
//				return true;
//			}
//		}
//
//		//Get error
//		string GetError() {
//
//			IsValid();
//
//			return error;
//		}
//	};
//
//
//	//Class for storing collection of images in images pyramid's format
//	class ImagePyramid {
//
//	private:
//		int octavesNum;				//Number of octaves
//		int layersNum;				//Number of layers
//		int layersAdditionalNum;	//Number of cross layers (not count in layersNum)
//		double sigmaA;				//Original sigma
//		double sigma0;				//Starter sigma
//
//		double sigmaInterval;		//Sigma interval between layers
//
//		Image* images;				//Images
//		int imagesNum;				//Number of Images
//
//		//public:
//		//double* sigmaArray;
//
//	public:
//
//#pragma region Constructors & Destructor
//
//		//Default constructor.
//		//Set all data as 0/null.
//		ImagePyramid();
//
//		//Basic constructor.
//		//Using ImagePyramidParams for set property values.
//		ImagePyramid(Image img, ImagePyramidParams params);
//
//		ImagePyramid(Image img, int octavesNum, int layersNum, double sigmaA, double sigma0, int layersAdditionalNum = 0);
//		~ImagePyramid();
//
//#pragma endregion
//
//
//		Image GetImage(int octave, int layer) const;
//		Image GetNearestImage(double sigma);
//		int GetOctavesNum() const;
//		int GetLayersNum() const;
//		int GetLayersAdditionalNum() const;
//		int GetImagesNum() const;
//		double GetSigma0() const;
//		double GetSigmaInterval() const;
//
//		uchar L(int x, int y, double sigma);
//
//		int FindNearestSigma(double sigma);
//		//Image L(double sigmaArray);
//	};
//}