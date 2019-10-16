#pragma once
#include "stdafx.h"

namespace CV_labs {

#pragma region Constants

	const int kRGB2GRAY_NTSC = 105;
	const int kRGB2GRAY_HDTV = 104;
	const int kRED = 103;
	const int kGREEN = 102;
	const int kBLUE = 101;

#pragma endregion

	//Class for storing image
	class Image {

	private:
		uchar * data;		//pixel value's matrix
		int rowsNum;		//number of rows
		int colsNum;		//number of columns

	public:

#pragma region Constructors & Destructor

		//Default constuctor.
		//Set data = NULL and rowsNum/colsNum = 0.
		Image();

		//Constructor.
		//Get OpenCV Mat as data source and int type as a way, how image must be read.
		Image(cv::Mat img, int type = kRGB2GRAY_HDTV);

		//Copy constructor.
		Image(const Image& other);

		//Constructor for one-color image with specific side sizes.
		Image(int rowsNum, int colsNum, uchar defValue);

		//Constructor from array.
		Image(int rowsNum, int colsNum, uchar* data);

		Image(int rowsNum, int colsNum, double* data);

		//Destructor.
		~Image();

#pragma endregion

#pragma region Operators

		//Assignment operator.
		Image& operator =(Image other);

		/*Возможно когда-нибудь доделаю
		Image& operator -(const Image& other);						//вычитание для ленивых
		Image& operator +(const Image& other);						//сложение для ленивых
		*/

#pragma endregion

#pragma region Geters

		//Get uchar* data(copy!).
		uchar* GetData() const;

		//Get double* data.
		double* GetDataD() const;

		//Get rows number.
		int GetRowsNumber() const;

		//Get columns number.
		int GetColsNumber() const;

		//Get rows number * columns number.
		int GetSize() const;

		//Get pixel value at [row][col].
		uchar GetValueAt(int row, int col) const;

		//Get min pixel value.
		uchar GetMinValue() const;

		//Get max pixel value.
		uchar GetMaxValue() const;

		//Convert image to Mat and return.
		cv::Mat GetMat() const;

		//Get copy of image with different size. (not recomended)
		Image GetScaledImage(double scale) const;

		//Get normalized copy of image data.
		uchar* GetNormalizedImageDataUC() const;

		//Get normalized copy of image data.
		double* GetNormalizedImageDataD() const;

		//Get normalized copy of image.
		Image GetNormalizedImage() const;

		//Get rotated copy of image.
		Image GetRotatedImage(double angle) const;

#pragma endregion

#pragma region Seters

		//Set pixel value by [row][col].
		void SetValueAt(int row, int col, uchar value);

#pragma endregion

#pragma region Methods

		//Print every pixel value in console.
		void Print() const;

		//Check is image set or empty.
		bool IsEmpty() const;

		//Change image size. (not recomended)
		void ScaleImage(double scale);

		//Normalize image.
		void NormalizeImage();

		//Rotate image by *** angle.
		void RotateImage(double angle);

		//Calculate absolute differense of two images and return result as image.
		Image CalcAbsoluteDiff(const Image&img2) const;

		//Insert other image into this image.
		void InsertImage(const Image &other, int posX, int posY) const;

		//cut out image's part
		Image CutOutImage(int row0, int col0, int row1, int col1) const;

		//Swap this image and other image.
		void Swap(Image& other);

#pragma endregion

	};

	//Get rows number.
	inline int Image::GetRowsNumber()  const {

		return this->rowsNum;
	}

	//Get columns number.
	inline int Image::GetColsNumber() const {

		return this->colsNum;
	}

	//Get rows number * columns number.
	inline int Image::GetSize() const {

		return rowsNum * colsNum;
	}

	//Get pixel value at [row][col].
	inline uchar Image::GetValueAt(int row, int col)  const {

		return this->data[row * colsNum + col];
	}

	//Set pixel value by [row][col].
	inline void Image::SetValueAt(int row, int col, uchar value) {

		this->data[row * colsNum + col] = value;
	}

	//Check is image set or empty.
	inline bool Image::IsEmpty() const {

		return data == NULL;
	}

}