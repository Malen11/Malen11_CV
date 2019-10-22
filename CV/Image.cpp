#include "stdafx.h"
#include "Image.hpp"

using namespace std;
using namespace cv;
using namespace CV_labs;

#pragma region Constructors & Destructor

//Default constuctor.
//Set data = NULL and rowsNum/colsNum = 0.
Image::Image() {

	this->data = NULL;
	this->colsNum = 0;
	this->rowsNum = 0;
}

//Constructor.
//Get OpenCV Mat as data source and int type as a way, how image must be read.
Image::Image(cv::Mat img, int type) {

	this->rowsNum = img.rows;
	this->colsNum = img.cols;

	this->data = new uchar[rowsNum*colsNum];

	int size = rowsNum * colsNum;
	uchar* matData= img.data;

	switch (type) {
	case kRED:
		for (int i = 0; i < size; i++) {
			this->data[i] = matData[3 * i + 2];
		}

		break;
	case kGREEN:
		for (int i = 0; i < size; i++) {
			this->data[i] = matData[3 * i + 1];
		}

		break;
	case kBLUE:
		for (int i = 0; i < size; i++) {
			this->data[i] = matData[3 * i];
		}

		break;
	case kRGB2GRAY_NTSC:
		for (int i = 0; i < size; i++) {
			this->data[i] = 0.299 * matData[3 * i + 2] + 0.587 * matData[3 * i + 1] + 0.114 * matData[3 * i];
		}

		break;
	default:
		for (int i = 0; i < size; i++) {
			this->data[i] = 0.213 * matData[3 * i + 2] + 0.715 * matData[3 * i + 1] + 0.072 * matData[3 * i];
		}

		break;
	}
}

//Copy constructor.
Image::Image(const Image& other) {

	this->rowsNum = other.rowsNum;
	this->colsNum = other.colsNum;
	this->data = new uchar[rowsNum*colsNum];

	copy(other.data, &(other.data[rowsNum * colsNum]), stdext::checked_array_iterator<uchar*>(this->data, rowsNum * colsNum));
}

//One-color image with specific side sizes' constructor.
Image::Image(int rowsNum, int colsNum, uchar defValue) {

	this->rowsNum = rowsNum;
	this->colsNum = colsNum;

	int size = this->rowsNum * this->colsNum;

	this->data = new uchar[size];

	for (int i = 0; i < size; i++) {

		data[i] = defValue;
	}
}

//Constructor from array.
Image::Image(int rowsNum, int colsNum, uchar* data) {

	this->rowsNum	= rowsNum;
	this->colsNum	= colsNum;
	this->data		= new uchar[rowsNum * colsNum];

	copy(data, &(data[rowsNum * colsNum]), stdext::checked_array_iterator<uchar*>(this->data, rowsNum * colsNum));
}

CV_labs::Image::Image(int rowsNum, int colsNum, double * data) {

	this->rowsNum = rowsNum;
	this->colsNum = colsNum;

	int size = this->rowsNum * this->colsNum;

	this->data = new uchar[size];

	for (int i = 0; i < size; i++) {

		this->data[i] = data[i] * 255;
	}
}

//Destructor.
Image::~Image() {

	if (!this->IsEmpty()) {

		delete[] data;
	}
}

#pragma endregion

#pragma region Operators

//Assignment operator.
Image& Image::operator=(Image other) {

	if (this != &other) {

		Swap(other);
	}

	return *this;
}

/*
Image & Image::operator-(const Image & other) {

	if (!(this->rows == other.rows && this->cols == other.cols)) {

		throw std::out_of_range("Images dimensions not match");
	}

	Image* temp = new Image();

	temp->rows = this->rows;
	temp->cols = this->cols;

	int size = rows * cols;
	for (int i = 0; i < size; i++) {
		temp->data[i] = std::max(0, this->data[i] - other.data[i]);
	}

	return *temp;
}

Image & Image::operator+(const Image & other) {

	if (!(this->rows == other.rows && this->cols == other.cols)) {

		throw std::out_of_range("Images dimensions not match");
	}

	Image* temp = new Image();

	temp->rows = this->rows;
	temp->cols = this->cols;

	int size = rows * cols;
	for (int i = 0; i < size; i++) {
		temp->data[i] = std::max(0, this->data[i] + other.data[i]);
	}

	return *temp;
}
*/

#pragma endregion

#pragma region Geters

//Get uchar* data(copy!).
uchar* Image::GetData() const {

	if (rowsNum * colsNum == 0) {

		return NULL;
	}
	else {

		uchar* temp = new uchar[rowsNum * colsNum];
		copy(data, &(data[rowsNum * colsNum]), stdext::checked_array_iterator<uchar*>(temp, rowsNum * colsNum));

		return temp;
	}
}

double * Image::GetDataD() const
{
	if (this->rowsNum * this->colsNum == 0) {

		return NULL;
	}
	else {

		int size = this->rowsNum * this->colsNum;
		double* temp = new double[this->rowsNum * this->colsNum];

		for(int i = 0; i < size; i++) {

			temp[i] = data[i] / 256.;
		}

		return temp;
	}
}

//Get min pixel value.
uchar Image::GetMinValue()  const {

	int size = this->GetSize();

	if (size > 0) {
		uchar min = data[0];

		for (int i = 1; i < size; i++) {

			if (data[i] < min)
				min = data[i];
		}

		return min;
	}
	else {
		throw std::out_of_range("Try to get value from empty array");
	}
}

//Get max pixel value.
uchar Image::GetMaxValue()  const {

	int size = this->GetSize();

	if (size > 0) {
		uchar max = data[0];

		for (int i = 1; i < size; i++) {

			if (data[i] > max)
				max = data[i];
		}

		return max;
	}
	else {
		throw std::out_of_range("Try to get value from empty array");
	}
}

//Convert image to Mat and return.
Mat Image::GetMat() const {

	if (this->IsEmpty()) {

		return Mat();
	}
	else {

		Mat temp(this->rowsNum, this->colsNum, CV_8UC1);
		copy(data, &(data[this->rowsNum * this->colsNum]), stdext::checked_array_iterator<uchar*>(temp.data, this->rowsNum * this->colsNum));

		return temp;
	}
}

//Get copy of image with lower size
Image CV_labs::Image::GetDownsampledImage(double scale) const{

	if (scale <= 1) {

		throw std::invalid_argument("Scale must be greater then 1");
	}
	if (2 * scale > this->GetColsNumber() || 2 * scale > this->GetRowsNumber()) {

		throw std::invalid_argument("Scale is bigger then half of one or many image dimensions");
	}

	Image result = Image();

	result.rowsNum = (int)round(this->rowsNum / scale);
	result.colsNum = (int)round(this->colsNum / scale);

	result.data = new uchar[result.rowsNum * result.colsNum];

	for (int i = 0; i < result.rowsNum; i++) {
		for (int j = 0; j < result.colsNum; j++) {

			result.data[i * result.colsNum + j] = this->data[(int)round(i * scale) * this->colsNum + (int)round(j * scale)];
		}
	}

	return result;
}

//Get copy of image with highter size (probably buged, image move to left)
Image CV_labs::Image::GetUpsampledImage(double scale) const {

	if (scale <= 1) {

		throw std::invalid_argument("Scale must be greater then 1");
	}

	Image result = Image();

	result.rowsNum = (int)round(this->rowsNum * scale);
	result.colsNum = (int)round(this->colsNum * scale);

	result.data = new uchar[result.rowsNum * result.colsNum];

	for (int i = 0; i < result.rowsNum; i++) {
		for (int j = 0; j < result.colsNum; j++) {

			result.data[i * result.colsNum + j] = this->data[std::min((int)round(i / scale), this->rowsNum - 1) * this->colsNum + std::min((int)round(j / scale), this->colsNum - 1)];
		}
	}

	return result;
}

//Get copy of image with different size. (not recomended)
Image Image::GetScaledImage(double scale) const {
	
	if (scale <= 0) {

		throw std::invalid_argument("Scale must be greater then 0");
	}

	Image result = Image();

	result.rowsNum = (int)round(this->rowsNum * scale);
	result.colsNum = (int)round(this->colsNum * scale);

	result.data = new uchar[result.rowsNum * result.colsNum];

	for (int i = 0; i < result.rowsNum; i++) {
		for (int j = 0; j < result.colsNum; j++) {

			result.data[i * result.colsNum + j] = this->data[(int)round(i / scale) * this->colsNum + (int)round(j / scale)];
		}
	}

	return result;
}

//Get normalized copy of image data.
uchar* Image::GetNormalizedImageDataUC() const {

	int dataSize = this->GetSize();
	double min, max, k;

	if (data != NULL) {

		min = max = data[0];

		for (int i = 0; i < dataSize; i++) {

			if (data[i] > max) {
				max = data[i];
			}

			if (data[i] < min) {
				min = data[i];
			}
		}

		k = 255 / (max - min);

		uchar* result = new uchar[dataSize];

		for (int i = 0; i < dataSize; i++) {

			result[i] = (uchar)(data[i] - min) * k;
		}

		return result;
	}
	else {

		throw std::invalid_argument("Image must not be empty");
	}
}

//Get normalized copy of image data.
double* Image::GetNormalizedImageDataD() const {

	int dataSize = this->GetSize();
	double min, max, k;

	if (data == NULL) {
		throw std::invalid_argument("Image must not be empty");
	}
		
	min = max = data[0];

	for (int i = 0; i < dataSize; i++) {

		if (data[i] > max) {
			max = data[i];
		}

		if (data[i] < min) {
			min = data[i];
		}
	}

	k = 1 / (max - min);

	double* result = new double[dataSize];

	for (int i = 0; i < dataSize; i++) {

		result[i] = (double)(data[i] - min) * k;
	}

	return result;
}

//Get normalized copy of image.
Image Image::GetNormalizedImage() const {

	return Image(this->rowsNum, this->colsNum, this->GetNormalizedImageDataUC());
}

//Get rotated copy of image.
Image Image::GetRotatedImage(double angle) const {

	int rowMid = this->rowsNum / 2;
	int colMid = this->colsNum / 2;
	int size = this->rowsNum * this->colsNum;

	uchar* dataRotated1 = new uchar[size];
	//uchar* dataRotated2 = new uchar[size];

	for (int i = 0; i < size; i++) {

		dataRotated1[i] = 255;
		//dataRotated2[i] = 255;
	}

	//int rowNew, colNew;
	int rowOld, colOld;

	for (int i = 0; i < this->rowsNum; i++) {
		for (int j = 0; j < this->colsNum; j++) {

			colOld = std::cos(angle)*(j - colMid) - std::sin(angle)*(i - rowMid) + colMid;
			rowOld = std::sin(angle)*(j - colMid) + std::cos(angle)*(i - rowMid) + rowMid;

			if (colOld >= 0 && colOld < this->colsNum && rowOld >= 0 && rowOld < this->rowsNum) {

				dataRotated1[i * this->colsNum + j] = data[rowOld * this->colsNum + colOld];
			}
		}
	}

	/*for (int i = -rowMid; i <= rowMid; i++) {
	for (int j = -colMid; j <= colMid; j++) {

	colNew = colMid+(std::cos(angle)*j - std::sin(angle)*i);
	rowNew = rowMid+(std::sin(angle)*j + std::cos(angle)*i);

	if (colNew >= 0 && colNew < cols && rowNew >= 0 && rowNew < rows) {

	dataRotated2[rowNew*cols + colNew] = data[(rowMid+i)*cols + (colMid+j)];
	}
	}
	}

	int pixel, n;
	for (int i = 0; i < rows; i++) {
	for (int j = 0; j < cols; j++) {

	pixel = 0;
	n = 0;

	if (dataRotated2[i*cols + j] == 255) {

	for (int u = -1; u <= 1; u++) {
	for (int v = -1; v <= 1; v++) {
	if (u + i >= 0 && u + i < rows && v + j>=0 && v + j < cols) {
	pixel += dataRotated2[(i + u)*cols + (j + v)];
	++n;
	}
	}
	}

	dataRotated2[i*cols + j] = pixel/n;
	}
	}
	}*/

	return Image(this->rowsNum, this->colsNum, dataRotated1);
}

#pragma endregion

#pragma region Methods

//Print every pixel value in console.
void Image::Print() const {

	int step;

	cout << "Rows number: " << this->rowsNum << endl;
	cout << "Columns number: " << this->colsNum << endl;
	cout << '[';

	for (int i = 0; i < rowsNum; i++) {

		step = i * colsNum;
		for (int j = 0; j < colsNum; j++) {

			cout << (int)data[step + j] << ", ";
		}
		if (i + 1 == rowsNum)
			cout << ']';
		else
			cout << endl;
	}

	cout << endl;
}

//Change image size. (not recomended)
void Image::ScaleImage(double scale) {

	(*this) = this->GetScaledImage(scale);
}

void CV_labs::Image::DownsampleImage(double scale) {
	(*this) = this->GetDownsampledImage(scale);
}

void CV_labs::Image::UpsampleImage(double scale) {
	(*this) = this->GetDownsampledImage(scale);
}

//Normalize image.
void Image::NormalizeImage() {

	(*this) = this->GetNormalizedImage();
}

//Rotate image by *** angle.
void Image::RotateImage(double angle) {

	(*this) = this->GetRotatedImage(angle);
}

//Calculate absolute differense of two images and return result as image.
Image Image::CalcAbsoluteDiff(const Image& img2) const {

	if (this->rowsNum == img2.rowsNum && this->colsNum == img2.colsNum) {

		Image result;

		result.rowsNum = img2.rowsNum;
		result.colsNum = img2.colsNum;
		result.data = new uchar[result.rowsNum * result.colsNum];

		int size = result.rowsNum * result.colsNum;
		for (int i = 0; i < size; i++) {

			result.data[i] = std::abs(this->data[i] - img2.data[i]);
		}

		return result;
	}
	else {

		throw std::out_of_range("Images dimensions not match");
	}
}

//Insert other image into this image.
void Image::InsertImage(const Image& other, int posX, int posY) const {

	if (this->rowsNum > posX && this->colsNum > posY) {

		for (int i = 0; i < other.rowsNum; i++) {
			for (int j = 0; j < other.colsNum; j++) {

				if (((i + posY) > 0 && (i + posY) <  this->rowsNum) && ((j + posX) > 0 && (j + posX) < this->colsNum)) {

					this->data[(i + posY) * this->colsNum + (j + posX)] = other.data[i * other.colsNum + j];
				}
			}
		}
	}
	else {

		throw std::invalid_argument("Cordinate [posX, posY] must be on image");
	}
}

Image CV_labs::Image::CutOutImage(int row0, int col0, int row1, int col1) const
{
	if (row0 < 0 || row0 > row1 || row1 > this->rowsNum || col0 < 0 || col0 > col1 || col1 > this->colsNum) {

		throw std::invalid_argument("Values must be in range [0 < row0 < row1 < rowsNum], [0 < col0 < col1 < colNum]");
	}

	int rowsNum = row1 - row0 + 1;
	int colsNum = col1 - col0 + 1;

	uchar* data = new uchar[rowsNum * colsNum];

	for (int rowInd = 0; rowInd < rowsNum; rowInd++) {
		for (int colInd = 0; colInd < colsNum; colInd++) {
			data[rowInd * colsNum + colInd] = this->data[(row0 + rowInd) * this->colsNum + (col0 + colInd)];
		}
	}

	return Image(rowsNum, colsNum, data);
}

//Swap this image and other image.
void Image::Swap(Image & other) {

	// swap all the members (and base subobject, if applicable) with other
	using std::swap; // because of ADL the compiler will use 
					 // custom swap for members if it exists
					 // falling back to std::swap
	swap(this->data, other.data);
	swap(this->rowsNum, other.rowsNum);
	swap(this->colsNum, other.colsNum);
}

#pragma endregion