#include "stdafx.h"
#include "Image.hpp"

using namespace std;
using namespace cv;

//Реализация

Image::Image() {

	this->data = NULL;
	this->cols = 0;
	this->rows = 0;
}

Image::Image(cv::Mat img, int type) {

	this->rows = img.rows;
	this->cols = img.cols;

	this->data = new uchar[cols*rows];

	int size = cols * rows;
	uchar* img_mat= img.data;

	switch (type) {
	case kRED:
		for (int i = 0; i < size; i++) {
			this->data[i] = img_mat[3 * i + 2];
		}

		break;
	case kGREEN:
		for (int i = 0; i < size; i++) {
			this->data[i] = img_mat[3 * i + 1];
		}

		break;
	case kBLUE:
		for (int i = 0; i < size; i++) {
			this->data[i] = img_mat[3 * i];
		}

		break;
	case kRGB2GRAY_NTSC:
		for (int i = 0; i < size; i++) {
			this->data[i] = 0.299*img_mat[3 * i + 2] + 0.587*img_mat[3 * i + 1] + 0.114*img_mat[3 * i];
		}

		break;
	default:
		for (int i = 0; i < size; i++) {
			this->data[i] = 0.213*img_mat[3 * i + 2] + 0.715*img_mat[3 * i + 1] + 0.072*img_mat[3 * i];
		}

		break;
	}
}

Image::Image(const Image& other) {

	this->cols = other.cols;
	this->rows = other.rows;
	this->data = new uchar[cols*rows];

	copy(other.data, &(other.data[rows * cols]), stdext::checked_array_iterator<uchar*>(this->data, rows*cols));
}

Mat Image::GetMat() const {

	if (this->IsEmpty()) {

		return Mat();
	}
	else {

		Mat temp(this->rows, this->cols, CV_8UC1);
		copy(data, &(data[rows * cols]), stdext::checked_array_iterator<uchar*>(temp.data, rows*cols));

		return temp;
	}
}

Image::~Image() {

	if (!this->IsEmpty())
		delete[] data;
}

void Image::Print() const {

	int step;

	cout << "Rows number: " << this->rows << endl;
	cout << "Columns number: " << this->cols << endl;
	cout << '[';

	for (int i = 0; i < rows; i++) {

		step = i * cols;
		for (int j = 0; j < cols; j++) {

			cout << (int)data[step + j] << ", ";
		}
		if (i + 1 == rows)
			cout << ']';
		else
			cout << endl;
	}

	cout << endl;
}

bool Image::IsEmpty() const {

	if (this->cols > 0 && this->rows > 0 && this->data != NULL)
		return false;
	else
		return true;
}

void Image::Swap(Image & other) {

	// swap all the members (and base subobject, if applicable) with other
	using std::swap; // because of ADL the compiler will use 
						// custom swap for members if it exists
						// falling back to std::swap
	swap(this->data, other.data);
	swap(this->rows, other.rows);
	swap(this->cols, other.cols);
}

int Image::GetRowsNumber()  const {

	return this->rows;
}

int Image::GetColsNumber() const {

	return this->cols;
}

uchar* Image::GetData() const {

	if (cols*rows == 0) {
	
		return NULL;
	}
	else {

		uchar* temp = new uchar[rows*cols];
		copy(data, &(data[rows * cols]), stdext::checked_array_iterator<uchar*>(temp, rows*cols));

		return temp;
	}
}

uchar* Image::GetNormalizeDataUC()  const {

	int max = this->GetMaxValue();
	int min = this->GetMinValue();
	int size = rows * cols;
	
	uchar* result = new uchar[size];

	for (int i = 0; i < size; i++) {

		result[i] = (255 * (data[i] - min)) / (max - min);
	}

	return result;
}

double* Image::GetNormalizeDataF()  const {

	int max = this->GetMaxValue();
	int min = this->GetMinValue();
	int size = rows * cols;

	double* result = new double[size];

	for (int i = 0; i < size; i++) {

		result[i] = (data[i] - min) / (double)(max - min);
	}

	return result;
}

uchar Image::GetValueAt(int row, int col)  const {

	return this->data[row*this->cols + col];
}

uchar Image::GetMinValue()  const {
	
	uchar min = 255;

	int size = rows * cols;
	for (int i = 0; i < size; i++) {

		if (data[i] < min)
			min = data[i];
	}

	return min;
}

uchar Image::GetMaxValue()  const {

	uchar max = 0;

	int size = rows * cols;
	for (int i = 0; i < size; i++) {

		if (data[i] > max)
			max = data[i];
	}

	return max;
}

void Image::SetValueAt(int row, int col, uchar value) {

	this->data[row*this->cols + col] = value;
}

Image& Image::operator=(Image other) {
	
	Swap(other);
	return *this;
}

//вычитает одно изображение из другого
//значения матрицы ограничены и приводятся к промежутку [0, 255]
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

//прибавляет одно изображение к другому
//значения матрицы ограничены и приводятся к промежутку [0, 255]
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