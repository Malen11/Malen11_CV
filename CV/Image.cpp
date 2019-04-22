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

	int step;

	switch (type) {
	case RED:
		for (int i = 0; i < rows; i++) {

			step = i * cols;

			for (int j = 0; j < cols; j++) {

				this->data[step + j] = img.at<cv::Vec3b>(cv::Point(j, i))[2];
			}
		}

		break;
	case GREEN:
		for (int i = 0; i < rows; i++) {

			step = i * cols;

			for (int j = 0; j < cols; j++) {

				this->data[step + j] = img.at<cv::Vec3b>(cv::Point(j, i))[1];
			}
		}

		break;
	case BLUE:
		for (int i = 0; i < rows; i++) {

			step = i * cols;

			for (int j = 0; j < cols; j++) {

				this->data[step + j] = img.at<cv::Vec3b>(cv::Point(j, i))[0];
			}
		}

		break;

	default:

		cv::Vec3b pixel;

		for (int i = 0; i < rows; i++) {

			step = i * cols;

			for (int j = 0; j < cols; j++) {

				this->data[step + j] = 0.213*img.at<cv::Vec3b>(cv::Point(j, i))[2] + 0.715*img.at<cv::Vec3b>(cv::Point(j, i))[1] + 0.072*img.at<cv::Vec3b>(cv::Point(j, i))[0];
			}
		}

		break;
	}
}

Image::Image(int rows, int cols, uchar** matrix) {

	this->rows = rows;
	this->cols = cols;

	if (cols > 0 && rows > 0 && matrix != NULL) {

		this->data = new uchar[cols*rows];

		int step;
		for (int i = 0; i < rows; i++) {

			step = i * cols;
			for (int j = 0; j < cols; j++) {

				this->data[step + j] = matrix[i][j];
			}
		}
	}
	else {

		std::cout << "Image(int rows, int cols, const uchar ** matrix): Empty array given, using default." << endl;
		this->data = NULL;
		this->cols = 0;
		this->rows = 0;
	}
}

Image::Image(const Image& other) {

	this->cols = other.cols;
	this->rows = other.rows;
	this->data = new uchar[cols*rows];

	copy(other.data, &(other.data[rows * cols]), stdext::checked_array_iterator<uchar*>(this->data, rows*cols));
}

Mat Image::GetMat() {

	if (this->IsEmpty()) {

		return Mat();
	}
	else {

		Mat temp(this->rows, this->cols, CV_8UC1);
		copy(data, &(data[rows * cols]), stdext::checked_array_iterator<uchar*>(temp.data,rows*cols));

		return temp;
	}
}

Image::~Image() {

	if (!this->IsEmpty())
		delete[] data;
}

void Image::Print() {

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

inline bool Image::IsEmpty() {

	if (this->cols > 0 && this->rows > 0 && this->data != NULL)
		return false;
	else
		return true;
}

inline int Image::GetRowsNumber() {

	return this->rows;
}

inline int Image::GetColsNumber() {

	return this->cols;
}

inline uchar Image::GetValueAt(int row, int col) {

	return this->data[row*this->cols + col];
}

inline void Image::SetValueAt(int row, int col, uchar value) {

	this->data[row*this->cols + col] = value;
}

Image& Image::operator=(const Image& other) {

	if (!this->IsEmpty())
		delete[] this->data;

	this->cols = other.cols;
	this->rows = other.rows;
	this->data = new uchar[cols*rows];

	copy(other.data, &(other.data[rows * cols]), stdext::checked_array_iterator<uchar*>(this->data, rows*cols));

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