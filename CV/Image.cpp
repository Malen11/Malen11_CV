#include "stdafx.h"
#include "Image.hpp"

using namespace std;
using namespace cv;

//Realization

//Default constuctor
Image::Image() {

	this->data = NULL;
	this->cols = 0;
	this->rows = 0;
}

//constructor from cv::Mat, type mean representation (ntfc,red, green, blue)
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

//copy constructor
Image::Image(const Image& other) {

	this->cols = other.cols;
	this->rows = other.rows;
	this->data = new uchar[cols*rows];

	copy(other.data, &(other.data[rows * cols]), stdext::checked_array_iterator<uchar*>(this->data, rows*cols));
}

Image::Image(int rows, int cols, uchar defValue) {

	int size = cols * rows;

	this->rows = rows;
	this->cols = cols;
	this->data = new uchar[size];

	for (int i = 0; i < size; i++) {
		data[i] = defValue;
	}
}

//create and get cv::Mat
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

/*vector<cv::Point> Image::GetPoints(int Threshold) const {
	vector<Point> points;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (data[i*cols + j] >= Threshold)
				points.push_back(Point(i, j));
		}
	}
	return points;
}*/

Image Image::GetDownsampleImage(int scale) const {

	Image result;
	result.cols = this->cols / scale;
	result.rows = this->rows / scale;
	result.data = new uchar[result.rows*result.cols];

	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			result.data[i*result.cols + j] = this->data[(i*scale)*this->cols + (j*scale)];
		}
	}

	return result;
}

//normalize image data
void Image::NormalizeImage() {

	uchar* temp = this->GetNormalizeDataUC();
	delete[] data;
	data = temp;
}

void Image::DownsampleImage(int scale) {

	Image temp = this->GetDownsampleImage(scale);
	Swap(temp);
}

Image Image::AbsoluteDiff(const Image& img1, const Image&img2) {
	
	if (img1.rows == img2.rows && img1.cols == img2.cols) {

		Image result;
		result.cols = img1.cols;
		result.rows = img1.rows;
		result.data = new uchar[result.rows*result.cols];

		int size = result.rows * result.cols;
		for (int i = 0; i < size; i++) {
			result.data[i] = std::abs(img1.data[i] - img2.data[i])%255;
		}

		return result;
	}
	else {
		throw std::out_of_range("Images dimensions not match");
	}
}

void Image::RotateImage(double angle) {

	int rowMid = this->rows / 2, colMid = this->cols/2;
	int size = rows*cols;
	uchar* dataRotated1 = new uchar[size];
	//uchar* dataRotated2 = new uchar[size];

	for (int i = 0; i < size; i++) {
		dataRotated1[i] = 255;
		//dataRotated2[i] = 255;
	}
	
	int rowNew, colNew;
	int rowOld, colOld;
	
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			colOld = std::cos(angle)*(j - colMid) - std::sin(angle)*(i - rowMid)+colMid;
			rowOld = std::sin(angle)*(j - colMid) + std::cos(angle)*(i - rowMid)+rowMid;

			if (colOld >= 0 && colOld < cols && rowOld >= 0 && rowOld < rows) {

				dataRotated1[i*cols + j] = data[rowOld*cols + colOld];
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

	delete[] data;
	//delete[] dataRotated2;
	data = dataRotated1;
}

//destructor
Image::~Image() {

	if (!this->IsEmpty())
		delete[] data;
}

//print data in console
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

//check, contain image correct data or not
bool Image::IsEmpty() const {

	if (this->cols > 0 && this->rows > 0 && this->data != NULL)
		return false;
	else
		return true;
}

Image Image::InsertImage(Image img, int posX, int posY) const {

	Image resImg(rows,cols, data);


	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {

			if ((i + posY)>0 && (i + posY) < rows && (j + posX)>0 && (j + posX) < cols) {
				resImg.data[(i + posY)*cols + (j + posX)] = img.data[i*img.cols + j];
			}
		}
	}
	return resImg;
}

//swap
void Image::Swap(Image & other) {

	// swap all the members (and base subobject, if applicable) with other
	using std::swap; // because of ADL the compiler will use 
						// custom swap for members if it exists
						// falling back to std::swap
	swap(this->data, other.data);
	swap(this->rows, other.rows);
	swap(this->cols, other.cols);
}

//get number rows
int Image::GetRowsNumber()  const {

	return this->rows;
}

//get number columns
int Image::GetColsNumber() const {

	return this->cols;
}

//get number pixels
int Image::GetSize() const {

	return rows*cols;
}

//get image (copy!) data 
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

//get image normalize data (uchar 0-255)
uchar* Image::GetNormalizeDataUC()  const {

	return LinearNormalization<uchar, uchar>(rows*cols, this->data, 0, 255);
}

//get image normalize data (double 0-1)
double* Image::GetNormalizeDataF()  const {

	/*int max = this->GetMaxValue();
	int min = this->GetMinValue();
	int size = rows * cols;

	double* result = new double[size];

	for (int i = 0; i < size; i++) {

		result[i] = (data[i] - min) / (double)(max - min);
	}*/

	return LinearNormalization<uchar, double>(rows*cols, this->data, 0, 1);
}

//get [row][col] value
uchar Image::GetValueAt(int row, int col)  const {

	return this->data[row*this->cols + col];
}

//get min pixel value
double Image::GetMinValue()  const {
	
	double min = data[0];

	int size = rows * cols;
	for (int i = 0; i < size; i++) {

		if (data[i] < min)
			min = data[i];
	}

	return min;
}

//get max pixel value
double Image::GetMaxValue()  const {

	double max = data[0];

	int size = rows * cols;
	for (int i = 0; i < size; i++) {

		if (data[i] > max)
			max = data[i];
	}

	return max;
}

//set [row][col] value
void Image::SetValueAt(int row, int col, uchar value) {

	this->data[row*this->cols + col] = value;
}

Image& Image::operator=(Image other) {
	
	if(this != &other)
		Swap(other);
	
	return *this;
}


//need modify!
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