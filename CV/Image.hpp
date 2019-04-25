#pragma once
#include <opencv2/opencv.hpp>

//Constant

const int kRGB2GRAY_NTSC = 105;
const int kRGB2GRAY_HDTV = 104;
const int kRED = 103;
const int kGREEN = 102;
const int kBLUE = 101;

//Definition

//for store, get info and modify grayscale image data
class Image {
private:
	uchar* data;		//������� ��������
	int rows;			//���������� �����
	int cols;			//���������� ��������
public:
	Image();														//����������� �� ��������� (0 ����� � ��������, ��������� ������� = NULL)
	Image(cv::Mat img, int type = kRGB2GRAY_HDTV);					//����������� �� ������ �������
	Image(const Image& other);										//����������� �����������
	template <typename T>									
	Image(int rows, int cols, T** matrix, bool normalize = false);	//����������� �� ������ 2-D �������
	template<typename T>									
	Image(int rows, int cols, T * matrix, bool normalize = false);	//����������� �� ������ 1-D �������

	Image& operator =(Image other);								//���������� ��� �������

	//get
	int GetRowsNumber() const;									//�������� ���������� �����
	int GetColsNumber() const;									//�������� ���������� ��������
	int GetSize() const;										//�������� ���������� pixels � �����������
	uchar* GetData() const;										//�������� ������
	uchar* GetNormalizeDataUC() const;							//�������� ��������������� ������ ���� uchar
	double* GetNormalizeDataF() const;							//�������� ��������������� ������ ���� double
	uchar GetValueAt(int row, int col) const;					//�������� �������� � ������� [row, col]
	double GetMinValue() const;									//�������� ����������� ��������
	double GetMaxValue() const;									//�������� ������������ ��������
	cv::Mat GetMat() const;										//�������������� ������� ������ � ������� ���� ������
	Image GetDownsampleImage(int scale) const;

	//set
	void SetValueAt(int row, int col, uchar value);				// ���������� �������� value � ������� [row, col]

	void Print() const;											//����� ������ �� �����
	bool IsEmpty() const;										//���������, ������ �� �����������
	void Swap(Image& other);
	template<typename srcT, typename dstT>
	static dstT* LinearNormalization(int dataSize, srcT* data, dstT newMin, dstT newMax);	//������� ��� �������� ������������
	void NormalizeImage();										//Normalize image to 0-255 uchar
	void DownsampleImage(int scale);
	Image AbsoluteDiff(const Image& img1, const Image&img2);

	~Image();

	//need modify!
	Image& operator -(const Image& other);						//��������� ��� �������
	Image& operator +(const Image& other);						//�������� ��� �������
};

//template function realization
template <typename T>
Image::Image(int rows, int cols, T** matrix, bool normalize) {

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

		if (normalize) {
			uchar* temp = LinearNormalization<T, uchar>(rows*cols, this->data, 0, 255);
			delete[] this->data;

			this->data = temp;
		}
	}
	else {

		std::cout << "Image(int rows, int cols, const uchar ** matrix): Empty array given, using default." << endl;
		this->data = NULL;
		this->cols = 0;
		this->rows = 0;
	}
}

template<typename T>
Image::Image(int rows, int cols, T * matrix, bool normalize) {

	this->rows = rows;
	this->cols = cols;

	if (cols > 0 && rows > 0 && matrix != NULL) {


		this->data = new uchar[cols*rows];
		int size = rows * cols;

		if (normalize) {

			this->data = LinearNormalization<T, uchar>(size, matrix, 0, 255);
		}
		else {
			for (int i = 0; i < size; i++) {

				this->data[i] = matrix[i];
			}
		}
	}
	else {

		std::cout << "Image(int rows, int cols, const uchar * matrix): Empty array given, using default." << endl;
		this->data = NULL;
		this->cols = 0;
		this->rows = 0;
	}
}

template<typename srcT, typename dstT>
static dstT * Image::LinearNormalization(int dataSize, srcT* data, dstT newMin, dstT newMax) {

	if (dataSize != 0 && data != NULL) {

		double min = data[0], max = data[0], k;

		for (int i = 0; i < dataSize; i++) {
			if (data[i] > max) {
				max = data[i];
			}
			
			if (data[i] < min) {
				min = data[i];
			}
		}

		k = (newMax - newMin) / (max - min);

		dstT* result = new dstT[dataSize];

		for (int i = 0; i < dataSize; i++) {

			result[i] = (dstT)((data[i] - min)*k + newMin);
		}

		return result;
	}
	else {
		return NULL;
	}
}