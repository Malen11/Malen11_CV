#pragma once
#include <opencv2/opencv.hpp>

//���������

const int kRGB2GRAY_NTSC = 105;
const int kRGB2GRAY_HDTV = 104;
const int kRED = 103;
const int kGREEN = 102;
const int kBLUE = 101;

//�����������

//����� ��� �������� �����������
class Image {
private:
	uchar* data;		//������� ��������
	int rows;			//���������� �����
	int cols;			//���������� ��������

	void Swap(Image& other);
public:
	Image();												//����������� �� ��������� (0 ����� � ��������, ��������� ������� = NULL)
	Image(cv::Mat img, int type = kRGB2GRAY_HDTV);			//����������� �� ������ �������
	Image(const Image& other);								//����������� �����������
	template <typename T>									
	Image(int rows, int cols, T** matrix);					//����������� �� ������ 2-D �������
	template<typename T>									
	Image(int rows, int cols, T * matrix);					//����������� �� ������ 1-D �������

	Image& operator =(Image other);					//���������� ��� �������
	Image& operator -(const Image& other);					//��������� ��� �������
	Image& operator +(const Image& other);					//�������� ��� �������

	//get
	int GetRowsNumber() const;									//�������� ���������� �����
	int GetColsNumber() const;									//�������� ���������� ��������
	uchar* GetData() const;										//�������� ������
	uchar* GetNormalizeDataUC() const;							//�������� ��������������� ������ ���� uchar
	double* GetNormalizeDataF() const;							//�������� ��������������� ������ ���� double
	uchar GetValueAt(int row, int col) const;					//�������� �������� � ������� [row, col]
	uchar GetMinValue() const;									//�������� ����������� ��������
	uchar GetMaxValue() const;									//�������� ������������ ��������
	cv::Mat GetMat() const;										//�������������� ������� ������ � ������� ���� ������

	//set
	void SetValueAt(int row, int col, uchar value);				// ���������� �������� value � ������� [row, col]

	void Print() const;											//����� ������ �� �����
	bool IsEmpty() const;										//���������, ������ �� �����������

	~Image();
};

//���������� template �������
template <typename T>
Image::Image(int rows, int cols, T** matrix) {

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

template<typename T>
Image::Image(int rows, int cols, T * matrix) {

	this->rows = rows;
	this->cols = cols;

	if (cols > 0 && rows > 0 && matrix != NULL) {

		this->data = new uchar[cols*rows];

		int size = rows * cols;
		for (int i = 0; i < size; i++) {

			this->data[i] = matrix[i];
		}
	}
	else {

		std::cout << "Image(int rows, int cols, const uchar ** matrix): Empty array given, using default." << endl;
		this->data = NULL;
		this->cols = 0;
		this->rows = 0;
	}
}
