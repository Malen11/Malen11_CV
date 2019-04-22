#pragma once
#include <opencv2/opencv.hpp>

//���������

const int RGB2GRAY = 100;
const int RED = 103;
const int GREEN = 102;
const int BLUE = 101;

//�����������

//����� ��� �������� �����������
class Image {
private:
	uchar* data;	//������� ��������
	int rows;			//���������� �����
	int cols;			//���������� ��������
public:
	Image();												//����������� �� ��������� (0 ����� � ��������, ��������� ������� = NULL)
	Image(cv::Mat img, int type = RGB2GRAY);				//����������� �� ������ �������
	Image(int rows, int cols, uchar** matrix);				//����������� �� ������ �������
	Image(const Image& other);								//����������� �����������

	Image& operator =(const Image& other);					//���������� ��� �������
	Image& operator -(const Image& other);					//��������� ��� �������
	Image& operator +(const Image& other);					//�������� ��� �������

	//get
	int GetRowsNumber();									//�������� ���������� �����
	int GetColsNumber();									//�������� ���������� ��������
	uchar* GetData();										//�������� ������
	uchar GetValueAt(int row, int col);						//�������� �������� � ������� [row, col]
	uchar GetMinValue();									//�������� ����������� ��������
	uchar GetMaxValue();									//�������� ������������ ��������
	cv::Mat GetMat();										//�������������� ������� ������ � ������� ���� ������

	//set
	void SetValueAt(int row, int col, uchar value);			// ���������� �������� value � ������� [row, col]

	void Print();											//����� ������ �� �����
	bool IsEmpty();											//���������, ������ �� �����������

	~Image();
};