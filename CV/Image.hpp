#pragma once
#include <opencv2/opencv.hpp>

//Константы

const int RGB2GRAY = 100;
const int RED = 103;
const int GREEN = 102;
const int BLUE = 101;

//Определение

//Класс для хранения изображений
class Image {
private:
	uchar* data;	//матрица значений
	int rows;			//количество строк
	int cols;			//количество столбцов
public:
	Image();												//конструктор по умолчанию (0 строк и столбцов, указатель матрицы = NULL)
	Image(cv::Mat img, int type = RGB2GRAY);				//конструктор на основе матрицы
	Image(int rows, int cols, uchar** matrix);				//конструктор на основе массива
	Image(const Image& other);								//конструктор копирования

	Image& operator =(const Image& other);					//присвоение для ленивых
	Image& operator -(const Image& other);					//вычитание для ленивых
	Image& operator +(const Image& other);					//сложение для ленивых

	//get
	int GetRowsNumber();									//получить количество строк
	int GetColsNumber();									//получить количество столбцов
	uchar* GetData();										//получить данные
	uchar GetValueAt(int row, int col);						//получить значение в позиции [row, col]
	uchar GetMinValue();									//получить минимальное значение
	uchar GetMaxValue();									//получить максимальное значение
	cv::Mat GetMat();										//преобразовывет матрицу данных к другому типу данных

	//set
	void SetValueAt(int row, int col, uchar value);			// установить значение value в позиции [row, col]

	void Print();											//вывод данных на экран
	bool IsEmpty();											//проверяет, пустое ли изображение

	~Image();
};