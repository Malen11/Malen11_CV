#pragma once
#include <opencv2/opencv.hpp>

//Константы

const int kRGB2GRAY_NTSC = 105;
const int kRGB2GRAY_HDTV = 104;
const int kRED = 103;
const int kGREEN = 102;
const int kBLUE = 101;

//Определение

//Класс для хранения изображений
class Image {
private:
	uchar* data;		//матрица значений
	int rows;			//количество строк
	int cols;			//количество столбцов

	void Swap(Image& other);
public:
	Image();												//конструктор по умолчанию (0 строк и столбцов, указатель матрицы = NULL)
	Image(cv::Mat img, int type = kRGB2GRAY_HDTV);			//конструктор на основе матрицы
	Image(const Image& other);								//конструктор копирования
	template <typename T>									
	Image(int rows, int cols, T** matrix);					//конструктор на основе 2-D массива
	template<typename T>									
	Image(int rows, int cols, T * matrix);					//конструктор на основе 1-D массива

	Image& operator =(Image other);					//присвоение для ленивых
	Image& operator -(const Image& other);					//вычитание для ленивых
	Image& operator +(const Image& other);					//сложение для ленивых

	//get
	int GetRowsNumber() const;									//получить количество строк
	int GetColsNumber() const;									//получить количество столбцов
	uchar* GetData() const;										//получить данные
	uchar* GetNormalizeDataUC() const;							//получить нормализованные данные типа uchar
	double* GetNormalizeDataF() const;							//получить нормализованные данные типа double
	uchar GetValueAt(int row, int col) const;					//получить значение в позиции [row, col]
	uchar GetMinValue() const;									//получить минимальное значение
	uchar GetMaxValue() const;									//получить максимальное значение
	cv::Mat GetMat() const;										//преобразовывет матрицу данных к другому типу данных

	//set
	void SetValueAt(int row, int col, uchar value);				// установить значение value в позиции [row, col]

	void Print() const;											//вывод данных на экран
	bool IsEmpty() const;										//проверяет, пустое ли изображение

	~Image();
};

//реализация template функций
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
