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
	uchar* data;		//матрица значений
	int rows;			//количество строк
	int cols;			//количество столбцов
public:
	Image();														//конструктор по умолчанию (0 строк и столбцов, указатель матрицы = NULL)
	Image(cv::Mat img, int type = kRGB2GRAY_HDTV);					//конструктор на основе матрицы
	Image(const Image& other);										//конструктор копирования
	template <typename T>									
	Image(int rows, int cols, T** matrix, bool normalize = false);	//конструктор на основе 2-D массива
	template<typename T>									
	Image(int rows, int cols, T * matrix, bool normalize = false);	//конструктор на основе 1-D массива

	Image& operator =(Image other);								//присвоение для ленивых

	//get
	int GetRowsNumber() const;									//получить количество строк
	int GetColsNumber() const;									//получить количество столбцов
	int GetSize() const;										//получить количество pixels в изображении
	uchar* GetData() const;										//получить данные
	uchar* GetNormalizeDataUC() const;							//получить нормализованные данные типа uchar
	double* GetNormalizeDataF() const;							//получить нормализованные данные типа double
	uchar GetValueAt(int row, int col) const;					//получить значение в позиции [row, col]
	double GetMinValue() const;									//получить минимальное значение
	double GetMaxValue() const;									//получить максимальное значение
	cv::Mat GetMat() const;										//преобразовывет матрицу данных к другому типу данных
	Image GetDownsampleImage(int scale) const;

	//set
	void SetValueAt(int row, int col, uchar value);				// установить значение value в позиции [row, col]

	void Print() const;											//вывод данных на экран
	bool IsEmpty() const;										//проверяет, пустое ли изображение
	void Swap(Image& other);
	template<typename srcT, typename dstT>
	static dstT* LinearNormalization(int dataSize, srcT* data, dstT newMin, dstT newMax);	//функция для линейной нормализации
	void NormalizeImage();										//Normalize image to 0-255 uchar
	void DownsampleImage(int scale);
	Image AbsoluteDiff(const Image& img1, const Image&img2);

	~Image();

	//need modify!
	Image& operator -(const Image& other);						//вычитание для ленивых
	Image& operator +(const Image& other);						//сложение для ленивых
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