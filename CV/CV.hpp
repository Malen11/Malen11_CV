#pragma once
#include "Image.hpp"

//определения 

//ядро (для фильтра)
struct Core {
	
	int xk = 0;
	int yk = 0;
	double* data = NULL;
};

/*struct SeparableCore {

	int k = 0;
	double* dataX = NULL;
	double* dataY = NULL;
};*/

//класс функций компьютерного зрения
class CV {
private:
public:
	static const int kInterpolateZero = 201;			//заполнить края 0
	static const int kInterpolateBorder = 202;			//использовать крайний пиксель
	static const int kInterpolateReflection = 203;		//отразить границу
	static const int kInterpolateWrap = 204;			//заворачивать границу

	template <typename T>
	static double* ApplyFilterRaw(int rows, int cols, T* data, Core core, int type);
	template <typename T>
	static T GetVirtualPixel(int row, int col, int interpolateType);	
	static Image Sobel(Image& img);
	template <typename T>
	static double* SobelRaw(int rows, int cols, T* data, int interpolateType = kInterpolateZero);
	static Image ApplyFilter(Image& img, Core core, int type);

	//устаревшее
	template <typename T>
	static T* ExpandImgByZero(int rows, int cols, T* data, int addX, int addY);
};

//реализация темплейтов

template <typename T>
double * CV::ApplyFilterRaw(int rows, int cols, T* data, Core core, int type) {

	double* result = new double[rows*cols];
	//int trows = rows + 2 * core.yk;
	//int tcols = cols + 2 * core.xk;
	//T* temp = NULL;
	int center = core.yk * (2 * core.xk + 1) + core.xk;;

	int pos, tpos;

	/*switch (type) {
	case kInterpolateZero:
	temp = ExpandImgByZero(rows, cols, data, core.xk, core.yk);
	break;
	default:
	break;
	}*/

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			result[pos] = 0;

			for (int u = -core.yk; u <= core.yk; u++) {
				for (int v = -core.xk; v <= core.xk; v++) {

					if ((i + u) > 0 && (i + u) < rows && (j + v) > 0 && (j + v) < cols) {
						result[pos] += data[pos + u * cols + v] * core.data[center + u * (2 * core.yk + 1) + v];
					}
					else {
						result[pos] += GetVirtualPixel<T>(i + u, j + v, type) * core.data[center + u * (2 * core.yk + 1) + v];
					}
				}
			}
		}
	}

	//delete[] temp;

	return result;
}

template<typename T>
T* CV::ExpandImgByZero(int rows, int cols, T* data, int addX, int addY) {

	int trows = rows + 2 * addX;
	int tcols = cols + 2 * addY;
	T* temp = new T[trows*tcols];

	//center
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			temp[(i + addY)*tcols + (j + addX)] = data[i*cols + j];
		}
	}

	//top&bottom + angles
	for (int i = 0; i < addX; i++) {
		for (int j = 0; j < tcols; j++) {

			temp[i*tcols + j] = 0;
			temp[(trows - 1 - i)*tcols + j] = 0;
		}
	}

	//left&right
	for (int i = addX; i < trows - addX; i++) {
		for (int j = 0; j < addY; j++) {

			temp[i*tcols + j] = 0;
			temp[i*tcols + (tcols - 1 - j)] = 0;
		}
	}

	return temp;
}

template<typename T>
T CV::GetVirtualPixel(int row, int col, int interpolateType) {

	switch (interpolateType) {
	case kInterpolateZero:
		return 0;
		break;
	default:
		break;
	}
	return NULL;
}

template <typename T>
double * CV::SobelRaw(int rows, int cols, T* data, int interpolateType) {

	Core Gx1, Gx2;//фильтры для градиента по x(т.к. сепарабельный, разбит на 2)
	Gx1.xk = 1;
	Gx1.yk = 0;
	Gx1.data = new double[3];
	Gx1.data[0] = 1; Gx1.data[1] = 0; Gx1.data[2] = -1;

	Gx2.xk = 0;
	Gx2.yk = 1;
	Gx2.data = new double[3];
	Gx2.data[0] = 1; Gx2.data[1] = 2; Gx2.data[2] = 1;

	Core Gy1, Gy2;//фильтры для градиента по y(т.к. сепарабельный, разбит на 2)
	Gy1.xk = 1;
	Gy1.yk = 0;
	Gy1.data = new double[3];
	Gy1.data[0] = 1; Gy1.data[1] = 2; Gy1.data[2] = 1;


	Gy2.xk = 0;
	Gy2.yk = 1;
	Gy2.data = new double[3];
	Gy2.data[0] = 1; Gy2.data[1] = 0; Gy2.data[2] = -1;

	double *temp, *temp1, *temp2, *result;

	temp = ApplyFilterRaw(rows, cols, data, Gx1, interpolateType);

	temp1 = ApplyFilterRaw(rows, cols, temp, Gx2, interpolateType);

	Image test1(rows, cols, temp1, true);
	imshow("Gx image", test1.GetMat());

	delete[] temp;

	temp = ApplyFilterRaw(rows, cols, data, Gy1, interpolateType);

	temp2 = ApplyFilterRaw(rows, cols, temp, Gy2, interpolateType);

	Image test2(rows, cols, temp2, true);
	imshow("Gy image", test2.GetMat());
	delete[] temp;

	int size = rows * cols;
	result = new double[size];
	for (int i = 0; i < size; i++) {

		result[i] = sqrt(temp1[i] * temp1[i] + temp2[i] * temp2[i]);
	}

	delete[] temp1;
	delete[] temp2;

	return result;
}
