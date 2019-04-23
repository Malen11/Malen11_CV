#include "stdafx.h"
#include "CV.hpp"

using namespace std;

template <typename T>
double * CV::ApplyFilterRaw(int rows, int cols, T* data, Core core, int type) {

	double* result = new double[rows*cols];
	int trows = rows + 2 * core.k;
	int tcols = cols + 2 * core.k;
	T* temp = NULL;
	int center = core.k * (2 * core.k + 1) + core.k;;
	
	int pos, tpos;
	
	switch (type) {
	case kZeroBorder:
		temp = ExpandImgByZero(rows, cols, data, core.k);
		break;
	default:
		break;
	}

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;
			tpos = (i + core.k) * tcols + (j + core.k);

			result[pos] = 0;

			for (int u = -core.k; u <= core.k; u++) {
				for (int v = -core.k; v <= core.k; v++) {

					result[pos] += temp[tpos + u * tcols + v] * core.data[center + u * (2*core.k+1) + v];
				}
			}
		}
	}

	delete[] temp;

	return result;
}

template<typename T>
T* CV::ExpandImgByZero(int rows, int cols, T* data, int numElements) {

	int trows = rows + 2 * numElements;
	int tcols = cols + 2 * numElements;
	T* temp = new T[trows*tcols];

	//center
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			temp[(i + numElements)*tcols + (j + numElements)] = data[i*cols + j];
		}
	}
	
	//top&bottom + angles
	for (int i = 0; i < numElements; i++) {
		for (int j = 0; j < tcols; j++) {

			temp[i*tcols + j] = 0;
			temp[(trows - 1 - i)*tcols + j] = 0;
		}
	}

	//left&right
	for (int i = numElements; i < trows - numElements; i++) {
		for (int j = 0; j < numElements; j++) {

			temp[i*tcols + j] = 0;
			temp[i*tcols + (tcols - 1 - j)] = 0;
		}
	}

	return temp;
}

Image CV::ApplyFilter(Image & img, Core core, int type) {

	uchar* data = img.GetData();

	double* temp = ApplyFilterRaw(img.GetRowsNumber(), img.GetColsNumber(), data, core, type);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp);

	delete[] temp;
	delete[] data;

	return result;
}
