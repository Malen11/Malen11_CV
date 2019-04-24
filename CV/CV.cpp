#include "stdafx.h"
#include "CV.hpp"

using namespace std;

Image CV::Sobel(Image & img) {

	double* data = img.GetNormalizeDataF();
	double* temp = SobelRaw(img.GetRowsNumber(), img.GetColsNumber(), data);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp, true);

	delete[] temp;
	delete[] data;

	return result;
}

Image CV::ApplyFilter(Image & img, Core core, int type) {

	uchar* data = img.GetData();

	double* temp = ApplyFilterRaw(img.GetRowsNumber(), img.GetColsNumber(), data, core, type);

	Image result(img.GetRowsNumber(), img.GetColsNumber(), temp);

	delete[] temp;
	delete[] data;

	return result;
}
