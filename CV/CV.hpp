#pragma once
#include "Image.hpp"

struct Core {
	
	int k = 0;
	double* data = NULL;
};

/*struct SeparableCore {

	int k = 0;
	double* dataX = NULL;
	double* dataY = NULL;
};*/

class CV {
private:
public:
	static const int kZeroBorder = 201;			//заполнить края 0
	static const int kExpandBorder = 202;		//использовать крайний пиксель
	static const int kReflectionBorder = 203;	//отразить границу
	static const int kWrapBorder = 204;			//заворачивать границу

	template <typename T>
	static double* ApplyFilterRaw(int rows, int cols, T* data, Core core, int type);
	template <typename T>
	static T* ExpandImgByZero(int rows, int cols, T* data, int k);
	static Image ApplyFilter(Image& img, Core core, int type);
};