#pragma once
#include "Image.hpp"

struct Core {
	
	int k;
	double* data;
};

const class CV {
private:
public:
	uchar* ApplyFilterRaw(int rows, int cols, uchar* data, Core core);
	Image ApplyFilter(Image& img, Core core);
};