#pragma once
#include "Image.hpp"

//definition 

struct Dot {

	int x;
	int y;
};

//Core (for filter's functions)
struct Core {
	
	int xk = 0;
	int yk = 0;
	double* data = NULL;
};

//class, that contains CV functions
class ComputerVision {
private:
public:
	static const int kInterpolateZero = 201;			//Interpolate with 0
	static const int kInterpolateBorder = 202;			//Interpolate with border value
	static const int kInterpolateReflection = 203;		//Interpolate with border reflection value
	static const int kInterpolateWrap = 204;			//Interpolate by wrap border

	static const int kPartDerivativeSobel = 301;		//Partial Derivative (Sobel)
	static const int kPartDerivativeXSobel = 302;		//Partial Derivative X (Sobel)
	static const int kPartDerivativeYSobel = 303;		//Partial Derivative Y (Sobel)
	static const int kPartDerivativePrewitt = 304;		//Partial Derivative (Prewitt)
	static const int kPartDerivativeXPrewitt = 305;		//Partial Derivative X (Prewitt)
	static const int kPartDerivativeYPrewitt = 306;		//Partial Derivative Y (Prewitt)
	static const int kPartDerivativeScharr = 307;		//Partial Derivative (Scharr)
	static const int kPartDerivativeXScharr = 308;		//Partial Derivative X (Scharr)
	static const int kPartDerivativeYScharr = 309;		//Partial Derivative Y (Scharr)

	static const int kHarrisResponseDirect = 401;		//response for Harris (Direct by calculate own numbers)
	static const int kHarrisResponseBase = 402;			//response for Harris (Base with k)
	static const int kHarrisResponseForstner = 403;		//response for Harris (Förstner and Gülch)
	
	static const double PI() { return std::atan(1.0) * 4;}

	//apply filter (calculate data*core) to Image
	static Image ApplyFilter(Image& img, Core core, int type);	

	//apply filter (calculate data*core) to Image data
	template <typename T>
	static double* ApplyFilterRaw(int rows, int cols, T* data, Core core, int interpolateType = kInterpolateZero);

	//apply Canny to Image
	static Image Canny(Image& img, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType = kInterpolateZero);
	
	//apply Canny to Image data (or array) 
	template <typename T>
	static double* CannyRaw(int rows, int cols, T* data, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType = kInterpolateZero);

	//create Gauss core
	static Core CreateGaussCore(double sigma, int k);

	//apply Gauss to Image
	static Image Gauss(Image& img, double sigma, int k, int interpolateType = kInterpolateZero);

	//apply Gauss to Image
	static Image GaussDefault(Image& img, double sigma, int interpolateType = kInterpolateZero);

	//apply Gauss to Image data (or array) 
	template <typename T>
	static double* GaussRaw(int rows, int cols, T* data, double sigma, int k, int interpolateType = kInterpolateZero);

	//calculate virtual pixel for interpolation 
	template <typename T>
	static T GetVirtualPixel(int row, int col, int rows, int cols, T* data, int interpolateType = kInterpolateZero);
	
	//apply Harris to Image
	static std::vector<Dot> Harris(Image& img, int wk, int localMinK, double Threshold, int ANMSNeeded = -1, int PartDerivativeType = kPartDerivativeSobel);

	//apply Harris to Image data (or array) 
	template <typename T>
	static std::vector<Dot> HarrisRaw(int rows, int cols, T* data, int wk, int localMinK, double Threshold, int ANMSNeeded = -1, int PartDerivativeType= kPartDerivativeSobel);

	//calculate response for Harris
	static double* HarrisResponse(int rows, int cols, double* A, double* B, double* C, int harrisResponseType = kHarrisResponseBase);

	template<typename srcT, typename dstT>
	static dstT* LinearNormalization(int dataSize, srcT* data, dstT newMin, dstT newMax);	//функция для линейной нормализации
	
	//apply Moravec to Image data (or array) //apply Sobel to Image
	static std::vector<Dot> Moravec(Image& img, int wk, int d, int p, double Threshold);

	//apply Moravec to Image data (or array) 
	template <typename T>
	static std::vector<Dot> MoravecRaw(int rows, int cols, T* data, int wk, int d, int p, double Threshold);

	//calculate Partial Derivative for Image data
	template <typename T>
	static double* PartDerivative(int rows, int cols, T* data, int PartDerivativeType, int interpolateType = kInterpolateZero);

	//apply Sobel to Image
	static Image Sobel(Image& img, int PartDerivativeType, int interpolateType = kInterpolateZero);
	
	//apply Sobel to Image data (or array) 
	template <typename T>
	static double* SobelRaw(int rows, int cols, T* data, int PartDerivativeType, int interpolateType = kInterpolateZero);

	///////////////////////////////////////////////

	//need modify!
	template <typename T>
	static T* ExpandImgByZero(int rows, int cols, T* data, int addX, int addY);
};

///////////////////////////////////
// template function realization //
///////////////////////////////////

template <typename T>
double * ComputerVision::ApplyFilterRaw(int rows, int cols, T* data, Core core, int type) {

	double* result = new double[rows*cols];

	int coreRowSize = 2 * core.xk + 1;
	int center = core.yk * coreRowSize + core.xk;

	int pos;


	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			result[pos] = 0;

			for (int u = -core.yk; u <= core.yk; u++) {
				for (int v = -core.xk; v <= core.xk; v++) {

					if (0 <= (i + u) && (i + u) < rows && 0 <= (j + v) && (j + v) < cols) {
						result[pos] += data[pos + u * cols + v] * core.data[center + u * coreRowSize + v];
					}
					else {
						result[pos] += GetVirtualPixel<T>(i + u, j + v, rows, cols, data, type) * core.data[center + u * coreRowSize + v];
					}
				}
			}
		}
	}

	return result;
}

/////////In work!!!
template<typename T>
double * ComputerVision::CannyRaw(int rows, int cols, T * data, double sigma, int k, double lowThreshhold, double highThreshhold, int interpolateType) {
	
	double* temp;
	double* gaussImage = GaussRaw(rows, cols, data, sigma, k, interpolateType);

	//double* G = SobelRaw(rows, cols, data, kPartDerivativeSobel);
	double* Gx = PartDerivative(rows, cols, gaussImage, kPartDerivativeXSobel);
	double* Gy = PartDerivative(rows, cols, gaussImage, kPartDerivativeYSobel);

	int size = rows * cols;
	double* G = new double[size];
	int* Th = new int[size];

	for (int i = 0; i < size; i++) {

		G[i] = sqrt(Gx[i] * Gx[i] + Gy[i] * Gy[i]);
		Th[i] = 45 * (int)round(4 * atan2(Gy[i], Gx[i]) / PI());
	}

	bool* extremum = new bool[size];
	int pos;
	double neighbourA, neighbourB;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {


			pos = i * cols + j;

			if (i == 0 || i == rows - 1 || j == 0 || j == cols - 1) {
				extremum[pos] = false;
			}
			else {
				pos = i * cols + j;
				neighbourA = 0;
				neighbourB = 0;

				switch (Th[pos]) {

				case 0:
				case 180:
				case -180:

					if (j != 0)
						neighbourA = G[pos - 1];
					if (j != cols - 1)
						neighbourB = G[pos + 1];
					break;

				case 90:
				case -90:

					if (i != 0)
						neighbourA = G[pos - cols];
					if (i != rows - 1)
						neighbourB = G[pos + cols];
					break;

				case 45:
				case -135:

					if (i != 0 && j != cols - 1)
						neighbourA = G[pos - cols + 1];
					if (i != rows - 1 && j != 0)
						neighbourA = G[pos + cols - 1];
					break;

				case -45:
				case 135:

					if (i != 0 && j != 0)
						neighbourA = G[pos - cols - 1];
					if (i != rows - 1 && j != cols - 1)
						neighbourA = G[pos + cols + 1];
					break;

				default:
					break;
				}

				extremum[pos] = (G[pos] > neighbourA && G[pos] > neighbourB);
			}
		}
	}

	temp = LinearNormalization<double, double>(size, G, 0, 100);

	double* result = new double[size];
	for (int i = 0; i < size; i++) {

		if (extremum[i] && temp[i] > highThreshhold)
			result[i] = 1;
		else
			result[i] = 0;
	}

	return result;
}

template<typename T>
double * ComputerVision::GaussRaw(int rows, int cols, T * data, double sigma, int k, int interpolateType) {
	
	double *temp = NULL, *result = NULL;

	Core gaussX;
	gaussX.xk = k;
	gaussX.yk = 0;
	gaussX.data = new double[2 * k + 1];

	Core gaussY;
	gaussY.xk = 0;
	gaussY.yk = k;
	gaussY.data = new double[2 * k + 1];

	double alpha = 1.0 / (sqrt(2 * PI())*sigma);
	//double val;
	
	for (int i = -k; i <= k; i++)
		gaussX.data[k + i] = alpha * std::exp(-(i*i) / (2.0 * sigma * sigma));

	for (int i = -k; i <= k; i++)
		gaussY.data[k + i] = alpha * std::exp(-(i*i) / (2.0 * sigma * sigma));

	temp = ApplyFilterRaw(rows, cols, data, gaussX, interpolateType);
	result = ApplyFilterRaw(rows, cols, temp, gaussY, interpolateType);


	delete[] temp;
	return result;
}

template<typename T>
T ComputerVision::GetVirtualPixel(int row, int col, int rows, int cols, T* data, int interpolateType) {

	if (row >= 0 && row < rows - 1 && col>=0 && col < cols - 1)
		return data[row*cols + col];

	int trow = row;
	int tcol = col;

	switch (interpolateType) {
	case kInterpolateZero:

		return 0;
		break;

	case kInterpolateBorder:

		return data[std::min(rows - 1, std::max(row,0)) *cols + std::min(cols- 1, std::max(col, 0))];
		break;

	case kInterpolateReflection:

		
		if (trow < 0)
			trow = 0 - row;
		else if (trow >= rows)
			trow = rows -1 - trow%(rows -1);

		if (tcol < 0)
			tcol = 0-col;
		else if (tcol >= cols)
			tcol = cols - 1 - tcol%(cols - 1);

		return data[trow*cols + tcol];

	case kInterpolateWrap:

		return data[((rows+row)%rows)*cols + (cols + col) % cols];

	default:
		break;
	}

	return NULL;
}

template<typename T>
std::vector<Dot> ComputerVision::HarrisRaw(int rows, int cols, T * data, int wk, int localMinK, double Threshold, int ANMSNeeded, int PartDerivativeType) {

	int size = rows * cols;
	int pos;

	double *partDerX = NULL, *partDerY = NULL;

	if (PartDerivativeType == kPartDerivativeSobel) {

		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXSobel, kInterpolateBorder);
		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYSobel, kInterpolateBorder);
	}
	else if (PartDerivativeType == kPartDerivativePrewitt) {

		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXPrewitt, kInterpolateBorder);
		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYPrewitt, kInterpolateBorder);
	}
	else if (PartDerivativeType == kPartDerivativeScharr) {

		partDerX = PartDerivative(rows, cols, data, kPartDerivativeXScharr, kInterpolateBorder);
		partDerY = PartDerivative(rows, cols, data, kPartDerivativeYScharr, kInterpolateBorder);
	}

	double* A = new double[size];
	double* B = new double[size];
	double* C = new double[size];

	double Ix, Iy;
	Core gaussCore = CreateGaussCore(wk / 3, wk);
	int gpos;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			A[pos] = 0;
			B[pos] = 0;
			C[pos] = 0;

			for (int u = -wk; u <= wk; u++) {
				for (int v = -wk; v <= wk; v++) {

					gpos = (wk + u)*(2 * wk + 1) + (wk + v);

					Ix = GetVirtualPixel<double>(i + u, j + v, rows, cols, partDerX, kInterpolateZero);
					Iy = GetVirtualPixel<double>(i + u, j + v, rows, cols, partDerY, kInterpolateZero);

					A[pos] += Ix * Ix*gaussCore.data[gpos];
					B[pos] += Ix * Iy*gaussCore.data[gpos];
					C[pos] += Iy * Iy*gaussCore.data[gpos];
				}
			}
		}
	}

	double* lambda = HarrisResponse(rows, cols, A, B, C, ComputerVision::kHarrisResponseForstner);

	Image test(rows, cols, lambda, true);
	test.GetMaxValue();
	cv::imshow("test", test.GetMat());

	std::vector<Dot> result;
	Dot point;
	bool isPoint;
	int r;

	if (ANMSNeeded == -1)
		r = localMinK;
	else
		r = 0;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			isPoint = true;
			pos = i * cols + j;

			if (lambda[pos] > Threshold) {

				for (int u = -r; u <= r && isPoint; u++) {
					for (int v = -r; v <= r && isPoint; v++) {

						if ((r >= sqrt(pow(u, 2) + pow(v, 2))) && !(u == 0 && v == 0) && lambda[pos] <= GetVirtualPixel<T>(i + u, j + v, rows, cols, lambda, kInterpolateZero))
							isPoint = false;
					}
				}
			}
			else {
				isPoint = false;
			}

			if (isPoint) {

				point.x = j;
				point.y = i;

				result.push_back(point);
			}
		}
	}

	/*do {

		result.clear();
		r++;

	} while (result.size() > ANMSNeeded && ANMSNeeded != -1);
	*/

	return result;
}

template<typename srcT, typename dstT>
dstT * ComputerVision::LinearNormalization(int dataSize, srcT * data, dstT newMin, dstT newMax) {

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

template <typename T>
double * ComputerVision::SobelRaw(int rows, int cols, T* data, int PartDerivativeType, int interpolateType) {

	double *temp1 = NULL, *temp2 = NULL, *result = NULL;

	if (PartDerivativeType == kPartDerivativeSobel) {

		temp1 = PartDerivative(rows, cols, data, kPartDerivativeXSobel, interpolateType);
		temp2 = PartDerivative(rows, cols, data, kPartDerivativeYSobel, interpolateType);
	}
	else if (PartDerivativeType == kPartDerivativePrewitt) {

		temp1 = PartDerivative(rows, cols, data, kPartDerivativeXPrewitt, interpolateType);
		temp2 = PartDerivative(rows, cols, data, kPartDerivativeYPrewitt, interpolateType);
	}
	else if (PartDerivativeType == kPartDerivativeScharr) {

		temp1 = PartDerivative(rows, cols, data, kPartDerivativeXScharr, interpolateType);
		temp2 = PartDerivative(rows, cols, data, kPartDerivativeYScharr, interpolateType);
	}

	//Image test1(rows, cols, temp1, true);
	//imshow("Gx image", test1.GetMat());

	//Image test2(rows, cols, temp2, true);
	//imshow("Gy image", test2.GetMat());

	int size = rows * cols;
	result = new double[size];

	for (int i = 0; i < size; i++) {

		result[i] = sqrt(temp1[i] * temp1[i] + temp2[i] * temp2[i]);
	}

	delete[] temp1;
	delete[] temp2;
	return result;
}

template<typename T>
std::vector<Dot> ComputerVision::MoravecRaw(int rows, int cols, T * data, int wk, int d, int p, double Threshold) {
	
	int size = rows*cols;
	int pos;
	
	double** shifts = new double*[8];

	//смещения, начиная с левого и против часовой
	for (int i = 0; i < 8; i++) {

		shifts[i] = new double[size];
	}

	//сначала считаем все смещения для каждого пикселя (их всего 8)
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			shifts[0][pos] = pow(data[pos] - GetVirtualPixel<T>(i, j + d, rows, cols, data, kInterpolateReflection), 2);
			shifts[1][pos] = pow(data[pos] - GetVirtualPixel<T>(i - d, j + d, rows, cols, data, kInterpolateReflection), 2);
			shifts[2][pos] = pow(data[pos] - GetVirtualPixel<T>(i - d, j, rows, cols, data, kInterpolateReflection), 2);
			shifts[3][pos] = pow(data[pos] - GetVirtualPixel<T>(i - d, j - d, rows, cols, data, kInterpolateReflection), 2);
			shifts[4][pos] = pow(data[pos] - GetVirtualPixel<T>(i, j - d, rows, cols, data, kInterpolateReflection), 2);
			shifts[5][pos] = pow(data[pos] - GetVirtualPixel<T>(i + d, j - d, rows, cols, data, kInterpolateReflection), 2);
			shifts[6][pos] = pow(data[pos] - GetVirtualPixel<T>(i + d, j, rows, cols, data, kInterpolateReflection), 2);
			shifts[7][pos] = pow(data[pos] - GetVirtualPixel<T>(i + d, j + d, rows, cols, data, kInterpolateReflection), 2);
		}
	}

	double temp;
	double* s = new double[size];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			pos = i * cols + j;

			for (int shift = 0; shift < 8; shift++) {

				temp = 0;

				for (int u = -wk; u <= wk; u++) {
					for (int v = -wk; v <= wk; v++) {

						if (u != 0 && v != 0)
							temp += GetVirtualPixel<T>(i + u, j + v, rows, cols, shifts[shift], kInterpolateBorder);
					}
				}

				if (shift == 0)
					s[pos] = temp;
				else
					s[pos] = std::min(s[pos], temp);
			}
		}
	}

	//Image test(rows, cols, s, true);
	//cv::imshow("test", test.GetMat());

	std::vector<Dot> result;
	Dot point;
	bool isPoint;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			isPoint = true;
			pos = i * cols + j;

			if (s[pos] > Threshold) {

				for (int u = -p; u <= p; u++) {
					for (int v = -p; v <= p; v++) {

						if (!(u == 0 && v == 0) && s[pos] <= GetVirtualPixel<T>(i + u, j + v, rows, cols, s, kInterpolateZero))
							isPoint = false;
					}
				}
			}
			else {
				isPoint = false;
			}

			if (isPoint) {

				point.x = j;
				point.y = i;

				result.push_back(point);
			}
		}
	}

	for (int shift = 0; shift < 8; shift++)
		delete[] shifts[shift];

	delete[] shifts;
	delete[] s;

	return result;
}

template<typename T>
double * ComputerVision::PartDerivative(int rows, int cols, T * data, int PartDerivativeType, int interpolateType) {

	Core Gp1, Gp2;// p(x or y) gradient's filters (separable)

	Gp1.xk = 1;
	Gp1.yk = 0;
	Gp1.data = new double[3];

	Gp2.xk = 0;
	Gp2.yk = 1;
	Gp2.data = new double[3];


	switch (PartDerivativeType) {
	case kPartDerivativeXSobel:

		Gp1.data[0] = 1; Gp1.data[1] = 0; Gp1.data[2] = -1;
		Gp2.data[0] = 1; Gp2.data[1] = 2; Gp2.data[2] = 1;
		break;

	case kPartDerivativeYSobel:

		Gp1.data[0] = 1; Gp1.data[1] = 2; Gp1.data[2] = 1;
		Gp2.data[0] = 1; Gp2.data[1] = 0; Gp2.data[2] = -1;
		break;

	case kPartDerivativeXPrewitt:

		Gp1.data[0] = 1; Gp1.data[1] = 0; Gp1.data[2] = -1;
		Gp2.data[0] = 1; Gp2.data[1] = 1; Gp2.data[2] = 1;
		break;

	case kPartDerivativeYPrewitt:

		Gp1.data[0] = 1; Gp1.data[1] = 1; Gp1.data[2] = 1;
		Gp2.data[0] = 1; Gp2.data[1] = 0; Gp2.data[2] = -1;
		break;

	case kPartDerivativeXScharr:

		Gp1.data[0] = sqrt(3); Gp1.data[1] = 0; Gp1.data[2] = -sqrt(3);
		Gp2.data[0] = sqrt(3); Gp2.data[1] = 10.0 / (3 * sqrt(3)); Gp2.data[2] = sqrt(3);
		break;

	case kPartDerivativeYScharr:

		Gp1.data[0] = sqrt(3); Gp1.data[1] = 10.0 / (3 * sqrt(3)); Gp1.data[2] = sqrt(3);
		Gp2.data[0] = sqrt(3); Gp2.data[1] = 0; Gp2.data[2] = -sqrt(3);
		break;

	default:
		break;
	}
	double *temp, *result;

	temp = ApplyFilterRaw(rows, cols, data, Gp1, interpolateType);
	result = ApplyFilterRaw(rows, cols, temp, Gp2, interpolateType);

	delete[] temp;

	return result;
}

template<typename T>
T* ComputerVision::ExpandImgByZero(int rows, int cols, T* data, int addX, int addY) {

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
