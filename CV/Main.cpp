// Main.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

using namespace std;
using namespace cv;

/** function StartCapture */
VideoCapture StartCapture() {

	string input;
	int selected;
	VideoCapture cap;

	while (true) {

		cout << "Use camera(1) or open existing file(2)? Input: ";
		cin >> input;

		selected = std::atoi(input.c_str());

		if (selected == 1) {

			cap.open(0); // open the default camera

			if (cap.isOpened()) {

				cout << "Success!" << endl;
				break;
			}
			else {

				cout << "Can't open the camera!" << endl;
				continue;
			}
		}
		else if (selected == 2) {

			while (true) {

				cout << "Enter path to file: ";
				cin >> input;
				cap.open(input);

				if (cap.isOpened()) {

					cout << "Success!" << endl;
					break;
				}
				else {

					cout << "Wrong file name!" << endl;
					continue;
				}
			}

			break;
		}
		else {
			cout << "Wrong input!" << endl;
		}
	}

	return cap;
}

/** function LoadImage */
Mat LoadImage() {

	string input;
	int selected;
	cv::Mat img;
	while (true) {

		cout << "Enter path to file: ";
		cin >> input;
		img = cv::imread(input);

		if (img.data != NULL) {

			cout << "Success!" << endl;
			break;
		}
		else {

			cout << "Wrong file name!" << endl;
			continue;
		}
	}

	return img;
}

//C:\Users\alist\Desktop\Bikesgray.jpg
//C:\Users\alist\Desktop\Valve_original_(1).PNG
int main() {

	Mat mat = LoadImage();

	Image test(mat);
	imshow("Base image", test.GetMat());
	
	double* ddata = test.GetNormalizeDataF();
	Image test1(test.GetRowsNumber(), test.GetColsNumber(), ddata,true);
	imshow("Normalize image", test1.GetMat());

	Image test2, test3;

	test2 = ComputerVision::Sobel(test1, ComputerVision::kPartDerivativeSobel,ComputerVision::kInterpolateZero);
	imshow("Sobel Before Gauss", test2.GetMat());

	test3 = ComputerVision::GaussDefault(test, 2, ComputerVision::kInterpolateZero);
	imshow("Gauss", test3.GetMat());

	test2 = ComputerVision::Sobel(test3, ComputerVision::kPartDerivativeScharr, ComputerVision::kInterpolateZero);
	imshow("Sobel After Gauss", test2.GetMat());

	//test2 = ComputerVision::Sobel(test1, ComputerVision::kPartDerivativePrewitt, ComputerVision::kInterpolateReflection);
	//imshow("Prewitt", test2.GetMat());

	//test2 = ComputerVision::Sobel(test1, ComputerVision::kPartDerivativeScharr, ComputerVision::kInterpolateWrap);
	//imshow("Scharr", test2.GetMat());


	/*int scale = 1; 
	int delta = 0;
	int ddepth = CV_16S;
	
	Mat src_gray;
	/// Convert it to gray
	cvtColor(mat, src_gray, CV_BGR2GRAY);
	imshow("Gray (OpenCV)", src_gray);

	/// Generate grad_x and grad_y
	Mat grad_x, grad_y;
	Mat abs_grad_x, abs_grad_y;

	/// Gradient X
	//Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
	Sobel(src_gray, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT);
	convertScaleAbs(grad_x, abs_grad_x);

	/// Gradient Y
	//Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
	Sobel(src_gray, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT);
	convertScaleAbs(grad_y, abs_grad_y);

	/// Total Gradient (approximate)
	addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, mat);
	imshow("Sobel (OpenCV)", mat);
	//test2 = ComputerVision::Sobel(test);
	//imshow("Sobel", test2.GetMat());*/

/*	
VideoCapture cap = StartCapture();

	while (true) {

		cap >> frame;
		test = Image(frame);

		imshow("test1", test.GetMat());
		test2 = ComputerVision::Sobel(test);
		imshow("test2", test2.GetMat());
		waitKey(30);
	}
	*/

	waitKey();

	return 0;
}
