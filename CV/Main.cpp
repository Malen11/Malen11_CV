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
	
	//double* ddata = test.GetNormalizeDataF();
	//Image test1(test.GetRowsNumber(), test.GetColsNumber(), ddata,true);
	//imshow("Normalize image", test1.GetMat());

	Image test2, test3, test4, test5;

	//lab 1
	/*
	test2 = ComputerVision::Sobel(test, ComputerVision::kPartDerivativeSobel,ComputerVision::kInterpolateZero);
	imshow("Sobel Before Gauss (Sobel)", test2.GetMat());

	test4 = ComputerVision::Sobel(test, ComputerVision::kPartDerivativeScharr, ComputerVision::kInterpolateZero);
	imshow("Sobel Before Gauss (Scharr)", test4.GetMat());

	test5 = test5.AbsoluteDiff(test2, test4);
	test5.NormalizeImage();
	imshow("Sobel Before Gauss (AbsDiff Sobel, Scharr)", test5.GetMat());

	//test3 = ComputerVision::GaussDefault(test, 2, ComputerVision::kInterpolateZero);
	//imshow("Gauss", test3.GetMat());

	//test2 = ComputerVision::Sobel(test3, ComputerVision::kPartDerivativeScharr, ComputerVision::kInterpolateZero);
	//imshow("Sobel After Gauss", test2.GetMat());
	*/
	
	//lab 2
	
	ImagePyramid imgPyr1(test, 4, 0.5, 1.5);

	ImagePyramid imgPyr2(test, 2, 0.5, 1.5);

	test2 = imgPyr1.GetImage(0, 3);

	test3 = imgPyr2.GetImage(0, 1);

	test4 = test2.AbsoluteDiff(test2, test3);
	test4.NormalizeImage();

	imshow("2-4 layers", test4.AbsoluteDiff(test2, test3).GetMat());

	/*for (int i = 0; i < imgPyr.GetOctavesNum(); i++) {
		for (int j = 0; j < imgPyr.GetLayersNum(); j++) {
			cv::imwrite("pyramid/Octave " + std::to_string(i) + " Layer " + std::to_string(j)+".png", imgPyr.GetImage(i, j).GetMat());
			//imshow("Octave: "+std::to_string(i)+" Layer: "+ std::to_string(j), imgPyr.GetImage(i, j).GetMat());
		}
	}
	*/

	//imgPyr.L()

	//lab 3
	//test2 = ComputerVision::Canny(test, 2, 6, 10, 10);
	//imshow("Canny", test2.GetMat());
	//video 
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
