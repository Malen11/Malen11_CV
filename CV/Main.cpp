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
int main() {

	Mat mat = LoadImage();

	Image test(mat);
	imshow("Base image", test.GetMat());
	
	double* ddata = test.GetNormalizeDataF();
	Image test1(test.GetRowsNumber(), test.GetColsNumber(), ddata,true);
	imshow("Normalize image", test1.GetMat());

	Image test2 = CV::Sobel(test1);
	imshow("Sobel", test2.GetMat());

	//test2 = CV::Sobel(test);
	//imshow("Sobel", test2.GetMat());
	waitKey();
/*	VideoCapture cap = StartCapture();

	while (true) {

		cap >> frame;
		test = Image(frame);

		imshow("test1", test.GetMat());
		test2 = CV::Sobel(test);
		imshow("test2", test2.GetMat());
		waitKey(30);
	}
	*/
	return 0;
}
