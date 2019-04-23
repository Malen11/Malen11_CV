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

int main() {

	Mat frame;

	VideoCapture cap = StartCapture();

	//namedWindow("test", CV_WINDOW_NORMAL);

	Core myCore;
	myCore.k = 1;
	myCore.data = new double[(2 * myCore.k + 1)*(2 * myCore.k + 1)];

	int center = myCore.k * (2 * myCore.k + 1) + myCore.k;
	double q= 1.0 / ((2 * myCore.k + 1)*(2 * myCore.k + 1));

	for (int i = -myCore.k; i <= myCore.k; i++) {
		for (int j = -myCore.k; j <= myCore.k; j++) {

			myCore.data[center + i*(2*myCore.k+1) + j] = q;
		}
	}

	Image test;
	Image test2;

	while (true) {

		cap >> frame;
		test = Image(frame);

		imshow("test1", test.GetMat());
		test2 = CV::ApplyFilter(test, myCore, CV::kZeroBorder);
		imshow("test2", test2.GetMat());
		waitKey(30);
	}

	return 0;
}
