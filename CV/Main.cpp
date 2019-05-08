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
	//int selected;
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

/** function MatPlotPoints */
cv::Mat MatPlotPoints(Image& img, vector<Dot> points, int color = 0) {

	cv::Mat temp = img.GetMat();
	cv::Mat colored;
	cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);
	
	for (vector<Dot>::iterator it = points.begin(); it != points.end(); it++) {

		circle(colored, Point((*it).x, (*it).y), 1, CV_RGB(255, 0, 0), 3);
	}

	return colored;
}

cv::Mat MatPlotLines(Image& img, vector<PairDot> lines, int color = 0) {

	cv::Mat temp = img.GetMat();
	cv::Mat colored;
	cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);
	RNG rng(12345);
	
	for (vector<PairDot>::iterator it = lines.begin(); it != lines.end(); it++) {
		Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
		line(colored, Point((*it).point0.x, (*it).point0.y), Point((*it).point1.x, (*it).point1.y), color, 2);
	}

	return colored;
}

//C:\Users\alist\Desktop\cv\Bikesgray.jpg
//C:\Users\alist\Desktop\cv\obj.jpg
//C:\Users\alist\Desktop\cv\Lena2.jpg
//C:\Users\alist\Desktop\cv\Lena2lheight.jpg
//C:\Users\alist\Desktop\cv\Lena2contrast.jpg
//C:\Users\alist\Desktop\cv\Lena2light.jpg
//C:\Users\alist\Desktop\cv\LenaRotate.jpg
//C:\Users\alist\Desktop\cv\LenaRotate2.jpg
//C:\Users\alist\Desktop\cv\Valve_original_(1).PNG
int main() {

	Mat mat = LoadImage();
	Image test(mat);
	imshow("Base image", test.GetMat());

	Mat mat1 = LoadImage();
	Image test1(mat1);
	imshow("Base image 1", test1.GetMat());

	Image emptyImg(test.GetRowsNumber() + test1.GetRowsNumber(), test.GetColsNumber() + test1.GetColsNumber(), 0);
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
/*
	ImagePyramid imgPyr1(test, 4, 0.5, 1.5);

	ImagePyramid imgPyr2(test, 2, 0.5, 1.5);

	test2 = imgPyr1.GetImage(0, 3);

	test3 = imgPyr2.GetImage(0, 1);

	test4 = test2.AbsoluteDiff(test2, test3);
	test4.NormalizeImage();

	imshow("2-4 layers", test4.AbsoluteDiff(test2, test3).GetMat());

	for (int i = 0; i < imgPyr.GetOctavesNum(); i++) {
		for (int j = 0; j < imgPyr.GetLayersNum(); j++) {
			cv::imwrite("pyramid/Octave " + std::to_string(i) + " Layer " + std::to_string(j)+".png", imgPyr.GetImage(i, j).GetMat());
			//imshow("Octave: "+std::to_string(i)+" Layer: "+ std::to_string(j), imgPyr.GetImage(i, j).GetMat());
		}
	}
	*/

	//lab 3
	/*
	//test2 = ComputerVision::Canny(test, 2, 6, 40, 30);
	//imshow("Canny", test2.GetMat());

	test2 = ComputerVision::GaussDefault(test, 1);
	vector<Dot> pointsThresh = ComputerVision::Harris(test,3,3,0.03);
	imshow("Harris Points (Thresh)", MatPlotPoints(test, pointsThresh));

	vector<Dot> pointsANMS = ComputerVision::Harris(test, 3, 3, 0.03, 40);
	imshow("Harris Points (ANMS)", MatPlotPoints(test, pointsANMS));


	test3 = ComputerVision::GaussDefault(test1, 1);
	pointsThresh = ComputerVision::Harris(test1, 3, 3, 0.03);
	imshow("Harris Points 1 (Thresh)", MatPlotPoints(test1, pointsThresh));

	pointsANMS = ComputerVision::Harris(test1, 3, 3, 0.03, 40);
	imshow("Harris Points 1 (ANMS)", MatPlotPoints(test1, pointsANMS));
	
	//test.RotateImage(-0.2);
	//imshow("rotate!", test.GetMat());
	
	//vector<Dot> points = ComputerVision::Moravec(test, 3, 1,3,0.02);
	//imshow("Moravec Points", MatPlotPoints(test, points));
	*/

	//lab4
	test3 = test;// ComputerVision::GaussDefault(test, 1.2, ComputerVision::kInterpolateReflection);
	vector<Dot> pointsANMS = ComputerVision::Harris(test3, 3, 3, 0.03, 50);
	imshow("Harris Points (ANMS)", MatPlotPoints(test3, pointsANMS));

	vector<Descriptor> descriptors = ComputerVision::CreateDescriptors(test3, pointsANMS, 16, 16, 2, 2, 8,ComputerVision::kDescriptorSimple, ComputerVision::kDescriptorNormalization2Times);

	test4 = test1;// ComputerVision::GaussDefault(test1, 1.2, ComputerVision::kInterpolateReflection);
	vector<Dot> pointsANMS1 = ComputerVision::Harris(test4, 3, 3, 0.03, 50);
	imshow("Harris Points 1 (ANMS)", MatPlotPoints(test4, pointsANMS1));

	vector<Descriptor> descriptors1 = ComputerVision::CreateDescriptors(test4, pointsANMS1, 16, 16, 2, 2, 8, ComputerVision::kDescriptorSimple, ComputerVision::kDescriptorNormalization2Times);
	emptyImg=emptyImg.InsertImage(test, 0, 0).InsertImage(test1, test.GetColsNumber(), test.GetRowsNumber());
	
	//vector<PairDot> matchedDesc = ComputerVision::DescriptorsMatching(descriptors, descriptors1);
	//vector<PairDot> matchedDesc = ComputerVision::DescriptorsMatching(descriptors, descriptors1, ComputerVision::kDescriptorsComparisonEuclid, ComputerVision::kDescriptorsMatchingNNDR,0.8);
	vector<PairDot> matchedDesc = ComputerVision::DescriptorsMatching(descriptors, descriptors1, ComputerVision::kDescriptorsComparisonEuclid, ComputerVision::kDescriptorsMatchingMutal); 
	vector<PairDot> lines(matchedDesc.size());

	PairDot line;
	for (int i = 0; i < matchedDesc.size(); i++) {
		line.point0 = matchedDesc[i].point0;

		line.point1.x = matchedDesc[i].point1.x + test.GetColsNumber();
		line.point1.y = matchedDesc[i].point1.y + test.GetRowsNumber();
		lines.push_back(line);
	}

	Mat colored = MatPlotLines(emptyImg, lines);
	/*
	int* matchedPoints = new int[descriptors.size()];
	double bestMatching, matching;

	for (int i = 0; i < descriptors.size(); i++) {

		matchedPoints[i]= 0;
		bestMatching = ComputerVision::DescriptorsDifference(descriptors[i], descriptors1[0]);

		for (int j = 1; j < descriptors1.size(); j++) {

			matching = ComputerVision::DescriptorsDifference(descriptors[i], descriptors1[j]);
			if (matching < bestMatching) {
				bestMatching = matching;
				matchedPoints[i] = j;
			}
		}
	}
	
	Mat temp = emptyImg.GetMat();
	cv::Mat colored;
	cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);

	for (int i = 0; i < descriptors.size(); i++) {

		line(colored, Point(descriptors[i].point.x, descriptors[i].point.y), Point(test.GetColsNumber()+descriptors1[matchedPoints[i]].point.x, test.GetRowsNumber() + descriptors1[matchedPoints[i]].point.y), CV_RGB(255, 0, 0),2);
	}*/
	
	double k = 1.2;
	resize(colored, colored, Size(k*test.GetColsNumber(), k*test.GetRowsNumber()));
	imshow("inserted image", colored);

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
