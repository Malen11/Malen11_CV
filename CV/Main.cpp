// Main.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

#include "Image.hpp"
#include "ImagePyramid.hpp"
#include "ComputerVision.hpp"
#include "ImageFilters.hpp"

using namespace std;
using namespace cv;
using namespace CV_labs;

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
//cv::Mat MatPlotPoints(Image& img, vector<Dot> points, int color = 0) {
//
//	cv::Mat temp = img.GetMat();
//	cv::Mat colored;
//	cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);
//	
//	for (vector<Dot>::iterator it = points.begin(); it != points.end(); it++) {
//
//		circle(colored, Point((*it).x, (*it).y), 1, CV_RGB(255, 0, 0), 3);
//	}
//
//	return colored;
//}

//cv::Mat MatPlotLines(Image& img, vector<PairDot> lines, int color = 0) {
//
//	cv::Mat temp = img.GetMat();
//	cv::Mat colored;
//	cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);
//	RNG rng(12345);
//	
//	for (vector<PairDot>::iterator it = lines.begin(); it != lines.end(); it++) {
//		Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
//		line(colored, Point((*it).point0.x, (*it).point0.y), Point((*it).point1.x, (*it).point1.y), color, 2);
//	}
//
//	return colored;
//}

void lab1() {

	Mat mat = LoadImage();
	Image test(mat);
	test.NormalizeImage();
	imshow("Base image", test.GetMat());

	Image gaussSepF = ImageFilters::ApplyFilter(test, ImageFilters::GenerateGaussSeparableCore(2.6), ImageFilters::kInterpolateReflection);
	Image gaussF = ImageFilters::ApplyFilter(test, ImageFilters::GenerateGaussCore(2.6), ImageFilters::kInterpolateZero);

	Image sobelX = ImageFilters::ApplyFilter(test, ImageFilters::GenerateSobelCore(ImageFilters::kPartDerivativeDirectionX));
	Image sobelY = ImageFilters::ApplyFilter(test, ImageFilters::GenerateSobelCore(ImageFilters::kPartDerivativeDirectionY));

	Image scharrSepX = ImageFilters::ApplyFilter(test, ImageFilters::GenerateScharrSeparableCore(ImageFilters::kPartDerivativeDirectionX));
	Image prewittSepY = ImageFilters::ApplyFilter(test, ImageFilters::GeneratePrewittSeparableCore(ImageFilters::kPartDerivativeDirectionY));

	Image gradientValSobel = ImageFilters::CalculateGradientValue(test, ImageFilters::kPartDerivativeTypeSobelCore, ImageFilters::kInterpolateZero);
	Image gradientValScharr = ImageFilters::CalculateGradientValue(test, ImageFilters::kPartDerivativeTypeScharrCore, ImageFilters::kInterpolateZero);

	Image absDif = gradientValSobel.CalcAbsoluteDiff(gradientValScharr);
	absDif.NormalizeImage();
	imshow("Gauss separable filter", gaussSepF.GetMat());
	imshow("Gauss filter", gaussF.GetMat());
	imshow("Sobel by X", sobelX.GetMat());
	imshow("Sobel by Y", sobelY.GetMat());
	imshow("Scharr by X", scharrSepX.GetMat());
	imshow("Prewitt by Y", prewittSepY.GetMat());
	imshow("Gradient Values Sobel", gradientValSobel.GetMat());
	imshow("Gradient Values Scharr", gradientValScharr.GetMat());
	imshow("Gradient Values Sobel abs dif Scharr", absDif.GetMat());
}

void lab2() {

	Mat mat = LoadImage();
	Image test(mat);
	test.NormalizeImage();
	imshow("Base image", test.GetMat());

	ImagePyramid imagePyramid(test, 0, 1.5, 3);

	//for(int i = 0; i < imagePyramid.GetOctavesNum(); i++) {
	//	for (int j = 0; j < imagePyramid.GetLayersNum(); j++) {

	//		string name = "octave " + std::to_string(i) + " layer " + std::to_string(j) + " sigma " + std::to_string(imagePyramid.GetImageSigma(i, j));
	//		//imshow(name, imagePyramid.GetImage(i, j).GetMat());
	//		cv::imwrite("pyramid/" + name + ".png", imagePyramid.GetImage(i, j).GetMat());
	//	}
	//}

	double sigma = 1.5 * 8;
	Image test2 = ImageFilters::ApplyFilter(test, ImageFilters::GenerateGaussSeparableCore(sigma), ImageFilters::kInterpolateReflection);

	int cols = test2.GetColsNumber();
	int rows = test2.GetRowsNumber();
	uchar* imgFromPyr = new uchar[cols *rows];


	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {

			imgFromPyr[i * cols + j] = imagePyramid.L(i, j, sigma);
		}
	}
	Image imgFromPyrRe(rows, cols, imgFromPyr);
	imshow("smothed image", test2.GetMat());
	imshow("rebuild image", imgFromPyrRe.GetMat());
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

	lab2();

	//Mat mat1 = LoadImage();
	//Image test1(mat1);
	//imshow("Base image 1", test1.GetMat());

	//Image emptyImg(test.GetRowsNumber() + test1.GetRowsNumber(), test.GetColsNumber() + test1.GetColsNumber(), 0);
	//
	////double* ddata = test.GetNormalizeDataF();
	////Image test1(test.GetRowsNumber(), test.GetColsNumber(), ddata,true);
	////imshow("Normalize image", test1.GetMat());

	//Image test2, test3, test4, test5;

	

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

	//lab 4-5

	
	//test.RotateImage(-ComputerVision::PI()/4);
	//imshow("rotate!", test.GetMat());


	//test3 =  test;//ComputerVision::GaussDefault(test, 1.2, ComputerVision::kInterpolateReflection);//
	//vector<Dot> pointsANMS = ComputerVision::Harris(test3, 4, 12, 0.03, 50);
	//imshow("Harris Points (ANMS)", MatPlotPoints(test3, pointsANMS));

	//vector<Descriptor> descriptors = ComputerVision::CreateDescriptors(test3, pointsANMS, 16, 16, 2, 2, 8,ComputerVision::kDescriptorSimple, ComputerVision::kDescriptorNormalization2Times);

	//test4 =  test1;//ComputerVision::GaussDefault(test1, 1.2, ComputerVision::kInterpolateReflection);//
	//vector<Dot> pointsANMS1 = ComputerVision::Harris(test4, 4, 4, 0.03, 50);
	//imshow("Harris Points 1 (ANMS)", MatPlotPoints(test4, pointsANMS1));

	//vector<Descriptor> descriptors1 = ComputerVision::CreateDescriptors(test4, pointsANMS1, 16, 16, 2, 2, 8, ComputerVision::kDescriptorSimple, ComputerVision::kDescriptorNormalization2Times);
	//emptyImg=emptyImg.InsertImage(test, 0, 0).InsertImage(test1, test.GetColsNumber(), test.GetRowsNumber());
	//
	////vector<PairDot> matchedDesc = ComputerVision::DescriptorsMatching(descriptors, descriptors1);
	////vector<PairDot> matchedDesc = ComputerVision::DescriptorsMatching(descriptors, descriptors1, ComputerVision::kDescriptorsComparisonEuclid, ComputerVision::kDescriptorsMatchingNNDR,0.8);
	//vector<PairDot> matchedDesc = ComputerVision::DescriptorsMatching(descriptors, descriptors1, ComputerVision::kDescriptorsComparisonEuclid, ComputerVision::kDescriptorsMatchingMutal); 
	//vector<PairDot> lines(matchedDesc.size());

	//PairDot line;
	//for (int i = 0; i < matchedDesc.size(); i++) {
	//	line.point0 = matchedDesc[i].point0;

	//	line.point1.x = matchedDesc[i].point1.x + test.GetColsNumber();
	//	line.point1.y = matchedDesc[i].point1.y + test.GetRowsNumber();
	//	lines.push_back(line);
	//}

	//Mat colored = MatPlotLines(emptyImg, lines);
	//
	//double k = 1.2;
	//resize(colored, colored, Size(k*test.GetColsNumber(), k*test.GetRowsNumber()));
	//imshow("inserted image", colored);
	

	//lab6 


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