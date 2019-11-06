// Main.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

#include "Image.hpp"
#include "ImagePyramid.hpp"
#include "ComputerVision.hpp"
#include "ImageFilters.hpp"
#include "ImageDetectors.hpp"
#include "ImageDescriptors.hpp"

using namespace std;
using namespace cv;
using namespace CV_labs;

double testAlpha = 0;

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
cv::Mat MatPlotPoints(Image& img, vector<CV_labs::Point> points, int color = 0) {

	cv::Mat temp = img.GetMat();
	cv::Mat colored;
	cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);
	
	for (vector<CV_labs::Point>::iterator it = points.begin(); it != points.end(); it++) {

		circle(colored, cv::Point((*it).x, (*it).y), 1, CV_RGB(255, 0, 0), 3);
	}

	return colored;
}

cv::Mat MatPlotLines(Image& img, vector<CV_labs::Line> lines, int color = 0) {

	cv::Mat temp = img.GetMat();
	cv::Mat colored;
	cv::cvtColor(temp, colored, cv::COLOR_GRAY2BGR);
	RNG rng(12345);
	
	for (vector<CV_labs::Line>::iterator it = lines.begin(); it != lines.end(); it++) {
		Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
		line(colored, cv::Point((*it).pointA.x, (*it).pointA.y), cv::Point((*it).pointB.x, (*it).pointB.y), color, 2);
	}

	return colored;
}

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

	ImagePyramid imagePyramid(test, 0, 1.5, 5);

	for(int i = 0; i < imagePyramid.GetOctavesNum(); i++) {
		for (int j = 0; j < imagePyramid.GetLayersNum(); j++) {

			string name = "octave " + std::to_string(i) + " layer " + std::to_string(j) + " sigma " + std::to_string(imagePyramid.GetImageSigma(i, j));
			//imshow(name, imagePyramid.GetImage(i, j).GetMat());
			cv::imwrite("pyramid/" + name + ".png", imagePyramid.GetImage(i, j).GetMat());
		}
	}

	double sigma = 1.5 * 4 * 1.4;
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

void lab3() {

	Mat mat = LoadImage();
	Image base(mat);
	base.NormalizeImage();
	imshow("Base image", base.GetMat());

	Image test1 = ImageFilters::ApplyFilter(base, ImageFilters::GenerateGaussSeparableCore(1.5), ImageFilters::kInterpolateReflection);
	imshow("Base image gauss", test1.GetMat());
	Image test2 = ImageFilters::ApplyFilter(base.GetRotatedImage(-0.2), ImageFilters::GenerateGaussSeparableCore(1.5), ImageFilters::kInterpolateReflection);
	imshow("Rotate gauss!", test2.GetMat());

	/*test2 = ComputerVision::Canny(test, 2, 6, 40, 30);
	imshow("Canny", test2.GetMat());*/

	//vector<CV_labs::Point> points = ImageDetectors::Moravec(base, 3, 1, 3, 0.02);
	//imshow("Moravec Points base", MatPlotPoints(base, points));

	//vector<CV_labs::Point> points1 = ImageDetectors::Moravec(test1, 3, 1, 3, 0.02);
	//imshow("Moravec Points gauss", MatPlotPoints(test1, points1));

	//vector<CV_labs::Point> points2 = ImageDetectors::Moravec(test2, 3, 1, 3, 0.02);
	//imshow("Moravec Points gauss rotated", MatPlotPoints(test2, points2));

	/*double* data = test1.GetNormalizedImageDataD();
	double* partDerX = ImageFilters::ApplyFilterRaw(test1.GetRowsNumber(), test1.GetColsNumber(), data, ImageFilters::GenerateSobelSeparableCore(ImageFilters::kPartDerivativeDirectionX), ImageFilters::kInterpolateBorder);
	partDerX = ImageFilters::NormalizeData(test1.GetRowsNumber() * test1.GetColsNumber(), partDerX);
	Image partDerXImg(test1.GetRowsNumber(), test1.GetColsNumber(), partDerX);
	imshow("partDerXImg", partDerXImg.GetMat());*/

	vector<CV_labs::Point> pointsThresh1 = ImageDetectors::Harris(test1, 5, 5, 0.03);
	imshow("Harris Points (Thresh)", MatPlotPoints(test1, pointsThresh1));

	vector<CV_labs::Point> pointsANMS1 = ImageDetectors::Harris(test1, 5, 5, 0.03, 40);
	imshow("Harris Points (ANMS)", MatPlotPoints(test1, pointsANMS1));


	vector<CV_labs::Point> pointsThresh2 = ImageDetectors::Harris(test2, 5, 5, 0.03);
	imshow("Harris Points rotate (Thresh)", MatPlotPoints(test2, pointsThresh2));

	vector<CV_labs::Point> pointsANMS2 = ImageDetectors::Harris(test2, 5, 5, 0.03, 40);
	imshow("Harris Points rotate (ANMS)", MatPlotPoints(test2, pointsANMS2));
}

void lab4() {

	//lab 4

	Mat mat1 = LoadImage();
	Image base1(mat1);
	base1.NormalizeImage();
	imshow("Base image 1", base1.GetMat());

	Image test1 = ImageFilters::ApplyFilter(base1, ImageFilters::GenerateGaussSeparableCore(1.5), ImageFilters::kInterpolateReflection);
	imshow("Base image 1 gauss", test1.GetMat());

	Mat mat2 = LoadImage();
	Image base2(mat2);
	base2.NormalizeImage();
	imshow("Base image 2", base2.GetMat());

	Image test2 = ImageFilters::ApplyFilter(base2, ImageFilters::GenerateGaussSeparableCore(1.5), ImageFilters::kInterpolateReflection);
	imshow("Base image 2 gauss", test2.GetMat());

	vector<CV_labs::Point> pointsANMS1 = ImageDetectors::Harris(test1, 5, 5, 0.03, 40);
	imshow("Harris Points image 1 (ANMS)", MatPlotPoints(test1, pointsANMS1));

	vector<CV_labs::Point> pointsANMS2 = ImageDetectors::Harris(test2, 5, 5, 0.03, 40);
	imshow("Harris Points image 2 (ANMS)", MatPlotPoints(test2, pointsANMS2));

	vector<CV_labs::Descriptor> descriptors1 = ImageDescriptorMethods::CreateDescriptors(test1, pointsANMS1, 16, 4, 16, ImageDescriptorMethods::kDescriptorTypeSquare, ImageDescriptorMethods::kDescriptorNormalization2Times);
	vector<CV_labs::Descriptor> descriptors2 = ImageDescriptorMethods::CreateDescriptors(test2, pointsANMS2, 16, 4, 16, ImageDescriptorMethods::kDescriptorTypeSquare, ImageDescriptorMethods::kDescriptorNormalization2Times);

	vector<CV_labs::Line> matchedDesc = ImageDescriptorMethods::DescriptorsMatching(descriptors1, descriptors2, ImageDescriptorMethods::kDescriptorsComparisonEuclid, ImageDescriptorMethods::kDescriptorsMatchingMutal);
	Image emptyImg(test1.GetColsNumber() + test2.GetColsNumber(), test1.GetRowsNumber() + test2.GetRowsNumber(), (uchar)0);
	emptyImg.InsertImage(test1, 0, 0); 
	emptyImg.InsertImage(test2, test1.GetColsNumber(), test1.GetRowsNumber());
	vector<CV_labs::Line> lines(matchedDesc.size());

	CV_labs::Line line;

	for (int i = 0; i < matchedDesc.size(); i++) {
		line.pointA = matchedDesc[i].pointA;

		line.pointB.x = matchedDesc[i].pointB.x + test1.GetColsNumber();
		line.pointB.y = matchedDesc[i].pointB.y + test1.GetRowsNumber();
		lines.push_back(line);
	}

	Mat colored = MatPlotLines(emptyImg, lines);
	
	double k = 1;
	resize(colored, colored, Size(k * test1.GetColsNumber(), k * test1.GetRowsNumber()));
	imshow("matched image", colored);
}

void lab5() {

	Mat mat1 = LoadImage();
	Image base1(mat1);
	base1.NormalizeImage();
	imshow("Base image 1", base1.GetMat());

	Image test1 = ImageFilters::ApplyFilter(base1, ImageFilters::GenerateGaussSeparableCore(1.6), ImageFilters::kInterpolateReflection);
	imshow("Base image 1 gauss", test1.GetMat());

	double angle;
	cout << "Enter angle: ";
	cin >> angle;

	int points;
	cout << "Enter points num: ";
	cin >> points;

	Image base2 = base1.GetRotatedImage(ImageFilters::PI() * angle / 180.0);
	base2.NormalizeImage();
	imshow("Base image 2", base2.GetMat());

	Image test2 = ImageFilters::ApplyFilter(base2, ImageFilters::GenerateGaussSeparableCore(1.6), ImageFilters::kInterpolateReflection);
	imshow("Base image 2 gauss", test2.GetMat());

	vector<CV_labs::Point> pointsANMS1 = ImageDetectors::Harris(test1, 5, 5, 0.03, points);
	imshow("Harris Points image 1 (ANMS)", MatPlotPoints(test1, pointsANMS1));

	vector<CV_labs::Point> pointsANMS2 = ImageDetectors::Harris(test2, 5, 5, 0.03, points);
	imshow("Harris Points image 2 (ANMS)", MatPlotPoints(test2, pointsANMS2));

	//testAlpha = 0;
	vector<CV_labs::Descriptor> descriptors1 = ImageDescriptorMethods::CreateDescriptors(test1, pointsANMS1, 16, 4, 8, ImageDescriptorMethods::kDescriptorTypeSquare, ImageDescriptorMethods::kDescriptorNormalization2Times);
	cout << "\n";
	//testAlpha = ImageFilters::PI() * angle / 180.0;
	vector<CV_labs::Descriptor> descriptors2 = ImageDescriptorMethods::CreateDescriptors(test2, pointsANMS2, 16, 4, 8, ImageDescriptorMethods::kDescriptorTypeSquare, ImageDescriptorMethods::kDescriptorNormalization2Times);
	cout << "\n";
	vector<CV_labs::Line> matchedDesc = ImageDescriptorMethods::DescriptorsMatching(descriptors1, descriptors2, ImageDescriptorMethods::kDescriptorsComparisonEuclid, ImageDescriptorMethods::kDescriptorsMatchingMutal);
	Image emptyImg(test1.GetColsNumber() + test2.GetColsNumber(), test1.GetRowsNumber() + test2.GetRowsNumber(), (uchar)0);
	emptyImg.InsertImage(test1, 0, 0);
	emptyImg.InsertImage(test2, test1.GetColsNumber(), test1.GetRowsNumber());
	vector<CV_labs::Line> lines(matchedDesc.size());

	CV_labs::Line line;

	for (int i = 0; i < matchedDesc.size(); i++) {
		line.pointA = matchedDesc[i].pointA;

		line.pointB.x = matchedDesc[i].pointB.x + test1.GetColsNumber();
		line.pointB.y = matchedDesc[i].pointB.y + test1.GetRowsNumber();
		lines.push_back(line);
	}

	Mat colored = MatPlotLines(emptyImg, lines);

	double k = 1;
	resize(colored, colored, Size(k * test1.GetColsNumber(), k * test1.GetRowsNumber()));
	imshow("matched image", colored);
}


//C:\Users\alist\Desktop\cv\Bikesgray.jpg
//C:\Users\alist\Desktop\cv\obj.jpg
//C:\Users\alist\Desktop\cv\star.jpg
//C:\Users\alist\Desktop\cv\Lena2.jpg
//C:\Users\alist\Desktop\cv\Lena2lheight.jpg
//C:\Users\alist\Desktop\cv\Lena2contrast.jpg
//C:\Users\alist\Desktop\cv\Lena2light.jpg
//C:\Users\alist\Desktop\cv\LenaRotate.jpg
//C:\Users\alist\Desktop\cv\LenaRotate2.jpg
//C:\Users\alist\Desktop\cv\Valve_original_(1).PNG
int main() {

	lab5();

	//Mat mat1 = LoadImage();
	//Image test1(mat1);
	//imshow("Base image 1", test1.GetMat());

	//Image emptyImg(test.GetRowsNumber() + test1.GetRowsNumber(), test.GetColsNumber() + test1.GetColsNumber(), 0);
	//
	////double* ddata = test.GetNormalizeDataF();
	////Image test1(test.GetRowsNumber(), test.GetColsNumber(), ddata,true);
	////imshow("Normalize image", test1.GetMat());

	//Image test2, test3, test4, test5;



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