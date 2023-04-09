#pragma once
#include<opencv2/opencv.hpp>
#include<vector>
#include<cmath>
#include<time.h>
using namespace cv;
using namespace std;
const double PI = acos(-1); // pi值
const int height = 256, width = 256;

struct Cpx // 定义一个复数结构体和复数运算法则
{
	double r, i;
	Cpx() : r(0), i(0) {} // 结构体的默认函数
	Cpx(double _r, double _i) : r(_r), i(_i) {} // 结构体的有参构造函数
};

//运算符重载


int ReverseBin(int a, int n);

Mat Resize_doubleline(Mat img);

void fft(vector<Cpx>& a, int lim, int opt);

void FFT2D(Cpx(*src)[width], Cpx(*dst)[width], int opt);

void ScaleMinMax(Mat src, Mat out);

void fftshift(Mat src, Mat out);

void fftshift(Cpx(*src)[width], Cpx(*dst)[width]);

void Mat2Cpx(Mat src, Cpx(*dst)[width]);

void Cpx2Mat(Cpx(*src)[width], Mat dst);

void get_AutoCor_from_img(Cpx(*src)[width], Cpx(*dst)[width]);

void get_FFTamp_from_AutoCor(Cpx(*src)[width], Cpx(*dst)[width]);

void sqrt(Cpx(*src)[width], Cpx(*dst)[width]);

void Cpx2MatDouble(Cpx(*src)[width], Mat dst);

void add_RandomPhase(Cpx(*src)[width], Cpx(*dst)[width]);

void ER_constraint(Cpx(*src)[width], Cpx(*dst)[width]);

void HIO_constraint(Cpx(*src_1)[width], Cpx(*src_2)[width], Cpx(*dst)[width], double beta);

void Update_Phase(Cpx(*FFTamp)[width], Cpx(*dst)[width]);

void CenterImg(Cpx(*src)[width], Cpx(*dst)[width]);