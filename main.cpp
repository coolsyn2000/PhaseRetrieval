#include<iostream>
#include"sig_process.h"


int main() {
	Mat img = imread("./test_sample/number4_256.jpg",0);
	imshow("read image", img);
	Cpx(*src)[width] = new Cpx[height][width];
	//src是指向含有width个Cpx结构体的数组的指针
	//用new创建二维数组，类型为结构体Cpx。 每个一维数组的类型为Cpx[width]

	Mat imgRez = Resize_doubleline(img);
	//imshow("resize", imgRez);
	Mat2Cpx(imgRez, src);

	//Cpx (*dst)[width] = new Cpx[height][width];
	//dst是傅里叶变换的内存空间

	Cpx(*xc)[width] = new Cpx[height][width];
	double t5 = getTickCount();
	get_AutoCor_from_img(src, xc);
	Cpx(*FFTamp)[width] = new Cpx[height][width];
	get_FFTamp_from_AutoCor(xc, FFTamp);
	Cpx(*FFTamp_rPha)[width] = new Cpx[height][width];
	//FFT2D(src, FFTamp, 1);
	add_RandomPhase(FFTamp,FFTamp_rPha);

	double t6 = getTickCount();
	double t5t6 = (t6 - t5) / getTickFrequency();
	std::cout << "autoCor耗时: " << t5t6 << "秒" << std::endl;

	Mat xc_mat = Mat::zeros(height, width, CV_64F);
	Mat xc_mat2 = Mat::zeros(height, width, CV_8UC1);
	Cpx2MatDouble(xc, xc_mat);
	ScaleMinMax(xc_mat, xc_mat2);
	imshow("xc", xc_mat2);
	
	const int N_HIO = 100;
	const int N_ER = 50;
	Cpx(*Obj)[width] = new Cpx[height][width];
	Cpx(*Obj_Cons)[width] = new Cpx[height][width];

	double ta;
	double tb;

	double beta = 2;
	//double ta = getTickCount();
	//HIO
	//ta = getTickCount();
	for (int i = 0; i < N_HIO; i++)
	{	
		FFT2D(FFTamp_rPha, Obj, -1);
		HIO_constraint(Obj, Obj_Cons,Obj_Cons,beta);
		FFT2D(Obj_Cons, FFTamp_rPha, 1);
		Update_Phase(FFTamp, FFTamp_rPha);
		beta = beta - (beta / N_HIO);
	}

	//ER
	
	for (int i = 0; i<N_ER; i++)
	{
		FFT2D(FFTamp_rPha, Obj, -1);
		ER_constraint(Obj, Obj_Cons);
		FFT2D(Obj_Cons, FFTamp_rPha, 1);
		Update_Phase(FFTamp, FFTamp_rPha);
	}
	
	Cpx(*Obj_Recons)[width] = new Cpx[height][width];
	ta = getTickCount();
	FFT2D(FFTamp_rPha, Obj, -1);
	tb = getTickCount();

	ER_constraint(Obj, Obj_Cons);
	CenterImg(Obj_Cons, Obj_Recons);
	double tatb = (tb - ta) / getTickFrequency();
	std::cout << "PRA耗时: " << tatb << "秒" << std::endl;

	//显示重建图像
	Mat out = Mat::zeros(height, width, CV_64F);
	Cpx2MatDouble(Obj_Recons, out);

	Mat out3 = Mat::zeros(height, width, CV_8UC1);
	ScaleMinMax(out, out3);
	imshow("相位恢复", out3);
	//imwrite("相位恢复.jpg", out3);

	waitKey(0);
}
