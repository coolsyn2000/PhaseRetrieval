#include<opencv2/opencv.hpp>
#include"sig_process.h"

Cpx operator + (Cpx a, Cpx b) { return Cpx(a.r + b.r, a.i + b.i); }
Cpx operator - (Cpx a, Cpx b) { return Cpx(a.r - b.r, a.i - b.i); }
Cpx operator * (Cpx a, Cpx b) { return Cpx(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
Cpx operator / (Cpx a, int b) { return Cpx(a.r * 1.0 / b, a.i * 1.0 / b); }


int ReverseBin(int a, int n)
	{
		int ret = 0;
		for (int i = 0; i < n; i++)
		{
			if (a & (1 << i)) ret |= (1 << (n - 1 - i));
		}
		return ret;
	}

void fft(vector<Cpx> &a, int lim, int opt) {
	int index;
	vector<Cpx> tempA(lim);
	for (int i = 0; i < lim; i++)
	{
		index = ReverseBin(i, log2(lim));
		tempA[i] = a[index];
	}

	vector<Cpx> WN(lim / 2);
	//生成WN表,避免重复计算
	for (int i = 0; i < lim / 2; i++)
	{
		WN[i] = Cpx(cos(2 * PI * i / lim), opt * -sin(2 * PI * i / lim));
	}

	//蝶形运算
	int Index0, Index1;
	Cpx temp;
	for (int steplenght = 2; steplenght <= lim; steplenght *= 2)
	{
		for (int step = 0; step < lim / steplenght; step++)
		{
			for (int i = 0; i < steplenght / 2; i++)
			{
				Index0 = steplenght * step + i;
				Index1 = steplenght * step + i + steplenght / 2;

				temp = tempA[Index1] * WN[lim / steplenght * i];
				tempA[Index1] = tempA[Index0] - temp;
				tempA[Index0] = tempA[Index0] + temp;
			}
		}
	}
	for (int i = 0; i < lim; i++)
	{
		if (opt == -1)
		{
			a[i] = tempA[i] / lim;
		}
		else
		{
			a[i] = tempA[i];
		}
	}
}

void FFT2D(Cpx(*src)[width], Cpx(*dst)[width], int opt)
{
	//第一遍fft
	for (int i = 0; i < height; i++)
	{
		vector<Cpx> tempData(src[i],src[i]+width);//行优化拷贝速度
		
		/*
		vector<Cpx> tempData(width);
		//获取每行数据
		for (int j = 0; j < width; j++)
		{
			tempData[j] = src[i][j];
		}
		*/
		
		//一维FFT
		fft(tempData, width, opt);
		//写入每行数据
		for (int j = 0; j < width; j++)
		{
			dst[i][j] = tempData[j];
		}
	}

	//第二遍fft
	for (int i = 0; i < width; i++)
	{
		//列拷贝优化
		vector<Cpx> tempData(height);
		//获取每列数据
		for (int j = 0; j < height; j++)
		{
			tempData[j] = dst[j][i];
		}
		//一维FFT
		fft(tempData, height, opt);
		//写入每列数据
		for (int j = 0; j < height; j++)
		{
			dst[j][i] = tempData[j];
		}
	}
}
void sqrt(Cpx(*src)[width], Cpx(*dst)[width])
{
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			dst[i][j].r = sqrt(src[i][j].r);
			dst[i][j].i = 0;

		}
	}
}

void get_AutoCor_from_img(Cpx(*src)[width], Cpx(*dst)[width])
{
	Cpx (*tempA)[width] = new Cpx[height][width];
	Cpx (*tempB)[width] = new Cpx[height][width];
	FFT2D(src, tempA, 1);
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			tempA[i][j].r= tempA[i][j].r * tempA[i][j].r + tempA[i][j].i * tempA[i][j].i;
			tempA[i][j].i = 0 ;

		}
	}
	FFT2D(tempA, tempB, -1);
	fftshift(tempB, dst);
	delete[] tempA;
	delete[] tempB;
}

void get_FFTamp_from_AutoCor(Cpx(*src)[width], Cpx(*dst)[width])
{
	Cpx(*tempA)[width] = new Cpx[height][width];
	Cpx(*tempB)[width] = new Cpx[height][width];
	fftshift(src, tempA);
	FFT2D(tempA, tempB,1);
	sqrt(tempB, dst);
	delete[]tempA;
	delete[]tempB;
}

void add_RandomPhase(Cpx(*src)[width], Cpx(*dst)[width])
{
	srand(time(NULL));
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			double N_rand = ((double)rand() / (RAND_MAX));
			//cout<<N_rand<<endl;
			dst[i][j].r = src[i][j].r*cos(2*PI*N_rand);
			dst[i][j].i = src[i][j].r*sin(2*PI*N_rand);
		}
	}
}

void ER_constraint(Cpx(*src)[width], Cpx(*dst)[width])
{
	double real;
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			dst[i][j].i = 0;
			real = src[i][j].r;
			if (real < 0) dst[i][j].r = 0;
			else dst[i][j].r = real;
		}
	}
}

void HIO_constraint(Cpx(*src_1)[width], Cpx(*src_2)[width], Cpx(*dst)[width], double beta)
{
	double real;
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			real = src_1[i][j].r;
			if (real< 0) dst[i][j].r = src_2[i][j].r - beta * real;
			else dst[i][j].r = real;
			dst[i][j].i = 0;
		}
	}
}

void Update_Phase(Cpx(*FFTamp)[width], Cpx(*dst)[width])
{
	double Amp, Pha;
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			Amp = sqrt(FFTamp[i][j].r * FFTamp[i][j].r + FFTamp[i][j].i * FFTamp[i][j].i);
			Pha = atan2(dst[i][j].i, dst[i][j].r);
			dst[i][j].r = Amp * cos(Pha);
			dst[i][j].i = Amp * sin(Pha);
		}
	}
}

void CenterImg(Cpx(*src)[width], Cpx(*dst)[width])
{
	int i_MAX = 0;
	int j_MAX = 0;
	//Cpx(*temp)[width] = new Cpx[height][width];
	//temp = src;
	double MAX=src[0][0].r;
	// 找最大值及索引
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			if (src[i][j].r > MAX)
			{
				MAX = src[i][j].r;
				i_MAX = i;
				j_MAX = j;
			}
		}
	}
	//矩阵shift操作
	int width_shift = (width / 2)- i_MAX;
	int height_shift = (height / 2)- j_MAX;

	//width_shift = ((width_shift % width) + width) % width;
	//height_shift = ((height_shift % width) + height) % height;
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			dst[i][j].r =src[(i - width_shift+width)%width][(j - (height_shift)+height)%height].r;
			//dst[i][j].r = src[(i + width_shift) % width][(j + height_shift) % height].r;
		}
	}

}

Mat Resize_doubleline(Mat img)
{
	int w = img.cols;
	int h = img.rows;
	Mat out(height, width, CV_8UC1);

	uchar* p = out.ptr<uchar>(0);//p是指向out第一行第一个元素的指针
	uchar* p2 = img.ptr<uchar>(0);//p2是指向img第一行第一个元素的指针

	int x_before, y_before;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++)
		{
			//寻找(x,y)在原图像中的对应像素
			float d_original_img_hnum = y * h / height;
			float d_original_img_wnum = x * w / width;

			//对应像素转整形，得到四领域中左上a点的像素坐标
			int i_original_img_hnum = d_original_img_hnum;
			int i_original_img_wnum = d_original_img_wnum;


			float distance_to_a_x = d_original_img_wnum - i_original_img_wnum;//在原图像中与a点的水平距离    
			float distance_to_a_y = d_original_img_hnum - i_original_img_hnum;//在原图像中与a点的垂直距离  

			//找到abcd四个坐标
			int point_a = i_original_img_hnum * w + i_original_img_wnum;
			int point_b = i_original_img_hnum * w + i_original_img_wnum + 1;
			int point_c = (i_original_img_hnum + 1) * w + i_original_img_wnum;
			int point_d = (i_original_img_hnum + 1) * w + i_original_img_wnum + 1;

			//边缘处理
			if (i_original_img_hnum == h - 1)
			{
				point_c = point_a;
				point_d = point_b;
			}
			if (i_original_img_wnum == w - 1)
			{
				point_b = point_a;
				point_d = point_c;
			}

			//计算双线性插值结果
			p[y * width + x] =
				p2[point_a] * (1 - distance_to_a_x) * (1 - distance_to_a_y) +
				p2[point_b] * distance_to_a_x * (1 - distance_to_a_y) +
				p2[point_c] * distance_to_a_y * (1 - distance_to_a_x) +
				p2[point_d] * distance_to_a_y * distance_to_a_x;

		}
	}
	return out;
}

void ScaleMinMax(Mat src, Mat out)
{
	int w = src.cols;
	int h = src.rows;
	double* p = src.ptr<double>(0);
	uchar* pOut = out.ptr<uchar>(0);
	double max = p[0];
	double min = p[0];

	//找minmax
	for (int i = 0; i < w * h; i++)
	{
		if (p[i] > max) max = p[i];
		if (p[i] < min) min = p[i];
	}

	double scale = 255.0 / (max - min);
	//fftshift
	for (int i = 0; i < w * h; i++)
	{
		pOut[i] = (uchar)((p[i] - min) * scale);
	}
}

void fftshift(Mat src, Mat out) {
	int w = src.cols;
	int h = src.rows;
	uchar* p = src.ptr<uchar>(0);
	uchar* pOut = out.ptr<uchar>(0);
	for (int i = 0; i < w * h; i++)
	{
		int j = i + w * h / 2 + w / 2;
		if (j > w * h) j = j - w * h;   //低频移至中间
		pOut[i] = (p[j] );
	}
}

void fftshift(Cpx(*src)[width], Cpx(*dst)[width]) {
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			if (i < width / 2 && j < height / 2) dst[i][j] = src[i + width / 2][j + height / 2];
			else if (i >= width / 2 && j >= height / 2) dst[i][j] = src[i - width / 2][j - height / 2];
			else if (i >= width / 2 && j < height / 2) dst[i][j] = src[i - width / 2][j + height / 2];
			else if (i < width / 2 && j >= height / 2) dst[i][j] = src[i + width / 2][j - height / 2];
			else break;
		}
	}
}

void Mat2Cpx(Mat src, Cpx(*dst)[width])
{
	//这里Mat里的数据得是unchar类型
	uchar* p = src.ptr<uchar>(0);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			dst[i][j] = Cpx(p[i * width + j], 0);
		}
	}
}

void Cpx2Mat(Cpx(*src)[width], Mat dst)
{
	//这里Mat里的数据得是unchar类型
	uchar* p = dst.ptr<uchar>(0);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			double g = sqrt(src[i][j].r * src[i][j].r);
			p[j + i * width] = (uchar)g;
		}
	}
}

void Cpx2MatDouble(Cpx(*src)[width], Mat dst)
{
	//这里Mat里的数据得是double类型
	double* p = dst.ptr<double>(0);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			double g = sqrt(src[i][j].r * src[i][j].r + src[i][j].i * src[i][j].i);
			//g = log(g + 1);  //转换为对数尺度
			p[j + i * width] = (double)g;
		}
	}
}