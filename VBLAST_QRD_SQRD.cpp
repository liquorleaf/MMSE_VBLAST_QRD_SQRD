#define _CRT_SECURE_NO_WARNINGS
#include"complex_matrix_operators.h"

//************************************辅助函数**********************************************

//获得一个均值为零，方差为1的高斯变量值；乘σ，加μ，可变为标准差σ，均值μ的高斯变量
//使用Box-Muller算法
//代码参考:https://blog.csdn.net/a616905919/article/details/45080539
double gaussianrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;
	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);	//由于涉及求对数，舍去S=0这个点情况，不影响连续随机变量的分布
		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);
	phase = 1 - phase;	//两个X都服从均值为0，方差为1的高斯分布，但每次选其中一个，下次选另一个
	return X;
}


//************************************产生信号、信道、噪声************************************

//本实验中的调制信号估算
Complex estimating_quantization_operation(Complex y)
{
	Complex result;
	//第一象限和正半实轴上的估算为(1+j)/sqrt(2)
	if ((y.re > 0) && (y.im >= 0)) { result.re = 1.0 / sqrt(2.0); result.im = 1.0 / sqrt(2.0); }
	//第二象限和正半虚轴上的估算为(-1+j)/sqrt(2)
	else if ((y.re <= 0) && (y.im > 0)) { result.re = -1.0 / sqrt(2.0); result.im = 1.0 / sqrt(2.0); }
	//第三象限和负半实轴上的估算为(-1-j)/sqrt(2)
	else if ((y.re < 0) && (y.im <= 0)) { result.re = -1.0 / sqrt(2.0); result.im = -1.0 / sqrt(2.0); }
	//第四象限和负半虚轴上的估算为(1-j)/sqrt(2)
	else if ((y.re >= 0) && (y.im < 0)) { result.re = 1.0 / sqrt(2.0); result.im = -1.0 / sqrt(2.0); }
	//原点上的估算为1+j
	else { result.re = 1.0 / sqrt(2.0); result.im = 1.0 / sqrt(2.0); }
	return result;
}

//产生均匀分布的随机01序列，注意输入参数nT是稍后调制后的发射天线数，故总比特数为nT的两倍
int* generate_bits(int nT)
{
	double* random_nums = new double[nT * 2];
	int* result = new int[nT * 2];
	for (int i = 0; i < nT * 2; i++)
	{
		random_nums[i] = ((double)rand() / RAND_MAX) * 2 - 1;
		if (random_nums[i] >= 0)result[i] = 1;
		else result[i] = 0;
	}
	delete[]random_nums;
	return result;
}

//把随机生成的bits调制，00到1+j，01到1-j，10到-1+j，11到-1-j，并且都/sqrt(2)
ComplexMatrix generate_signal(int* bits, int nT)
{
	ComplexMatrix signal(nT, 1, false);
	for (int i = 0; i < nT * 2; i++)
	{
		signal.c[i/2][0].re = bits[i] ? (-1 / sqrt(2.0)) : 1 / sqrt(2.0);
		i++;
		signal.c[i/2][0].im = bits[i] ? (-1 / sqrt(2.0)) : 1 / sqrt(2.0);
	}
	return signal;
}

//把还原出的信号翻译回bits
int* signal_to_bits(ComplexMatrix& c_hat)
{
	int* bits = new int[c_hat.lr * 2];
	for (int i = 0; i < c_hat.lr * 2; i++)
	{
		bits[i] = (abs(c_hat.c[i / 2][0].re + 1 / sqrt(2.0)) < 1e-16) ? 1 : 0;
		i++;
		bits[i] = (abs(c_hat.c[i / 2][0].im + 1 / sqrt(2.0)) < 1e-16) ? 1 : 0;
	}
	return bits;
}

//产生一个高斯随机变量信道矩阵（均值0，方差1，互相独立）
//根据本处使用的定义，实部虚部分别为（均值0，方差1）的独立高斯变量的复数为（均值0，方差1）的复高斯变量
ComplexMatrix generate_H(int nR, int nT)
{
	ComplexMatrix H(nR, nT, false);
	for(int i = 0; i < H.lr; i++)
		for (int j = 0; j < H.lc; j++)
		{
			H.c[i][j].re = gaussianrand();
			H.c[i][j].im = gaussianrand();
		}
	return H;
}

//产生一个高斯随机变量噪声向量（均值0，方差1，互相独立）
//根据本处使用的定义，实部虚部分别为（均值0，方差1）的独立高斯变量的复数为（均值0，方差1）的复高斯变量
ComplexMatrix generate_noise(int nR)
{
	ComplexMatrix niu(nR, 1, false);
	for (int i = 0; i < niu.lr; i++)
	{
		niu.c[i][0].re = gaussianrand();
		niu.c[i][0].im = gaussianrand();
	}
	return niu;
}

//************************************还原算法**********************************************

//直接乘伪逆算法
ComplexMatrix Pseudo_inverse(ComplexMatrix& x, ComplexMatrix& H)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix x_temp(x.lr, x.lc, x.is_real, x.c);	//本次信号处理中x的临时副本，防止原x被改
	ComplexMatrix c_hat(H.lc, 1, x.is_real);			//估算的还原出的发射信号
	ComplexMatrix G = H.Moore_Penrose_pseudo_inverse();	//H的伪逆
	x_temp = G * x_temp;
	for (int i = 0; i < H.lc; i++)c_hat.c[i][0] = estimating_quantization_operation(x_temp.c[i][0]);
	//释放内存
	//x_temp.clear(); G.clear();;
	return c_hat;
}

//V-BLAST算法
ComplexMatrix V_BLAST(ComplexMatrix& x, ComplexMatrix& H)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix x_temp(x.lr, x.lc, x.is_real, x.c);	//本次信号处理中x的临时副本，防止原x被改
	ComplexMatrix c_hat(H.lc, 1, x.is_real);			//估算的还原出的发射信号
	ComplexMatrix y(1, 1);								//还原出的发射信号和噪声之和
	Complex y_ = 0.0;									//还原出的发射信号和噪声之和
	ComplexMatrix G;									//H的伪逆
	ComplexMatrix g;									//存储g的行
	double g_2_norm;									//存储g的2-范数
	double g_2_norm_min = 0.0;							//存储G的2-范数最小的一行的2-范数
	int g_2_norm_min_row = 0;							//存储G的2-范数最小的一行的存储下标
	ComplexMatrix H_temp;								//存储H在迭代过程中每次减一列后的结果
	ComplexMatrix h;									//存储H_temp的列
	ComplexMatrix temp;									//计算中介
	//记录H的列是否已被处理，以便在新一轮中H去掉对应列
	int** flag1 = new int* [2];
	flag1[0] = new int[H.lc];
	flag1[1] = new int[H.lc];
	for (int i = 0; i < H.lc; i++)
	{
		flag1[0][i] = i;	//flag1第一行记录列下标
		flag1[1][i] = 0;	//flag1第二行记录当前列是否已被处理，初始为0(未处理)
	}
	int k = 0;	//记录flag2赋值到哪一列了

	//共循环nT次
	for (int i = 0; i < H.lc; i++)
	{
		//根据flag1拼出去掉已处理列的新H，和记录当前剩下的列对应的原下标的flag2
		int* flag2 = new int[H.lc - i];
		k = 0;
		for (int j = 0; j < H.lc; j++)
			if (!flag1[1][j])
			{
				flag2[k] = flag1[0][j];
				k++;
			}
		H_temp = H.get_column(flag2[0]);
		for (int j = 1; j < H.lc - i; j++)
		{
			h = H.get_column(flag2[j]);
			H_temp = H_temp.combine_columns(H_temp, h);
		}
		//求H_temp的伪逆G
		G = H_temp.Moore_Penrose_pseudo_inverse();
		//找出G中2-范数最小的一行（初始为第一行）
		g_2_norm_min_row = 0;
		g = G.get_row(0);
		g_2_norm_min = g.vector_2_norm();
		for (int j = 1; j < G.lr; j++)
		{
			g = G.get_row(j);
			g_2_norm = g.vector_2_norm();
			if (g_2_norm < g_2_norm_min)
			{
				g_2_norm_min_row = j;
				g_2_norm_min = g_2_norm;
			}
		}
		//取得G中2-范数最小的一行
		g = G.get_row(g_2_norm_min_row);
		//接收信号还原
		y = g * x_temp;
		//还原的信号用量化算子估算（经过设计，本实验测试运行范围内不会越界）
		c_hat.c[flag2[g_2_norm_min_row]][0] = estimating_quantization_operation(y.c[0][0]);
		//用估算的还原信号修正接收信号x_temp
		h = H_temp.get_column(g_2_norm_min_row);
		temp = h * c_hat.c[flag2[g_2_norm_min_row]][0];
		x_temp = x_temp - temp;
		//flag1中当前处理的一列标记为已处理
		flag1[1][flag2[g_2_norm_min_row]] = 1;
		//释放内存
		delete[]flag2;
	}
	//释放内存
	delete[]flag1[1];	delete[]flag1[0];	delete[]flag1;
	//x_temp.clear();	y.clear(); G.clear(); g.clear(); H_temp.clear(); h.clear(); temp.clear();
	return c_hat;
}

//V-BLAST分层误码率检测
ComplexMatrix V_BLAST_layer(ComplexMatrix& x, ComplexMatrix& H, int* bits)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix x_temp(x.lr, x.lc, x.is_real, x.c);	//本次信号处理中x的临时副本，防止原x被改
	ComplexMatrix c_hat(H.lc, 1, x.is_real);			//估算的还原出的发射信号
	ComplexMatrix y(1, 1);								//还原出的发射信号和噪声之和
	Complex y_ = 0.0;									//还原出的发射信号和噪声之和
	ComplexMatrix G;									//H的伪逆
	ComplexMatrix g;									//存储g的行
	double g_2_norm;									//存储g的2-范数
	double g_2_norm_min = 0.0;							//存储G的2-范数最小的一行的2-范数
	int g_2_norm_min_row = 0;							//存储G的2-范数最小的一行的存储下标
	ComplexMatrix H_temp;								//存储H在迭代过程中每次减一列后的结果
	ComplexMatrix h;									//存储H_temp的列
	ComplexMatrix temp;									//计算中介

	//记录H的列是否已被处理，以便在新一轮中H去掉对应列
	int** flag1 = new int* [3];
	flag1[0] = new int[H.lc];
	flag1[1] = new int[H.lc];
	flag1[2] = new int[H.lc];
	for (int i = 0; i < H.lc; i++)
	{
		flag1[0][i] = i;	//flag1第一行记录列下标
		flag1[1][i] = 0;	//flag1第二行记录当前列是否已被处理，初始为0(未处理)
		flag1[2][i] = 0;	//flag1第三行记录当前列是第几个被处理的
	}
	int k = 0;	//记录flag2赋值到哪一列了

	//共循环nT次
	for (int i = 0; i < H.lc; i++)
	{
		//根据flag1拼出去掉已处理列的新H，和记录当前剩下的列对应的原下标的flag2
		int* flag2 = new int[H.lc - i];
		k = 0;
		for (int j = 0; j < H.lc; j++)
			if (!flag1[1][j])
			{
				flag2[k] = flag1[0][j];
				k++;
			}
		H_temp = H.get_column(flag2[0]);
		for (int j = 1; j < H.lc - i; j++)
		{
			h = H.get_column(flag2[j]);
			H_temp = H_temp.combine_columns(H_temp, h);
		}
		//求H_temp的伪逆G
		G = H_temp.Moore_Penrose_pseudo_inverse();
		//找出G中2-范数最小的一行（初始为第一行）
		g_2_norm_min_row = 0;
		g = G.get_row(0);
		g_2_norm_min = g.vector_2_norm();
		for (int j = 1; j < G.lr; j++)
		{
			g = G.get_row(j);
			g_2_norm = g.vector_2_norm();
			if (g_2_norm < g_2_norm_min)
			{
				g_2_norm_min_row = j;
				g_2_norm_min = g_2_norm;
			}
		}
		//取得G中2-范数最小的一行
		g = G.get_row(g_2_norm_min_row);
		//接收信号还原
		y = g * x_temp;
		//还原的信号用量化算子估算（经过设计，本实验测试运行范围内不会越界）
		c_hat.c[flag2[g_2_norm_min_row]][0] = estimating_quantization_operation(y.c[0][0]);
		//用估算的还原信号修正接收信号x_temp
		h = H_temp.get_column(g_2_norm_min_row);
		temp = h * c_hat.c[flag2[g_2_norm_min_row]][0];
		x_temp = x_temp - temp;
		//flag1中当前处理的一列标记为已处理，并记录被处理顺序
		flag1[1][flag2[g_2_norm_min_row]] = 1;
		flag1[2][flag2[g_2_norm_min_row]] = i;
		//释放内存
		delete[]flag2;
	}
	//还原信号翻译回比特流
	int* bits_hat = signal_to_bits(c_hat);
	//按照记录的处理顺序执行分层误码率检测，最先为最高层
	ComplexMatrix BER_per_layer(1, H.lc, true);			//初始化分层误码率记录矩阵
	for (int i = 0; i < H.lc; i++)						//共nT层
	{
		int order = flag1[2][i];		//第i列的处理顺位(0开始)
		int layer = H.lc - order;		//第i列的层数(第0个处理对应第nT层，之后处理顺位每+1，层数-1)
		//获取第i列原始的2位比特
		int bit0 = bits[2 * i];
		int bit1 = bits[2 * i + 1];
		//获取第i列翻译回的2位比特
		int bit0_hat = bits_hat[2 * i];
		int bit1_hat = bits_hat[2 * i + 1];
		//第i列是第layer层，储存在误码率检测的第(layer-1)列
		//误码率：错误一个码加1，最后除以2
		BER_per_layer.c[0][layer - 1].re += (bit0 != bit0_hat);
		BER_per_layer.c[0][layer - 1].re += (bit1 != bit1_hat);
		BER_per_layer.c[0][layer - 1].re /= 2.0;
	}
	//释放内存
	delete[] flag1[2]; delete[]flag1[1]; delete[]flag1[0];	delete[]flag1;	delete[] bits_hat;
	//x_temp.clear();	y.clear(); G.clear(); g.clear(); H_temp.clear(); h.clear(); temp.clear();
	return BER_per_layer;
}

//无排序QRD算法
ComplexMatrix QRD(ComplexMatrix& x, ComplexMatrix& H)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix c_hat(H.lc, 1, x.is_real);		//估算的还原出的发射信号
	ComplexMatrix H_temp(H.lr, H.lc, H.is_real, H.c);	//存储H在本次V-BLAST计算中的临时复制
	ComplexMatrix y;		//存储接收信号x被Q修正后的信号
	ComplexMatrix Q_H;		//存储Q的共轭转置
	Complex temp;		//计算中介
	//对H_temp进行QR分解，之后H_temp中存储了Q
	ComplexMatrix R = H_temp.Gram_Schmidt_QR_modified__Complex(H_temp);
	//x根据Q变换为y
	Q_H = !H_temp;
	y = Q_H * x;
	//利用R的nT行和y还原出发射信号
	for (int i = c_hat.lr - 1; i >= 0; i--)
	{
		c_hat.c[i][0] = y.c[i][0];
		for (int j = c_hat.lr - 1; j > i; j--)
		{
			temp = R.c[i][j] * c_hat.c[j][0];
			c_hat.c[i][0] = c_hat.c[i][0] - temp;
		}
		c_hat.c[i][0] = c_hat.c[i][0] / R.c[i][i];
		c_hat.c[i][0] = estimating_quantization_operation(c_hat.c[i][0]);
	}
	//释放内存
	//H_temp.clear();	y.clear(); Q_H.clear(); R.clear();
	return c_hat;
}

//SQRD算法
ComplexMatrix SQRD(ComplexMatrix& x, ComplexMatrix& H)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix c_hat(H.lc, 1, x.is_real);		//估算的还原出的发射信号
	ComplexMatrix c_hat_temp(H.lc, 1, x.is_real);	//估算的还原出的发射信号的临时副本
	ComplexMatrix H_temp(H.lr, H.lc, H.is_real, H.c);	//存储H在本次V-BLAST计算中的临时复制
	ComplexMatrix y;		//存储接收信号x被Q修正后的信号
	ComplexMatrix Q_H;		//存储Q的共轭转置
	Complex temp;		//计算中介
	//初始化存储交换记录的向量
	int* S = new int[H.lc];
	for (int j = 0; j < H.lc; j++) S[j] = j;
	//对H_temp进行有排序的QR分解，之后H_temp中存储了Q
	ComplexMatrix R = H_temp.sorted_Gram_Schmidt_QR_modified__Complex(H_temp, S);
	//x根据Q变换为y
	Q_H = !H_temp;
	y = Q_H * x;
	//利用R的nT行和y还原出发射信号
	for (int i = H.lc - 1; i >= 0; i--)
	{
		c_hat.c[i][0] = y.c[i][0];
		for (int j = H.lc - 1; j > i; j--)
		{
			temp = R.c[i][j] * c_hat.c[j][0];
			c_hat.c[i][0] = c_hat.c[i][0] - temp;
		}
		c_hat.c[i][0] = c_hat.c[i][0] / R.c[i][i];
		c_hat.c[i][0] = estimating_quantization_operation(c_hat.c[i][0]);
	}
	//从S中还原发射信号c的原顺序（先复制一份，再按S顺序将复制品中的元素向结果向量中塞）
	for (int i = 0; i < H.lc; i++) c_hat_temp.c[i][0] = c_hat.c[i][0];
	for (int i = 0; i < H.lc; i++) c_hat.c[S[i]][0] = c_hat_temp.c[i][0];
	//释放内存
	delete[]S;
	//c_hat_temp.clear(); H_temp.clear();	y.clear(); Q_H.clear(); R.clear();
	return c_hat;
}

//************************************MMSE还原算法******************************************

//MMSE-线性(直接乘伪逆)算法
ComplexMatrix MMSE_Pseudo_inverse(ComplexMatrix& x, ComplexMatrix& H)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix x_temp(x.lr, x.lc, x.is_real, x.c);	//本次信号处理中x的临时副本，防止原x被改
	ComplexMatrix x0(H.lc, 1, x.is_real);				//用于拼接的临时零向量
	x_temp = x_temp.combine_rows(x_temp, x0);			//MMSE检测器所用接收向量
	ComplexMatrix c_hat(H.lc, 1, x.is_real);			//储存估算的还原出的发射信号

	ComplexMatrix H_temp(H);							//H的临时副本
	ComplexMatrix I0 = H.make_eyes(H.lc);				//单位阵
	H_temp = H_temp.combine_rows(H_temp, I0);			//MMSE检测器所用H矩阵
	ComplexMatrix G = H_temp.Moore_Penrose_pseudo_inverse();	//MMSE-H的伪逆

	x_temp = G * x_temp;
	for (int i = 0; i < H.lc; i++)c_hat.c[i][0] = estimating_quantization_operation(x_temp.c[i][0]);
	//释放内存
	//x_temp.clear(); G.clear();;
	return c_hat;
}

//MMSE-V-BLAST算法
ComplexMatrix MMSE_V_BLAST(ComplexMatrix& x, ComplexMatrix& H)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix x_temp(x.lr, x.lc, x.is_real, x.c);	//本次信号处理中x的临时副本，防止原x被改
	ComplexMatrix x0(H.lc, 1, x.is_real);				//用于拼接的临时零向量
	x_temp = x_temp.combine_rows(x_temp, x0);			//MMSE检测器所用接收向量
	ComplexMatrix c_hat(H.lc, 1, x.is_real);			//估算的还原出的发射信号

	ComplexMatrix I0 = H.make_eyes(H.lc);				//单位阵
	ComplexMatrix H_temp0 = H.combine_rows(H, I0);		//MMSE检测器所用信道矩阵

	ComplexMatrix y(1, 1);								//还原出的发射信号和噪声之和
	Complex y_ = 0.0;									//还原出的发射信号和噪声之和
	ComplexMatrix G;									//MMSE-H的伪逆
	ComplexMatrix g;									//存储g的行
	double g_2_norm;									//存储g的2-范数
	double g_2_norm_min = 0.0;							//存储G的2-范数最小的一行的2-范数
	int g_2_norm_min_row = 0;							//存储G的2-范数最小的一行的存储下标
	ComplexMatrix H_temp;								//存储H在迭代过程中每次减一列后的结果
	ComplexMatrix h;									//存储H_temp的列
	ComplexMatrix temp;									//计算中介
	//记录H的列是否已被处理，以便在新一轮中H去掉对应列
	int** flag1 = new int* [2];
	flag1[0] = new int[H_temp0.lc];
	flag1[1] = new int[H_temp0.lc];
	for (int i = 0; i < H_temp0.lc; i++)
	{
		flag1[0][i] = i;	//flag1第一行记录列下标
		flag1[1][i] = 0;	//flag1第二行记录当前列是否已被处理，初始为0(未处理)
	}
	int k = 0;	//记录flag2赋值到哪一列了

	//共循环nT次
	for (int i = 0; i < H_temp0.lc; i++)
	{
		//根据flag1拼出去掉已处理列的新H和去相应行的x，和记录当前剩下的列对应的原下标的flag2
		int* flag2 = new int[H_temp0.lc - i];
		k = 0;
		for (int j = 0; j < H_temp0.lc; j++)
			if (!flag1[1][j])
			{
				flag2[k] = flag1[0][j];
				k++;
			}
		H_temp = H_temp0.get_column(flag2[0]);
		for (int j = 1; j < H_temp0.lc - i; j++)
		{
			h = H_temp0.get_column(flag2[j]);
			H_temp = H_temp.combine_columns(H_temp, h);
		}
		//求H_temp的伪逆G
		G = H_temp.Moore_Penrose_pseudo_inverse();
		//找出G中2-范数最小的一行（初始为第一行）
		g_2_norm_min_row = 0;
		g = G.get_row(0);
		g_2_norm_min = g.vector_2_norm();
		for (int j = 1; j < G.lr; j++)
		{
			g = G.get_row(j);
			g_2_norm = g.vector_2_norm();
			if (g_2_norm < g_2_norm_min)
			{
				g_2_norm_min_row = j;
				g_2_norm_min = g_2_norm;
			}
		}
		//取得G中2-范数最小的一行
		g = G.get_row(g_2_norm_min_row);
		//接收信号还原
		y = g * x_temp;
		//还原的信号用量化算子估算（经过设计，本实验测试运行范围内不会越界）
		c_hat.c[flag2[g_2_norm_min_row]][0] = estimating_quantization_operation(y.c[0][0]);
		//用估算的还原信号修正接收信号x_temp
		h = H_temp.get_column(g_2_norm_min_row);
		temp = h * c_hat.c[flag2[g_2_norm_min_row]][0];
		x_temp = x_temp - temp;
		//flag1中当前处理的一列标记为已处理
		flag1[1][flag2[g_2_norm_min_row]] = 1;
		//释放内存
		delete[]flag2;
	}
	//释放内存
	delete[]flag1[1];	delete[]flag1[0];	delete[]flag1;
	//x_temp.clear();	y.clear(); G.clear(); g.clear(); H_temp.clear(); h.clear(); temp.clear();
	return c_hat;
}

//MMSE-V-BLAST分层误码率检测
ComplexMatrix MMSE_V_BLAST_layer(ComplexMatrix& x, ComplexMatrix& H, int* bits)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix x_temp(x.lr, x.lc, x.is_real, x.c);	//本次信号处理中x的临时副本，防止原x被改
	ComplexMatrix x0(H.lc, 1, x.is_real);				//用于拼接的临时零向量
	x_temp = x_temp.combine_rows(x_temp, x0);			//MMSE检测器所用接收向量
	ComplexMatrix c_hat(H.lc, 1, x.is_real);			//估算的还原出的发射信号

	ComplexMatrix I0 = H.make_eyes(H.lc);				//单位阵
	ComplexMatrix H_temp0 = H.combine_rows(H, I0);		//MMSE检测器所用信道矩阵

	ComplexMatrix y(1, 1);								//还原出的发射信号和噪声之和
	Complex y_ = 0.0;									//还原出的发射信号和噪声之和
	ComplexMatrix G;									//MMSE-H的伪逆
	ComplexMatrix g;									//存储g的行
	double g_2_norm;									//存储g的2-范数
	double g_2_norm_min = 0.0;							//存储G的2-范数最小的一行的2-范数
	int g_2_norm_min_row = 0;							//存储G的2-范数最小的一行的存储下标
	ComplexMatrix H_temp;								//存储H在迭代过程中每次减一列后的结果
	ComplexMatrix h;									//存储H_temp的列
	ComplexMatrix temp;									//计算中介
	//记录H的列是否已被处理，以便在新一轮中H去掉对应列
	int** flag1 = new int* [3];
	flag1[0] = new int[H_temp0.lc];
	flag1[1] = new int[H_temp0.lc];
	flag1[2] = new int[H_temp0.lc];
	for (int i = 0; i < H_temp0.lc; i++)
	{
		flag1[0][i] = i;	//flag1第一行记录列下标
		flag1[1][i] = 0;	//flag1第二行记录当前列是否已被处理，初始为0(未处理)
		flag1[2][i] = 0;	//flag1第三行记录当前列是第几个被处理的
	}
	int k = 0;	//记录flag2赋值到哪一列了

	//共循环nT次
	for (int i = 0; i < H_temp0.lc; i++)
	{
		//根据flag1拼出去掉已处理列的新H和去相应行的x，和记录当前剩下的列对应的原下标的flag2
		int* flag2 = new int[H_temp0.lc - i];
		k = 0;
		for (int j = 0; j < H_temp0.lc; j++)
			if (!flag1[1][j])
			{
				flag2[k] = flag1[0][j];
				k++;
			}
		H_temp = H_temp0.get_column(flag2[0]);
		for (int j = 1; j < H_temp0.lc - i; j++)
		{
			h = H_temp0.get_column(flag2[j]);
			H_temp = H_temp.combine_columns(H_temp, h);
		}
		//求H_temp的伪逆G
		G = H_temp.Moore_Penrose_pseudo_inverse();
		//找出G中2-范数最小的一行（初始为第一行）
		g_2_norm_min_row = 0;
		g = G.get_row(0);
		g_2_norm_min = g.vector_2_norm();
		for (int j = 1; j < G.lr; j++)
		{
			g = G.get_row(j);
			g_2_norm = g.vector_2_norm();
			if (g_2_norm < g_2_norm_min)
			{
				g_2_norm_min_row = j;
				g_2_norm_min = g_2_norm;
			}
		}
		//取得G中2-范数最小的一行
		g = G.get_row(g_2_norm_min_row);
		//接收信号还原
		y = g * x_temp;
		//还原的信号用量化算子估算（经过设计，本实验测试运行范围内不会越界）
		c_hat.c[flag2[g_2_norm_min_row]][0] = estimating_quantization_operation(y.c[0][0]);
		//用估算的还原信号修正接收信号x_temp
		h = H_temp.get_column(g_2_norm_min_row);
		temp = h * c_hat.c[flag2[g_2_norm_min_row]][0];
		x_temp = x_temp - temp;
		//flag1中当前处理的一列标记为已处理
		flag1[1][flag2[g_2_norm_min_row]] = 1;
		flag1[2][flag2[g_2_norm_min_row]] = i;
		//释放内存
		delete[]flag2;
	}
	//还原信号翻译回比特流
	int* bits_hat = signal_to_bits(c_hat);
	//按照记录的处理顺序执行分层误码率检测，最先为最高层
	ComplexMatrix BER_per_layer(1, H_temp0.lc, true);	//初始化分层误码率记录矩阵
	for (int i = 0; i < H_temp0.lc; i++)				//共nT层
	{
		int order = flag1[2][i];		//第i列的处理顺位(0开始)
		int layer = H.lc - order;		//第i列的层数(第0个处理对应第nT层，之后处理顺位每+1，层数-1)
		//获取第i列原始的2位比特
		int bit0 = bits[2 * i];
		int bit1 = bits[2 * i + 1];
		//获取第i列翻译回的2位比特
		int bit0_hat = bits_hat[2 * i];
		int bit1_hat = bits_hat[2 * i + 1];
		//第i列是第layer层，储存在误码率检测的第(layer-1)列
		//误码率：错误一个码加1，最后除以2
		BER_per_layer.c[0][layer - 1].re += (bit0 != bit0_hat);
		BER_per_layer.c[0][layer - 1].re += (bit1 != bit1_hat);
		BER_per_layer.c[0][layer - 1].re /= 2.0;
	}
	//释放内存
	delete[] flag1[2]; delete[]flag1[1]; delete[]flag1[0];	delete[]flag1;	delete[] bits_hat;
	return BER_per_layer;
}

//MMSE-无排序QRD算法
ComplexMatrix MMSE_QRD(ComplexMatrix& x, ComplexMatrix& H, double sigma_n)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix c_hat(H.lc, 1, x.is_real);		//估算的还原出的发射信号
	ComplexMatrix H_temp(H.lr, H.lc, H.is_real, H.c);	//存储H在本次V-BLAST计算中的临时复制
	ComplexMatrix I0 = H.make_eyes(H.lc);				//单位阵
	I0 = I0 * sigma_n;
	H_temp = H_temp.combine_rows(H_temp, I0);			//拼成扩充H矩阵
	ComplexMatrix y;		//存储接收信号x被Q修正后的信号
	ComplexMatrix Q_H;		//存储Q的共轭转置
	Complex temp;		//计算中介
	//对H_temp进行QR分解，之后H_temp中存储了Q
	ComplexMatrix R = H_temp.Gram_Schmidt_QR_modified__Complex(H_temp);
	//取出Q1，x根据Q变换为y
	H_temp = H_temp.get_rows(0, H.lr - 1);
	Q_H = !H_temp;
	y = Q_H * x;
	//利用R的nT行和y还原出发射信号
	for (int i = c_hat.lr - 1; i >= 0; i--)
	{
		c_hat.c[i][0] = y.c[i][0];
		for (int j = c_hat.lr - 1; j > i; j--)
		{
			temp = R.c[i][j] * c_hat.c[j][0];
			c_hat.c[i][0] = c_hat.c[i][0] - temp;
		}
		c_hat.c[i][0] = c_hat.c[i][0] / R.c[i][i];
		c_hat.c[i][0] = estimating_quantization_operation(c_hat.c[i][0]);
	}
	//释放内存
	//H_temp.clear();	y.clear(); Q_H.clear(); R.clear();
	return c_hat;
}

//MMSE-SQRD算法
ComplexMatrix MMSE_SQRD(ComplexMatrix& x, ComplexMatrix& H, double sigma_n)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix c_hat(H.lc, 1, x.is_real);		//估算的还原出的发射信号
	ComplexMatrix qi;								//提取Q的列
	ComplexMatrix qk;								//提取Q的列
	ComplexMatrix qiH;								//储存Q的列的共轭转置

	double q_2_norm_min;			//存储Q的2-范数最小的一列的2-范数
	int q_2_norm_min_column;		//存储Q的2-范数最小的一列的存储下标
	double r_ii;					//临时存储r_ii相关计算
	ComplexMatrix rik;				//临时存储r_ik相关计算

	ComplexMatrix temp;				//计算中间变量
															//QR分解有问题
	//初始化QR分解的R
	ComplexMatrix R(H.lc, H.lc, x.is_real);
	//初始化QR分解的扩充Q为扩充H
	ComplexMatrix Q(H);
	ComplexMatrix I0 = Q.make_eyes(H.lc);
	I0 = I0 * sigma_n;
	Q = Q.combine_rows(Q, I0);
	//初始化存储交换记录的向量
	int* p = new int[Q.lc];
	//初始化存储Q的列的2-范数的向量
	double* norms = new double[Q.lc];
	for (int i = 0; i < Q.lc; i++)
	{
		p[i] = i;
		norms[i] = 0.0;
	}
	for (int i = 0; i < Q.lc; i++)
	{
		qi = Q.get_column(i);
		norms[i] = qi.vector_2_norm();
		norms[i] *= norms[i];
	}
	for (int i = 0; i < Q.lc; i++)
	{
		//找出Q剩下列中2-范数最小的一列下标
		q_2_norm_min_column = i;
		q_2_norm_min = norms[i];
		for (int k = i + 1; k < Q.lc; k++)
			if (norms[k] < q_2_norm_min)
				q_2_norm_min_column = k;
		//交换第i列和第q_2_norm_min_column列
		R.exchange_column(i, q_2_norm_min_column);
		std::swap(p[i], p[q_2_norm_min_column]);
		std::swap(norms[i], norms[q_2_norm_min_column]);
		//交换扩充Q的第i列和第q_2_norm_min_column列的前(nR+i-1)行
		Q.exchange_some_rows_of_column(i, q_2_norm_min_column, 0, x.lr + i - 1);
		//R的对角元
		R.c[i][i] = sqrt(norms[i]);
		//Q的列标准化
		r_ii = R.c[i][i].re;
		for (int j = 0; j < Q.lr; j++)
			Q.c[j][i] = Q.c[j][i] / r_ii;
		qi = Q.get_column(i);
		qiH = !qi;
		//从后续中减去分量，以正交化
		for (int k = i + 1; k < Q.lc; k++)
		{
			qk = Q.get_column(k);
			//R的非对角上三角元
			rik = qiH * qk;
			R.c[i][k] = rik.c[0][0];
			//更新Q的第k列
			temp = rik.c[0][0] * qi;
			qk = qk - temp;
			for (int j = 0; j < Q.lr; j++)
				Q.c[j][k] = qk.c[j][0];
			//更新norms的第k个
			norms[k] -= pow(rik.c[0][0].modulus(), 2);
		}
	}

	//取出Q1，x根据Q变换为y
	Q = Q.get_rows(0, x.lr - 1);
	ComplexMatrix Q_H = !Q;
	ComplexMatrix y = Q_H * x;
	//利用R的nT行和y还原出发射信号
	Complex temp1;						//计算中间变量
	for (int i = c_hat.lr - 1; i >= 0; i--)
	{
		c_hat.c[i][0] = y.c[i][0];
		for (int j = c_hat.lr - 1; j > i; j--)
		{
			temp1 = R.c[i][j] * c_hat.c[j][0];
			c_hat.c[i][0] = c_hat.c[i][0] - temp1;
		}
		c_hat.c[i][0] = c_hat.c[i][0] / R.c[i][i];
		c_hat.c[i][0] = estimating_quantization_operation(c_hat.c[i][0]);
	}
	//从S中还原发射信号c的原顺序（先复制一份，再按S顺序将复制品中的元素向结果向量中塞）
	//（经过设计，本实验测试运行范围内不会越界）
	ComplexMatrix c_hat_temp(c_hat);
	for (int i = 0; i < Q.lc; i++) c_hat.c[p[i]][0] = c_hat_temp.c[i][0];
	//释放内存
	delete[]p; delete[]norms;
	return c_hat;
}

//MMSE-SQRD-PSA算法
ComplexMatrix MMSE_SQRD_PSA(ComplexMatrix& x, ComplexMatrix& H, double sigma_n)
{
	//接收x:nR * 1		//信道H:nR * nT		//还原出c_hat:nT * 1
	ComplexMatrix c_hat(H.lc, 1, x.is_real);		//估算的还原出的发射信号
	ComplexMatrix qi;								//提取Q的列
	ComplexMatrix qk;								//提取Q的列
	ComplexMatrix qiH;								//储存Q的列的共轭转置

	double q_2_norm_min;			//存储Q的2-范数最小的一列的2-范数
	int q_2_norm_min_column;		//存储Q的2-范数最小的一列的存储下标
	double r_ii;					//临时存储r_ii相关计算
	ComplexMatrix rik;				//临时存储r_ik相关计算

	ComplexMatrix temp;				//计算中间变量
															//QR分解有问题
	//初始化QR分解的R
	ComplexMatrix R(H.lc, H.lc, x.is_real);
	//初始化QR分解的扩充Q为扩充H
	ComplexMatrix Q(H);
	ComplexMatrix I0 = Q.make_eyes(H.lc);
	I0 = I0 * sigma_n;
	Q = Q.combine_rows(Q, I0);
	//初始化存储交换记录的向量
	int* p = new int[Q.lc];
	//初始化存储Q的列的2-范数的向量
	double* norms = new double[Q.lc];
	for (int i = 0; i < Q.lc; i++)
	{
		p[i] = i;
		norms[i] = 0.0;
	}
	for (int i = 0; i < Q.lc; i++)
	{
		qi = Q.get_column(i);
		norms[i] = qi.vector_2_norm();
		norms[i] *= norms[i];
	}
	for (int i = 0; i < Q.lc; i++)
	{
		//找出Q剩下列中2-范数最小的一列下标
		q_2_norm_min_column = i;
		q_2_norm_min = norms[i];
		for (int k = i + 1; k < Q.lc; k++)
			if (norms[k] < q_2_norm_min)
				q_2_norm_min_column = k;
		//交换第i列和第q_2_norm_min_column列
		R.exchange_column(i, q_2_norm_min_column);
		std::swap(p[i], p[q_2_norm_min_column]);
		std::swap(norms[i], norms[q_2_norm_min_column]);
		//交换扩充Q的第i列和第q_2_norm_min_column列的前(nR+i-1)行
		Q.exchange_some_rows_of_column(i, q_2_norm_min_column, 0, x.lr + i - 1);
		//R的对角元
		R.c[i][i] = sqrt(norms[i]);
		//Q的列标准化
		r_ii = R.c[i][i].re;
		for (int j = 0; j < Q.lr; j++)
			Q.c[j][i] = Q.c[j][i] / r_ii;
		qi = Q.get_column(i);
		qiH = !qi;
		//从后续中减去分量，以正交化
		for (int k = i + 1; k < Q.lc; k++)
		{
			qk = Q.get_column(k);
			//R的非对角上三角元
			rik = qiH * qk;
			R.c[i][k] = rik.c[0][0];
			//更新Q的第k列
			temp = rik.c[0][0] * qi;
			qk = qk - temp;
			for (int j = 0; j < Q.lr; j++)
				Q.c[j][k] = qk.c[j][0];
			//更新norms的第k个
			norms[k] -= pow(rik.c[0][0].modulus(), 2);
		}
	}

	//*********************************PSA*********************************
	ComplexMatrix Q1 = Q.get_rows(0, x.lr - 1);		//取出Q1
	ComplexMatrix Q2 = Q.get_rows(x.lr, Q.lr - 1);	//取出Q2
	ComplexMatrix Householder;						//储存Householder变换阵
	ComplexMatrix a, u, aH, uH, w1, w2;				//计算Householder变换阵的中间变量
	int k_min = H.lc;								//寻找误差最小行
	int ki = 0;										//寻找误差最小行
	double* error_l = new double[H.lc];				//储存误差error

	for (int i = H.lc - 1; i >= 1; i--)
	{
		//Q2每行属于下三角的部分的2-范数视为误差
		for (int l = 0; l <= i; l++)
		{
			temp = Q2.get_sub_matrix(l, l, 0, i);
			error_l[l] = pow(temp.vector_2_norm(), 2);
		}
		//找到上述误差最小行放到当前子矩阵最后
		ki = i;
		for (int k = i - 1; k >= 0; k--)
			if (error_l[k] < error_l[ki])
				ki = k;
		if (k_min > ki)k_min = ki;
		//如果误差最小行不在当前子矩阵最后一行，则和最后交换，并用Householder矩阵把Q2变回上三角
		if (ki < i)
		{
			Q2.exchange_row(i, ki);
			std::swap(p[i], p[ki]);
		}
		if (k_min < i)
		{
			//计算Q2的第i行的第k_min到第(i-1)列的Householder变换阵
			a = Q2.get_sub_matrix(i, i, k_min, i);
			ComplexMatrix en(1, a.lc, a.is_real);
			en.c[0][en.lc - 1].re = a.vector_2_norm();
			u = a - en;
			double temp1 = 1 / u.vector_2_norm();
			u = u * temp1;
			uH = !u;
			aH = !a;
			w1 = u * aH;
			w2 = a * uH;
			Complex w = w1.c[0][0] / w2.c[0][0];
			w.re += 1.0;
			temp = uH * u;
			temp = w * temp;
			Householder = Householder.make_eyes(a.lc);
			Householder = Householder - temp;
			//Q2的第0到i行的第k_min到第(i-1)列右乘Householder变换阵
			//Q1的第k_min到第(i-1)列右乘Householder变换阵
			//先分别挖出子矩阵，拼起来右乘，再填回去
			ComplexMatrix Q1__ = Q1.get_sub_matrix(0, Q1.lr - 1, k_min, i);
			ComplexMatrix Q2__ = Q2.get_sub_matrix(0, i, k_min, i);
			temp = temp.combine_rows(Q1__, Q2__);
			temp = temp * Householder;
			for (int jj = 0; jj < temp.lc; jj++)
			{
				for (int ii = 0; ii < Q1.lr; ii++)
					Q1.c[ii][jj + k_min] = temp.c[ii][jj];
				for (int ii = 0; ii <= i; ii++)
					Q2.c[ii][jj + k_min] = temp.c[ii + Q1.lr][jj];
			}
		}
	}
	double temp2 = 1 / sigma_n;
	ComplexMatrix Q2_inverse = Q2.square_inverse();
	R = temp2 * Q2_inverse;

	//*********************************还原信号*****************************
	//x根据Q变换为y
	ComplexMatrix Q1_H = !Q1;
	ComplexMatrix y = Q1_H * x;
	//利用R的nT行和y还原出发射信号
	Complex temp1;						//计算中间变量
	for (int i = c_hat.lr - 1; i >= 0; i--)
	{
		c_hat.c[i][0] = y.c[i][0];
		for (int j = c_hat.lr - 1; j > i; j--)
		{
			temp1 = R.c[i][j] * c_hat.c[j][0];
			c_hat.c[i][0] = c_hat.c[i][0] - temp1;
		}
		c_hat.c[i][0] = c_hat.c[i][0] / R.c[i][i];
		c_hat.c[i][0] = estimating_quantization_operation(c_hat.c[i][0]);
	}
	//从S中还原发射信号c的原顺序（先复制一份，再按S顺序将复制品中的元素向结果向量中塞）
	//（经过设计，本实验测试运行范围内不会越界）
	ComplexMatrix c_hat_temp(c_hat);
	for (int i = 0; i < Q.lc; i++) c_hat.c[p[i]][0] = c_hat_temp.c[i][0];

	//释放内存
	delete[]p; delete[]norms; delete[]error_l;
	return c_hat;
}

//************************************统计并向文件输出结果************************************

//定义文件输出流
ofstream dataOut;
//比较两个比特流，得出误码率
double BER(int* bits, int* bits_hat, int nT)
{
	int error = 0;
	for (int i = 0; i < nT * 2; i++)
		if (bits[i] != bits_hat[i])
			error++;
	double BER = double(error) / (double(nT) * 2.0);
	return BER;
}

//总体误码率测试：指定发射和接收天线数（输入参数先接收，后发射），产生n个H信道
//每次用这n个H跑一个采样点，随后改变信噪比，再跑下一次same_SNR_n_H
//存储误码率采用一个行向量，分别记录每种方法，由一个信噪比和一种方法可确定该方法在此信噪比下的平均误码率
//采样点的数据都输出到txt文件中以便MATLAB画图处理
//一个采样点：指定发射和接收天线数，输入一个信噪比（功率比）
//n个H以该信噪比跑一次one_H_nn，返回n个H用每种方法时的平均误码率，即每种方法在该点的误码率BER
//一个采样点的一部分：指定发射和接收天线数（在H的维度信息里体现），同一个H
//每种方法“同步”跑nn次（源比特流和噪声都不同），返回每种方法的平均误码率
void BER_versus_SNR(int nR, int nT)
{
	int n = 100;						//H个数
	cout << "H个数："; cin >> n;
	double SNR = 1.0;					//信噪功率比
	double SNR_sqrt = sqrt(SNR);		//信噪幅度比
	int nn = 100;						//运行次数
	cout << "每个H运行次数："; cin >> nn;

	string filename;						//根据nR和nT命名当前测试数据文件
	char nR_1[5] = { '_','n','R','_','\0' }, nR_2[16] = { 0 }, nT_1[4] = { 'n','T','_','\0' }, nT_2[16] = { 0 };
	_itoa(nR, nR_2, 10); _itoa(nT, nT_2, 10);
	filename = string(nT_1) + string(nT_2) + string(nR_1) + string(nR_2) + ".txt";
	dataOut.open(filename.c_str());

	//产生信道
	ComplexMatrix** H_n = new ComplexMatrix * [n];
	for (int i = 0; i < n; i++)
		H_n[i] = new ComplexMatrix(generate_H(nR, nT));
	//运行仿真
	for (int k = 0; k < 10; k++)		//前十个采样点SNR每次增大1，最小1开始
	{
		cout << "开始运行第" << k + 1 << "个采样点...\n";
		ComplexMatrix AVG_BER(1, 4, true);	//平均误码率初始化
		SNR = 1.0 + double(k);
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix BERs(1, 4, true);	//储存四种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//四种方法从接收信号还原，并翻译出比特流
				c_hat_0 = Pseudo_inverse(x, H_temp);	bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = V_BLAST(x, H_temp);			bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = QRD(x, H_temp);				bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = SQRD(x, H_temp);				bits_hat_3 = signal_to_bits(c_hat_3);
				//四个翻译比特流和源比特流比较，更新误码率
				BERs.c[0][0].re += BER(bits, bits_hat_0, H_temp.lc);
				BERs.c[0][1].re += BER(bits, bits_hat_1, H_temp.lc);
				BERs.c[0][2].re += BER(bits, bits_hat_2, H_temp.lc);
				BERs.c[0][3].re += BER(bits, bits_hat_3, H_temp.lc);
				//释放内存
				delete[]bits;
				delete[]bits_hat_0; delete[]bits_hat_1;
				delete[]bits_hat_2; delete[]bits_hat_3;
				c.clear(); niu.clear(); x.clear();
				c_hat_0.clear(); c_hat_1.clear(); c_hat_2.clear(); c_hat_3.clear();
			}
			//四个误码率各求平均值
			for (int i = 0; i < 4; i++)BERs.c[0][i].re /= nn;
			for (int i = 0; i < 4; i++)AVG_BER.c[0][i].re += BERs.c[0][i].re;
			//释放内存
			//H_temp.clear();
		}
		//误码率求平均值并输出
		for (int i = 0; i < 4; i++)AVG_BER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_BER;
		//释放内存
		//AVG_BER.clear();
	}
	for (int k = 0; k < 2; k++)			//后2个采样点SNR每次增大10，最大20结束,因为发现测试精度内后面四种方法都为0
	{
		cout << "开始运行第" << k + 11 << "个采样点...\n";
		ComplexMatrix AVG_BER(1, 4, true);	//平均误码率初始化
		SNR = 20.0 + double(k) * 10;
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix BERs(1, 4, true);	//储存四种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//四种方法从接收信号还原，并翻译出比特流
				c_hat_0 = Pseudo_inverse(x, H_temp);	bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = V_BLAST(x, H_temp);			bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = QRD(x, H_temp);				bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = SQRD(x, H_temp);				bits_hat_3 = signal_to_bits(c_hat_3);
				//四个翻译比特流和源比特流比较，更新误码率
				BERs.c[0][0].re += BER(bits, bits_hat_0, H_temp.lc);
				BERs.c[0][1].re += BER(bits, bits_hat_1, H_temp.lc);
				BERs.c[0][2].re += BER(bits, bits_hat_2, H_temp.lc);
				BERs.c[0][3].re += BER(bits, bits_hat_3, H_temp.lc);
				//释放内存
				delete[]bits;
				delete[]bits_hat_0; delete[]bits_hat_1;
				delete[]bits_hat_2; delete[]bits_hat_3;
				c.clear(); niu.clear(); x.clear();
				c_hat_0.clear(); c_hat_1.clear(); c_hat_2.clear(); c_hat_3.clear();
			}
			//四个误码率各求平均值
			for (int i = 0; i < 4; i++)BERs.c[0][i].re /= nn;
			for (int i = 0; i < 4; i++)AVG_BER.c[0][i].re += BERs.c[0][i].re;
			//释放内存
			//H_temp.clear();
		}
		//误码率求平均值并输出
		for (int i = 0; i < 4; i++)AVG_BER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_BER;
		//释放内存
		//AVG_BER.clear();
	}
}

//MMSE总体误码率测试
void MMSE_BER_versus_SNR(int nR, int nT)
{
	int n = 100;						//H个数
	cout << "H个数："; cin >> n;
	double SNR = 1.0;					//信噪功率比
	double SNR_sqrt = sqrt(SNR);		//信噪幅度比
	int nn = 100;						//运行次数
	cout << "每个H运行次数："; cin >> nn;

	string filename;						//根据nR和nT命名当前测试数据文件
	char nR_1[5] = { '_','n','R','_','\0' }, nR_2[16] = { 0 }, nT_1[4] = { 'n','T','_','\0' }, nT_2[16] = { 0 };
	_itoa(nR, nR_2, 10); _itoa(nT, nT_2, 10);
	filename = string(nT_1) + string(nT_2) + string(nR_1) + string(nR_2) + "_MMSE.txt";
	dataOut.open(filename.c_str());

	//产生信道
	ComplexMatrix** H_n = new ComplexMatrix * [n];
	for (int i = 0; i < n; i++)
		H_n[i] = new ComplexMatrix(generate_H(nR, nT));
	//运行仿真
	for (int k = 0; k < 10; k++)		//前10个采样点SNR每次增大1，最小1开始
	{
		cout << "开始运行第" << k + 1 << "个采样点...\n";
		ComplexMatrix AVG_BER(1, 5, true);	//平均误码率初始化
		SNR = 1.0 + double(k);
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix BERs(1, 5, true);	//储存每种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_4;			//储存SQRD_PSA从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits_hat_4;				//储存SQRD_PSA从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//从接收信号还原，并翻译出比特流
				c_hat_0 = MMSE_Pseudo_inverse(x, H_temp);	bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = MMSE_V_BLAST(x, H_temp);			bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = MMSE_QRD(x, H_temp, 1.0);			bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = MMSE_SQRD(x, H_temp, 1.0);		bits_hat_3 = signal_to_bits(c_hat_3);
				c_hat_4 = MMSE_SQRD_PSA(x, H_temp, 1.0);	bits_hat_4 = signal_to_bits(c_hat_4);
				//翻译比特流和源比特流比较，更新误码率
				BERs.c[0][0].re += BER(bits, bits_hat_0, H_temp.lc);
				BERs.c[0][1].re += BER(bits, bits_hat_1, H_temp.lc);
				BERs.c[0][2].re += BER(bits, bits_hat_2, H_temp.lc);
				BERs.c[0][3].re += BER(bits, bits_hat_3, H_temp.lc);
				BERs.c[0][4].re += BER(bits, bits_hat_4, H_temp.lc);
				//释放内存
				delete[]bits;		delete[]bits_hat_4;
				delete[]bits_hat_0; delete[]bits_hat_1;
				delete[]bits_hat_2; delete[]bits_hat_3;
			}
			//误码率各求平均值
			for (int i = 0; i < 5; i++)BERs.c[0][i].re /= nn;
			for (int i = 0; i < 5; i++)AVG_BER.c[0][i].re += BERs.c[0][i].re;
		}
		//误码率求平均值并输出
		for (int i = 0; i < 5; i++)AVG_BER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_BER;
	}
	for (int k = 0; k < 9; k++)			//后9个采样点SNR每次增大10，最大100结束
	{
		cout << "开始运行第" << k + 11 << "个采样点...\n";
		ComplexMatrix AVG_BER(1, 5, true);	//平均误码率初始化
		SNR = 20.0 + double(k) * 10;
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix BERs(1, 5, true);	//储存每种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_4;			//储存SQRD_PSA从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits_hat_4;				//储存SQRD_PSA从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//从接收信号还原，并翻译出比特流
				c_hat_0 = MMSE_Pseudo_inverse(x, H_temp);	bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = MMSE_V_BLAST(x, H_temp);			bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = MMSE_QRD(x, H_temp, 1.0);			bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = MMSE_SQRD(x, H_temp, 1.0);		bits_hat_3 = signal_to_bits(c_hat_3);
				c_hat_4 = MMSE_SQRD_PSA(x, H_temp, 1.0);	bits_hat_4 = signal_to_bits(c_hat_4);
				//翻译比特流和源比特流比较，更新误码率
				BERs.c[0][0].re += BER(bits, bits_hat_0, H_temp.lc);
				BERs.c[0][1].re += BER(bits, bits_hat_1, H_temp.lc);
				BERs.c[0][2].re += BER(bits, bits_hat_2, H_temp.lc);
				BERs.c[0][3].re += BER(bits, bits_hat_3, H_temp.lc);
				BERs.c[0][4].re += BER(bits, bits_hat_4, H_temp.lc);
				//释放内存
				delete[]bits;		delete[]bits_hat_4;
				delete[]bits_hat_0; delete[]bits_hat_1;
				delete[]bits_hat_2; delete[]bits_hat_3;
			}
			//误码率各求平均值
			for (int i = 0; i < 5; i++)BERs.c[0][i].re /= nn;
			for (int i = 0; i < 5; i++)AVG_BER.c[0][i].re += BERs.c[0][i].re;
		}
		//误码率求平均值并输出
		for (int i = 0; i < 5; i++)AVG_BER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_BER;
	}
}

//ZF和MMSE总体误码率测试
void ZF_MMSE_BER_versus_SNR(int nR, int nT)
{
	int n = 100;						//H个数
	cout << "H个数："; cin >> n;
	double SNR = 1.0;					//信噪功率比
	double SNR_sqrt = sqrt(SNR);		//信噪幅度比
	int nn = 100;						//运行次数
	cout << "每个H运行次数："; cin >> nn;

	string filename;						//根据nR和nT命名当前测试数据文件
	char nR_1[5] = { '_','n','R','_','\0' }, nR_2[16] = { 0 }, nT_1[4] = { 'n','T','_','\0' }, nT_2[16] = { 0 };
	_itoa(nR, nR_2, 10); _itoa(nT, nT_2, 10);
	filename = string(nT_1) + string(nT_2) + string(nR_1) + string(nR_2) + "_ZF_MMSE.txt";
	dataOut.open(filename.c_str());

	//产生信道
	ComplexMatrix** H_n = new ComplexMatrix * [n];
	for (int i = 0; i < n; i++)
		H_n[i] = new ComplexMatrix(generate_H(nR, nT));
	//运行仿真
	for (int k = 0; k < 10; k++)		//前10个采样点SNR每次增大1，最小1开始
	{
		cout << "开始运行第" << k + 1 << "个采样点...\n";
		ComplexMatrix AVG_BER(1, 9, true);	//平均误码率初始化
		SNR = 1.0 + double(k);
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix BERs(1, 9, true);		//储存每种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_4;			//储存MMSE线性(直接乘伪逆)从接收信号中还原的发射信号
				ComplexMatrix c_hat_5;			//储存MMSE-V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_6;			//储存MMSE-无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_7;			//储存MMSE-SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_8;			//储存MMSE-SQRD_PSA从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits_hat_4;				//储存MMSE线性(直接乘伪逆)从还原信号翻译出的比特流
				int* bits_hat_5;				//储存MMSE-V-BLAST从还原信号翻译出的比特流
				int* bits_hat_6;				//储存MMSE-无排序QRD从还原信号翻译出的比特流
				int* bits_hat_7;				//储存MMSE-SQRD从还原信号翻译出的比特流
				int* bits_hat_8;				//储存MMSE-SQRD_PSA从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//从接收信号还原，并翻译出比特流
				c_hat_0 = Pseudo_inverse(x, H_temp);		bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = V_BLAST(x, H_temp);				bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = QRD(x, H_temp);					bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = SQRD(x, H_temp);					bits_hat_3 = signal_to_bits(c_hat_3);
				c_hat_4 = MMSE_Pseudo_inverse(x, H_temp);	bits_hat_4 = signal_to_bits(c_hat_4);
				c_hat_5 = MMSE_V_BLAST(x, H_temp);			bits_hat_5 = signal_to_bits(c_hat_5);
				c_hat_6 = MMSE_QRD(x, H_temp, 1.0);			bits_hat_6 = signal_to_bits(c_hat_6);
				c_hat_7 = MMSE_SQRD(x, H_temp, 1.0);		bits_hat_7 = signal_to_bits(c_hat_7);
				c_hat_8 = MMSE_SQRD_PSA(x, H_temp, 1.0);	bits_hat_8 = signal_to_bits(c_hat_8);
				//翻译比特流和源比特流比较，更新误码率
				BERs.c[0][0].re += BER(bits, bits_hat_0, H_temp.lc);
				BERs.c[0][1].re += BER(bits, bits_hat_1, H_temp.lc);
				BERs.c[0][2].re += BER(bits, bits_hat_2, H_temp.lc);
				BERs.c[0][3].re += BER(bits, bits_hat_3, H_temp.lc);
				BERs.c[0][4].re += BER(bits, bits_hat_4, H_temp.lc);
				BERs.c[0][5].re += BER(bits, bits_hat_5, H_temp.lc);
				BERs.c[0][6].re += BER(bits, bits_hat_6, H_temp.lc);
				BERs.c[0][7].re += BER(bits, bits_hat_7, H_temp.lc);
				BERs.c[0][8].re += BER(bits, bits_hat_8, H_temp.lc);
				//释放内存
				delete[]bits;		delete[]bits_hat_0;
				delete[]bits_hat_1; delete[]bits_hat_2;
				delete[]bits_hat_3; delete[]bits_hat_4;
				delete[]bits_hat_5; delete[]bits_hat_6;
				delete[]bits_hat_7; delete[]bits_hat_8;
			}
			//误码率各求平均值
			for (int i = 0; i < 9; i++)BERs.c[0][i].re /= nn;
			for (int i = 0; i < 9; i++)AVG_BER.c[0][i].re += BERs.c[0][i].re;
		}
		//误码率求平均值并输出
		for (int i = 0; i < 9; i++)AVG_BER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_BER;
	}
	for (int k = 0; k < 9; k++)			//后9个采样点SNR每次增大10，最大100结束
	{
		cout << "开始运行第" << k + 11 << "个采样点...\n";
		ComplexMatrix AVG_BER(1, 9, true);	//平均误码率初始化
		SNR = 20.0 + double(k) * 10;
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix BERs(1, 9, true);		//储存每种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_4;			//储存MMSE线性(直接乘伪逆)从接收信号中还原的发射信号
				ComplexMatrix c_hat_5;			//储存MMSE-V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_6;			//储存MMSE-无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_7;			//储存MMSE-SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_8;			//储存MMSE-SQRD_PSA从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits_hat_4;				//储存MMSE线性(直接乘伪逆)从还原信号翻译出的比特流
				int* bits_hat_5;				//储存MMSE-V-BLAST从还原信号翻译出的比特流
				int* bits_hat_6;				//储存MMSE-无排序QRD从还原信号翻译出的比特流
				int* bits_hat_7;				//储存MMSE-SQRD从还原信号翻译出的比特流
				int* bits_hat_8;				//储存MMSE-SQRD_PSA从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//从接收信号还原，并翻译出比特流
				c_hat_0 = Pseudo_inverse(x, H_temp);		bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = V_BLAST(x, H_temp);				bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = QRD(x, H_temp);					bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = SQRD(x, H_temp);					bits_hat_3 = signal_to_bits(c_hat_3);
				c_hat_4 = MMSE_Pseudo_inverse(x, H_temp);	bits_hat_4 = signal_to_bits(c_hat_4);
				c_hat_5 = MMSE_V_BLAST(x, H_temp);			bits_hat_5 = signal_to_bits(c_hat_5);
				c_hat_6 = MMSE_QRD(x, H_temp, 1.0);			bits_hat_6 = signal_to_bits(c_hat_6);
				c_hat_7 = MMSE_SQRD(x, H_temp, 1.0);		bits_hat_7 = signal_to_bits(c_hat_7);
				c_hat_8 = MMSE_SQRD_PSA(x, H_temp, 1.0);	bits_hat_8 = signal_to_bits(c_hat_8);
				//翻译比特流和源比特流比较，更新误码率
				BERs.c[0][0].re += BER(bits, bits_hat_0, H_temp.lc);
				BERs.c[0][1].re += BER(bits, bits_hat_1, H_temp.lc);
				BERs.c[0][2].re += BER(bits, bits_hat_2, H_temp.lc);
				BERs.c[0][3].re += BER(bits, bits_hat_3, H_temp.lc);
				BERs.c[0][4].re += BER(bits, bits_hat_4, H_temp.lc);
				BERs.c[0][5].re += BER(bits, bits_hat_5, H_temp.lc);
				BERs.c[0][6].re += BER(bits, bits_hat_6, H_temp.lc);
				BERs.c[0][7].re += BER(bits, bits_hat_7, H_temp.lc);
				BERs.c[0][8].re += BER(bits, bits_hat_8, H_temp.lc);
				//释放内存
				delete[]bits;		delete[]bits_hat_0;
				delete[]bits_hat_1; delete[]bits_hat_2;
				delete[]bits_hat_3; delete[]bits_hat_4;
				delete[]bits_hat_5; delete[]bits_hat_6;
				delete[]bits_hat_7; delete[]bits_hat_8;
			}
			//误码率各求平均值
			for (int i = 0; i < 9; i++)BERs.c[0][i].re /= nn;
			for (int i = 0; i < 9; i++)AVG_BER.c[0][i].re += BERs.c[0][i].re;
		}
		//误码率求平均值并输出
		for (int i = 0; i < 9; i++)AVG_BER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_BER;
	}
}

//V-BLAST和MMSE-V-BLAST分层误码率测试
void BER_versus_SNR_per_layer(int nR, int nT)
{
	int n = 100;						//H个数
	cout << "H个数："; cin >> n;
	double SNR = 1.0;					//信噪功率比
	double SNR_sqrt = sqrt(SNR);		//信噪幅度比
	int nn = 100;						//运行次数
	cout << "每个H运行次数："; cin >> nn;

	string filename;						//根据nR和nT命名当前测试数据文件
	char nR_1[5] = { '_','n','R','_','\0' }, nR_2[16] = { 0 }, nT_1[4] = { 'n','T','_','\0' }, nT_2[16] = { 0 };
	_itoa(nR, nR_2, 10); _itoa(nT, nT_2, 10);
	filename = string(nT_1) + string(nT_2) + string(nR_1) + string(nR_2) + "_V-BLAST_layers.txt";
	dataOut.open(filename.c_str());

	//产生信道
	ComplexMatrix** H_n = new ComplexMatrix * [n];
	for (int i = 0; i < n; i++)
		H_n[i] = new ComplexMatrix(generate_H(nR, nT));
	//运行仿真
	for (int k = 0; k < 10; k++)		//前10个采样点SNR每次增大1，最小1开始
	{
		cout << "开始运行第" << k + 1 << "个采样点...\n";
		ComplexMatrix AVG_BER1(1, 4, true);	//V-BLAST平均分层误码率初始化
		ComplexMatrix AVG_BER2(1, 4, true);	//MMSE-V-BLAST平均分层误码率初始化
		SNR = 1.0 + double(k);
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix BER_layers1;		//储存V-BLAST分层误码率
				ComplexMatrix BER_layers2;		//储存MMSE-V-BLAST分层误码率
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//获取新的分层误码率
				BER_layers2 = MMSE_V_BLAST_layer(x, H_temp, bits);
				BER_layers1 = V_BLAST_layer(x, H_temp, bits);
				//更新平均误码率
				AVG_BER1 = AVG_BER1 + BER_layers1;
				AVG_BER2 = AVG_BER2 + BER_layers2;
				//释放内存
				delete[]bits;
			}
		}
		//误码率求平均值并输出
		for (int i = 0; i < AVG_BER1.lc; i++) AVG_BER1.c[0][i].re /= (double(n) * double(nn));
		for (int i = 0; i < AVG_BER2.lc; i++) AVG_BER2.c[0][i].re /= (double(n) * double(nn));
		ComplexMatrix AVG_BER = AVG_BER1.combine_columns(AVG_BER1, AVG_BER2);
		dataOut << SNR << " " << AVG_BER;
	}
	for (int k = 0; k < 9; k++)			//后九个采样点SNR每次增大10，最大100结束
	{
		cout << "开始运行第" << k + 11 << "个采样点...\n";
		ComplexMatrix AVG_BER1(1, 4, true);	//V-BLAST平均分层误码率初始化
		ComplexMatrix AVG_BER2(1, 4, true);	//MMSE-V-BLAST平均分层误码率初始化
		SNR = 20.0 + double(k) * 10;
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix BER_layers1;		//储存V-BLAST分层误码率
				ComplexMatrix BER_layers2;		//储存MMSE-V-BLAST分层误码率
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//获取新的分层误码率
				BER_layers2 = MMSE_V_BLAST_layer(x, H_temp, bits);
				BER_layers1 = V_BLAST_layer(x, H_temp, bits);
				//更新平均误码率
				AVG_BER1 = AVG_BER1 + BER_layers1;
				AVG_BER2 = AVG_BER2 + BER_layers2;
				//释放内存
				delete[]bits;
			}
		}
		//误码率求平均值并输出
		for (int i = 0; i < AVG_BER1.lc; i++) AVG_BER1.c[0][i].re /= (double(n) * double(nn));
		for (int i = 0; i < AVG_BER2.lc; i++) AVG_BER2.c[0][i].re /= (double(n) * double(nn));
		ComplexMatrix AVG_BER = AVG_BER1.combine_columns(AVG_BER1, AVG_BER2);
		dataOut << SNR << " " << AVG_BER;
	}
}

//MMSE总体误帧率测试
void MMSE_FER_versus_SNR(int nR, int nT)
{
	int n = 100;						//H个数
	cout << "H个数："; cin >> n;
	double SNR = 1.0;					//信噪功率比
	double SNR_sqrt = sqrt(SNR);		//信噪幅度比
	int nn = 100;						//运行次数
	cout << "每个H运行次数："; cin >> nn;

	string filename;						//根据nR和nT命名当前测试数据文件
	char nR_1[5] = { '_','n','R','_','\0' }, nR_2[16] = { 0 }, nT_1[4] = { 'n','T','_','\0' }, nT_2[16] = { 0 };
	_itoa(nR, nR_2, 10); _itoa(nT, nT_2, 10);
	filename = string(nT_1) + string(nT_2) + string(nR_1) + string(nR_2) + "_MMSE_FER.txt";
	dataOut.open(filename.c_str());

	//产生信道
	ComplexMatrix** H_n = new ComplexMatrix * [n];
	for (int i = 0; i < n; i++)
		H_n[i] = new ComplexMatrix(generate_H(nR, nT));
	//运行仿真
	for (int k = 0; k < 10; k++)		//前10个采样点SNR每次增大1，最小1开始
	{
		cout << "开始运行第" << k + 1 << "个采样点...\n";
		ComplexMatrix AVG_FER(1, 5, true);	//平均误帧率初始化
		SNR = 1.0 + double(k);
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix FERs(1, 5, true);	//储存每种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_4;			//储存SQRD_PSA从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits_hat_4;				//储存SQRD_PSA从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//从接收信号还原，并翻译出比特流
				c_hat_0 = MMSE_Pseudo_inverse(x, H_temp);	bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = MMSE_V_BLAST(x, H_temp);			bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = MMSE_QRD(x, H_temp, 1.0);			bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = MMSE_SQRD(x, H_temp, 1.0);		bits_hat_3 = signal_to_bits(c_hat_3);
				c_hat_4 = MMSE_SQRD_PSA(x, H_temp, 1.0);	bits_hat_4 = signal_to_bits(c_hat_4);
				//翻译比特流和源比特流比较，更新误帧率
				bool* Frame_error = new bool[5];
				for (int i = 0; i < 5; i++)Frame_error[i] = false;
				for (int i = 0; i < nT * 2; i++)
				{
					if (bits[i] != bits_hat_0[i]) Frame_error[0] = true;
					if (bits[i] != bits_hat_1[i]) Frame_error[1] = true;
					if (bits[i] != bits_hat_2[i]) Frame_error[2] = true;
					if (bits[i] != bits_hat_3[i]) Frame_error[3] = true;
					if (bits[i] != bits_hat_4[i]) Frame_error[4] = true;
				}
				for (int i = 0; i < 5; i++)
					if (Frame_error[i]) FERs.c[0][i].re += 1.0;
				//释放内存
				delete[]bits;		delete[]Frame_error;
				delete[]bits_hat_0; delete[]bits_hat_1;
				delete[]bits_hat_2; delete[]bits_hat_3;
				delete[]bits_hat_4;
			}
			//误帧率各求平均值
			for (int i = 0; i < 5; i++)FERs.c[0][i].re /= nn;
			for (int i = 0; i < 5; i++)AVG_FER.c[0][i].re += FERs.c[0][i].re;
		}
		//误帧率求平均值并输出
		for (int i = 0; i < 5; i++)AVG_FER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_FER;
	}
	for (int k = 0; k < 9; k++)			//后9个采样点SNR每次增大10，最大100结束
	{
		cout << "开始运行第" << k + 11 << "个采样点...\n";
		ComplexMatrix AVG_FER(1, 5, true);	//平均误帧率初始化
		SNR = 20.0 + double(k) * 10;
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix FERs(1, 5, true);	//储存每种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_4;			//储存SQRD_PSA从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits_hat_4;				//储存SQRD_PSA从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//从接收信号还原，并翻译出比特流
				c_hat_0 = MMSE_Pseudo_inverse(x, H_temp);	bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = MMSE_V_BLAST(x, H_temp);			bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = MMSE_QRD(x, H_temp, 1.0);			bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = MMSE_SQRD(x, H_temp, 1.0);		bits_hat_3 = signal_to_bits(c_hat_3);
				c_hat_4 = MMSE_SQRD_PSA(x, H_temp, 1.0);	bits_hat_4 = signal_to_bits(c_hat_4);
				//翻译比特流和源比特流比较，更新误帧率
				bool* Frame_error = new bool[5];
				for (int i = 0; i < 5; i++)Frame_error[i] = false;
				for (int i = 0; i < nT * 2; i++)
				{
					if (bits[i] != bits_hat_0[i]) Frame_error[0] = true;
					if (bits[i] != bits_hat_1[i]) Frame_error[1] = true;
					if (bits[i] != bits_hat_2[i]) Frame_error[2] = true;
					if (bits[i] != bits_hat_3[i]) Frame_error[3] = true;
					if (bits[i] != bits_hat_4[i]) Frame_error[4] = true;
				}
				for (int i = 0; i < 5; i++)
					if (Frame_error[i]) FERs.c[0][i].re += 1.0;
				//释放内存
				delete[]bits;		delete[]Frame_error;
				delete[]bits_hat_0; delete[]bits_hat_1;
				delete[]bits_hat_2; delete[]bits_hat_3;
				delete[]bits_hat_4;
			}
			//误帧率各求平均值
			for (int i = 0; i < 5; i++)FERs.c[0][i].re /= nn;
			for (int i = 0; i < 5; i++)AVG_FER.c[0][i].re += FERs.c[0][i].re;
		}
		//误帧率求平均值并输出
		for (int i = 0; i < 5; i++)AVG_FER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_FER;
	}
	for (int k = 0; k < 9; k++)			//后9个采样点SNR每次增大100，最大1000结束
	{
		cout << "开始运行第" << k + 20 << "个采样点...\n";
		ComplexMatrix AVG_FER(1, 5, true);	//平均误帧率初始化
		SNR = 200.0 + double(k) * 100;
		for (int j = 0; j < n; j++)
		{
			SNR_sqrt = sqrt(SNR);
			ComplexMatrix H_temp(*(H_n[j]));
			H_temp = H_temp * SNR_sqrt;
			ComplexMatrix FERs(1, 5, true);	//储存每种方法对一个H的误码率，初始化
			for (int l = 0; l < nn; l++)
			{
				ComplexMatrix c_hat_0;			//储存直接乘伪逆从接收信号中还原的发射信号
				ComplexMatrix c_hat_1;			//储存V-BLAST从接收信号中还原的发射信号
				ComplexMatrix c_hat_2;			//储存无排序QRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_3;			//储存SQRD从接收信号中还原的发射信号
				ComplexMatrix c_hat_4;			//储存SQRD_PSA从接收信号中还原的发射信号
				int* bits_hat_0;				//储存直接乘伪逆从还原信号翻译出的比特流
				int* bits_hat_1;				//储存V-BLAST从还原信号翻译出的比特流
				int* bits_hat_2;				//储存无排序QRD从还原信号翻译出的比特流
				int* bits_hat_3;				//储存SQRD从还原信号翻译出的比特流
				int* bits_hat_4;				//储存SQRD_PSA从还原信号翻译出的比特流
				int* bits = generate_bits(H_temp.lc);					//产生源比特流
				ComplexMatrix c = generate_signal(bits, H_temp.lc);		//产生发射信号
				ComplexMatrix niu = generate_noise(H_temp.lr);			//产生噪声
				ComplexMatrix x = H_temp * c; x = x + niu;				//产生接收信号
				//从接收信号还原，并翻译出比特流
				c_hat_0 = MMSE_Pseudo_inverse(x, H_temp);	bits_hat_0 = signal_to_bits(c_hat_0);
				c_hat_1 = MMSE_V_BLAST(x, H_temp);			bits_hat_1 = signal_to_bits(c_hat_1);
				c_hat_2 = MMSE_QRD(x, H_temp, 1.0);			bits_hat_2 = signal_to_bits(c_hat_2);
				c_hat_3 = MMSE_SQRD(x, H_temp, 1.0);		bits_hat_3 = signal_to_bits(c_hat_3);
				c_hat_4 = MMSE_SQRD_PSA(x, H_temp, 1.0);	bits_hat_4 = signal_to_bits(c_hat_4);
				//翻译比特流和源比特流比较，更新误帧率
				bool* Frame_error = new bool[5];
				for (int i = 0; i < 5; i++)Frame_error[i] = false;
				for (int i = 0; i < nT * 2; i++)
				{
					if (bits[i] != bits_hat_0[i]) Frame_error[0] = true;
					if (bits[i] != bits_hat_1[i]) Frame_error[1] = true;
					if (bits[i] != bits_hat_2[i]) Frame_error[2] = true;
					if (bits[i] != bits_hat_3[i]) Frame_error[3] = true;
					if (bits[i] != bits_hat_4[i]) Frame_error[4] = true;
				}
				for (int i = 0; i < 5; i++)
					if (Frame_error[i]) FERs.c[0][i].re += 1.0;
				//释放内存
				delete[]bits;		delete[]Frame_error;
				delete[]bits_hat_0; delete[]bits_hat_1;
				delete[]bits_hat_2; delete[]bits_hat_3;
				delete[]bits_hat_4;
			}
			//误帧率各求平均值
			for (int i = 0; i < 5; i++)FERs.c[0][i].re /= nn;
			for (int i = 0; i < 5; i++)AVG_FER.c[0][i].re += FERs.c[0][i].re;
		}
		//误帧率求平均值并输出
		for (int i = 0; i < 5; i++)AVG_FER.c[0][i].re /= n;
		dataOut << SNR << " " << AVG_FER;
	}
}

//************************************主程序************************************************

int main()
{
	srand(int(time(NULL)));	//更新随机数种子
	int nT, nR;
	
	//*******************************(1)个别测试1*******************************************
	nT = 4;		//发射天线数
	nR = 4;		//接收天线数
	cout << "(1)个别测试1\n";
	//产生bits并显示
	int* bits = generate_bits(nT);
	cout << "比特流 ="; for (int i = 0; i < nT * 2; i++) cout << bits[i]; cout << endl;
	//调制成发射信号，产生信道，产生噪声，全部显示
	ComplexMatrix c = generate_signal(bits, nT);
	cout << "发射信号c = \n" << c << endl;
	ComplexMatrix H = generate_H(nR, nT);
	cout << "信道矩阵H = \n" << H << endl;
	ComplexMatrix niu = generate_noise(nR);
	cout << "噪声ν = \n" << niu << endl;
	//产生接收信号
	ComplexMatrix x = H * c;
	x = x + niu;
	//九种算法尝试恢复接收信号为发射信号，显示结果；并还原为bits，显示结果
	ComplexMatrix result = Pseudo_inverse(x, H);  int* result_bits = signal_to_bits(result);
	cout << "直接乘伪逆算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;
	
	result = V_BLAST(x, H);  result_bits = signal_to_bits(result); 
	cout << "V-BLAST算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;

	result = QRD(x, H);	result_bits = signal_to_bits(result); 
	cout << "无排序QRD算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;

	result = SQRD(x, H); result_bits = signal_to_bits(result); 
	cout << "SQRD算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;

	result = MMSE_Pseudo_inverse(x, H); result_bits = signal_to_bits(result);
	cout << "MMSE_直接乘伪逆算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;

	result = MMSE_V_BLAST(x, H); result_bits = signal_to_bits(result);
	cout << "MMSE_V-BLAST算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;

	result = MMSE_QRD(x, H, 1.0); result_bits = signal_to_bits(result);
	cout << "MMSE_QRD算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;

	result = MMSE_SQRD(x, H, 1.0); result_bits = signal_to_bits(result);
	cout << "MMSE_SQRD算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;

	result = MMSE_SQRD_PSA(x, H, 1.0); result_bits = signal_to_bits(result);
	cout << "MMSE_SQRD_PSA算法结果 = \n" << result << endl << "结果比特流 =";
	for (int i = 0; i < nT * 2; i++) cout << result_bits[i]; cout << endl << endl;

	//*******************************(2)蒙特卡洛仿真测试************************************
	int input;
	cout << "选择：\n2.ZF总体误码率测试；3.V-BLAST和MMSE-V-BLAST分层误码率测试；4.MMSE总体误帧率测试；5.ZF和MMSE总体误码率测试：";
	cin >> input;
	cout << "输入接收天线数nR："; cin >> nR;
	cout << "输入发射天线数nT："; cin >> nT;
	if (input == 2)
	{
		BER_versus_SNR(nR, nT);
		cout << "测试完成\n";
	}
	else if (input == 3)
	{
		BER_versus_SNR_per_layer(nR, nT);
		cout << "测试完成\n";
	}
	else if (input == 4)
	{
		MMSE_FER_versus_SNR(nR, nT);
		cout << "测试完成\n";
	}
	else if (input == 5)
	{
		ZF_MMSE_BER_versus_SNR(nR, nT);
		cout << "测试完成\n";
	}
	return 0;
}
