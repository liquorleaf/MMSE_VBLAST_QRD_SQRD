#include<iostream>
#include<iomanip>
#include"complex_matrix_operators.h"
#include"complex_operators.h"
//ofstream dataOut;			//输出到文件，不用时注释掉

//构造函数
ComplexMatrix::ComplexMatrix(int m, int n, bool flag, Complex** incm)
{
	is_real = flag;
	lr = m;		//行数
	lc = n;		//列数
	//创建矩阵内存空间
	c = new Complex*[m];
	for (int i = 0; i < m; i++)
		c[i] = new Complex[n];
	//为矩阵赋值
	if (incm == NULL)	//没有输入矩阵时
	{
	
	
	}
	else				//输入矩阵时
	{
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				c[i][j] = incm[i][j];
	}
}
//构造函数，从二维double数组转换来
ComplexMatrix::ComplexMatrix(double** a, int m, int n)
{
	is_real = 1;//一定是实数矩阵
	lr = m;		//行数
	lc = n;		//列数
	//创建矩阵内存空间
	c = new Complex * [m];
	for (int i = 0; i < m; i++)
		c[i] = new Complex[n];
	//为矩阵赋值
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c[i][j].re = a[i][j];
}
//复制构造函数
ComplexMatrix::ComplexMatrix(const ComplexMatrix& A)
{
	lr = A.lr;
	lc = A.lc;
	is_real = A.is_real;
	//创建矩阵内存空间
	c = new Complex * [lr];
	for (int i = 0; i < lr; i++)
		c[i] = new Complex[lc];
	//为矩阵赋值
	if (A.c == NULL)	//输入矩阵为空时
		c = NULL;
	else				//输入矩阵时
	{
		for (int i = 0; i < lr; i++)
			for (int j = 0; j < lc; j++)
				c[i][j] = A.c[i][j];
	}
}
//析构函数
ComplexMatrix::~ComplexMatrix() { clear(); }
//清除函数
void ComplexMatrix::clear()
{
	if (c != NULL)
	{
		for (int i = 0; i < lr; i++)
			delete[] c[i];
		delete[] c;
	}
	c = NULL;
	lr = lc = 0;
}
//初始化为特殊矩阵
//单位阵
ComplexMatrix ComplexMatrix::make_eyes(int n)
{
	ComplexMatrix result(n, n);
	for (int i = 0; i < n; i++) result.c[i][i].re = 1.0;
	return result;
}

//矩阵加号+
ComplexMatrix operator + (ComplexMatrix& cm1, ComplexMatrix& cm2)
{
	if (cm1.lr != cm2.lr || cm1.lc != cm2.lc)	//维度不匹配则不能运算
	{
		cout << "维度不匹配，不能执行矩阵加法" << endl;
		system("pause");
	}
	ComplexMatrix cm3(cm1.lr, cm1.lc, (cm1.is_real&&cm2.is_real));//两个都为实数矩阵则返回实数矩阵
	for (int i = 0; i < cm3.lr; i++)
		for (int j = 0; j < cm3.lc; j++)
			cm3.c[i][j] = cm1.c[i][j] + cm2.c[i][j];
	return cm3;
}
//矩阵和二维数组加号+
ComplexMatrix operator + (ComplexMatrix& cm1, double**& d1)
{
	ComplexMatrix cm2(d1, cm1.lr, cm1.lc);
	return cm1 + cm2;
}
ComplexMatrix operator + (double**& d1, ComplexMatrix& cm1)
{
	ComplexMatrix cm2(d1, cm1.lr, cm1.lc);
	return cm1 + cm2;
}

//矩阵减号-
ComplexMatrix operator - (ComplexMatrix& cm1, ComplexMatrix& cm2)
{
	if (cm1.lr != cm2.lr || cm1.lc != cm2.lc)	//维度不匹配则不能运算
	{
		cout << "维度不匹配，不能执行矩阵减法" << endl;
		system("pause");
	}
	ComplexMatrix cm3(cm1.lr, cm1.lc, (cm1.is_real && cm2.is_real));//两个都为实数矩阵则返回实数矩阵
	for (int i = 0; i < cm3.lr; i++)
		for (int j = 0; j < cm3.lc; j++)
			cm3.c[i][j] = cm1.c[i][j] - cm2.c[i][j];
	return cm3;
}
//矩阵和二维数组减号-
ComplexMatrix operator - (ComplexMatrix& cm1, double**& d1)
{
	ComplexMatrix cm2(d1, cm1.lr, cm1.lc);
	return cm1 - cm2;
}
ComplexMatrix operator - (double**& d1, ComplexMatrix& cm1)
{
	ComplexMatrix cm2(d1, cm1.lr, cm1.lc);
	return cm1 - cm2;
}

//矩阵乘号*
ComplexMatrix operator * (ComplexMatrix& cm1, ComplexMatrix& cm2)
{
	if (cm1.lc != cm2.lr)	//维度不匹配则不能运算
	{
		cout << "维度不匹配，不能执行矩阵乘法" << endl;
		system("pause");
	}
	ComplexMatrix cm3(cm1.lr, cm2.lc, (cm1.is_real && cm2.is_real));//两个都为实数矩阵则返回实数矩阵
	for (int i = 0; i < cm3.lr; i++)
		for (int j = 0; j < cm3.lc; j++)
			for (int k = 0; k < cm1.lc; k++)
			{
				Complex temp_c;
				temp_c = cm1.c[i][k] * cm2.c[k][j];
				cm3.c[i][j] = cm3.c[i][j] + temp_c;
			}
	return cm3;
}
//矩阵和二维数组乘号*
ComplexMatrix operator * (ComplexMatrix& cm1, double**& d1)
{
	ComplexMatrix cm2(d1, cm1.lr, cm1.lc);
	return cm1 * cm2;
}
ComplexMatrix operator * (double**& d1, ComplexMatrix& cm1)
{
	ComplexMatrix cm2(d1, cm1.lr, cm1.lc);
	return cm1 * cm2;
}

//矩阵数乘*
ComplexMatrix operator * (double& d1, ComplexMatrix& cm1)
{
	ComplexMatrix cm2(cm1.lr, cm1.lc, cm1.is_real, cm1.c);
	for (int i = 0; i < cm1.lr; i++)
		for (int j = 0; j < cm1.lc; j++)
			cm2.c[i][j] = d1 * cm2.c[i][j];
	return cm2;
}
ComplexMatrix operator * (ComplexMatrix& cm1, double& d1)
{
	ComplexMatrix cm2(cm1.lr, cm1.lc, cm1.is_real, cm1.c);
	for (int i = 0; i < cm1.lr; i++)
		for (int j = 0; j < cm1.lc; j++)
			cm2.c[i][j] = cm2.c[i][j] * d1;
	return cm2;
}
ComplexMatrix operator * (Complex& d1, ComplexMatrix& cm1)
{
	ComplexMatrix cm2(cm1.lr, cm1.lc, cm1.is_real, cm1.c);
	for (int i = 0; i < cm1.lr; i++)
		for (int j = 0; j < cm1.lc; j++)
			cm2.c[i][j] = d1 * cm2.c[i][j];
	return cm2;
}
ComplexMatrix operator * (ComplexMatrix& cm1, Complex& d1)
{
	ComplexMatrix cm2(cm1.lr, cm1.lc, cm1.is_real, cm1.c);
	for (int i = 0; i < cm1.lr; i++)
		for (int j = 0; j < cm1.lc; j++)
			cm2.c[i][j] = cm2.c[i][j] * d1;
	return cm2;
}

//矩阵输出<<
ostream & operator << (ostream& out, ComplexMatrix& cm1)
{
	if (cm1.is_real == false)		//复数矩阵输出
		for (int i = 0; i < cm1.lr; i++)
		{
			for (int j = 0; j < cm1.lc; j++)
				out << cm1.c[i][j] << " ";
			out << "\n";
		}
	else if (cm1.is_real == true)	//实数矩阵输出
		for (int i = 0; i < cm1.lr; i++)
		{
			for (int j = 0; j < cm1.lc; j++)
				out << setw(6) << cm1.c[i][j].re << " ";
			out << "\n";
		}
	return out;
}
//矩阵输入>>
istream & operator >> (istream& in, ComplexMatrix& cm1)
{
	if (cm1.is_real == false)		//复数矩阵输入
		//先输入一行所有的"实部 虚部",再输入下一行
		for (int i = 0; i < cm1.lr; i++)
			for (int j = 0; j < cm1.lc; j++)
				in >> cm1.c[i][j];
	else if (cm1.is_real == true)	//实数矩阵输入
		for (int i = 0; i < cm1.lr; i++)
			for (int j = 0; j < cm1.lc; j++)
				in >> cm1.c[i][j].re;
	return in;
}

//矩阵共轭转置!
ComplexMatrix ComplexMatrix::operator !()
{
	ComplexMatrix cm1(this->lc, this->lr, this->is_real);
	for (int i = 0; i < cm1.lr; i++)
		for (int j = 0; j < cm1.lc; j++)
			cm1.c[i][j] = !((this->c)[j][i]);
	return cm1;
}

//重载赋值运算符
ComplexMatrix& ComplexMatrix::operator = (const ComplexMatrix A)
{
	if (this->c != A.c)
	{
		this->clear();		//先清除当前等号左边对释放内存非常重要
		this->lr = A.lr;
		this->lc = A.lc;
		this->is_real = A.is_real;
		//创建矩阵内存空间
		this->c = new Complex * [A.lr];
		for (int i = 0; i < A.lr; i++)
			c[i] = new Complex[A.lc];
		//为矩阵赋值
		if (A.c == NULL)	//输入矩阵为空时
			this->c = NULL;
		else				//输入矩阵时
		{
			for (int i = 0; i < A.lr; i++)
				for (int j = 0; j < A.lc; j++)
					this->c[i][j] = A.c[i][j];
		}
	}
	return *this;
}

//换行
void ComplexMatrix::exchange_row(int i1, int i2)
{
	Complex temp;
	for (int j = 0; j < lc; j++)
	{
		temp = c[i1][j];
		c[i1][j] = c[i2][j];
		c[i2][j] = temp;
	}
}
//换列
void ComplexMatrix::exchange_column(int j1, int j2)
{	//换列
	Complex temp;
	for (int i = 0; i < lr; i++)
	{
		temp = c[i][j1];
		c[i][j1] = c[i][j2];
		c[i][j2] = temp;
	}
}
//换列的一个范围内的行
void ComplexMatrix::exchange_some_rows_of_column(int j1, int j2, int i1, int i2)
{
	Complex temp;
	for (int i = i1; i <= i2; i++)
	{
		temp = c[i][j1];
		c[i][j1] = c[i][j2];
		c[i][j2] = temp;
	}
}
//得到行（参数为存储下标）
ComplexMatrix ComplexMatrix::get_row(int i)
{
	ComplexMatrix target_row(1, lc, is_real);
	for (int k = 0; k < lc; k++) target_row.c[0][k] = c[i][k];
	return target_row;
}
//得到列（参数为存储下标）
ComplexMatrix ComplexMatrix::get_column(int j)
{
	ComplexMatrix target_column(lr, 1, is_real);
	for (int k = 0; k < lr; k++) target_column.c[k][0] = c[k][j];
	return target_column;
}
//得到连续的许多行（参数为存储下标范围）
ComplexMatrix ComplexMatrix::get_rows(int i1, int i2)
{
	//if (i1 > i2)  NULL
	ComplexMatrix target_rows(i2 - i1 + 1, lc, is_real);
	for (int i = i1; i <= i2; i++)
		for (int j = 0; j < lc; j++)
			target_rows.c[i - i1][j] = c[i][j];
	return target_rows;
}
//得到子矩阵（参数为存储下标范围）
ComplexMatrix ComplexMatrix::get_sub_matrix(int i1, int i2, int j1, int j2)
{
	//if (i1 > i2 || j1 > j2)  NULL
	ComplexMatrix target_sub_matrix(i2 - i1 + 1, j2 - j1 + 1, is_real);
	for (int i = i1; i <= i2; i++)
		for (int j = j1; j <= j2; j++)
			target_sub_matrix.c[i - i1][j - j1] = c[i][j];
	return target_sub_matrix;
}

//列合并：行数相同，拼成更多列的一个矩阵
ComplexMatrix ComplexMatrix::combine_columns(ComplexMatrix& A, ComplexMatrix& B)
{
	ComplexMatrix result(A.lr, A.lc + B.lc, (A.is_real && B.is_real));
	for (int i = 0; i < A.lr; i++)
	{
		for (int j = 0; j < A.lc; j++)result.c[i][j] = A.c[i][j];
		for (int j = A.lc; j < result.lc; j++)result.c[i][j] = B.c[i][j - A.lc];
	}
	return result;
}
//行合并：列数相同，拼成更多行的一个矩阵
ComplexMatrix ComplexMatrix::combine_rows(ComplexMatrix& A, ComplexMatrix& B)
{
	ComplexMatrix result(A.lr + B.lr, A.lc, (A.is_real && B.is_real));
	for (int j = 0; j < A.lc; j++)
	{
		for (int i = 0; i < A.lr; i++)result.c[i][j] = A.c[i][j];
		for (int i = A.lr; i < result.lr; i++)result.c[i][j] = B.c[i - A.lr][j];
	}
	return result;
}
//去除指定的一列（参数为存储下标）
ComplexMatrix ComplexMatrix::column_delete(int k)
{
	ComplexMatrix result(lr, lc - 1, is_real);
	for (int i = 0; i < lr; i++)
	{
		for (int j = 0; j < k; j++)result.c[i][j] = c[i][j];
		for (int j = k + 1; j < lc; j++)result.c[i][j - 1] = c[i][j];
	}
	return result;
}

//前向带入消元得到解x（实数）
void ComplexMatrix::forward_substitution(ComplexMatrix& A_b, ComplexMatrix& x)
{
	for (int i = 0; i < A_b.lr; i++)
	{
		for (int j = 0; j < i; j++)
			A_b.c[i][A_b.lc - 1].re -= (A_b.c[i][j].re * x.c[j][0].re);
		x.c[i][0].re = A_b.c[i][A_b.lc - 1].re / A_b.c[i][i].re;
	}
}
//前向带入消元得到解x（复数）
void ComplexMatrix::forward_substitution__Complex(ComplexMatrix& A_b, ComplexMatrix& x)
{
	Complex temp;
	for (int i = 0; i < A_b.lr; i++)
	{
		for (int j = 0; j < i; j++)
		{
			temp = A_b.c[i][j] * x.c[j][0];
			A_b.c[i][A_b.lc - 1] = A_b.c[i][A_b.lc - 1] - temp;
		}
		x.c[i][0] = A_b.c[i][A_b.lc - 1] / A_b.c[i][i];
	}
}
//后向带入消元得到解x（实数）
void ComplexMatrix::backward_substitution(ComplexMatrix& A_b, ComplexMatrix& x)
{
	for (int i = A_b.lr - 1; i >= 0; i--)
	{
		for (int j = A_b.lr - 1; j > i; j--)
			A_b.c[i][A_b.lc - 1].re -= (A_b.c[i][j].re * x.c[j][0].re);
		x.c[i][0].re = A_b.c[i][A_b.lc - 1].re / A_b.c[i][i].re;
	}
}
//后向带入消元得到解x（复数）
void ComplexMatrix::backward_substitution__Complex(ComplexMatrix& A_b, ComplexMatrix& x)
{
	Complex temp;
	for (int i = A_b.lr - 1; i >= 0; i--)
	{
		for (int j = A_b.lr - 1; j > i; j--)
		{
			temp = A_b.c[i][j] * x.c[j][0];
			A_b.c[i][A_b.lc - 1] = A_b.c[i][A_b.lc - 1] - temp;
		}
		x.c[i][0] = A_b.c[i][A_b.lc - 1] / A_b.c[i][i];
	}
}

//方阵求逆（复数）
ComplexMatrix ComplexMatrix::square_inverse()
{
	//解线程方程组Ax=b的方式：b分别取单位阵的各个列向量，得x即逆矩阵的对应列向量，拼成逆矩阵即可
	//if (lr != lc)return NULL;
	ComplexMatrix A(lr, lc, is_real, c);		//当前矩阵
	ComplexMatrix b(lr, 1, is_real);			//储存每次的b向量
	ComplexMatrix Ab(lr, lc + 1, is_real);		//储存每次的增广矩阵
	ComplexMatrix x(lr, 1, is_real);			//储存每次的解向量
	ComplexMatrix result(lr, lc, is_real);		//结果矩阵
	int ii = 0;									//迭代优化次数计数
	ComplexMatrix r;							//残差向量
	ComplexMatrix z(lr, 1, is_real);			//解修正向量
	for (int i = 0; i < lc; i++)	//外层循环，次数与列数相同
	{
		//b更新为第i行元素为1的列向量
		for (int k = 0; k < b.lr; k++)
		{
			b.c[k][0].re = 0.0;
			b.c[k][0].im = 0.0;
		}
		b.c[i][0].re = 1.0;
		//A和b拼成增广矩阵，并由部分选主元高斯解x
		Ab = combine_columns(A, b);
		Gaussian_elimination_partial_pivoting__Complex(Ab, x);
		ii = 0;
		//迭代优化10次
		do {
			r = A * x; r = b - r;
			Ab = Ab.column_delete(Ab.lc - 1);
			Ab = Ab.combine_columns(Ab, r);
			backward_substitution__Complex(Ab, z);
			x = x + z;
			ii++;
		} while (ii < 10);
		//解出的x是结果矩阵的第i列
		for (int k = 0; k < x.lr; k++) result.c[k][i] = x.c[k][0];
	}
	//释放内存
	//A.clear(); b.clear(); Ab.clear(); x.clear(); r.clear(); z.clear();
	return result;
}
//求伪逆（复数）
ComplexMatrix ComplexMatrix::Moore_Penrose_pseudo_inverse()
{
	ComplexMatrix A(lr, lc, is_real, c);
	ComplexMatrix A_H = !A;
	ComplexMatrix temp = A_H * A;
	temp = temp.square_inverse();
	temp = temp * A_H;
	return temp;
}

//向量2-范数（欧几里得范数）（复数）
double ComplexMatrix::vector_2_norm()
{
	double result = 0.0;
	if (lr == 1)		//行向量
	{
		for (int j = 0; j < lc; j++)
			result = result + c[0][j].re * c[0][j].re + c[0][j].im * c[0][j].im;
		result = sqrt(result);
		return result;
	}
	else if (lc == 1)	//列向量
	{
		for (int i = 0; i < lr; i++)
			result = result + c[i][0].re * c[i][0].re + c[i][0].im * c[i][0].im;
		result = sqrt(result);
		return result;
	}
	else return 0.0;	//都不是（控制输入以使该情况不要出现）
}

//部分选主元高斯消元（实数）
void ComplexMatrix::Gaussian_elimination_partial_pivoting(ComplexMatrix& A_b, ComplexMatrix& x)
{
	//化系数矩阵为上三角阵
	for (int i = 0; i < A_b.lr; i++)
	{
		//选出最大的主元
		double temp_max = A_b.c[i][i].re;
		int temp_max_row = i;
		for (int j = i + 1; j < A_b.lr; j++)
			if (abs(A_b.c[j][i].re) > abs(temp_max))
			{
				temp_max = A_b.c[j][i].re;
				temp_max_row = j;
			}
		A_b.exchange_row(i, temp_max_row);

		//用主元消元
		for (int k = i + 1; k < A_b.lr; k++)
		{
			double temp = -A_b.c[k][i].re / A_b.c[i][i].re;
			for (int j = i; j < A_b.lc; j++)
				A_b.c[k][j].re = A_b.c[k][j].re + temp * A_b.c[i][j].re;
		}
	}
	//后向带入消元得到解x
	backward_substitution(A_b, x);
}
//部分选主元高斯消元（复数）
void ComplexMatrix::Gaussian_elimination_partial_pivoting__Complex(ComplexMatrix& A_b, ComplexMatrix& x)
{
	//化系数矩阵为上三角阵
	for (int i = 0; i < A_b.lr; i++)
	{
		//选出最大的主元
		Complex temp_max = A_b.c[i][i];
		int temp_max_row = i;
		for (int j = i + 1; j < A_b.lr; j++)
			if (A_b.c[j][i].modulus() > temp_max.modulus())
			{
				temp_max = A_b.c[j][i];
				temp_max_row = j;
			}
		A_b.exchange_row(i, temp_max_row);

		//用主元消元
		for (int k = i + 1; k < A_b.lr; k++)
		{
			Complex temp = A_b.c[k][i] / A_b.c[i][i];
			Complex temptemp;
			for (int j = i; j < A_b.lc; j++)
			{
				temptemp = temp * A_b.c[i][j];
				A_b.c[k][j] = A_b.c[k][j] - temptemp;
			}
		}
	}
	//后向带入消元得到解x
	backward_substitution__Complex(A_b, x);
}


//解QR分解后的增广矩阵表示的超定方程（实数）
void ComplexMatrix::solution_of_augmentedMatrix_after_QR(ComplexMatrix& Ab, ComplexMatrix& x)
{
	//后向带入解系统的最小二乘解
	ComplexMatrix Rb(x.lr, x.lr + 1, true);
	for (int i = 0; i < Rb.lr; i++)
		for (int j = 0; j < Rb.lc; j++)
			Rb.c[i][j].re = Ab.c[i][j].re;
	backward_substitution(Rb, x);
}
//增广矩阵的Householder变换（实数）
void ComplexMatrix::Householder_QR_augmented(ComplexMatrix& A)
{
	ComplexMatrix v_k(A.lr, 1, true);		//储存豪斯霍尔德向量
	ComplexMatrix v_k_T(1, A.lr, true);		//储存豪斯霍尔德向量的转置
	ComplexMatrix v_k_j(A.lr, 1, true);		//储存参与运算的豪斯霍尔德向量
	ComplexMatrix vkT_vk(1, 1, true);		//储存豪斯霍尔德向量的自身内积
	double beta_k = 0.0;					//储存豪斯霍尔德向量的自身内积
	ComplexMatrix A_j(A.lr, 1, true);		//储存A矩阵提取出来的列
	ComplexMatrix gamma_j_CM(1, 1, true);	//储存对剩余子矩阵做变换时的系数
	double gamma_j = 0.0;					//储存对剩余子矩阵做变换时的系数
	for (int k = 0; k < A.lc - 1; k++)	//for增广矩阵
	{
		//计算当前列的豪斯霍尔德向量
		double square_sum = 0.0;
		for (int i = k; i < A.lr; i++)square_sum += pow(A.c[i][k].re, 2);
		double alpha_k = ((A.c[k][k].re >= 0) ? (-1.0) : 1.0) * sqrt(square_sum);
		for (int i = 0; i < k; i++)v_k.c[i][0].re = 0.0;
		for (int i = k; i < A.lr; i++)v_k.c[i][0].re = A.c[i][k].re;
		v_k.c[k][0].re = v_k.c[k][0].re - alpha_k;
		v_k_T = !v_k;
		vkT_vk = v_k_T * v_k;
		beta_k = vkT_vk.c[0][0].re;
		//如果当前列已经为零，跳过
		if (abs(beta_k) < 1e-323) continue;
		//对剩余的子矩阵做变换
		for (int j = k; j < A.lc; j++)
		{
			for (int i = 0; i < A.lr; i++)A_j.c[i][0].re = A.c[i][j].re;
			gamma_j_CM = v_k_T * A_j;
			gamma_j = gamma_j_CM.c[0][0].re;
			gamma_j = 2.0 * gamma_j / beta_k;
			for (int i = 0; i < A.lr; i++)v_k_j.c[i][0].re = v_k.c[i][0].re;
			v_k_j = gamma_j * v_k;
			A_j = A_j - v_k_j;
			for (int i = 0; i < A.lr; i++)A.c[i][j].re = A_j.c[i][0].re;
		}
	}
}
//增广矩阵的吉文斯旋转变换（实数）
void ComplexMatrix::Givens_QR(ComplexMatrix& A)
{	//对系数矩阵主对角线下元素循环
	for (int i = 0; i < A.lc - 1; i++)
		for (int j = i + 1; j < A.lr; j++)
			if (abs(A.c[j][i].re) >= 1e-323)	//非零元素才用消除
			{
				//制造吉文斯旋转矩阵
				ComplexMatrix Givens_rotation(A.lr, A.lr, true);
				for (int k = 0; k < Givens_rotation.lr; k++)Givens_rotation.c[k][k].re = 1.0;
				double a1 = A.c[i][i].re;	//利用所在列主对角元完成旋转，防止前一列又出现非零元
				double a2 = A.c[j][i].re;
				double c = a1 / sqrt(pow(a1, 2) + pow(a2, 2));
				double s = a2 / sqrt(pow(a1, 2) + pow(a2, 2));
				Givens_rotation.c[i][i].re = c;
				Givens_rotation.c[i][j].re = s;
				Givens_rotation.c[j][i].re = -s;
				Givens_rotation.c[j][j].re = c;
				//吉文斯旋转
				A = Givens_rotation * A;
			}
}
//古典格拉姆-施密特正交化的QR分解（实数）
ComplexMatrix ComplexMatrix::Gram_Schmidt_QR_classical(ComplexMatrix& Q)
{
	//Q将被分解，原始数据存入A
	ComplexMatrix A(Q.lr, Q.lc, true);
	for (int i = 0; i < Q.lr; i++)for (int j = 0; j < Q.lc; j++)A.c[i][j].re = Q.c[i][j].re;
	//QR分解的R
	ComplexMatrix R(Q.lc, Q.lc, true);

	ComplexMatrix qk(A.lr, 1, true);	//对列操作
	ComplexMatrix qkT(1, A.lr, true);	//对列操作
	ComplexMatrix qj(A.lr, 1, true);	//对列操作
	ComplexMatrix qjT(1, A.lr, true);	//对列操作
	ComplexMatrix ak(A.lr, 1, true);	//提取列

	ComplexMatrix rjk(1, 1, true);		//内积临时变量
	double r_jk = 0.0;					//内积临时变量
	ComplexMatrix rkk(1, 1, true);		//二范数临时变量
	double r_kk = 0.0;					//二范数临时变量
	for (int k = 0; k < A.lc; k++)	//对列循环
	{
		for (int i = 0; i < A.lr; i++)
		{
			qk.c[i][0].re = A.c[i][k].re;
			ak.c[i][0].re = A.c[i][k].re;
		}
		//从当前列中减去它在前面列中的分量
		for (int j = 0; j < k; j++)
		{
			for (int i = 0; i < Q.lr; i++)
				qj.c[i][0].re = Q.c[i][j].re;
			qjT = !qj;
			rjk = qjT * ak;
			r_jk = rjk.c[0][0].re;
			R.c[j][k].re = r_jk;
			qj = r_jk * qj;
			qk = qk - qj;
		}
		qkT = !qk;
		rkk = qkT * qk;
		r_kk = sqrt(rkk.c[0][0].re);
		R.c[k][k].re = r_kk;
		//如果线性相关，则中断
		if (abs(r_kk) < 1e-323)break;
		//将当前列标准化
		r_kk = 1 / r_kk;
		qk = r_kk * qk;
		for (int i = 0; i < Q.lr; i++)Q.c[i][k].re = qk.c[i][0].re;
	}
	return R;
}
//改进格拉姆-施密特正交化的QR分解（实数）
ComplexMatrix ComplexMatrix::Gram_Schmidt_QR_modified(ComplexMatrix& Q)
{
	//Q将被分解，原始数据存入A
	ComplexMatrix A(Q.lr, Q.lc, true);
	for (int i = 0; i < Q.lr; i++)for (int j = 0; j < Q.lc; j++)A.c[i][j].re = Q.c[i][j].re;
	//QR分解的R
	ComplexMatrix R(Q.lc, Q.lc, true);

	ComplexMatrix ak(A.lr, 1, true);		//提取列
	ComplexMatrix akT(1, A.lr, true);		//列转置
	ComplexMatrix aj(A.lr, 1, true);		//对列操作
	ComplexMatrix qk(A.lr, 1, true);		//对列操作
	ComplexMatrix r_kj_qk(A.lr, 1, true);	//对列操作
	ComplexMatrix qkT(1, A.lr, true);		//对列操作

	ComplexMatrix rkj(1, 1, true);		//内积临时变量
	double r_kj = 0.0;					//内积临时变量
	ComplexMatrix rkk(1, 1, true);		//二范数临时变量
	double r_kk = 0.0;					//二范数临时变量

	for (int k = 0; k < A.lc; k++)	//对列循环
	{
		for (int i = 0; i < A.lr; i++)
			ak.c[i][0].re = A.c[i][k].re;
		akT = !ak;
		rkk = akT * ak;
		r_kk = sqrt(rkk.c[0][0].re);
		R.c[k][k] = r_kk;
		//如果线性相关，则中断
		if (abs(r_kk) < 1e-323)break;
		//将当前列标准化
		r_kk = 1 / r_kk;
		qk = r_kk * ak;
		qkT = !qk;
		for (int i = 0; i < Q.lr; i++)Q.c[i][k].re = qk.c[i][0].re;
		//减去后续列在当前列上的分量
		for (int j = k + 1; j < A.lc; j++)
		{
			for (int i = 0; i < A.lr; i++) aj.c[i][0].re = A.c[i][j].re;
			rkj = qkT * aj;
			R.c[k][j] = r_kj = rkj.c[0][0].re;
			r_kj_qk = r_kj * qk;
			aj = aj - r_kj_qk;
			for (int i = 0; i < A.lr; i++) A.c[i][j].re = aj.c[i][0].re;
		}
	}
	return R;
}
//改进格拉姆-施密特正交化的QR分解（复数）
ComplexMatrix ComplexMatrix::Gram_Schmidt_QR_modified__Complex(ComplexMatrix& Q)
{
	//QR分解的R
	ComplexMatrix R(Q.lc, Q.lc, Q.is_real);

	ComplexMatrix ak(Q.lr, 1);		//提取列
	ComplexMatrix aj(Q.lr, 1);		//对列操作
	ComplexMatrix qk(Q.lr, 1);		//对列操作
	ComplexMatrix r_kj_qk(Q.lr, 1);	//对列操作
	ComplexMatrix qkH(1, Q.lr);		//对列操作

	ComplexMatrix r_kj(1, 1);		//内积临时变量
	double r_kk = 0.0;				//二范数临时变量

	for (int k = 0; k < Q.lc; k++)	//对列循环
	{
		ak = Q.get_column(k);
		r_kk = ak.vector_2_norm();
		R.c[k][k].re = r_kk;
		//如果线性相关，则中断
		if (abs(r_kk) < 1e-323)break;
		//将当前列标准化
		r_kk = 1 / r_kk;
		qk = r_kk * ak;
		qkH = !qk;
		for (int i = 0; i < Q.lr; i++)Q.c[i][k] = qk.c[i][0];
		//减去后续列在当前列上的分量
		for (int j = k + 1; j < Q.lc; j++)
		{
			aj = Q.get_column(j);
			r_kj = qkH * aj;
			R.c[k][j] = r_kj.c[0][0];
			r_kj_qk = r_kj.c[0][0] * qk;
			aj = aj - r_kj_qk;
			for (int i = 0; i < Q.lr; i++) Q.c[i][j] = aj.c[i][0];
		}
	}
	//释放内存
	//ak.clear(); aj.clear(); qk.clear(); r_kj_qk.clear(); qkH.clear(); r_kj.clear();
	return R;
}
//有选择的改进格拉姆-施密特正交化的QR分解（复数）
//算法来源：大二上《数值计算方法》涉及的论文 Efficient Algorithm for Detecting Layered Space-Time Codes
ComplexMatrix ComplexMatrix::sorted_Gram_Schmidt_QR_modified__Complex(ComplexMatrix& Q, int* S)
{
	//QR分解的R
	ComplexMatrix R(Q.lc, Q.lc, Q.is_real);

	ComplexMatrix ai(Q.lr, 1);		//提取列
	ComplexMatrix aj(Q.lr, 1);		//对列操作
	ComplexMatrix qi(Q.lr, 1);		//对列操作
	ComplexMatrix r_ij_qi(Q.lr, 1);	//对列操作
	ComplexMatrix qiH(1, Q.lr);		//对列操作

	ComplexMatrix q;				//存储Q的列
	double q_2_norm;				//存储q的2-范数
	double q_2_norm_min;			//存储Q的2-范数最小的一列的2-范数
	int q_2_norm_min_column;		//存储Q的2-范数最小的一列的存储下标
	int temp;						//交换中介

	ComplexMatrix r_ij(1, 1);		//内积临时变量
	double r_ii = 0.0;				//二范数临时变量

	for (int i = 0; i < Q.lc; i++)	//对列循环
	{
		//找出Q剩下列中2-范数最小的一列
		q_2_norm_min_column = i;
		q = Q.get_column(i);
		q_2_norm_min = q.vector_2_norm();
		for (int k = i + 1; k < Q.lc; k++)
		{
			q = Q.get_column(k);
			q_2_norm = q.vector_2_norm();
			if (q_2_norm < q_2_norm_min)
			{
				q_2_norm_min_column = k;
				q_2_norm_min = q_2_norm;
			}
		}
		//交换Q,R,S中刚求出的2-范数最小列对应的列和当前的第i列
		Q.exchange_column(i, q_2_norm_min_column);
		R.exchange_column(i, q_2_norm_min_column);
		temp = S[i];
		S[i] = S[q_2_norm_min_column];
		S[q_2_norm_min_column] = temp;
		//开始正交化
		ai = Q.get_column(i);
		r_ii = ai.vector_2_norm();
		R.c[i][i].re = r_ii;
		//如果线性相关，则中断
		if (abs(r_ii) < 1e-323)break;
		//将当前列标准化
		r_ii = 1 / r_ii;
		qi = r_ii * ai;
		qiH = !qi;
		for (int k = 0; k < Q.lr; k++)Q.c[k][i] = qi.c[k][0];
		//减去后续列在当前列上的分量
		for (int j = i + 1; j < Q.lc; j++)
		{
			aj = Q.get_column(j);
			r_ij = qiH * aj;
			R.c[i][j] = r_ij.c[0][0];
			r_ij_qi = r_ij.c[0][0] * qi;
			aj = aj - r_ij_qi;
			for (int k = 0; k < Q.lr; k++) Q.c[k][j] = aj.c[k][0];
		}
	}
	//释放内存
	//ai.clear(); aj.clear(); qi.clear(); r_ij_qi.clear(); qiH.clear(); q.clear(); r_ij.clear();
	return R;
}

//牛顿插值获得多项式（实数）
ComplexMatrix ComplexMatrix::Newton_interpolation_get_polynomial(ComplexMatrix& t, ComplexMatrix& y)
{
	ComplexMatrix pi_(t.lr, 1, true);		//插值多项式系数
	ComplexMatrix Ab(t.lr, t.lr + 1, true);	//增广矩阵
	//给增广矩阵赋值
	for (int i = 0; i < Ab.lr; i++)
	{
		Ab.c[i][Ab.lc - 1].re = y.c[i][0].re;	//等式右边的y
		for (int j = Ab.lc - 2; j > i; j--)Ab.c[i][j].re = 0.0;	//严格上三角的零
		for (int j = 0; j <= i; j++)	//下三角的系数
		{
			Ab.c[i][j].re = 1.0;
			for (int k = 0; k < j; k++)
				Ab.c[i][j].re *= (t.c[i][0].re - t.c[k][0].re);
		}
	}
	//解插值多项式系数
	forward_substitution(Ab, pi_);
	return pi_;
}
//牛顿插值多项式秦九韶算法算值（实数）
double ComplexMatrix::Newton_interpolation_get_value(double t, ComplexMatrix& t_x, ComplexMatrix& pi_)
{
	double p = pi_.c[pi_.lr-1][0].re;
	for (int i = pi_.lr - 2; i >= 0; i--)
	{
		p *= (t - t_x.c[i][0].re);
		p += pi_.c[i][0].re;
	}
	return p;
}
//牛顿插值多一个点得到新的多项式（实数）
ComplexMatrix ComplexMatrix::Newton_interpolation_add_one_point(double& x_new, double& y_new, ComplexMatrix& t, ComplexMatrix& pi_)
{
	//t是老的自变量表，pi_是老的插值多项式系数
	ComplexMatrix pi_new(pi_.lr + 1, 1, true);
	for (int i = 0; i < pi_.lr; i++)pi_new.c[i][0].re = pi_.c[i][0].re;
	ComplexMatrix new_equation(1, pi_.lr + 1, true);	//新加一个点导致求插值系数新加一个方程，之前的系数不变
	//新方程的老自变量个数个系数
	for (int j = 0; j < pi_.lr; j++)
	{
		new_equation.c[0][j].re = 1.0;
		for (int k = 0; k < j; k++)
			new_equation.c[0][j].re *= (x_new - t.c[k][0].re);
	}
	//新方程的最后一个系数
	new_equation.c[0][pi_.lr].re = 1.0;
	for (int k = 0; k < pi_.lr; k++)
		new_equation.c[0][pi_.lr].re *= (x_new - t.c[k][0].re);
	//代入解新的系数
	double pn = y_new;
	for (int i = 0; i < pi_.lr; i++)
		pn -= (new_equation.c[0][i].re * pi_.c[i][0].re);
	pn /= new_equation.c[0][pi_.lr].re;
	//返回结果
	pi_new.c[pi_new.lr - 1][0].re = pn;
	return pi_new;
}
//牛顿插值获得多项式（递归方法）（输入参数为最大下标）（实数）
ComplexMatrix ComplexMatrix::Newton_interpolation_get_polynomial_recursive(ComplexMatrix& t, ComplexMatrix& y, int count)
{
	if (count == 0)
	{
		ComplexMatrix pi_(1, 1, true);
		pi_.c[0][0].re = y.c[0][0].re;
		return pi_;
	}
	else
	{
		ComplexMatrix pi_0 = Newton_interpolation_get_polynomial_recursive(t, y, count - 1);
		ComplexMatrix pi_ = Newton_interpolation_add_one_point(t.c[count][0].re, y.c[count][0].re, t, pi_0);
		return pi_;
	}
}

//优化：最速下降法（实数）
double ComplexMatrix::Steepest_Descent(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x))
{
	ComplexMatrix xk(x0.lr, x0.lc, x0.is_real, x0.c);	//储存迭代解
	ComplexMatrix sk(xk.lr, xk.lc, xk.is_real);			//储存负梯度
	ComplexMatrix sk0;									//储存修正用的负梯度
	double alpha0;										//线性搜索中上一次步长
	double alpha;										//线性搜索步长
	double fx0;											//线性搜索步长中上一次的函数值
	double fx;											//线性搜索步长中本次的函数值
	//dataOut.open("Steepest_Descent.txt");							//输出到文件，不用时注释掉
	//dataOut << xk.c[0][0].re << " " << xk.c[1][0].re << endl;		//输出x到文件，不用时注释掉
	//求初始负梯度
	for (int i = 0; i < xk.lr; i++)
		sk.c[i][0].re = -((grad[i])(xk));
	do {
		//线性搜索求步长（alpha每次增加0.001）					//考虑第一次直接用sk，之后每次/2的搜索方法
		ComplexMatrix sk_temp(sk.lr, sk.lc, sk.is_real, sk.c);	//储存线性搜索步长时临时s
		ComplexMatrix xk_temp(xk.lr, xk.lc, xk.is_real, xk.c);	//储存线性搜索步长时临时x
		alpha = 0.0;	//初始化步长
		do {
			alpha0 = alpha;
			fx0 = f(xk_temp);
			alpha += 0.001;
			sk_temp = alpha * sk;
			xk_temp = xk_temp + sk_temp;
			fx = f(xk_temp);
		} while (fx < fx0);	//脱出时fx >= fx0，fx0对应的alpha0为所求找到最小值的步长
		//修正解
		sk = alpha0 * sk;
		sk0 = sk;
		xk = xk + sk;		//cout << xk << endl;		//输出迭代路径，不用时注释掉
		//求新负梯度
		for (int i = 0; i < xk.lr; i++)
			sk.c[i][0].re = -((grad[i])(xk));
	//	dataOut << xk.c[0][0].re << " " << xk.c[1][0].re << endl;		//输出x到文件，不用时注释掉
	} while (sk0.vector_2_norm() >= 1e-15);	//负梯度为零时脱出								//尝试增加try计尝试次数
	//dataOut.close();													//输出到文件，不用时注释掉
	cout << xk << endl;													//输出最终x，不用时注释掉
	return f(xk);
}
//优化：牛顿法（实数）
double ComplexMatrix::Newton_unconstrained_optimization
(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x), double(*Hessian[])(ComplexMatrix& x))
{
	ComplexMatrix xk(x0.lr, x0.lc, x0.is_real, x0.c);	//储存迭代解向量
	ComplexMatrix gradk(xk.lr, xk.lc, xk.is_real);		//储存负梯度向量
	ComplexMatrix Hk(xk.lr, xk.lr, xk.is_real);			//储存黑塞矩阵
	ComplexMatrix Hk_gradk;								//储存增广矩阵
	ComplexMatrix sk(xk.lr, xk.lc, xk.is_real);			//储存牛顿步长向量
	//dataOut.open("Newton_optimization.txt");						//输出到文件，不用时注释掉
	//dataOut << xk.c[0][0].re << " " << xk.c[1][0].re << endl;		//输出x到文件，不用时注释掉
	//求初始负梯度
	for (int i = 0; i < xk.lr; i++)
		gradk.c[i][0].re = -((grad[i])(xk));
	do {
		//求牛顿步长
			//求黑塞矩阵
		for (int i = 0; i < xk.lr; i++)
			for (int j = 0; j < xk.lr; j++)
				Hk.c[i][j].re = (Hessian[i * x0.lr + j])(xk);
			//产生增广矩阵，用部分选主元高斯消元解方程组
		Hk_gradk = Hk.combine_columns(Hk, gradk);
		Gaussian_elimination_partial_pivoting(Hk_gradk, sk);
		//修正解
		xk = xk + sk;						//cout << xk << endl;		//输出迭代路径，不用时注释掉
		//求新负梯度
		for (int i = 0; i < xk.lr; i++)
			gradk.c[i][0].re = -((grad[i])(xk));
	//	dataOut << xk.c[0][0].re << " " << xk.c[1][0].re << endl;		//输出x到文件，不用时注释掉
	} while (gradk.vector_2_norm() >= 1e-15);	//负梯度为零时脱出								//尝试增加try计尝试次数
	//dataOut.close();													//输出到文件，不用时注释掉
	cout << xk << endl;													//输出最终x，不用时注释掉
	return f(xk);
}
//优化：阻尼牛顿法（实数）
double ComplexMatrix::damped_Newton_unconstrained_optimization
(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x), double(*Hessian[])(ComplexMatrix& x))
{
	ComplexMatrix xk(x0.lr, x0.lc, x0.is_real, x0.c);	//储存迭代解向量
	ComplexMatrix gradk(xk.lr, xk.lc, xk.is_real);		//储存负梯度向量
	ComplexMatrix Hk(xk.lr, xk.lr, xk.is_real);			//储存黑塞矩阵
	ComplexMatrix Hk_gradk;								//储存增广矩阵
	ComplexMatrix sk(xk.lr, xk.lc, xk.is_real);			//储存牛顿步长向量
	double alpha0;										//线性搜索中上一次步长
	double alpha;										//线性搜索步长
	double fx0;											//线性搜索步长中上一次的函数值
	double fx;											//线性搜索步长中本次的函数值
	//dataOut.open("Damped_Newton.txt");								//输出到文件，不用时注释掉
	//dataOut << xk.c[0][0].re << " " << xk.c[1][0].re << endl;			//输出x到文件，不用时注释掉
	//求初始负梯度
	for (int i = 0; i < xk.lr; i++)
		gradk.c[i][0].re = -((grad[i])(xk));
	do {
		//求牛顿步长
			//求黑塞矩阵
		for (int i = 0; i < xk.lr; i++)
			for (int j = 0; j < xk.lr; j++)
				Hk.c[i][j].re = (Hessian[i * x0.lr + j])(xk);
			//产生增广矩阵，用部分选主元高斯消元解方程组
		Hk_gradk = Hk.combine_columns(Hk, gradk);
		Gaussian_elimination_partial_pivoting(Hk_gradk, sk);
		//充分接近负梯度为零的点前线性搜索步长，充分接近后步长为1					//考虑第一次直接用sk，之后每次/2的搜索方法
		if (gradk.vector_2_norm() < 1e-3)	//负梯度很小，充分接近
			alpha0 = 1;
		else								//线性搜索求步长（alpha每次增加0.001）
		{
			ComplexMatrix sk_temp(sk.lr, sk.lc, sk.is_real, sk.c);	//储存线性搜索步长时临时s
			ComplexMatrix xk_temp(xk.lr, xk.lc, xk.is_real, xk.c);	//储存线性搜索步长时临时x
			alpha = 0.0;	//初始化步长
			do {
				alpha0 = alpha;
				fx0 = f(xk_temp);
				alpha += 0.001;
				sk_temp = alpha * sk;
				xk_temp = xk_temp + sk_temp;
				fx = f(xk_temp);
			} while (fx < fx0);	//脱出时fx >= fx0，fx0对应的alpha0为所求找到最小值的步长
		}
		//修正解
		sk = alpha0 * sk;
		xk = xk + sk;					//cout << xk << endl;		//输出迭代路径，不用时注释掉
		//求新负梯度
		for (int i = 0; i < xk.lr; i++)
			gradk.c[i][0].re = -((grad[i])(xk));
	//	dataOut << xk.c[0][0].re << " " << xk.c[1][0].re << endl;		//输出x到文件，不用时注释掉
	} while (gradk.vector_2_norm() >= 1e-15);	//负梯度为零时脱出								//尝试增加try计尝试次数
	//dataOut.close();													//输出到文件，不用时注释掉
	cout << xk << endl;													//输出最终x，不用时注释掉
	return f(xk);
}
//无约束优化的BFGS割线修正法（实数）
double ComplexMatrix::BFGS_optimization(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x))
{
	ComplexMatrix Bk = make_eyes(x0.lr);				//初始黑塞近似为单位阵
	ComplexMatrix xk(x0);								//初始化迭代解
	ComplexMatrix yk(x0);								//修正黑塞近似中的临时变量
	ComplexMatrix sk(x0.lr, x0.lc, x0.is_real);				//储存步长
	ComplexMatrix gradk(x0.lr, x0.lc, x0.is_real);			//储存负梯度
	ComplexMatrix gradk0(gradk);							//临时储存上一次的负梯度
	ComplexMatrix Bk_gradk = combine_columns(Bk, gradk);	//储存增广矩阵
	//初始负梯度
	for (int i = 0; i < xk.lr; i++)
		gradk.c[i][0].re = -((grad[i])(xk));
	do {

		//产生增广矩阵，用部分选主元高斯消元解方程组
		Bk_gradk = combine_columns(Bk, gradk);
		Gaussian_elimination_partial_pivoting(Bk_gradk, sk);
		//修正解
		xk = xk + sk;					cout << xk << endl;		//输出迭代路径，不用时注释掉
		//修正黑塞近似矩阵
		//求新的负梯度
		gradk0 = gradk;
		for (int i = 0; i < xk.lr; i++)
			gradk.c[i][0].re = -((grad[i])(xk));
		yk = gradk0 - gradk;

		ComplexMatrix ykT = (!yk);
		ComplexMatrix yk_ykT = yk * ykT;
		ComplexMatrix ykT_sk = ykT * sk;
		double ykT_sk_val = ykT_sk.c[0][0].re; ykT_sk_val = 1 / ykT_sk_val;
		yk_ykT = yk_ykT * ykT_sk_val;

		ComplexMatrix skT = (!sk);
		ComplexMatrix temp1 = Bk * sk;
		temp1 = temp1 * skT;
		temp1 = temp1 * Bk;
		ComplexMatrix temp2 = skT * Bk;
		temp2 = temp2 * sk;
		double temp2_val = temp2.c[0][0].re; temp2_val = 1 / temp2_val;
		temp1 = temp1 * temp2_val;

		Bk = Bk + yk_ykT;
		Bk = Bk - temp1;
	} while (gradk.vector_2_norm() >= 1e-15);			//负梯度为零时脱出								//尝试增加try计尝试次数
	cout << xk << endl;													//输出最终x，不用时注释掉
	return f(xk);
}
//无约束优化的共轭梯度法（实数）beta由Fletcher-Reeves公式给出
double ComplexMatrix::Conjugate_Gradient_F_R(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x))
{
	ComplexMatrix temp_0(x0.lr, x0.lc, x0.is_real);	//零向量
	ComplexMatrix temp1, temp2;						//计算中间变量
	double beta;									//迭代修正解的梯度修正因子
	double alpha0;										//线性搜索中上一次步长
	double alpha;										//线性搜索步长
	double fx0;											//线性搜索步长中上一次的函数值
	double fx;											//线性搜索步长中本次的函数值
	ComplexMatrix xk(x0);							//初始迭代解
	ComplexMatrix gk(x0.lr, x0.lc, x0.is_real);		//初始梯度
	for (int i = 0; i < xk.lr; i++) gk.c[i][0].re = ((grad[i])(xk));
	ComplexMatrix gk0(gk);							//上次迭代中的梯度
	ComplexMatrix gkT;								//梯度转置
	ComplexMatrix gk0T;								//上次迭代中的梯度转置
	ComplexMatrix sk(x0.lr, x0.lc, x0.is_real);		//初始化步长
	sk = temp_0 - gk;
	ComplexMatrix sk0(sk);							//上一次迭代中的步长
	do {
		//线性搜索求步长
		ComplexMatrix sk_temp(sk);	//储存线性搜索步长时临时s
		ComplexMatrix xk_temp(xk);	//储存线性搜索步长时临时x
		alpha = 0.0;	//初始化步长
		do {
			alpha0 = alpha;
			fx0 = f(xk_temp);
			alpha += 0.001;
			sk_temp = alpha * sk;
			xk_temp = xk_temp + sk_temp;
			fx = f(xk_temp);
		} while (fx < fx0);	//脱出时fx >= fx0，fx0对应的alpha0为所求找到最小值的步长
		//修正解
		sk0 = sk;
		sk = alpha0 * sk;
		xk = xk + sk;		//cout << xk << endl;		//输出迭代路径，不用时注释掉
		//计算新梯度
		gk0 = gk; gk0T = (!gk0);
		for (int i = 0; i < xk.lr; i++) gk.c[i][0].re = (grad[i])(xk);
		gkT = (!gk);
		//beta由Fletcher-Reeves公式给出
		temp1 = gkT * gk;
		temp2 = gk0T * gk0;
		beta = temp1.c[0][0].re / temp2.c[0][0].re;
		//修正梯度
		sk = beta * sk;
		sk = sk - gk;							//static int i = 0; i++; cout << i << " ";//输出迭代次数，不用时注释掉
	}while (gk.vector_2_norm() >= 1e-15);				//梯度为零时脱出								//尝试增加try计尝试次数
	cout << xk << endl;													//输出最终x，不用时注释掉
	return f(xk);
}
//无约束优化的共轭梯度法（实数）beta由Polak-Ribiere公式给出
double ComplexMatrix::Conjugate_Gradient_P_R(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x))
{
	ComplexMatrix temp_0(x0.lr, x0.lc, x0.is_real);	//零向量
	ComplexMatrix temp1, temp2;						//计算中间变量
	double beta;									//迭代修正解的梯度修正因子
	double alpha0;										//线性搜索中上一次步长
	double alpha;										//线性搜索步长
	double fx0;											//线性搜索步长中上一次的函数值
	double fx;											//线性搜索步长中本次的函数值
	ComplexMatrix xk(x0);							//初始迭代解
	ComplexMatrix gk(x0.lr, x0.lc, x0.is_real);		//初始梯度
	for (int i = 0; i < xk.lr; i++) gk.c[i][0].re = ((grad[i])(xk));
	ComplexMatrix gk0(gk);							//上次迭代中的梯度
	ComplexMatrix gk0T;								//上次迭代中的梯度转置
	ComplexMatrix sk(x0.lr, x0.lc, x0.is_real);		//初始化步长
	sk = temp_0 - gk;
	ComplexMatrix sk0(sk);							//上一次迭代中的步长
	do {
		//线性搜索求步长
		ComplexMatrix sk_temp(sk.lr, sk.lc, sk.is_real, sk.c);	//储存线性搜索步长时临时s
		ComplexMatrix xk_temp(xk.lr, xk.lc, xk.is_real, xk.c);	//储存线性搜索步长时临时x
		alpha = 0.0;	//初始化步长
		do {
			alpha0 = alpha;
			fx0 = f(xk_temp);
			alpha += 0.001;
			sk_temp = alpha * sk;
			xk_temp = xk_temp + sk_temp;
			fx = f(xk_temp);
		} while (fx < fx0);	//脱出时fx >= fx0，fx0对应的alpha0为所求找到最小值的步长
		//修正解
		sk0 = sk;
		sk = alpha0 * sk;
		xk = xk + sk;		//cout << xk << endl;		//输出迭代路径，不用时注释掉
		//计算新梯度
		gk0 = gk; gk0T = (!gk0);
		for (int i = 0; i < xk.lr; i++) gk.c[i][0].re = (grad[i])(xk);
		//beta由Polak-Ribiere公式给出
		temp1 = gk - gk0;
		temp1 = (!temp1);
		temp1 = temp1 * gk;
		temp2 = gk0T * gk0;
		beta = temp1.c[0][0].re / temp2.c[0][0].re;
		//修正梯度
		sk = beta * sk0;
		sk = sk - gk;			//static int i = 0; i++; cout << i << " ";//输出迭代次数，不用时注释掉
	} while (sk.vector_2_norm() >= 1e-15);				//梯度为零时脱出								//尝试增加try计尝试次数
	cout << xk << endl;													//输出最终x，不用时注释掉
	return f(xk);
}
//非线性最小二乘：高斯-牛顿法（实数）
ComplexMatrix ComplexMatrix::Gauss_Newton_nonlinear_least_squares
(ComplexMatrix& x0, const ComplexMatrix& t, const ComplexMatrix& y,
	double(*f)(double, ComplexMatrix&), double(*Jacobian[])(double, ComplexMatrix&))
{
	ComplexMatrix xk(x0);						//初始化迭代最小二乘解
	ComplexMatrix sk(x0.lr, 1, x0.is_real);		//修正步长
	ComplexMatrix r(t.lr, 1, t.is_real);		//负残差向量
	ComplexMatrix J(t.lr, x0.lr, t.is_real);	//雅可比矩阵
	ComplexMatrix Jr;							//储存增广矩阵
	for (int i = 0; i < t.lr; i++)				//计算初始负残差
		r.c[i][0].re = -y.c[i][0].re + f(t.c[i][0].re, xk);
	ComplexMatrix r0(r);						//前一个负残差向量
	ComplexMatrix r_r0 = r - r0;				//负残差向量的变化量
	do {
		for (int i = 0; i < t.lr; i++)		//计算雅可比矩阵
			for (int j = 0; j < xk.lr; j++)
				J.c[i][j].re = Jacobian[j](t.c[i][0].re, xk);
		Jr = combine_columns(J, r);			//计算新步长的最小二乘解
		Householder_QR_augmented(Jr);
		solution_of_augmentedMatrix_after_QR(Jr, sk);
		xk = xk + sk;			cout << xk << endl;//显示迭代路径，不用时注释掉
		r0 = r;
		for (int i = 0; i < t.lr; i++)		//计算新的负残差
			r.c[i][0].re = -y.c[i][0].re + f(t.c[i][0].re, xk);
		r_r0 = r - r0;
	} while (r_r0.vector_2_norm() >= 1e-4);	//残差在误差限内不再变化时脱出
	return xk;
}
//等式约束优化：拉格朗日乘子法，逐次二次规划（实数）
double ComplexMatrix::Lagrange_multipliers_sequential_quadratic_programming(ComplexMatrix& x0, ComplexMatrix& lambda0,
	double(*f)(ComplexMatrix& x), double(*gradf[])(ComplexMatrix& x), double(*Hf[])(ComplexMatrix& x),
	double(*g[])(ComplexMatrix& x), double(*Jg[])(ComplexMatrix& x), double(*Hg[])(ComplexMatrix& x))
{
	ComplexMatrix xk(x0);							//初始化迭代解x
	ComplexMatrix lambdak(lambda0);					//初始化迭代解lambda
	ComplexMatrix xlk = combine_rows(xk, lambdak);	//初始化非线性方程组的迭代解
	ComplexMatrix sk(xlk.lr, 1, x0.is_real);		//初始化修正步长为零
	//关于牛顿步长的线性方程组
	ComplexMatrix B(x0.lr, x0.lr, x0.is_real);					//储存等式左边的B
	ComplexMatrix J(lambda0.lr, x0.lr, x0.is_real);				//储存g的雅可比矩阵
	ComplexMatrix JT = !J;										//储存g雅可比矩阵的转置
	ComplexMatrix A0(lambda0.lr, lambda0.lr, x0.is_real);		//等式左边矩阵右下零部分
	ComplexMatrix A;											//储存等式左边的系数矩阵
	ComplexMatrix b1(x0.lr, 1, x0.is_real);			//储存等式右边的向量上半部分
	ComplexMatrix b2(lambda0.lr, 1, x0.is_real);	//储存等式右边的向量下半部分
	ComplexMatrix b;								//储存等式右边的向量
	ComplexMatrix Ab;						//储存增广矩阵
	ComplexMatrix temp0(x0.lr, 1);			//临时零向量
	ComplexMatrix temp1;					//计算中间变量
	do {
		//修正解
		xlk = xlk + sk;
		xk = xlk.get_rows(0, x0.lr - 1);
		lambdak = xlk.get_rows(x0.lr, xlk.lr - 1);		cout << xlk << endl;
		//求等式左边的系数矩阵
			//求等式左边左上的B
		for (int i = 0; i < B.lr; i++)
			for (int j = 0; j < B.lc; j++)
			{
				B.c[i][j].re = (Hf[i * x0.lr + j])(xk);
				for (int m = 0; m < lambdak.lr; m++)
					B.c[i][j].re += (lambdak.c[m][0].re * (Hg[m * x0.lr * x0.lr + i * x0.lr + j])(xk));
			}
			//求g的雅可比矩阵及其转置
		for (int i = 0; i < J.lr; i++)
			for (int j = 0; j < J.lc; j++)
				J.c[i][j].re = (Jg[i * x0.lr + j])(xk);
		JT = !J;
			//拼成等式左边的系数矩阵
		A = combine_columns(B, JT);
		temp1 = combine_columns(J, A0);
		A = combine_rows(A, temp1);
		//求等式右边的向量
			//求上半部分
		temp1 = JT * lambdak;
		for (int i = 0; i < b1.lr; i++)
			b1.c[i][0].re = (gradf[i])(xk);
		b1 = b1 + temp1;
		b1 = temp0 - b1;
			//求下半部分
		for (int i = 0; i < b2.lr; i++)
			b2.c[i][0].re = -((g[i])(xk));
			//拼成右边的向量
		b = combine_rows(b1, b2);
		//拼成增广矩阵并解步长
		Ab = combine_columns(A, b);
		Gaussian_elimination_partial_pivoting(Ab, sk);
	} while (sk.vector_2_norm() >= 1e-15);		//步长在误差限内为零时脱出
	xk = xlk.get_rows(0, x0.lr - 1);
	return f(xk);
}