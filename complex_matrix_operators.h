#pragma once

#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<fstream>
#include"complex_operators.h"
//复数矩阵类
class ComplexMatrix
{
public:
	//构造函数
	explicit ComplexMatrix(int m = 1, int n = 1, bool flag = 0, Complex** incm = NULL);
	ComplexMatrix(double** a, int m, int n);
	ComplexMatrix(const ComplexMatrix& A);
	//析构函数
	~ComplexMatrix();
	//清除
	void clear();
	//*************************************初始化为特殊矩阵*************************************
	//单位阵
	ComplexMatrix make_eyes(int n);

	//*************************************运算符重载******************************************
	friend ComplexMatrix operator + (ComplexMatrix& cm1, ComplexMatrix& cm2);
	friend ComplexMatrix operator + (ComplexMatrix& cm1, double**& d1);
	friend ComplexMatrix operator + (double**& d1, ComplexMatrix& cm1);

	friend ComplexMatrix operator - (ComplexMatrix& cm1, ComplexMatrix& cm2);
	friend ComplexMatrix operator - (ComplexMatrix& cm1, double**& d1);
	friend ComplexMatrix operator - (double**& d1, ComplexMatrix& cm1);

	friend ComplexMatrix operator * (ComplexMatrix& cm1, ComplexMatrix& cm2);
	friend ComplexMatrix operator * (ComplexMatrix& cm1, double**& d1);
	friend ComplexMatrix operator * (double**& d1, ComplexMatrix& cm1);

	friend ComplexMatrix operator * (double& d1, ComplexMatrix& cm1);
	friend ComplexMatrix operator * (ComplexMatrix& cm1, double& d1);
	friend ComplexMatrix operator * (Complex& d1, ComplexMatrix& cm1);
	friend ComplexMatrix operator * (ComplexMatrix& cm1, Complex& d1);

	friend ostream & operator << (ostream& out, ComplexMatrix& cm1);
	friend istream & operator >> (istream& in, ComplexMatrix& cm1);

	ComplexMatrix operator !();
	ComplexMatrix& operator = (const ComplexMatrix A);

	//*************************************常用行列操作****************************************
	//换行
	void exchange_row(int i1, int i2);
	//换列
	void exchange_column(int j1, int j2);
	//换列的一个范围内的行
	void exchange_some_rows_of_column(int j1, int j2, int i1, int i2);
	//得到行（参数为存储下标）
	ComplexMatrix get_row(int i);
	//得到列（参数为存储下标）
	ComplexMatrix get_column(int j);
	//得到连续的许多行（参数为存储下标范围）
	ComplexMatrix get_rows(int i1, int i2);
	//得到子矩阵（参数为存储下标范围）
	ComplexMatrix get_sub_matrix(int i1, int i2, int j1, int j2);

	//列合并：行数相同，拼成更多列的一个矩阵
	ComplexMatrix combine_columns(ComplexMatrix& A, ComplexMatrix& B);
	//行合并：列数相同，拼成更多行的一个矩阵
	ComplexMatrix combine_rows(ComplexMatrix& A, ComplexMatrix& B);
	//去除指定的一列（参数为存储下标）
	ComplexMatrix column_delete(int j);

	//*****************************************解三角阵****************************************
	//前向带入消元得到解x（实数）
	void forward_substitution(ComplexMatrix& A_b, ComplexMatrix& x);
	//前向带入消元得到解x（复数）
	void forward_substitution__Complex(ComplexMatrix& A_b, ComplexMatrix& x);
	//后向带入消元得到解x（实数）
	void backward_substitution(ComplexMatrix& A_b, ComplexMatrix& x);
	//后向带入消元得到解x（复数）
	void backward_substitution__Complex(ComplexMatrix& A_b, ComplexMatrix& x);

	//*****************************************求逆********************************************
	//方阵求逆（复数）
	ComplexMatrix square_inverse();
	//求伪逆（复数）
	ComplexMatrix Moore_Penrose_pseudo_inverse();

	//*****************************************求范数******************************************
	//向量2-范数（欧几里得范数）（复数）
	double vector_2_norm();

	//*****************************************高斯消元****************************************
	//部分选主元高斯消元（实数）
	void Gaussian_elimination_partial_pivoting(ComplexMatrix& A_b, ComplexMatrix& x);
	//部分选主元高斯消元（复数）
	void Gaussian_elimination_partial_pivoting__Complex(ComplexMatrix& A_b, ComplexMatrix& x);

	//*****************************************QR分解*****************************************
	//解QR分解后的增广矩阵表示的超定方程（实数）
	void solution_of_augmentedMatrix_after_QR(ComplexMatrix& Ab, ComplexMatrix& x);
	//增广矩阵的Householder变换（实数）
	void Householder_QR_augmented(ComplexMatrix& A);
	//增广矩阵的吉文斯旋转变换（实数）
	void Givens_QR(ComplexMatrix& A);
	//古典格拉姆-施密特正交化的QR分解（实数）
	ComplexMatrix Gram_Schmidt_QR_classical(ComplexMatrix& Q);
	//改进格拉姆-施密特正交化的QR分解（实数）
	ComplexMatrix Gram_Schmidt_QR_modified(ComplexMatrix& Q);
	//改进格拉姆-施密特正交化的QR分解（复数）
	ComplexMatrix Gram_Schmidt_QR_modified__Complex(ComplexMatrix& Q);
	//有选择的改进格拉姆-施密特正交化的QR分解（复数）
	//算法来源：大二上《数值计算方法》涉及的论文 Efficient Algorithm for Detecting Layered Space-Time Codes
	ComplexMatrix sorted_Gram_Schmidt_QR_modified__Complex(ComplexMatrix& Q, int* S);

	//*****************************************插值********************************************
	//牛顿插值获得多项式（实数）
	ComplexMatrix Newton_interpolation_get_polynomial(ComplexMatrix& t, ComplexMatrix& y);
	//牛顿插值多项式秦九韶算法算值（实数）
	double Newton_interpolation_get_value(double t, ComplexMatrix& t_x, ComplexMatrix& pi);
	//牛顿插值多一个点得到新的多项式（实数）
	ComplexMatrix Newton_interpolation_add_one_point(double& x_new, double& y_new, ComplexMatrix& t, ComplexMatrix& pi);
	//牛顿插值获得多项式（递归方法）（输入参数为最大下标）（实数）
	ComplexMatrix Newton_interpolation_get_polynomial_recursive(ComplexMatrix& t, ComplexMatrix& y, int count);

	//*****************************************优化********************************************
	//优化：最速下降法（实数）
	double Steepest_Descent(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x));
	//优化：牛顿法（实数）
	double Newton_unconstrained_optimization
	(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x), double(*Hessian[])(ComplexMatrix& x));
	//优化：阻尼牛顿法（实数）
	double damped_Newton_unconstrained_optimization
	(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x), double(*Hessian[])(ComplexMatrix& x));
	//无约束优化的BFGS割线修正法（实数）
	double BFGS_optimization(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x));
	//无约束优化的共轭梯度法（实数）beta由Fletcher-Reeves公式给出
	double Conjugate_Gradient_F_R(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x));
	//无约束优化的共轭梯度法（实数）beta由Polak-Ribiere公式给出
	double Conjugate_Gradient_P_R(ComplexMatrix& x0, double(*f)(ComplexMatrix& x), double(*grad[])(ComplexMatrix& x));
	//非线性最小二乘：高斯-牛顿法（实数）
	ComplexMatrix Gauss_Newton_nonlinear_least_squares
	(ComplexMatrix& x0, const ComplexMatrix& t, const ComplexMatrix& y, 
		double(*f)(double, ComplexMatrix&), double(*Jacobian[])(double, ComplexMatrix&));
	//等式约束优化：拉格朗日乘子法，逐次二次规划（实数）
	double Lagrange_multipliers_sequential_quadratic_programming(ComplexMatrix& x0, ComplexMatrix& lambda0,
		double(*f)(ComplexMatrix& x), double(*gradf[])(ComplexMatrix& x),double(*Hf[])(ComplexMatrix& x),
		double(*g[])(ComplexMatrix& x), double(*Jg[])(ComplexMatrix& x),  double(*Hg[])(ComplexMatrix& x));
public:
	int lr;		//行数
	int lc;		//列数
	Complex** c;	//元素
	bool is_real;	//标记是否是实数矩阵
};