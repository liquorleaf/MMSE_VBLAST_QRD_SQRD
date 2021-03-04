#pragma once

#include<iostream>
using namespace std;
//复数类
class Complex
{
public:
	//构造函数
	Complex(double r = 0.0, double i = 0.0) { re = r; im = i; }

	//运算符重载
	friend Complex operator + (Complex& c1, Complex& c2);
	friend Complex operator + (double& d1, Complex& c2);
	friend Complex operator + (Complex& c1, double& d2);

	friend Complex operator - (Complex& c1, Complex& c2);
	friend Complex operator - (double& d1, Complex& c2);
	friend Complex operator - (Complex& c1, double& d2);

	friend Complex operator * (Complex& c1, Complex& c2);
	friend Complex operator * (double& d1, Complex& c2);
	friend Complex operator * (Complex& c1, double& d2);

	friend Complex operator / (Complex& c1, Complex& c2);
	friend Complex operator / (double& d1, Complex& c2);
	friend Complex operator / (Complex& c1, double& d2);

	//虚数单位用j表示
	friend ostream& operator << (ostream& out, Complex& c1);
	friend istream& operator >> (istream& in, Complex& c1);

	Complex operator !();

	//常用操作
	//求模值
	double modulus();
public:		//为了单独访问实部
	double re, im;
};