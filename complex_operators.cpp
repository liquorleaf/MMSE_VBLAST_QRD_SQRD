#include<iostream>
#include<iomanip>
#include<cmath>
#include"complex_operators.h"

//单个复数的运算符重载
//加号+
Complex operator + (Complex& c1, Complex& c2)
{
	Complex c3;
	c3.re = c1.re + c2.re;
	c3.im = c1.im + c2.im;
	return c3;
}
Complex operator + (double& d1, Complex& c2)
{
	Complex c3;
	c3.re = d1 + c2.re;
	c3.im = c2.im;
	return c3;
}
Complex operator + (Complex& c1, double& d2)
{
	Complex c3;
	c3.re = c1.re + d2;
	c3.im = c1.im;
	return c3;
}
//减号-
Complex operator - (Complex& c1, Complex& c2)
{
	Complex c3;
	c3.re = c1.re - c2.re;
	c3.im = c1.im - c2.im;
	return c3;
}
Complex operator - (double& d1, Complex& c2)
{
	Complex c3;
	c3.re = d1 - c2.re;
	c3.im = -c2.im;
	return c3;
}
Complex operator - (Complex& c1, double& d2)
{
	Complex c3;
	c3.re = c1.re - d2;
	c3.im = c1.im;
	return c3;
}
//乘号*
Complex operator * (Complex& c1, Complex& c2)
{
	Complex c3;
	c3.re = c1.re * c2.re - c1.im * c2.im;
	c3.im = c1.re * c2.im + c1.im * c2.re;
	return c3;
}
Complex operator * (double& d1, Complex& c2)
{
	Complex c3;
	c3.re = d1 * c2.re;
	c3.im = d1 * c2.im;
	return c3;
}
Complex operator * (Complex& c1, double& d2)
{
	Complex c3;
	c3.re = c1.re * d2;
	c3.im = c1.im * d2;
	return c3;
}
//除号/
Complex operator / (Complex& c1, Complex& c2)
{
	Complex c3;
	double temp = c2.re * c2.re + c2.im * c2.im;
	c3.re = (c1.re * c2.re + c1.im * c2.im) / temp;
	c3.im = (c1.im * c2.re - c1.re * c2.im) / temp;
	return c3;
}
Complex operator / (double& d1, Complex& c2)
{
	Complex c3;
	double temp = c2.re * c2.re + c2.im * c2.im;
	c3.re = d1 * c2.re / temp;
	c3.im = -d1 * c2.im / temp;
	return c3;
}
Complex operator / (Complex& c1, double& d2)
{
	Complex c3(c1.re / d2, c1.im / d2);
	return c3;
}

//输出<<,输入>>
ostream& operator << (ostream& out, Complex& c1)
{
	if (c1.im >= 0)
		out << c1.re << "+" << c1.im << "j";
	else if (c1.im < 0)
		out << c1.re << c1.im << "j";
	return out;
}
istream& operator >> (istream& in, Complex& c1)	//先输入实部，后输入虚部
{
	in >> c1.re >> c1.im;
	return in;
}
//共轭!
Complex Complex::operator !()
{
	Complex c1;
	c1.re = this->re;
	c1.im = -(this->im);
	return c1;
}

//求模值
double Complex::modulus()
{
	return sqrt(re * re + im * im);
}