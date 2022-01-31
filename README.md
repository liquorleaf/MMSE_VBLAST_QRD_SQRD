# MMSE_VBLAST_QRD_SQRD
MMSE_VBLAST_QRD_SQRD Space-Time Code Detecting with Complex Matrix.

A simple simulation for:
D. Wubben, R. Bohnke, V. Kuhn and K. -. Kammeyer, "MMSE extension of V-BLAST based on sorted QR decomposition," 2003 IEEE 58th Vehicular Technology Conference. VTC 2003-Fall (IEEE Cat. No.03CH37484), 2003, pp. 508-512 Vol.1, doi: 10.1109/VETECF.2003.1285069.

# Complex Matirx Methods

## 运算符重载

矩阵加、减、数乘、乘法、转置(!)

## 常用行列操作

交换两行、交换两列、交换两列中的一部分行、提取行、提取列、提取连续的一些行、提取子矩阵、
行数相同的矩阵沿列合并、列数相同的矩阵沿行合并、去除指定的一列

## 解三角阵

前向/后向带入消元得到解x（实数/复数）

## 求逆

方阵求逆（复数）、求Moore Penrose伪逆（复数）

## 求范数

向量2-范数（欧几里得范数）（复数）

## 高斯消元

部分选主元高斯消元（实数/复数）

## QR分解

解QR分解后的增广矩阵表示的超定方程（实数）、增广矩阵的Householder变换（实数）、增广矩阵的吉文斯旋转变换（实数）、
古典格拉姆-施密特正交化的QR分解（实数）、改进格拉姆-施密特正交化的QR分解（实数/复数）、有选择的改进格拉姆-施密特正交化的QR分解（复数）

## 插值

牛顿插值获得多项式（实数）、牛顿插值多项式秦九韶算法算值（实数）、牛顿插值多一个点得到新的多项式（实数）、牛顿插值获得多项式（递归方法）（输入参数为最大下标）（实数）

## 优化

* 使用函数指针接收导函数

最速下降法（实数）、牛顿法（实数）、阻尼牛顿法（实数）、无约束优化的BFGS割线修正法（实数）、无约束优化的共轭梯度法（实数）beta由Fletcher-Reeves公式给出、无约束优化的共轭梯度法（实数）beta由Polak-Ribiere公式给出、非线性最小二乘：高斯-牛顿法（实数）、等式约束优化：拉格朗日乘子法，逐次二次规划（实数）
