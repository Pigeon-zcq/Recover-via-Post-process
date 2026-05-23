描述：
Gegenbauer文件夹中用 Fortran 语言写的代码，主要功能是：
由
\[
f(x) = a(x) + b(x)(1 + x)^s, \quad x \in [-1, 1] \tag{1.1}
\]
Fourier 系数或配点值，应用坐标变换到 y，然后计算 Gegenbauer 多项式投影系数从而重构 f 
相关参考文献在 Chen和Shu - 2014 - Recovering Exponential Accuracy from Collocation Point Values of Smooth Functions with End-Point Sin.md 和 Chen和Shu - 2015 - Recovering Exponential Accuracy in Fourier Spectral Methods Involving Piecewise Smooth Functions wit.md 中。
要求：
根据当前文件的 Fortran 代码，转成 matlab 代码，按照参考文献中的重构算法
