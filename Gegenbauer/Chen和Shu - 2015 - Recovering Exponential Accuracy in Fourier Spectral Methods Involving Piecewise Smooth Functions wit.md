# Recovering Exponential Accuracy in Fourier Spectral Methods Involving Piecewise Smooth Functions with Unbounded Derivative Singularities

**Zheng Chen**¹ · **Chi-Wang Shu**²

Received: 18 July 2014 / Revised: 24 February 2015 / Accepted: 26 February 2015  
Published online: 26 March 2015  
© Springer Science+Business Media New York 2015

## Abstract
Fourier spectral methods achieve exponential accuracy both on the approximation level and for solving partial differential equations, if the solution is analytic. If the solution is discontinuous but piecewise analytic up to the discontinuities, Fourier spectral methods produce pointwise accuracy, but still maintain exponential accuracy after post-processing (Gottlieb and Shu in SIAM Rev 30:644–668, 1997). In Chen and Shu (J Comput Appl Math 265:83–95, 2014), an extended technique is provided to recover exponential accuracy for functions which have end-point singularities, from the knowledge of point values on standard collocation points. In this paper, we develop a technique to recover exponential accuracy from the first \( N \) Fourier coefficients of functions which are analytic in the open interval but have unbounded derivative singularities at end points. With this post-processing method, we are able to obtain exponential accuracy of spectral methods applied to linear transport equations involving such functions.

**Keywords**: Spectral method · Exponential accuracy · Fourier coefficients · Gegenbauer expansion · Transport equation · Variable coefficients · Singular initial conditions · Noise

## 1 Introduction
In this paper, we are concerned with the accuracy of spectral methods when applied to problems involving piecewise smooth functions with unbounded derivative singularities. We

Research supported by NSF Grants DMS-1112700 and DMS-1418750, and AFOSR Grant F49550-12-1-0399.

**Chi-Wang Shu**  
shu@dam.brown.edu

**Zheng Chen**  
zchen@iastate.edu

1. Department of Mathematics, Iowa State University, Ames, IA 50011, USA  
2. Division of Applied Mathematics, Brown University, Providence, RI 02912, USA

investigate the issue of overcoming the Gibbs phenomenon, which describes how a global spectral approximation of a piecewise analytic function behaves at the jump discontinuity. A prototype is to use Fourier series to approximate an analytic but non-periodic function \( u(x) \) on interval \([-1, 1]\), which has discontinuities at the boundaries of the interval when extended periodically with period 2. The Fourier partial sum using the first \( 2N + 1 \) modes

\[
u_N(x) = \sum_{|k| \leq N} \hat{u}_k e^{ik\pi x},
\]

with Fourier coefficients \( \hat{u}_k \) defined by

\[
\hat{u}_k = \frac{1}{2} \int_{-1}^{1} u(x)e^{-ik\pi x} \, dx,
\]

has large oscillations near the jumps, which are not improved as the number of terms in the partial sum increases. In smooth regions away from the discontinuities, convergence is only first order. Therefore, there is no convergence in the maximum norm. This is the so-called Gibbs phenomenon.

In [9–12, 14], Gottlieb et al. developed a general framework to overcome this difficulty. This technique recovers exponential accuracy in the maximum norm for any (sub-)interval of analyticity (up to and including the boundaries of this interval), from the knowledge of either the first \( N \) spectral expansion coefficients, or the point values at \( N \) standard collocation points. This means that exponential accuracy is recovered at all points, including at the actual discontinuity points (the left and right limits at these points), if the locations of these discontinuity points are known. If the locations of the discontinuity points are not known exactly but are known to be within certain fixed intervals, then exponential accuracy can be recovered from any interval which does not overlap with these fixed intervals containing the discontinuities. In this framework, an important tool is the set of Gegenbauer polynomials, which are orthogonal in the interval \([-1, 1]\] with the weight \( (1 - x^2)^{\lambda - \frac{1}{2}} \). The key in this technique is that the parameter \( \lambda \) in the weight function as well as the number of terms \( m \) retained in the Gegenbauer expansion should both be chosen proportional to \( N \). For a review, we refer to [13]. These methods are widely used as "reconstruction" or "post-processing" techniques to recover exponential accuracy for point values or coefficients based on the spectral approximation, such as in recovering high order information of the discontinuous solutions of scalar nonlinear hyperbolic PDEs [18], and in the simulation of sophisticated problems [8]. These techniques have been successfully applied to the field of image reconstruction [2–5] as well. The Gegenbauer basis has also been successful in recovering lost order of accuracy in other types of approximations, such as weighted essentially non-oscillatory (WENO) solutions of hyperbolic PDEs [15], and in radial basis functions approximations of linear and nonlinear hyperbolic PDEs [17].

Besides piecewise analytic functions, functions with end-point singularities exist in many applications. Most fractional differential equations have singular solutions. Many standard numerical methods solving fractional differential problems give poor accuracy, due to the lack of regularities. Therefore, it is important to provide a way to recover accuracy of solutions as well as to obtain high order accuracy at the approximation level. In [1], Adcock et al. focused on the approximation of functions which are analytic on a compact interval except at the end-points, and utilized variable transform methods. They introduced two new mappings from the original interval to either semi-infinite or infinite interval, and provided approximation procedure from sampling information on the new region. The two new mappings,

J Sci Comput (2015) 65:1145–1165

compared with the standard transformations, vastly improve resolution power, and achieve root exponential decays, with proper choice of parameters in the mappings.

We are interested in developing post-processors to recover high order accuracy for such functions as well. This task is significantly more difficult than the recovery of accuracy for piecewise analytic functions. In particular, the extension of the technique in [13] to functions with end-point singularities is highly non-trivial. In [7], we made this extension to handle spectral collocation methods for such functions. The reconstruction procedure is performed on functions of the following form

$$ f(x) = a(x) + b(x)(1 + x)^s, \quad x \in [-1, 1] \tag{1.3} $$

where \( s \) is a given fractional constant

$$ 0 < s = \frac{p}{q} < 1 \tag{1.4} $$

with relatively prime integers \( p \) and \( q \), and \( a(x) \) and \( b(x) \) are both analytic but unknown functions. Such functions lack regularities, and the derivatives blow up at the left end point. With this extension, exponential accuracy can be obtained from standard collocation point values of such functions, which is different from sampling on the mapped region as in [1], by properly choosing the parameters \( \lambda \) and \( m \) to be linearly dependent on \( N \). A crucial modification of the choice of parameters and a more refined estimate were necessary to balance the terms in the truncation error for this analysis.

In this study, we are interested in recovering high order accuracy from the first \( 2N + 1 \) Fourier coefficients for functions in the form (1.3). The objective is to extract the hidden information from the truncated Fourier series (1.1) and recover exponentially accurate point values at every point including at the singularities.

As in [7], we assume that the analytic functions \( a(x) \) and \( b(x) \), denoted generically as \( c(x) \), satisfy the following condition.

**Assumption 1.1** There exists a constant \( \rho \geq 1 \) and a constant \( C(\rho) \) such that, for every \( k \geq 0 \),

$$ \max_{-1 \leq x \leq 1} \left| \frac{d^k c(x)}{dx^k} \right| \leq C(\rho) \frac{k!}{\rho^k}. $$

This is a standard assumption for analytic functions, where \( \rho \) is the distance from the interval \( [-1, 1] \) to the nearest singularity of the function \( c(x) \) in the complex plane.

We will use the following one-to-one transformation between \( x \in [-1, 1] \) and \( y \in [-1, 1] \):

$$ (2^{q-1}(1 + x))^{\frac{1}{q}} = 1 + y \tag{1.5} $$

where \( q \) is defined in (1.4).

The function \( F(y) = f(x(y)) \) of the variable \( y \) has its usual Gegenbauer expansion under the basis \( \{C_l^q(y)\} \):

$$ f(x(y)) = F(y) = \sum_{l=0}^{\infty} \hat{f}^q(l) C_l^q(y) $$

with the Gegenbauer coefficients \( \hat{f}^q(l) \) given by

$$ \hat{f}^q(l) = \frac{1}{\hat{h}_l^q} \int_{-1}^{1} (1 - y^2)^{\lambda - \frac{1}{2}} F(y) C_l^q(y) dy \tag{1.6} $$

where the precise value of the normalization constant \( h^1 \) will be given later by (2.4) in Definition 2.1.

Our goal is to find a good approximation to the first \( m \sim N \) Gegenbauer coefficients \( \hat{f}^k(l) \) in (1.6), denoted as \( \hat{g}^k(l) \) [defined later in (3.1)], from the given Fourier coefficients. We will then obtain the approximation of \( f(x) \) using these \( m \sim N \) terms of its Gegenbauer expansion:

\[
f^{m, \lambda}_N(x) = \sum_{l=0}^{m} \hat{g}^k(l) C^k_l(\varphi(x)).
\]

With proper choice of the parameters \( \lambda \) and \( m \), the error between the reconstructed approximation \( f^{m, \lambda}_N(x) \) and the function \( f(x) \), measured in the maximum norm, decays exponentially as \( N \) increases. Therefore, the reconstruction method provides a way to post-process functions with such singularities from the first \( 2N + 1 \) accurate Fourier coefficients.

Another major concern in the reconstruction methods is the possible existence of noise in the data. The Gegenbauer reconstruction techniques in [13] work well with noise. This is because the noise is projected to be very small in the Gegenbauer basis. Similar robustness can also be observed in the new reconstruction methods for piecewise smooth functions with end-point singularities.

The analysis of the error from reconstruction is divided into two parts: the truncation error and the regularization error. The truncation error measures the difference between the exact Gegenbauer coefficients of \( f(x(\varphi)) \) with \( \lambda \sim N \), and the approximate Gegenbauer coefficients \( \hat{g}^k(l) \) obtained by using the truncated Fourier series. This will be investigated in Sect. 3. In relation to the collocation case as in [7], the Galerkin case is more difficult for such singularity cases, which involves different and more complicated analysis. The regularization error measures the difference between the Gegenbauer expansion with \( \lambda \sim N \), using the first \( m \sim N \) Gegenbauer coefficients, and the function itself. This error is estimated in Sect. 4. The results for the reconstruction are summarized in Theorem 4.3 in Sect. 4. We also give an analysis of the Fourier Galerkin methods for solving initial problems of linear hyperbolic time-dependent partial differential equations with reconstructions in Sect. 5. Section 6 contains several numerical examples to illustrate our results and robustness to noise. In Sect. 2, we shall give several useful preliminary properties and estimates. Concluding remarks are given in Sect. 7.

Throughout this paper, we will use \( C \) to denote a generic constant either independent of the growing parameters, or depending on them at most in polynomial growth. The details will be indicated clearly in the text. These constants may not take the same value at different places.

## 2 Preliminaries

In this section, we will introduce the Gegenbauer polynomials and discuss some of their asymptotic behavior (see Bateman [6]). Then we will give some estimates as preparation for the error estimates in Sect. 3.

### Definition 2.1
The Gegenbauer polynomial \( C^k_n(x) \), for \( \lambda \geq 0 \), is defined by

\[
(1 - x^2)^{\frac{\lambda - \frac{1}{2}}{2}} C^k_n(x) = \frac{(-1)^n}{2^n n!} G(\lambda, n) \frac{d^n}{dx^n} \left[ (1 - x^2)^{n + \frac{\lambda - \frac{1}{2}}{2}} \right]
\]

where \( G(\lambda, n) \) is given by

\[
G(\lambda, n) = \frac{\Gamma(\lambda + \frac{1}{2})\Gamma(n + 2\lambda)}{\Gamma(2\lambda)\Gamma(n + \lambda + \frac{1}{2})} \tag{2.1}
\]

for \( \lambda > 0 \), by

\[
G(0, n) = \frac{2\sqrt{\pi(n - 1)!}}{\Gamma(n + \frac{1}{2})} \tag{2.2}
\]

for \( \lambda = 0 \) and \( n \geq 1 \), and by 

\[
G(0, 0) = 1 \tag{2.3}
\]

for \( \lambda = 0 \) and \( n = 0 \). Notice that by this standardization, \( C^0_n(x) \) is defined by:

\[
C^0_n(x) = \lim_{\lambda \to 0^+} \frac{1}{\lambda} C^1_n(x) = \frac{2}{n} T_n(x), \quad n > 0; \quad C^0_0(x) = 1, 
\]

where \( T_n(x) \) are the Chebyshev polynomials.

Under this definition we have, for \( \lambda > 0 \),

\[
C^\lambda_n(1) = \frac{\Gamma(n + 2\lambda)}{n!\Gamma(2\lambda)}; \tag{2.4}
\]

for \( \lambda = 0 \) and \( n \geq 1 \),

\[
C^0_n(1) = \frac{2}{n}; 
\]

for \( \lambda = 0 \) and \( n = 0 \),

\[
C^0_0(1) = 1; 
\]

and 

\[
|C^\lambda_n(x)| \leq C^\lambda_n(1), \quad -1 \leq x \leq 1. \tag{2.5}
\]

The Gegenbauer polynomials are orthogonal under their weight function \( (1 - x^2)^{\lambda - \frac{1}{2}} \):

\[
\int_{-1}^{1} (1 - x^2)^{\lambda - \frac{1}{2}} C^k_\lambda(x) C^\lambda_n(x) \, dx = \delta_{k,n} h^\lambda_n, 
\]

where, for \( \lambda > 0 \),

\[
h^\lambda_n = \pi \frac{1}{2} C^1_n(1) \frac{\Gamma(\lambda + \frac{1}{2})}{\Gamma(\lambda)(n + \lambda)}; \tag{2.6}
\]

for \( \lambda = 0 \) and \( n \geq 1 \),

\[
h^0_n = \frac{2\pi}{n^2}; 
\]

for \( \lambda = 0 \) and \( n = 0 \),

\[
h^0_0 = \pi. 
\]

We will need to use the Stirling's formula and the estimate of \( h^\lambda_n \) for the asymptotics of the Gegenbauer polynomials for large \( n \) and \( \lambda \).

**Lemma 2.2** We have the Stirling's formula

\[
(2\pi)^{\frac{1}{2}} x^{x + \frac{1}{2}} e^{-x} \leq \Gamma(x + 1) \leq (2\pi)^{\frac{1}{2}} x^{x + \frac{1}{2}} e^{-x + \frac{1}{1} x}, \quad x \geq 1. \tag{2.7}
\]

# Lemma 2.3
There exists a constant \( C \) independent of \( \lambda \) and \( n \) such that

\[
C^{-1} \frac{\lambda^{\frac{1}{2}}}{n + \lambda} C_n^{\hat{n}} \leq h_n \leq C \frac{\lambda^{\frac{1}{2}}}{n + \lambda} C_n^{\hat{n}}(1) \tag{2.6}
\]

We would need to estimate \( \left\| \frac{d^{t}}{dx^{t}} \left\{ \frac{d}{dy} (1 - y^2)^{\frac{1}{2} + \frac{1}{2}} \frac{dy}{dx} \right\} \right\|_{\infty} \), therefore we need the following preliminaries first.

## Remark 2.4
\[
\frac{d}{dy} (1 - y^2)^{\frac{1}{2} + \frac{1}{2}} \frac{dy}{dx}
\]
has up to \( t \)-th derivatives in \( x \), where \( t = \left[ \frac{\lambda + \frac{1}{2}}{q} \right] - 1 \), the largest integer below \( \frac{\lambda + \frac{1}{2}}{q} - 1 \).

It is easy to observe that

\[
\frac{d^{m}}{dx^{n}} \left( (1 - y(x)^2)^{\frac{1}{2} + \frac{1}{2}} \right) = A^{n} y_1^{\hat{n}} y_2^{\hat{n}} y_3, \quad 0 \leq n \leq t + 1 \tag{2.7}
\]

where 

\[
A = \frac{2q}{2q} , \quad y_1 = (1 - y(x)^2)^{\frac{1}{2} - qn}, \quad y_2 = (1 - y(x))^{n(q-1)},
\]

and \( y_3 \) satisfies the following recursive relation:

\[
Y_3^0 = 1
\]

\[
Y_3^{n+1} = -\left[ (2l + 2\lambda - qn - n - 1)y + n(q - 1) \right] Y_3^{n} + (1 - y^2) \frac{d}{dy} Y_3^{n}, \quad 0 \leq n \leq t.
\]

It is easy to show that \( Y_3^{n} \) is an \( n \)-th degree polynomial of \( y \). We have the following estimate on \( Y_3^{n} \).

# Lemma 2.5
We have, for \( 0 \leq n \leq t + 1 \),

\[
|Y_3^{n}| \leq (2l + 2\lambda)^{n}, \quad y \in [-1, 1] \tag{2.8}
\]

**Proof** The proof can be found in [7]. In the proof, \( Y_3^{n} \) is rewritten as \( Y_3^{n} = \sum_{i=0}^{n} a_i y^i \) and we denote \( S_n = \sum_{i=0}^{n} |a_i| \). The proof also provides an estimate for \( S_n \):

\[
S_n \leq (2l + 2\lambda)^{n} \tag{2.9}
\]

This will be used later in the estimation.

For \( l \geq 1 \),

\[
\frac{d^{t}}{dx^{t}} \left\{ \frac{d^{t}}{dx^{t}} \left( (1 - y^2(x))^{\frac{1}{2} + \frac{1}{2}} \right) \right\} = \frac{d^{t-1}}{dy^{t-1}} \left[ \left( (1 - y^2(x))^{\frac{1}{2} + \frac{1}{2}} \right) \right] \tag{2.10}
\]

# Lemma 2.6
We have the following estimate for \( l \geq 1 \),

\[
\left| \frac{d^{t-1}}{dy^{t-1}} \left\{ \left( (1 - y^2(x))^{\frac{1}{2} + \frac{1}{2}} \right)^{t + 1} \right\} \right| \leq C A^{t + 1} l^{t + 1} (l + \lambda)^{t}, \quad y \in [-1, 1] \tag{2.11}
\]

**Proof** 

\[
\frac{d^{t-1}}{dy^{t-1}} \left\{ \left( (1 - y^2(x))^{\frac{1}{2} + \frac{1}{2}} \right)^{t + 1} \right\} = A^{t + 1} \frac{d^{t-1}}{dy^{t-1}} \left\{ Y_1^{t + 1} Y_2^{t + 1} Y_3^{t + 1} \right\}
\]

For simplicity, we denote 

$$
\frac{d^i}{dy^i} \left( Y_{1}^{l+1} Y_{2}^{l+1} Y_{3}^{l+1} \right) = X_{1} X_{2}^{l}, \quad 0 \leq i < l - 1
$$

where, 

$$
X_{1}^{i} = (1 - y^2)^{l+\lambda - \frac{1}{2} - q(t+1)-i}
$$

and $X_{2}^{i}$ satisfies the following recursive relationship:

$$
X_{2}^{0} = Y_{2}^{l+1} Y_{3}^{l+1}
$$

$$
X_{2}^{i+1} = -\left[ 2(l + \lambda - q(t + 1) - 1) - \frac{1}{2} - i \right] X_{2}^{i} + (1 - y^2) \frac{d}{dy} X_{2}^{i}, \quad 0 \leq i < l - 2.
$$

It is easy to find out that $X_{2}^{i}$ is a polynomial of degree $(t + 1)q + i$. We only need to prove that 

$$
|X_{2}^{l-1}| \leq C^{2} l^{+\lambda}(l + \lambda)^{t} \tag{2.12}
$$

We can rewrite it as 

$$
X_{2}^{i} = \sum_{j=(t+1)q+i}^{j=0} \beta_{j} y^{j}. 
$$

In order to get the upper bound for $X_{2}^{l-1}$ as in (2.12), we need to estimate $X_{2}^{0}$ first. Let us rewrite $Y_{3}^{l+1} = \sum_{i=0}^{r+1} a_{i} y^{i}$ with its coefficients $a_{i}$.

$$
X_{2}^{0} = (1 - y)(t+1)(q-1) Y_{3}^{l+1}
$$

$$
= \sum_{j=0}^{(r+1)(q-1)} \left( \binom{(t + 1)(q - 1)}{j} (-y)^{j} \right) \sum_{i=0}^{r+1} a_{i} y^{i}
$$

$$
= \sum_{j=0}^{(r+1)(q-1)} \sum_{i=0}^{r+1} \left( \binom{(t + 1)(q - 1)}{j} (-1)^{j} a_{i} y^{i+j} \right)
$$

Then we measure the sum of the coefficients of $X_{2}^{0}$,

$$
\sum_{k=0}^{(t+1)q} |\beta_{k}| \leq \sum_{j=0}^{(t+1)(q-1)} \sum_{i=0}^{r+1} \left( \binom{(t + 1)(q - 1)}{j} |a_{i}| \right)
$$

$$
= 2^{(t+1)(q-1)} S_{t+1} \leq 2^{(t+1)q}(l + \lambda)^{t+1}
$$

Using induction, we can easily get similar result for $X_{2}^{l-1}$,

$$
\sum_{k=0}^{(t+1)q+t-1} |\beta_{k}| \leq 2^{(t+1)q+t-1}(l + \lambda)^{t} 
$$

which implies that 

$$
|X_{2}^{l-1}| \leq 2^{(t+1)q+t-1}(l + \lambda)^{t} \leq C^{2} l^{+t}(l + \lambda)^{t}
$$

Thus, we complete the proof.

Remark 2.7 The result in Lemma 2.6 is also true for \( l = 0 \). Hence we have

\[
\left\| \frac{d^l}{dx^l} \left[ \frac{d}{dy} \left( 1 - y^2(x) \right)^{\lambda - \frac{1}{2}} \frac{dy}{dx} \right] \right\|_{L^\infty} \leq C A^{\tau + 1} 2^{\lambda}(l + \lambda)^{\tau}, \quad l \geq 0 \tag{2.13}
\]

### 3 Truncation Error

In this section, we will establish the error estimate for replacing the Gegenbauer coefficients \( \hat{f}^{\lambda}(l) \) by the new approximated coefficients \( \hat{g}^{\lambda}(l) \), defined later in (3.1), in the Gegenbauer expansion.

Consider the function in the form of

\[
f(x) = a(x) + b(x)(1 + x)^s
\]

where \( s \) is a given constant \( 0 < s < 1 \) with relatively prime positive integers \( p \) and \( q \), and \( a(x) \) and \( b(x) \) are analytic functions satisfying Assumption 1.1.

We assume that the Fourier coefficients \( \hat{f}_n \) \( (-N \leq n \leq N) \) are given. Thus we have its truncated Fourier series

\[
f_N(x) = \sum_{|n| \leq N} \hat{f}_n e^{in\pi x}.
\]

We are interested in recovering the first \( m \) coefficients in the Gegenbauer expansion of \( f(x) \). For the function \( F(y) = f(x(y)) \), we have the usual Gegenbauer expansion with the basis \( \{C_i^{\lambda}\} \):

\[
f(x(y)) = F(y) = \sum_{l=0}^{\infty} \hat{f}^{\lambda}(l) C_i^{\lambda}(y).
\]

where the Gegenbauer coefficients \( \hat{f}^{\lambda}(l) \) are given by (1.6).

The candidate for approximating the Gegenbauer coefficients \( \hat{f}^{\lambda}(l) \) is:

\[
\hat{g}^{\lambda}(l) = \frac{1}{h_l^{\lambda}} \int_{-1}^{1} (1 - y^2)^{\lambda - \frac{1}{2}} f_N(x(y)) C_i^{\lambda}(y) dy. \tag{3.1}
\]

### Definition 3.1

The truncation error is defined as

\[
TE(\lambda, m, N) = \max_{-1 \leq y \leq 1} \left| \sum_{l=0}^{m} \left( \hat{f}^{\lambda}(l) - \hat{g}^{\lambda}(l) \right) C_i^{\lambda}(y) \right| \tag{3.2}
\]

where \( \hat{f}^{\lambda}(l), \hat{g}^{\lambda}(l) \) are defined in (1.6) and (3.1).

In the next two lemmas, we bound the truncation error using the regularity of the function \( M(y) \) as in (3.3), and then in terms of the number of given Fourier coefficients \( N \), the number of Gegenbauer polynomials \( m \) and parameter of Gegenbauer polynomial \( \lambda \).

### Lemma 3.2

The truncation error is bounded by

\[
TE(\lambda, m, N) \leq \frac{C}{(N\pi)^{r-1}} \sum_{l=0}^{m} C_i^{\lambda}(1) \left\| \frac{d^l}{dx^l} M(y) \right\|_{L^\infty}
\]

J Sci Comput (2015) 65:1145–1165

where \( t = \left[ \frac{\lambda + \frac{2}{q}}{q} \right] - 1 \), and  
\( M(y) = (1 - y^2)^{\lambda - \frac{1}{2}} C^1(y) \frac{dy}{dx} \) (3.3)

**Proof**  
We have  
\[
TE(\lambda, m, N) = \max_{-1 \leq y \leq 1} \left| \sum_{l=0}^{m} \frac{C^l(y)}{h_l^1} \int_{-1}^{1} (1 - y^2)^{\lambda - \frac{1}{2}} (f - f_N) \circ x(y) C^l(y) dy \right|
\]
\[
\leq \sum_{l=0}^{m} \frac{C^l(1)}{h_l^1} \int_{-1}^{1} (1 - y^2)^{\lambda - \frac{1}{2}} (f - f_N) \circ x(y) C^l(y) dy
\]
\[
= \sum_{l=0}^{m} \frac{C^l(1)}{h_l^1} \int_{-1}^{1} (1 - y^2)^{\lambda - \frac{1}{2}} \sum_{|n| > N} \hat{f}_n e^{inx} C^l(y) dy
\]
\[
= \sum_{l=0}^{m} \frac{C^l(1)}{h_l^1} \sum_{|n| > N} \hat{f}_n \int_{-1}^{1} e^{inx} M(y) dy
\]
\[
= \sum_{l=0}^{m} \frac{C^l(1)}{h_l^1} \sum_{|n| > N} \hat{f}_n \int_{-1}^{1} e^{inx} \frac{d}{dx} M(y) dx
\]
\[
= \sum_{l=0}^{m} \frac{C^l(1)}{h_l^1} \sum_{|n| > N} \hat{f}_n \int_{-1}^{1} e^{inx} \frac{d}{dx} M(y) dx
\]
\[
\leq C \sum_{|n| > N} \frac{C^l(1)}{h_l^1} \sum_{|n| > N} 1 \cdot \left\| \frac{d}{dx} M(y) \right\|_{L^\infty}
\]
\[
\leq \frac{C}{(N\pi)^{t-1}} \sum_{l=0}^{m} \frac{C^l(1)}{h_l^1} \left\| \frac{d}{dx} M(y) \right\|_{L^\infty}
\]

The definitions of \( \hat{f}^\lambda \) in (1.6) and \( \hat{g}^\lambda(l) \) in (3.1) are used in the first equality; (2.3) is used in the second inequality; in the third equality, the error of the Fourier partial sum \( f(x) - f_N(x) = \sum_{|n| > N} \hat{f}_n e^{inx} \) is used; the substitution (1.5) has been made in the integral in the fourth equality; in the fifth equality, we use integration by parts \( t \) times, and the fact that \( \frac{d}{dx} M(y) \) vanishes at \( y = \pm 1 \) for \( 0 \leq i \leq t - 1 \); since \( f(x) \) is an \( L^2 \)-function, its Fourier coefficients \( \hat{f}_n \) are uniformly bounded, i.e., \( |\hat{f}_n| \leq C \), which is used in the sixth inequality. □

**Lemma 3.3**  
\[
TE(\lambda, m, N) \leq C \frac{(m + 1)2^\lambda A^{t+1}(\lambda + m)^{m+t+1} \Gamma(\lambda) \Gamma(m + 2\lambda)}{m! \Gamma(2\lambda) \Gamma(m + \frac{1}{2})(N\pi)^{t-1}}
\]

**Proof**  
From Lemma 3.2, we have  
\[
TE(\lambda, m, N) \leq C \frac{(N\pi)^{t-1}}{m} \sum_{l=0}^{m} \frac{C^l(1)}{h_l^1} \left\| \frac{d}{dx'} \left[ (1 - y^2)^{\lambda - \frac{1}{2}} C^l(y) \frac{dy}{dx} \right] \right\|_{L^\infty}
\]
\[
= C \frac{(N\pi)^{t-1}}{m} \sum_{l=0}^{m} \frac{C^l(1)}{h_l^1} \frac{1}{2l/h_l^1} \left\| \frac{d}{dx} \frac{d}{dy} (1 - y^2)^{\lambda - \frac{1}{2}} \frac{dy}{dx} \right\|_{L^\infty}
\]
\[
\leq C \frac{A^{t+1}2^{(N\pi)^{t-1}}}{l h_l^1} \sum_{l=0}^{m} C^l(1) G(\lambda, l)(\lambda + t)^{t}
\]

1154
J Sci Comput (2015) 65:1145–1165

\[
= C \frac{A^{t+1}2^{t}\Gamma(\lambda)}{\Gamma(2)(N\pi)^{t-1}} \sum_{l=0}^{m} \Gamma(l + 2\lambda)(l + \lambda)^{t+1} \frac{l!}{l + \frac{1}{2}}
\]
\[
\leq C \frac{A^{t+1}2^{t}\Gamma(\lambda)(m + 2\lambda)(m + \lambda)^{m+t+1}}{m! \Gamma(2\lambda)(m + \frac{1}{2})(N\pi)^{t-1}}
\]

where the definition of Gegenbauer polynomial in (2.1) is used in the second step; Lemma 2.6 and Remark 2.7 are used in the third inequality; (2.1) and (2.4) are used in fourth equality; in the last step, we use the fact that \(\Gamma(l + 2\lambda)(l + \lambda)^{t+1}/l!(l + \lambda + \frac{1}{2})\) is an increasing function of \(l\).

If both \(\alpha\) and \(m\) grow linearly with \(N\), the truncation error can be made to decay exponentially. This result is stated in the following theorem.

**Theorem 3.4** (The exponential decay of the truncation error) Let \(\lambda = \alpha N\), \(m = \beta N\) with \(0 < \alpha, \beta < 1\), then

\[
TE(\alpha N, \beta N, N) \leq C N^{2} q_{T}
\]

with

\[
q_{T} = \frac{A^{\frac{\alpha}{q}} A^{\frac{1}{q}} e^{\beta(2\alpha + \beta)^{2\alpha + \beta}}}{\pi^{\frac{1}{q}} (2\alpha)^{\alpha} \beta^{\beta} (\alpha + \beta)^{(1 - \frac{1}{q})\alpha}}.
\]

When we choose \(\beta = \gamma\), i.e. \(m = \gamma N\), we have

\[
q_{T} = \left( \frac{1}{4} e^{\frac{2}{(2 + \gamma)^{2}}} \alpha^{\frac{1}{q}} \right)^{\alpha}.
\]

If we choose \(\alpha\) to satisfy

\[
\alpha < \frac{\pi}{A} \left( \frac{2\gamma(1 + \gamma)^{1 - \frac{1}{q}}}{e^{\frac{2}{(2 + \gamma)^{2}}}} \right)^{q},
\]

then \(q_{T} < 1\).

**Proof** From lemma 3.3, we use Stirling's formula (2.5) and get

\[
TE(\lambda, m, N) \leq C (m + 1)2^{A^{t + 1}}(m + \lambda)^{m + t + 1}\Gamma(t + 2\lambda) \frac{1}{m!(2\lambda)^{t}\Gamma(m + \lambda + \frac{1}{2})(N\pi)^{t - 1}}
\]
\[
\leq C \frac{e^{m A^{t + 1}(m + 2\lambda)^{m + 2 - \frac{1}{2}}}}{2^{\lambda} \lambda^{m + \lambda - 1} m^{m - \frac{1}{2}}(N\pi)^{t - 1}}
\]

hence

\[
TE(\alpha N, \beta N, N) \leq CN^{2} \frac{A(\alpha + \beta)}{\pi} \left( \frac{e^{\beta(2\alpha + \beta)^{2\alpha + \beta}}}{(2\alpha)^{\alpha} \beta^{\beta}(\alpha + \beta)^{N}} \right)^{\frac{1}{q}}.
\]

\[
\leq CN^{2} \frac{A(\alpha + \beta)}{\pi} \left( \frac{e^{\beta(2\alpha + \beta)^{2\alpha + \beta}}}{(2\alpha)^{\alpha} \beta^{\beta}(\alpha + \beta)^{N}} \right)^{\frac{1}{q}}.
\]

\[
= CN^{2} q_{T}
\]

where we have used \(\lambda^{\frac{2}{q}} - 2 \leq t = \left| \frac{\lambda + \frac{1}{2}}{q} - 1 \right| \leq \frac{1}{q}\) in the second step.

Remark 3.5 With proper choice of the ratios between \( \lambda, m \) and \( N, q_T \) can be a positive number less than 1, and then the upper bound of truncation error is exponentially small. In the proof, there is no effort to optimize the parameters. In practice, \( \lambda \) and \( m \) should be chosen to balance the truncation error as well as the regularization error, which will be analyzed in the next section.

# 4 Regularization Error

Regularization error, the second part of the post-processing error, is caused by approximating \( f(x) = a(x) + b(x)(1+x)^s \) on \([-1, 1]\) by its Gegenbauer expansion based on the Gegenbauer polynomials \( C^l_1(x) \). This has been studied in [7]. We will just quote the result. The function \( F(y) = f(x(y)) = a \cdot x(y) + b \cdot x(y)(1 + y)^2s - p \) of the variable \( y \) is analytic on \( y \in [-1, 1] \). Let us consider the Gegenbauer partial sum of the first \( m \) terms for the function \( f(x) \) given by:

\[
f^{\lambda, m}(x) = \sum_{l=0}^{m} \hat{f}^{\lambda}(l) C^l_1(y(x)), \tag{4.1}
\]

with the Gegenbauer coefficients \( \hat{f}^{\lambda}(l) \) defined by (1.6).

## Definition 4.1

The regularization error is defined by

\[
RE(\lambda, m) = \max_{-1 \leq y \leq 1} \left| f(x) - \sum_{l=0}^{m} \hat{f}^{\lambda}(l) C^l_1(y(x)) \right|
\]

\[
= \max_{-1 \leq y \leq 1} \left| f(x(y)) - \sum_{l=0}^{m} \hat{f}^{\lambda}(l) C^l_1(y) \right|.
\]

When both \( \lambda \) and \( m \) grow linearly with \( N \), we have the following result for the estimate of regularization error [7].

## Theorem 4.2

(The exponential decay of the regularization error) For the function \( f(x) = a(x) + b(x)(1 + x)^s \), with analytic functions \( a(x) \) and \( b(x) \) satisfying Assumption 1.1, if we assume \( \lambda = \alpha N \) and \( m = \gamma N \), then

\[
\max_{-1 \leq y \leq 1} \left| f(x) - \sum_{l=0}^{m} \hat{f}^{\lambda}(l) C^l_1(y(x)) \right| \leq C q_R^N
\]

where

\[
q_R = \left( \frac{(y + 2)^{y + 2}}{\rho y^2 + 2(y + 1)^{y + 1}} \right)^{\alpha}
\]

which is always less than 1.

**Proof** The proof of this theorem can be found in [7]. 

We can combine the estimates for the truncation error and the regularization error and obtain the following theorem about the reconstruction error:

Theorem 4.3 (Removal of the Gibbs Phenomenon) Consider a function in the form of 

$$
f(x) = a(x) + b(x)(1 + x)^s,
$$ 

with given fractional constant \(0 < s < \frac{2}{q} < 1\), and \(a(x)\) and \(b(x)\) are analytic functions satisfying Assumption 1.1. Assume that the Fourier coefficient 

$$
\tilde{f}_n = \frac{1}{2} \int_{-1}^{1} f(x)e^{-i k \pi x} dx
$$ 

are known for \(-N \leq n \leq N\). Let \( \hat{g}^{(l)} \), with \(0 \leq l \leq m\), be the Gegenbauer expansion coefficients of \(f_N(x) = \sum_{|n| \leq N} \hat{f}_n e^{i n \pi x}\), defined in (3.1). Then for \( \lambda = \alpha N \) and \( m = \gamma \) with 

$$
\alpha < \frac{\pi}{A} \left( \frac{2 \gamma (1 + \gamma)^{\frac{1 - l}{q}}}{e^{\gamma}(2 + \gamma)^{2 + \gamma}} \right)^{q},
$$ 

we have 

$$
\left| \max_{-1 \leq x \leq 1} f(x) - \sum_{l=0}^{m} \hat{g}^{(l)} C_l^{\gamma}(y(x)) \right| \leq C \left( q_N \frac{1}{q_T} + q_R^N \right)
$$ 

where 

$$
q_T = \left( \frac{1}{4} e^{\gamma}(2 + \gamma)^{2 + \gamma} \frac{1}{2 \pi^q \gamma (1 + \gamma)^{\frac{1 - 1/q}}} \right)^{\alpha} < 1, \quad q_R = \left( \frac{(\gamma + 2) y^{2} }{ \rho \gamma^{2} + (2 + 1) y^{1}} \right)^{\alpha} < 1.
$$ 

Remark 4.4 In the proof, no attempt has been made to optimize the parameters. We will discuss about the choice of parameters in numerical examples in Sect. 6.

Remark 4.5 For functions with multiple singularities, the reconstructions can be built piece by piece, and on each piece the function has at most one end-point singularity. The technique only needs to be modified by scaling. It is easy to prove that the reconstruction for each piece can also be exponentially accurate with proper parameters.

5 Error Analysis for Fourier Galerkin Methods with Reconstructions

In this section, we will consider using the new reconstruction method to post-process the solution obtained by using Fourier Galerkin methods to solve linear, hyperbolic time-dependent partial differential equations (PDEs) [16]. Consider the problem 

$$
\begin{cases}
u_t - \mathcal{L} u = 0, & x \in [-1, 1], \ t > 0 \\
u(x, 0) = g(x), & x \in [-1, 1]
\end{cases}
$$ 

(5.1)

where \( \mathcal{L} \) is a differential operator. In the Fourier Galerkin method, we seek solution \(u_N(x, t)\) from the space 

$$
\mathcal{B}_N = \text{span}\{e^{ikx}\}_{|k| \leq N}, \text{ i.e. }:
$$ 

$$
u_N(x, t) = \sum_{|k| \leq N} a_k(t) e^{ik \pi x}
$$ 

(5.2)

where \(a_k(t)\) are unknown coefficients which will be determined by the method. In general, the coefficients \(a_k(t)\) are not equal to the Fourier coefficients \(\hat{u}_k\); they will be equal only

if we obtain the exact solution of the problem. In the Fourier Galerkin method, \( a_k(t) \) are determined by requiring that the residual 

\[
R_N(x, t) = \frac{\partial u_N(x, t)}{\partial t} - \mathcal{L} u_N(x, t) \tag{5.3}
\]

is orthogonal to \( \hat{B}_N \). The method is defined by the requirement that the orthogonal projection of the residual onto the space \( \hat{B}_N \) is zero. If the residual is smooth enough, this requirement implies that the residual itself is small.

When spectral method is used to solve PDEs, we will estimate the errors in the maximum norm. If the solution lacks regularity,

\[
\|u(x, t) - u_N(x, t)\|_{L^\infty}
\]

cannot be small even for \( t = 0 \), due to the Gibbs phenomenon. We would like to reconstruct the solution from \( u_N(x, t) \) at the final time \( t = T \), and recover exponentially decaying error

\[
\|u(x, T) - Q_N u_N(x, T)\|_{L^\infty}
\]

where \( Q_N \) is a post-processing operator.

In order to illustrate the post-processing methods clearly, we consider the simple linear transport equations involving piecewise smooth functions with unbounded derivative singularities:

\[
\begin{cases}
u_t + c u_x = 0, & x \in [-1, 1], \, t > 0 \\
u(x, 0) = a(x) + b(x)(1 + x)^s, & x \in [-1, 1]
\end{cases} \tag{5.4}
\]

with periodic boundary conditions, where \( c \) is the phase speed (for simplicity, assume that \( c > 0 \)), \( s \) is a given fractional constant 

\[
0 < s = \frac{p}{q} < 1 \tag{5.5}
\]

with relatively prime integers \( p \) and \( q \), and \( a(x) \) and \( b(x) \) are both analytic functions. The initial condition is \( C^{0,\alpha} \) Hölder continuous with \( \alpha \leq s \) and 2-periodic. It has a singularity in \( [-1, 1] \), which means its first derivative blows up at \( x = -1 \). The solution to this problem is 

\[
u(x, t) = u(x - ct, 0) \tag{5.6}
\]

The singularities will move along the lines through points \( (x, t) = (2n - 1, 0) \) (here, \( n \in \mathbb{Z} \)) with direction \( (c, 1) \). Even though we consider functions with one singularity in this paper, our technique can be easily applied to functions with finitely many such singularities in \( [-1, 1] \) of known locations.

The coefficients \( a_k(t) \) are determined by requiring that the residual 

\[
R_N(x, t) = \frac{\partial u_N(x, t)}{\partial t} + c \frac{\partial u_N(x, t)}{\partial x} \tag{5.7}
\]

is orthogonal to \( \hat{B}_N \). If we rewrite the residual in terms of the Fourier series,

\[
R_N(x, t) = \sum_{|k| \leq \infty} \hat{R}_k(t)e^{ikx}, \tag{5.8}
\]

the orthogonality requirement yields 

\[
\hat{R}_k(t) = \frac{1}{2} \int_{-1}^{1} R_N(x, t)e^{-ik\pi x} dx = 0, \quad -N \leq k \leq N. \tag{5.9}
\]

1158
J Sci Comput (2015) 65:1145–1165

This requirement provides \( (2N + 1) \) ordinary differential equations to determine the \( (2N + 1) \) unknowns \( a_k(t) \),

\[
\frac{d}{dt} a_k(t) + i k \pi c a_k(t) = 0, \quad -N \leq k \leq N, \tag{5.10}
\]

and the corresponding initial conditions are

\[
u_N(x, 0) = \sum_{|k| \leq N} a_k(0) e^{i k \pi x}, \tag{5.11}
\]

\[
a_k(0) = \frac{1}{2} \int_{-1}^{1} u(x, 0) e^{-i k \pi x} dx \tag{5.12}
\]

For this linear constant coefficient problem, \( R_N(x, t) \) itself is in the space \( \mathcal{B}_N \), therefore the orthogonal complement must be zero, i.e.,

\[
R_N(x, t) = 0.
\]

It means that we obtain the exact solution of this linear problem by Fourier Galerkin methods. Therefore, the coefficients \( a_k(t) \) are exactly the Fourier coefficients for the solution \( u(x, t) \), which has singularities at \( (2n - 1 + ct, t) \) (here, \( n \in \mathbb{Z} \)).

In order to recover the exponential accuracy of the solution, we need to take post-processing on its translation \( v(x, t) = u(x + ct, t) \), which has a singularity at \( x = -1 \) in the interval \( [-1, 1] \). The function \( v(x, t) \) behaves like \( (1 + x)^s \) near the singularity. Under this translation, the Fourier coefficients \( b_k(t) \) of \( v(x, t) \) are

\[
b_k(t) = e^{i k \pi c t} a_k(t). \tag{5.13}
\]

From the truncated Fourier series \( v_N(x, t) = \sum_{|k| \leq N} b_k(t) e^{i k \pi x} \), we get a reconstructed approximation

\[
v^j_N(x, t) = \sum_{l=0}^{m} \hat{v}^j_l C^j_l(y(x)), \quad x \in [-1, 1] \tag{5.14}
\]

with

\[
\hat{v}^j_l = \frac{1}{h^j_t} \int_{-1}^{1} (1 - y^2)^{\frac{l-1}{2}} (1 - y^2)^{\frac{1}{2}} v_N(y(x), t) C^j_l(y) dy. \tag{5.15}
\]

Then the numerical solution after processing \( v^m_N(x, t) \) is obtained by

\[
u^m_N(x, t) = Q_N u_N(x, t) = v^m_N(x - ct, t). \tag{5.16}
\]

From Theorem 4.3, we obtain the following theorem on the accuracy in spectral methods solving the problem (5.4).

**Theorem 5.1** Using Fourier Galerkin method to solve problem (5.4), the solution \( u_N(x, t) \) is defined by

\[
u_N(x, t) = \sum_{|k| \leq N} a_k(t) e^{i k \pi x},
\]

where the coefficients \( a_k(t) \) are determined by (5.10) and (5.12). After post-processing, we get

\[
Q_N u_N(x, t) = v^m_N(x - ct, t).
\]

where, \( v^{m,N}_\lambda(x, t) \) is defined in (5.14) and (5.15). The error
\[
\max_{-1 \leq s \leq 1} |u(x, t) - Q_N u_N(x, t)|
\]
is exponentially small, when the parameters \( \lambda \) and \( m \) chosen in the post-processing are proportional to \( N \) and satisfy (3.4).

**Remark 5.2** For linear problems (5.4) with variable coefficients \( c(x) \), the computed coefficients \( a_k(t) \) are not exact for \( t > 0 \). But they are still computed with exponential accuracy, even with nonsmooth \( L^2 \) initial data. Therefore, we can still observe improved accuracy after post-processing, and numerical tests are shown in Sect. 6.

# 6 Numerical Results

In this section, we give numerical examples to illustrate the results of the reconstruction technique. In order to testify that exponential accuracy could be recovered from the knowledge of the first \( 2N + 1 \) Fourier coefficients, we first test on the same two examples as in [7].

## Example 6.1
We take the first \( 2N + 1 \) Fourier coefficients of function
\[
u(x) = \cos(x) + \sin(x)\sqrt{1+x}, \quad x \in [-1, 1] \tag{6.1}
\]
and try to recover the pointwise values of the solution over \([-1, 1]\).

First, the parameters are chosen according to
\[
\lambda = \frac{1}{16} N, \quad m = \frac{3}{80} N. \tag{6.2}
\]

For the sake of easy computation, we choose \( \lambda \) and \( m \) to be the biggest integer satisfying (6.2).

In Fig. 1a, we show the errors in the logarithm scale, for \( N = 40, 80, 160, 320 \) and 640. In Table 1 "linear choice", we show the maximum errors of the reconstruction over \( x \in [-1, 1] \) for each \( N \) and the orders of convergence. We can clearly find out that the error is exponentially decaying.

The other way is to find the optimal choice of the parameters such that the reconstruction achieves the smallest error. To choose the optimal parameters, we compare the errors of all the possible pairs \( 0 \leq \lambda \leq \frac{N}{2} \) and \( 0 \leq m \leq \frac{N}{2} \). In Table 1 "optimal choice", we also show the \( L^\infty \) errors of the reconstruction with the optimal choice for each \( N \).

## Example 6.2
We take the first \( 2N + 1 \) Fourier coefficients of function
\[
u(x) = \cos(x) + \sin(x)(1 + x)^{\frac{1}{3}}, \quad x \in [-1, 1] \tag{6.3}
\]
and try to recover the pointwise values of the solution over \([-1, 1]\).

First, the parameters are chosen linearly dependent on \( N \):
\[
\lambda = \frac{1}{8} N, \quad m = \frac{1}{32} N. \tag{6.4}
\]

In Fig. 1b, we show the errors in the logarithm scale, for \( N = 40, 80, 160, 320 \) and 640. In Table 2, we show the maximum errors of the reconstruction over \( x \in [-1, 1] \) with both ways to choose parameters, for each \( N \). We again clearly find out that the error is exponentially decaying with linear choice.

\[
\text{Fig. 1 Pointwise errors in the logarithm scale.} \quad \text{a Example 1} \quad \left( \lambda = \frac{1}{8} N, m = \frac{1}{32} N \right).
\]
\[
\text{b Example 2} \quad \left( \lambda = \frac{1}{16} N, m = \frac{3}{80} N \right).
\]

Next, we reconstruct the same two functions with different levels of noise in data to test the robustness to noise. The way we add noise to Fourier coefficients is to generate two random and independent sequences, which are uniformly distributed in \(\left[-\frac{\delta}{2}, \frac{\delta}{2}\right]\) with noise level \(\delta\), and then add the sequences to the real parts and imaginary parts of the Fourier coefficients separately.

\textbf{Example 6.3 (Noise)} For the function 
\[
u(x) = \cos(x) + \sin(x) \sqrt{1 + x}, \quad x \in [-1, 1] \tag{6.5}
\]
we try to recover it with given information of the first \(2N + 1\) Fourier coefficients with a certain level of noise. We compare the maximum error of the reconstruction with different levels of noise, for either the linear choice or the optimal choice of parameters \(\lambda\) and \(m\) as in Example 6.1. The results are listed in Tables 3 and 4.

# J Sci Comput (2015) 65:1145–1165

## Table 1 Maximum error table

| N   | Linear choice  | \( L^\infty \) error | Order | \( \lambda \) | m | Optimal choice | \( L^\infty \) error | \( \lambda \) | m |
|-----|----------------|----------------------|-------|---------------|---|----------------|----------------------|---------------|---|
| 40  | 5.91E–001      | 2                    | 1     | 1             | 1 | 1.21E–003      | 8                    | 18            |
| 80  | 2.95E–001      | 1.00                 | 5     | 3             | 3 | 5.25E–005      | 12                   | 40            |
| 160 | 1.56E–002      | 4.25                 | 10    | 6             | 6 | 5.57E–010      | 20                   | 80            |
| 320 | 1.33E–005      | 10.19                | 20    | 12            | 12| 2.75E–011      | 22                   | 51            |
| 640 | 6.46E–012      | 20.98                | 40    | 24            | 24| 1.11E–012      | 23                   | 44            |

## Table 2 Maximum error table

| N   | Linear choice  | \( L^\infty \) error | Order | \( \lambda \) | m | Optimal choice | \( L^\infty \) error | \( \lambda \) | m |
|-----|----------------|----------------------|-------|---------------|---|----------------|----------------------|---------------|---|
| 40  | 7.13E–001      | 5                    | 1     | 1             | 5 | 1.83E–001      | 5                    | 12            |
| 80  | 3.64E–001      | 0.97                 | 10    | 2             | 2 | 1.98E–002      | 7                    | 33            |
| 160 | 2.91E–001      | 0.33                 | 20    | 5             | 5 | 1.24E–003      | 12                   | 76            |
| 320 | 1.17E–002      | 4.64                 | 40    | 10            | 10| 5.50E–006      | 17                   | 97            |
| 640 | 1.83E–006      | 12.64                | 80    | 20            | 20| 3.15E–007      | 21                   | 89            |

## Table 3 Maximum error with different levels of noise in data (linear choice)

| N   | \( L^\infty \) (no noise) | \( L^\infty \) Noise | \( L^\infty \) Noise | \( L^\infty \) Noise |
|-----|----------------------------|----------------------|----------------------|----------------------|
|     | \( L^\infty \) Noise      | \( L^\infty \) Noise | \( L^\infty \) Noise |                      |
| 40  | 5.91E–01                  | 6.07E–01            | 1.00E–02            | 5.53E–01            |
| 80  | 2.95E–01                  | 1.82E–01            | 1.00E–02            | 3.25E–01            |
| 160 | 1.56E–02                  | 1.53E–02            | 1.00E–05            | 8.16E–03            |
| 320 | 1.33E–05                  | 1.36E–05            | 1.00E–11            | 1.12E–05            |

## Table 4 Maximum error with different levels of noise in data (optimal choice)

| N   | \( L^\infty \) (no noise) | \( L^\infty \) Noise | \( L^\infty \) Noise | \( L^\infty \) Noise |
|-----|----------------------------|----------------------|----------------------|----------------------|
|     | \( L^\infty \) Noise      | \( L^\infty \) Noise | \( L^\infty \) Noise |                      |
| 40  | 1.21E–03                  | 1.16E–03            | 1.00E–07            | 3.30E–03            |
| 80  | 5.25E–05                  | 4.74E–05            | 1.00E–12            | 1.20E–04            |
| 160 | 5.57E–11                  | 4.97E–10            | 1.00E–22            | 4.81E–09            |
| 320 | 2.75E–11                  | 2.38E–11            | 1.00E–23            | 3.25E–11            |

## Example 6.4 (Noise) For the function

\[
u(x) = \cos(x) + \sin(x)(1 + x)^{\frac{1}{3}}, \quad x \in [-1, 1]
\]

we try to recover it with given information of the first \( 2N + 1 \) Fourier coefficients with certain level of noise. We compare the maximum error of the reconstruction with different levels of noise, for either the linear choice or the optimal choice of parameters \( \lambda \) and \( m \) as in Example 6.2. The results are listed in Tables 5 and 6.

Table 5 Maximum error with different levels of noise in data (linear choice)

| N   | L∞ (no noise) | L∞ Noise | L∞ (no noise) | L∞ Noise |
|-----|----------------|----------|----------------|----------|
| 40  | 7.13E-01       | 6.82E-01 | 1.00E-02       | 7.21E-01 |
| 80  | 3.64E-01       | 3.43E-01 | 1.00E-02       | 8.04E-01 |
| 160 | 2.91E-01       | 2.81E-01 | 1.00E-04       | 8.87E-01 |
| 320 | 1.17E-02       | 1.65E-02 | 1.00E-08       | 4.93E-02 |

Table 6 Maximum error with different levels of noise in data (optimal choice)

| N   | L∞ (no noise) | L∞ Noise | L∞ (no noise) | L∞ Noise |
|-----|----------------|----------|----------------|----------|
| 40  | 1.83E-01       | 1.80E-01 | 1.00E-04       | 1.77E-01 |
| 80  | 1.98E-02       | 1.01E-02 | 1.00E-06       | 8.46E-02 |
| 160 | 1.24E-03       | 1.21E-03 | 1.00E-12       | 7.29E-03 |
| 320 | 5.50E-06       | 5.16E-06 | 1.00E-18       | 2.73E-06 |

From the results in Example 6.3 and Example 6.4, we observe that the reconstruction methods work pretty well with noise. The more accurate the reconstruction is, the more sensitive to noise it is. The linear choice of parameters is less sensitive, since it is patterned. The optimal one is numerically chosen for best errors, thus is more sensitive to any perturbation. And the allowable noise level to actual error level is about \( \text{error}^2 \) for Example 6.3 and \( \text{error}^3 \) for Example 6.4.

At last, we have numerical tests on Fourier Galerkin methods with post-processing for initial boundary problems of linear transport equation with variable coefficients.

Example 6.5 We solve the following initial boundary problem
\[
\begin{cases}
u_t - xu_{xx} = 0, & x \in [-1, 1], \ t > 0 \\
u(x, 0) = g(x) = \sqrt{1+x}, & x \in [-1, 1]
\end{cases}
\]
with periodic boundary condition.

The exact solution is \( u(x, t) = \sqrt{1 + x e^t} \mod 2 \). If we look at the solution at \( T = \log 2 \), it has two deformed “copies” of the initial condition, and has singularities at two places -0.5 and 0.5. Therefore, the filtered Fourier Galerkin method with Runge-Kutta time stepping (as in [18]) provides solutions with poor accuracy and oscillations near the singularities, as in Fig. 2. Here, we would like to recover the solution on \( [-0.5, 0] \), since the solution behaves like square root at the left end of this sub-interval.

\[
u(x, \log 2) = \sqrt{1 + 2x}, \quad x \in [-0.5, 0].
\]

We choose the parameters \( \lambda \) and \( m \) in two different ways. First, we choose \( \lambda \) and \( m \) linearly dependent on \( N \):
\[
\lambda = \frac{N}{80}, \quad m = \frac{N}{4}.
\]

The other way is to choose the pair, which gives best recovery, in the range \( 0 \leq \lambda, m \leq \frac{N}{2} \). The results are listed in Table 7 and shown in Fig. 3. With both choices, the solutions after...

J Sci Comput (2015) 65:1145–1165

**Fig. 2** Numerical solution \( u_N(x) \) (N = 10, 20, 40, and 80) on \([-0.5, 0]\), using Fourier Galerkin method with filters. 
a Numerical solutions \( u_N(x) \) and exact solution \( u(x) \). 
b Pointwise errors in logarithm scale

post-processing provide improved accuracy, even up to the singularities. We also observe exponential accuracy from the reconstruction with the linear choice.

### 7 Concluding Remarks

We have built the Gegenbauer polynomial based technique to reconstruct approximations with exponential accuracy in the maximum norm, based on the Fourier coefficients, for piecewise analytic functions where in each piece the function is analytic only in the open interval with unbounded derivative end-point singularities. Such functions are significantly more difficult to handle than piecewise analytic functions which are analytic up to the discontinuity points. Our technique provides a post-processing method for spectral methods to solve transport

| N  | Linear choice         | Optimal choice       | Without PP           |
|----|-----------------------|----------------------|----------------------|
|    | L∞ error | Order | λ | m | L∞ error | λ | m | L∞ error |
|----|-----------|-------|---|---|-----------|---|---|-----------|
| 10 | 6.93E-01 | 0     | 3 | 5.81E-01 | 1 | 2 | 6.57E-01 |
| 20 | 5.13E-01 | 0.43  | 5 | 1.72E-01 | 1 | 8 | 5.04E-01 |
| 40 | 1.01E-02 | 5.67  | 1 | 10 | 7.34E-03 | 1 | 6 | 4.45E-01 |
| 80 | 2.61E-05 | 8.60  | 1 | 20 | 2.61E-05 | 1 | 20 | 4.07E-01 |

### Fig. 3
Pointwise errors of the reconstructions with two different strategies for choosing parameters \( \lambda \) and \( m \), in logarithm scale. 
- a) Linear choices. 
- b) Optimal pairs.

equations involving such functions to achieve exponential accuracy. Preliminary numerical results are provided to demonstrate the theory. This technique is easy to be generated to higher dimensional cases, and has potential applications for solving other partial differential equations whose solutions have such singularities.

## References

1. Adcock, B., Richardson, M.: New exponential variable transform methods for functions with endpoint singularities. SIAM J. Numer. Anal. 52, 1887–1912 (2014)
2. Archibald, R., Chen, K., Gelb, A., Renaut, R.: Improving tissue segmentation of human brain MRI through preprocessing by the Gegenbauer reconstruction method. NeuroImage 20, 489–502 (2003)
3. Archibald, R., Gelb, A.: A method to reduce the Gibbs ringing artifact in MRI scans while keeping tissue boundary integrity. IEEE Med. Imaging 21, 305–319 (2002)
4. Archibald, R., Gelb, A.: Reducing the effects of noise in image reconstruction. J. Sci. Comput. 17, 167–180 (2002)
5. Archibald, R., Hu, J., Gelb, A., Farin, G.: Improving the accuracy of volumetric segmentation using pre-processing boundary detection and image reconstruction. IEEE Trans. Med. Imaging 13, 459–466 (2004)
6. Bateman, H.: Higher Transcendental Functions, v2. McGraw-Hill, New York (1953)
7. Chen, Z., Shu, C.-W.: Recovering exponential accuracy from collocation point values of smooth functions with end-point singularities. J. Comput. Appl. Math. 265, 83–95 (2014)
8. Gottlieb, D., Gottlieb, S.: Spectral methods for compressible reactive flows. C. R. Mec. 333, 3–16 (2005)
9. Gottlieb, D., Shu, C.-W.: Resolution properties of the Fourier method for discontinuous waves. Comput. Methods Appl. Mech. Eng. 116, 27–37 (1994)
10. Gottlieb, D., Shu, C.-W.: On the Gibbs phenomenon IV: recovering exponential accuracy in a sub-interval from a Gegenbauer partial sum of a piecewise analytic function. Math. Comput. 64, 1081–1095 (1995)
11. Gottlieb, D., Shu, C.-W.: On the Gibbs phenomenon V: recovering exponential accuracy from collocation point values of a piecewise analytic function. Numer. Math. 71, 511–526 (1995)
12. Gottlieb, D., Shu, C.-W.: On the Gibbs phenomenon III: recovering exponential accuracy in a sub-interval from a spectral partial sum of a piecewise analytic function. SIAM J. Numer. Anal. 33, 280–290 (1996)
13. Gottlieb, D., Shu, C.-W.: On the Gibbs phenomenon and its resolution. SIAM Rev. 30, 644–668 (1997)
14. Gottlieb, D., Shu, C.-W., Solomonoff, A., Vandeven, H.: On the Gibbs phenomenon I: recovering exponential accuracy from the Fourier partial sum of a non-periodic analytic function. J. Comput. Appl. Math. 43, 81–98 (1992)
15. Gottlieb, S., Gottlieb, D., Shu, C.-W.: Recovering high order accuracy in WENO computations of steady state hyperbolic systems. J. Sci. Comput. 28, 307–318 (2006)
16. Hesthaven, J.S., Gottlieb, S., Gottlieb, D.: Spectral Methods for Time-Dependent Problems. Cambridge Monographs on Applied and Computational Mathematics, Vol. 21. Cambridge University Press, Cambridge (2007)
17. Jung, J.-H., Gottlieb, S., Kim, S.O., Bresten, C.L., Higgs, D.: Recovery of high order accuracy in radial basis function approximations of discontinuous problems. J. Sci. Comput. 45, 359–381 (2010)
18. Shu, C.-W., Wong, P.S.: A note on the accuracy of spectral method applied to nonlinear conservation laws. J. Sci. Comput. 10, 357–369 (1995)