# Recovering exponential accuracy from collocation point values of smooth functions with end-point singularities

**Zheng Chen, Chi-Wang Shu**  
Division of Applied Mathematics, Brown University, Providence, RI 02912, United States

---

## A R T I C L E  I N F O

**Article history:**  
Received 18 November 2012  
Received in revised form 5 April 2013  

Dedicated to Professor Ben-yu Guo on his seventieth birthday

**Keywords:**  
Spectral approximation  
Exponential accuracy  
Gegenbauer expansion  
Collocation  
Gaussian points  

---

## A B S T R A C T

Gibbs phenomenon is the particular manner how a global spectral approximation of a piecewise analytic function behaves at the jump discontinuity. The truncated spectral series has large oscillations near the jump, and the overshoot does not decay as the number of terms in the truncated series increases. There is therefore no convergence in the maximum norm, and convergence in smooth regions away from the discontinuity is also slow. In Gottlieb and Shu (1995) [5], a methodology is proposed to completely overcome this difficulty in the context of spectral collocation methods, resulting in the recovery of exponential accuracy from collocation point values of a piecewise analytic function. In this paper, we extend this methodology to handle spectral collocation methods for functions which are analytic in the open interval but have singularities at end-points. With this extension, we are able to obtain exponential accuracy from collocation point values of such functions. Similar to Gottlieb and Shu (1995) [5], the proof is constructive and uses the Gegenbauer polynomials \( C_n^{\lambda}(x) \). The result implies that the Gibbs phenomenon can be overcome for smooth functions with endpoint singularities.

© 2013 Elsevier B.V. All rights reserved.

---

## 1. Introduction

In this paper, we investigate the issue of overcoming the Gibbs phenomenon. Gibbs phenomenon is the particular manner how a global spectral approximation of a piecewise analytic function behaves at the jump discontinuity. A prototype is the Fourier approximation of an analytic but non-periodic function, which, when viewed as a periodic function, has a discontinuity at the boundaries of the interval. The truncated Fourier series has large oscillations near this jump, and the overshoot does not decay as the number of terms retained in the series increases. There is therefore no convergence in the maximum norm. Also, convergence in smooth regions away from the discontinuity is slow. In a series of papers [1–5], Gottlieb et al. developed a general framework to overcome this difficulty, in the sense that exponential accuracy is recovered in the maximum norm for any sub-interval (including the whole interval if there is only one discontinuity which is located at the interval boundaries), from the knowledge of either the first \( N \) spectral expansion coefficients, or the point values at \( N \) standard collocation points. This means that exponential accuracy is recovered at all points, including at the actual discontinuity points (the left and right limits at these points), if the locations of these discontinuity points are not known exactly but are known to be within certain fixed intervals, then exponential accuracy can be recovered from any interval which does not overlap with these fixed intervals containing the discontinuities. An important tool used in this framework is the set of Gegenbauer polynomials, which are orthogonal on the interval \([-1, 1]\ with the weight \( (1 - x^2)^{\lambda - \frac{1}{2}} \). It turns out that the main ingredient in this technique is that the parameter

\[
\lambda \text{ in the weight function as well as the number of terms m retained in the Gegenbauer expansion should both be chosen proportional to N. For a review of this series of work, we refer to [6].}
\]

In this study, we are interested in collocation methods. In [5], a methodology is proposed to completely overcome the difficulty of Gibbs phenomenon, in the context of spectral collocation methods, resulting in the recovery of exponential accuracy from collocation point values of a piecewise analytic function. In this paper, we extend this methodology to handle spectral collocation methods for functions which are analytic in the open interval but have singularities at end-points. Such functions appear often in applications, for example, as solutions of certain partial differential equations (PDEs). We assume that we are given the point values \( f(x_i) \), where \( x_i \) are the standard collocation points (Gaussian points of the orthogonal basis \( \{p_k(x)\} \), which may be the Fourier trigonometric polynomials, the Legendre polynomials, the Chebyshev polynomials, or the general Gegenbauer polynomials). The objective is to recover exponentially accurate point values at every point including at the singularities.

The reconstruction procedure is performed on functions of the following form

\[
f(x) = a(x) + b(x)(1 + x)^s, \quad x \in [-1, 1] \tag{1.1}
\]

where \( s \) is a given fractional constant

\[
0 < s < \frac{p}{q} < 1 \tag{1.2}
\]

in which \( p \) and \( q \) are relatively prime positive integers, and \( a(x) \) and \( b(x) \) are both analytic functions. Unfortunately, our current approach cannot treat functions of the form (1.1) with an irrational power \( s \).

Since we assume that the value of \( s \) is known, we are not losing generality by assuming (1.1) with \( 0 < s < 1 \). All functions same from with \( 0 < s < 1 \) multiplying with an analytic function \( (1 + x)^k \) with a positive integer \( k \), or by absorbing an analytic function \( (1 + x)^k \) with a positive integer \( k \) into \( b(x) \). Since we are handling collocation methods and are given values of the original functions at the collocation points, we also have access to the collocation point values of the modified functions which are the original functions multiplied with \( (1 + x)^k \). Singularity at the right end or at both ends can be handled in the same fashion. Also, the result can be easily generalized to the situation of finitely many singularities (of the form \( (x - z_j)^k \) at finitely many points \( z_j \) inside the interval \( [-1, 1] \) with fraction \( s_k \), along the lines of [3,4] and using the techniques in this paper.

As in [1], we assume that the analytic functions \( a(x) \) and \( b(x) \), denoted generically as \( c(x) \), satisfy the following condition.

**Assumption 1.1.** There exists a constant \( \rho \geq 1 \) and a constant \( C(\rho) \) such that, for every \( k \geq 0 \),

\[
\max_{-1 \leq x \leq 1} \left| \frac{d^k c(x)}{dx^k} \right| \leq C(\rho) \frac{k!}{\rho^k}.
\]

This is a standard assumption for analytic functions, where \( \rho \) is the distance from the interval \( [-1, 1] \) to the nearest singularity of the function \( c(x) \) in the complex plane.

We will use the following one to transformation between \( x \in [-1, 1] \) and \( y \in [-1, 1] \):

\[
(2^{q-1}(1 + x))^j = 1 + y
\]

where \( q \) is defined in (1.2) in which \( p \) and \( q \) are assumed to be relatively prime.

The function \( f(y) = f(x(y)) \) of the variable \( y \) has its usual Gegenbauer expansion under the basis \( \{C_n^q(y)\} \):

\[
f(x(y)) = \bar{f}(y) = \sum_{l=0}^{\infty} \hat{f}^l (C_l^q(y))
\]

with the Gegenbauer coefficients \( \hat{f}^l \) given by

\[
\hat{f}^l = \frac{1}{h_l^q} \int_{-1}^{1} (1 - y^2)^{l-1} f(y) C_l^q(y) dy \tag{1.3}
\]

where the precise value of the normalization constant \( h_l^q \) will be given by (2.3) in Definition 2.1.

Our goal is to find a good approximation to the first \( m \sim N \) Gegenbauer coefficients \( \hat{f}^l \) in (1.3), denoted as \( g^l \), from the known point values \( f(x_i) \) at the standard Gaussian collocation points. We will then obtain the approximation of \( f(x) \) using these \( m \sim N \) terms of its Gegenbauer expansion:

\[
f_m^N(x) = \sum_{l=0}^{m} g^l (C_l^q(y)(x)).
\]

As in [5], we will separate the analysis of the error into two parts: the truncation error and the regularization error. Unlike [5], the analysis for the truncation error must be completely changed because of the end-point singularity. The truncation error measures the difference between the exact Gegenbauer coefficients of \( f(x) \) with \( \sim N \), and the approximate

Z. Chen, C.-W. Shu / Journal of Computational and Applied Mathematics 265 (2014) 83–95

Gegenbauer coefficients \( g^k(1) \) obtained by using the collocation point values. This will be investigated in Section 3. The regularization error measures the difference between the Gegenbauer expansion with \( \lambda \sim N \), using the first \( m \sim N \) Gegenbauer coefficients, and the function itself. This error is estimated in Section 4. The results are summarized in Theorem 4.3 in Section 4. Section 5 contains several numerical examples to illustrate our results. In Section 2, we shall give several useful preliminary properties, estimates, and approximation results. Concluding remarks are given in Section 6.

## 2. Preliminaries

In this section, we will first introduce the Gegenbauer polynomials and discuss some of their asymptotic behavior. We rely heavily on the standardization in Bateman [7].

### Definition 2.1
The Gegenbauer polynomial \( C_n^\lambda(x) \), for \( \lambda \geq 0 \), is defined by

\[
(1 - x^2)^{\frac{1}{2} - \frac{1}{2} n} C_n^\lambda(x) = \frac{(-1)^n}{2^n n!} G(\lambda, n) \frac{d^n}{dx^n} \left[ (1 - x^2)^{n + \frac{\lambda}{2} - \frac{1}{2}} \right]
\]

where \( G(\lambda, n) \) is given by

\[
G(\lambda, n) = \frac{\Gamma\left(\lambda + \frac{1}{2}\right) \Gamma(\lambda + 2)}{\Gamma(2\lambda) \Gamma\left(n + \frac{1}{2}\right)}
\]

for \( \lambda > 0 \), by

\[
G(0, n) = \frac{2}{\sqrt{\pi} (n - 1)!} \quad \text{for } n = 0 \text{ and } n \geq 1,
\]

and by \( G(0, 0) = 1 \).

For \( \lambda = 0 \) and \( n = 0 \), notice that this standardization, \( C_n^0(x) \) is defined by (see [7]):

\[
C_0^0(x) = \lim_{n \to 0^+} \frac{1}{n} C_n^0(x) = T_n(x), \quad n > 0; \quad C_0^0(x) = 1,
\]

where \( T_n(x) \) are the Chebyshev polynomials.

Under this definition we have, for \( \lambda > 0 \),

\[
C_n^\lambda(1) = \frac{\Gamma(n + 2\lambda)}{n! (2\lambda)}; \quad \text{for } \lambda = 0 \text{ and } n = 2,
\]

\[
C_0^\lambda(1) = \frac{2}{n}; \quad \text{for } \lambda = 0 \text{ and } n = 0,
\]

and

\[
C_0^\lambda(1) = 1;
\]

and 

\[
|C_n^\lambda(x)| \leq C_n^\lambda(1), \quad -1 \leq x \leq 1. \tag{2.2}
\]

The Gegenbauer polynomials are orthogonal under their weight function \( (1 - x^2)^{\lambda - \frac{1}{2}} \):

\[
\int_{-1}^{1} (1 - x^2)^{\frac{1}{2} - \frac{1}{2} n} C_n^\lambda(x) C_n^\lambda(x) dx = \delta_{k,n} h_n^\lambda
\]

where, for \( \lambda > 0 \),

\[
h_n^\lambda = \pi \frac{C_n^\lambda(1)}{\Gamma\left(\lambda + \frac{1}{2}\right) \Gamma(\lambda + n)}; \quad \text{for } \lambda = 0 \text{ and } n \geq 1,
\]

\[
h_0^\lambda = \frac{2\pi}{n^2}; \quad \text{for } \lambda = 0 \text{ and } n = 0,
\]

\[
h_0^0 = \pi.
\]

We will need to use the Stirling's formula and the estimate of \( h_k^\lambda \) for the asymptotics of the Gegenbauer polynomials for large \( n \) and \( \lambda \).

## Lemma 2.2
We have the Stirling's formula
\[
(2\pi)^{\frac{1}{2}} x^{x+\frac{1}{2}} e^{-x} \leq \Gamma(x + 1) \leq (2\pi)^{\frac{1}{2}} x^{x+\frac{1}{2}} e^{-x}, \quad x \geq 1. \tag{2.4}
\]

## Lemma 2.3
There exists a constant \( C \) independent of \( \lambda \) and \( n \) such that
\[
C^{-1} \frac{\lambda}{n + \lambda} C_k^\lambda(1) \leq h_k^\lambda \leq C \frac{\lambda}{n + \lambda} C_k^\lambda(1). \tag{2.5}
\]

We also need to quote the approximation results for the interpolation polynomials. Given the point values \( f(x_i) \) of the function \( f(x) \) at the \( N \) Gaussian points \( \{x_i\} \) of the basis functions \( \{\phi_k(x)\} \), we denote the unique interpolation polynomial by
\[
I_N(f)(x) = \sum_{k=1}^{N} \hat{f}_k \phi_k(x) \tag{2.6}
\]
satisfying
\[
I_N(f)(x_i) = f(x_i), \quad i = 1, 2, \ldots, N. \tag{2.7}
\]

For the interpolation polynomial \( I_N(f)(x) \), we have the following error estimates:

## Lemma 2.4
If \( f(x) \) has continuous derivatives in \( [-1, 1] \), then we have the following estimate for the interpolation polynomial \( I_N(f) \) defined by (2.6) and (2.7), where the basis functions \( \{\phi_k(x)\} \) are either trigonometric polynomials \( e^{ikx} \) (with \( k = -N, \ldots, N \)) or the Gegenbauer polynomials \( C_k^\lambda(x) \) with \( \mu > -\frac{1}{2} \) (including the Chebyshev case for \( \mu = 0 \) and the Legendre case for \( \mu = \frac{1}{2} \)):
\[
\| f - I_N(f) \|_2 \leq \frac{C}{N^m} \| f \|_{m, \infty}, 
\]
where the weighted \( L^2 \) norm is defined as
\[
\| f \|_{L^2}^2 = \int_{-1}^{1} \omega(x) |f(x)|^2 \, dx
\]
with the weight function \( \omega(x) = 1 \) for the trigonometric or Legendre polynomials bases and \( \omega(x) = (1 - x^2)^{\mu - \frac{1}{2}} \) for the general Gegenbauer polynomial bases. Here \( C \) is a constant independent of \( N \) and \( m \).

## Proof
The proof for the trigonometric polynomial case can be found in [8]. For the Gegenbauer polynomial case, the proof can be found in [9]. Although the authors do not explicitly point out, the constant \( C \) is independent of \( m \) in both situations.

## Remark 2.5
\[
(1 - y(x))^{\lambda - \frac{1}{2}} \left( (1 - y(x))^{\frac{1}{2}} (1 + y) \right)^{-\frac{2}{3}} f(x) 
\]
has up to \( t \)-th derivatives in \( x \), where \( t = \left\lfloor \frac{a}{2} \right\rfloor - 1 \leq a_2 - 1 \), the largest integer below \( \frac{a}{2} - 1 \).

It is easy to observe that
\[
\frac{d^m}{dx^m} (1 - y(x)^2)^{\lambda - \frac{1}{2}} = A y^{\gamma_1} y^{\gamma_2} y^{\gamma_3}, \quad 0 \leq n \leq t
\]
where
\[
A = \frac{2^n}{2q}, \quad y^{\gamma_1} = (1 - y(x))^{\lambda - \gamma_n}, \quad y^{\gamma_2} = (1 - y(x))^{\gamma - 1},
\]
and \( y^{\gamma_3} \) satisfies the following recursive relation:
\[
y_0^{\gamma} = 1
\]
\[
y_{n+1}^{\gamma} = -[2A - qn - n(1 - y)(n - qn - 1)] y_n^{\gamma}, \quad 0 \leq n < t.
\]

It is easy to show that \( y^{\gamma_3} \) is an \( n \)-th degree polynomial of \( y \). We have the following estimate on \( y^{\gamma_3} \).

# Z. Chen, C.-W. Shu / Journal of Computational and Applied Mathematics 265 (2014) 83–95

## Lemma 2.6
We have, for \( 0 \leq n \leq t \),
\[
|y_n^3| \leq (2\lambda)^n, \quad y \in [-1, 1].
\]

**Proof.** Let us rewrite \( y_n^3 = \sum_{i=0}^{n} a_i y^i \) and denote \( S_n = \sum_{i=0}^{n} |a_i| \).
If \( |S_n| \leq 2^{-n} \), then \( |y| \leq \sum_{i=0}^{n} |a_i| \leq (2\lambda)^n \).
Therefore, to prove \( |y_n^3| \leq (2\lambda)^n \), we only need to prove \( |S_n| \leq (2\lambda)^n \), for which we will use induction.
1. For \( n = 0 \): \( S_0 = 1 \leq (2\lambda)^0 \).
2. Suppose \( S_n \leq (2\lambda)^n \) for \( n \).

Then for \( n + 1 \),
\[
|S_{n+1}| \leq |[(2\lambda - qn - n - 1)y + n(q - 1)]S_n + 2n S_n| \leq |[(2\lambda - qn - n + n(q - 1) + 2n](2\lambda)^n|.
\]

Thus, we have proved \( |y_n^3| \leq (2\lambda)^n \).

## Lemma 2.7
We have, for \( 0 \leq n \leq t \),
\[
\frac{d^m}{dx^m}[(1 - y^2)^{-\frac{1}{2}}] \leq (2A)^{\lambda}(1 - y)^{n(q-1)}(1 - y^2)^{-\frac{1}{2} - \alpha_n}.
\]

**Proof.** We only need to use (2.8).
Now let us rewrite the function
\[
(1 - y^2)^{-\frac{1}{2}}(x)(1 + y) = C(1 - y^2)^{-\frac{1}{2}} \{ a(x)(1 + x)^s_1 + b(x)(1 + x)^s_2 \}
\]
where
\[
s_1 = \left( -\frac{q - 1}{2q}, \left( -\frac{1}{2}, 0 \right) \right), \quad s_2 = \left( \frac{p - q - 1}{2q}, \left( -\frac{1}{2}, 1 \right) \right).
\]

## Lemma 2.8
We have, for \( 0 \leq n \leq t \) and \( j = 1, 2 \),
\[
\frac{d^m}{dx^m}[(1 + x)^j(1 - y^2)^{-\frac{1}{2}}] \leq C(2A)^{\lambda}(1 + y)^{j}(1 - y)^{n(q-1)}(1 - y^2)^{-\frac{1}{2} - \alpha_n} \left( 1 + \frac{1}{e} \right)^n,
\]

where \( C \) is at most a constant multiplied by \( \sqrt{n} \).

**Proof.** We have
\[
\frac{d^m}{dx^m}[(1 + x)^j(1 - y^2)^{-\frac{1}{2}}] = \sum_{i=0}^{n} \binom{n}{i} \frac{d^i}{dx^i}(1 + x)^j \frac{d^{n-i}}{dx^{n-i}}(1 - y^2)^{-\frac{1}{2}}.
\]
\[
= \sum_{i=1}^{n} \binom{n}{i} s_j \cdots (s_j + 1 - i)(1 + x)^s_j \frac{d^{n-i}}{dx^{n-i}}(1 - y^2)^{-\frac{1}{2}} + |(1 + x)^j| \frac{d^m}{dx^m}(1 - y^2)^{-\frac{1}{2}}.
\]
\[
\leq \sum_{i=1}^{n} \binom{n}{i} iA(q - s)^{(s)}(1 + y)^{s}(2A)^{\lambda}(1 + y)^{s}(1 - y)^{n(q-1)} \left[ \sum_{j=1}^{n} \frac{(1 - y^q)}{2\lambda} \right].
\]
\[
= C(2A)^{\lambda}(1 - y)^{n(q-1)}(1 + y)^{j}(1 - y^2)^{-\frac{1}{2} - \alpha_n} \leq C(2A)^{\lambda}(1 - y)^{n(q-1)}(1 + y)^{j}\sum_{i=0}^{n} \binom{n}{i} \left( \frac{(1 - y)q}{2e} \right)^i.
\]

88
Z. Chen, C.-W. Shu / Journal of Computational and Applied Mathematics 265 (2014) 83–95

\[
\leq C(2\lambda)^{n}(1 - y)^{\gamma - 1}(1 + y)^{\delta}(1 - y)^{2 - \eta n} \sum_{i=0}^{n} \binom{n}{i} \left( \frac{1}{e} \right)^{i}
\]

\[
= C(2\lambda)^{n}(1 - y)^{\gamma - 1}(1 + y)^{\delta}(1 - y)^{2 - \eta n} \left( 1 + \frac{1}{e} \right)^{n}
\]

where we have used Lemma 2.7 in the third inequality, (2.4) for the fifth inequality, and the binomial formula for the last equality.

**Lemma 2.9.** We have, for \( 0 \leq s \leq t \) and \( j = 1, 2 \)

\[
\frac{d}{dx} \left[ (1 - y)^{2 - \frac{1}{2} (1 + x)} \alpha(x) \right] \leq C(0)(1 + y)^{\delta}(1 - y)^{2 - \eta}(2\lambda)^{n} \left[ 1 + \frac{1}{e} \left( \frac{1}{2\lambda q} \right) \right]
\]

where \( \alpha(x) \) is an analytic function satisfying Assumption 1.1, and \( s_j \in \left( -\frac{1}{2}, 1 \right) \) are defined in (2.12). \( C(0) \) is at most a constant multiplied by \( n \).

**Proof.** We have

\[
\frac{d}{dx} \left[ (1 - y)^{2 - \frac{1}{2} (1 + x)} \alpha(x) \right] = \sum_{i=0}^{n} \binom{n}{i} \frac{d}{dx} \left[ (1 - y)^{2 - \frac{1}{2} (1 + x)} \right] \frac{d^{n-i}}{dx^{n-i}} \alpha(x)
\]

\[
\leq \sum_{i=0}^{n} (2\lambda)^{i}(1 + y)^{\delta}(1 - y)^{\gamma - 1}(1 - y)^{2 - \eta n} \left( 1 + \frac{1}{e} \right)^{i} C(0)(n - i)! \frac{d^{n-i}}{dx^{n-i}} \alpha(x)
\]

\[
\leq C(0)(1 + y)^{\delta}(1 - y)^{2 - \eta} \sum_{i=0}^{n} \binom{n}{i} 2\lambda(1 - y)^{\gamma - 1}(1 - y)^{2 - \eta} \left( 1 + \frac{1}{e} \right)^{i} \left( \frac{n}{e \rho^{i}} \right)^{n-i}
\]

\[
= C(0)(1 + y)^{\delta}(1 - y)^{2 - \eta}(2\lambda)^{n} \left[ (1 - y)^{\gamma - 1} \left( 1 + \frac{1}{e} \right) \right]^{n}
\]

\[
\leq C(0)(1 + y)^{\delta}(1 - y)^{2 - \eta}(2\lambda)^{n} \left[ 2^{-1} \left( 1 + \frac{1}{e} \right) + \frac{n}{2\lambda e} \right]^{n}
\]

\[
= C(0)(1 + y)^{\delta}(1 - y)^{2 - \eta}(2\lambda)^{n} \left[ 1 + \frac{1}{e} + \frac{n}{2\lambda e} \right]^{n}
\]

\[
\leq C(0)(1 + y)^{\delta}(1 - y)^{2 - \eta}(2\lambda)^{n} \left[ 1 + \frac{1}{e} \left( \frac{1}{2\lambda q} \right) \right]^{n}
\]

where we have used Lemma 2.8 and Assumption 1.1 in the second inequality, (2.4) for the third inequality, and the binomial formula for the fourth equality.

We have the following properties for the Gegenbauer expansion of analytic functions, which will be used in Section 4.

**Lemma 2.10 (The Exponential Decay of the Regularization Error).** Let \( f(x) \) be an analytic function on \( [-1, 1] \) satisfying Assumption 1.1. Let \( \hat{f}(l) \), with \( 0 \leq l \leq m \), be its Gegenbauer coefficients defined by

\[
\hat{f}(l) = \frac{1}{h_{l}^{2}} \int_{-1}^{1} (1 - x^{2})^{-\frac{1}{2}} f(x) C_{l}^{\gamma}(x) \, dx.
\]

Furthermore, assume \( \lambda = \gamma m \), where \( \gamma \) is a positive constant, then the regularization error defined by

\[
RE(\lambda, m) = \max_{|x| \leq 1} \left| f(x) - \sum_{l=0}^{m} \hat{f}(l) C_{l}^{\gamma}(x) \right|
\]

satisfies

\[
RE(\gamma, m) \leq Cq^{m}
\]

Z. Chen, C.-W. Shu / Journal of Computational and Applied Mathematics 265 (2014) 83–95

where \( q \) is given by

\[
q = \frac{(1 + 2\gamma)^{1 + 2\alpha}}{p^{2(1 + 2\gamma)}\gamma^{1 + \gamma}}
\]

which is always less than 1. In particular, if \( \gamma = 1 \) and \( m = \beta N \) where \( \beta \) is a positive constant, then

\[
\text{RE}(\beta N, \beta N) \leq C q^N
\]

with

\[
q = \left( \frac{27}{32p} \right)^{\beta}.
\]

**Proof.** The proof of this lemma can be found in [1].

### 3. Truncation error

Consider the function in the form of

\[
f(x) = a(x) + b(x)(1 + x)^s
\]

where \( s \) is a given constant \( 0 < s < \frac{p}{q} < 1 \) with relatively prime positive integers \( p \) and \( q \), and \( a(x) \) and \( b(x) \) are analytic functions satisfying Assumption 1.1.

We assume that the point values \( f(x_i) \) on \( N \) Gaussian points are given. We are interested in recovering the first \( m \) coefficients in the Gegenbauer expansion of \( f(x) \). For the function \( f(y) = f(x(y)) \), we have the usual Gegenbauer expansion with the basis \( \{ C_l^s(y) \} \):

\[
f(x(y)) = \bar{f}(y) = \sum_{l=0}^{\infty} \hat{f}^s(l) C_l^s(y)
\]

where the Gegenbauer coefficients \( \hat{f}^s(l) \) are given by

\[
\hat{f}^s(l) = \frac{1}{h_l^s} \int_{-1}^{1} (1 - y^2)^{\frac{l}{2}} f(x(y)) C_l^s(y) dy. \tag{3.1}
\]

Our goal here is to find a good approximation to the first \( m \sim N \) Gegenbauer coefficients \( \hat{f}^s(l) \) in this expansion. Based on the known point values \( \{ f(x_i) \} \), we define not the usual interpolation polynomial \( I_N f(x) \), but

\[
(I^s_N f)(y) = I_N \left( (1 - y^2)^{\frac{l}{2}} f(x) \right) = I_N \left( (1 - y^2)^{\frac{l}{2}} (1 + y)^{\frac{l - 1}{2}} f(x) \right). \tag{3.2}
\]

Intuitively, the function being interpolated has about \( \frac{\hat{a}}{q} \sim N \) continuous derivatives, hence the interpolation would produce nice error estimates.

Our candidate for approximating the Gegenbauer coefficients \( \hat{f}^s(l) \) is:

\[
\hat{g}^s(l) = \frac{1}{h_l^s} \int_{-1}^{1} (1 + y)^{\frac{q - 1}{2}} f(x(y)) C_l^s(y) dy. \tag{3.3}
\]

**Definition 3.1.** The truncation error is defined as

\[
TE(\lambda, m, N) = \max_{-1 \leq y \leq 1} \left| \hat{f}^s(l) - \hat{g}^s(l) C_l^s(y) \right|. \tag{3.4}
\]

Let us start to estimate the truncation error.

**Lemma 3.2.** The truncation error is bounded by

\[
TE(\lambda, m, N) \leq \frac{(m + 1)(m + \lambda)(m + 2\lambda)}{\sqrt{m! r(2A)}} \left\| \frac{d^r}{dx^r} \left\{ (1 - y^2)^{\frac{l}{2}}(1 + y)^{\frac{q - 1}{2}} f(x) \right\} \right\|_{L^\infty}
\]

where \( r = \left\lceil \frac{q}{2} \right\rceil - 1.

Proof. We have

$$
TE(\lambda, m, N) = \max_{|s| \leq 1} \sum_{l=0}^{m} \left( \hat{f}(l) - \hat{g}(l) \right) C^l \hat{y}
$$

$$
\leq \sum_{l=0}^{m} C^l \frac{1}{h_l^2} \int_{-1}^{1} \left\{ (1 - y^2)^{\frac{l}{2}} f(x) - (1 + y)^{\frac{g+1}{2}} I_N \left( \frac{(1 - y)(x)^{\frac{g-1}{2}}}{(1 + y)(x)} \right) \right\} C^l \hat{y} dy
$$

$$
\leq \sum_{l=0}^{m} C^l \frac{(1 + y)^{\frac{g-1}{2}}}{h_l^2} \left\{ (1 - y^2)^{\frac{l}{2}} - I_N \left( \frac{(1 - y)(x)^{\frac{g-1}{2}}}{(1 + y)(x)} \right) \right\} dy
$$

$$
\leq C \sum_{l=0}^{m} C^l \frac{1}{h_l^2} \int_{-1}^{1} \left( 1 - y^2 \right)^{\frac{l}{2}} f(x) dx dy
$$

$$
= C \sum_{l=0}^{m} C^l \frac{(1 - y^2)^{\frac{l}{2}} f(x)}{(1 + y)^{\frac{g-1}{2}}} \left( \int_{-1}^{1} \frac{1}{(1 + y)(x)} \right) \| \cdot \|_{l^2}
$$

$$
\leq C \left( \cdot \left( \int_{-1}^{1} \alpha(x) dx \right) \right)^{\frac{1}{2}} \left( \int_{-1}^{1} \alpha(x) dx \right)^{\frac{1}{2}} 
$$

where $\alpha(x)$ is the weight function of the basis $\{\phi_k\}$, and the definitions of $\hat{f}(l)$ and $\hat{g}(l)$ in (3.3) are used in the second inequality; (2.2) is used in the second and third inequalities, the Cauchy-Schwarz inequality and $\frac{dx}{dy} = \frac{(1+y)^{g-1}}{(1+y)(x)}$ are used in the fourth inequality, Lemma 2.4 is used in the sixth inequality, and (2.1) and (2.5) are used in the last inequality.

From the estimates in Lemma 2.9, we obtain

**Theorem 3.3 (The Exponential Decay of the Truncation Error).** Let $\lambda = \alpha N, m = \beta N$ with $0 < \alpha, \beta < 1$, then

$$
TE(\alpha N, \beta N, N) \leq C \rho_q N_T
$$

with

$$
q_T = \left( \frac{\beta + 2\alpha}{\beta \alpha} \right)^{\frac{q}{\alpha}} \cdot \frac{\alpha}{2q(1 + \delta)} \quad , \quad \delta = \frac{1}{e} \left( 1 + \frac{1}{2qAqP} \right).
$$

When we choose $\beta = \gamma \alpha$, i.e. $m = \gamma N$, we have

$$
q_T = \left( \frac{y + 2y^2 + 2}{y^2} \right)^{\frac{q}{2}} \cdot \frac{\alpha}{2q(1 + \delta)}.
$$

If we choose $\alpha$ to satisfy

$$
\alpha < \frac{2q}{(1 + \delta)} \left( \frac{y}{(2 + y)(2 + y)} \right)^{q},
$$

then $q_T < 1$.

Proof. We have

$$
TE(\lambda, m, N) = \max_{|s| \leq 1} \sum_{l=0}^{m} \left( \hat{f}(l) - \hat{g}(l) \right) C^l \hat{y}
$$

$$
\leq C \left( m + 1 \right)(m + \lambda) \frac{r(m + 2\lambda)}{\sqrt{m! r(2\lambda)}} \cdot \frac{1}{N} \left\| \frac{d}{dx} \left( (1 - y^2)^{\frac{l}{2}} f(x) \right) \frac{1}{(1 + y)} \right\|_{l^\infty}
$$

Z. Chen, C.-W. Shu / Journal of Computational and Applied Mathematics 265 (2014) 83–95

\[
\sum_{j=1}^{\frac{2}{q}} \leq C(\rho)(m + 1)(m + \lambda j \Gamma(m + 2\lambda)) \frac{1}{N^t (1 + y)^{\frac{q}{2} - 1} (2\lambda A \xi^{t})}
\]

where we have used Lemma 3.2 in the second inequality, and we need Lemma 2.9 in the third inequality and \( t = \left\lceil \frac{\lambda - \frac{1}{2}}{q} \right\rceil - 1 \leq \frac{\lambda - \frac{1}{2}}{q} \) in the fourth inequality. Therefore, we have

\[
TE(\alpha, m, N) \leq C(\rho) \left( (m + 1)(m + \lambda j) \Gamma(m + 2\lambda) \frac{(2\lambda A (1 + \delta))^t}{\sqrt{m!} \Gamma(2\lambda)} \frac{(2\lambda A (1 + \delta))^t}{N} \right)
\]

hence

\[
TE(\alpha, \beta, N) \leq C(\rho) \left( \frac{(2\alpha + \beta)2^{\alpha}}{\beta^2(2\alpha)2\alpha} \right)^{N} \left( (2\lambda A (1 + \delta))^t \right)
\]

\[
\leq C(\rho) \left( \frac{(2\alpha + \beta)2^{\alpha}}{\beta^2(2\alpha)2\alpha} \right)^{N} \left( (2\lambda A (1 + \delta))^t \right)^{N}
\]

\[
= C(\rho) q^{N}
\]

with 

\[
q = \left( \frac{(2\alpha + \beta)2^{\alpha}}{\beta^2(2\alpha)2\alpha} \right) \left( \frac{\alpha}{2q(1 + \delta)} \right)^{\frac{q}{2}}
\]

where we have used Stirling's formula (2.4) in the first inequality and \( \frac{\lambda}{q} - 2 \leq t \leq t = \left\lceil \frac{\lambda}{q} \right\rceil - 1 \leq \frac{\lambda}{q} \) in the second inequality.

4. Regularization error and the main theorem

In this section we would like to establish the error estimates for approximating \( f(x) = a(x) + b(x)(1 + x)^s \) on \([-1, 1]\) by its Gegenbauer expansion based on the Gegenbauer polynomials \( C_n^{\lambda}(y(x)) \). Since our goal is to remove the Gibbs phenomenon, we will use the maximum norm. As shown in the last section, we would need both \( \lambda \) and \( m \) to grow linearly with \( N \) in order to obtain exponential accuracy in the maximum norm. Thus, we will consider the case of large \( \lambda \) and \( m \) in this section.

We will assume that \( a(x) \) and \( b(x) \) are analytic functions on \([-1, 1]\) satisfying Assumption 1.1. We would like to estimate the regularization error in the maximum norm. We consider the Gegenbauer partial sum of the first \( m \) terms for the function \( f(x) \) given by:

\[
f^{m}(x) = \sum_{j=0}^{m} \hat{f}^{(j)} C_{j}^{\lambda}(y(x)),
\]

with the Gegenbauer coefficients \( \hat{f}^{(j)} \) defined by (3.1).

Definition 4.1. The regularization error is defined by

\[
RE(\alpha, m) = \max_{|x| \leq 1} \left| f(x) - \sum_{j=0}^{m} \hat{f}^{(j)} C_{j}^{\lambda}(y(x)) \right|
\]

\[
= \max_{|x| \leq 1} \left| f(y(x)) - \sum_{j=0}^{m} \hat{f}^{(j)} (C_{j}^{\lambda}(y) \right) .
\]

Since the function \( f(x) \) is analytic with respect to the variable \( y \), we can get the following result from Lemma 2.10.

Theorem 4.2 (The Exponential Decay of the Regularization Error). For the function \( f(x) = a(x) + b(x)(1 + x)^s \), with analytic functions \( a(x) \) and \( b(x) \) satisfying Assumption 1.1, if we assume \( \lambda = \alpha N \) and \( m = \gamma N \), then

\[
\max_{|x| \leq 1} \left| f(x) - \sum_{j=0}^{m} \hat{f}^{(j)} C_{j}^{\lambda}(y(x)) \right| \leq C N^{-\gamma q}
\]

92  
Z. Chen, C.-W. Shu / Journal of Computational and Applied Mathematics 265 (2014) 83–95  

where  
\[ 
q_r = \left( \frac{(y + 2y^2 + 2)}{(y^2 + 2y + 1)} \right)^{\alpha} 
\]  
which is always less than 1.  

We can then combine the estimates for the truncation error and the regularization error to obtain the main theorem:  

### Theorem 4.3 (Removal of the Gibbs Phenomenon).  
Consider a function in the form of \( f(x) = a(x) + b(x)(1 + x)^{\delta} \), with given fractional constant \( 0 < s < \frac{1}{\alpha} < 1 \), and \( a(x) \) and \( b(x) \) are analytic functions satisfying Assumption 1.1. Assume that the point values \( f(x_j) \) at \( N \) Gaussian points are given. Let \( g^2(I) \), with \( 0 \leq l \leq m \), be the Gegenbauer expansion coefficients, defined in (3.2) and (3.3). Then, for \( x = \xi \) and \( m = \gamma \) with  
\[
\alpha < \frac{2q}{(1 + \delta)(y^2 + 2y)} \cdot \left( \frac{y^2}{(2 + y)(2 + y)} \right)^{\frac{q}{8}} \cdot \delta = \frac{1}{e} \left( 1 + \frac{1}{24q\alpha p} \right), 
\]  
we have  
\[
\max_{-1 \leq s \leq 1} \left| f(x) - \sum_{n=0}^{m} g^2(\xi_j)(y(x)) \right| \leq C \left( q_N^T + q_N^R \right) 
\]  
where  
\[ 
q_T = \left( \frac{(y + 2y^2 + 2) \left( \frac{\alpha}{2q(1 + \delta)} \right)^{\frac{1}{\alpha}}}{y^r} \right) < 1, 
\]  
\[ 
q_R = \left( \frac{(y + 2y^2 + 2)}{(y^2 + 2y + 1)} \right)^{\alpha} < 1. 
\]  

**Proof.** We can get the proof by combining the results of Theorems 3.3 and 4.2.  

**Remark 4.4.** In the proof, no attempt has been made to optimize the parameters.  

### 5. Numerical results  
In this section, we give two numerical examples to illustrate the result. We have tested both the Fourier collocation and the Chebyshev collocation.  

**Example 5.1.** We take the function  
\[ 
f(x) = \cos(x) + \sin(x)\sqrt{1 + x} 
\]  
and try to recover the pointwise values of this function over \([-1, 1]\).  

First, we consider the Fourier collocation. We assume the point values \( \{f(x_j)\} \) on \( 2N + 1 \) uniform points:  
\[ 
x_i = \frac{2i}{2N + 1}, \quad i = -N, \ldots, N 
\]  
are given. The parameters are chosen as  
\[
\lambda = 0.2N, \quad m = 0.075N. 
\]  

The choice of the parameters here and below are guided by the theory in previous sections. Second, we consider the Chebyshev collocation. We assume the point values \( \{f(x_i)\} \) on the \( N + 1 \) Chebyshev collocation points:  
\[ 
x_i = \cos\left( \frac{\pi(2i + 1)}{2N + 2} \right), \quad i = 0, \ldots, N 
\]  
are given. The parameters are chosen as  
\[
\lambda = 0.2N, \quad m = 0.1N. 
\]  

In Fig. 5.1 we show the errors, in the logarithmic scale, for \( N = 40, 80, 160, 320 \) and \( 640 \). In Table 5.1, we show the maximum errors for each \( N \) and the orders of convergence. We can clearly find out that both the Fourier case and the Chebyshev case provide exponentially decaying errors, and the errors are smaller in the Chebyshev case than in the Fourier case.

Z. Chen, C.-W. Shu / Journal of Computational and Applied Mathematics 265 (2014) 83–95

93

(a) Fourier case \( (\lambda = 0.2N, m = 0.075N) \)

\[
\text{Pointwise Error}
\]

(b) Chebyshev case \( (\lambda = 0.2N, m = 0.1N) \)

Fig. 5.1. Pointwise errors in the logarithm scale.

Table 5.1  
Maximum error table.

\[
\begin{array}{|c|c|c|}
\hline
N & \text{Fourier case} & \text{Chebyshev case} \\
\hline
 & \lambda = 0.2N, m = 0.075N & \lambda = 0.2N, m = 0.1N \\
\hline
L^1 \text{ error} & \text{Order} & L^1 \text{ error} & \text{Order} \\
\hline
40 & 3.46E-01 & 1.34E-01 \\
80 & 2.01E-02 & 4.11 & 2.49E-03 & 5.75 \\
160 & 2.19E-05 & 9.84 & 5.17E-08 & 15.56 \\
320 & 5.93E-11 & 18.49 & 3.58E-19 & 37.07 \\
640 & 1.82E-18 & 24.96 & 8.48E-43 & 78.48 \\
\hline
\end{array}
\]

Example 5.2. We take the function 

\[
f(x) = \cos(x) + \sin(x)(1 + x)^{\frac{1}{3}}
\]

and again try to recover the pointwise values of this function over \([-1, 1]\).

First, we consider the Fourier collocation. We assume the point values \( f(x_j) \) on the \( 2N + 1 \) uniform points \( x_j \) given by (5.2) are given. The parameters are chosen as

\[
\lambda = 0.2N, \quad m = 0.025N.
\]

94  
Z. Chen, C.-W. Shu / Journal of Computational and Applied Mathematics 265 (2014) 83–95  

![Figure 5.2](image_url)  
(a) Fourier case \((\lambda = 0.2N, m = 0.025N)\).  
(b) Chebyshev case \((\lambda = 0.2N, m = 0.05N)\).  

**Fig. 5.2.** Pointwise errors in the logarithm scale.  

### Table 5.2  
**Maximum error table.**  

| N   | Fourier case                          | Chebyshev Case                     |  
|-----|--------------------------------------|------------------------------------|  
|     | \(\lambda = 0.2N, m = 0.025N\)      | \(\lambda = 0.2N, m = 0.05N\)     |  
|     | \(L^\infty\) error | Order        | \(L^\infty\) error | Order        |  
| 40  | 7.956E–01          | 3.49E–01    |                      |              |  
| 80  | 3.906E–01          | 0.93        | 2.55E–01            | –0.59        |  
| 160 | 6.00E–01           | –0.60       | 4.25E–02            | 3.63         |  
| 320 | 5.16E–02           | 3.54        | 3.95E–05            | 10.07        |  
| 640 | 6.36E–05           | 9.66        | 1.15E–12            | 25.03        |  

Second, we consider the Chebyshev collocation. We assume the point values \(f(\xi_i)\) on the \(N + 1\) Chebyshev collocation points given by (5.4) are given. The parameters are chosen as  
\[
\lambda = 0.2N, \quad m = 0.05N.
\]  

In Fig. 5.2 we show the errors, in the logarithm scale, for \(N = 40, 80, 160, 320\) and \(640\). In Table 5.2, we show the maximum errors for each \(N\) and the orders of convergence. We again clearly see that both the Fourier case and the Chebyshev case provide exponentially decaying errors, and errors are smaller in the Chebyshev case than in the Fourier case.

Remark 5.3. No attempt is made to optimize the parameters for each individual case in the numerical examples.

6. Concluding remarks

We have generalized the Gegenbauer polynomial based technique [5] to reconstruct approximations with exponential accuracy in the maximum norm, based on collocation point values of the function, from the case of piecewise analytic functions where each piece the function is analytic in the closed interval, to piecewise functions where in each piece the function is analytic only on the open interval with end point singularities. Numerical results are provided to demonstrate the theory. This technique has potential applications for solving certain partial differential equations whose solutions have such end point singularities, and also in accurate numerical integration and differentiation of such functions.

References

[1] D. Gottlieb, C.-W. Shu, A. Solomonoff, H. Vandeven, On the Gibbs phenomenon I: recovering exponential accuracy from the Fourier partial sum of a non-periodic analytic function, Journal of Computational and Applied Mathematics 43 (1992) 81–98.  
[2] D. Gottlieb, C.-W. Shu, Resolution properties of the Fourier method for discontinuous waves, Computer Methods in Applied Mechanics and Engineering 116 (1994) 27–37.  
[3] D. Gottlieb, C.-W. Shu, On the Gibbs phenomenon III: recovering exponential accuracy in a sub-interval from a spectral partial sum of a piecewise analytic function, SIAM Journal on Numerical Analysis 33 (1996) 280–290.  
[4] D. Gottlieb, C.-W. Shu, On the Gibbs phenomenon V: recovering exponential accuracy in a sub-interval from a Gegenbauer partial sum of a piecewise analytic function, Mathematics of Computation 64 (1995) 1081–1095.  
[5] D. Gottlieb, C.-W. Shu, On the Gibbs phenomenon V: recovering exponential accuracy from collocation point values of a piecewise analytic function, Numerische Mathematik 71 (1995) 511–526.  
[6] D. Gottlieb, C.-W. Shu, On the Gibbs phenomenon and its resolution, SIAM Review 30 (1997) 644–668.  
[7] H. Bateman, Higher Transcendental Functions, Vol. 2, McGraw-Hill, 1953.  
[8] E. Tadmor, The exponential accuracy of Fourier and Chebyshev differencing methods, SIAM Journal on Numerical Analysis 23 (1986) 1–10.  
[9] D. Funaro, Polynomial Approximation of Differential Equations, in: Lecture Notes in Physics, vol. 8, Springer-Verlag, 1992.