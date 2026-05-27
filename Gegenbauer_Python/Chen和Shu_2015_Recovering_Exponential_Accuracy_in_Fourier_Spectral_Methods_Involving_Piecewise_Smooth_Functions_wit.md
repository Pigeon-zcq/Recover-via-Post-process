# Recovering Exponential Accuracy in Fourier Spectral Methods Involving Piecewise Smooth Functions with Unbounded Derivative Singularities

Zheng Chen1 · Chi-Wang Shu2 

Received: 18 July 2014 / Revised: 24 February 2015 / Accepted: 26 February 2015 

Published online: 26 March 2015 

© Springer Science+Business Media New York 2015 

Abstract Fourier spectral methods achieve exponential accuracy both on the approximation level and for solving partial differential equations, if the solution is analytic. If the solution is discontinuous but piecewise analytic up to the discontinuities, Fourier spectral methods produce poor pointwise accuracy, but still maintain exponential accuracy after post-processing (Gottlieb and Shu in SIAM Rev 30:644–668, 1997) . In Chen and Shu (J Comput Appl Math 265:83–95, 2014), an extended technique is provided to recover exponential accuracy for functions which have end-point singularities, from the knowledge of point values on standard collocation points. In this paper, we develop a technique to recover exponential accuracy from the first N Fourier coefficients of functions which are analytic in the open interval but have unbounded derivative singularities at end points. With this post-processing method, we are able to obtain exponential accuracy of spectral methods applied to linear transport equations involving such functions. 

Keywords Spectral method  Exponential accuracy  Fourier coefficients  Gegenbauer expansion · Transport equation · Variable coefficients · Singular initial conditions · Noise 

# 1 Introduction

In this paper, we are concerned with the accuracy of spectral methods when applied to problems involving piecewise smooth functions with unbounded derivative singularities. We 

Research supported by NSF Grants DMS-1112700 and DMS-1418750, and AFOSR Grant F49550-12-1-0399. 

investigate the issue of overcoming the Gibbs phenomenon, which describes how a global spectral approximation of a piecewise analytic function behaves at the jump discontinuity. A prototype is to use Fourier series to approximate an analytic but non-periodic function u(x) on interval 1, 1 , which has discontinuities at the boundaries of the interval when extended periodically with period 2. The Fourier partial sum using the first 2N 1 modes 

$$
u _ {N} (x) = \sum_ {| k | \leq N} \tilde {u} _ {k} e ^ {i k \pi x}, \tag {1.1}
$$

with Fourier coefficients $\tilde { u } _ { k }$ defined by 

$$
\tilde {u} _ {k} = \frac {1}{2} \int_ {- 1} ^ {1} u (x) e ^ {- i k \pi x} d x, \tag {1.2}
$$

has large oscillations near the jumps, which are not improved as the number of terms in the partial sum increases. In smooth regions away from the discontinuities, convergence is only first order. Therefore, there is no convergence in the maximum norm. This is the so-called Gibbs phenomenon. 

In [9–12,14], Gottlieb et al. developed a general framework to overcome this difficulty. This technique recovers exponential accuracy in the maximum norm for any (sub-)interval of analyticity (up to and including the boundaries of this interval), from the knowledge of either the first N spectral expansion coefficients, or the point values at N standard collocation points. This means that exponential accuracy is recovered at all points, including at the actual discontinuity points (the left and right limits at these points), if the locations of these discontinuity points are known. If the locations of these discontinuity points are not known exactly but are known to be within certain fixed intervals, then exponential accuracy can be recovered from any interval which does not overlap with these fixed intervals containing the discontinuities. In this framework, an important tool is the set of Gegenbauer polynomials, which are orthogonal in the interval [ 1, 1] with the weight $( 1 - x ^ { 2 } ) ^ { \lambda - \frac { 1 } { 2 } }$ . The key in this technique is that the parameter λ in the weight function as well as the number of terms m retained in the Gegenbauer expansion should both be chosen proportional to N . For a review, we refer to [13]. These methods are widely used as “reconstruction” or “post-processing” techniques to recover exponential accuracy for point values or coefficients based on the spectral approximation, such as in recovering high order information of the discontinuous solutions of scalar nonlinear hyperbolic PDEs [18], and in the simulation of sophisticated problems [8]. These techniques have been successfully applied to the field of image reconstruction [2–5] as well. The Gegenbauer basis has also been successful in recovering lost order of accuracy in other types of approximations, such as weighted essentially non-oscillatory (WENO) solutions of hyperbolic PDEs [15], and in radial basis functions approximations of linear and nonlinear hyperbolic PDEs [17]. 

Besides piecewise analytic functions, functions with end-point singularities exist in many applications. Most fractional differential equations have singular solutions. Many standard numerical methods solving fractional differential problems give poor accuracy, due to the lack of regularities. Therefore, it is important to provide a way to recover accuracy of solutions as well as to obtain high order accuracy at the approximation level. In [1], Adcock et al. focused on the approximation of functions which are analytic on a compact interval except at the end-points, and utilized variable transform methods. They introduced two new mappings from the original interval to either semi-infinite or infinite interval, and provided approximation procedure from sampling information on the new region. The two new mappings, compared with the standard transformations, vastly improve resolution power, and achieve root exponential decays, with proper choice of parameters in the mappings. 

We are interested in developing post-processors to recover high order accuracy for such functions as well. This task is significantly more difficult than the recovery of accuracy for piecewise analytic functions. In particular, the extension of the technique in [13] to functions with end-point singularities is highly non-trivial. In [7], we made this extension to handle spectral collocation methods for such functions. The reconstruction procedure is performed on functions of the following form 

$$
f (x) = a (x) + b (x) (1 + x) ^ {s}, \quad x \in [ - 1, 1 ] \tag {1.3}
$$

where s is a given fractional constant 

$$
0 <   s = \frac {p}{q} <   1 \tag {1.4}
$$

with relatively prime integers p and q, and $a ( x )$ and $b ( x )$ are both analytic but unknown functions. Such functions lack regularities, and the derivatives blow up at the left end point. With this extension, exponential accuracy can be obtained from standard collocation point values of such functions, which is different from sampling on the mapped region as in [1], by properly choosing the parameters λ and m to be linearly dependent on N . A crucial modification of the choice of parameters and a more refined estimate were necessary to balance the terms in the truncation error for this analysis. 

In this study, we are interested in recovering high order accuracy from the first $2 N + 1$ Fourier coefficients for functions in the form (1.3). The objective is to extract the hidden information from the truncated Fourier series (1.1) and recover exponentially accurate point values at every point including at the singularities. 

As in [7], we assume that the analytic functions $a ( x )$ and $b ( x )$ , denoted generically as $c ( x )$ , satisfy the following condition. 

Assumption 1.1 There exists a constant $\rho \ge 1$ and a constant $C ( \rho )$ such that, for every $k \geq 0$ , 

$$
\max _ {- 1 \leq x \leq 1} \left| \frac {d ^ {k} c (x)}{d x ^ {k}} \right| \leq C (\rho) \frac {k !}{\rho^ {k}}.
$$

This is a standard assumption for analytic functions, where ρ is the distance from the interval $[ - 1 , 1 ]$ to the nearest singularity of the function $c ( x )$ in the complex plane. 

We will use the following one-to-one transformation between $x \in [ - 1 , 1 ]$ and $y \in [ - 1 , 1 ]$ : 

$$
(2 ^ {q - 1} (1 + x)) ^ {\frac {1}{q}} = 1 + y \tag {1.5}
$$

where $q$ is defined in (1.4). 

The function $F ( y ) = f ( x ( y ) )$ of the variable y has its usual Gegenbauer expansion under the basis $\{ C _ { l } ^ { \lambda } ( y ) \}$ : 

$$
f (x (y)) = F (y) = \sum_ {l = 0} ^ {\infty} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y)
$$

with the Gegenbauer coefficients ${ \hat { f } } ^ { \lambda } ( l )$ given by 

$$
\hat {f} ^ {\lambda} (l) = \frac {1}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} F (y) C _ {l} ^ {\lambda} (y) d y \tag {1.6}
$$

where the precise value of the normalization constant $h _ { l } ^ { \lambda }$ will be given later by (2.4) in Definition 2.1. 

Our goal is to find a good approximation to the first m ∼ N Gegenbauer coefficients ${ \hat { f } } ^ { \lambda } ( l )$ in (1.6), denoted as $\hat { g } ^ { \lambda } ( l )$ [defined later in (3.1)], from the given Fourier coefficients. 

We will then obtain the approximation of f (x) using these m $\sim N$ terms of its Gegenbauer expansion: 

$$
f _ {N} ^ {m, \lambda} (x) = \sum_ {l = 0} ^ {m} \hat {g} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)).
$$

With proper choice of the parameters λ and m, the error between the reconstructed approximation $f _ { N } ^ { m , \lambda } ( x )$ and the function $f ( x )$ , measured in the maximum norm, decays exponentially as N increases. Therefore, the reconstruction method provides a way to post-process functions with such singularities from the first $2 N + 1$ accurate Fourier coefficients. 

Another major concern in the reconstruction methods is the possible existence of noise in the data. The Gegenbauer reconstruction techniques in [13] work well with noise. This is because the noise is projected to be very small in the Gegenbauer basis. Similar robustness can also be observed in the new reconstruction methods for piecewise smooth functions with end-point singularities. 

The analysis of the error from reconstruction is divided into two parts: the truncation error and the regularization error. The truncation error measures the difference between the exact Gegenbauer coefficients of $f ( x ( y ) )$ ) with $\lambda \sim N$ , and the approximate Gegenbauer coefficients $\hat { g } ^ { \lambda } ( l )$ obtained by using the truncated Fourier series. This will be investigated in Sect. 3. In relation to the collocation case as in [7], the Galerkin case is more difficult for such singularity cases, which involves different and more complicated analysis. The regularization error measures the difference between the Gegenbauer expansion with $\lambda \sim N$ , using the first m N Gegenbauer coefficients, and the function itself. This error is estimated in Sect. 4. The results for the reconstruction are summarized in Theorem 4.3 in Sect. 4. We also give an analysis of the Fourier Galerkin methods for solving initial problems of linear hyperbolic time-dependent partial differential equations with reconstructions in Sect. 5. Section 6 contains several numerical examples to illustrate our results and robustness to noise. In Sect. 2, we shall give several useful preliminary properties and estimates. Concluding remarks are given in Sect. 7. 

Throughout this paper, we will use C to denote a generic constant either independent of the growing parameters, or depending on them at most in polynomial growth. The details will be indicated clearly in the text. These constants may not take the same value at different places. 

# 2 Preliminaries

In this section, we will introduce the Gegenbauer polynomials and discuss some of their asymptotic behavior (see Bateman [6]). Then we will give some estimates as preparation for the error estimates in Sect. 3. 

Definition 2.1 The Gegenbauer polynomial $C _ { n } ^ { \lambda } ( x )$ , for $\lambda \geq 0$ , is defined by 

$$
(1 - x ^ {2}) ^ {\lambda - \frac {1}{2}} C _ {n} ^ {\lambda} (x) = \frac {(- 1) ^ {n}}{2 ^ {n} n !} G (\lambda , n) \frac {d ^ {n}}{d x ^ {n}} \left[ (1 - x ^ {2}) ^ {n + \lambda - \frac {1}{2}} \right]
$$

where G(λ, n) is given by 

$$
G (\lambda , n) = \frac {\Gamma (\lambda + \frac {1}{2}) \Gamma (n + 2 \lambda)}{\Gamma (2 \lambda) \Gamma (n + \lambda + \frac {1}{2})} \tag {2.1}
$$

for λ > 0, by 

$$
G (0, n) = \frac {2 \sqrt {\pi} (n - 1) !}{\Gamma (n + \frac {1}{2})}
$$

for λ  0 and n  1, and by 

$$
G (0, 0) = 1
$$

for λ = 0 and n = 0. Notice that by this standardization, $C _ { n } ^ { 0 } ( x )$ is defined by: 

$$
C _ {n} ^ {0} (x) = \lim _ {\lambda \rightarrow 0 ^ {+}} \frac {1}{\lambda} C _ {n} ^ {\lambda} (x) = \frac {2}{n} T _ {n} (x), \quad n > 0; \quad C _ {0} ^ {0} (x) = 1,
$$

where $T _ { n } ( x )$ are the Chebyshev polynomials. 

Under this definition we have, for λ > 0, 

$$
C _ {n} ^ {\lambda} (1) = \frac {\Gamma (n + 2 \lambda)}{n ! \Gamma (2 \lambda)}; \tag {2.2}
$$

for λ = 0 and n ≥ 1, 

$$
C _ {n} ^ {0} (1) = \frac {2}{n};
$$

for λ 0 and n 0, 

$$
C _ {n} ^ {0} (1) = 1;
$$

and 

$$
\left| C _ {n} ^ {\lambda} (x) \right| \leq C _ {n} ^ {\lambda} (1), - 1 \leq x \leq 1. \tag {2.3}
$$

The Gegenbauer polynomials are orthogonal under their weight function $( 1 - x ^ { 2 } ) ^ { \lambda - \frac { 1 } { 2 } }$ : 

$$
\int_ {- 1} ^ {1} (1 - x ^ {2}) ^ {\lambda - \frac {1}{2}} C _ {k} ^ {\lambda} (x) C _ {n} ^ {\lambda} (x) d x = \delta_ {k, n} h _ {n} ^ {\lambda}
$$

where, for λ > 0, 

$$
h _ {n} ^ {\lambda} = \pi^ {\frac {1}{2}} C _ {n} ^ {\lambda} (1) \frac {\Gamma (\lambda + \frac {1}{2})}{\Gamma (\lambda) (n + \lambda)}; \tag {2.4}
$$

for λ 0 and n 1, 

$$
h _ {n} ^ {0} = \frac {2 \pi}{n ^ {2}};
$$

for λ 0 and n 0, 

$$
h _ {0} ^ {0} = \pi .
$$

We will need to use the Stirling’s formula and the estimate of $h _ { n } ^ { \lambda }$ for the asymptotics of the Gegenbauer polynomials for large n and λ. 

Lemma 2.2 We have the Stirling’s formula 

$$
(2 \pi) ^ {\frac {1}{2}} x ^ {x + \frac {1}{2}} e ^ {- x} \leq \Gamma (x + 1) \leq (2 \pi) ^ {\frac {1}{2}} x ^ {x + \frac {1}{2}} e ^ {- x + \frac {1}{1 2 x}}, \quad x \geq 1. \tag {2.5}
$$

Lemma 2.3 There exists a constant C independent of λ and n such that 

$$
C ^ {- 1} \frac {\lambda^ {\frac {1}{2}}}{n + \lambda} C _ {n} ^ {\lambda} (1) \leq h _ {n} ^ {\lambda} \leq C \frac {\lambda^ {\frac {1}{2}}}{n + \lambda} C _ {n} ^ {\lambda} (1). \tag {2.6}
$$

We would need to estimate $| | \frac { d ^ { t } } { d x ^ { t } } \left\{ \frac { d ^ { l } } { d y ^ { l } } ( 1 - y ^ { 2 } ) ^ { l + \lambda - \frac { 1 } { 2 } } \frac { d y } { d x } \right\} | | _ { L ^ { \infty } }$ , therefore we need the following preliminaries first. 

Remark 2.4 $\begin{array} { r } { \frac { d ^ { l } } { d y ^ { l } } ( 1 - y ^ { 2 } ) ^ { l + \lambda - \frac { 1 } { 2 } } \frac { d y } { d x } } \end{array}$ has up to t-th derivatives in x, where $t = \lfloor \frac { \lambda + \frac { 1 } { 2 } } { q } \rfloor - 1$ λ 1 , the largest integer below λ+ 12 1. $\frac { \lambda + \frac { 1 } { 2 } } { q } - 1$ q 

It is easy to observe that 

$$
\frac {d ^ {n}}{d x ^ {n}} \left(1 - y (x) ^ {2}\right) ^ {l + \lambda - \frac {1}{2}} = A ^ {n} Y _ {1} ^ {n} Y _ {2} ^ {n} Y _ {3} ^ {n}, \quad 0 \leq n \leq t + 1 \tag {2.7}
$$

where 

$$
A = \frac {2 ^ {q}}{2 q}, \quad Y _ {1} ^ {n} = \left(1 - y (x) ^ {2}\right) ^ {l + \lambda - \frac {1}{2} - q n}, \quad Y _ {2} ^ {n} = (1 - y (x)) ^ {n (q - 1)},
$$

and $Y _ { 3 } ^ { n }$ satisfies the following recursive relation: 

$$
Y _ {3} ^ {0} = 1
$$

$$
Y _ {3} ^ {n + 1} = - \left[ (2 l + 2 \lambda - q n - n - 1) y + n (q - 1) \right] Y _ {3} ^ {n} + (1 - y ^ {2}) \frac {d}{d y} Y _ {3} ^ {n}, \quad 0 \leq n \leq t.
$$

It is easy to show that $Y _ { 3 } ^ { n }$ is an n-th degree polynomial of y. We have the following estimate on $Y _ { 3 } ^ { n }$ . 

Lemma 2.5 We have, for $0 \leq n \leq t + 1$ , 

$$
\left| Y _ {3} ^ {n} \right| \leq (2 l + 2 \lambda) ^ {n}, \quad y \in [ - 1, 1 ] \tag {2.8}
$$

Proof The proof can be found in [7]. In the proof, $Y _ { 3 } ^ { n }$ is rewritten as $\begin{array} { r } { Y _ { 3 } ^ { n } = \sum _ { i = 0 } ^ { n } a _ { i } y ^ { i } } \end{array}$ and we denote $\begin{array} { r } { S _ { n } = \sum _ { i = 0 } ^ { n } \left| a _ { i } \right| } \end{array}$ . The proof also provides an estimate for $S _ { n }$ : 

$$
S _ {n} \leq (2 l + 2 \lambda) ^ {n} \tag {2.9}
$$

This will be used later in the estimation. 

For $l \geq 1$ , 

$$
\frac {d ^ {t}}{d x ^ {t}} \left\{\frac {d ^ {l}}{d y ^ {l}} \left(1 - y ^ {2} (x)\right) ^ {l + \lambda - \frac {1}{2}} \frac {d y}{d x} \right\} = \frac {d ^ {l - 1}}{d y ^ {l - 1}} \frac {d ^ {t + 1}}{d x ^ {t + 1}} \left\{(1 - y ^ {2} (x)) ^ {l + \lambda - \frac {1}{2}} \right\} \tag {2.10}
$$

Lemma 2.6 We have the following estimate for $l \geq 1$ , 

$$
\left| \frac {d ^ {l - 1}}{d y ^ {l - 1}} \frac {d ^ {t + 1}}{d x ^ {t + 1}} \left\{\left(1 - y ^ {2} (x)\right) ^ {l + \lambda - \frac {1}{2}} \right\} \right| \leq C A ^ {t + 1} 2 ^ {l + \lambda} (l + \lambda) ^ {l + t}, \quad y \in [ - 1, 1 ] \tag {2.11}
$$

$$
\frac {d ^ {l - 1}}{d y ^ {l - 1}} \frac {d ^ {t + 1}}{d x ^ {t + 1}} \left\{\left(1 - y ^ {2} (x)\right) ^ {l + \lambda - \frac {1}{2}} \right\} = A ^ {t + 1} \frac {d ^ {l - 1}}{d y ^ {l - 1}} \left\{Y _ {1} ^ {t + 1} Y _ {2} ^ {t + 1} Y _ {3} ^ {t + 1} \right\}
$$

For simplicity, we denote 

$$
\frac {d ^ {i}}{d y ^ {i}} \{Y _ {1} ^ {t + 1} Y _ {2} ^ {t + 1} Y _ {3} ^ {t + 1} \} = X _ {1} ^ {i} X _ {2} ^ {i}, 0 \leq i \leq l - 1
$$

where, 

$$
X _ {1} ^ {i} = (1 - y ^ {2}) ^ {l + \lambda - \frac {1}{2} - q (t + 1) - i}
$$

and $X _ { 2 } ^ { i }$ satisfies the following recursive relationship: 

$$
X _ {2} ^ {0} = Y _ {2} ^ {t + 1} Y _ {3} ^ {t + 1}
$$

$$
X _ {2} ^ {i + 1} = - \left[ 2 (l + \lambda - q (t + 1) - \frac {1}{2} - i) y \right] X _ {2} ^ {i} + (1 - y ^ {2}) \frac {d}{d y} X _ {2} ^ {i}, \quad 0 \leq i <   l - 2.
$$

It is easy to find out that $X _ { 2 } ^ { i }$ is a polynomial of degree $( t + 1 ) q + i$ . We only need to prove that 

$$
\left| X _ {2} ^ {l - 1} \right| \leq C 2 ^ {l + \lambda} (l + \lambda) ^ {l + t} \tag {2.12}
$$

We can rewrite it as $\begin{array} { r } { X _ { 2 } ^ { i } = \sum _ { j = 0 } ^ { ( t + 1 ) q + i } \beta _ { j } ^ { i } y ^ { j } } \end{array}$ . In order to get the upper bound for $X _ { 2 } ^ { l - 1 }$ as in 2 (2.12), we need to estimate $X _ { 2 } ^ { 0 }$ =0 j first. Let us rewrite $\begin{array} { r } { Y _ { 3 } ^ { t + 1 } = \sum _ { i = 0 } ^ { t + 1 } a _ { i } y ^ { i } } \end{array}$ 2 with its coefficients $a _ { i }$ . 

$$
\begin{array}{l} X _ {2} ^ {0} = (1 - y) ^ {(t + 1) (q - 1)} Y _ {3} ^ {t + 1} \\ = \sum_ {j = 0} ^ {(t + 1) (q - 1)} \binom {(t + 1) (q - 1)} {j} (- y) ^ {j} \sum_ {i = 0} ^ {t + 1} a _ {i} y ^ {i} \\ = \sum_ {j = 0} ^ {(t + 1) (q - 1)} \sum_ {i = 0} ^ {t + 1} \binom {(t + 1) (q - 1)} {j} (- 1) ^ {j} a _ {i} y ^ {i + j} \\ \end{array}
$$

Then we measure the sum of the coefficients of $X _ { 2 } ^ { 0 }$ , 

$$
\begin{array}{l} \sum_ {k = 0} ^ {(t + 1) q} | \beta_ {k} ^ {0} | \leq \sum_ {j = 0} ^ {(t + 1) (q - 1)} \sum_ {i = 0} ^ {t + 1} \binom {(t + 1) (q - 1)} {j} | a _ {i} | \\ = \sum_ {j = 0} ^ {(t + 1) (q - 1)} \binom {(t + 1) (q - 1)} {j} \sum_ {i = 0} ^ {t + 1} | a _ {i} | \\ = 2 ^ {(t + 1) (q - 1)} S _ {t + 1} \\ \leq 2 ^ {(t + 1) q} (l + \lambda) ^ {t + 1} \\ \end{array}
$$

Using induction, we can easily get similar result for $X _ { 2 } ^ { l - 1 }$ 

$$
\sum_ {k = 0} ^ {(t + 1) q + l - 1} | \beta_ {k} ^ {l - 1} | \leq 2 ^ {(t + 1) q + l - 1} (l + \lambda) ^ {t + l}
$$

which implies that 

$$
\begin{array}{l} | X _ {2} ^ {l - 1} | \leq 2 ^ {(t + 1) q + l - 1} (l + \lambda) ^ {t + l} \\ \leq C 2 ^ {\lambda + l} (l + \lambda) ^ {t + l} \\ \end{array}
$$

Thus, we complete the proof. 

Remark 2.7 The result in Lemma 2.6 is also true for $l = 0 ,$ . Hence we have 

$$
\left\| \frac {d ^ {t}}{d x ^ {t}} \left\{\frac {d ^ {l}}{d y ^ {l}} (1 - y ^ {2} (x)) ^ {l + \lambda - \frac {1}{2}} \frac {d y}{d x} \right\} \right\| _ {L ^ {\infty}} \leq C A ^ {t + 1} 2 ^ {l + \lambda} (l + \lambda) ^ {l + t}, \quad l \geq 0 \tag {2.13}
$$

# 3 Truncation Error

In this section, we will establish the error estimate for replacing the Gegenbauer coefficients ${ \hat { f } } ^ { \lambda } ( l )$ by the new approximated coefficients $\hat { g } ^ { \lambda } ( l )$ , defined later in (3.1), in the Gegenbauer expansion. 

Consider the function in the form of 

$$
f (x) = a (x) + b (x) (1 + x) ^ {s}
$$

where s is a given constant $\textstyle 0 < s = { \frac { p } { q } } < 1$ with relatively prime positive integers $p$ and $q _ { \mathrm { : } }$ , and $a ( x )$ and $b ( x )$ are analytic functions satisfying Assumption 1.1. 

We assume that the Fourier coefficients $\tilde { f } _ { n } \ ( - N \leq n \leq N )$ are given. Thus we have its truncated Fourier series 

$$
f _ {N} (x) = \sum_ {| n | \leq N} \tilde {f} _ {n} e ^ {i n \pi x}.
$$

We are interested in recovering the first m coefficients in the Gegenbauer expansion of $f ( x )$ . For the function $F ( y ) = f ( x ( y ) )$ ), we have the usual Gegenbauer expansion with the basis $\{ C _ { l } ^ { \lambda } ( y ) \}$ : 

$$
f (x (y)) = F (y) = \sum_ {l = 0} ^ {\infty} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x))
$$

where the Gegenbauer coefficients ${ \hat { f } } ^ { \lambda } ( l )$ are given by (1.6). 

The candidate for approximating the Gegenbauer coefficients ${ \hat { f } } ^ { \lambda } ( l )$ is: 

$$
\hat {g} ^ {\lambda} (l) = \frac {1}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f _ {N} \circ x (y) C _ {l} ^ {\lambda} (y) d y. \tag {3.1}
$$

Definition 3.1 The truncation error is defined as 

$$
T E (\lambda , m, N) = \max _ {- 1 \leq y \leq 1} \left| \sum_ {l = 0} ^ {m} (\hat {f} ^ {\lambda} (l) - \hat {g} ^ {\lambda} (l)) C _ {l} ^ {\lambda} (y) \right| \tag {3.2}
$$

where, ${ \hat { f } } ^ { \lambda } ( l ) , { \hat { g } } ^ { \lambda } ( l )$ are defined in (1.6) and (3.1) 

In the next two lemmas, we bound the truncation error using the regularity of the function M(y) as in (3.3), and then in terms of the number of given Fourier coefficients N , the number of Gegenbauer polynomials m and parameter of Gegenbauer polynomial λ. 

Lemma 3.2 The truncation error is bounded by 

$$
T E (\lambda , m, N) \leq \frac {C}{(N \pi) ^ {t - 1}} \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \left\| \frac {d ^ {t}}{d x ^ {t}} M (y) \right\| _ {L ^ {\infty}}
$$

where $t = \lfloor \frac { \lambda + \frac { 1 } { 2 } } { q } \rfloor - 1$ λ+ 1 , and 

$$
M (y) = (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} C _ {l} ^ {\lambda} (y) \frac {d y}{d x} \tag {3.3}
$$

Proof We have 

$$
\begin{array}{l} T E (\lambda , m, N) = \max _ {- 1 \leq y \leq 1} \left| \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (y)}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} (f - f _ {N}) \circ x (y) C _ {l} ^ {\lambda} (y) d y \right| \\ \leq \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \left| \int_ {- 1} ^ {1} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} (f - f _ {N}) \circ x (y) C _ {l} ^ {\lambda} (y) d y \right| \\ = \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \left| \int_ {- 1} ^ {1} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \sum_ {| n | > N} \tilde {f} _ {n} e ^ {i n \pi x} C _ {l} ^ {\lambda} (y) d y \right| \\ = \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \left| \sum_ {| n | > N} \tilde {f} _ {n} \int_ {- 1} ^ {1} e ^ {i n \pi x} M (y) d x \right| \\ = \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \left| \sum_ {| n | > N} \frac {\tilde {f} _ {n}}{(i n \pi) ^ {t}} \int_ {- 1} ^ {1} e ^ {i n \pi x} \frac {d ^ {t}}{d x ^ {t}} M (y) d x \right| \\ \leq C \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \sum_ {| n | > N} \frac {1}{(| n | \pi) ^ {t}} \left\| \frac {d ^ {t}}{d x ^ {t}} M (y) \right\| _ {L ^ {\infty}} \\ \leq \frac {C}{(N \pi) ^ {t - 1}} \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \left\| \frac {d ^ {t}}{d x ^ {t}} M (y) \right\| _ {L ^ {\infty}} \\ \end{array}
$$

The definitions of ${ \hat { f } } ^ { \lambda } ( l )$ in (1.6) and $\hat { g } ^ { \lambda } ( l )$ in (3.1) are used in the first equality; (2.3) is used in the second inequality; in the third equality, the error of the Fourier partial sum $\begin{array} { r } { f ( x ) - f _ { N } ( x ) = \sum _ { | n | > N } \bar { f } _ { n } e ^ { i n \pi x } } \end{array}$ is used; the substitution (1.5) has been made in the integral in the fourth equality; in the fifth equality, we use integration by parts t times, and the fact that ${ \frac { d ^ { i } } { d x ^ { i } } } M ( y )$ vanishes at $y = \pm 1$ for $0 \leq i \leq t - 1$ ; since $f ( x )$ is an $L ^ { 2 } \mathrm { - f u n c t i o n }$ , its Fourier coefficients $\tilde { f _ { n } }$ are uniformly bounded, i.e., $| \tilde { f _ { n } } | \leq C$ , which is used in the sixth inequality. 

# Lemma 3.3

$$
T E (\lambda , m, N) \leq C \frac {(m + 1) 2 ^ {\lambda} A ^ {t + 1} (m + \lambda) ^ {m + t + 1} \Gamma (\lambda) \Gamma (m + 2 \lambda)}{m ! \Gamma (2 \lambda) \Gamma (m + \lambda + \frac {1}{2}) (N \pi) ^ {t - 1}}
$$

Proof From Lemma 3.2, we have 

$$
\begin{array}{l} T E (\lambda , m, N) \leq \frac {C}{(N \pi) ^ {t - 1}} \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \left\| \frac {d ^ {t}}{d x ^ {t}} \left[ (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} C _ {l} ^ {\lambda} (y) \frac {d y}{d x} \right] \right\| _ {L ^ {\infty}} \\ = \frac {C}{(N \pi) ^ {t - 1}} \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1) G (\lambda , l)}{2 ^ {l} l ! h _ {l} ^ {\lambda}} \left\| \frac {d ^ {t}}{d x ^ {t}} \left\{\frac {d ^ {l}}{d y ^ {l}} (1 - y ^ {2} (x)) ^ {l + \lambda - \frac {1}{2}} \frac {d y}{d x} \right\} \right\| _ {L ^ {\infty}} \\ \leq C \frac {A ^ {t + 1} 2 ^ {\lambda}}{(N \pi) ^ {t - 1}} \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1) G (\lambda , l) (l + \lambda) ^ {l + t}}{l ! h _ {l} ^ {\lambda}} \\ \end{array}
$$

$$
\begin{array}{l} = C \frac {A ^ {t + 1} 2 ^ {\lambda} \Gamma (\lambda)}{\Gamma (2 \lambda) (N \pi) ^ {t - 1}} \sum_ {l = 0} ^ {m} \frac {\Gamma (l + 2 \lambda) (l + \lambda) ^ {l + t + 1}}{l ! \Gamma (l + \lambda + \frac {1}{2})} \\ \leq C \frac {A ^ {t + 1} 2 ^ {\lambda} \Gamma (\lambda) (m + 1) \Gamma (m + 2 \lambda) (m + \lambda) ^ {m + t + 1}}{m ! \Gamma (2 \lambda) \Gamma (m + \lambda + \frac {1}{2}) (N \pi) ^ {t - 1}} \\ \end{array}
$$

where the definition of Gegenbauer polynomial in (2.1) is used in the second step; Lemma 2.6 and Remark 2.7 are used in the third inequality; (2.1) and (2.4) are used in fourth equality; in the last step, we use the fact that $\Gamma ( l + 2 \dot { \lambda } ) ( \dot { l } + \lambda ) ^ { l + t + 1 } / \{ l ! \Gamma ( l + \lambda + \frac { 1 } { 2 } ) \}$ is an increasing function of l. □ 

If both λ and m grow linearly with N , the truncation error can be made to decay exponentially. This result is stated in the following theorem. 

Theorem 3.4 (The exponential decay of the truncation error) $L e t \lambda = \alpha N , m = \beta N$ with $0 < \alpha , \beta < 1$ , then 

$$
T E (\alpha N, \beta N, N) \leq C N ^ {2} q _ {T} ^ {N}
$$

with 

$$
q _ {T} = \frac {A ^ {\frac {\alpha}{q}} e ^ {\beta} (2 \alpha + \beta) ^ {2 \alpha + \beta}}{\pi^ {\frac {\alpha}{q}} (2 \alpha) ^ {\alpha} \beta^ {\beta} (\alpha + \beta) ^ {(1 - \frac {1}{q}) \alpha}}.
$$

When we choose $\beta = \gamma \alpha ,$ , i.e. $m = \gamma \lambda$ , we have 

$$
q _ {T} = \left(\frac {A ^ {\frac {1}{q}} e ^ {\gamma} (2 + \gamma) ^ {2 + \gamma} \alpha^ {\frac {1}{q}}}{2 \pi^ {\frac {1}{q}} \gamma^ {\gamma} (1 + \gamma) ^ {1 - \frac {1}{q}}}\right) ^ {\alpha}.
$$

If we choose α to satisfy 

$$
\alpha <   \frac {\pi}{A} \left(\frac {2 \gamma^ {\gamma} (1 + \gamma) ^ {1 - \frac {1}{q}}}{e ^ {\gamma} (2 + \gamma) ^ {2 + \gamma}}\right) ^ {q}, \tag {3.4}
$$

then $q _ { T } < 1$ . 

Proof From lemma 3.3, we use Stirling’s formula (2.5) and get 

$$
\begin{array}{l} T E (\lambda , m, N) \leq C \frac {(m + 1) 2 ^ {\lambda} A ^ {t + 1} (m + \lambda) ^ {m + t + 1} \Gamma (\lambda) \Gamma (m + 2 \lambda)}{m ! \Gamma (2 \lambda) \Gamma (m + \lambda + \frac {1}{2}) (N \pi) ^ {t - 1}} \\ \leq C \frac {e ^ {m} A ^ {t + 1} (m + 2 \lambda) ^ {m + 2 \lambda - \frac {1}{2}}}{2 ^ {\lambda} \lambda^ {\lambda} (m + \lambda) ^ {\lambda - t - 1} m ^ {m - \frac {1}{2}} (N \pi) ^ {t - 1}} \\ \end{array}
$$

hence 

$$
\begin{array}{l} T E (\alpha N, \beta N, N) \leq C N ^ {2} \left(\frac {A (\alpha + \beta)}{\pi}\right) ^ {t} \left\{\frac {e ^ {\beta} (2 \alpha + \beta) ^ {2 \alpha + \beta}}{(2 \alpha) ^ {\alpha} \beta^ {\beta} (\alpha + \beta) ^ {\alpha}} \right\} ^ {N} \\ \leq C N ^ {2} \left(\frac {A (\alpha + \beta)}{\pi}\right) ^ {\frac {\alpha}{q} N} \left\{\frac {e ^ {\beta} (2 \alpha + \beta) ^ {2 \alpha + \beta}}{(2 \alpha) ^ {\alpha} \beta^ {\beta} (\alpha + \beta) ^ {\alpha}} \right\} ^ {N} \\ = C N ^ {2} q _ {T} ^ {N} \\ \end{array}
$$

where we have used $\begin{array} { r } { \frac { \lambda } { q } - 2 \leq t = \lfloor \frac { \lambda + \frac { 1 } { 2 } } { q } \rfloor - 1 \leq \frac { \lambda } { q } } \end{array}$ λ 1 in the second step. q 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/051a6ce5-4596-4416-8f09-ffb550db0f85/5fc2702f6f08f94db63734373972df0bbb4a8742bbbadb1d91fb1b818fb856dc.jpg)


Remark 3.5 With proper choice of the ratios between λ, m and N , $q _ { T }$ can be a positive number less than 1, and then the upper bound of truncation error is exponentially small. In the proof, there is no effort to optimize the parameters. In practice, λ and m should be chosen to balance the truncation error as well as the regularization error, which will be analyzed in the next section. 

# 4 Regularization Error

Regularization error, the second part of the post-processing error, is caused by approximating $f ( x ) = a ( x ) + b ( x ) ( 1 + x ) ^ { s } \mathrm { o n } [ - 1 , 1 ]$ by its Gegenbauer expansion based on the Gegenbauer polynomials $C _ { l } ^ { \lambda } ( x ( y ) )$ . This has been studied in [7]. We will just quote the result. 

The function $F ( y ) = f ( x ( y ) ) = a \circ x ( y ) + b \circ x ( y ) ( 1 + y ) ^ { p } 2 ^ { s - p }$ of the variable y is analytic on $y \in [ - 1$ , 1 . Let us consider the Gegenbauer partial sum of the first m terms for the function $f ( x )$ given by: 

$$
f ^ {\lambda , m} (x) = \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)), \tag {4.1}
$$

with the Gegenbauer coefficients ${ \hat { f } } ^ { \lambda } ( l )$ defined by (1.6). 

Definition 4.1 The regularization error is defined by 

$$
\begin{array}{l} R E (\lambda , m) = \max _ {- 1 \leq x \leq 1} \left| f (x) - \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)) \right| \\ = \max _ {- 1 \leq y \leq 1} \left| f (x (y)) - \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y) \right|. \\ \end{array}
$$

When both λ and m grow linearly with N , we have the following result for the estimate of regularization error [7]. 

Theorem 4.2 (The exponential decay of the regularization error) For the function $f ( x ) =$ $a ( x ) + b ( x ) ( 1 + x ) ^ { s }$ , with analytic functions a(x) and $b ( x )$ satisfying Assumption 1.1, if we assume $\lambda = \alpha N$ and $m = \gamma \lambda$ , then 

$$
\max _ {- 1 \leq x \leq 1} \left| f (x) - \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)) \right| \leq C q _ {R} ^ {N}
$$

where 

$$
q _ {R} = \left(\frac {(\gamma + 2) ^ {\gamma + 2}}{\rho^ {\gamma} 2 ^ {\gamma + 2} (\gamma + 1) ^ {\gamma + 1}}\right) ^ {\alpha}
$$

which is always less than 1. 

Proof The proof of this theorem can be found in [7]. 

□ 

We can combine the estimates for the truncation error and the regularization error and obtain the following theorem about the reconstruction error: 

Theorem 4.3 (Removal of the Gibbs Phenomenon) Consider a function in the form of $f ( x ) = a ( x ) + b ( x ) ( 1 + x ) ^ { s }$ , with given fractional constant $\textstyle 0 < s = { \frac { p } { q } } < 1$ , and a(x) and b(x) are analytic functions satisfying Assumption 1.1. Assume that the Fourier coefficient 

$$
\tilde {f} _ {n} = \frac {1}{2} \int_ {- 1} ^ {1} f (x) e ^ {- i k \pi x} d x
$$

are known for $- N \leq n \leq N .$ . Let $\hat { g } ^ { \lambda } ( l )$ , with $0 \leq l \leq m$ , be the Gegenbauer expansion coefficients of $\begin{array} { r } { f _ { N } ( x ) = \sum _ { | n | \leq N } \hat { f } _ { n } e ^ { i n \pi x } } \end{array}$ , defined in (3.1). Then for $\lambda = \alpha N$ and $m = \gamma \lambda$ with 

$$
\alpha <   \frac {\pi}{A} \left(\frac {2 \gamma^ {\gamma} (1 + \gamma) ^ {1 - \frac {1}{q}}}{e ^ {\gamma} (2 + \gamma) ^ {2 + \gamma}}\right) ^ {q},
$$

we have 

$$
\max _ {- 1 \leq x \leq 1} \left| f (x) - \sum_ {l = 0} ^ {m} \hat {g} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)) \right| \leq C \left(q _ {T} ^ {N} + q _ {R} ^ {N}\right)
$$

where 

$$
q _ {T} = \left(\frac {A ^ {\frac {1}{q}} e ^ {\gamma} (2 + \gamma) ^ {2 + \gamma} \alpha^ {\frac {1}{q}}}{2 \pi^ {\frac {1}{q}} \gamma^ {\gamma} (1 + \gamma) ^ {1 - \frac {1}{q}}}\right) ^ {\alpha} <   1, q _ {R} = \left(\frac {(\gamma + 2) ^ {\gamma + 2}}{\rho^ {\gamma} 2 ^ {\gamma + 2} (\gamma + 1) ^ {\gamma + 1}}\right) ^ {\alpha} <   1.
$$

Remark 4.4 In the proof, no attempt has been made to optimize the parameters. We will discuss about the choice of parameters in numerical examples in Sect. 6. 

Remark 4.5 For functions with multiple singularities, the reconstructions can be built piece by piece, and on each piece the function has at most one end-point singularity. The technique only needs to be modified by scaling. It is easy to prove that the reconstruction for each piece can also be exponentially accurate with proper parameters. 

# 5 Error Analysis for Fourier Galerkin Methods with Reconstructions

In this section, we will consider using the new reconstruction method to post-process the solution obtained by using Fourier Galerkin methods to solve linear, hyperbolic time-dependent partial differential equations (PDEs) [16]. Consider the problem 

$$
\left\{ \begin{array}{l l} u _ {t} - \mathscr {L} u = 0, & x \in [ - 1, 1 ], \quad t > 0 \\ u (x, 0) = g (x), & x \in [ - 1, 1 ] \end{array} \right. \tag {5.1}
$$

where $\mathcal { L }$ is a differential operator. In the Fourier Galerkin method, we seek solution $u _ { N } ( x , t )$ from the space $\hat { B } _ { N } = \mathrm { s p a n } \{ e ^ { i k \pi x } \} _ { | k | \leq N }$ , i.e.: 

$$
u _ {N} (x, t) = \sum_ {| k | \leq N} a _ {k} (t) e ^ {i k \pi x} \tag {5.2}
$$

where, $a _ { k } ( t )$ are unknown coefficients which will be determined by the method. In general, the coefficients $a _ { k } ( t )$ are not equal to the Fourier coefficients $\tilde { u } _ { k } .$ ; they will be equal only if we obtain the exact solution of the problem. In the Fourier Galerkin method, $a _ { k } ( t )$ are determined by requiring that the residual 

$$
R _ {N} (x, t) = \frac {\partial u _ {N} (x , t)}{\partial t} - \mathscr {L} u _ {N} (x, t) \tag {5.3}
$$

is orthogonal to ${ \hat { B } } _ { N }$ . The method is defined by the requirement that the orthogonal projection of the residual onto the space ${ \hat { B } } _ { N }$ is zero. If the residual is smooth enough, this requirement implies that the residual itself is small. 

When spectral method is used to solve PDEs, we will estimate the errors in the maximum norm. If the solution lacks regularity, 

$$
\| u (x, t) - u _ {N} (x, t) \| _ {L ^ {\infty}}
$$

can not be small even for $t = 0 ,$ , due to the Gibbs phenomenon. We would like to reconstruct the solution from $u _ { N } ( x , t )$ at the final time $t = T$ , and recover exponentially decaying error 

$$
\| u (x, T) - Q _ {N} u _ {N} (x, T) \| _ {L ^ {\infty}}
$$

where $Q _ { N }$ is a post-processing operator. 

In order to illustrate the post-processing methods clearly, we consider the simple linear transport equations involving piecewise smooth functions with unbounded derivative singularities: 

$$
\left\{ \begin{array}{l l} u _ {t} + c u _ {x} = 0, & x \in [ - 1, 1 ], \quad t > 0 \\ u (x, 0) = a (x) + b (x) (1 + x) ^ {s}, & x \in [ - 1, 1 ] \end{array} \right. \tag {5.4}
$$

with periodic boundary conditions, where $c$ is the phase speed (for simplicity, assume that $c > 0 )$ , s is a given fractional constant 

$$
0 <   s = \frac {p}{q} <   1 \tag {5.5}
$$

with relatively prime integers $p$ and $q ,$ , and $a ( x )$ and $b ( x )$ are both analytic functions. The initial condition is $C ^ { 0 , \alpha }$ Hölder continuous with $\alpha \leq s$ and 2-periodic. It has a singularity in 1, 1 , which means its first derivative blows up at $x = - 1$ . The solution to this problem is 

$$
u (x, t) = u (x - c t, 0) \tag {5.6}
$$

The singularities will move along the lines through points $( x , t ) = ( 2 n - 1 , 0 )$ (here, $n \in Z )$ with direction $( c , 1 )$ . Even though we consider functions with one singularity in this paper, our technique can be easily applied to functions with finitely many such singularities in 1, 1 of known locations. 

The coefficients $a _ { k } ( t )$ are determined by requiring that the residual 

$$
R _ {N} (x, t) = \frac {\partial u _ {N} (x , t)}{\partial t} + c \frac {\partial u _ {N} (x , t)}{\partial x} \tag {5.7}
$$

is orthogonal to ${ \hat { B } } _ { N }$ . If we rewrite the residual in terms of the Fourier series, 

$$
R _ {N} (x, t) = \sum_ {| k | \leq \infty} \hat {R} _ {k} (t) e ^ {i k \pi x}, \tag {5.8}
$$

the orthogonality requirement yields 

$$
\hat {R} _ {k} (t) = \frac {1}{2} \int_ {- 1} ^ {1} R _ {N} (x, t) e ^ {- i k \pi x} d x = 0, - N \leq k \leq N. \tag {5.9}
$$

This requirement provides $( 2 N + 1 )$ ordinary differential equations to determine the $( 2 N + 1 )$ unknowns $a _ { k } ( t )$ , 

$$
\frac {d}{d t} a _ {k} (t) + i k \pi c a _ {k} (t) = 0, \quad - N \leq k \leq N, \tag {5.10}
$$

and the corresponding initial conditions are 

$$
u _ {N} (x, 0) = \sum_ {| k | \leq N} a _ {k} (0) e ^ {i k \pi x}, \tag {5.11}
$$

$$
a _ {k} (0) = \frac {1}{2} \int_ {- 1} ^ {1} u (x, 0) e ^ {- i k \pi x} d x \tag {5.12}
$$

For this linear constant coefficient problem, $R _ { N } ( x , t )$ itself is in the space $\hat { B } _ { N }$ , therefore the orthogonal complement must be zero, i.e., 

$$
R _ {N} (x, t) = 0.
$$

It means that we obtain the exact solution of this linear problem by Fourier Galerkin methods. Therefore, the coefficients $a _ { k } ( t )$ are exactly the Fourier coefficients for the solution $u ( x , t )$ , which has singularities at $( 2 n - 1 + c t , t )$ (here, $n \in Z )$ . 

In order to recover the exponential accuracy of the solution, we need to take postprocessing on its translation $v ( x , t ) = u ( x + c t , t )$ , which has a singularity $x \ = \ - 1$ in the interval 1, 1 . The function v(x, t) behaves like $( 1 + x ) ^ { s }$ near the singularity. Under this translation, the Fourier coefficients $b _ { k } ( t )$ of $v ( x , t )$ are 

$$
b _ {k} (t) = e ^ {i k \pi c t} a _ {k} (t). \tag {5.13}
$$

From the truncated Fourier series $\begin{array} { r } { v _ { N } ( x , t ) = \sum _ { | k | \leq N } b _ { k } ( t ) e ^ { i k \pi x } } \end{array}$ , we get a reconstructed approximation 

$$
v _ {N} ^ {\lambda , m} (x, t) = \sum_ {l = 0} ^ {m} \hat {v} _ {l} ^ {\lambda} C _ {l} ^ {\lambda} (y (x)), \quad x \in [ - 1, 1 ] \tag {5.14}
$$

with 

$$
\hat {v} _ {l} ^ {\lambda} = \frac {1}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} v _ {N} (x (y), t) C _ {l} ^ {\lambda} (y) d y. \tag {5.15}
$$

Then the numerical solution after processing $u _ { N } ^ { \lambda , m } ( x , t )$ 

$$
u _ {N} ^ {\lambda , m} (x, t) = Q _ {N} u _ {N} (x, t) = v _ {N} ^ {\lambda , m} (x - c t, t). \tag {5.16}
$$

From Theorem 4.3, we obtain the following theorem on the accuracy in spectral methods solving the problem (5.4). 

Theorem 5.1 Using Fourier Galerkin method to solve problem (5.4), the solution $u _ { N } ( x , t )$ is defined by 

$$
u _ {N} (x, t) = \sum_ {| k | \leq N} a _ {k} (t) e ^ {i k \pi x},
$$

where the coefficients $a _ { k } ( t )$ are determined by (5.10) and (5.12). After post-processing, we get 

$$
Q _ {N} u _ {N} (x, t) = v _ {N} ^ {\lambda , m} (x - c t, t)
$$

where, v $v _ { N } ^ { \lambda , m } ( x , t )$ λ,m is defined in (5.14) and (5.15). The error 

$$
\max _ {- 1 \leq x \leq 1} | u (x, t) - Q _ {N} u _ {N} (x, t) |
$$

is exponentially small, when the parameters λ and m chosen in the post-processing are proportional to N and satisfy (3.4). 

Remark 5.2 For linear problems (5.4) with variable coefficients c(x), the computed coefficients $a _ { k } ( t )$ are not exact for $t > 0$ . But they are still computed with exponential accuracy, even with nonsmooth $L ^ { 2 }$ initial data. Therefore, we can still observe improved accuracy after post-processing, and numerical tests are shown in Sect. 6. 

# 6 Numerical Results

In this section, we give numerical examples to illustrate the results of the reconstruction technique. In order to testify that exponential accuracy could be recovered from the knowledge of the first $2 N + 1$ Fourier coefficients, we first test on the same two examples as in [7]. 

Example 6.1 We take the first 2N + 1 Fourier coefficients of function 

$$
u (x) = \cos (x) + \sin (x) \sqrt {1 + x}, \quad x \in [ - 1, 1 ] \tag {6.1}
$$

and try to recover the pointwise values of the solution over 1, 1 . 

First, the parameters are chosen according to 

$$
\lambda = \frac {1}{1 6} N, \quad m = \frac {3}{8 0} N. \tag {6.2}
$$

For the sake of easy computation, we choose λ and m to be the biggest integer satisfying (6.2). 

In Fig. 1a, we show the errors in the logarithm scale, for N 40, 80, 160, 320 and 640. In Table 1 “linear choice”, we show the maximum errors of the reconstruction over x 1, 1 , for each N and the orders of convergence. We can clearly find out that the error is exponentially decaying. 

The other way is to find the optimal choice of the parameters such that the reconstruction achieves the smallest error. To choose the optimal parameters, we compare the errors of all the possible pairs $\begin{array} { r } { 0 \leq \lambda \leq \frac { N } { 2 } } \end{array}$ and $\begin{array} { r } { 0 \leq m \leq \frac { N } { 2 } } \end{array}$ . In Table 1 “optimal choice”, we also show the $L ^ { \infty }$ errors of the reconstruction with the optimal choice for each N . 

Example 6.2 We take the first 2N 1 Fourier coefficients of function 

$$
u (x) = \cos (x) + \sin (x) (1 + x) ^ {\frac {1}{3}}, \quad x \in [ - 1, 1 ] \tag {6.3}
$$

and try to recover the pointwise values of the solution over [−1, 1]. 

First, the parameters are chosen linearly dependent on N : 

$$
\lambda = \frac {1}{8} N, \quad m = \frac {1}{3 2} N. \tag {6.4}
$$

In Fig. 1b, we show the errors in the logarithm scale, for N 40, 80, 160, 320 and 640. In Table 2, we show the maximum errors of the reconstruction over $x \in [ - 1,1]$ ,  with both ways to choose parameters, for each N . We again clearly find out that the error is exponentially decaying with linear choice. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/051a6ce5-4596-4416-8f09-ffb550db0f85/bdfa0b669cb6e4dae3425291b1ff38a61b6946e630a4bb7d186717afdf2efc11.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/051a6ce5-4596-4416-8f09-ffb550db0f85/5f76eed0a13bf56583f158a1a0850b31500c1a4ca6a5d23af600597e1812eea1.jpg)



Fig. 1 Pointwise errors in the logarithm scale. a Example 1 $\begin{array} { r } { \left( \lambda = \frac { 1 } { 1 6 } N , m = \frac { 3 } { 8 0 } N \right) } \end{array}$ . b Example 2 $\begin{array} { r } { \left( \lambda = \frac { 1 } { 8 } N , m = \frac { 1 } { 3 2 } N \right) } \end{array}$


Next, we reconstruct the same two functions with different levels of noise in data to test the robustness to noise. The way we add noise to Fourier coefficients is to generate two random and independent sequences, which are uniformly distributed in $[ - \frac { \delta } { 2 } , \frac { \delta } { 2 } ]$ with noise level δ, and then add the sequences to the real parts and imaginary parts of the Fourier coefficients separately. 

Example 6.3 (Noise) For the function 

$$
u (x) = \cos (x) + \sin (x) \sqrt {1 + x}, \quad x \in [ - 1, 1 ] \tag {6.5}
$$

we try to recover it with given information of the first 2N 1 Fourier coefficients with certain level of noise. We compare the maximum error of the reconstruction with different levels of noise, for either the linear choice or the optimal choice of parameters λ and m as in Example 6.1. The results are listed in Tables 3 and 4. 


Table 1 Maximum error table


<table><tr><td>N</td><td>Linear choice<eq>L^{\infty}</eq> error</td><td>Order</td><td><eq>\lambda</eq></td><td>m</td><td>Optimal choice<eq>L^{\infty}</eq> error</td><td><eq>\lambda</eq></td><td>m</td></tr><tr><td>40</td><td>5.91E-001</td><td></td><td>2</td><td>1</td><td>1.21E-003</td><td>8</td><td>18</td></tr><tr><td>80</td><td>2.95E-001</td><td>1.00</td><td>5</td><td>3</td><td>5.25E-005</td><td>12</td><td>40</td></tr><tr><td>160</td><td>1.56E-002</td><td>4.25</td><td>10</td><td>6</td><td>5.57E-010</td><td>20</td><td>80</td></tr><tr><td>320</td><td>1.33E-005</td><td>10.19</td><td>20</td><td>12</td><td>2.75E-011</td><td>22</td><td>51</td></tr><tr><td>640</td><td>6.46E-012</td><td>20.98</td><td>40</td><td>24</td><td>1.11E-012</td><td>23</td><td>44</td></tr></table>


Table 2 Maximum error table


<table><tr><td>N</td><td>Linear choice<eq>L^{\infty}</eq> error</td><td>Order</td><td><eq>\lambda</eq></td><td>m</td><td>Optimal choice<eq>L^{\infty}</eq> error</td><td><eq>\lambda</eq></td><td>m</td></tr><tr><td>40</td><td>7.13E-001</td><td></td><td>5</td><td>1</td><td>1.83E-001</td><td>5</td><td>12</td></tr><tr><td>80</td><td>3.64E-001</td><td>0.97</td><td>10</td><td>2</td><td>1.98E-002</td><td>7</td><td>33</td></tr><tr><td>160</td><td>2.91E-001</td><td>0.33</td><td>20</td><td>5</td><td>1.24E-003</td><td>12</td><td>76</td></tr><tr><td>320</td><td>1.17E-002</td><td>4.64</td><td>40</td><td>10</td><td>5.50E-006</td><td>17</td><td>97</td></tr><tr><td>640</td><td>1.83E-006</td><td>12.64</td><td>80</td><td>20</td><td>3.15E-007</td><td>21</td><td>89</td></tr></table>


Table 3 Maximum error with different levels of noise in data (linear choice)


<table><tr><td>N</td><td><eq>L^{\infty}</eq> (no noise)</td><td><eq>L^{\infty}</eq></td><td>Noise</td><td><eq>L^{\infty}</eq></td><td>Noise</td><td><eq>L^{\infty}</eq></td><td>Noise</td></tr><tr><td>40</td><td>5.91E-01</td><td>6.07E-01</td><td>1.00E-02</td><td>5.53E-01</td><td>1.00E-01</td><td>8.54E-01</td><td>1.00E+00</td></tr><tr><td>80</td><td>2.95E-01</td><td>1.82E-01</td><td>1.00E-02</td><td>3.25E-01</td><td>1.00E-01</td><td>8.59E+00</td><td>1.00E+00</td></tr><tr><td>160</td><td>1.56E-02</td><td>1.53E-02</td><td>1.00E-05</td><td>8.16E-03</td><td>1.00E-04</td><td>9.44E-02</td><td>1.00E-03</td></tr><tr><td>320</td><td>1.33E-05</td><td>1.36E-05</td><td>1.00E-11</td><td>1.12E-05</td><td>1.00E-10</td><td>8.37E-05</td><td>1.00E-09</td></tr></table>


Table 4 Maximum error with different levels of noise in data (optimal choice)


<table><tr><td>N</td><td><eq>L^{\infty}</eq> (no noise)</td><td><eq>L^{\infty}</eq></td><td>Noise</td><td><eq>L^{\infty}</eq></td><td>Noise</td><td><eq>L^{\infty}</eq></td><td>Noise</td></tr><tr><td>40</td><td>1.21E-03</td><td>1.16E-03</td><td>1.00E-07</td><td>3.30E-03</td><td>1.00E-06</td><td>2.66E-02</td><td>1.00E-05</td></tr><tr><td>80</td><td>5.25E-05</td><td>4.74E-05</td><td>1.00E-12</td><td>1.20E-04</td><td>1.00E-11</td><td>4.19E-04</td><td>1.00E-10</td></tr><tr><td>160</td><td>5.57E-10</td><td>4.97E-10</td><td>1.00E-22</td><td>4.81E-09</td><td>1.00E-21</td><td>6.24E-08</td><td>1.00E-20</td></tr><tr><td>320</td><td>2.75E-11</td><td>2.38E-11</td><td>1.00E-23</td><td>3.25E-11</td><td>1.00E-22</td><td>1.76E-10</td><td>1.00E-21</td></tr></table>

Example 6.4 (Noise) For the function 

$$
u (x) = \cos (x) + \sin (x) (1 + x) ^ {\frac {1}{3}}, \quad x \in [ - 1, 1 ] \tag {6.6}
$$

we try to recover it with given information of the first 2N 1 Fourier coefficients with certain level of noise. We compare the maximum error of the reconstruction with different levels of noise, for either the linear choice or the optimal choice of parameters λ and m as in Example 6.2. The results are listed in Tables 5 and 6. 


Table 5 Maximum error with different levels of noise in data (linear choice)


<table><tr><td>N</td><td><eq>L^{\infty}</eq> (no noise)</td><td><eq>L^{\infty}</eq></td><td>Noise</td><td><eq>L^{\infty}</eq></td><td>Noise</td><td><eq>L^{\infty}</eq></td><td>Noise</td></tr><tr><td>40</td><td>7.13E-01</td><td>6.82E-01</td><td>1.00E-02</td><td>7.21E-01</td><td>1.00E-01</td><td>1.25E+00</td><td>1.00E+00</td></tr><tr><td>80</td><td>3.64E-01</td><td>3.43E-01</td><td>1.00E-02</td><td>8.04E-01</td><td>1.00E-01</td><td>8.30E+00</td><td>1.00E+00</td></tr><tr><td>160</td><td>2.91E-01</td><td>2.81E-01</td><td>1.00E-04</td><td>8.87E-01</td><td>1.00E-03</td><td>6.06E+00</td><td>1.00E-02</td></tr><tr><td>320</td><td>1.17E-02</td><td>1.65E-02</td><td>1.00E-08</td><td>4.93E-02</td><td>1.00E-07</td><td>1.13E+00</td><td>1.00E-06</td></tr></table>


Table 6 Maximum error with different levels of noise in data (optimal choice)


<table><tr><td>N</td><td><eq>L^{\infty}</eq> (no noise)</td><td><eq>L^{\infty}</eq></td><td>Noise</td><td><eq>L^{\infty}</eq></td><td>Noise</td><td><eq>L^{\infty}</eq></td><td>Noise</td></tr><tr><td>40</td><td>1.83E-01</td><td>1.80E-01</td><td>1.00E-04</td><td>1.77E-01</td><td>1.00E-03</td><td>2.14E+00</td><td>1.00E-02</td></tr><tr><td>80</td><td>1.98E-02</td><td>1.01E-02</td><td>1.00E-06</td><td>8.46E-02</td><td>1.00E-05</td><td>4.89E-01</td><td>1.00E-04</td></tr><tr><td>160</td><td>1.24E-03</td><td>1.21E-03</td><td>1.00E-12</td><td>7.29E-03</td><td>1.00E-11</td><td>1.53E-02</td><td>1.00E-10</td></tr><tr><td>320</td><td>5.50E-06</td><td>5.16E-06</td><td>1.00E-18</td><td>2.73E-06</td><td>1.00E-17</td><td>1.13E-04</td><td>1.00E-16</td></tr></table>

From the results in Example 6.3 and Example 6.4, we observe that the reconstruction method works pretty well with noise. The more accurate the reconstruction is, the more sensitive to noise it is. The linear choice of parameters is less sensitive, since it is patterned. The optimal one is numerically chosen for best errors, thus is more sensitive to any perturbation. And the allowable noise level to actual error level is about $e r r o r ^ { 2 }$ for Example 6.3 and $e r r o r ^ { 3 }$ for Example 6.4. 

At last, we have numerical tests on Fourier Galerkin methods with post-processing for initial boundary problems of linear transport equation with variable coefficients. 

Example 6.5 We solve the following initial boundary problem 

$$
\left\{ \begin{array}{l l} u _ {t} - x u _ {x} = 0, & x \in [ - 1, 1 ], \quad t > 0 \\ u (x, 0) = g (x) = \sqrt {1 + x}, & x \in [ - 1, 1 ] \end{array} \right. \tag {6.7}
$$

with periodic boundary condition. 

The exact solution is $u ( x , t ) = \sqrt { 1 + x e ^ { t } \ ( \mathrm { m o d } \ 2 ) }$ . If we look at the solution at T log 2, it has two deformed “copies” of the initial condition, and has singularities at two places 0.5 and 0.5. Therefore, the filtered Fourier Galerkin method with Runge–Kutta time stepping (as in [18]) provides solutions with poor accuracy and oscillations near the singularities, as in Fig. 2. Here, we would like to recover the solution on [−0.5, 0], since the solution behaves like square root at the left end of this sub-interval. 

$$
u (x, \log 2) = \sqrt {1 + 2 x}, \quad x \in [ - 0. 5, 0 ].
$$

We choose the parameters λ and m in two different ways. First, we choose λ and m linearly dependent on N : 

$$
\lambda = \frac {N}{8 0}, \quad m = \frac {N}{4}.
$$

The other way is to choose the pair, which gives best recovery, in the range $\begin{array} { r } { 0 \leq \lambda , m \leq \frac { N } { 2 } } \end{array}$ . The results are listed in Table 7 and shown in Fig. 3. With both choices, the solutions after post-processing provide improved accuracy, even up to the singularities. We also observe exponential accuracy from the reconstruction with the linear choice. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/051a6ce5-4596-4416-8f09-ffb550db0f85/d9c042217e22acfa7118e105a911dfd40b0355e9a372c3989dae0393f262b1da.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/051a6ce5-4596-4416-8f09-ffb550db0f85/1700f3d97c1f92bf95c38557ff654cec99341b0bb5c74265d0c776d8b3347747.jpg)



Fig. 2 Numerical solution u (x) (N = 10, 20, 40, and 80) on 0.5, 0 , using Fourier Galerkin method with filters. a Numerical solutions u (x) and exact solution u(x). b Pointwise errors in logarithm scale


# 7 Concluding Remarks

We have built the Gegenbauer polynomial based technique to reconstruct approximations with exponential accuracy in the maximum norm, based on the Fourier coefficients, for piecewise analytic functions where in each piece the function is analytic only in the open interval with unbounded derivative end-point singularities. Such functions are significantly more difficult to handle than piecewise analytic functions which are analytic up to the discontinuity points. Our technique provides a post-processing method for spectral methods to solve transport 


Table 7 Maximum error table


<table><tr><td>N</td><td>Linear choice<eq>L^{\infty}</eq> error</td><td>Order</td><td><eq>\lambda</eq></td><td>m</td><td>Optimal choice<eq>L^{\infty}</eq> error</td><td><eq>\lambda</eq></td><td>m</td><td>Without PP<eq>L^{\infty}</eq> error</td></tr><tr><td>10</td><td>6.93E-01</td><td></td><td>0</td><td>3</td><td>5.81E-01</td><td>1</td><td>2</td><td>6.57E-01</td></tr><tr><td>20</td><td>5.13E-01</td><td>0.43</td><td>0</td><td>5</td><td>1.72E-01</td><td>1</td><td>8</td><td>5.04E-01</td></tr><tr><td>40</td><td>1.01E-02</td><td>5.67</td><td>1</td><td>10</td><td>7.34E-03</td><td>1</td><td>6</td><td>4.45E-01</td></tr><tr><td>80</td><td>2.61E-05</td><td>8.60</td><td>1</td><td>20</td><td>2.61E-05</td><td>1</td><td>20</td><td>4.07E-01</td></tr></table>

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/051a6ce5-4596-4416-8f09-ffb550db0f85/3de5682840b852c837ae27c5c60ea0248465e39d3dd506af068438ad47783c13.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/051a6ce5-4596-4416-8f09-ffb550db0f85/b52ad41aadd84f7946b2ae582e672268b4e7da5aaaec8362d8fc28e72b84f803.jpg)



Fig. 3 Pointwise errors of the reconstructions with two different strategies for choosing parameters λ and m, in logarithm scale. a Linear choises. b Optimal pairs


equations involving such functions to achieve exponential accuracy. Preliminary numerical results are provided to demonstrate the theory. This technique is easy to be generated to higher dimensional cases, and has potential applications for solving other partial differential equations whose solutions have such singularities. 

# References



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

