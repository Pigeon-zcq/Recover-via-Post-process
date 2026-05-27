# Recovering exponential accuracy from collocation point values of smooth functions with end-point singularities✩

Zheng Chen, Chi-Wang Shu ∗ 

Division of Applied Mathematics, Brown University, Providence, RI 02912, United States 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/21786fb1-7f2e-43f0-b753-a48f89c56ea3/6f76c7803e5c61293cf7c364c362db353ec1720dd4f13268816b0397791a460f.jpg)


CrossMark 

# a r t i c l e i n f o

Article history: 

Received 18 November 2012 

Received in revised form 5 April 2013 

Dedicated to Professor Ben-yu Guo on his seventieth birthday 

Keywords: 

Spectral approximation 

Exponential accuracy 

Gegenbauer expansion 

Collocation 

Gaussian points 

# a b s t r a c t

Gibbs phenomenon is the particular manner how a global spectral approximation of a piecewise analytic function behaves at the jump discontinuity. The truncated spectral series has large oscillations near the jump, and the overshoot does not decay as the number of terms in the truncated series increases. There is therefore no convergence in the maximum norm, and convergence in smooth regions away from the discontinuity is also slow. In Gottlieb and Shu (1995) [5], a methodology is proposed to completely overcome this difficulty in the context of spectral collocation methods, resulting in the recovery of exponential accuracy from collocation point values of a piecewise analytic function. In this paper, we extend this methodology to handle spectral collocation methods for functions which are analytic in the open interval but have singularities at end-points. With this extension, we are able to obtain exponential accuracy from collocation point values of such functions. Similar to Gottlieb and Shu (1995) [5], the proof is constructive and uses the Gegenbauer polynomials $C _ { n } ^ { \lambda } ( x )$ . The result implies that the Gibbs phenomenon can be overcome for smooth functions with endpoint singularities. 

© 2013 Elsevier B.V. All rights reserved. 

# 1. Introduction

In this paper, we investigate the issue of overcoming the Gibbs phenomenon. Gibbs phenomenon is the particular manner how a global spectral approximation of a piecewise analytic function behaves at the jump discontinuity. A prototype is the Fourier approximation of an analytic but non-periodic function, which, when viewed as a periodic function, has a discontinuity at the boundaries of the interval. The truncated Fourier series has large oscillations near this jump, and the overshoot does not decay as the number of terms retained in the series increases. There is therefore no convergence in the maximum norm. Also, convergence in smooth regions away from the discontinuity is slow. In a series of papers [1–5], Gottlieb et al. developed a general framework to overcome this difficulty, in the sense that exponential accuracy is recovered in the maximum norm for any sub-interval (including the whole interval if there is only one discontinuity which is located at the interval boundaries), from the knowledge of either the first N spectral expansion coefficients, or the point values at N standard collocation points. This means that exponential accuracy is recovered at all points, including at the actual discontinuity points (the left and right limits at these points), if the locations of these discontinuity points are known. If the locations of these discontinuity points are not known exactly but are known to be within certain fixed intervals, then exponential accuracy can be recovered from any interval which does not overlap with these fixed intervals containing the discontinuities. An important tool used in this framework is the set of Gegenbauer polynomials, which are orthogonal in the interval [−1, 1] with the weight $( 1 - x ^ { 2 } ) ^ { \lambda - \frac { 1 } { 2 } }$ . It turns out that the main ingredient in this technique is that the parameter λ in the weight function as well as the number of terms m retained in the Gegenbauer expansion should both be chosen proportional to N. For a review of this series of work, we refer to [6]. 

In this study, we are interested in collocation methods. In [5], a methodology is proposed to completely overcome the difficulty of Gibbs phenomenon, in the context of spectral collocation methods, resulting in the recovery of exponential accuracy from collocation point values of a piecewise analytic function. In this paper, we extend this methodology to handle spectral collocation methods for functions which are analytic in the open interval but have singularities at end-points. Such functions appear often in applications, for example, as solutions of certain partial differential equations (PDEs). We assume that we are given the point values {f (xi)}, where $x _ { i }$ are the standard collocation points (Gaussian points of the orthogonal basis $\{ \phi _ { k } ( x ) \}$ , which may be the Fourier trigonometric polynomials, the Legendre polynomials, the Chebyshev polynomials, or the general Gegenbauer polynomials). The objective is to recover exponentially accurate point values at every point including at the singularities. 

The reconstruction procedure is performed on functions of the following form 

$$
f (x) = a (x) + b (x) (1 + x) ^ {s}, \quad x \in [ - 1, 1 ] \tag {1.1}
$$

where s is a given fractional constant 

$$
0 <   s = \frac {p}{q} <   1, \tag {1.2}
$$

in which p and q are relatively prime positive integers, and a(x) and $b ( \boldsymbol { x } )$ are both analytic functions. Unfortunately, our current approach cannot treat functions of the form (1.1) with an irrational power s. 

Since we assume that the value of s is known, we are not losing generality by assuming (1.1) with $0 < s < 1$ . All functions with end-point singularity at the left end, of the form (1.1) with an arbitrary positive or negative $s ,$ can be rewritten to the same form with $0 < s < 1$ by multiplying with an analytic function $( 1 + x ) ^ { k }$ with a positive integer k, or by absorbing an analytic function $( 1 + x ) ^ { k }$ with a positive integer k into b(x). Since we are handling collocation methods and are given values of the original functions at the collocation points, we also have access to the collocation point values of the modified functions which are the original functions multiplied with $( 1 + x ) ^ { k }$ . Singularity at the right end or at both ends can be handled in the same fashion. Also, the result can be easily generalized to the situation of finitely many singularities (of the form $( x - z _ { k } ) ^ { s _ { k } }$ at finitely many points $z _ { k }$ inside the interval [−1, 1] with fractional $s _ { k } ) ,$ , along the lines of [3,4] and using the techniques in this paper. 

As in [1], we assume that the analytic functions a(x) and b(x), denoted generically as $c ( \boldsymbol x )$ , satisfy the following condition. 

Assumption 1.1. There exists a constant $\rho \ge 1$ and a constant $C ( \rho )$ such that, for every $k \geq 0$ , 

$$
\max _ {- 1 \leq x \leq 1} \left| \frac {d ^ {k} c (x)}{d x ^ {k}} \right| \leq C (\rho) \frac {k !}{\rho^ {k}}.
$$

This is a standard assumption for analytic functions, where $\rho$ is the distance from the interval $[ - 1 , 1 ]$ to the nearest singularity of the function c(x) in the complex plane. 

We will use the following one to one transformation between $x \in [ - 1 , 1 ]$ and $y \in [ - 1 , 1 ]$ : 

$$
(2 ^ {q - 1} (1 + x)) ^ {\frac {1}{q}} = 1 + y
$$

where q is defined in (1.2) in which p and q are assumed to be relatively prime. 

The function ${ \bar { f } } ( y ) = f ( x ( y ) )$ of the variable y has its usual Gegenbauer expansion under the basis $\{ C _ { l } ^ { \lambda } ( y ) \}$ : 

$$
f (x (y)) = \bar {f} (y) = \sum_ {l = 0} ^ {\infty} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y)
$$

with the Gegenbauer coefficients $\hat { f } ^ { \lambda } ( l )$ given by 

$$
\hat {f} ^ {\lambda} (l) = \frac {1}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \bar {f} (y) C _ {l} ^ {\lambda} (y) d y \tag {1.3}
$$

where the precise value of the normalization constant $h _ { \scriptscriptstyle { l } } ^ { \lambda }$ will be given by (2.3) in Definition 2.1. 

Our goal is to find a good approximation to the first m ∼ N Gegenbauer coefficients $\hat { f } ^ { \lambda } ( l )$ in (1.3), denoted as $\hat { g } ^ { \lambda } ( l )$ , from the known point values $\{ f ( x _ { i } ) \}$ at the standard Gaussian collocation points. We will then obtain the approximation of f (x) using these $m \sim N$ terms of its Gegenbauer expansion: 

$$
f _ {N} ^ {m, \lambda} (x) = \sum_ {l = 0} ^ {m} \hat {g} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)).
$$

As in [5], we will separate the analysis of the error into two parts: the truncation error and the regularization error. Unlike [5], the analysis for the truncation error must be completely changed because of the end-point singularity. The truncation error measures the difference between the exact Gegenbauer coefficients of $f ( x ( y ) )$ with $\lambda \sim N$ , and the approximate 

Gegenbauer coefficients $\hat { g } ^ { \lambda } ( l )$ obtained by using the collocation point values. This will be investigated in Section 3. The regularization error measures the difference between the Gegenbauer expansion with $\lambda \sim N$ , using the first m ∼ N Gegenbauer coefficients, and the function itself. This error is estimated in Section 4. The results are summarized in Theorem 4.3 in Section 4. Section 5 contains several numerical examples to illustrate our results. In Section 2, we shall give several useful preliminary properties, estimates, and approximation results. Concluding remarks are given in Section 6. 

Throughout this paper, we will use C to denote a generic constant either independent of the growing parameters, or depending on them at most in polynomial growth. The details will be indicated clearly in the text. These constants may not take the same value at different places. 

# 2. Preliminaries

In this section, we will first introduce the Gegenbauer polynomials and discuss some of their asymptotic behavior. We rely heavily on the standardization in Bateman [7]. 

Definition 2.1. The Gegenbauer polynomial $C _ { n } ^ { \lambda } ( x ) , \mathrm { f o r } \lambda \ge 0$ , is defined by 

$$
(1 - x ^ {2}) ^ {\lambda - \frac {1}{2}} C _ {n} ^ {\lambda} (x) = \frac {(- 1) ^ {n}}{2 ^ {n} n !} G (\lambda , n) \frac {d ^ {n}}{d x ^ {n}} \left[ (1 - x ^ {2}) ^ {n + \lambda - \frac {1}{2}} \right]
$$

where G(λ, n) is given by 

$$
G (\lambda , n) = \frac {\Gamma \left(\lambda + \frac {1}{2}\right) \Gamma (n + 2 \lambda)}{\Gamma (2 \lambda) \Gamma \left(n + \lambda + \frac {1}{2}\right)}
$$

for λ > 0, by 

$$
G (0, n) = \frac {2 \sqrt {\pi} (n - 1) !}{\Gamma \left(n + \frac {1}{2}\right)}
$$

for λ = 0 and n ≥ 1, and by 

$$
G (0, 0) = 1
$$

for λ = 0 and n = 0. Notice that by this standardization, $C _ { n } ^ { 0 } ( x )$ is defined by (see [7]): 

$$
C _ {n} ^ {0} (x) = \lim _ {\lambda \rightarrow 0 ^ {+}} \frac {1}{\lambda} C _ {n} ^ {\lambda} (x) = \frac {2}{n} T _ {n} (x), \quad n > 0; \qquad C _ {0} ^ {0} (x) = 1,
$$

where $T _ { n } ( x )$ are the Chebyshev polynomials. 

Under this definition we have, for λ > 0, 

$$
C _ {n} ^ {\lambda} (1) = \frac {\Gamma (n + 2 \lambda)}{n ! \Gamma (2 \lambda)}; \tag {2.1}
$$

for λ = 0 and n ≥ 1, 

$$
C _ {n} ^ {0} (1) = \frac {2}{n};
$$

for λ = 0 and n = 0, 

$$
C _ {n} ^ {0} (1) = 1;
$$

and 

$$
\left| C _ {n} ^ {\lambda} (x) \right| \leq C _ {n} ^ {\lambda} (1), \quad - 1 \leq x \leq 1. \tag {2.2}
$$

The Gegenbauer polynomials are orthogonal under their weight function $( 1 - x ^ { 2 } ) ^ { \lambda - \frac { 1 } { 2 } }$ : 

$$
\int_ {- 1} ^ {1} (1 - x ^ {2}) ^ {\lambda - \frac {1}{2}} C _ {k} ^ {\lambda} (x) C _ {n} ^ {\lambda} (x) d x = \delta_ {k, n} h _ {n} ^ {\lambda}
$$

where, for λ > 0, 

$$
h _ {n} ^ {\lambda} = \pi^ {\frac {1}{2}} C _ {n} ^ {\lambda} (1) \frac {\Gamma \left(\lambda + \frac {1}{2}\right)}{\Gamma (\lambda) (n + \lambda)}; \tag {2.3}
$$

for λ = 0 and n ≥ 1, 

$$
h _ {n} ^ {0} = \frac {2 \pi}{n ^ {2}};
$$

for λ = 0 and n = 0, 

$$
h _ {0} ^ {0} = \pi .
$$

We will need to use the Stirling’s formula and the estimate of $h _ { n } ^ { \lambda }$ for the asymptotics of the Gegenbauer polynomials for large n and λ. 

Lemma 2.2. We have the Stirling’s formula 

$$
(2 \pi) ^ {\frac {1}{2}} x ^ {x + \frac {1}{2}} e ^ {- x} \leq \Gamma (x + 1) \leq (2 \pi) ^ {\frac {1}{2}} x ^ {x + \frac {1}{2}} e ^ {- x + \frac {1}{1 2 x}}, \quad x \geq 1. \tag {2.4}
$$

Lemma 2.3. There exists a constant C independent of λ and n such that 

$$
C ^ {- 1} \frac {\lambda^ {\frac {1}{2}}}{n + \lambda} C _ {n} ^ {\lambda} (1) \leq h _ {n} ^ {\lambda} \leq C \frac {\lambda^ {\frac {1}{2}}}{n + \lambda} C _ {n} ^ {\lambda} (1). \tag {2.5}
$$

We also need to quote the approximation results for the interpolation polynomials. Given the point values $\{ f ( x _ { i } ) \}$ of the function $f ( x )$ at the N Gaussian points $\left\{ x _ { i } \right\}$ of the basis functions $\{ \phi _ { k } ( x ) \}$ , we denote the unique interpolation polynomial by 

$$
I _ {N} f (x) = \sum_ {k = 1} ^ {N} \hat {f} _ {k} \phi_ {k} (x) \tag {2.6}
$$

satisfying 

$$
I _ {N} f (x _ {i}) = f (x _ {i}), \quad i = 1, 2, \dots , N. \tag {2.7}
$$

For the interpolation polynomial $I _ { N } f ( x )$ , we have the following error estimates: 

Lemma 2.4. $I f f ( x )$ has m continuous derivatives in [−1, 1], then we have the following estimate for the interpolation polynomial $I _ { N } f ( x )$ defined by (2.6) and (2.7), where the basis functions $\{ \phi _ { k } ( x ) \}$ } are either trigonometric polynomials $\{ e ^ { i k \pi x } \}$ (with $k =$ $- N , \ldots , N )$ or the Gegenbauer polynomials $\{ C _ { k } ^ { \mu } ( x ) \}$ with $\mu > - \frac { 1 } { 2 }$ (including the Chebyshev case for $\mu = 0$ and the Legendre case for $\begin{array} { r } { \mu = \frac { 1 } { 2 } ) } \end{array}$ 

$$
\| f - I _ {N} f \| _ {L _ {\omega} ^ {2}} \leq \frac {C}{N ^ {m}} \| f ^ {(m)} \| _ {L ^ {\infty}},
$$

where the weighted $L ^ { 2 }$ norm is defined as 

$$
\| f \| _ {L _ {\omega} ^ {2}} ^ {2} = \int_ {- 1} ^ {1} \omega (x) | f (x) | ^ {2} d x
$$

with the weight function ω $\ b ( x ) = 1 f o r$ the trigonometric or Legendre polynomials bases and $\omega ( x ) = ( 1 - x ^ { 2 } ) ^ { \mu - 1 / 2 }$ for the general Gegenbauer polynomial bases. Here C is a constant independent of N and m. 

Proof. The proof for the trigonometric polynomial case can be found in [8]. For the Gegenbauer polynomial case, the proof can be found in [9]. Although the authors do not explicitly point out, the constant C is independent of m in both situations.  

We would need to estimate $\begin{array} { r } { \parallel \frac { d ^ { t } } { d x ^ { t } } \{ ( 1 - y ^ { 2 } ) ^ { \lambda - \frac { 1 } { 2 } } ( 1 + y ) ^ { - \frac { q - 1 } { 2 } } f ( x ) \} \parallel _ { L ^ { \infty } } } \end{array}$ , therefore we need the following preliminaries first. 

Remark 2.5. $( 1 - y ( x ) ^ { 2 } ) ^ { \lambda - \frac { 1 } { 2 } } ( 1 + y ( x ) ) ^ { - \frac { q - 1 } { 2 } } f ( x )$ has up to t-th derivatives in x, where $\begin{array} { r } { t = \lfloor \frac { \lambda } { q } \rfloor - 1 \leq \frac { \lambda } { q } } \end{array}$ , the largest integer below $\begin{array} { r } { \frac { \lambda } { q } - 1 } \end{array}$ . 

It is easy to observe that 

$$
\frac {d ^ {n}}{d x ^ {n}} (1 - y (x) ^ {2}) ^ {\lambda - \frac {1}{2}} = A ^ {n} Y _ {1} ^ {n} Y _ {2} ^ {n} Y _ {3} ^ {n}, \quad 0 \leq n \leq t \tag {2.8}
$$

where 

$$
A = \frac {2 ^ {q}}{2 q}, \qquad Y _ {1} ^ {n} = (1 - y (x) ^ {2}) ^ {\lambda - \frac {1}{2} - q n}, \qquad Y _ {2} ^ {n} = (1 - y (x)) ^ {n (q - 1)},
$$

and $Y _ { 3 } ^ { n }$ satisfies the following recursive relation: 

$$
Y _ {3} ^ {0} = 1
$$

$$
Y _ {3} ^ {n + 1} = - \left[ (2 \lambda - q n - n - 1) y (x) + n (q - 1) \right] Y _ {3} ^ {n} + (1 - y (x) ^ {2}) \frac {d}{d y} Y _ {3} ^ {n}, \quad 0 \leq n <   t.
$$

It is easy to show that $Y _ { 3 } ^ { n }$ is an n-th degree polynomial of y. We have the following estimate on $Y _ { 3 } ^ { n }$ . 

Lemma 2.6. We have, for $0 \leq n \leq t ,$ 

$$
\left| Y _ {3} ^ {n} \right| \leq (2 \lambda) ^ {n}, \quad y \in [ - 1, 1 ]. \tag {2.9}
$$

Proof. Let us rewrite $\begin{array} { r } { Y _ { 3 } ^ { n } = \sum _ { i = 0 } ^ { n } a _ { i } y ^ { i } } \end{array}$ and denote $\begin{array} { r } { S _ { n } = \sum _ { i = 0 } ^ { n } | a _ { i } | . } \end{array}$ 

If $\vert S _ { n } \vert \leq ( 2 \lambda ) ^ { n }$ , then $\begin{array} { r } { | Y _ { 3 } ^ { n } | \leq \sum _ { i = 0 } ^ { n } | a _ { i } | | y | ^ { i } \leq \sum _ { i = 0 } ^ { n } | a _ { i } | \leq ( 2 \lambda ) ^ { n } . } \end{array}$ 

Therefore, to prove $| Y _ { 3 } ^ { \check { n } } | \leq ( \overline { { 2 \lambda } } ) ^ { n }$ , we only need to prove $| S _ { n } | \leq ( 2 \lambda ) ^ { n }$ , for which we will use induction. 

$1 . n = 0 \colon Y _ { 3 } ^ { 0 } = 1 , s 0 \bar { S _ { 0 } } = 1 \leq ( 2 \lambda ) ^ { 0 } .$ 

2. Suppose $S _ { n } \leq ( 2 \lambda ) ^ { n }$ . Then for n + 1, 

$$
\left| Y _ {3} ^ {n + 1} \right| \leq | - [ (2 \lambda - q n - n - 1) y + n (q - 1) ] Y _ {3} ^ {n} | + \left| (1 - y ^ {2}) \frac {d}{d y} Y _ {3} ^ {n} \right|.
$$

Therefore, 

$$
\begin{array}{l} \left| S _ {n + 1} \right| \leq [ (2 \lambda - q n - n - 1) + n (q - 1) ] S _ {n} + 2 n S _ {n} \\ \leq [ (2 \lambda - q n - n - 1) + n (q - 1) + 2 n ] (2 \lambda) ^ {n} \\ \leq (2 \lambda) ^ {n + 1}. \\ \end{array}
$$

Thus, we have proved $| Y _ { 3 } ^ { n } | \leq ( 2 \lambda ) ^ { n }$ .  

Lemma 2.7. We have, for $0 \leq n \leq t ,$ , 

$$
\left| \frac {d ^ {n}}{d x ^ {n}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \right| \leq (2 A \lambda) ^ {n} (1 - y) ^ {n (q - 1)} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - q n}. \tag {2.10}
$$

Proof. We only need to use (2.8). 

Now let us rewrite the function 

$$
(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x) (1 + y) ^ {- \frac {q - 1}{2}} = C (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \left\{a (x) (1 + x) ^ {s _ {1}} + b (x) (1 + x) ^ {s _ {2}} \right\} \tag {2.11}
$$

where 

$$
s _ {1} = - \frac {q - 1}{2 q} \in \left(- \frac {1}{2}, 0\right), \quad s _ {2} = \frac {p}{q} - \frac {q - 1}{2 q} \in \left(- \frac {1}{2}, 1\right). \tag {2.12}
$$

Lemma 2.8. We have, for $0 \leq n \leq t a n d j = 1 , 2$ , 

$$
\left| \frac {d ^ {n}}{d x ^ {n}} \left[ (1 + x) ^ {s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \right] \right| \leq C (2 A \lambda) ^ {n} (1 + y) ^ {q s _ {j}} (1 - y) ^ {n (q - 1)} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - q n} \left(1 + \frac {1}{e}\right) ^ {n},
$$

where C is at most a constant multiplied by ${ \sqrt { n } } .$ 

Proof. We have 

$$
\begin{array}{l} \left| \frac {d ^ {n}}{d x ^ {n}} \left[ (1 + x) ^ {s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \right] \right| = \left| \sum_ {i = 0} ^ {n} \binom {n} {i} \frac {d ^ {i}}{d x ^ {i}} (1 + x) ^ {s _ {j}} \frac {d ^ {n - i}}{d x ^ {n - i}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \right| \\ = \left| \sum_ {i = 1} ^ {n} \binom {n} {i} s _ {j} \cdot (s _ {j} - 1) \cdot (s _ {j} - 2) \dots (s _ {j} + 1 - i) (1 + x) ^ {s _ {j} - i} \frac {d ^ {n - i}}{d x ^ {n - i}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \right| \\ + \left| (1 + x) ^ {s _ {j}} \frac {d ^ {n}}{d x ^ {n}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \right| \\ \leq \sum_ {i = 1} ^ {n} \binom {n} {i} i! (A q) ^ {(i - s _ {j})} (1 + y) ^ {q (s _ {j} - i)} (2 A \lambda) ^ {n - i} (1 - y) ^ {(q - 1) (n - i)} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - q (n - i)} \\ + (A q) ^ {- s _ {j}} (1 + y) ^ {q s _ {j}} (2 A \lambda) ^ {n} (1 - y) ^ {(q - 1) n} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - q n} \\ = C (2 A \lambda) ^ {n} (1 - y) ^ {(q - 1) n} (1 + y) ^ {q s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - n q} \left[ \sum_ {i = 1} ^ {n} \binom {n} {i} i! \left(\frac {(1 - y) q}{2 \lambda}\right) ^ {i} + 1 \right] \\ \leq C (2 A \lambda) ^ {n} (1 - y) ^ {(q - 1) n} (1 + y) ^ {q s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - n q} \sum_ {i = 0} ^ {n} \binom {n} {i} \left(\frac {(1 - y) n q}{2 e \lambda}\right) ^ {i} \\ \end{array}
$$

$$
\begin{array}{l} \leq C (2 A \lambda) ^ {n} (1 - y) ^ {(q - 1) n} (1 + y) ^ {q s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - n q} \sum_ {i = 0} ^ {n} \binom {n} {i} \left(\frac {1}{e}\right) ^ {i} \\ = C (2 A \lambda) ^ {n} (1 - y) ^ {(q - 1) n} (1 + y) ^ {q s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - n q} \left(1 + \frac {1}{e}\right) ^ {n} \\ \end{array}
$$

where we have used Lemma 2.7 in the third inequality, (2.4) for the fifth inequality, and the binomial formula for the last equality.  

Lemma 2.9. We have, for $0 \leq n \leq t a n d j = 1 , 2$ 

$$
\left| \frac {d ^ {n}}{d x ^ {n}} \left[ (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} (1 + x) ^ {s _ {j}} a (x) \right] \right| \leq C (\rho) (1 + y) ^ {q s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - q n} (2 ^ {q} A \lambda) ^ {n} \left[ 1 + \frac {1}{e} \left(1 + \frac {1}{2 ^ {q} A q \rho}\right) \right] ^ {n},
$$

where a(x) is an analytic function satisfying Assumption 1.1, and $s _ { j } \in ( - \frac { 1 } { 2 } , 1 )$ are defined in (2.12). $C ( \rho )$ is at most a constant multiplied by n. 

Proof. We have 

$$
\begin{array}{l} \left| \frac {d ^ {n}}{d x ^ {n}} \left[ (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} (1 + x) ^ {s _ {j}} a (x) \right] \right| = \left| \sum_ {i = 0} ^ {n} \binom {n} {i} \frac {d ^ {i}}{d x ^ {i}} \left[ (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} (1 + x) ^ {s _ {j}} \right] \frac {d ^ {n - i}}{d x ^ {n - i}} a (x) \right| \\ \leq \sum_ {i = 0} ^ {n} \binom {n} {i} (2 A \lambda) ^ {i} (1 + y) ^ {q s _ {j}} (1 - y) ^ {(q - 1) i} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - i q} \left(1 + \frac {1}{e}\right) ^ {i} \frac {C (\rho) (n - i) !}{\rho^ {n - i}} \\ \leq C (\rho) (1 + y) ^ {q s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} \sum_ {i = 0} ^ {n} \binom {n} {i} \left[ 2 A \lambda (1 - y) ^ {q - 1} (1 - y ^ {2}) ^ {- q} \left(1 + \frac {1}{e}\right) \right] ^ {i} \left(\frac {n}{e \rho}\right) ^ {n - i} \\ = C (\rho) (1 + y) ^ {q s _ {j}} \left(1 - y ^ {2}\right) ^ {\lambda - \frac {1}{2}} \left[ 2 A \lambda (1 - y) ^ {q - 1} \left(1 - y ^ {2}\right) ^ {- q} \left(1 + \frac {1}{e}\right) + \frac {n}{e \rho} \right] ^ {n} \\ = C (\rho) (1 + y) ^ {q s _ {j}} \left(1 - y ^ {2}\right) ^ {\lambda - \frac {1}{2} - q n} (2 A \lambda) ^ {n} \left[ (1 - y) ^ {q - 1} \left(1 + \frac {1}{e}\right) + \frac {n}{2 A \lambda e \rho} \right] ^ {n} \\ \leq C (\rho) (1 + y) ^ {q s _ {j}} \left(1 - y ^ {2}\right) ^ {\lambda - \frac {1}{2} - q n} (2 A \lambda) ^ {n} \left[ 2 ^ {q - 1} \left(1 + \frac {1}{e}\right) + \frac {n}{2 A \lambda e \rho} \right] ^ {n} \\ = C (\rho) (1 + y) ^ {q s _ {j}} \left(1 - y ^ {2}\right) ^ {\lambda - \frac {1}{2} - q n} \left(2 ^ {q} A \lambda\right) ^ {n} \left[ 1 + \frac {1}{e} + \frac {n}{2 ^ {q} A \lambda e \rho} \right] ^ {n} \\ \leq C (\rho) (1 + y) ^ {q s _ {j}} \left(1 - y ^ {2}\right) ^ {\lambda - \frac {1}{2} - q n} \left(2 ^ {q} A \lambda\right) ^ {n} \left[ 1 + \frac {1}{e} \left(1 + \frac {1}{2 ^ {q} A q \rho}\right) \right] ^ {n} \\ \end{array}
$$

where we have used Lemma 2.8 and Assumption 1.1 in the second inequality, (2.4) for the third inequality, and the binomial formula for the fourth equality. 

We have the following properties for the Gegenbauer expansion of analytic functions, which will be used in Section 4. 

Lemma 2.10 (The Exponential Decay of the Regularization Error). Let f (x) be an analytic function on [−1, 1] satisfying Assumption 1.1. Let $\hat { f } ^ { \lambda } ( l )$ , with $0 \leq l \leq m$ , be its Gegenbauer coefficients defined by 

$$
\hat {f} ^ {\lambda} (l) = \frac {1}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} (1 - x ^ {2}) ^ {\lambda - \frac {1}{2}} f (x) C _ {l} ^ {\lambda} (x) d x. \tag {2.13}
$$

Furthermore, assume $\lambda = \gamma m ,$ , where γ is a positive constant, then the regularization error defined by 

$$
R E (\lambda , m) = \max _ {- 1 \leq x \leq 1} \left| f (x) - \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (x) \right| \tag {2.14}
$$

satisfies 

$$
R E (\gamma m, m) \leq C q ^ {m} \tag {2.15}
$$

where q is given by 

$$
q = \frac {(1 + 2 \gamma) ^ {1 + 2 \gamma}}{\rho 2 ^ {1 + 2 \gamma} \gamma^ {\gamma} (1 + \gamma) ^ {1 + \gamma}} \tag {2.16}
$$

which is always less than 1. In particular, if $\gamma = 1$ and $m = \beta N$ where $\beta$ is a positive constant, then 

$$
R E (\beta N, \beta N) \leq C q ^ {N} \tag {2.17}
$$

with 

$$
q = \left(\frac {2 7}{3 2 \rho}\right) ^ {\beta}. \tag {2.18}
$$

Proof. The proof of this lemma can be found in [1]. 

# 3. Truncation error

Consider the function in the form of 

$$
f (x) = a (x) + b (x) (1 + x) ^ {s}
$$

where s is a given constant $\textstyle 0 < s = { \frac { p } { q } } < 1$ with relatively prime positive integers p and q, and $a ( x )$ and $b ( \boldsymbol x )$ are analytic functions satisfying Assumption 1.1. 

We assume that the point values $\{ f ( x _ { i } ) \}$ on N Gaussian points are given. We are interested in recovering the first m coefficients in the Gegenbauer expansion of $\dot { \boldsymbol { f } } ( \boldsymbol { x } )$ . For the function ${ \bar { f } } ( y ) = f ( x ( y ) )$ ), we have the usual Gegenbauer expansion with the basis $\{ C _ { l } ^ { \lambda } ( y ) \}$ : 

$$
f (x (y)) = \bar {f} (y) = \sum_ {l = 0} ^ {\infty} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x))
$$

where the Gegenbauer coefficients ${ \hat { f } } ^ { \lambda } ( l )$ are given by 

$$
\hat {f} ^ {\lambda} (l) = \frac {1}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x (y)) C _ {l} ^ {\lambda} (y) d y. \tag {3.1}
$$

Our goal here is to find a good approximation to the first $m \sim N$ Gegenbauer coefficients $\hat { f } ^ { \lambda } ( l )$ in this expansion. Based on the known point values $\{ f ( x _ { i } ) \}$ , we define not the usual interpolation polynomial $I _ { N } f ( x )$ , but 

$$
(f ^ {\lambda}) _ {N} (x) = I _ {N} \left(\frac {(1 - y (x) ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{\sqrt {A _ {d y} ^ {d x}}}\right) = I _ {N} \left((1 - y (x) ^ {2}) ^ {\lambda - \frac {1}{2}} (1 + y (x)) ^ {- \frac {q - 1}{2}} f (x)\right). \tag {3.2}
$$

Intuitively, the function being interpolated has about $\begin{array} { l } { \displaystyle { \frac { \lambda } { q } } \ \sim \ N } \end{array}$ continuous derivatives, hence the interpolation would produce nice error estimates. 

Our candidate for approximating the Gegenbauer coefficients $\hat { f } ^ { \lambda } ( l )$ is: 

$$
\hat {g} ^ {\lambda} (l) = \frac {1}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} (1 + y) ^ {\frac {q - 1}{2}} (f ^ {\lambda}) _ {N} (x (y)) C _ {l} ^ {\lambda} (y) d y. \tag {3.3}
$$

Definition 3.1. The truncation error is defined as 

$$
T E (\lambda , m, N) = \max _ {- 1 \leq y \leq 1} \left| \sum_ {l = 0} ^ {m} (\hat {f} ^ {\lambda} (l) - \hat {g} ^ {\lambda} (l)) C _ {l} ^ {\lambda} (y) \right|. \tag {3.4}
$$

Let us start to estimate the truncation error. 

Lemma 3.2. The truncation error is bounded by 

$$
T E (\lambda , m, N) \leq \frac {C (m + 1) (m + \lambda) \Gamma (m + 2 \lambda)}{\sqrt {\lambda} m ! \Gamma (2 \lambda)} \frac {1}{N ^ {t}} \left\| \frac {d ^ {t}}{d x ^ {t}} \left\{(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} (1 + y) ^ {- \frac {q - 1}{2}} f (x) \right\} \right\| _ {L ^ {\infty}}
$$

where $t = \lfloor \textstyle { \frac { \lambda } { q } } \rfloor - 1$ 

Proof. We have 

$$
\begin{array}{l} T E (\lambda , m, N) = \max _ {- 1 \leq y \leq 1} \left| \sum_ {l = 0} ^ {m} (\hat {f} ^ {\lambda} (l) - \hat {g} ^ {\lambda} (l)) C _ {l} ^ {\lambda} (y) \right| \\ \leq \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1)}{h _ {l} ^ {\lambda}} \left| \int_ {- 1} ^ {1} \left\{(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x) - (1 + y) ^ {\frac {q - 1}{2}} I _ {N} \left(\frac {(1 - y (x) ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y (x)) ^ {\frac {q - 1}{2}}}\right) \right\} C _ {l} ^ {\lambda} (y) d y \right| \\ \leq \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1) ^ {2}}{h _ {l} ^ {\lambda}} \int_ {- 1} ^ {1} \left| (1 + y) ^ {\frac {q - 1}{2}} \left\{\frac {(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y) ^ {\frac {q - 1}{2}}} - I _ {N} \left(\frac {(1 - y (x) ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y (x)) ^ {\frac {q - 1}{2}}}\right) \right\} \right| d y \\ \leq C \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1) ^ {2}}{h _ {l} ^ {\lambda}} \left(\int_ {- 1} ^ {1} \omega (x) \frac {d x}{d y} \left[ \frac {(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y) ^ {\frac {q - 1}{2}}} - I _ {N} \left(\frac {(1 - y (x) ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y (x)) ^ {\frac {q - 1}{2}}}\right) \right] ^ {2} d y\right) ^ {\frac {1}{2}} \\ \cdot \left(\int_ {- 1} ^ {1} \frac {1}{\omega (x)} d y\right) ^ {\frac {1}{2}} \\ = C \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1) ^ {2}}{h _ {l} ^ {\lambda}} \left\| \frac {(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y) ^ {\frac {q - 1}{2}}} - I _ {N} \left(\frac {(1 - y (x) ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y (x)) ^ {\frac {q - 1}{2}}}\right) \right\| _ {L _ {\omega} ^ {2}} \\ \leq C \sum_ {l = 0} ^ {m} \frac {C _ {l} ^ {\lambda} (1) ^ {2}}{h _ {l} ^ {\lambda}} \frac {1}{N ^ {t}} \left\| \frac {d ^ {t}}{d x ^ {t}} \left\{\frac {(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y) ^ {\frac {q - 1}{2}}} \right\} \right\| _ {L _ {\infty}} \\ \leq \frac {C (m + 1) (m + \lambda) \Gamma (m + 2 \lambda)}{\sqrt {\lambda} m ! \Gamma (2 \lambda)} \frac {1}{N ^ {t}} \left\| \frac {d ^ {t}}{d x ^ {t}} \left\{\frac {(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y) ^ {\frac {q - 1}{2}}} \right\} \right\| _ {L ^ {\infty}} \\ \end{array}
$$

where $\omega ( x )$ is the weight function of the basis $\{ \phi _ { k } \}$ , and the definitions of ${ \hat { f } } ^ { \lambda } ( l )$ in (3.1) and $\hat { g } ^ { \lambda } ( l )$ in (3.3) are used in the second inequality, (2.2) is used in the second and third inequalities, the Cauchy–Schwarz inequality and dxdy = (1+y)q−1A $\begin{array} { r } { \frac { d x } { d y } = \frac { ( 1 + y ) ^ { q - 1 } } { A } } \end{array}$ are used in the fourth inequality, Lemma 2.4 is used in the sixth inequality, and (2.1) and (2.5) are used in the last inequality. 

From the estimates in Lemma 2.9, we obtain 

Theorem 3.3 (The Exponential Decay of the Truncation Error). Let λ = αN, m = βN with $0 < \alpha , \beta < 1$ , then 

$$
T E (\alpha N, \beta N, N) \leq C (\rho) q _ {T} ^ {N}
$$

with 

$$
q _ {T} = \frac {(\beta + 2 \alpha) ^ {\beta + 2 \alpha}}{\beta^ {\beta} \alpha^ {2 \alpha}} \left[ \frac {\alpha}{2 q} (1 + \delta) \right] ^ {\frac {\alpha}{q}}, \quad \delta = \frac {1}{e} \left(1 + \frac {1}{2 ^ {q} A q \rho}\right).
$$

When we choose $\beta = \gamma \alpha ,$ , i.e. $m = \gamma \lambda$ , we have 

$$
q _ {T} = \left(\frac {(\gamma + 2) ^ {\gamma + 2}}{\gamma^ {\gamma}} \left[ \frac {\alpha}{2 q} (1 + \delta) \right] ^ {\frac {1}{q}}\right) ^ {\alpha}.
$$

If we choose α to satisfy 

$$
\alpha <   \frac {2 q}{(1 + \delta)} \left(\frac {\gamma^ {\gamma}}{(2 + \gamma) ^ {(2 + \gamma)}}\right) ^ {q},
$$

then $q _ { T } < 1$ 

Proof. We have 

$$
\begin{array}{l} T E (\lambda , m, N) = \max _ {- 1 \leq y \leq 1} \left| \sum_ {l = 0} ^ {m} (\hat {f} ^ {\lambda} (l) - \hat {g} ^ {\lambda} (l)) C _ {l} ^ {\lambda} (y) \right| \\ \leq \frac {C (m + 1) (m + \lambda) \Gamma (m + 2 \lambda)}{\sqrt {\lambda} m ! \Gamma (2 \lambda)} \frac {1}{N ^ {t}} \left\| \frac {d ^ {t}}{d x ^ {t}} \left\{\frac {(1 - y ^ {2}) ^ {\lambda - \frac {1}{2}} f (x)}{(1 + y) ^ {\frac {q - 1}{2}}} \right\} \right\| _ {L ^ {\infty}} \\ \end{array}
$$

$$
\begin{array}{l} \leq \sum_ {j = 1} ^ {2} \frac {C (\rho) (m + 1) (m + \lambda) \Gamma (m + 2 \lambda)}{\sqrt {\lambda} m ! \Gamma (2 \lambda)} \frac {1}{N ^ {t}} (1 + y) ^ {q s _ {j}} (1 - y ^ {2}) ^ {\lambda - \frac {1}{2} - q t} (2 ^ {q} A \lambda) ^ {t} (1 + \delta) ^ {t} \\ \leq \left(\sum_ {j = 1} ^ {2} 2 ^ {q | s _ {j} |} C (\rho)\right) \frac {(m + 1) (m + \lambda) \Gamma (m + 2 \lambda)}{\sqrt {\lambda} m ! \Gamma (2 \lambda)} \left(\frac {2 ^ {q} A \lambda (1 + \delta)}{N}\right) ^ {t} \\ \end{array}
$$

where we have used Lemma 3.2 in the second inequality, and we need Lemma 2.9 in the third inequality and $t = \lfloor \textstyle { \frac { \lambda } { q } } \rfloor - 1 \leq$ $\begin{array} { r } { \frac { \lambda - \frac { 1 } { 2 } } { q } - \frac { 1 } { 2 } \leq \frac { \lambda - \frac { 1 } { 2 } } { q } } \end{array}$ λ− 1 q in the fourth inequality. Therefore, we have 

$$
T E (\lambda , m, N) \leq C (\rho) \frac {(m + 1) (m + \lambda) \Gamma (m + 2 \lambda)}{\sqrt {\lambda} m ! \Gamma (2 \lambda)} \left(\frac {2 ^ {q} A \lambda (1 + \delta)}{N}\right) ^ {t},
$$

hence 

$$
\begin{array}{l} T E (\alpha N, \beta N, N) \leq C (\rho) \left(\frac {(2 \alpha + \beta) ^ {2 \alpha + \beta}}{\beta^ {\beta} (2 \alpha) ^ {2 \alpha}}\right) ^ {N} \left(2 ^ {q} A \alpha (1 + \delta)\right) ^ {t} \\ \leq C (\rho) \left(\frac {(2 \alpha + \beta) ^ {2 \alpha + \beta}}{\beta^ {\beta} (2 \alpha) ^ {2 \alpha}}\right) ^ {N} \left(2 ^ {q} A \alpha (1 + \delta)\right) ^ {\frac {\alpha}{q} N} \\ = C (\rho) q ^ {N} \\ \end{array}
$$

with 

$$
q = \left(\frac {(2 \alpha + \beta) ^ {2 \alpha + \beta}}{\beta^ {\beta} \alpha^ {2 \alpha}}\right) \left(\frac {\alpha}{2 q} (1 + \delta)\right) ^ {\frac {\alpha}{q}}
$$

where we have used Stirling’s formula (2.4) in the first inequality and $\begin{array} { r } { \frac { \lambda } { q } - 2 \leq t = \lfloor \frac { \lambda } { q } \rfloor - 1 \leq \frac { \lambda } { q } } \end{array}$ in the second inequality. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/21786fb1-7f2e-43f0-b753-a48f89c56ea3/17ffc6cb1686cfe36d654be9b0e871308ab8220abe8992c13a7c3f02f37d55dd.jpg)


# 4. Regularization error and the main theorem

In this section we would like to establish the error estimates for approximating $f ( x ) = a ( x ) + b ( x ) ( 1 + x ) ^ { s } 0 \mathrm { n } [ - 1$ , 1] by its Gegenbauer expansion based on the Gegenbauer polynomials $C _ { n } ^ { \lambda } ( y ( x ) )$ ). Since our goal is to remove the Gibbs phenomenon, we will use the maximum norm. As shown in the last section, we would need both λ and m to grow linearly with N in order to obtain exponential accuracy in the maximum norm. Thus, we will consider the case of large λ and m in this section. 

We will assume that a(x) and b(x) are analytic functions on [−1, 1] satisfying Assumption 1.1. We would like to estimate the regularization error in the maximum norm. We consider the Gegenbauer partial sum of the first m terms for the function f (x) given by: 

$$
f ^ {\lambda , m} (x) = \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)), \tag {4.1}
$$

with the Gegenbauer coefficients $\hat { f } ^ { \lambda } ( l )$ defined by (3.1). 

Definition 4.1. The regularization error is defined by 

$$
\begin{array}{l} R E (\lambda , m) = \max _ {- 1 \leq x \leq 1} \left| f (x) - \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)) \right| \\ = \max _ {- 1 \leq y \leq 1} \left| f (x (y)) - \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y) \right|. \\ \end{array}
$$

Since the function $f ( x ( y ) )$ is analytic with respect to the variable y, we can get the following result from Lemma 2.10. 

Theorem 4.2 (The Exponential Decay of the Regularization Error). For the function $f ( x ) = a ( x ) + b ( x ) ( 1 + x ) ^ { s }$ , with analytic functions a(x) and b(x) satisfying Assumption 1.1, if we assume $\lambda = \alpha N$ and $m = \gamma \lambda ,$ , then 

$$
\max _ {- 1 \leq x \leq 1} \left| f (x) - \sum_ {l = 0} ^ {m} \hat {f} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)) \right| \leq C q _ {R} ^ {N}
$$

where 

$$
q _ {R} = \left(\frac {(\gamma + 2) ^ {\gamma + 2}}{\rho^ {\gamma} 2 ^ {\gamma + 2} (\gamma + 1) ^ {\gamma + 1}}\right) ^ {\alpha}
$$

which is always less than 1. 

We can then combine the estimates for the truncation error and the regularization error to obtain the main theorem: 

Theorem 4.3 (Removal of the Gibbs Phenomenon). Consider a function in the form of $f ( x ) = a ( x ) + b ( x ) ( 1 + x ) ^ { s }$ , with given fractional constant $0 < \bar { s ^ { } } = \frac { p } { q } < 1$ , and a(x) and b(x) are analytic functions satisfying Assumption 1.1. Assume that the point values {f (xi)} at N Gaussian points are given. Let $\hat { g } ^ { \lambda } ( l )$ , with $0 \leq l \leq m$ , be the Gegenbauer expansion coefficients, defined in (3.2) and (3.3). Then for λ = αN and m = γ λ with 

$$
\alpha <   \frac {2 q}{(1 + \delta)} \left(\frac {\gamma^ {\gamma}}{(2 + \gamma) ^ {(2 + \gamma)}}\right) ^ {q}, \quad \delta = \frac {1}{e} \left(1 + \frac {1}{2 ^ {q} A q \rho}\right),
$$

we have 

$$
\max _ {- 1 \leq x \leq 1} \left| f (x) - \sum_ {l = 0} ^ {m} \hat {g} ^ {\lambda} (l) C _ {l} ^ {\lambda} (y (x)) \right| \leq C \left(q _ {T} ^ {N} + q _ {R} ^ {N}\right)
$$

where 

$$
q _ {T} = \left(\frac {(\gamma + 2) ^ {\gamma + 2}}{\gamma^ {\gamma}} \left[ \frac {\alpha}{2 q} (1 + \delta) \right] ^ {\frac {1}{q}}\right) ^ {\alpha} <   1,
$$

$$
q _ {R} = \left(\frac {(\gamma + 2) ^ {\gamma + 2}}{\rho^ {\gamma} 2 ^ {\gamma + 2} (\gamma + 1) ^ {\gamma + 1}}\right) ^ {\alpha} <   1.
$$

Proof. We can get the proof by combining the results of Theorems 3.3 and 4.2. 

Remark 4.4. In the proof, no attempt has been made to optimize the parameters. 

# 5. Numerical results

In this section, we give two numerical examples to illustrate the result. We have tested both the Fourier collocation and the Chebyshev collocation. 

Example 5.1. We take the function 

$$
f (x) = \cos (x) + \sin (x) \sqrt {1 + x} \tag {5.1}
$$

and try to recover the pointwise values of this function over [−1, 1]. 

First, we consider the Fourier collocation. We assume the point values {f (xi)} on the 2N + 1 uniform points: 

$$
x _ {i} = \frac {2 i}{2 N + 1}, \quad i = - N, \dots , N \tag {5.2}
$$

are given. The parameters are chosen as 

$$
\lambda = 0. 2 N, \quad m = 0. 0 7 5 N. \tag {5.3}
$$

The choice of the parameters here and below are guided by the theory in previous sections. 

Second, we consider the Chebyshev collocation. We assume the point values {f (x )} on the N + 1 Chebyshev collocation points: 

$$
x _ {i} = \cos \left(\frac {\pi (2 i + 1)}{2 N + 2}\right), \quad i = 0, \dots , N \tag {5.4}
$$

are given. The parameters are chosen as 

$$
\lambda = 0. 2 N, \quad m = 0. 1 N. \tag {5.5}
$$

In Fig. 5.1 we show the errors, in the logarithm scale, for N = 40, 80, 160, 320 and 640. In Table 5.1, we show the maximum errors for each N and the orders of convergence. We can clearly find out that both the Fourier case and the Chebyshev case provide exponentially decaying errors, and the errors are smaller in the Chebyshev case than in the Fourier case. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/21786fb1-7f2e-43f0-b753-a48f89c56ea3/b94f0dd5fc404af8270fb2a9c0d242998d6d5c4313f013b86301baf3ab6a4253.jpg)



(a) Fourier case (λ = 0.2N, m = 0.075N).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/21786fb1-7f2e-43f0-b753-a48f89c56ea3/1fcc794eef948ee5baf66c12e26348d0c8dfacb4313d885a41ee2d5137901daf.jpg)



(b) Chebyshev case (λ = 0.2N, m = 0.1N).



Fig. 5.1. Pointwise errors in the logarithm scale.



Table 5.1 Maximum error table.


<table><tr><td rowspan="2">N</td><td colspan="2">Fourier case<eq>\lambda = 0.2N, m = 0.075N</eq></td><td colspan="2">Chebyshev case<eq>\lambda = 0.2N, m = 0.1N</eq></td></tr><tr><td><eq>L^{\infty}</eq> error</td><td>Order</td><td><eq>L^{\infty}</eq> error</td><td>Order</td></tr><tr><td>40</td><td>3.46E-01</td><td></td><td>1.34E-01</td><td></td></tr><tr><td>80</td><td>2.01E-02</td><td>4.11</td><td>2.49E-03</td><td>5.75</td></tr><tr><td>160</td><td>2.19E-05</td><td>9.84</td><td>5.17E-08</td><td>15.56</td></tr><tr><td>320</td><td>5.93E-11</td><td>18.49</td><td>3.58E-19</td><td>37.07</td></tr><tr><td>640</td><td>1.82E-18</td><td>24.96</td><td>8.48E-43</td><td>78.48</td></tr></table>

# Example 5.2. We take the function

$$
f (x) = \cos (x) + \sin (x) (1 + x) ^ {\frac {1}{3}} \tag {5.6}
$$

and again try to recover the pointwise values of this function over [−1, 1]. 

First, we consider the Fourier collocation. We assume the point values {f (xi)} on the 2N + 1 uniform points xi given by (5.2) are given. The parameters are chosen as 

$$
\lambda = 0. 2 N, \quad m = 0. 0 2 5 N. \tag {5.7}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/21786fb1-7f2e-43f0-b753-a48f89c56ea3/82560c9c71ca19e54fe084161517a540732e6df368892b61a1e88cdcc0f6c17c.jpg)



(a) Fourier case (λ = 0.2N, m = 0.025N).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-27/21786fb1-7f2e-43f0-b753-a48f89c56ea3/e3b3c095fbdef8ac4fdae2a91dbe7230916dc7873041081a7039ac84e1dc82a7.jpg)



(b) Chebyshev case $( \lambda = 0 . 2 N , m = 0 . 0 5 N )$ .



Fig. 5.2. Pointwise errors in the logarithm scale.



Table 5.2 Maximum error table.


<table><tr><td rowspan="2">N</td><td colspan="2">Fourier case<eq>\lambda = 0.2N, m = 0.025N</eq></td><td colspan="2">Chebyshev Case<eq>\lambda = 0.2N, m = 0.05N</eq></td></tr><tr><td><eq>L^{\infty}</eq> error</td><td>Order</td><td><eq>L^{\infty}</eq> error</td><td>Order</td></tr><tr><td>40</td><td>7.56E-01</td><td></td><td>3.49E-01</td><td></td></tr><tr><td>80</td><td>3.96E-01</td><td>0.93</td><td>5.25E-01</td><td>-0.59</td></tr><tr><td>160</td><td>6.00E-01</td><td>-0.60</td><td>4.25E-02</td><td>3.63</td></tr><tr><td>320</td><td>5.16E-02</td><td>3.54</td><td>3.95E-05</td><td>10.07</td></tr><tr><td>640</td><td>6.36E-05</td><td>9.66</td><td>1.15E-12</td><td>25.03</td></tr></table>

Second, we consider the Chebyshev collocation. We assume the point values {f (xi)} on the N + 1 Chebyshev collocation points given by (5.4) are given. The parameters are chosen as 

$$
\lambda = 0. 2 N, \quad m = 0. 0 5 N. \tag {5.8}
$$

In Fig. 5.2 we show the errors, in the logarithm scale, for N = 40, 80, 160, 320 and 640. In Table 5.2, we show the maximum errors for each N and the orders of convergence. We again clearly see that both the Fourier case and the Chebyshev case provide exponentially decaying errors, and errors are smaller in the Chebyshev case than in the Fourier case. 

Remark 5.3. No attempt is made to optimize the parameters for each individual case in the numerical examples. 

# 6. Concluding remarks

We have generalized the Gegenbauer polynomial based technique [5] to reconstruct approximations with exponential accuracy in the maximum norm, based on collocation point values of the function, from the case of piecewise analytic functions where in each piece the function is analytic in the closed interval, to piecewise functions where in each piece the function is analytic only in the open interval with end point singularities. Numerical results are provided to demonstrate the theory. This technique has potential applications for solving certain partial differential equations whose solutions have such end point singularities, and also in accurate numerical integration and differentiation of such functions. 

# References



[1] D. Gottlieb, C.-W. Shu, A. Solomonoff, H. Vandeven, On the Gibbs phenomenon I: recovering exponential accuracy from the Fourier partial sum of a non-periodic analytic function, Journal of Computational and Applied Mathematics 43 (1992) 81–98. 





[2] D. Gottlieb, C.-W. Shu, Resolution properties of the Fourier method for discontinuous waves, Computer Methods in Applied Mechanics and Engineering 116 (1994) 27–37. 





[3] D. Gottlieb, C.-W. Shu, On the Gibbs phenomenon III: recovering exponential accuracy in a sub-interval from a spectral partial sum of a piecewise analytic function, SIAM Journal on Numerical Analysis 33 (1996) 280–290. 





[4] D. Gottlieb, C.-W. Shu, On the Gibbs phenomenon IV: recovering exponential accuracy in a sub-interval from a Gegenbauer partial sum of a piecewise analytic function, Mathematics of Computation 64 (1995) 1081–1095. 





[5] D. Gottlieb, C.-W. Shu, On the Gibbs phenomenon V: recovering exponential accuracy from collocation point values of a piecewise analytic function, Numerische Mathematik 71 (1995) 511–526. 





[6] D. Gottlieb, C.-W. Shu, On the Gibbs phenomenon and its resolution, SIAM Review 30 (1997) 644–668. 





[7] H. Bateman, Higher Transcendental Functions, Vol. 2, McGraw-Hill, 1953. 





[8] E. Tadmor, The exponential accuracy of Fourier and Chebyshev differencing methods, SIAM Journal on Numerical Analysis 23 (1986) 1–10. 





[9] D. Funaro, Polynomial Approximation of Differential Equations, in: Lecture Notes in Physics, vol. 8, Springer-Verlag, 1992. 

