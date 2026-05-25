# Optimal asymptotic analyses on Laguerre and Hermite orthogonal approximation for functions of algebraic and logarithmic regularities⋆

Yali Zhanga,b, Guidong Liuc, Shuhuang Xiangb,∗ 

aSchool of Mathematics and Statistics, Yancheng Teachers University, Yancheng 224000, P. R. China bSchool of Mathematics and Statistics, Central South University, Changsha 410083, P. R. China cSchool of Mathematics, Nanjing Audit University, Nanjing 211815, P. R. China 

# Abstract

Based on the Hilb-type formula and van der Corput-type lemmas, we present optimal asymptotic estimates for the decay of the Laguerre and Hermite coefficients for functions with algebraic and logarithmic singularities, which in turn yield the convergence rates of the corresponding spectral orthogonal projections. Numerous examples are provided to verify the optimality of these asymptotic results. 

Keywords: van der Corput-type lemmas, Laguerre approximation, Hermite approximation, Asymptotic analysis, Convergence rate 2010 MSC: 41A10, 41A25 

# 1. Introduction

The approximation of functions using orthogonal polynomials plays a fundamental role in numerical analysis and scientific computing, particularly in spectral methods for solving differential equations, Gaussian quadrature rules, and expansions in special functions [21, 20, 6, 12, 27, 26, 11]. Among the various families of orthogonal polynomials, Laguerre and Hermite polynomials are of particularly significance for representing functions over semi-infinite and infinite domains, respectively. Their effectiveness in function 

⋆This work was supported by the National Natural Science Foundation of China (No. 12271528), the Natural Science Foundation of Jiangsu Province (BK20241087). 

∗Corresponding author 

Email addresses: ylzhang_math@163.com (Yali Zhang), liugd@nau.edu.cn (Guidong Liu), xiangsh@csu.edu.cn (Shuhuang Xiang) 

approximation has led to applications in various fields, including quantum mechanics, signal processing and numerical solutions of partial differential equations [10, 22, 16, 15, 7, 14, 9, 18]. Research on these polynomials has yielded many new and important results [2, 17, 4, 13, 3]. 

A prominent advantage of orthogonal polynomial approximations is their well-defined orthogonality properties, which facilitate stable numerical computations and rapid convergence. The rate of convergence is primarily determined by the regularity of the function being approximated. Classical results in spectral approximation theory focus on analytic or highly smooth functions, for which the expansion coefficients decay rapidly, leading to exponential convergence in appropriate norms. For example, Elliott and Tuan [8] established contour integral representations for the Laguerre and Hermite coefficients of analytic functions. Wang in [23], and Wang and Zhang in [24] investigated the rootexponential convergence of Laguerre and Hermite approximations for such functions. For functions with lower regularity, Shen and Wang [19] and Xiang [26] analyzed the algebraic convergence rates of order (n−s), where s is a parameter that depends on the regularity of the underlying functions. However, these results become suboptimal when the regularity parameter s is non-integer. More recently, Zhang et al. [28, 29] investigated optimal pointwise estimates for Laguerre approximations by a Hilb-type asymptotic formula for Laguerre polynomials. 

In this paper, we study the asymptotic behavior of Laguerre and Hermite expansion coefficients for functions with endpoint or interior singularities. We first derive several asymptotic estimates for integrals involving highly oscillatory Bessel functions. Following Hilb-type formulas, these estimates are subsequently applied to integrals involving Laguerre and Hermite polynomials. By employing Rodrigues’ formulas, we then establish the decay rates of the Laguerre and Hermite coefficients for functions with endpoint or interior singularities. It is important to emphasize that all of the derived estimates are optimal in the sense that the convergence rates cannot be improved, which will be confirmed by ample numerical experiments. Finally, these results are used to characterize the decay of the corresponding Laguerre and Hermite spectral orthogonal projections. 

The structure of this paper is as follows. In Section 2, we review the fundamental properties of Laguerre polynomials and Hermite polynomials. In Section 3, we derive the generalized van der Corput-type Lemmas for Bessel transforms. Sections 4 and 5 present the optimal decay rates of the Laguerre and Hermite coefficients, respectively, as well as the convergence rates of the corresponding spectral orthogonal projections for functions with algebraic and logarithmic singularities. Finally, Section 6 provides a brief summary of the main results and conclusions. 

# 2. Preliminaries

# 2.1. Gneralized Laguerre polynomials

The generalized Laguerre polynomials $L _ { n } ^ { ( \alpha ) } ( x )$ form an orthogonal basis on the semiinfinite interval $[ 0 , + \infty )$ with respect to the weight function $\omega _ { \alpha } ( x ) = x ^ { \alpha } \mathbf { e } ^ { - x }$ , where $\alpha > - 1$ . Specifically, they satisfy the orthogonality property 

$$
\int_ {0} ^ {+ \infty} \omega_ {\alpha} (x) L _ {n} ^ {(\alpha)} (x) L _ {m} ^ {(\alpha)} (x) \mathrm{d} x = \sigma_ {n} ^ {(\alpha)} \delta_ {m, n}, \quad \sigma_ {n} ^ {(\alpha)} = \frac {\Gamma (n + \alpha + 1)}{n !}, \tag {1}
$$

where $\delta _ { m , n }$ denotes the Kronecker delta symbol. 

Let f (x) be a suitably smooth function on $[ 0 , \infty )$ . Then it can be expanded in a Laguerre series [20, p. 265] as 

$$
f (x) = \sum_ {n = 0} ^ {\infty} a _ {n} (\alpha) L _ {n} ^ {(\alpha)} (x),
$$

where the coefficients $a _ { n } ( \alpha )$ are given by 

$$
a _ {n} (\alpha) = \frac {1}{\sigma_ {n} ^ {(\alpha)}} \int_ {0} ^ {\infty} f (x) L _ {n} ^ {(\alpha)} (x) \omega_ {\alpha} (x) \mathrm{d} x. \tag {2}
$$

A practical approximation of $f ( x )$ is obtained by truncating the Laguerre series after the first $N + 1 \mathrm { - t e r m s }$ 

$$
S _ {N} ^ {(\alpha)} [ f ] (x) = \sum_ {n = 0} ^ {N} a _ {n} (\alpha) L _ {n} ^ {(\alpha)} (x).
$$

The corresponding approximation error in the weighted $L ^ { 2 } { \mathrm { - n o r m } }$ with respect to $\omega _ { \alpha } ( x )$ is given by 

$$
\| f (x) - S _ {N} ^ {(\alpha)} [ f ] (x) \| _ {L _ {\omega \alpha} ^ {2} [ 0, \infty)} = \sqrt {\sum_ {n = N + 1} ^ {\infty} a _ {n} ^ {2} (\alpha) \sigma_ {n} ^ {(\alpha)}}. \tag {3}
$$

# 2.2. Hermite orthogonal polynomials

The Hermite polynomials $H _ { n } ( x )$ form an orthogonal basis on the entire real line $( - \infty , + \infty )$ with respect to the Gaussian weight function $\omega ( x ) = \mathtt { e } ^ { - x ^ { 2 } }$ . They are widely used in probability theory, quantum mechanics and signal processing. For a function $g ( x )$ defined on $( - \infty , + \infty )$ , the Hermite series expansion takes the form [20, p. 270] 

$$
g (x) = \sum_ {n = 0} ^ {\infty} h _ {n} H _ {n} (x),
$$

where the Hermite coefficients $h _ { n }$ are defined as 

$$
h _ {n} = \frac {1}{\gamma_ {n}} \int_ {- \infty} ^ {+ \infty} g (x) H _ {n} (x) \omega (x) \mathrm{d} x, \quad \gamma_ {n} = \sqrt {\pi} 2 ^ {n} n!. \tag {4}
$$

The Nth partial sum of the Hermite series provides an approximation to $g ( x )$ 

$$
S _ {N} [ g ] (x) = \sum_ {n = 0} ^ {N} h _ {n} H _ {n} (x).
$$

with the corresponding approximation error in the weighted $L ^ { 2 }$ -norm with respect to $\omega ( x )$ given by 

$$
\| g (x) - S _ {N} [ g ] (x) \| _ {L _ {\omega} ^ {2} (- \infty , \infty)} = \sqrt {\sum_ {n = N + 1} ^ {\infty} h _ {n} ^ {2} \gamma_ {n}}. \tag {5}
$$

Moreover, Hermite polynomials are closely connected to generalized Laguerre polynomials. In particular, they satisfy the following identities (see [21, 20]) 

$$
\begin{array}{l} H _ {2 n} (x) = (- 1) ^ {n} 2 ^ {2 n} n! L _ {n} ^ {(- 1 / 2)} \left(x ^ {2}\right), \\ H _ {2 n - 1} (x) = (- 1) ^ {n} 2 ^ {2 n + 1} \dots \tau^ {(1 / 2)} \left(x ^ {2}\right). \end{array} \tag {6}
$$

$$
H _ {2 n + 1} (x) = (- 1) ^ {n} 2 ^ {2 n + 1} n! x L _ {n} ^ {(1 / 2)} (x ^ {2}),
$$

where $L _ { n } ^ { ( - 1 / 2 ) } ( x ^ { 2 } )$ and $L _ { n } ^ { ( 1 / 2 ) } ( x ^ { 2 } )$ are generalized Laguerre polynomials of degree n with parameters $\alpha = - 1 / 2$ and $1 / 2$ , respectively. These relations illustrate the structural interplay between Hermite and Laguerre systems, which is of particular interest in asymptotic and approximation analyses. 

# 2.3. The asymptotic properties of Laguerre and Hermite polynomials

From an asymptotic point of view as $n  \infty$ , the generalized Laguerre polynomials exhibit a close connection with the Bessel function of the first kind. This relationship is formalized by the following Hilb-type asymptotic formula. 

Lemma 2.1. (Asymptotic formula of Hilb’s type [21, Theorem 8.22.4]) Let $\alpha > - 1$ and define $\tilde { n } = n + ( \alpha + 1 ) / 2$ . As $n  \infty$ , the generalized Laguerre polynomial satisfies 

$$
\mathrm{e} ^ {- x / 2} x ^ {\alpha / 2} L _ {n} ^ {(\alpha)} (x) = \frac {\tilde {n} ^ {- \alpha / 2} \Gamma (n + \alpha + 1)}{n !} J _ {\alpha} (2 \sqrt {\tilde {n} x}) + \left\{ \begin{array}{l l} x ^ {5 / 4} O \left(n ^ {\frac {\alpha}{2} - \frac {3}{4}}\right), & c n ^ {- 1} \leq x \leq \omega , \\ x ^ {\alpha / 2 + 2} O \left(n ^ {\alpha}\right), & 0 <   x \leq c n ^ {- 1}, \end{array} \right. \tag {7}
$$

where c and ω are fixed positive constants, and $J _ { \alpha } ( z )$ denotes the Bessel function of the first kind of order α. The bounds above hold uniformly for $0 < x \leq \omega .$ . In the special case $\alpha = 0$ , the last bound in (7) is to be replaced by $x ^ { 2 }$ log $( x ^ { - 1 } n ^ { - 1 } )$ . 

In addition to their asymptotic representations, the maximum behavior of Laguerre and Hermite polynomials is of significant interest in approximation theory. The following results provide sharp estimates for their maximum magnitudes in weighted forms. 

Lemma 2.2. ([21, Theorem 8.91.2]) Let α, $\lambda \in \mathbb { R } , a > 0$ and $0 < \eta < 4 .$ . As $n  \infty ,$ , the generalized Laguerre polynomial satisfies 

$$
\max _ {x \geq a} \mathrm{e} ^ {- x / 2} x ^ {\lambda} \left| L _ {n} ^ {(\alpha)} (x) \right| = \left\{ \begin{array}{l l} O \left(n ^ {\max \{\lambda - 1 / 2, \alpha / 2 - 1 / 4 \}}\right), & \text {   if   } a \leq x \leq (4 - \eta) n, \\ O \left(n ^ {\max \{\lambda - 1 / 3, \alpha / 2 - 1 / 4 \}}\right), & \text {   if   } x \geq a. \end{array} \right. \tag {8}
$$

Lemma 2.3. ([21, Theorem 8.91.3]) Let $\lambda \in \mathbb { R } , a > 0$ and $0 < \eta < 2$ . As $n  \infty$ , the Hermite polynomial satisfies 

$$
\max _ {| x | \geq a} \mathrm{e} ^ {- x ^ {2} / 2} x ^ {\lambda} \left| H _ {n} (x) \right| = \sqrt {2 ^ {n} n !} \left\{ \begin{array}{l l} O \left(n ^ {\max \{\lambda / 2 - 1 / 4, - 1 / 4 \}}\right), & \text {   if   } a \leq | x | \leq \sqrt {(2 - \eta) n}, \\ O \left(n ^ {\max \{\lambda / 2 - 1 / 1 2, - 1 / 4 \}}\right), & \text {   if   } | x | \geq a. \end{array} \right. \tag {9}
$$

These results play a fundamental role in establishing precise convergence rates and asymptotic estimates for Laguerre and Hermite approximations, particularly for functions exhibiting singularities. 

# 3. Some useful lemmas

In this section, we review several asymptotic estimates for integrals involving highly oscillatory Bessel functions, which arise through the Hilb-type formula (7) in connection with Laguerre polynomials. These results extend the classical van der Corput lemma for Fourier transforms to the following Bessel transforms. 

Lemma 3.1 ([25]). Let $b > a > 0 , \psi ( x ) \in C [ a , b ]$ and $\psi ^ { \prime } ( x ) \in L ^ { 1 } [ a , b ]$ . Then as $\omega  \infty$ , it holds that 

$$
\int_ {a} ^ {b} J _ {v} (\omega x) \psi (x) \mathrm{d} x = O \left(\omega^ {- 3 / 2}\right).
$$

Lemma 3.2 ([27]). Let $\alpha + \nu > - 1 , \beta > - 1$ and $\mu \in \mathbb { N } .$ . Suppose that $\psi ( x ) \in C [ 0 , b ]$ and $\psi ^ { \prime } ( x ) \in L ^ { 1 } [ 0 , b ]$ for a constant $b > 0$ . Then, as $\omega  \infty ,$ , the following asymptotic estimates hold 

$$
\int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {v} (\omega x) \psi (x) \mathrm{d} x = O \left(\max \left\{\frac {\ln^ {\mu} \omega}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\min \{\beta + 3 / 2 , 3 / 2 \}}} \right\}\right), \tag {10}
$$

and 

$$
\int_ {0} ^ {b} \ln^ {\mu} (b - x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \psi (x) \mathrm{d} x = O \left(\max \left\{\frac {1}{\omega^ {\alpha + 1}}, \frac {\ln^ {\mu} \omega}{\omega^ {\min \{\beta + 3 / 2 , 3 / 2 \}}} \right\}\right). \tag {11}
$$

Lemma 3.3. Let $\beta > - 1$ and $\mu \in \mathbb { N } .$ For any fixed constants a and b satisfying $0 < a <$ < $b < \infty ,$ , the following estimates hold uniformly for $t \in [ a , b ]$ 

$$
\int_ {a} ^ {t} \ln^ {\mu} (x - a) (x - a) ^ {\beta} J _ {v} (\omega x) \mathrm{d} x = O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right), \tag {12}
$$

and 

$$
\int_ {t} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} J _ {v} (\omega x) \mathrm{d} x = O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right). \tag {13}
$$

as $\omega  \infty .$ . 

Proof. Recall that for $z \in [ a , a + 1 / \omega ]$ , the Bessel function satisfies $J _ { \nu } ( \omega z ) = O ( \omega ^ { - 1 / 2 } )$ $( \sec , \mathsf { e . g . } , [ 2 1 , ( 1 . 7 1 . 7 ) ] )$ . Thus, for $t \in [ a , a + 1 / \omega ]$ , it follows that 

$$
\begin{array}{l} \int_ {a} ^ {t} \ln^ {\mu} (x - a) (x - a) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x \\ = \mathcal {O} \left(\omega^ {- 1 / 2}\right) \int_ {a} ^ {t} \ln^ {\mu} (x - a) (x - a) ^ {\beta} d x \\ = \mathcal {O} (\omega^ {- 1 / 2}) \left(\frac {\ln^ {\mu} (x - a) (x - a) ^ {\beta + 1}}{\beta + 1} \Big | _ {a} ^ {t} - \frac {\mu}{\beta + 1} \int_ {a} ^ {t} \ln^ {\mu - 1} (x - a) (x - a) ^ {\beta} \mathrm{d} x\right) \\ = \mathcal {O} \left(\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega\right) + \mathcal {O} \left(\omega^ {- \beta - 3 / 2} \ln^ {\mu - 1} \omega\right) + \dots + \mathcal {O} \left(\omega^ {- \beta - 3 / 2}\right) \\ = \mathcal {O} (\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega). \\ \end{array}
$$

For $t \in ( a + 1 / \omega , b ]$ , we split the integral and derive that 

$$
\begin{array}{l} \int_ {a} ^ {t} \ln^ {\mu} (x - a) (x - a) ^ {\beta} J _ {v} (\omega x) \mathrm{d} x \\ = \int_ {a} ^ {a + 1 / \omega} \ln^ {\mu} (x - a) (x - a) ^ {\beta} J _ {v} (\omega x) d + \int_ {a + 1 / \omega} ^ {t} \ln^ {\mu} (x - a) (x - a) ^ {\beta} J _ {v} (\omega x) d x \\ = O \left(\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega\right) + O \left(\omega^ {- 3 / 2}\right) \\ = \mathcal {O} \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right), \\ \end{array}
$$

where the estimate of the second integral follows from Lemma 3.1. Then we get Equation (12). A similar argument applies to the second integral leads to the estimate (13), by the symmetry of the analysis.  

Lemma 3.4. Suppose $\beta > - 1 , 0 < a < b < \infty , \mu \in \mathbb { N } .$ . Suppose that $\psi ( x ) \in C [ a , b ]$ and $\psi ^ { \prime } ( x ) \in L ^ { 1 } [ a , b ]$ . Then the following asymptotic results hold as $\omega \to \infty$ 

$$
\int_ {a} ^ {b} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {v} (\omega x) \mathrm{d} x = O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right) \tag {14}
$$

and 

$$
\int_ {a} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \psi (x) J _ {v} (\omega x) \mathrm{d} x = O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right). \tag {15}
$$

Proof. According to Equation (12), we have 

$$
\begin{array}{l} \int_ {a} ^ {b} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {v} (\omega x) d x \\ = \int_ {a} ^ {b} \psi (x) \left[ \int_ {a} ^ {x} \ln^ {\mu} (u - a) (u - a) ^ {\beta} J _ {v} (\omega u) d u \right] ^ {\prime} d x \\ = \psi (b) \int_ {a} ^ {b} \ln^ {\mu} (x - a) (x - a) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x - \int_ {a} ^ {b} \psi^ {\prime} (x) \left[ \int_ {a} ^ {x} \ln^ {\mu} (u - a) (u - a) ^ {\beta} J _ {\nu} (\omega u) \mathrm{d} u \right] \mathrm{d} x \\ \leq C _ {0} \left(| \psi (b) | + \int_ {a} ^ {b} | \psi^ {\prime} (x) | \mathrm{d} x\right) \cdot \max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\} \\ \leq C \cdot \max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}, \\ \end{array}
$$

where $C _ { 0 } , C$ are some constants independent of $x .$ 

The proof of (15) follows analogously by applying the same argument to Equation (13).  

Lemma 3.5. Let a and b be fixed finite constants, $\mu \in \mathbb { N } , \psi ( x ) \in C [ a , b ]$ and $\psi ^ { \prime } ( x ) \ \in$ $L ^ { 1 } [ a , b ]$ . Then the following asymptotic results hold as $\omega  \infty ;$ 

• Suppose the following conditions are satisfied to ensure the integrability of the integral: (i) $H 0 < a < b ,$ , assume $\beta > - 1 . \ ( i i ) \ I f 0 = a < b ,$ assume $\alpha + \delta + \beta + \nu > - 1$ . (iii) $H a < 0 \leq b ,$ , assume $\alpha + \delta + \nu > - 1$ and $\beta > - 1$ . (iv) ${ \mathit { I f a } } < { \mathit { b } } < 0$ , assume 

$\beta > - 1$ . Then we have 

$$
\begin{array}{l} \int_ {a} ^ {b} x ^ {\alpha} | x | ^ {\delta} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {v} (\omega | x |) \mathrm{d} x \\ = \left\{ \begin{array}{l l} O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right), & 0 <   a <   b, \\ O \left(\max \left\{\omega^ {- \alpha - \delta - \beta - 1} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right), & 0 = a <   b, \\ O \left(\max \left\{\frac {1}{\omega^ {\alpha + \delta + 1}}, \frac {\ln^ {\mu} \omega}{\omega^ {\min \left\{\beta + 3 / 2 , 3 / 2 \right\}}} \right\}\right), & a <   0 \leq b, \\ O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right), & a <   b <   0. \end{array} \right. \tag {16} \\ \end{array}
$$

• Suppose the following conditions hold: (i) $H 0 < a < b ,$ assume $\beta > - 1$ . (ii) If $a \ \leq \ 0 \ < \ b ,$ assume $\alpha + \delta + \nu > - 1$ and $\beta > - 1$ . (iii) $f \ * a \ < \ b \ = \ 0 .$ , assume $\alpha + \delta + \beta + \nu > - 1$ . (iv) ${ \cal I } f a < b < 0$ , assume $\beta > - 1$ . Then we have 

$$
\begin{array}{l} \int_ {a} ^ {b} x ^ {\alpha} | x | ^ {\delta} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \psi (x) J _ {v} (\omega | x |) \mathrm{d} x \\ = \left\{ \begin{array}{l l} O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right), & 0 <   a <   b, \\ O \left(\max \left\{\frac {1}{\omega^ {\alpha + \delta + 1}}, \frac {\ln^ {\mu} \omega}{\omega^ {\min \left\{\beta + 3 / 2 , 3 / 2 \right\}}} \right\}\right), & a \leq 0 <   b, \\ O \left(\max \left\{\omega^ {- \alpha - \delta - \beta - 1} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right), & a <   b = 0, \\ O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right), & a <   b <   0. \end{array} \right. \tag {17} \\ \end{array}
$$

Proof. The proof is divided into four cases: 

Case $( \mathbf { i } ) \colon 0 < a < b$ . According to Equation (14), we have 

$$
\int_ {a} ^ {b} x ^ {\alpha + \delta} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {\nu} (\omega x) \mathrm{d} x = O \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right),
$$

where we used the fact that $x ^ { \alpha + \delta } \psi ( x ) \in C [ a , b ]$ and $( x ^ { \alpha + \delta } \psi ( x ) ) ^ { \prime } \in L ^ { 1 } [ a , b ]$ 

Case (ii): $0 = a < b .$ . Applying Equation (10), we obtain 

$$
\int_ {0} ^ {b} x ^ {\alpha + \delta + \beta} \ln^ {\mu} (x) \psi (x) J _ {v} (\omega x) \mathrm{d} x = O \left(\max \left\{\omega^ {- \alpha - \delta - \beta - 1} \ln^ {\mu} \omega , \omega^ {- 3 / 2} \right\}\right).
$$

Case (iii): $a < 0 \leq b .$ By Equation (11), we estimate 

$$
\begin{array}{l} \int_ {a} ^ {b} x ^ {\alpha} | x | ^ {\delta} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {\nu} (\omega | x |) \mathrm{d} x \\ = (- 1) ^ {\delta} \int_ {a} ^ {0} x ^ {\alpha + \delta} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {v} (- \omega x) d x \\ + \int_ {0} ^ {b} x ^ {\alpha + \delta} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {\nu} (\omega x) d x \\ = (- 1) ^ {\alpha} \int_ {0} ^ {- a} x ^ {\alpha + \delta} \ln^ {\mu} (- a - x) (- a - x) ^ {\beta} \psi (- x) J _ {v} (\omega x) d x + O \left(\max \left\{\omega^ {- \alpha - \delta - 1}, \omega^ {- 3 / 2} \right\}\right) \\ = \mathcal {O} \left(\max \left\{\frac {1}{\omega^ {\alpha + \delta + 1}}, \frac {\ln^ {\mu} (\omega)}{\omega^ {\min \{\beta + 3 / 2 , 3 / 2 \}}} \right\}\right) + \mathcal {O} \left(\max \left\{\omega^ {- \alpha - \delta - 1}, \omega^ {- 3 / 2} \right\}\right) \\ = O \left(\max \left\{\frac {1}{\omega^ {\alpha + \delta + 1}}, \frac {\ln^ {\mu} \omega}{\omega^ {\min \{\beta + 3 / 2 , 3 / 2 \}}} \right\}\right), \\ \end{array}
$$

where we used the fact that $\ln ^ { \mu } ( x - a ) ( x - a ) ^ { \beta } \psi ( x ) \in C [ 0 , b ]$ and $( \ln ^ { \mu } ( x - a ) ( x - a ) ^ { \beta } \psi ( x ) ) ^ { \prime } \in$ $L ^ { 1 } [ 0 , b ]$ for $a < 0$ in the second equality. In particular, if $b = 0$ , the same estimate holds: 

$$
\begin{array}{l} \int_ {a} ^ {0} x ^ {\alpha} | x | ^ {\delta} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {v} (\omega | x |) d x \\ = (- 1) ^ {\alpha} \int_ {0} ^ {- a} x ^ {\alpha + \delta} \ln^ {\mu} (- a - x) (- a - x) ^ {\beta} \psi (- x) J _ {v} (\omega x) d x \\ = \mathcal {O} \left(\max \left\{\frac {1}{\omega^ {\alpha + \delta + 1}}, \frac {\ln^ {\mu} (\omega)}{\omega^ {\min \{\beta + 3 / 2 , 3 / 2 \}}} \right\}\right). \\ \end{array}
$$

Case (iv): $a < b < 0$ . Using Equation (15), we get 

$$
\begin{array}{l} \int_ {a} ^ {b} x ^ {\alpha} | x | ^ {\delta} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {v} (\omega | x |) \mathrm{d} x \\ = (- 1) ^ {\delta} \int_ {a} ^ {b} x ^ {\alpha + \delta} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) J _ {v} (- \omega x) d x \\ = (- 1) ^ {\alpha} \int_ {- b} ^ {- a} x ^ {\alpha + \delta} \ln^ {\mu} (- a - x) (- a - x) ^ {\beta} \psi (- x) J _ {v} (\omega x) \mathrm{d} x \\ = \mathcal {O} \left(\max \left\{\omega^ {- \beta - 3 / 2} \ln^ {\mu} (\omega), \omega^ {- 3 / 2} \right\}\right). \\ \end{array}
$$

The result in Equation (17) follows analogously by symmetry. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/33b7038c67e67de42bdab133589f98102389faddaf8f86e78442e6b4b39a862f.jpg)


With Lemma 3.2, Lemma 3.4 and Lemma 3.5 at disposal, we now establish a connection between the Bessel function and the Laguerre or Hermite polynomial by using the Hilb-type formula in (7). This connection enables us to derive the following asymptotic estimates presented in Theorems 3.1 and 3.2. 

Theorem 3.1. Let $\mu \in \mathbb { N } , \psi ( x ) \in C [ a , b ] , \psi ^ { \prime } ( x ) \in L ^ { 1 } [ a , b ]$ and $\alpha , \beta , \tau$ are selected such that the following integrals are integrable. Then the following asymptotic results hold for $n \to \infty$ . 

(i) $I f 0 = a < b < \infty ,$ , we have 

$$
\int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\tau} (b - x) ^ {\beta} \mathrm{e} ^ {- x} L _ {n} ^ {(\alpha)} (x) \psi (x) \mathrm{d} x = \mathcal {O} \left(\max \left\{\ln^ {\mu} (2 \sqrt {n}) n ^ {\alpha - \tau - 1}, n ^ {\max \left\{\frac {\alpha - \beta}{2} - \frac {3}{4}, \frac {\alpha}{2} - \frac {3}{4} \right\}} \right\}\right), \tag {18}
$$

and 

$$
\int_ {0} ^ {b} \ln^ {\mu} (b - x) x ^ {\tau} (b - x) ^ {\beta} \mathrm{e} ^ {- x} L _ {n} ^ {(\alpha)} (x) \psi (x) \mathrm{d} x = O \left(\max \left\{n ^ {\alpha - \tau - 1}, \ln^ {\mu} (2 \sqrt {n}) n ^ {\max \left\{\frac {\alpha - \beta}{2} - \frac {3}{4}, \frac {\alpha}{2} - \frac {3}{4} \right\}} \right\}\right). \tag {19}
$$

(ii) $I f 0 < a < b < \infty ,$ , we have 

$$
\int_ {a} ^ {b} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \psi (x) L _ {n} ^ {(\alpha)} (x) \omega_ {\alpha} (x) \mathrm{d} x = O \left(\max \left\{\ln^ {\mu} (2 \sqrt {n}) n ^ {\frac {\alpha - \beta}{2} - \frac {3}{4}}, n ^ {\frac {\alpha}{2} - \frac {3}{4}} \right\}\right), \tag {20}
$$

and 

$$
\int_ {a} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \psi (x) L _ {n} ^ {(\alpha)} (x) \omega_ {\alpha} (x) \mathrm{d} x = O \left(\max \left\{\ln^ {\mu} (2 \sqrt {n}) n ^ {\frac {\alpha - \beta}{2} - \frac {3}{4}}, n ^ {\frac {\alpha}{2} - \frac {3}{4}} \right\}\right). \tag {21}
$$

Proof. We begin by applying the substitution $x = t ^ { 2 }$ with $t > 0$ , it yields 

$$
\ln^ {\mu} x = \ln^ {\mu} t ^ {2} = 2 ^ {\mu} \ln^ {\mu} t.
$$

By the Hilb-type formula (7), the integral (18) becomes 

$$
\begin{array}{l} \int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\tau} (b - x) ^ {\beta} \mathrm{e} ^ {- x} L _ {n} ^ {(\alpha)} (x) \psi (x) \mathrm{d} x \\ = \frac {\Gamma (n + \alpha + 1)}{n ! \tilde {n} ^ {\frac {\alpha}{2}}} \int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\frac {2 \tau - \alpha}{2}} (b - x) ^ {\beta} J _ {\alpha} \left\{2 (\tilde {n} x) ^ {1 / 2} \right\} \mathrm{e} ^ {- \frac {x}{2}} \psi (x) \mathrm{d} x + O \left(n ^ {\frac {\alpha}{2} - \frac {3}{4}}\right) \\ = \frac {2 ^ {\mu + 1} \Gamma (n + \alpha + 1)}{n ! \tilde {n} ^ {\frac {\alpha}{2}}} \int_ {0} ^ {\sqrt {b}} \ln^ {\mu} (t) t ^ {2 \tau - \alpha + 1} (\sqrt {b} - t) ^ {\beta} J _ {\alpha} (2 \tilde {n} ^ {1 / 2} t) \Psi (t) \mathrm{d} t + O \left(n ^ {\frac {\alpha}{2} - \frac {3}{4}}\right), \\ \end{array}
$$

where $\tilde { n } = n + ( \alpha + 1 ) / 2$ and $\Psi ( t ) = \mathrm { e } ^ { - t ^ { 2 } / 2 } ( \sqrt { b } + x ) ^ { \beta } \psi ( t ^ { 2 } )$ . It is easy to verity that $\Psi ( t ) \in$ $C [ 0 , \sqrt { b } ]$ and $\Psi ( t ) \in L ^ { 1 } [ 0 , \sqrt { b } ]$ . Therefore, by Lemma 3.2, we have 

$$
\int_ {0} ^ {\sqrt {b}} \ln^ {\mu} (t) t ^ {2 \tau - \alpha + 1} (\sqrt {b} - t) ^ {\beta} J _ {\alpha} (2 \tilde {n} ^ {1 / 2} t) \Psi (t) \mathrm{d} t = O \left(\max \left\{\ln^ {\mu} (2 \sqrt {n}) n ^ {\alpha / 2 - \tau - 1}, n ^ {\max \left\{- \frac {\beta}{2} - \frac {3}{4}, - \frac {3}{4} \right\}} \right\}\right).
$$

In view of the asymptotic result [1, (6.1.2)] 

$$
\sigma_ {n} ^ {(\alpha)} = \frac {\Gamma (n + \alpha + 1)}{n !} = O (n ^ {\alpha}). \tag {22}
$$

Then we have 

$$
\int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\tau} (b - x) ^ {\beta} \mathrm{e} ^ {- x} L _ {n} ^ {(\alpha)} (x) \psi (x) \mathrm{d} x = O \left(\max \left\{\ln^ {\mu} (2 \sqrt {n}) n ^ {\alpha - \tau - 1}, n ^ {\max \left\{\frac {\alpha - \beta}{2} - \frac {3}{4}, \frac {\alpha}{2} - \frac {3}{4} \right\}} \right\}\right).
$$

By employing similar techniques, along with Lemmas 3.2 and 3.4, one can derive the remaining estimates (19), (20) and (21). 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/2d7dd5d982d06d0cfde99d4e0f608c85c9f3b199c553b899338cd97a3eca97cd.jpg)


Theorem 3.2. Let $\beta > - 1 , \mu \in \mathbb { N } , \psi ( x ) \in C [ a , b ]$ and $\psi ^ { \prime } ( x ) \in L ^ { 1 } [ a , b ]$ . Then, as $n  \infty$ , we have 

$$
\begin{array}{l} \int_ {a} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \mathrm{e} ^ {- x ^ {2}} H _ {n} (x) \psi (x) \mathrm{d} x \\ = \left\{ \begin{array}{l l} O \left(\ln^ {\mu} (2 \sqrt {n}) 2 ^ {n} \left(\frac {n}{2}\right)! \cdot n ^ {\max \left\{- \frac {\beta}{2} - 1, - 1 \right\}}\right), & n \text {   is   even } \\ O \left(\ln^ {\mu} (2 \sqrt {n}) 2 ^ {n} \left(\frac {n - 1}{2}\right)! \cdot n ^ {\max \left\{- \frac {\beta}{2} - \frac {1}{2}, - \frac {1}{2} \right\}}\right), & n \text {   is   odd } \end{array} \right. \tag {23} \\ \end{array}
$$

and 

$$
\begin{array}{l} \int_ {a} ^ {b} \ln^ {\mu} (x - a) (x - a) ^ {\beta} \mathrm{e} ^ {- x ^ {2}} H _ {n} (x) \psi (x) \mathrm{d} x \\ = \left\{ \begin{array}{l l} O \left(\ln^ {\mu} (2 \sqrt {n}) 2 ^ {n} \left(\frac {n}{2}\right)! \cdot n ^ {\max \left\{- \frac {\beta}{2} - 1, - 1 \right\}}\right), & n \text {   is   even } \\ O \left(\ln^ {\mu} (2 \sqrt {n}) 2 ^ {n} \left(\frac {n - 1}{2}\right)! \cdot n ^ {\max \left\{- \frac {\beta}{2} - \frac {1}{2}, - \frac {1}{2} \right\}}\right), & n \text {   is   odd } \end{array} . \right. \tag {24} \\ \end{array}
$$

Proof. Using Equation (6), the integral can be expressed as follows: 

$$
\begin{array}{l} \int_ {a} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \mathrm{e} ^ {- x ^ {2}} H _ {n} (x) \psi (x) \mathrm{d} x \\ = \left\{ \begin{array}{l l} (- 1) ^ {\frac {n}{2}} 2 ^ {n} \Big (\frac {n}{2} \Big)! \int_ {a} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \mathrm{e} ^ {- x ^ {2}} L _ {n / 2} ^ {(- 1 / 2)} (x ^ {2}) \psi (x) \mathrm{d} x, & n \text {is even}, \\ (- 1) ^ {\frac {n - 1}{2}} 2 ^ {n} \Big (\frac {n - 1}{2} \Big)! \int_ {a} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \mathrm{e} ^ {- x ^ {2}} x L _ {(n - 1) / 2} ^ {(1 / 2)} (x ^ {2}) \psi (x) \mathrm{d} x, & n \text {is odd}. \end{array} \right. \\ \end{array}
$$

For even values of n, applying Hilb-type asymptotic formula (7) and Equation (17), we obtain: 

$$
\begin{array}{l} \int_ {a} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \mathrm{e} ^ {- x ^ {2}} L _ {n / 2} ^ {(- 1 / 2)} (x ^ {2}) \psi (x) \mathrm{d} x \\ = \mathcal {O} (n ^ {- \frac {1}{4}}) \int_ {a} ^ {b} | x | ^ {\frac {1}{2}} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \mathrm{e} ^ {- \frac {x ^ {2}}{2}} \psi (x) J _ {- 1 / 2} (2 \bar {n} ^ {1 / 2} | x |) \mathrm{d} x + \mathcal {O} (n ^ {- 1}) \\ = O \left(\ln^ {\mu} (2 \sqrt {n}) n ^ {\max \left\{- \frac {\beta}{2} - 1, - 1 \right\}}\right), \\ \end{array}
$$

where $\bar { n } = n / 2 + 1 / 4$ . For odd values of n, we similarly obtain: 

$$
\begin{array}{l} \int_ {a} ^ {b} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \mathrm{e} ^ {- x ^ {2}} x L _ {(n - 1) / 2} ^ {(1 / 2)} (x ^ {2}) \psi (x) \mathrm{d} x \\ = O \left(n ^ {\frac {1}{4}}\right) \int_ {a} ^ {b} x | x | ^ {- \frac {1}{2}} \ln^ {\mu} (b - x) (b - x) ^ {\beta} \mathrm{e} ^ {- \frac {x ^ {2}}{2}} \psi (x) J _ {1 / 2} \left(2 \bar {n} ^ {1 / 2} | x |\right) \mathrm{d} x + O \left(n ^ {- \frac {1}{2}}\right) \\ = \mathcal {O} \left(\ln^ {\mu} (2 \sqrt {n}) n ^ {\max \left\{- \frac {\beta}{2} - \frac {1}{2}, - \frac {1}{2} \right\}}\right). \\ \end{array}
$$

Thus, the desired result in (23) is derived. A similar approach leads to the estimate in (24) by formula (7) and Equation (16).  

# 4. Asymptotics on Laguerre coefficients and convergence rates on Laguerre orthogonal projections for functions with algebraic and logarithmic regularities

# 4.1. Function with boundary regularities

Consider the function 

$$
f (x) = x ^ {\delta} \ln^ {\mu} (x) g (x), \tag {25}
$$

where $\mu \in \mathbb { N }$ and $g ( x ) \in C ^ { \infty } [ 0 , \infty )$ , such that $\int _ { 0 } ^ { \infty } \mathrm { e } ^ { - x } x ^ { \alpha } f ( x ) \mathrm { d } x < \infty$ with $\alpha > - 1$ . 

Theorem 4.1. Let $f ( x )$ be defined by (25). Suppose that $f ( x )$ satisfies the following conditions for $i = 0 , 1 , . . . , k$ with k being the smallest integer greater than or equal to $\alpha + 2 \delta + 1 / 2 \delta$ : 

$$
\lim _ {x \rightarrow \infty} \mathrm{e} ^ {- x / 2} x ^ {\alpha + i + 1} f ^ {(i)} (x) = 0, \int_ {b} ^ {\infty} \mathrm{e} ^ {- x / 2} x ^ {\alpha + k + 1} | f ^ {(k + 1)} (x) | \mathrm{d} x <   \infty , \tag {26}
$$

where $b > 0$ is a fixed positive constant. Then, for $\alpha + \delta > - 1$ and $n  \infty ,$ , the Laguerre coefficient $a _ { n } ( \alpha )$ given by (2) satisfies the following asymptotic estimate 

$$
\left| a _ {n} (\alpha) \right| = O \left(n ^ {- \alpha - \delta - 1} \ln^ {\mu} (2 \sqrt {n})\right). \tag {27}
$$

Proof. Recalling the Rodrigues formula for the generalized Laguerre polynomials (see, e.g., [21, (5.1.5)]) 

$$
\mathrm{e} ^ {- x} x ^ {\alpha} L _ {n} ^ {(\alpha)} (x) = \frac {1}{(n) _ {k}} \frac {\mathrm{d} ^ {k}}{\mathrm{d} x ^ {k}} \left\{\mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) \right\}, \tag {28}
$$

where $( n ) _ { k } = n ( n - 1 ) \cdots ( n - k + 1 )$ denotes the falling Pochhammer symbol. Then, the Laguerre coefficients can be expressed as 

$$
\begin{array}{l} a _ {n} (\alpha) = \frac {1}{\sigma_ {n} ^ {(\alpha)} n} \int_ {0} ^ {\infty} f (x) \mathrm{d} \left\{\mathrm{e} ^ {- x} x ^ {\alpha + 1} L _ {n - 1} ^ {(\alpha + 1)} (x) \right\} \\ = \frac {1}{\sigma_ {n} ^ {(\alpha)} n} \left(\mathrm{e} ^ {- x} x ^ {\alpha + 1} L _ {n - 1} ^ {(\alpha + 1)} (x) f (x) \Big | _ {0} ^ {\infty} - \int_ {0} ^ {\infty} f ^ {\prime} (x) \mathrm{e} ^ {- x} x ^ {\alpha + 1} L _ {n - 1} ^ {(\alpha + 1)} (x) \mathrm{d} x\right) \\ = \dots \tag {29} \\ = \frac {(- 1) ^ {k}}{\sigma_ {n} ^ {(\alpha)} (n) _ {k}} \int_ {0} ^ {\infty} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) f ^ {(k)} (x) \mathrm{d} x \\ = \frac {(- 1) ^ {k}}{\sigma_ {n} ^ {(\alpha)} (n) _ {k}} \left(\int_ {0} ^ {b} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) f ^ {(k)} (x) \mathrm{d} x + \int_ {b} ^ {\infty} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) f ^ {(k)} (x) \mathrm{d} x\right), \\ \end{array}
$$

where repeated integration by parts is applied, and the maximum estimate of the Laguerre polynomials (8), together with condition (26) ensures that 

$$
\lim _ {x \to + \infty} \mathrm{e} ^ {- x} x ^ {\alpha + i + 1} L _ {n - i - 1} ^ {(\alpha + i + 1)} (x) f ^ {(i)} (x) = 0, \quad i = 0, 1, \dots , k.
$$

If $\delta > 0$ is non-integer, applying Leibniz rule and Faà di Bruno’s formula yields 

$$
f ^ {(k)} (x) = \sum_ {i = 0} ^ {k} {\binom {k} {i}} (\delta) _ {i} x ^ {\delta - i} \sum_ {j = 0} ^ {k - i} {\binom {k - i} {j}} \Bigl (\ln^ {\mu} x \Bigr) ^ {(j)} g ^ {(k - i - j)} (x),
$$

which, for simplicity, can be written as 

$$
f ^ {(k)} (x) = x ^ {\delta - k} \sum_ {j = 0} ^ {\min \{\mu , k \}} h _ {\mu - j} (x) \ln^ {\mu - j} x,
$$

where $h _ { \mu - j } ( x ) \in C ^ { \infty } [ 0 , \infty )$ are smooth functions corresponding to the derivatives of $g ( x )$ . 

Using Equation (18), the finite integral part of $a _ { n } ( \alpha )$ can be estimated by 

$$
\begin{array}{l} \int_ {0} ^ {b} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) f ^ {(k)} (x) \mathrm{d} x \\ = \sum_ {j = 0} ^ {\min \{\mu , k \}} \int_ {0} ^ {b} \mathrm{e} ^ {- x} x ^ {\alpha + \delta} L _ {n - k} ^ {(\alpha + k)} (x) h _ {\mu - j} (x) \ln^ {\mu - j} x \mathrm{d} x \\ = \sum_ {j = 0} ^ {\min \{\mu , k \}} O \left(n ^ {k - \delta - 1} \ln^ {\mu - j} (2 \sqrt {n})\right) \\ = \mathcal {O} \left(n ^ {k - \delta - 1} \ln^ {\mu} (2 \sqrt {n})\right), \\ \end{array}
$$

where the estimate follows from the assumption $\begin{array} { r } { k - \delta - 1 \ge \frac { \alpha + k } { 2 } - \frac { 3 } { 4 } } \end{array}$ 2 . For the second term in Equation (29), integrating by parts once more gives 

$$
\begin{array}{l} \int_ {b} ^ {\infty} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) f ^ {(k)} (x) \mathrm{d} x \\ = - \frac {1}{n - k} \mathrm{e} ^ {- x} x ^ {\alpha + k + 1} L _ {n - k - 1} ^ {(\alpha + k + 1)} (x) f ^ {(k)} (x) \big | _ {x = b} \tag {30} \\ - \frac {1}{n - k} \int_ {b} ^ {\infty} \mathrm{e} ^ {- x} x ^ {\alpha + k + 1} L _ {n - k - 1} ^ {(\alpha + k + 1)} (x) f ^ {(k + 1)} (x) \mathrm{d} x. \\ \end{array}
$$

According to the estimate in Lemma 2.2, the integral in (30) satisfies 

$$
\int_ {b} ^ {\infty} \mathrm{e} ^ {- x} x ^ {\alpha + k + 1} L _ {n - k - 1} ^ {(\alpha + k + 1)} (x) f ^ {(k + 1)} (x) \mathrm{d} x = O \left(n ^ {\frac {\alpha + k}{2} + \frac {1}{4}}\right).
$$

Combining this with the asymptotic result (22), we conclude that 

$$
a _ {n} (\alpha) = \mathcal {O} \left(n ^ {- \alpha - \delta - 1} \ln^ {\mu} (2 \sqrt {n})\right).
$$

If $\delta$ is a nonnegative integer, we require an additional condition for $k ,$ , namely $k >$ max $\{ \delta , \mu \}$ , then 

$$
f ^ {(k)} (x) = \sum_ {i = 0} ^ {k} x ^ {\delta - k + i} \sum_ {j = 0} ^ {\mu - 1} \ln^ {\mu} (x) \phi_ {i, j} (x) + x ^ {\delta} \ln^ {\mu} (x) \phi_ {0} (x),
$$

where $\phi _ { 0 } ( x ) , \phi _ { i , j } ( x ) \in C ^ { \infty } [ 0 , + \infty )$ . By the similar above argument, we get the desired result (27). 

If $\delta < 0$ , by a similar proof as above, without integrating by parts, we can also obtain the desired result (27).  

Example 4.1. Figure 1 illustrates the optimal decay of the Laguerre expansion coefficients $| a _ { n } ( \alpha ) |$ for the function $f ( x ) = x ^ { \delta } \ln ^ { \mu } ( x )$ , with various values of $\alpha , \delta ,$ and $\mu .$ As observed in the figure, the asymptotic orders are consistent with the results established in Theorem 4.1. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/de236bf265930f9241cadaff9902c62d7b443e3d1f15faf20c2f0cb9f2db8562.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/ef47a2bc567cfde1998ca10c803ffea1a13f20d2e33b4b9f43923eeef3fbfaa7.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/b2aed349e969d4dce90369c5a8b6eddd7c36747532119e9edb441aee86f8deed.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/ca3640fb8e5cadc4bb03db23af5fb3b323d6c77adde0e3c11546afa58131037c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/9a1b547bd2cf26d6c02e003268c65943b596659b9ec9260b393fbd4691fe2670.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/13f5f9a9dfd4db07855151ef1f3eb0a006c2fd77a798d33377da6a2d1f7c1b09.jpg)



Figure 1: The asymptotic estimates of $\left| a _ { n } ( \alpha ) \right|$ for $f ( x ) = x ^ { \delta } \ln ^ { \mu } ( x ) \mathrm { : } \delta = 1 . 2 , \mu = 3$ (first row); δ = 3, µ = 3 (second row). In all figures, the blue dashed lines represent the asymptotic orders provided in Theorem 4.1.


# 4.2. Function with interior regularities

Consider the function 

$$
f (x) = | x - x _ {0} | ^ {\gamma} \ln^ {\mu} | x - x _ {0} | g (x), | x _ {0} | <   \infty , \tag {31}
$$

where $\mu \in \mathbb { N } , \gamma > 0 , g ( x ) \in C ^ { \infty } [ 0 , \infty )$ such that $f ( x )$ satisfies $\textstyle \int _ { 0 } ^ { \infty } \mathrm { e } ^ { - x } x ^ { \alpha } f ( x ) \mathrm { d } x < \infty$ for all $\alpha > - 1$ . 

Theorem 4.2. Let $f ( x )$ be defined by (31). Suppose that $f ( x )$ satisfies the following conditions for $i = 0 , 1 , \ldots , m$ with k being the largest integer less than or equal to γ 

$$
\lim _ {x \rightarrow \infty} \mathrm{e} ^ {- x / 2} x ^ {\alpha + i + 1} f ^ {(i)} (x) = 0, \int_ {b} ^ {\infty} \mathrm{e} ^ {- x / 2} x ^ {\alpha + k + 1} | f ^ {(k + 1)} (x) | \mathrm{d} x <   \infty ,
$$

where $b > x _ { 0 }$ is a fixed constant. Then, as n → ∞, the Laguerre coefficient (2) satisfies 

$$
\left| a _ {n} (\alpha) \right| = O \left(n ^ {- \frac {\alpha + \gamma}{2} - \frac {3}{4}} \ln^ {\mu} (2 \sqrt {n})\right). \tag {32}
$$

Proof. Observe that 

$$
a _ {n} (\alpha) = \frac {1}{\sigma_ {n} ^ {(\alpha)}} \left[ \int_ {0} ^ {x _ {0}} + \int_ {x _ {0}} ^ {+ \infty} \right] \mathrm{e} ^ {- x} x ^ {\alpha} f (x) L _ {n} ^ {(\alpha)} (x) \mathrm{d} x,
$$

then we estimate the asymptotic behavior of each integral term separately. 

Case (i): For $x \in [ 0 , x _ { 0 } ]$ , we have 

$$
f ^ {(k)} (x) = (x _ {0} - x) ^ {\gamma - k} \sum_ {j = 0} ^ {\mu} \ln^ {\mu - j} (x _ {0} - x) e _ {\mu - j} (x),
$$

where $e _ { \mu - j } ( x ) \in C ^ { \infty }$ and vanish when $j > k$ . Applying Rodrigues’ formula (28), we obtain 

$$
I _ {1} = \frac {1}{\sigma_ {n} ^ {(\alpha)}} \int_ {0} ^ {x _ {0}} \mathrm{e} ^ {- x} x ^ {\alpha} f (x) L _ {n} ^ {(\alpha)} (x) \mathrm{d} x = \frac {(- 1) ^ {k}}{\sigma_ {n} ^ {(\alpha)} (n) _ {k}} \int_ {0} ^ {x _ {0}} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) f ^ {(k)} (x) \mathrm{d} x. \tag {33}
$$

If γ is a positive integer, then $k = \gamma$ . According to the assumptions, we have $\textstyle { \frac { \alpha + \gamma } { 2 } } - { \frac { 3 } { 4 } } >$ −1. By using Equation (19), we get 

$$
\begin{array}{l} \int_ {0} ^ {x _ {0}} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) f ^ {(k)} (x) \mathrm{d} x \\ = \sum_ {j = 0} ^ {\mu} \int_ {0} ^ {x _ {0}} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) (x _ {0} - x) ^ {\gamma - k} \ln^ {\mu - j} (x _ {0} - x) e _ {\mu - j} (x) \mathrm{d} x \\ = \sum_ {j = 0} ^ {\mu} \mathcal {O} \left(n ^ {\frac {\alpha + \gamma}{2} - \frac {3}{4}} \ln^ {\mu - j} (2 \sqrt {n})\right) \\ = \mathcal {O} \left(n ^ {\frac {\alpha + \gamma}{2} - \frac {3}{4}} \ln^ {\mu} (2 \sqrt {n})\right). \\ \end{array}
$$

Combined with the normalization factor in (22), this yields 

$$
I _ {1} = \mathcal {O} \left(n ^ {- \frac {\alpha + \gamma}{2} - \frac {3}{4}} \ln^ {\mu} (2 \sqrt {n})\right). \tag {34}
$$

If γ is not an integer, then $k < \gamma < k + 1$ . By assumption, we have 

$$
\frac {\alpha + 2 k - \gamma}{2} - \frac {3}{4} > \frac {\alpha + k}{2} - \frac {3}{4} > - 1.
$$

Applying integration by parts to (33) yields 

$$
I _ {1} = \frac {(- 1) ^ {k + 1}}{\sigma_ {n} ^ {(\alpha)} (n) _ {k + 1}} \int_ {0} ^ {x _ {0}} e ^ {- x} x ^ {\alpha + k + 1} L _ {n - k - 1} ^ {(\alpha + k + 1)} (x) f ^ {(k + 1)} (x) d x,
$$

where 

$$
f ^ {(k + 1)} (x) = (x _ {0} - x) ^ {\gamma - k - 1} \sum_ {j = 0} ^ {\mu} \ln^ {\mu - j} (x _ {0} - x) \tilde {e} _ {\mu - j} (x),
$$

with $\tilde { e } _ { \mu - j } ( x ) \in C ^ { \infty } [ 0 , x _ { 0 } ]$ . Using estimate (19) again, we obtain 

$$
\begin{array}{l} \sum_ {j = 0} ^ {\mu} \int_ {0} ^ {x _ {0}} \mathrm{e} ^ {- x} x ^ {\alpha + k + 1} L _ {n - k - 1} ^ {(\alpha + k + 1)} (x) (x _ {0} - x) ^ {\gamma - k - 1} \ln^ {\mu - j} (x _ {0} - x) \tilde {e} _ {\mu - j} (x) \mathrm{d} x \\ = \sum_ {j = 0} ^ {\mu} O \left(n ^ {\frac {\alpha - \gamma}{2} + k + \frac {1}{4}} \ln^ {\mu - j} (2 \sqrt {n})\right) \\ = O \left(n ^ {\frac {\alpha - \gamma}{2} + k + \frac {1}{4}} \ln^ {\mu} (2 \sqrt {n})\right). \\ \end{array}
$$

Together with (22), it implies 

$$
I _ {1} = \mathcal {O} \left(n ^ {- \frac {\alpha + \gamma}{2} - \frac {3}{4}} \ln^ {\mu} (2 \sqrt {n})\right). \tag {35}
$$

Combining (34) and (35), we conclude that 

$$
\frac {1}{\sigma_ {n} ^ {(\alpha)}} \int_ {0} ^ {x _ {0}} \mathrm{e} ^ {- x} x ^ {\alpha} f (x) L _ {n} ^ {(\alpha)} (x) \mathrm{d} x = O \left(n ^ {- \frac {\alpha + \gamma}{2} - \frac {3}{4}} \ln^ {\mu} (2 \sqrt {n})\right). \tag {36}
$$

Case (ii): For $x \in [ x _ { 0 } , + \infty )$ , we write 

$$
f ^ {(k)} (x) = (x - x _ {0}) ^ {\gamma - k} \sum_ {j = 0} ^ {\mu} \ln^ {\mu - j} (x - x _ {0}) p _ {\mu - j} (x),
$$

where $p _ { \mu - j } \in C ^ { \infty } [ x _ { 0 } , \infty )$ and vanish when $j > k$ . By applying Rodrigues’ formula again, we obtain 

$$
I _ {2} = \frac {1}{\sigma_ {n} ^ {(\alpha)}} \int_ {x _ {0}} ^ {\infty} \mathrm{e} ^ {- x} x ^ {\alpha} f (x) L _ {n} ^ {(\alpha)} (x) \mathrm{d} x = \frac {(- 1) ^ {k}}{\sigma_ {n} ^ {(\alpha)} (n) _ {k}} \int_ {x _ {0}} ^ {\infty} \mathrm{e} ^ {- x} x ^ {\alpha + k} L _ {n - k} ^ {(\alpha + k)} (x) f ^ {(k)} (x) \mathrm{d} x.
$$

The asymptotic estimate of $I _ { 2 }$ follows by applying the same technique used in the proof of Theorem 4.1. 

Combining this with (36), we conclude that the asymptotic behavior of $a _ { n } ( \alpha )$ satisfies 

$$
a _ {n} (\alpha) = \mathcal {O} \left(n ^ {- \frac {\alpha + \gamma}{2} - \frac {3}{4}} \ln^ {\mu} (2 \sqrt {n})\right),
$$

which is the desired result stated in (32). 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/0d345711348be88dd8baf6793b8ea501c29eb5be366a5df70c27debf1f7ed57c.jpg)


Remark 4.1. For $- 1 < \gamma \leq 0$ and $\alpha > - 1$ , a similar argument to that used in Theorem 4.2, but without integration by parts, yields the following estimate 

$$
a _ {n} (\alpha) = \mathcal {O} \left(n ^ {- \frac {\alpha + \gamma}{2} - \frac {3}{4}} \ln^ {\mu} (2 \sqrt {n})\right).
$$

Example 4.2. Figure 2 illustrates the optimal decay of the Laguerre expansion coefficients $| a _ { n } ( \alpha ) |$ of the function $f ( x ) = | x - 0 . 3 | ^ { \gamma } \ln ^ { \mu } | x - 0 . 3 |$ for different values of $\alpha , \gamma ,$ and $\mu ,$ respectively. As observed from the graph, the asymptotic behavior of the coefficients agrees with the estimates established in Theorem 4.2. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/c9301f114141119a967f307e94a23b910386b4d1514100bfbea0a11e1fbad493.jpg)



n= 10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/1ddc878ed587e54e28c6ed170c35bf3dd2b81e339fce0646dd54bf4bfbe17e1d.jpg)



n = 10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/c11b0ab238f7248322cdb7b76ebcc075ac4e66c72ae9666839802dcd0c76eb9f.jpg)



n = 10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/de8fcfc36a25077e21e2d19cf093cefef6e3dc0f4494dcfb2b4ab204eddcc15d.jpg)



n =10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/abab25c1a856f549a35e1d4bf1458ac48056f1d6cb6a18bbd99c169f59cd0c19.jpg)



n=10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/d6d48f49dc8f3470a8be591776d2358d08f751c427908f61d9c820df19530006.jpg)



n=10:10:1000



Figure 2: The asymptotic estimates of $| a _ { n } ( \alpha ) |$ for $f ( x ) = | x - 0 . 3 | ^ { \gamma } \ln ^ { \mu } | x - 0 . 3 |$ with $\gamma = 1 . 2 , \mu = 2$ (first row) and $\gamma = 3 , \mu = 1$ (second row). In all figures, the blue dashed lines represent the asymptotic orders provided in Theorem 4.2.


# 4.3. The convergence rates on the Laguerre orthogonal projections

For the function $f ( x ) = x ^ { \delta } \ln ^ { \mu } ( x ) g ( x )$ , it is straightforward to verify that $f \in L _ { \omega _ { \alpha } } ^ { 2 } [ 0 , + \infty )$ if $\alpha + 2 \delta > - 1$ . Similarly, for the function $f ( x ) = | x - x _ { 0 } | ^ { \gamma } \ln ^ { \mu } | x - x _ { 0 } | g ( x )$ , it can be observed that $f \in L _ { \omega _ { \alpha } } ^ { 2 } [ 0 , + \infty )$ if $\gamma > - 1 / 2$ . By applying the asymptotic estimates for the Laguerre coefficients of functions with algebraic and logarithmic singularities, as established in Theorems 4.1 and 4.2, we will obtain the convergence rates of the Laguerre orthogonal projections for these respective functions. 

Theorem 4.3. Let $f ( x )$ be defined as in (25) and satisfy the assumptions in Theorem $4 . l .$ . Then for $\alpha + \delta > - 1$ , the Laguerre expansion follows as $N $ ∞ that 

$$
\begin{array}{l} \| f (x) - S _ {N} ^ {(\alpha)} [ f ] (x) \| _ {L _ {\omega_ {\alpha}} ^ {2} [ 0, + \infty)} = \mathcal {O} \left(N ^ {\frac {- \alpha - 2 \delta - 1}{2}} \ln^ {\mu} (2 \sqrt {N})\right), \quad \alpha + 2 \delta > - 1, \\ \left\| \mathrm{e} ^ {- x / 2} x ^ {\alpha / 2} (f - S _ {N} ^ {(\alpha)} [ f ]) \right\| _ {L ^ {\infty} (0, \infty)} = \mathcal {O} (N ^ {- \frac {\alpha}{2} - \delta - \frac {1}{4}} \ln^ {\mu} (2 \sqrt {N})), \quad \alpha + 2 \delta > - \frac {1}{2}. \\ \end{array}
$$

Proof. By Equation (3), we have 

$$
\| f (x) - S _ {N} ^ {(\alpha)} [ f ] (x) \| _ {L _ {\omega_ {\alpha}} ^ {2} [ 0, + \infty)} = \left[ \sum_ {n = N + 1} ^ {\infty} a _ {n} ^ {2} (\alpha) \sigma_ {n} ^ {(\alpha)} \right] ^ {\frac {1}{2}},
$$

which, together with the estimate of $a _ { n } ( \alpha )$ in Theorem 4.1. 

Similarly, by substituting (3), we get 

$$
\mathrm{e} ^ {- x / 2} x ^ {\alpha / 2} \left(f - S _ {N} ^ {(\alpha)} [ f ]\right) = \sum_ {n = N + 1} ^ {\infty} a _ {n} (\alpha) \left(\mathrm{e} ^ {- x / 2} x ^ {\alpha / 2} L _ {n} ^ {(\alpha)}\right),
$$

from which the second estimate follows directly by applying Lemma 2.2, Theorem 4.1, and integration by parts.  

Example 4.3. Figure 3 shows the convergence of the Laguerre orthogonal projection $\parallel f -$ $S _ { N } ^ { ( \alpha ) } [ f ] \vert _ { \mathrm { L } _ { w _ { \alpha } } ^ { 2 } [ 0 , + \infty ) }$ for the function $f ( x ) \ = \ x ^ { \delta } \ln ^ { \mu } ( x )$ with different values of $\alpha , \delta ,$ , and $\mu ,$ respectively. As can be seen from the figure, these asymptotic orders are consistent with the results in Theorem 4.3. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/7fbc30f10f2ca85e4cdae89f8c04175794a3973a11e28230782a5c694a5d1603.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/159de4ccae3fe1ccd4f6ba408cc6c331d11d1b48f9fb5d7ca842c90d561e9bc6.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/90d2554ec7e8914769617f79162c78a20782adc13137641145326bc8ea175807.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/3f1aa320ceaeb2a66d1f4cc70dd25681fffb9b591ba3e91adabb93bf1ab51bbc.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/e3ebef9c6cab2ee5df799253f232109a99232471d1414d0bb4ee093b78ef74ff.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/8d27a08513367088bae4f3d3a93cce8f4a892789559418be36639a74efcd15fa.jpg)



Figure 3: The asymptotic estimates of $\lVert f - S _ { N } ^ { ( \alpha ) } [ f ] \rVert _ { \mathrm { L } _ { w _ { \alpha } } ^ { 2 } [ 0 , + \infty ) }$ for $f ( x ) = x ^ { \delta }$ lnµ(x): δ = 1.2, µ = 3 (first row); $\delta = 4 , \mu = 1$ α (second row). In all figures, the blue dashed lines represent the asymptotic orders provided in Theorem 4.3.


Theorem 4.4. Let f (x) be defined by (31) and satisfy the assumptions in Theorem 4.2. Then, as $N  \infty$ , the Laguerre expansion follows that 

$$
\begin{array}{l} \| f (x) - S _ {N} ^ {(\alpha)} [ f ] (x) \| _ {L _ {\omega_ {\alpha}} ^ {2} [ 0, + \infty)} = \mathcal {O} \left(N ^ {- \frac {\gamma}{2} - \frac {1}{4}} \ln^ {\mu} (2 \sqrt {N})\right), \quad \gamma > - 1 / 2, \\ \left\| \mathrm{e} ^ {- x / 2} x ^ {\alpha / 2} (f - S _ {N} ^ {(\alpha)} [ f ]) \right\| _ {L ^ {\infty} (0, \infty)} = O \Big (N ^ {- \frac {\gamma}{2}} \ln^ {\mu} (2 \sqrt {N}) \Big), \quad \gamma > 0. \\ \end{array}
$$

Proof. By combining Theorem 4.2, the results can be obtained using the same proof method as in Theorem 4.3.  

Example 4.4. Figure 4 shows the convergence of the Laguerre orthogonal projection $\parallel f -$ $S _ { N } ^ { ( \alpha ) } [ f ] \vert _ { \mathrm { L } _ { w _ { \alpha } } ^ { 2 } [ 0 , + \infty ) }$ for the function $f ( x ) = | x - 0 . 3 | ^ { \gamma } \ln ^ { \mu } | x - 0 . 3 |$ with different values of α, γ, and $\mu ,$ respectively. As can be seen from the figure, these asymptotic orders are consistent with the results in Theorem 4.4. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/707c7480d325376131aaea4bbb3ead1d6983d11f2ae8c971d8e7313bdd15dd02.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/fb68f2a35c58ac98320a9e170d32630c09268f51bfeda78a279d93e7e0996fa9.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/f2c835b7cfaf5950eb9cf24fa9500488f655838cb2342d2d887518e12fa7aa49.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/b2bceef39d127cd4b9f7f75b593317ad125de64f0bb63addc068e9de9c971537.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/0e5c15c871d62724597a23e33609a5050274f09602b9aacf7f284c088d88d6ec.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/62494fd73d21df6a372ead8c87569805b5f314b2be57fcf87b04436c90267e91.jpg)



Figure 4: The asymptotic estimates of $\lVert f - S _ { N } ^ { ( \alpha ) } [ f ] \rVert _ { \mathrm { L } _ { w _ { \alpha } } ^ { 2 } [ 0 , + \infty ) }$ for $f ( x ) = | x - 0 . 3 | ^ { \gamma } \ln ^ { \mu } | x - 0 . 3 | \colon \gamma = 1 . 2 , \mu = 2$ (first row); $\gamma = 3 , \mu = 1$ α (second row). In all figures, the blue dashed lines represent the asymptotic orders provided in Theorem 4.4.


Moreover, for the non-uniformly Laguerre-weighted Sobolev space $H ^ { m , \alpha } ( \Omega ) \left[ 1 2 \right]$ with any integer $m \ge 0 , \alpha > - 1 , \Omega = ( 0 , + \infty )$ , the weighted norm of $H ^ { m , \alpha } ( \Omega )$ is defined by 

$$
\| u \| _ {H ^ {m, \alpha} (\Omega)} = \left\{\sum_ {q = 0} ^ {m} \int_ {0} ^ {+ \infty} \mathrm{e} ^ {- x} x ^ {\alpha + q} [ u ^ {(q)} (x) ] ^ {2} \mathrm{d} x \right\} ^ {\frac {1}{2}},
$$

define 

$$
f ^ {(q)} (x) = \sum_ {n = 0} ^ {\infty} a _ {n} ^ {(q)} (\alpha + q) L _ {n} ^ {(\alpha + q)} (x), q = 0, 1, \dots , m
$$

with 

$$
a _ {n} ^ {(q)} (\alpha + q) = \frac {1}{\sigma_ {n} ^ {(\alpha + q)}} \int_ {0} ^ {+ \infty} \mathrm{e} ^ {- x} x ^ {\alpha + q} f ^ {(q)} (x) L _ {n} ^ {(\alpha + q)} (x) \mathrm{d} x.
$$

From (1) and (29), we get 

$$
\begin{array}{l} a _ {n} ^ {(q)} (\alpha + q) = \frac {1}{\sigma_ {n} ^ {(\alpha + q)}} \int_ {0} ^ {+ \infty} \mathrm{e} ^ {- x} x ^ {\alpha + q} f ^ {(q)} (x) L _ {n} ^ {(\alpha + q)} (x) \mathrm{d} x \\ = \frac {(- 1) ^ {q}}{\sigma_ {n} ^ {(\alpha + q)}} (n + q) (n + q - 1) \dots (n + 1) \int_ {0} ^ {+ \infty} \mathrm{e} ^ {- x} x ^ {\alpha} f (x) L _ {n + q} ^ {(\alpha)} (x) \mathrm{d} x \\ = \frac {(- 1) ^ {q} \sigma_ {n + q} ^ {(\alpha)}}{\sigma_ {n} ^ {(\alpha + q)}} (n + q) (n + q - 1) \dots (n + 1) a _ {n + q} (\alpha) \\ = a _ {n + q} (\alpha) \mathcal {O} (1), \\ \end{array}
$$

which together with Theorems 4.1 and 4.2, we can obtain the following convergence rates. 

For the function $f ( x ) = x ^ { \delta } \ln ^ { \mu } ( x ) g ( x )$ , it is easy to get $f \in H ^ { m , \alpha } ( \Omega )$ if $\alpha + 2 \delta > m - 1$ . For the function $f ( x ) = | x - x _ { 0 } | ^ { \gamma } \ln ^ { \mu } | x - x _ { 0 } | g ( x )$ , it is realize that $f \in H ^ { m , \alpha } ( \Omega ) \mathrm { i f } \gamma > m - \frac { 1 } { 2 }$ 

Corollary 4.1. Let f (x) be defined by (25) and satisfy the assumptions in Theorem 4.1. Then for $\alpha + \delta > - 1$ and $N \gg 1$ , the Laguerre expansion follows that 

$$
\| f (x) - S _ {N} ^ {(\alpha)} [ f ] (x) \| _ {H ^ {m, \alpha} (\Omega)} = O \left(N ^ {\frac {m - \alpha - 2 \delta - 1}{2}} \ln^ {\mu} (2 \sqrt {N})\right), \alpha + 2 \delta > m - 1.
$$

Corollary 4.2. Let f (x) be defined by (31) and satisfy the assumptions in Theorem 4.2. Then for $N \gg 1$ , the Laguerre expansion follows that 

$$
\| f (x) - S _ {N} ^ {(\alpha)} [ f ] (x) \| _ {H ^ {m, \alpha} (\Omega)} = \mathcal {O} \left(N ^ {\frac {m - \gamma}{2} - \frac {1}{4}} \ln^ {\mu} (2 \sqrt {N})\right), \gamma > m - \frac {1}{2}.
$$

# 5. Asymptotics on Hermite coefficients and convergence rates on Hermite orthogonal projections for the functions with algebraic and logarithmic regularities

# 5.1. Function with interior regularities

Consider the function 

$$
f (x) = \left| x - z _ {0} \right| ^ {s} \ln^ {\mu} | x - z _ {0} | g (x), z _ {0} \in (- \infty , + \infty), \tag {37}
$$

where $\mu$ is a positive integer, $s > 0 , g ( x ) \in C ^ { \infty } ( - \infty , + \infty )$ such that $\textstyle \int _ { - \infty } ^ { + \infty } \mathbf { e } ^ { - x ^ { 2 } } f ( x ) \mathrm { d } x < \infty$ 

Theorem 5.1. Let f (x) be defined by (37). Suppose that f (x) satisfies the following conditions for $i = 0 , 1 , \ldots , k ,$ , with k being largest integer less than or equal to s: 

$$
\begin{array}{l} \lim _ {x \to - \infty} \mathrm{e} ^ {- \frac {x ^ {2}}{2}} f ^ {(i)} (x) = 0, \lim _ {x \to + \infty} \mathrm{e} ^ {- \frac {x ^ {2}}{2}} f ^ {(i)} (x) = 0, \\ \int_ {- \infty} ^ {- b} \mathrm{e} ^ {- \frac {x ^ {2}}{2}} | f ^ {(k + 1)} (x) | \mathrm{d} x <   \infty , \int_ {b} ^ {\infty} \mathrm{e} ^ {- \frac {x ^ {2}}{2}} | f ^ {(k + 1)} (x) | \mathrm{d} x <   \infty , \\ \end{array}
$$

where b is a positive constant and $b > \left| \boldsymbol { z } _ { 0 } \right|$ . Then, as $n  \infty ,$ , the Hermite coefficients in (4) satisfy 

$$
\left| h _ {n} \right| = O \left(n ^ {- \frac {n + s}{2} - 1} \ln^ {\mu} (2 \sqrt {n})\right). \tag {38}
$$

Proof. We begin by splitting $h _ { n }$ as follows 

$$
h _ {n} = \frac {1}{\gamma_ {n}} \left[ \int_ {- \infty} ^ {z _ {0}} + \int_ {z _ {0}} ^ {+ \infty} \right] \mathrm{e} ^ {- x ^ {2}} f (x) H _ {n} (x) \mathrm{d} x =: I _ {1} + I _ {2}.
$$

We now estimate the asymptotic behavior of these two integrals. 

Case (i): For $x \in ( - \infty , z _ { 0 } ]$ , the k-th derivative of $f ( x )$ can be written as 

$$
f ^ {(k)} (x) = (z _ {0} - x) ^ {s - k} \sum_ {j = 0} ^ {\mu} \ln^ {\mu - j} (z _ {0} - x) q _ {\mu - j} (x),
$$

where $q _ { \mu - j } ( x ) \in C ^ { \infty } ( - \infty , z _ { 0 } ]$ . Recalling the Rodrigues’ formula for Hermite polynomials 

$$
\mathrm{e} ^ {- x ^ {2}} H _ {n} (x) = (- 1) ^ {k} \frac {\mathrm{d} ^ {k}}{\mathrm{d} x ^ {k}} \left\{\mathrm{e} ^ {- x ^ {2}} H _ {n - k} (x) \right\},
$$

and performing repeated integration by parts, we obtain 

$$
\begin{array}{l} I _ {1} = \frac {1}{\gamma_ {n}} \int_ {- \infty} ^ {z _ {0}} f (x) \mathrm{e} ^ {- x ^ {2}} H _ {n} (x) \mathrm{d} x \\ = \frac {1}{\gamma_ {n}} \int_ {- \infty} ^ {z _ {0}} f ^ {\prime} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - 1} (x) \mathrm{d} x \tag {39} \\ = \dots \\ = \frac {1}{\gamma_ {n}} \int_ {- \infty} ^ {z _ {0}} f ^ {(k)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k} (x) \mathrm{d} x. \\ \end{array}
$$

If $s \in \mathbb { N } ^ { + }$ , then $k = s .$ . The above integral can be split as 

$$
\int_ {- \infty} ^ {z _ {0}} f ^ {(k)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k} (x) \mathrm{d} x = \int_ {- b} ^ {z _ {0}} f ^ {(k)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k} (x) \mathrm{d} x + \int_ {- \infty} ^ {- b} f ^ {(k)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k} (x) \mathrm{d} x, \tag {40}
$$

where $b > \left| z _ { 0 } \right|$ is a constant. Substituting the expression for $f ^ { ( k ) } ( x )$ , the first integral in 

Equation (40) can be estimated as 

$$
\begin{array}{l} \int_ {- b} ^ {z _ {0}} f ^ {(k)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k} (x) \mathrm{d} x \\ = \sum_ {j = 0} ^ {\mu} \int_ {- b} ^ {z _ {0}} \ln^ {\mu - j} (z _ {0} - x) e ^ {- x ^ {2}} H _ {n - k} (x) q _ {\mu - j} (x) d x \\ = \left\{ \begin{array}{l l} O \bigg (2 ^ {n} \big (\frac {n - k}{2} \big)! n ^ {- 1} \ln^ {\mu} (2 \sqrt {n}) \bigg), & n - k \text {is even}, \\ O \bigg (2 ^ {n} \big (\frac {n - k - 1}{2} \big)! n ^ {- \frac {1}{2}} \ln^ {\mu} (2 \sqrt {n}) \bigg), & n - k \text {is odd}. \end{array} \right. \\ \end{array}
$$

For the second integral in (40), we apply integration by parts and obtain 

$$
\int_ {- \infty} ^ {- b} f ^ {(k)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k} (x) \mathrm{d} x \tag {41}
$$

$$
= - f ^ {(k)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \big | _ {x = - b} + \int_ {- \infty} ^ {- b} f ^ {(k + 1)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \mathrm{d} x.
$$

According to Equation (9), the two terms on the right-hand side of Equation (41) can be estimated as 

$$
\left| f ^ {(k)} (b) \mathrm{e} ^ {- b ^ {2}} H _ {n - k - 1} (b) \right| = O \left(\left(\frac {2}{\mathrm{e}}\right) ^ {\frac {n}{2}} n ^ {\frac {n - k}{2} - \frac {1}{2}}\right),
$$

and 

$$
\int_ {- \infty} ^ {- b} f ^ {(k + 1)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \mathrm{d} x = O \left(\left(\frac {2}{\mathrm{e}}\right) ^ {\frac {n}{2}} n ^ {\frac {n - k}{2} - \frac {1}{2}}\right).
$$

Using Stirling’s formula, and applying Equation (4), we obtain 

$$
I _ {1} = \mathcal {O} \left(n ^ {- \frac {n + s}{2} - 1} \ln^ {\mu} (2 \sqrt {n})\right).
$$

If $s > 0$ is not an integer, we apply integration by parts again to $I _ { 1 } .$ , yielding 

$$
I _ {1} = \frac {1}{\gamma_ {n}} \int_ {- \infty} ^ {z _ {0}} f ^ {(k + 1)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \mathrm{d} x.
$$

Here, $f ^ { ( k + 1 ) } ( x )$ can be expressed as 

$$
f ^ {(k + 1)} (x) = (z _ {0} - x) ^ {s - k - 1} \sum_ {j = 0} ^ {\mu} \ln^ {\mu - j} (z _ {0} - x) \tilde {q} _ {\mu - j} (x),
$$

where $\tilde { q } _ { \mu - j } \in C ^ { \infty } ( - \infty , z _ { 0 } ]$ . The integral can then be split again as 

$$
\begin{array}{l} \int_ {- \infty} ^ {z _ {0}} f ^ {(k + 1)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \mathrm{d} x \tag {42} \\ = \int_ {- b} ^ {z _ {0}} f ^ {(k + 1)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \mathrm{d} x + \int_ {- \infty} ^ {- b} f ^ {(k + 1)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \mathrm{d} x. \\ \end{array}
$$

Substituting the expression for $f ^ { ( k + 1 ) } ( x )$ , we obtain 

$$
\begin{array}{l} \int_ {- b} ^ {z _ {0}} f ^ {(k + 1)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \mathrm{d} x \\ = \sum_ {j = 0} ^ {\mu} \int_ {- b} ^ {z _ {0}} (z _ {0} - x) ^ {s - k - 1} \ln^ {\mu - j} (z _ {0} - x) e ^ {- x ^ {2}} H _ {n - k - 1} (x) \tilde {q} _ {\mu - j} (x) d x \\ = \left\{ \begin{array}{l l} O \bigg (2 ^ {n} \Big (\frac {n - k - 1}{2} \Big)! n ^ {\frac {k - s - 1}{2}} \ln^ {\mu} (2   \sqrt {n}) \bigg), & n - k - 1 \text { is even}, \\ O \bigg (2 ^ {n} \Big (\frac {n - k - 2}{2} \Big)! n ^ {\frac {k - s}{2}} \ln^ {\mu} (2   \sqrt {n}) \bigg), & n - k - 1 \text { is odd}. \end{array} \right. \\ \end{array}
$$

Similar to Equation (41), the second integral in (42) can be estimated as 

$$
\int_ {- \infty} ^ {- b} f ^ {(k + 1)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k - 1} (x) \mathrm{d} x = O \left(\left(\frac {2}{\mathrm{e}}\right) ^ {\frac {n}{2}} n ^ {\frac {n - k}{2} - 1}\right).
$$

Using the Stirling formula and Equation (4), we conclude that 

$$
I _ {1} = \mathcal {O} \left(n ^ {- \frac {n + s}{2} - 1} \ln^ {\mu} (2 \sqrt {n})\right).
$$

Case (ii): For x  [z0, + ), then 

$$
f ^ {(k)} (x) = (x - z _ {0}) ^ {s - k} \sum_ {j = 0} ^ {\mu} \ln^ {\mu - j} (x - z _ {0}) r _ {\mu - j} (x),
$$

where $r _ { \mu - j } ( x ) \in C ^ { \infty } [ z _ { 0 } , \infty )$ . By Rodrigues’ formula, we get 

$$
I _ {2} = \frac {1}{\gamma_ {n}} \int_ {z _ {0}} ^ {\infty} f ^ {(k)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n - k} (x) \mathrm{d} x.
$$

Following a similar approach to Case (i), we derive the desired result in Equation (38).  

Example 5.1. Figure 5 illustrates the decay rates of $\log _ { n } | h _ { n } |$ for the function $f ( x ) \ =$ $| x - z _ { 0 } | ^ { s } \ln ^ { \mu } | x - z _ { 0 } |$ with various values of s and $\mu .$ As shown in the graph, the observed asymptotic behavior aligns with the predictions made in Theorem 5.1. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/508f87e3a51ab86c09a05b09bbbfb40ec057bd9de688cbdcc53aefe952d5450c.jpg)



n = 10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/549a049efbca1bb5e35e132ae810683859ad3eeab29ab50fd03bd007ac051457.jpg)



= 10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/6a13bd9343c30fba8c3e83549f6be4f71db7d11a17b2fc64726e1dae22297e49.jpg)



n = 10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/3ee856a72a83a153ba571831777348eced62445c4c3f237363c26fca14330a1a.jpg)



n=10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/5715648c48a2645a746d803b4ddb67c1bfd2e6bbe0667189aa6d1dfa4dfd6be2.jpg)



n =10:10:1000


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/e9a55ccb216e7ad26bddc79b7cae0f14ea6ae1e2737144fc8bf4d5884a15aadf.jpg)



n=10:10:1000



Figure 5: The asymptotic estimates of $\log _ { n } | h _ { n } |$ for $f ( x ) = | x - 3 | ^ { s } \ln ^ { \mu } | x - 3 |$ . In all figures, the blue dashed lines represent the asymptotic orders provided in Theorem 5.1.


# 5.2. The convergence rates on the Hermite orthogonal projections

The asymptotic behavior of the Hermite spectral expansion coefficients for functions with algebraic and logarithmic regularities at the interior, as described in Theorem 5.1, enables us to determine the convergence rate of the Hermite orthogonal projection for a given function. 

Theorem 5.2. Let $f ( x )$ be defined by (37) and satisfies the assumptions given in Theorem 5.1. Then, as $N \gg 1$ , the Hermite expansion satisfies that 

$$
\left\| f (x) - S _ {N} [ f ] (x) \right\| _ {L _ {\omega} ^ {2} (R)} = \mathcal {O} \left(N ^ {- \frac {s}{2} - \frac {1}{4}} \ln^ {\mu} (2 \sqrt {N})\right).
$$

$$
\left\| \mathrm{e} ^ {- x ^ {2} / 2} (f - S _ {N} [ f ]) \right\| _ {L ^ {\infty} (R)} = O \left(N ^ {- \frac {s}{2}} \ln^ {\mu} (2 \sqrt {N})\right).
$$

Proof. By Equation (5), we have 

$$
\| f (x) - S _ {N} [ f ] (x) \| _ {L _ {\omega} ^ {2} (R)} = \left[ \sum_ {n = N + 1} ^ {\infty} h _ {n} ^ {2} \gamma_ {n} \right] ^ {\frac {1}{2}},
$$

which directly leads to the desired result by Theorem 5.1 and the method of integration by parts. 

Similarly, from (5), we obtain 

$$
\mathrm{e} ^ {- x ^ {2} / 2} (f - S _ {N} [ f ]) = \sum_ {n = N + 1} ^ {\infty} h _ {n} \left(\mathrm{e} ^ {- x ^ {2} / 2} H _ {n}\right),
$$

which again leads to the desired result by applying Lemma 2.3, Theorem 5.1, and integration by parts.  

Example 5.2. Figure 6 shows the convergence rates of the weighted Hermite truncation errors $\lVert f ( x ) - S _ { N } [ f ] ( x ) \rVert _ { L _ { \omega } ^ { 2 } ( R ) }$ for the function $f ( x ) = | x - z _ { 0 } | ^ { s } \ln ^ { \mu } | x - z _ { 0 } |$ with various values of s and $\mu .$ As can be seen from the graph, these asymptotic orders are consistent with the results of Theorem 5.2. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/ed0b3d65a60c09fac021a050edb696f45741203f96a8b2ba8015949a43999f50.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/6a16ef0a9e643170117efa9fa0375b197e73f3ab3be2d4d3c6d71da25a0e25e1.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/9565aa1f9b5546da711e8cc462f411155471b35ea1ddc2cafa9631c534ae0747.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/fde5fa2ae31f72417e0b0869ff6ebeb234619530b9bca679b34cfd17c08087cd.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/310874383818426b6b20e8e55bdf7d15671dce9175640a6c95aac864aa26774c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/b045178e-6d67-421f-a222-cbee7803e1ab/f11db98d62a5f4b865987eee44babecd88e8d8244f67454f5af2504db29aea8f.jpg)



Figure 6: The asymptotic estimates of $\lVert f ( x ) - S _ { N } [ f ] ( x ) \rVert _ { L _ { \omega } ^ { 2 } ( R ) }$ for $f ( x ) = | x - 3 | ^ { s } \ln ^ { \mu } | x - 3 |$ . In all figures, the blue dashed lines represent the asymptotic orders provided in Theorem 5.2.


Moreover, for the Hermite-weighted Sobolev space $W ^ { m } ( R ) \left[ 5 \right]$ with any integer $m \geq 0$ , the weighted norm of $W ^ { m } ( R )$ is defined by 

$$
\| v \| _ {W ^ {m} (R)} = \left\{\sum_ {p = 0} ^ {m} \int_ {- \infty} ^ {+ \infty} \mathrm{e} ^ {- x ^ {2}} [ v ^ {(p)} (x) ] ^ {2} \mathrm{d} x \right\} ^ {\frac {1}{2}},
$$

define 

$$
f ^ {(p)} (x) = \sum_ {n = 0} ^ {\infty} h _ {n} ^ {(p)} H _ {n} (x), p = 0, 1,..., m
$$

with 

$$
h _ {n} ^ {(p)} = \frac {1}{\gamma_ {n}} \int_ {- \infty} ^ {+ \infty} f ^ {(p)} (x) H _ {n} (x) \mathrm{e} ^ {- x ^ {2}} \mathrm{d} x.
$$

From (4) and (39), we get 

$$
\begin{array}{l} h _ {n} ^ {(p)} = \frac {1}{\gamma_ {n}} \int_ {- \infty} ^ {+ \infty} f ^ {(p)} (x) \mathrm{e} ^ {- x ^ {2}} H _ {n} (x) \mathrm{d} x \\ = \frac {1}{\gamma_ {n}} \int_ {- \infty} ^ {+ \infty} f (x) \mathrm{e} ^ {- x ^ {2}} H _ {n + p} (x) \mathrm{d} x \\ = h _ {n + p} \mathcal {O} (1), \\ \end{array}
$$

which together with Theorems 5.1 and 5.2, we can obtain the following convergence rates. 

Corollary 5.1. Let f (x) be defined by (37) and satisfy the assumptions in Theorem 5.1. Then for N  1, the Hermite expansion follows that 

$$
\left\| f (x) - S _ {N} [ f ] (x) \right\| _ {W ^ {m} (R)} = O \left(N ^ {\frac {m - s}{2} - \frac {1}{4}} \ln^ {\mu} (2 \sqrt {N})\right).
$$

# 6. Conclusions

This paper investigates the optimal asymptotic behavior of the Laguerre and Hermite spectral expansion coefficients for functions with algebraic and logarithmic regularities. By applying the Hilb-type formula and deriving precise asymptotic estimates for integrals involving Laguerre and Hermite polynomials, we establish the optimal decay rates of the coefficients. These results are then used to characterize the decay of the corresponding Laguerre and Hermite spectral orthogonal projections. Numerical experiments are performed to validate the theoretical findings, demonstrating their accuracy and applicability. These results provide valuable insights into the spectral approximation of functions with regularities, paving the way for more efficient numerical methods in computational mathematics. 

# Conflict of Interest

The authors declare that they have no conflict of interest. 

# Data Availability

The code used in this work will be made available upon request to the authors. 

# References

# References



[1] M. Abramowitz, I.A. Stegun, Handbook of mathematical functions with formulas, graphs, and mathematical tables, Dover Publications, Inc., New York, 1966. 





[2] P. Agarwal, F. Qi, M. Chand, S. Jain, Certain integrals involving the generalized hypergeometric function and the Laguerre polynomials, J. Comput. Appl. Math. 313 (2017) 307–317 . https://doi.org/10.1016/j.cam.2016.09.034 





[3] Z. Avazzadeh, H. Hassani, P. Agarwal, S. Mehrabi, M.J. Ebadi, M.S. Dahaghin, An optimization method for studying fractional-order tuberculosis disease model via generalized Laguerre polynomials, Soft comput. 27(14) (2023) 9519–9531 . https://doi.org/10.1007/s00500-023-08086-z 





[4] Z. Avazzadeh, H. Hassani, M.J. Ebadi, P. Agarwal, M. Poursadeghfard, E. Naraghirad, Optimal Approximation of Fractional Order Brain Tumor Model Using Generalized Laguerre Polynomials, Iran. J. Sci. 47(2) (2023) 501–513 . https://doi.org/10.1007/s40995-022-01388-1 





[5] H. Brezis, Functional analysis, Sobolev spaces and partial differential equations, Vol. 2. Springer, 2011. 





[6] S. Chen, J. Shen, Log orthogonal functions in semi-infinite intervals: approximation results and applications, SIAM J. Numer. Anal. 61(1) (2023) 110–134. https://doi.org/10.1137/21M1466840 





[7] D. Colton, J. Wimp, Analytic solutions of the heat equation and some formulas for Laguerre and Hermite polynomials, Complex Var. Elliptic Equations 3(4) (1984) 397– 412 . https://doi.org/10.1080/17476938408814079 





[8] D. Elliott, P. Tuan, Asymptotic estimates of Fourier coefficients, SIAM J. Math. Anal. 5(1) (1974) 1–10. https://doi.org/10.1137/0505001 





[9] D. Funaro, O. Kavian, Approximation of some diffusion evolution equations in unbounded domains by Hermite functions, Math. Comp. 57(196) (1991) 597–619. https://doi.org/10.1090/S0025-5718-1991-1094949-X 





[10] A.G. Fatyanov, A.V. Terekhov, High-performance modeling acoustic and elastic waves using the Parallel Dichotomy Algorithm, J. Comput. Phys. 230(5) (2011) 1992– 2003. https://doi.org/10.1016/j.jcp.2010.11.046 





[11] A. Gil, J. Segura, N.M. Temme, Asymptotic approximations to the nodes and weights of Gauss-Hermite and Gauss-Laguerre quadratures, Stud. Appl. Math. 140(3) (2018) 298–332. https://doi.org/10.1111/sapm.12201 





[12] B. Guo, J. Shen, C. Xu, Generalized Laguerre approximation and its applications to exterior problems, J. Comput. Math. 23 (2005) 113–130. https://www.jstor.org/stable/43693219 





[13] H. Hassani, Z. Avazzadeh, P. Agarwal, S. Mehrabi, M. Ebadi, M.S. Dahaghin, E. Naraghirad, A study on fractional tumor-immune interaction model related to lung cancer via generalized Laguerre polynomials, BMC Med. Res. Methodol. 23(1) (2023) 189. https://doi.org/10.1186/s12874-023-02006-3 





[14] J.A. Jo, Q. Fang, T. Papaioannou, J.D. Baker, A.H. Dorafshar, T. Reil, J.H. Qiao, M.C. Fishbein, J.A. Freischlag, L. Marcu, Laguerre-based method for analysis of time-resolved fluorescence data: application to in-vivo characterization and diagnosis of atherosclerotic lesions, J. Biomed. Opt. 11(2) (2006) 021004. https://doi.org/10.1117/1.2186045 





[15] A.F. Mastryukov, B.G. Mikhailenko, Numerical solution of Maxwell’s equations for anisotropic media using the Laguerre transform, Russ. Geol. Geophys. 49(8) (2008) 621–627. https://doi.org/10.1016/j.rgg.2007.12.011 





[16] B.G. Mikhailenko, Spectral Laguerre method for the approximate solution of time dependent problems, Appl. Math. Lett. 12(4) (1999) 105–110. https://doi.org/10.1016/S0893-9659(99)00043-9 





[17] J.C. Prajapati, N.K. Ajudia, P. Agarwal, Some results due to Konhauser polynomials of first kind and Laguerre polynomials, Appl. Math. Comput. 247 (2014) 639–650. https://doi.org/10.1016/j.amc.2014.09.020 





[18] J.W. Schumer, J.P. Holloway, Vlasov simulations using velocityscaled Hermite representations, J. Comput. Phys. 144(2) (1998) 626–661. https://doi.org/10.1006/jcph.1998.5925 





[19] J. Shen, L.-L. Wang, Some recent advances on spectral methods for unbounded domains, Commun. Comput. Phys. 5(2-4) (2009) 195–241. 





[20] J. Shen, T. Tang, L.-L. Wang, Spectral methods: algorithms, analysis and applications, Springer Science & Business Media, 2011. 





[21] G. Szego, Orthogonal polynomials, American Mathematical Society ˝ , 1939. 





[22] A.V. Terekhov, A fast parallel algorithm for solving block-tridiagonal systems of linear equations including the domain decomposition method, Parallel Comput. 39(6- 7) (2013) 245–258. https://doi.org/10.1016/j.parco.2013.03.003 





[23] H. Wang, Convergence analysis of Laguerre approximations for analytic functions, Math. Comput. 93(350) (2024) 2861–2884. https://doi.org/10.48550/arXiv.2304.05744 





[24] H. Wang, L. Zhang, Convergence analysis of Hermite approximations for analytic functions, Math. Comput.. https://doi.org/10.48550/arXiv.2312.07940 





[25] S. Xiang, Numerical analysis of a fast integration method for highly oscillatory functions, BIT Numer. Math. 47 (2007) 469–482. https://doi.org/10.1007/s10543-007- 0127-y 





[26] S. Xiang, Asymptotics on Laguerre or Hermite polynomial expansions and their applications in Gauss quadrature, J. Math. Anal. Appl. 393(2) (2012) 434–444. https://doi.org/10.1016/j.jmaa.2012.03.056 





[27] S. Xiang, Convergence rates on spectral orthogonal projection approximation for functions of algebraic and logarithmatic regularities, SIAM J. Numer. Anal. 59(3) (2021) 1374–1398. https://doi.org/10.1137/20M134407X 





[28] Y. Zhang, S. Xiang, D. Kong, On optimal convergence rates of Laguerre polynomial expansions for piecewise functions, J. Comput. Appl. Math. 425 (2023) 115053. https://doi.org/10.1016/j.cam.2022.115053 





[29] Y. Zhang, S. Xiang, D. Kong, Optimal pointwise error estimates for piecewise functions expanded with Laguerre polynomials, J. Comput. Appl. Math. 443 (2024) 115749. https://doi.org/10.1016/j.cam.2023.115749 

