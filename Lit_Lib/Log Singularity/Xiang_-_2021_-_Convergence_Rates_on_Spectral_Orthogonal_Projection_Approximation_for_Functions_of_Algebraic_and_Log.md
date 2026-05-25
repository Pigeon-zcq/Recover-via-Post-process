# CONVERGENCE RATES ON SPECTRAL ORTHOGONAL PROJECTION APPROXIMATION FOR FUNCTIONS OF ALGEBRAIC AND LOGARITHMATIC REGULARITIES∗

SHUHUANG XIANG† 

Abstract. Based on the Hilb type formula between Jacobi polynomials and Bessel functions, optimal decay rates on the Jacobi expansion coefficients are derived by applying van der Corput type lemmas for functions of algebraic and logarithmatic singularities, which leads to the optimal convergence rates on the Jacobi, Gegenbauer, and Chebyshev orthogonal projections. It is interesting to see that for boundary singularities, one may get faster convergence rate on the Jacobi or Gegenbauer projection as $( \alpha , \beta )$ and λ increases. The larger values of parameters, the higher convergence rates can be achieved. In particular, the truncated error of Legendre projection has one half order higher than that of Chebyshev projection. Moreover, if min $\left\{ \alpha , \beta \right\} > 0$ and $\lambda > { \frac { 1 } { 2 } } ,$ , the Jacobi and Gegenbauer orthogonal projections have higher convergence orders compared with Legendre. While for interior singularity, the convergence order is independent of $( \alpha , \beta )$ and λ. 

Key words. asymptotic, coefficient, convergence rate, Gegenbauer polynomial, Jacobi polynomial, orthogonal expansion, truncated spectral expansion 

AMS subject classifications. 41A10, 41A25, 41A50, 65N35, 65M70 

DOI. 10.1137/20M134407X 

1. Introduction. The p and hp versions of the finite element methods, or spectral and spectral-element methods have attracted large interest both in theory and computational practice. To deal with corner singularities $( 1 \pm x ) ^ { \gamma } \ln ^ { \mu } ( 1 \pm x )$ (µ a nonnegative integer), non-uniformly Jacobi-weighted Sobolev spaces, $H ^ { m , \beta } ( \Omega )$ with integer $m \geq 0 { \mathrm { ~ a n d ~ } } \beta > - 1$ , are introduced to instead of the standard weighted Sobolev space $H _ { w } ^ { m } ( \Omega )$ with $w ( x )$ a weight on Ω ([2, 3, 4, 8, 9, 10, 12, 14]), which is applied to estimate the orthogonal projection 

$$
\left\| f - \mathcal {P} _ {N} ^ {f} \right\| _ {W} \leq \rho (N) \| f \| _ {H ^ {m, \beta} (\Omega)},
$$

where $\mathcal { P } _ { N } ^ { f } ( x )$ is the truncated polynomial of the Jacobi (Gegenbauer or Chebyshev) expansion of $f ( x )$ , W is a related Sobolev or Besov space, $H ^ { m , \beta } ( \Omega )$ defined as a closure of $C ^ { \infty }$ -functions endowed with the weighted norm. The convergence rate $\rho ( N )$ depends on the regularity exponentials of $H ^ { m , \beta } ( \Omega ) \ [ 1 2 ]$ . 

It is worth noticing that even for $H ^ { m , \beta } ( \Omega )$ with $\Omega = ( - 1 , 1 )$ and weighted norm 

$$
\| u \| _ {H ^ {m, \beta} (\Omega)} = \left\{\sum_ {j = 0} ^ {m} \int_ {- 1} ^ {1} (1 - x ^ {2}) ^ {\beta + j} \left[ u ^ {(j)} (x) \right] ^ {2} d x \right\} ^ {\frac {1}{2}}, \tag {1.1}
$$

it could not lead to the optimal order for (1 + x)γ-type singular functions with noninteger $\gamma > 0$ (see [6, p. 474] and [12]). Indeed, for $f ( x ) = ( 1 + x ) ^ { \gamma } \in H ^ { m , - \frac { 1 } { 2 } } ( \Omega )$ with $\begin{array} { r } { m < 2 \gamma + \frac { 1 } { 2 } } \end{array}$ , the Chebyshev approximation $\Vert f - \mathcal { P } _ { N } ^ { f } \Vert _ { L _ { w } ^ { 2 } ( \Omega ) }$ with $w ( x ) = ( 1 - x ^ { 2 } ) ^ { - \frac { 1 } { 2 } }$ loss an order of the fractional part of $2 \gamma + { \textstyle \frac { 1 } { 2 } }$ , or one order when $\begin{array} { r } { 2 \gamma = k + \frac { 1 } { 2 } } \end{array}$ with nonnegative integer k. For more details, see Liu, Wang, and Li [12]. 

To overcome the above deficiency, Liu, Wang, and Li [12] introduced a new framework of fractional Sobolev-type spaces: generalized Gegenbauer functions of fractional degree (GGF-Fs). Under this framework, the optimal decay rate of Chebyshev expansion coefficients for a large class of functions with interior and endpoint singularities are presented. In addition, Hilb type estimates of GGF-Fs are derived in [11]. However, for interior and logarithmic singularities, it yields 

$$
| x - \theta | ^ {s} \ln | x - \theta | \in W _ {\theta} ^ {s + 1 - \varepsilon} (\Omega) \quad \forall \varepsilon \in (0, 1)
$$

for $\theta \in \Omega = ( - 1 , 1 )$ , and the Chebyshev expansion coefficients $c _ { n }$ satisfy 

$$
| c _ {n} | = O (n ^ {- (s + 1 - \varepsilon)})
$$

(see [12]). But the estimate is still suboptimal. How to modify the fractional space to best characterize this type singularity appears nontrivial and is still open [12, Remark 4.5]. 

For common Gegenbauer expansion, Wang [18] proposed an alternative derivation of the contour integral representation. With this representation, optimal estimates for the Gegenbauer expansion coefficients are derived for analytic functions or functions with endpoint algebraic singularities. However, a precise result on the asymptotic behavior of interior singularities or endpoint logarithmic singularities is still open [18]. 

The main results obtained in this paper are given in Theorems 3.1–3.2, Corollaries 3.3–3.4, Theorem 3.5, Theorems 4.1 and 4.3, and Corollary 4.4. These results show that for the endpoint singularity $( 1 \pm x ) ^ { \gamma } \ln ^ { \mu } ( 1 \pm x ) g ( x )$ with $g \in C ^ { \infty } [ - 1 , 1 ]$ , the Jacobi or Gegenbauer orthogonal projection can achieve faster convergence rate as $( \alpha , \beta )$ and λ increases. The larger values of $( \alpha , \beta )$ and λ, the higher convergence rates can be obtained. In particular, if min $\{ \alpha , \beta \} > - \frac { 1 } { 2 }$ and $\lambda > 0$ , the Jacobi and Gegenbauer orthogonal projections have higher convergence orders than Chebyshev. While for interior singularity $| x - z _ { 0 } | ^ { s } \ln ^ { \mu } | x - z _ { 0 } | g ( x ) ~ ( z _ { 0 } \in ( - 1 , 1 ) )$ , the convergence order is independent of $( \alpha , \beta )$ and λ. Specifically, for the Jacobi expansion, it establishes that (1.2) 

$$
\| f - \mathcal {P} _ {N} ^ {f} \| _ {H ^ {m, \alpha , \beta} (\Omega)} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {m - \alpha - 2 \gamma - 1} \ln^ {\mu} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {m - s - 1 / 2} \ln^ {\mu} (N)), & f (x) = | x - z _ {0} | ^ {s} \ln^ {\mu} | x - z _ {0} | g (x), \\ \mathcal {O} (N ^ {m - \beta - 2 \delta - 1} \ln^ {\mu} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \end{array} \right.
$$

where $g \in C ^ { \infty } [ - 1 , 1 ]$ , µ is a nonnegative integer, $z _ { 0 } \in ( - 1 , 1 )$ , min $\{ \alpha + \gamma , \beta + \delta , \alpha +$ $2 \gamma , \beta + 2 \delta \} > m - 1$ for the boundary singularities, $\begin{array} { r } { s > m - \frac { 1 } { 2 } } \end{array}$ and min $\{ \alpha , \beta \} \geq - \frac { 1 } { 2 }$ for the interior singularity, and 

$$
\| u \| _ {H ^ {m, \alpha , \beta} (\Omega)} = \left\{\sum_ {j = 0} ^ {m} \int_ {- 1} ^ {1} (1 - x) ^ {\alpha + j} (1 + x) ^ {\beta + j} \left[ u ^ {(j)} (x) \right] ^ {2} d x \right\} ^ {\frac {1}{2}}. \tag {1.3}
$$

Furthermore, if $\mathit { \Psi } \gamma , \delta$ are integers and $\mu$ is a positive integer, then 

(1.4) 

$$
\| f - \mathcal {P} _ {N} ^ {f} \| _ {H ^ {m, \alpha , \beta} (\Omega)} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {m - \alpha - 2 \gamma - 1} \ln^ {\mu - 1} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {m - \beta - 2 \delta - 1} \ln^ {\mu - 1} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x). \end{array} \right.
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/e71f5c115d43ee0510c0d2c952d657e2c50a163b221df5c7dc28e70a3f3ef1f6.jpg)



Fig. 1.1. $P _ { 1 0 0 } ^ { ( - 0 . 5 , - 0 . 5 ) } ( x ) ~ ( l e f t ) , ~ P _ { 1 0 0 } ^ { ( 0 , 0 ) } ( x ) ~ ( m i d d l e ) , ~ a n d ~ P _ { 1 0 0 } ^ { ( 1 , 1 ) } ( x ) ~ ( r i g h t )$


In addition, the above optimal estimates on the decay rates of the Jacobi expansion can easily lead to the optimal estimates on Gegenbauer expansion for $\begin{array} { r } { \alpha = \beta = \lambda - \frac { 1 } { 2 } } \end{array}$ and the Chebyshev expansion for $\alpha = \beta = - \frac { 1 } { 2 }$ 

These results and the optimal convergence rates can be deduced from the optimal estimates on the decayed rates of the Jacobi expansion coefficients. To avoid the deficiency of the above frameworks and get the optimal asymptotic orders on the expansion coefficients, we will apply van der Corput type lemmas on highly oscillatory Bessel transforms with a large frequency. 

Note that Jacobi polynomial $P _ { n } ^ { ( \alpha , \beta ) } ( x )$ is an oscillatory function around $x = 0$ particularly, a highly oscillatory function as n $\gg 1$ (see Figure 1.1). The Hilb type formula given in Darboux [7] and Szeg¨o [15, Theorem 8.21.12] establishes the relation between the Jacobi polynomial and a highly oscillatory Bessel function with a larger frequency when the degree of the Jacobi polynomial increases. 

The paper is organized as follows. In section 2, we first present the relationship among the Jacobi, Gegenbauer, and Chebyshev expansion coefficients and derive van der Corput type lemmas for Bessel transforms. Based on the Hilb type formula, the optimal decay rates on the expansion coefficients for logarithmatic singularities are presented in section $^ { 3 , }$ and the convergence rates on the spectral orthogonal projections are given in section 4. Final remarks are included in section 5. 

2. Preliminaries. Assume $f ( x )$ is a suitably smooth function on $[ - 1 , 1 ]$ . Consider the continuous polynomial expansion 

$$
f (x) = \sum_ {n = 0} ^ {\infty} a _ {n} (\alpha , \beta) P _ {n} ^ {(\alpha , \beta)} (x), \quad \alpha , \beta > - 1 \tag {2.1}
$$

with the expansion coefficients 

$$
a _ {n} (\alpha , \beta) = \frac {1}{\sigma_ {n} ^ {\alpha , \beta}} \int_ {- 1} ^ {1} (1 - x) ^ {\alpha} (1 + x) ^ {\beta} f (x) P _ {n} ^ {(\alpha , \beta)} (x)   \mathrm{d} x, \tag {2.2}
$$

where $P _ { n } ^ { ( \alpha , \beta ) } ( x )$ is the Jacobi polynomial of degree n and 

$$
\sigma_ {n} ^ {\alpha , \beta} = 2 ^ {\alpha + \beta + 1} \frac {\Gamma (n + \alpha + 1) \Gamma (n + \beta + 1)}{n ! (2 n + \alpha + \beta + 1) \Gamma (n + \alpha + \beta + 1)} \tag {2.3}
$$

(see [1, p. 774]), which leads to the error bound on the orthogonal projection with respect to the weight function $w ( \alpha , \beta ) = ( 1 - x ) ^ { \alpha } ( 1 + x ) ^ { \beta }$ 

$$
\| f - \mathcal {P} _ {N} ^ {f} \| _ {L _ {w (\alpha , \beta)} ^ {2} [ - 1, 1 ]} = \sqrt {\sum_ {n = N + 1} ^ {\infty} a _ {n} ^ {2} (\alpha , \beta) \sigma_ {n} ^ {\alpha , \beta}} \tag {2.4}
$$

and implies that the convergence is only on the decay of the expansion coefficients. 

In particular, the optimal estimates on the decay rates of the Jacobi expansion can easily lead to the optimal estimates on Gegenbauer expansion, then derives the estimates on Chebyshev expansion. 

The Gegenbauer expansion with respect to the weight function $w ( \lambda ) = ( 1 - x ^ { 2 } ) ^ { \lambda - \frac { 1 } { 2 } }$ is 

$$
f (x) = \sum_ {n = 0} ^ {\infty} a _ {n} (\lambda) C _ {n} ^ {(\lambda)} (x), \qquad \lambda > - \frac {1}{2}, \quad \lambda \neq 0 \tag {2.5}
$$

with the expansion coefficients 

$$
a _ {n} (\lambda) = \frac {1}{\hbar_ {n}} \int_ {- 1} ^ {1} (1 - x ^ {2}) ^ {\lambda - \frac {1}{2}} C _ {n} ^ {(\lambda)} (x) f (x) \mathrm{d} x, \quad \hbar_ {n} = \frac {2 ^ {1 - 2 \lambda} \pi}{\Gamma^ {2} (\lambda)} \frac {\Gamma (n + 2 \lambda)}{n ! (n + \lambda)}
$$

(see $\left[ 1 0 , \mathrm { p } . \ 7 9 \right] )$ , where $C _ { n } ^ { ( \lambda ) } ( x )$ is the Gegenbauer polynomial of degree n 

$$
C _ {n} ^ {(\lambda)} (x) = \frac {\Gamma \left(\lambda + \frac {1}{2}\right) \Gamma (n + 2 \lambda)}{\Gamma (2 \lambda) \Gamma \left(n + \lambda + \frac {1}{2}\right)} P _ {n} ^ {(\lambda - \frac {1}{2}, \lambda - \frac {1}{2})} (x), \tag {2.6}
$$

which is related to Legendre and Chebyshev polynomials as follows [10, p. 76] 

$$
P _ {n} (x) = C _ {n} ^ {(\frac {1}{2})} (x), \quad T _ {n} (x) = n \lim _ {\lambda \to 0} \Gamma (2 \lambda) C _ {n} ^ {(\lambda)} (x). \tag {2.7}
$$

From (2.1) and (2.6) it follows 

$$
\begin{array}{l} f (x) = \sum_ {n = 0} ^ {\infty} a _ {n} \left(\lambda - \frac {1}{2}, \lambda - \frac {1}{2}\right) P _ {n} ^ {\left(\lambda - \frac {1}{2}, \lambda - \frac {1}{2}\right)} (x) \\ = \sum_ {n = 0} ^ {\infty} a _ {n} \left(\lambda - \frac {1}{2}, \lambda - \frac {1}{2}\right) \frac {\Gamma (2 \lambda) \Gamma \left(n + \lambda + \frac {1}{2}\right)}{\Gamma \left(\lambda + \frac {1}{2}\right) \Gamma (n + 2 \lambda)} C _ {n} ^ {(\lambda)} (x) \\ = \sum_ {n = 0} ^ {\infty} a _ {n} (\lambda) C _ {n} ^ {(\lambda)} (x), \\ \end{array}
$$

which derives 

$$
a _ {n} (\lambda) = a _ {n} \left(\lambda - \frac {1}{2}, \lambda - \frac {1}{2}\right) \frac {\Gamma (2 \lambda) \Gamma \left(n + \lambda + \frac {1}{2}\right)}{\Gamma \left(\lambda + \frac {1}{2}\right) \Gamma (n + 2 \lambda)}. \tag {2.8}
$$

Moreover, note that 

$$
f (x) = \sum_ {n = 0} ^ {\infty} a _ {n} (\lambda) C _ {n} ^ {(\lambda)} (x) = \sum_ {n = 0} ^ {\infty} \frac {a _ {n} (\lambda)}{n \Gamma (2 \lambda)} n \Gamma (2 \lambda) C _ {n} ^ {(\lambda)} (x),
$$

which, together with 

$$
f (x) = \sum_ {n = 0} ^ {\infty} c _ {n} T _ {n} (x) \tag {2.9}
$$

and (2.7)–(2.8), derives $c _ { 0 } = a _ { 0 } ( - \textstyle \frac { 1 } { 2 } , - \textstyle \frac { 1 } { 2 } )$ and for $n \geq 1$ 

$$
\begin{array}{l} c _ {n} = \lim _ {\lambda \to 0} \frac {a _ {n} (\lambda)}{n \Gamma (2 \lambda)} = \lim _ {\lambda \to 0} \frac {a _ {n} \left(\lambda - \frac {1}{2} , \lambda - \frac {1}{2}\right) \frac {\Gamma (2 \lambda) \Gamma \left(n + \lambda + \frac {1}{2}\right)}{\Gamma (\lambda + \frac {1}{2}) \Gamma (n + 2 \lambda)}}{n \Gamma (2 \lambda)} \tag {2.10} \\ = \frac {\Gamma \left(n + \frac {1}{2}\right)}{n \Gamma \left(\frac {1}{2}\right) \Gamma (n)} a _ {n} \left(- \frac {1}{2}, - \frac {1}{2}\right). \\ \end{array}
$$

The following asymptotic Hilb type formula for Jacobi polynomials, related to a highly oscillatory Bessel function with a large frequency, is introduced in Darboux [7] and Szeg¨o [15, Theorem 8.21.12]. 

Lemma 2.1 ([7, 15]). Let $\alpha , \beta > - 1 ,$ then as $n \to \infty$ 

$$
\begin{array}{l} \theta^ {- \frac {1}{2}} \sin^ {\alpha + \frac {1}{2}} \left(\frac {\theta}{2}\right) \cos^ {\beta + \frac {1}{2}} \left(\frac {\theta}{2}\right) P _ {n} ^ {(\alpha , \beta)} (\cos \theta) \\ = \frac {\Gamma (n + \alpha + 1)}{\sqrt {2} n ! \tilde {N} ^ {\alpha}} J _ {\alpha} (\tilde {N} \theta) + \left\{ \begin{array}{l l} \theta^ {\frac {1}{2}} \mathcal {O} \left(\tilde {N} ^ {- \frac {3}{2}}\right), & c n ^ {- 1} \leq \theta \leq \pi - \epsilon , \\ \theta^ {\alpha + 2} \mathcal {O} \left(\tilde {N} ^ {\alpha}\right), & 0 <   \theta \leq c n ^ {- 1}, \end{array} \right. \tag {2.11} \\ \end{array}
$$

where $\tilde { N } = n + ( \alpha + \beta + 1 ) / 2 , c ,$ , and  are fixed positive numbers, and $J _ { \alpha } ( z )$ is the first kind of Bessel function of order α. The constants in the O-terms depend on $\alpha ,$ $\beta , c ,$ and . 

Lemma 2.2. Suppose $\alpha + \nu > 0 , b > 0$ , and $\mu$ is a nonnegative integer; then it is satisfied $f o r \ x \in [ 0 , b ]$ that 

$$
\left| \ln^ {\mu} (\omega x) (\omega x) ^ {\alpha} J _ {\nu} (\omega x) \right| = \left\{ \begin{array}{l l} \mathcal {O} \left((\ln^ {\mu} (\omega)) \omega^ {\alpha - \frac {1}{2}}\right), & \alpha \geq \frac {1}{2}, \\ \mathcal {O} (1), & \alpha <   \frac {1}{2}, \end{array} \right. \quad \omega \gg 1. \tag {2.12}
$$

Proof. Define $F _ { 1 } ( z ) = \ln ^ { \mu } ( z ) z ^ { \alpha } J _ { \nu } ( z )$ for $z \in [ 0 , + \infty )$ , where $F _ { 1 } ( 0 )$ is defined by its limit as z tends to 0. From the definition of $J _ { \nu }$ (Abramowitz and Stegun [1, Eq. 9.1.10]) 

$$
J _ {\nu} (z) = \left(\frac {z}{2}\right) ^ {\nu} \sum_ {n = 0} ^ {\infty} \frac {\left(- \frac {1}{4} z ^ {2}\right) ^ {n}}{n ! \Gamma (\nu + n + 1)}, \tag {2.13}
$$

we see that $F _ { 1 } ( 0 ) = 0 , F _ { 1 }$ is continuous for $z \in [ 0 , + \infty )$ . In addition, from [1, Eq. (9.1.30), Eq. (9.2.1)] and [20, p. 199] 

$$
J _ {\nu} (z) = \sqrt {\frac {2}{\pi z}} \cos \left(z - \frac {1}{2} \nu \pi - \frac {\pi}{4}\right) + \mathcal {O} \left(z ^ {- \frac {3}{2}}\right), \quad z \to + \infty ,
$$

there exists a $z _ { 0 } \geq 1$ such that for $z \geq z _ { 0 } , | J _ { \nu } ( z ) | \leq C _ { 1 } z ^ { - \frac { 1 } { 2 } }$ for some positive constant $C _ { 1 }$ independent of $z ,$ which implies for $\omega x \geq z _ { 0 }$ 

$$
\left| F _ {1} (\omega x) \right| \leq C _ {1} \ln^ {\mu} (\omega x) (\omega x) ^ {\alpha - \frac {1}{2}}, \quad \alpha \geq \frac {1}{2}; \quad \left| F _ {1} (\omega x) \right| = \mathcal {O} (1), \quad \alpha <   \frac {1}{2};
$$

then (2.12) is satisfied for ωx $\geq z _ { 0 }$ . Notice that $F _ { 1 } ( z )$ is uniformly bounded on $[ 0 , z _ { 0 } ]$ , which yields for $\omega x < z _ { 0 }$ 

$$
| F _ {1} (\omega x) | = \mathcal {O} (1).
$$

These together complete the proof. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/28fec4e2a9c6a1f61b5c17e5d1f7a9695cb3baacc7465f9378e62a7773685e9d.jpg)


Lemma 2.3 (van der Corput lemma for Bessel transform I). Suppose $\alpha + \nu > - 1$ , $b > 0 , \mu$ is a nonnegative integer, $\psi \in C [ 0 , b ]$ and $\psi ^ { \prime } \in L ^ { 1 } [ 0 , b ]$ ; then it is satisfied for $\omega \gg 1$ , 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu} (\omega x)   \mathrm{d} x = \left\{ \begin{array}{l l} \mathcal {O} \left(\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}\right), & \alpha \leq \frac {1}{2}, \\ \mathcal {O} \left(\omega^ {- \frac {3}{2}}\right), & \alpha > \frac {1}{2}, \end{array} \right. \tag {2.14}
$$

and 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} \psi (x) J _ {\nu} (\omega x)   \mathrm{d} x = \left\{ \begin{array}{l l} \mathcal {O} \left(\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}\right), & \alpha \leq \frac {1}{2}, \\ \mathcal {O} \left(\omega^ {- \frac {3}{2}}\right), & \alpha > \frac {1}{2}, \end{array} \right. \tag {2.15}
$$

uniformly $f o r t \in [ 0 , b ]$ . 

Proof. In the case $\alpha < \textstyle { \frac { 1 } { 2 } }$ : Let $u = \omega x$ . It follows 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu} (\omega x) \mathrm{d} x = \frac {1}{\omega^ {\alpha + 1}} \int_ {0} ^ {\omega t} \sum_ {j = 0} ^ {\mu} (- 1) ^ {j} C _ {\mu} ^ {j} \ln^ {\mu - j} (u) \ln^ {j} (\omega) u ^ {\alpha} J _ {\nu} (u) \mathrm{d} u, \tag {2.16}
$$

where $\begin{array} { r } { C _ { \mu } ^ { j } = \frac { \Gamma \left( \mu + 1 \right) } { \Gamma \left( j + 1 \right) \Gamma \left( \mu - j + 1 \right) } } \end{array}$ . It is worthy of noting that 

$$
\int_ {0} ^ {+ \infty} u ^ {\alpha} J _ {\nu} (u) d u = \left\{ \begin{array}{l l} \frac {2 ^ {\alpha} \Gamma \left(\frac {\alpha + \nu + 1}{2}\right)}{\Gamma \left(\frac {\nu - \alpha + 1}{2}\right)} <   + \infty , & \nu - \alpha + 1 \neq 0, \\ x ^ {\nu + 1} J _ {\nu + 1} (x) \Big | _ {0} ^ {+ \infty} = 0, & \nu - \alpha + 1 = 0, \end{array} \right. \Re (\alpha + \nu) > - 1, \Re (\alpha) <   \frac {1}{2},
$$

([1, Eq. (9.1.30), Eq. (11.4.16)]), which implies that $\begin{array} { r l } { \int _ { 0 } ^ { \omega t } u ^ { \alpha } J _ { \nu } ( u ) } \end{array}$ du is uniformly bounded for $t \in [ 0 , + \infty )$ since this integral with a variable upper bound is continuous and convergent as $t \to + \infty$ . 

Then from the assumption $\alpha + \nu > - 1$ , we may choose η with $0 < \eta <$ < min $\bigl \{ \frac { 1 } { 2 } -$ $\alpha , \alpha + \nu + 1 \}$ such that $\alpha \pm \eta < \frac { 1 } { 2 } , \alpha + \nu \pm \eta > - 1$ , which implies for $j = 0 , 1 , \ldots , \mu$ that $\begin{array} { r l } { \int _ { 0 } ^ { 1 } [ u ^ { \eta } \ln ^ { \mu - j } ( u ) ] u ^ { \alpha - \eta } J _ { \nu } ( u ) } \end{array}$ du is convergent, and $\begin{array} { r } { \int _ { 1 } ^ { + \infty } [ u ^ { - \eta } \ln ^ { \mu - j } ( u ) ] u ^ { \alpha + \eta } J _ { \nu } ( u ) } \end{array}$ du is also convergent by Abel criterion for improper integrals [13] due to that $\int _ { 1 } ^ { + \infty } u ^ { \alpha + \eta }$ $\boldsymbol { J } _ { \nu } ( u ) \mathrm { d } u$ is convergent and $[ u ^ { - \eta } \ln ^ { \mu - j } ( u ) ]$ is decreasing and tends to 0 as $u  + \infty$ . Thus, $\begin{array} { r } { \int _ { 0 } ^ { + \infty } \ln ^ { \mu - j } ( u ) u ^ { \alpha } J _ { \nu } ( u ) } \end{array}$ du is convergent and $\begin{array} { r l } { \int _ { 0 } ^ { \omega t } \ln ^ { \mu - j } ( u ) u ^ { \alpha } J _ { \nu } ( u ) } \end{array}$ du is uniformly bounded for $t \in [ 0 , + \infty )$ too. These together with (2.16) lead to the desired result (2.14) in the case $\alpha < \frac { 1 } { 2 }$ . 

In the case $\alpha = { \textstyle { \frac { 1 } { 2 } } } \colon$ Noting that $[ z ^ { \nu + 1 } J _ { \nu + 1 } ( z ) ] ^ { \prime } = z ^ { \nu + 1 } J _ { \nu } ( z )$ , the integral can be represented for $\mu \geq 1$ as 

$$
\begin{array}{l} \int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu} (\omega x) \mathrm{d} x = \omega^ {- 1} \int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha - \nu - 1} \mathrm{d} [ x ^ {\nu + 1} J _ {\nu + 1} (\omega x) ] \\ = \omega^ {- 1} \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu + 1} (\omega x) \Big | _ {0} ^ {t} \\ - \omega^ {- 1} \int_ {0} ^ {t} [ \mu x ^ {\alpha - 1} \ln^ {\mu - 1} (u) + (\alpha - \nu - 1) \ln^ {\mu} (x) x ^ {\alpha - 1} ] J _ {\nu + 1} (\omega x) d x \\ = \omega^ {- 1 - \alpha} [ \ln (\omega t) - \ln (\omega) ] ^ {\mu} (\omega t) ^ {\alpha} J _ {\nu + 1} (\omega t) \\ - \omega^ {- 1} \int_ {0} ^ {t} [ \mu x ^ {\alpha - 1} \ln^ {\mu - 1} (u) + (\alpha - \nu - 1) \ln^ {\mu} (x) x ^ {\alpha - 1} ] J _ {\nu + 1} (\omega x) d x \\ = \mathcal {O} \left(\ln^ {\mu} (\omega) \omega^ {- 1 - \alpha}\right) \\ \end{array}
$$

by Lemma 2.2 and the above proof since $\alpha - 1 < { \frac { 1 } { 2 } }$ . Similarly, the above estimate is also satisfied for $\mu = 0$ . 

In the case $\begin{array} { r } { \alpha > \frac { 1 } { 2 } ; } \end{array}$ : It yields 

$$
\begin{array}{l} \int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu} (\omega x) d x \\ = \int_ {0} ^ {t} [ x ^ {\alpha - \frac {1}{2}} \ln^ {\mu} (x) ] x ^ {\frac {1}{2}} J _ {\nu} (\omega x) d x \\ = \int_ {0} ^ {t} \left[ x ^ {\alpha - \frac {1}{2}} \ln^ {\mu} (x) \right] d \left[ \int_ {0} ^ {x} u ^ {\frac {1}{2}} J _ {\nu} (\omega u) d u \right] \\ = \left[ t ^ {\alpha - \frac {1}{2}} \ln^ {\mu} (t) \right] \int_ {0} ^ {t} x ^ {\frac {1}{2}} J _ {\nu} (\omega x) d x - \int_ {0} ^ {t} \left[ x ^ {\alpha - \frac {1}{2}} \ln^ {\mu} (x) \right] ^ {\prime} \left(\int_ {0} ^ {x} u ^ {\frac {1}{2}} J _ {\nu} (\omega u) d u\right) d x \\ = \mathcal {O} \left(\omega^ {- \frac {3}{2}}\right) + \mathcal {O} \left(\omega^ {- \frac {3}{2}}\right) \int_ {0} ^ {t} [ x ^ {\alpha - \frac {1}{2}} \ln^ {\mu} (x) ] ^ {\prime} d x \\ = \mathcal {O} \left(\omega^ {- \frac {3}{2}}\right). \\ \end{array}
$$

Expression (2.15) follows from (2.14) and 

$$
\begin{array}{l} \int_ {0} ^ {t} \psi (x) \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu} (\omega x) d x \\ = \int_ {0} ^ {t} \psi (x) \left[ \int_ {0} ^ {x} \ln^ {\mu} (u) u ^ {\alpha} J _ {\nu} (\omega u) d u \right] ^ {\prime} d x \\ = \psi (t) \int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu} (\omega x) d x - \int_ {0} ^ {t} \psi^ {\prime} (x) \left[ \int_ {0} ^ {x} \ln^ {\mu} (u) u ^ {\alpha} J _ {\nu} (\omega u) d u \right] d x. \\ \end{array}
$$

Remark 1. From the proof of Lemma 2.3, we see that for $\alpha + \nu > - 1$ , 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu} (\omega x) d x = \mathcal {O} \left(\ln^ {\mu} (\omega) \omega^ {- 1 - \alpha}\right)
$$

is satisfied uniformly for $t \in [ 0 , + \infty )$ in the case $\alpha < \textstyle { \frac { 1 } { 2 } }$ , while 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\frac {1}{2}} J _ {\nu} (\omega x) d x = \mathcal {O} \left(\ln^ {\mu} (\omega) \omega^ {- \frac {3}{2}}\right), \int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} J _ {\nu} (\omega x) d x = \mathcal {O} \left(\omega^ {- \frac {3}{2}}\right) \left(\alpha > \frac {1}{2}\right)
$$

uniformly for $t \in [ 0 , b ]$ . In particular, for $\mu = 0 \ [ 2 4 ]$ 

(2.17) 

$$
\int_ {0} ^ {t} x ^ {\alpha} J _ {\nu} (\omega x) d x = \mathcal {O} \left(\omega^ {- \min \{1 + \alpha , \frac {3}{2} \}}\right) \left\{ \begin{array}{l l} \text { uniformly   for } t \in [ 0, + \infty) \text { for } \alpha <   \frac {1}{2}, \\ \text { uniformly   for } t \in [ 0, b ] \text { for } \alpha \geq \frac {1}{2}. \end{array} \right.
$$

Lemma 2.4 ([21]). For $0 < a < b , \psi \in C [ a , b ]$ and $\psi ^ { \prime } \in L ^ { 1 } [ a , b ]$ , 

$$
\int_ {a} ^ {b} \psi (t) J _ {\nu} (\omega t) d t = \mathcal {O} \left(\omega^ {- \frac {3}{2}}\right).
$$

Lemma 2.5 (van der Corput lemma for Bessel transform II). Suppose $\alpha { + } \nu > { - } 1$ , $\beta > - 1 , b > 0 ;$ , µ is a nonnegative integer, $\psi \in C [ 0 , b ]$ and $\hat { \psi } ^ { \prime } \in L ^ { 1 } [ 0 , b ]$ ; then it is satisfied $f o r \omega \gg 1$ , 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x)   \mathrm{d} x = \mathcal {O} \left(\max \left\{\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\min \left\{\beta + \frac {3}{2} , \frac {3}{2} \right\}}} \right\}\right) \tag {2.18}
$$

uniformly for $t \in [ 0 , b ]$ . In particular, $i f b = 1$ and $\mu \geq 1$ , 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (1 - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x)   \mathrm{d} x = \mathcal {O} \left(\max \left\{\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\frac {3}{2}}} \right\}\right). \tag {2.19}
$$

Proof. For $t \in [ 0 , \frac { b } { 2 } ]$ , from Lemma 2.3, it establishes 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x = \mathcal {O} \left(\max \left\{\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\frac {3}{2}}} \right\}\right).
$$

For $\textstyle t \in ( { \frac { b } { 2 } } , b ]$ , the integral can be written as 

$$
\begin{array}{l} \int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x)   \mathrm{d} x \tag {2.20} \\ = \int_ {0} ^ {\frac {b}{2}} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x + \int_ {\frac {b}{2}} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x \\ = \mathcal {O} \left(\max \left\{\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\frac {3}{2}}} \right\}\right) + \int_ {\frac {b}{2}} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x. \\ \end{array}
$$

For the second term on the right-hand side of the first identity (2.20), setting $F ( z ) =$ $\begin{array} { r } { \int _ { \frac { b } { 2 } } ^ { z } \ln ^ { \mu } ( x ) x ^ { \alpha } J _ { \nu } ( \omega x ) d x } \end{array}$ , from Lemma 2.4 it yields $F ( z ) ~ = ~ \mathcal { O } ( \omega ^ { - \frac { 3 } { 2 } } )$ . Moreover, for $\begin{array} { r } { t \in [ \frac { b } { 2 } , b - \frac { 1 } { \omega } ] } \end{array}$ it follows that 

$$
\begin{array}{l} \left| \int_ {\frac {b}{2}} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x \right| \\ = \left| \int_ {\frac {b}{2}} ^ {t} (b - x) ^ {\beta} F ^ {\prime} (x) \mathrm{d} x \right| \\ \leq \max \left\{\left(\frac {b}{2}\right) ^ {\beta}, | b - t | ^ {\beta} \right\} | F (t) | + | \beta | \| F \| _ {\infty} \int_ {\frac {b}{2}} ^ {b - \frac {1}{\omega}} (b - x) ^ {\beta - 1} \mathrm{d} x \\ = \mathcal {O} \left(\omega^ {- \min \left\{\beta + \frac {3}{2}, \frac {3}{2} \right\}}\right). \\ \end{array}
$$

For $\begin{array} { r } { t \in ( b - \frac 1 { \omega } , b ] } \end{array}$ it follows 

$$
\begin{array}{l} \int_ {\frac {b}{2}} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x \\ = \int_ {\frac {b}{2}} ^ {b - \frac {1}{\omega}} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x + \int_ {b - \frac {1}{\omega}} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x \\ = \mathcal {O} \left(\omega^ {- \min \left\{\beta + \frac {3}{2}, \frac {3}{2} \right\}}\right) + \int_ {b - \frac {1}{\omega}} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x. \\ \end{array}
$$

From $J _ { \nu } ( \omega z ) = \mathcal { O } ( \omega ^ { - \frac { 1 } { 2 } } )$ for $z \in [ \frac { b } { 2 } , b ] \ ( [ 1 ] )$ , it implies 

$$
\int_ {b - \frac {1}{\omega}} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) \mathrm{d} x = \mathcal {O} \left(\omega^ {- \frac {1}{2}}\right) \int_ {b - \frac {1}{\omega}} ^ {t} (b - x) ^ {\beta} \mathrm{d} x = \mathcal {O} \left(\omega^ {- \beta - \frac {3}{2}}\right).
$$

These together lead to 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x)   \mathrm{d} x = \mathcal {O} \left(\max \left\{\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\min \{\beta + \frac {3}{2} , \frac {3}{2} \}}} \right\}\right). \tag {2.21}
$$

Similarly, (2.18) follows from (2.21) and 

$$
\begin{array}{l} \int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x) \mathrm{d} x \\ = \int_ {0} ^ {t} \psi (x) \left[ \int_ {0} ^ {x} \ln^ {\mu} (u) u ^ {\alpha} (b - u) ^ {\beta} J _ {\nu} (\omega u) d u \right] ^ {\prime} d x \\ = \psi (t) \int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} J _ {\nu} (\omega x) d x - \int_ {0} ^ {t} \psi^ {\prime} (x) \left[ \int_ {0} ^ {x} \ln^ {\mu} (u) u ^ {\alpha} (b - u) ^ {\beta} J _ {\nu} (\omega u) d u \right] d x. \\ \end{array}
$$

Particularly, for b = 1 and $\mu \geq 1$ , rewriting 

$$
\int_ {0} ^ {t} \ln^ {\mu} (x) x ^ {\alpha} (1 - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x) \mathrm{d} x = \int_ {0} ^ {t} [ \ln (x) (1 - x) ^ {- 1} ] ^ {\mu} x ^ {\alpha} (1 - x) ^ {\beta + \mu} \psi (x) J _ {\nu} (\omega x) \mathrm{d} x
$$

yields (2.19) due to that $\begin{array} { r } { \beta + \mu + \frac { 3 } { 2 } > \frac { 3 } { 2 } . } \end{array}$ 

Lemma 2.6 (Generalized van der Corput lemma [23]). Suppose ω $\gg 1$ and $\phi ( x ) \in C ^ { \infty } [ 0 , b ]$ , then it is satisfied that 

$$
\int_ {0} ^ {b} x ^ {\alpha} (b - x) ^ {\delta - 1} \phi (x) J _ {\nu} (\omega x)   \mathrm{d} x = O \left(\omega^ {- \min \left\{\alpha + 1, \delta + \frac {1}{2} \right\}}\right) \tag {2.22}
$$

for $\alpha > - 1 , \alpha + \nu > - 1 , \delta > 0 .$ 

Lemma 2.7 (van der Corput lemma for Bessel transform III). Suppose $\alpha + \nu >$ $- 1 , \ : \beta > - 1 , \ : b > 0 , \ : \mu$ is a nonnegative integer, and $\psi \in C ^ { \infty } [ 0 , b ]$ ; then it is satisfied $f o r \omega \gg 1$ , 

$$
\int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x)   \mathrm{d} x = \mathcal {O} \left(\max \left\{\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\beta + \frac {3}{2}}} \right\}\right). \tag {2.23}
$$

In particular, if b = 1 and $\mu \geq 1$ , 

$$
\int_ {0} ^ {1} \ln^ {\mu} (x) x ^ {\alpha} (1 - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x)   \mathrm{d} x = \mathcal {O} \left(\max \left\{\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\beta + \mu + \frac {3}{2}}} \right\}\right). \tag {2.24}
$$

Proof. If min $\textstyle \left\{ \alpha + 1 , \beta + { \frac { 3 } { 2 } } \right\} \leq { \frac { 3 } { 2 } }$ , setting $\begin{array} { r } { F ( x ) = \int _ { 0 } ^ { x } \mathrm { l n } ^ { \mu } ( t ) t ^ { \alpha } ( b - t ) ^ { \beta } J _ { \nu } ( \omega t ) } \end{array}$ dt, it follows 

$$
\int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x) \mathrm{d} x = \int_ {0} ^ {b} \psi (x) \mathrm{d} F (x) = \psi (b) F (b) - \int_ {0} ^ {b} \psi^ {\prime} (x) F (x) \mathrm{d} x,
$$

which implies the desired result by Lemma 2.5. 

$\mathrm { I f ~ } \textstyle { \frac { 3 } { 2 } } < \operatorname* { m i n } \left\{ \alpha + 1 , \beta + \frac { 3 } { 2 } \right\} \le \frac { 5 } { 2 }$ , by integrating by parts, it follows for $\mu \geq 1$ that 

$$
\begin{array}{l} \int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\alpha} (b - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x) \mathrm{d} x \\ = \omega^ {- 1} \int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\alpha - \nu - 1} (b - x) ^ {\beta} \psi (x) \mathrm{d} [ x ^ {\nu + 1} J _ {\nu + 1} (\omega x) ] \\ = - \omega^ {- 1} \int_ {0} ^ {b} \ln^ {\mu} (x) x ^ {\alpha - 1} (b - x) ^ {\beta - 1} \phi (x) J _ {\nu + 1} (\omega x) d x \\ - \mu \omega^ {- 1} \int_ {0} ^ {b} x ^ {\alpha - 1} (b - x) ^ {\beta - 1} \ln^ {\mu - 1} (x) (b - x) \psi (x) J _ {\nu + 1} (\omega x) d x \\ = \mathcal {O} \left(\max \left\{\frac {\ln^ {\mu} (\omega)}{\omega^ {\alpha + 1}}, \frac {1}{\omega^ {\beta + \frac {3}{2}}} \right\}\right), \\ \end{array}
$$

where $\phi ( x ) = ( \alpha - \nu - 1 ) ( b - x ) \psi ( x ) - \beta x \psi ( x ) + x ( b - x ) \psi ^ { \prime } ( x )$ . In the case $\mu = 0$ , it directly follows from Lemma 2.6. 

By induction on k for $k + \textstyle { \frac { 3 } { 2 } } <$ < min $\textstyle \left\{ \alpha + 1 , \beta + { \frac { 3 } { 2 } } \right\} \leq k + { \frac { 5 } { 2 } }$ with $k = 0 , 1 , . . . ,$ by a similar proof to the above, it is easy to derive the desired result (2.23). 

The special case $b = 1$ and $\mu \geq 1$ follows from (2.23) together with 

$$
\int_ {0} ^ {1} \ln^ {\mu} (x) x ^ {\alpha} (1 - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x) \mathrm{d} x = \int_ {0} ^ {1} [ \ln (x) (1 - x) ^ {- 1} ] ^ {\mu} x ^ {\alpha} (1 - x) ^ {\beta + \mu} \psi (x) J _ {\nu} (\omega x) \mathrm{d} x.
$$

□ 

By a similar proof, the following are satisfied. 

Lemma 2.8 (van der Corput lemma for Bessel transform IV). Suppose $\alpha + \nu >$ $- 1 , \ : \beta > - 1 , \ : b > 0 , \ : \mu$ is a nonnegative integer, $\psi \in C [ 0 , b ]$ and $\psi ^ { \prime } \in L ^ { 1 } [ 0 , b ] ,$ then it is satisfied $f o r \omega \gg 1$ , 

$$
\int_ {0} ^ {t} \ln^ {\mu} (b - x) x ^ {\alpha} (b - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x)   \mathrm{d} x = \mathcal {O} \left(\max \left\{\frac {1}{\omega^ {\alpha + 1}}, \frac {\ln^ {\mu} (\omega)}{\omega^ {\min \{\beta + \frac {3}{2} , \frac {3}{2} \}}} \right\}\right) \tag {2.25}
$$

uniformly for $t \in [ 0 , b ]$ . In particular, $i f b = 1$ and $\mu \geq 1$ , α in (2.25) can be replaced by $\alpha + \mu .$ . 

Lemma 2.9 (van der Corput lemma for Bessel transform V). Suppose $\alpha { + } \nu > { - } 1$ , $\beta > - 1 , b > 0 ,$ µ is a nonnegative integer, and $\psi \in C ^ { \infty } [ 0 , b ] .$ ; then it is satisfied for $\omega \gg 1$ , 

$$
\int_ {0} ^ {b} \ln^ {\mu} (b - x) x ^ {\alpha} (b - x) ^ {\beta} \psi (x) J _ {\nu} (\omega x)   \mathrm{d} x = \mathcal {O} \left(\max \left\{\frac {1}{\omega^ {\alpha + 1}}, \frac {\ln^ {\mu} (\omega)}{\omega^ {\beta + \frac {3}{2}}} \right\}\right). \tag {2.26}
$$

In particular, if $b = 1$ and $\mu \geq 1$ , α in (2.26) can be replaced by $\alpha + \mu$ . 

The asymptotics on Lemma 2.7 and Lemma 2.9 are illustrated with different values of $( \alpha , \beta )$ by $\begin{array} { r } { \int _ { 0 } ^ { b } \ln ^ { \mu } ( x ) x ^ { \alpha } ( b - x ) ^ { \beta } \psi ( x ) J _ { \nu } ( \omega x ) d x } \end{array}$ with $b = { \textstyle { \frac { 1 } { 2 } } } , \mu = 1$ , and $\psi ( x ) = $ cos(x) for $\omega = 1 : 1$ , 000 in Figure 2.1, which shows the estimates are tight. 

# 3. On functions of limited regularity at endpoints or interior points.

# 3.1. Functions with boundary regularities. Now we consider

$$
f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \tag {3.1}
$$

where µ is a positive integer and $g \in C ^ { \infty } [ - 1 , 1 ]$ . 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/13d9398d855656d56b9450b1815d41f5065350c310464eb50d281f47cc9e07df.jpg)



Fig. 2.1. The asymptotics on Lemma 2.7 and Lemma 2.9 are illustrated with different values of $( \alpha , \beta )$ .


Theorem 3.1. Suppose $f ( x )$ is defined by (3.1); then the Jacobi coefficients (2.2) satisfy for $\alpha + \gamma > - 1$ that 

$$
| a _ {n} (\alpha , \beta) | = \left\{ \begin{array}{l l} \mathcal {O} \left(n ^ {- \alpha - 2 \gamma - 1} \ln^ {\mu} (n)\right), & \gamma \text {   is   not   an   integer }, \\ \mathcal {O} \left(n ^ {- \alpha - 2 \gamma - 1} \ln^ {\mu - 1} (n)\right), & \gamma \text {   is   an   integer }. \end{array} \right. \tag {3.2}
$$

Proof. Let $k _ { 0 }$ be a positive integer such that $\alpha + 2 \gamma - k _ { 0 } + 2 \leq \frac { 3 } { 2 } \leq \beta + k _ { 0 } + 2$ . Then it is easy to verify that $f ^ { ( k ) } ( x )$ can be represented for $\gamma > 0$ not an integer as for $k = 1 , 2 , \ldots , k _ { 0 }$ 

$$
\begin{array}{l} f ^ {(k)} (x) = (1 - x) ^ {\gamma - k} \left[ \ln^ {\mu} (1 - x) h _ {\mu , k} (x) + \ln^ {\mu - 1} (1 - x) h _ {\mu - 1, k} (x) \right. \\ \left. + \dots + \ln (1 - x) h _ {1, k} (x) + h _ {0, k} (x) \right] \psi (x) \\ \end{array}
$$

with $h _ { j , k } \in C ^ { \infty } [ - 1 , 1 ] , j = 0 , 1 , \ldots , \mu$ . In particular, let 

$$
f ^ {(k _ {0})} (x) = (1 - x) ^ {\gamma - k _ {0}} \left[ \ln^ {\mu} (1 - x) h _ {\mu} (x) + \ln^ {\mu - 1} (1 - x) h _ {\mu - 1} (x) \right.
$$

$$
\left. + \dots + \ln (1 - x) h _ {1} (x) + h _ {0} (x) \right]
$$

$$
=: (1 - x) ^ {\gamma - k _ {0}} \psi (x),
$$

which implies from $\alpha + \gamma > - 1$ that 

$$
(1 - x) ^ {k + \alpha} (1 + x) ^ {k + \beta} P _ {n - k} ^ {(k + \alpha , k + \beta)} (x) f ^ {(k - 1)} (x) \Big | _ {- 1} ^ {1} = 0, \quad k = 1, 2, \ldots , k _ {0}.
$$

Then by Rodrigues’ formula [15, p. 94, (4.10.1)] 

$$
\begin{array}{l} (1 - x) ^ {\alpha} (1 + x) ^ {\beta} P _ {n} ^ {(\alpha , \beta)} (x) \\ = \frac {(- 1) ^ {k}}{2 ^ {k} n (n - 1) \cdots (n - k + 1)} \frac {\mathrm{d} ^ {k}}{\mathrm{d} x ^ {k}} \left\{(1 - x) ^ {k + \alpha} (1 + x) ^ {k + \beta} P _ {n - k} ^ {(k + \alpha , k + \beta)} (x) \right\}, \\ \end{array}
$$

and integration by parts $k _ { 0 }$ times, we have 

$$
\begin{array}{l} a _ {n} (\alpha , \beta) = \frac {1}{\sigma_ {n} ^ {\alpha , \beta}} \int_ {- 1} ^ {1} (1 - x) ^ {\alpha} (1 + x) ^ {\beta} f (x) P _ {n} ^ {(\alpha , \beta)} (x)   \mathrm{d} x \tag {3.3} \\ = \frac {\int_ {- 1} ^ {1} (1 - x) ^ {\alpha + 1} (1 + x) ^ {\beta + 1} P _ {n - 1} ^ {(\alpha + 1 , \beta + 1)} (x) f ^ {\prime} (x) \mathrm{d} x}{2 n \sigma_ {n} ^ {\alpha , \beta}} \\ = \dots \\ = \frac {\int_ {- 1} ^ {1} (1 - x) ^ {\alpha + k _ {0}} (1 + x) ^ {\beta + k _ {0}} P _ {n - k _ {0}} ^ {(\alpha + k _ {0} , \beta + k _ {0})} (x) f ^ {(k _ {0})} (x) \mathrm{d} x}{2 ^ {k _ {0}} \sigma_ {n} ^ {\alpha , \beta} n (n - 1) \cdots (n - k _ {0} + 1)} \\ = \frac {\int_ {- 1} ^ {1} (1 - x) ^ {\alpha + \gamma} (1 + x) ^ {\beta + k _ {0}} P _ {n - k _ {0}} ^ {(\alpha + k _ {0} , \beta + k _ {0})} (x) \psi (x) \mathrm{d} x}{2 ^ {k _ {0}} \sigma_ {n} ^ {\alpha , \beta} n (n - 1) \cdots (n - k _ {0} + 1)}. \\ \end{array}
$$

Setting $x = \cos \theta$ and noting that $\begin{array} { r } { \ln ( 1 - x ) = 2 \ln ( \theta ) + 2 \ln \frac { \sin \theta / 2 } { \theta / 2 } - \ln 2 \mathrm { ~ a n d ~ } \ln \frac { \sin \theta / 2 } { \theta / 2 } \in \left( \frac { 1 } { \theta } , \frac { 1 } { \theta } \right) . } \end{array}$ n sin θ/2 − ln 2 and ln s $C ^ { \infty } [ 0 , \pi ] ;$ ; then $\psi ( x )$ can be represented as 

$$
\psi (\cos \theta) = \ln^ {\mu} (\theta) \hat {h} _ {\mu} (\cos \theta) + \ln^ {\mu - 1} (\theta) \hat {h} _ {\mu - 1} (\cos \theta) + \dots + \ln (\theta) \hat {h} _ {1} (\cos \theta) + \hat {h} _ {0} (\cos \theta)
$$

with $\hat { h } _ { j } \in C ^ { \infty } [ 0 , \pi ] , j = 0 , 1 , \ldots , \mu$ . Then by Lemma 2.1, Lemma 2.3, and (2.17), the numerator in (3.3) can be estimated as 

(3.4) 

$$
\int_ {- 1} ^ {1} (1 - x) ^ {\alpha + \gamma} (1 + x) ^ {\beta + k _ {0}} P _ {n - k _ {0}} ^ {(\alpha + k _ {0}, \beta + k _ {0})} (x) \psi (x) \mathrm{d} x
$$

$$
= 2 ^ {\alpha + \beta + \gamma + k _ {0} + 1} \int_ {0} ^ {\pi} \sin^ {2 \alpha + 2 \gamma + 1} \left(\frac {\theta}{2}\right) \cos^ {2 \beta + 2 k _ {0} + 1} \left(\frac {\theta}{2}\right) P _ {n - k _ {0}} ^ {(\alpha + k _ {0}, \beta + k _ {0})} (\cos \theta) \psi (\cos \theta) d \theta
$$

$$
= 2 ^ {\alpha + \beta + \gamma + k _ {0} + 1} \int_ {0} ^ {\pi / 2} \sin^ {2 \alpha + 2 \gamma + 1} \left(\frac {\theta}{2}\right) \cos^ {2 \beta + 2 k _ {0} + 1} \left(\frac {\theta}{2}\right) P _ {n - k _ {0}} ^ {(\alpha + k _ {0}, \beta + k _ {0})} (\cos \theta) \psi (\cos \theta) d \theta
$$

$$
+ 2 ^ {\alpha + \beta + \gamma + k _ {0} + 1} \int_ {\pi / 2} ^ {\pi} \sin^ {2 \alpha + 2 \gamma + 1} \left(\frac {\theta}{2}\right) \cos^ {2 \beta + 2 k _ {0} + 1} \left(\frac {\theta}{2}\right) P _ {n - k _ {0}} ^ {(\alpha + k _ {0}, \beta + k _ {0})} (\cos \theta) \psi (\cos \theta) d \theta
$$

$$
= 2 ^ {\alpha + \beta + \gamma + k _ {0} + 1} \left\{\int_ {0} ^ {\pi / 2} \theta^ {- \frac {1}{2}} \sin^ {\alpha + k _ {0} + \frac {1}{2}} \left(\frac {\theta}{2}\right) \cos^ {\beta + k _ {0} + \frac {1}{2}} \left(\frac {\theta}{2}\right) P _ {n - k _ {0}} ^ {(\alpha + k _ {0}, \beta + k _ {0})} (\cos \theta) \right.
$$

$$
\cdot \theta^ {\frac {1}{2}} \sin^ {\alpha + 2 \gamma - k _ {0} + \frac {1}{2}} \left(\frac {\theta}{2}\right) \cos^ {\beta + k _ {0} + \frac {1}{2}} \left(\frac {\theta}{2}\right) \psi (\cos \theta) d \theta
$$

$$
+ (- 1) ^ {n - k _ {0} + 1} \int_ {0} ^ {\pi / 2} \cos^ {2 \alpha + 2 \gamma + 1} \left(\frac {\theta}{2}\right) \sin^ {2 \beta + 2 k _ {0} + 1} \left(\frac {\theta}{2}\right) P _ {n - k _ {0}} ^ {(\beta + k _ {0}, \alpha + k _ {0})} (\cos \theta) \psi (- \cos \theta) \mathrm{d} \theta \Bigg \}
$$

$$
= \frac {\Gamma (n + \alpha + 1)}{(n - k _ {0}) ! \tilde {N} ^ {\alpha + k _ {0}}} \int_ {0} ^ {\pi / 2} \theta^ {\alpha + 2 \gamma - k _ {0} + 1} [ \ln^ {\mu} (\theta) \tilde {h} _ {\mu} (\theta) + \dots + \ln (\theta) \tilde {h} _ {1} (\theta) + \tilde {h} _ {0} (\theta) ] J _ {\alpha + k _ {0}} (\tilde {N} \theta) \mathrm{d} \theta
$$

$$
+ \mathcal {O} (\tilde {N} ^ {- \frac {3}{2}}) + \frac {(- 1) ^ {n - k _ {0} + 1} \Gamma (n + \beta + 1)}{(n - k _ {0}) ! \tilde {N} ^ {\beta + k _ {0}}} \int_ {0} ^ {\pi / 2} \theta^ {\beta + k _ {0} + 1} J _ {\beta + k _ {0}} (\tilde {N} \theta) \hat {\psi} _ {1} (\theta) \mathrm{d} \theta + \mathcal {O} (\tilde {N} ^ {- \frac {3}{2}})
$$

$$
= \mathcal {O} \left(\tilde {N} ^ {- \min \{\alpha + 2 \gamma - k _ {0} + 2, \frac {3}{2} \}} \ln^ {\mu} (\tilde {N})\right) + \mathcal {O} \left(\tilde {N} ^ {- \min \{\alpha + 2 \gamma - k _ {0} + 2, \frac {3}{2} \}} \ln^ {\mu - 1} (\tilde {N})\right)
$$

$$
+ \dots + \mathcal {O} \left(\tilde {N} ^ {- \min \{\alpha + 2 \gamma - k _ {0} + 2, \frac {3}{2} \}}\right) + \mathcal {O} \left(\tilde {N} ^ {- \min \{\beta + k _ {0} + 2, \frac {3}{2} \}}\right) + \mathcal {O} (\tilde {N} ^ {- \frac {3}{2}})
$$

$$
= \mathcal {O} \left(\tilde {N} ^ {- \alpha - 2 \gamma + k _ {0} - 2} \ln^ {\mu} (\tilde {N})\right),
$$

where 

$$
\left\{ \begin{array}{l} \tilde {h} _ {j} (\theta) = 2 ^ {\alpha + \beta + \gamma + k _ {0} + \frac {1}{2}} \left(\frac {\sin \theta / 2}{\theta}\right) ^ {\alpha + 2 \gamma - k _ {0} + \frac {1}{2}} \cos^ {\beta + k _ {0} + \frac {1}{2}} \left(\frac {\theta}{2}\right) \hat {h} _ {j} (\theta), \\ \hat {\psi} _ {1} (\theta) = 2 ^ {\alpha + \beta + \gamma + k _ {0} + \frac {1}{2}} \left(\frac {\sin \theta / 2}{\theta}\right) ^ {\beta + k _ {0} + \frac {1}{2}} \cos^ {\alpha + 2 \gamma - k _ {0} + \frac {1}{2}} \left(\frac {\theta}{2}\right) \psi (- \cos \theta), \end{array} \right.
$$

and $\tilde { N } = n + ( \alpha + \beta + 1 ) / 2$ . It is easy to verify that $\begin{array} { r } { \hat { \psi } _ { 1 } ( \theta ) , \tilde { h } _ { j } ( \theta ) \in C ^ { \infty } [ 0 , \frac { \pi } { 2 } ] } \end{array}$ for $j = 0 , \ldots , \mu ;$ then together with $\sigma _ { n } ^ { \alpha , \beta } = \mathcal { O } ( n ^ { - 1 } )$ , it derives 

$$
\begin{array}{l} a _ {n} (\alpha , \beta) = \frac {\int_ {- 1} ^ {1} (1 - x) ^ {\alpha + \gamma} (1 + x) ^ {\beta + k _ {0}} P _ {n - k _ {0}} ^ {(\alpha + k _ {0} , \beta + k _ {0})} (x) \psi (x) \mathrm{d} x}{2 ^ {k _ {0}} \sigma_ {n} ^ {\alpha , \beta} n (n - 1) \cdots (n - k _ {0} + 1)} \\ = \mathcal {O} \left(n ^ {- \alpha - 2 \gamma - 1} \ln^ {\mu} (n)\right). \\ \end{array}
$$

In the case γ is a nonnegative integer and $\mu > 1$ , we select $k _ { 0 }$ further satisfied $k _ { 0 } > \operatorname* { m a x } \{ \gamma , \mu \}$ ; then $f ^ { ( k _ { 0 } ) } ( x )$ can be represented as 

$$
f ^ {(k _ {0})} (x) = \sum_ {j = 0} ^ {k _ {0}} (1 - x) ^ {\gamma - k _ {0} + j} \sum_ {i = 0} ^ {\mu - 1} \ln^ {i} (1 - x) \phi_ {j, i} (x) + (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) \phi_ {0} (x)
$$

with $\phi _ { 0 } \in C ^ { \infty } [ - 1 , 1 ]$ and $\phi _ { j , i } \in C ^ { \infty } [ - 1 , 1 ]$ . In analogy to the above proof, it also leads to the desired result. 

In the case γ is a nonnegative integer and $\mu = 1 , f ^ { ( k _ { 0 } ) } ( x )$ can be represented as 

$$
f ^ {(k _ {0})} (x) = \ln (1 - x) \psi_ {1} (x) + (1 - x) ^ {\gamma - k _ {0}} \psi_ {2} (x), \quad \psi_ {1}, \psi_ {2} \in C ^ {\infty} [ - 1, 1 ].
$$

In analogy to the above proof, it leads to 

$$
\int_ {- 1} ^ {1} (1 - x) ^ {\alpha + k _ {0}} (1 + x) ^ {\beta + k _ {0}} P _ {n - k _ {0}} ^ {(\alpha + k _ {0}, \beta + k _ {0})} (x) \ln (1 - x) \psi_ {1} (x) \mathrm{d} x = \mathcal {O} (n ^ {- \frac {3}{2}})
$$

by $\alpha + k _ { 0 } + 1 > \frac { 1 } { 2 }$ , and 

$$
\int_ {- 1} ^ {1} (1 - x) ^ {\alpha + \gamma} (1 + x) ^ {\beta + k _ {0}} P _ {n - k _ {0}} ^ {(\alpha + k _ {0}, \beta + k _ {0})} (x) \psi_ {2} (x) \mathrm{d} x = \mathcal {O} \left(n ^ {- \alpha - 2 \gamma + k _ {0} - 2}\right),
$$

which together deduces the desired result (3.2). 

Similar results can be obtained for 

$$
f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \tag {3.5}
$$

where µ is a positive integer and $g \in C ^ { \infty } [ - 1 , 1 ]$ . 

Theorem 3.2. Suppose f(x) is defined by (3.5); then the Jacobi coefficients (2.2) satisfy that for $\beta + \delta > - 1$ 

$$
| a _ {n} (\alpha , \beta) | = \left\{ \begin{array}{l l} \mathcal {O} \left(n ^ {- \beta - 2 \delta - 1} \ln^ {\mu} (n)\right), & \delta \text {   is   not   an   integer }, \\ \mathcal {O} \left(n ^ {- \beta - 2 \delta - 1} \ln^ {\mu - 1} (n)\right), & \delta \text {   is   an   integer }. \end{array} \right. \tag {3.6}
$$

Now we consider 

$$
f (x) = (1 - x) ^ {\gamma} (1 + x) ^ {\delta} \ln^ {\mu} (1 - x ^ {2}) g (x), \tag {3.7}
$$

where $\mu$ is a positive integer and $g \in C ^ { \infty } [ - 1 , 1 ]$ . Then from Tuan and Elloiit [17], $f ( x )$ can be rewritten as 

$$
f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g _ {1} (x) + (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g _ {2} (x)
$$

with $g _ { 1 } , g _ { 2 } \in C ^ { \infty } [ - 1 , 1 ]$ . 

Corollary 3.3. Suppose $f ( x )$ is defined as (3.7); then the Jacobi coefficients (2.2) satisfy that 

(3.8) 

$$
| a _ {n} (\alpha , \beta) | = \left\{ \begin{array}{l l} \mathcal {O} \left(n ^ {- \min \{\alpha + 2 \gamma + 1, \beta + 2 \delta + 1 \}} \ln^ {\mu - 1} (n)\right), & \gamma , \delta \in \mathbf {Z} _ {+}, \\ \mathcal {O} \left(\max \left\{n ^ {- \alpha - 2 \gamma - 1}, n ^ {- \beta - 2 \delta - 1} \ln (n) \right\} \ln^ {\mu - 1} (n)\right), & \gamma \in \mathbf {Z} _ {+}, \delta \notin \mathbf {Z} _ {+}, \\ \mathcal {O} \left(\max \left\{n ^ {- \beta - 2 \delta - 1} \ln (n), n ^ {- \alpha - 2 \gamma - 1} \right\} \ln^ {\mu - 1} (n)\right), & \gamma \notin \mathbf {Z} _ {+}, \delta \in \mathbf {Z} _ {+}, \\ \mathcal {O} \left(n ^ {- \min \{\alpha + 2 \gamma + 1, \beta + 2 \delta + 1 \}} \ln^ {\mu} (n)\right), & \gamma \notin \mathbf {Z} _ {+}, \delta \notin \mathbf {Z} _ {+}, \end{array} \right.
$$

where min $\{ \alpha + \gamma , \beta + \delta \} > - 1$ and $\mathbf { Z } _ { + }$ is the set of nonnegative integers. 

We illustrate the decay rates by $f ( x ) = ( 1 - x ) ^ { \gamma } \ln ( 1 - x )$ with different values of α and $\beta$ (see Figure 3.1). These numerical results are in accordance with the estimates. The asymptotic orders of the decay on the coefficients are sharp. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/d62db15e68714268e6a18585bf894db4fa318eabea6fdab5fb2d800e1c6bf848.jpg)



Fig. 3.1. The asymptotic decay of the Jacobi coefficients $\left| a _ { n } ( \alpha , \beta ) \right|$ in (2.1) for $f ( x ) \ =$ $( 1 - x ) ^ { \gamma } \ln ( 1 - x )$ with $\gamma = 0$ (first row), $\gamma = 0 . 5$ (second row), $\gamma = 1$ (third row), and $d i f f e r _ { - }$ ent values of $( \alpha , \beta )$ .


One can expect that sharp bounds for the Gegenbauer coefficients (2.5) and Chebyshev coefficients (2.9) from (2.8) and (2.10), respectively. 

Corollary 3.4. The Gegenbauer and Chebyshev expansion coefficients for f(x) are satisfied 

$\bullet f ( x ) = ( 1 - x ) ^ { \gamma } \ln ^ { \mu } ( 1 - x ) g ( x ) \ a n d \ g \in C ^ { \infty } [ - 1 , 1 ] .$ 

$$
| a _ {n} (\lambda) | = \left\{ \begin{array}{l l} \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- 2 \lambda - 2 \gamma}\right), & \lambda + \gamma > - \frac {1}{2} a n d \gamma \notin \mathbf {Z} _ {+}, \\ \mathcal {O} \left(\ln^ {\mu - 1} (n) n ^ {- 2 \lambda - 2 \gamma}\right), & \lambda + \gamma > - \frac {1}{2} a n d \gamma \in \mathbf {Z} _ {+}, \end{array} \right. \tag {3.9}
$$

$$
| c _ {n} | = \left\{ \begin{array}{l l} \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- 1 - 2 \gamma}\right), & \gamma > - \frac {1}{2} a n d \gamma \notin \mathbf {Z} _ {+}, \\ \mathcal {O} \left(\ln^ {\mu - 1} (n) n ^ {- 1 - 2 \gamma}\right), & \gamma > - \frac {1}{2} a n d \gamma \in \mathbf {Z} _ {+}. \end{array} \right. \tag {3.10}
$$

$\bullet f ( x ) = ( 1 + x ) ^ { \delta } \ln ^ { \mu } ( 1 + x ) g ( x ) \ a n d \ g \in C ^ { \infty } [ - 1 , 1 ] .$ 

$$
(3. 1 1) \quad | a _ {n} (\lambda) | = \left\{ \begin{array}{l l} \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- 2 \lambda - 2 \delta}\right), & \lambda + \delta > - \frac {1}{2} a n d \delta \notin \mathbf {Z} _ {+}, \\ \mathcal {O} \left(\ln^ {\mu - 1} (n) n ^ {- 2 \lambda - 2 \delta}\right), & \lambda + \delta > - \frac {1}{2} a n d \delta \in \mathbf {Z} _ {+}, \end{array} \right.
$$

$$
| c _ {n} | = \left\{ \begin{array}{l l} \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- 1 - 2 \delta}\right), & \delta > - \frac {1}{2} a n d \delta \notin \mathbf {Z} _ {+}, \\ \mathcal {O} \left(\ln^ {\mu - 1} (n) n ^ {- 1 - 2 \delta}\right), & \delta > - \frac {1}{2} a n d \delta \in \mathbf {Z} _ {+}. \end{array} \right. \tag {3.12}
$$

Remark 2. (i) Theorems 3.1–3.2 and Corollaries 3.3–3.4 show that for functions of endpoint singularities, the decay of the coefficients in a Legendre polynomial series has the same asymptotic order as the Chebyshev. 

(ii) One may get faster convergence rate by increasing $( \alpha , \beta )$ for endpoint singularities. In the case min $\{ \alpha , \beta \} > 0 \mathrm { o r } \lambda > { \frac { 1 } { 2 } }$ , the decay order of coefficients in a Jacobi or Gegenbauer series is faster than that in a Legendre or Chebyshev series. 

Numerical results for these estimates on the two end points about the Chebyshev, Legendre and Jacobi expansion coefficients are illustrated in Figures 3.2 and 3.3, which indicate the optimal orders of the estimates. To get more clearly, in the second row of Figure 3.2, we consider $n = 1 0 : 1 0 : 2 , 0 0 0$ . 

# 3.2. Functions with interior regularities. Let us consider

$$
f (x) = | x - z _ {0} | ^ {s} \ln^ {\mu} | x - z _ {0} | g (x), \quad z _ {0} \in (- 1, 1) \tag {3.13}
$$

with $g \in C ^ { \infty } [ - 1 , 1 ] , s > 0$ a real number, and µ a positive integer. 

Theorem 3.5. Suppose that $f ( x )$ is defined as (3.13); then the Jacobi coefficients (2.2) satisfy that 

$$
| a _ {n} (\alpha , \beta) | = \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- s - \frac {1}{2}}\right) \quad \text { as   } n \to \infty . \tag {3.14}
$$

Proof. Note that 

$$
a _ {n} (\alpha , \beta) = \frac {1}{\sigma_ {n} ^ {\alpha , \beta}} \left[ \int_ {- 1} ^ {z _ {0}} + \int_ {z _ {0}} ^ {1} \right] (1 - x) ^ {\alpha} (1 + x) ^ {\beta} P _ {n} ^ {(\alpha , \beta)} (x) f (x) \mathrm{d} x. \tag {3.15}
$$

Without loss of generality, here we consider only the second integral in (3.15). 

Let $k _ { 0 } = \lfloor s \rfloor$ , the greatest integer less than or equals to s. Since 

$$
\begin{array}{l} f ^ {(k _ {0})} (x) = (x - z _ {0}) ^ {s - k _ {0}} \left[ \ln^ {\mu} (x - z _ {0}) h _ {\mu} (x) + \ln^ {\mu - 1} (x - z _ {0}) h _ {\mu - 1} (x) \right. \\ \left. + \dots + \ln (x - z _ {0}) h _ {1} (x) + h _ {0} (x) \right] \\ =: (x - z _ {0}) ^ {s - k _ {0}} \psi (x) \\ \end{array}
$$

with $h _ { j } \in C ^ { \infty } [ z _ { 0 } , 1 ] , j = 0 , 1 , \ldots , \mu .$ . 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/581151672cdace2fa9da84ae46002ddaead2a6303713cb5ff31e26dd96b14de8.jpg)



Fig. 3.2. The asymptotic decay of the Chebyshev, Legendre, and Jacobi expansion coefficients for $f ( x ) = ( 1 - x ) ^ { 0 . 3 } \stackrel {  } { ( 1 + x ) } ^ { 0 . 7 } \ln ^ { \mu } \stackrel {  } { ( 1 - x ^ { 2 } ) }$ sin x with different values of $( \alpha , \beta )$ , respectively: $\mu = 1$ (first row) and $\mu = 2$ (second row).


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/73cc2bf60f7182db1b59164e846ff2ab27ff40db7a82f273afdfb85bab1f62cf.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/3e92111a6584165883ff7dc6e95a39053777156a63eaef91b94800a758b3d507.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/8a8b5b187222e3cc2d24484887ad000ba64e6354dd4e8017a7516ba43dac037e.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/9f7cde984251f17013fb7c1b2bee7fce2713b3fc22143459ba07c6d4058e457f.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/5b441d92b0d2ecf313deef5617b046a36eb990ceb75ccdf4b07453116f5e019c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/dc9a7bbc900a5d7f9c3c69c70c966355aa368634685f5235d4ed6f6b2b3357f9.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/21f6284066c7e29ef787f7218746cb9c39aeed017741da7a21a9ecf380da483e.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/205d943ffd4198e4ffd30fd18e31f4da74c601ac9792c5f7173a2974c7ad7111.jpg)



Fig. 3.3. The asymptotic decay of the Chebyshev, Legendre, and Jacobi expansion coefficients $f o r \ f ( x ) = ( 1 - x ) ^ { 1 } ( 1 + \bar { x ) } ^ { 2 } \ln ^ { \mu } ( 1 - \bar { x } ^ { 2 } )$ ) sin x with different values of $( \alpha , \beta )$ , respectively: $\mu = 1$ (first row) and $\mu = 2$ (second row).


Applying Rodrigues’ formula follows 

$$
I _ {1} = \frac {\int_ {- 1} ^ {1} (1 - x) ^ {\alpha + k _ {0}} (1 + x) ^ {\beta + k _ {0}} P _ {n - k _ {0}} ^ {(\alpha + k _ {0} , \beta + k _ {0})} (x) f ^ {(k _ {0})} (x)   \mathrm{d} x}{2 ^ {k _ {0}} \sigma_ {n} ^ {\alpha , \beta} n (n - 1) \cdots (n - k _ {0} + 1)}.
$$

For simplicity, we only consider the first term in $f ^ { ( k _ { 0 } ) } ( x )$ denoted by $\varphi _ { 1 } ( x ) = ( x -$ $z _ { 0 } ) ^ { s - k _ { 0 } } \ln ^ { \mu } ( x - z _ { 0 } ) \psi _ { 1 } ( x )$ . Similar proof can be directly applied to estimate the other terms. Without loss of generality, assume $f ^ { ( k _ { 0 } ) } ( x ) = \bar { \varphi _ { 1 } ( x ) } = ( x - z _ { 0 } ) ^ { s - k _ { 0 } } \ln ^ { \mu } ( x -$ $z _ { 0 } ) \psi _ { 1 } ( x )$ . 

In the case $s = k _ { 0 }$ is a positive integer: Setting $x = \cos \theta$ and $\theta _ { 0 } = \operatorname { a r c c o s } z _ { 0 }$ together with Lemma 2.1, analogously to the proofs of Theorem 3.1, it indicates that by Lemma 2.9 

$$
\begin{array}{l} I _ {1} = \frac {\int_ {0} ^ {\theta_ {0}} (1 - \cos \theta) ^ {\alpha + k _ {0}} (1 + \cos \theta) ^ {\beta + k _ {0}} P _ {n - k _ {0}} ^ {(\alpha + k _ {0} , \beta + k _ {0})} (\cos \theta) \varphi_ {1} (\cos \theta) \sin \theta \mathrm{d} \theta}{2 ^ {k _ {0}} \sigma_ {n} ^ {\alpha , \beta} n (n - 1) \cdots (n - k _ {0} + 1)} \\ = \frac {\Gamma (n + \alpha + 1) \int_ {0} ^ {\theta_ {0}} \left(\frac {\theta}{2}\right) ^ {\frac {1}{2}} \sin^ {\alpha + k _ {0} + \frac {1}{2}} \frac {\theta}{2} \cos^ {\beta + k _ {0} + \frac {1}{2}} \frac {\theta}{2} J _ {\alpha + k _ {0}} (\tilde {N} \theta) \varphi_ {1} (\cos \theta) d \theta}{2 ^ {- \alpha - \beta - 1 - k _ {0}} (n - k _ {0}) ! \tilde {N} ^ {\alpha + k _ {0}} \sigma_ {n} ^ {\alpha , \beta} n (n - 1) \cdots (n - k _ {0} + 1)} \\ + \mathcal {O} (\tilde {N} ^ {- 3 / 2 - k _ {0} + 1}) \\ = \frac {2 ^ {\alpha + \beta + 1 + k _ {0}} \Gamma (n + \alpha + 1)}{n ! \tilde {N} ^ {\alpha + k _ {0}} \sigma_ {n} ^ {\alpha , \beta}} \int_ {0} ^ {\theta_ {0}} \theta^ {\alpha + k _ {0} + 1} (\theta_ {0} - \theta) ^ {s - k _ {0}} \left\{\ln^ {\mu} (\theta_ {0} - \theta) \hat {\varphi} _ {\mu} (\theta) \right. \\ + \dots + \ln (\theta_ {0} - \theta) \hat {\varphi} _ {1} (\theta) + \hat {\varphi} _ {0} (\theta) \} J _ {\alpha + k _ {0}} (\tilde {N} \theta) d \theta + \mathcal {O} \left(\tilde {N} ^ {- \frac {1}{2} - k _ {0}}\right) \\ = \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- s - \frac {1}{2}}\right) + \dots + \mathcal {O} \left(\ln (n) n ^ {- s - \frac {1}{2}}\right) + \mathcal {O} \left(n ^ {- s - \frac {1}{2}}\right) + \mathcal {O} \left(n ^ {- \frac {1}{2} - k _ {0}}\right) \\ = \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- s - \frac {1}{2}}\right) \\ \end{array}
$$

due to $\alpha + k _ { 0 } + 2 > { \frac { 3 } { 2 } }$ , where $\tilde { N } = n + ( \alpha + \beta + 1 ) / 2$ and $\hat { \varphi } _ { j } ( \theta ) \in C ^ { \infty } [ 0 , \theta _ { 0 } ]$ for $j = 0 , 1 , \dotsc , \mu .$ 

In the case that $s > 0$ is not an integer: By integrating by parts once again, it follows 

$$
\begin{array}{l} I _ {1} = \frac {\int_ {0} ^ {\theta_ {0}} (1 - \cos \theta) ^ {\alpha + k _ {0} + 1} (1 + \cos \theta) ^ {\beta + k _ {0} + 1} P _ {n - k _ {0} - 1} ^ {(\alpha + k _ {0} + 1 , \beta + k _ {0} + 1)} (\cos \theta) \varphi_ {1} ^ {\prime} (\cos \theta) \sin \theta \mathrm{d} \theta}{2 ^ {k _ {0} + 1} \sigma_ {n} ^ {\alpha , \beta} n (n - 1) \cdots (n - k _ {0})} \\ = \frac {\Gamma (n + \alpha + 1) \int_ {0} ^ {\theta_ {0}} \left(\frac {\theta}{2}\right) ^ {\frac {1}{2}} \sin^ {\alpha + k _ {0} + \frac {3}{2}} \frac {\theta}{2} \cos^ {\beta + k _ {0} + \frac {3}{2}} \frac {\theta}{2} J _ {\alpha + k _ {0} + 1} (\tilde {N} \theta) \varphi_ {1} ^ {\prime} (\cos \theta) d \theta}{2 ^ {- \alpha - \beta - 2 - k _ {0}} (n - k _ {0} - 1) ! \tilde {N} ^ {\alpha + k _ {0}} \sigma_ {n} ^ {\alpha , \beta} n (n - 1) \cdots (n - k _ {0})} \\ + \mathcal {O} (\tilde {N} ^ {- 3 / 2 - k _ {0}}) \\ = \frac {2 ^ {\alpha + \beta + 2 + k _ {0}} \Gamma (n + \alpha + 1)}{n ! \tilde {N} ^ {\alpha + k _ {0} + 1} \sigma_ {n} ^ {\alpha , \beta}} \int_ {0} ^ {\theta_ {0}} \theta^ {\alpha + k _ {0} + 2} (\theta_ {0} - \theta) ^ {s - k _ {0} - 1} \{\ln^ {\mu} (\theta_ {0} - \theta) \hat {\varphi} _ {\mu} (\theta) \\ + \dots + \ln (\theta_ {0} - \theta) \hat {\varphi} _ {1} (\theta) + \hat {\varphi} _ {0} (\theta) \} J _ {\alpha + k _ {0} + 1} (\tilde {N} \theta) \mathrm{d} \theta + \mathcal {O} (\tilde {N} ^ {- \frac {3}{2} - k _ {0}}) \\ = \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- s - \frac {1}{2}}\right) + \dots + \mathcal {O} \left(\ln (n) n ^ {- s - \frac {1}{2}}\right) + \mathcal {O} \left(n ^ {- s - \frac {1}{2}}\right) + \mathcal {O} \left(n ^ {- \frac {3}{2} - k _ {0}}\right) \\ = \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- s - \frac {1}{2}}\right) \\ \end{array}
$$

by Lemma 2.9 and $\alpha + 2 k _ { 0 } + 3 > k _ { 0 } + \textstyle \frac { 3 } { 2 } > s + \textstyle \frac { 1 } { 2 }$ 

Similar results can be derived for the first integral in the right-hand side of (3.15). These lead to the desired result (3.14). □ 

Remark 3. In the case $- 1 < s \le 0$ , in analogy to the above proof without integrating by parts, from Lemma 2.9 we get that 

$$
| a _ {n} (\alpha , \beta) | = \mathcal {O} \left(\max \left\{\ln^ {\mu} (n) n ^ {- s - \frac {1}{2}}, n ^ {- \min \{1 + \alpha , 1 + \beta \}} \right\}\right) \quad \text { as } n \to \infty . \tag {3.16}
$$

Corollary 3.6. Suppose that f(x) is defined as (3.13); then for the Gegenbauer and Chebyshev expansion coefficients, the following hold 

$$
| a _ {n} (\lambda) | = \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- \lambda - s}\right), \tag {3.17}
$$

$$
| c _ {n} | = \mathcal {O} \left(\ln^ {\mu} (n) n ^ {- 1 - s}\right). \tag {3.18}
$$

Remark 4. (i) For functions of interior singularities, the decay order of the coefficients in a Jacobi polynomial series is independent of $( \alpha , \beta )$ if $s > 0$ . 

(ii) The decay of the coefficients in a Chebyshev polynomial series has the highest asymptotic order. Figure 3.4 illustrates this phenomenon exactly. However, in section 4, we will see that all these spectral expansions have the same convergence order. 

Remark 5. From the proofs, we see that all the above thereoms and corollaries still hold for $\mu = 0 ~ [ 2 5 ]$ except for s even. In the case $\mu = 0$ and s even, $a _ { n } ( \alpha , \beta )$ , $a _ { n } ( \lambda )$ , and $c _ { n }$ are exponentially decayed [5, 16, 18, 19, 22, 26, 27]. 

4. The convergence rates on the spectral orthogonal projections. For $f ( x ) = ( 1 - x ) ^ { \gamma } \ln ^ { \mu } ( 1 - x ) g ( x )$ or $f ( x ) = ( 1 + x ) ^ { \delta } \ln ^ { \mu } ( 1 + x ) g ( x )$ , it is easy to verify that $f \in L _ { w ( \alpha , \beta ) } ^ { 2 } [ - 1 , 1 ]$ if min $\{ \alpha + 2 \gamma , \beta + 2 \delta \} > - 1$ . For $f ( x ) = | x - z _ { 0 } | ^ { s } \ln ^ { \mu } | x - z _ { 0 } | g ( x )$ $f \in L _ { w ( \alpha , \beta ) } ^ { 2 } [ - 1 , 1 ] \mathrm { ~ i f ~ } s > - \frac { 1 } { 2 }$ . 

Theorem 4.1. Suppose that f is of algebraic and logarithmatic regularity at an endpoint or interior point; then for the Jacobi expansion, it follows that 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/ecd1abe6905bf8ef9e8602561837dbeb4629bcdcce43cf337aa3e8d10b278a6a.jpg)



Fig. 3.4. The asymptotic decay of the Chebyshev, Legendre, Jacobi coefficients for $f ( x ) =$ $| x - 0 . 5 | ^ { s } \ln | x - 0 . 5 |$ cos x with $s = 0 . 5$ and $s = 3 ,$ , respectively: $n = 1 0 : 1 0 : 1 , 0 0 0$ .


(4.1) 

$$
\| f - \mathcal {P} _ {N} ^ {f, J a} \| _ {L _ {w (\alpha , \beta)} ^ {2} [ - 1, 1 ]} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {- \alpha - 2 \gamma - 1} \ln^ {\mu} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {- s - 1 / 2} \ln^ {\mu} (N)), & f (x) = | x - z _ {0} | ^ {s} \ln^ {\mu} | x - z _ {0} | g (x), \\ \mathcal {O} (N ^ {- \beta - 2 \delta - 1} \ln^ {\mu} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \end{array} \right.
$$

where $g \in C ^ { \infty } [ - 1 , 1 ]$ , µ is a nonnegative integer, $z _ { 0 } \in ( - 1 , 1 )$ , min $\{ \alpha + \gamma , \beta + \delta , \alpha +$ $2 \gamma , \beta + 2 \delta \} > - 1$ for the boundary singularities, and $s > - \frac { 1 } { 2 }$ and m $\begin{array} { r } { \operatorname { i n } \{ \alpha , \beta \} \ge - \frac { 1 } { 2 } ~ f o r } \end{array}$ the interior singularity. In particular, $i f \gamma , \delta$ are integers and $\mu \geq 1$ , then 

(4.2) 

$$
\| f - \mathcal {P} _ {N} ^ {f, J a} \| _ {L _ {w (\alpha , \beta)} ^ {2} [ - 1, 1 ]} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {- \alpha - 2 \gamma - 1} \ln^ {\mu - 1} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {- \beta - 2 \delta - 1} \ln^ {\mu - 1} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x). \end{array} \right.
$$

While the optimal convergence rate on the orthogonal projection can be obtained from the sharp bounds for the Gegenbauer coefficients and 

$$
\| f - \mathcal {P} _ {N} ^ {f, G e} \| _ {L _ {w (\lambda)} ^ {2} [ - 1, 1 ]} = \sqrt {\sum_ {n = N + 1} ^ {\infty} a _ {n} ^ {2} (\lambda) \hbar_ {n}}
$$

together with 

$$
\hbar_ {n} = \frac {2 ^ {1 - 2 \lambda} \pi}{\Gamma^ {2} (\lambda)} \frac {\Gamma (n + 2 \lambda)}{n ! (n + \lambda)} = \mathcal {O} (n ^ {2 \lambda - 2})
$$

(see [10, p. 79] and ([15, (7.33.1) p. 171]). 

Corollary 4.2. Suppose that f is of algebraic and logorithmatic regularity at an endpoint or interior point; then for the Gegenbauer and Chebyshev expansions, it follows that 

(4.3) 

$$
\| f - \mathcal {P} _ {N} ^ {f, G e} \| _ {L _ {w (\lambda)} ^ {2} [ - 1, 1 ]} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {- \lambda - 2 \gamma - 1 / 2} \ln^ {\mu} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {- s - 1 / 2} \ln^ {\mu} (N)), & f (x) = | x - z _ {0} | ^ {s} \ln^ {\mu} | x - z _ {0} | g (x), \\ \mathcal {O} (N ^ {- \lambda - 2 \delta - 1 / 2} \ln^ {\mu} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \end{array} \right.
$$

(4.4) 

$$
\| f - \mathcal {P} _ {N} ^ {f, C h} \| _ {L _ {w (- \frac {1}{2}, - \frac {1}{2})} ^ {2} [ - 1, 1 ]} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {- 1 / 2 - 2 \gamma} \ln^ {\mu} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {- s - 1 / 2} \ln^ {\mu} (N)), & f (x) = | x - z _ {0} | ^ {s} \ln^ {\mu} | x - z _ {0} | g (x), \\ \mathcal {O} (N ^ {- 1 / 2 - 2 \delta} \ln^ {\mu} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \end{array} \right.
$$

where $\begin{array} { r } { g \in C ^ { \infty } [ - 1 , 1 ] , \ z _ { 0 } \in ( - 1 , 1 ) , \operatorname* { m i n } \{ \lambda + \gamma , \lambda + 2 \gamma , \lambda + \delta , , \lambda + 2 \delta \} > - \frac { 1 } { 2 } , \ s > - \frac { 1 } { 2 } } \end{array}$ with $\lambda \geq 0$ . In particular, $i \dot { f } \gamma , \delta$ are integers and $\mu \geq 1$ , then 

(4.5) 

$$
\| f - \mathcal {P} _ {N} ^ {f, G e} \| _ {L _ {w (\lambda)} ^ {2} [ - 1, 1 ]} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {- \lambda - 2 \gamma - 1 / 2} \ln^ {\mu - 1} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {- \lambda - 2 \delta - 1 / 2} \ln^ {\mu - 1} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \end{array} \right.
$$

$$
\| f - \mathcal {P} _ {N} ^ {f, C h} \| _ {L _ {w (- \frac {1}{2}, - \frac {1}{2})} ^ {2} [ - 1, 1 ]} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {- 1 / 2 - 2 \gamma} \ln^ {\mu - 1} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {- 1 / 2 - 2 \delta} \ln^ {\mu - 1} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x). \end{array} \right.
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/0ce1762cc5f577cd6035330aa3caff4e8d122e9b8593b2c051e30e19a4ec5ae1.jpg)



Fig. 4.1. The weighted norm errors of the truncated Chebyshev, Legendre, and Jacobi $( \alpha = 3 . 6 ,$ $\beta = 3 . 7 )$ expansions $\| f - \mathcal { P } _ { N } ^ { f } \| _ { L _ { \omega ( \alpha , \beta ) } ^ { 2 } [ - 1 , 1 ] }$ L2 [−1,1] for f (x) = (1 − x)γ (1 + x)δ ln(1 − x2): γ = 0.6 and $f ( x ) = ( 1 - x ) ^ { \gamma } ( 1 + x ) ^ { \delta } \ln ( 1 - x ^ { 2 } ) \colon \gamma = 0 . 6$ $\delta = 0 . 4$ (first row), and $\gamma = 1$ and $\delta = 2$ (second row), respectively.


Remark 6. From Theorem 4.1 and Corollary 4.2, we see that for functions with limited regularities at the endpoints $f ( x ) = ( \mathrm { { 1 } } - x ) ^ { \gamma } ( \mathrm { { 1 } } + x ) ^ { \delta } \ln ^ { \mu } ( \mathrm { { 1 } } - x ^ { 2 } ) g ( x )$ , the Jacobi or Gegenbauer expansion can achieve faster convergence order than Chebyshev expansion if min $\displaystyle \{ \alpha , \beta \} > - \frac { 1 } { 2 }$ or $\lambda > 0$ , which is illustrated by Figures 4.1 and 4.2. 

Remark 7. It is obvious from Theorem 4.1 and Corollary 4.2 that for the functions with interior singularity, all these spectral expansions have the same convergence order in the case $s > - { \frac { 1 } { 2 } }$ with min $\{ \alpha , \beta \} \ge - { \textstyle \frac { 1 } { 2 } } \mathrm { o r } \lambda \ge 0$ . All these estimates are attainable. 

In particular, for the nonuniformly Jacobi-weighted Sobolev space $H ^ { m , \alpha , \beta } ( \Omega )$ with integer $m \ge 0 , \alpha > - 1 , \beta > - 1 , \Omega = [ - 1 , 1 ]$ and weighted norm 

$$
\| u \| _ {H ^ {m, \alpha , \beta} (\Omega)} = \left\{\sum_ {q = 0} ^ {m} \int_ {- 1} ^ {1} (1 - x) ^ {\alpha + q} (1 + x) ^ {\beta + q} [ u ^ {(q)} (x) ] ^ {2} d x \right\} ^ {\frac {1}{2}}, \tag {4.7}
$$

define 

$$
f ^ {(q)} (x) = \sum_ {n = 0} ^ {\infty} a _ {n} ^ {(q)} (\alpha + q, \beta + q) P _ {n} ^ {(\alpha + q, \beta + q)} (x), \quad q = 0, 1, \dots , m \tag {4.8}
$$

with 

$$
a _ {n} ^ {(q)} (\alpha + q, \beta + q) = \frac {1}{\sigma_ {n} ^ {\alpha + q , \beta + q}} \int_ {- 1} ^ {1} (1 - x) ^ {\alpha + q} (1 + x) ^ {\beta + q} f ^ {(q)} (x) P _ {n} ^ {(\alpha + q, \beta + q)} (x) d x.
$$

From (3.3), we see that 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/500d1a00e9b528c855d46217bd81168920cc4ac487ee2715c6b25180e820f666.jpg)



Fig. 4.2. The weighted norm errors of the truncated Chebyshev, Legendre, and Jacobi $( \alpha = 3 . 6 ,$ $\beta = 3 . 7 )$ expansions $\| f - \mathcal { P } _ { N } ^ { f } \| _ { L _ { \omega ( \alpha , \beta ) } ^ { 2 } [ - 1 , 1 ] }$ L2 [−1,1] for f (x) = |x − 12 |s ln |x − 12 |: s = 1 (first row) and $f ( x ) = | x - { \textstyle \frac { 1 } { 2 } } | ^ { s } \ln | x - { \textstyle \frac { 1 } { 2 } } | \colon s = 1$ $s = 2 . 5$ (second row), respectively.


$$
\begin{array}{l} a _ {n} ^ {(q)} (\alpha + q, \beta + q) \\ = \frac {1}{\sigma_ {n} ^ {\alpha + q , \beta + q}} \int_ {- 1} ^ {1} (1 - x) ^ {\alpha + q} (1 + x) ^ {\beta + q} f ^ {(q)} (x) P _ {n} ^ {(\alpha + q, \beta + q)} (x) d x \\ = \frac {2 ^ {q}}{\sigma_ {n} ^ {\alpha + q , \beta + q}} (n + q) (n + q - 1) \dots (n + 1) \int_ {- 1} ^ {1} (1 - x) ^ {\alpha} (1 + x) ^ {\beta} f (x) P _ {n + q} ^ {(\alpha , \beta)} (x) d x \\ = \frac {\sigma_ {n + q} ^ {\alpha , \beta}}{\sigma_ {n} ^ {\alpha + q , \beta + q}} 2 ^ {q} (n + q) (n + q - 1) \dots (n + 1) a _ {n + q} (\alpha , \beta) \\ = a _ {n + q} (\alpha , \beta) \mathcal {O} (n ^ {q}), \\ \end{array}
$$

which together with Theorems 3.1–3.2 and 3.5 yields the following convergence rate. 

Notice that for $f ( x ) = ( 1 - x ) ^ { \gamma } \ln ^ { \mu } ( 1 - x ) g ( x ) { \mathrm { ~ o r ~ } } f ( x ) = ( 1 + x ) ^ { \delta } \ln ^ { \mu } ( 1 + x ) g ( x )$ , it is easy to verify that $f \in H ^ { m , \tilde { \alpha } , \beta } ( \Omega )$ if min $\{ \alpha + 2 \gamma - m , \beta + 2 \delta - m \} > - 1$ . For $\begin{array} { r } { f ( x ) = | x - z _ { 0 } | ^ { s } \ln ^ { \mu } | x - z _ { 0 } | g ( x ) , f \in H ^ { m , \alpha , \beta } ( \Omega ) { \mathrm { ~ i f ~ } } s > m - { \frac { 1 } { 2 } } } \end{array}$ . 

Theorem 4.3. For the Jacobi expansion, (1.2) and (1.4) are satisfied. 

Proof. Note that 

$$
\| f - \mathcal {P} _ {N} ^ {f, J a} \| _ {H ^ {m, \alpha , \beta} (\Omega)} ^ {2} = \sum_ {q = 0} ^ {m} \sum_ {n = N + 1} ^ {\infty} [ a _ {n} ^ {(q)} (\alpha + q, \beta + q) ] ^ {2} \sigma_ {n} ^ {\alpha + q, \beta + q},
$$

which directly leads to the desired result by Theorems 3.1–3.2 and 3.5 together with (4.10). □ 

Similar results can be obtained for Gegenbauer and Chebyshev projections from 

$$
\operatorname{span} \left\{T _ {n} (x) \right\} _ {n = 0} ^ {N} = \operatorname{span} \left\{P _ {n} ^ {\left(- \frac {1}{2}, - \frac {1}{2}\right)} (x) \right\} _ {n = 0} ^ {N}
$$

and 

$$
\operatorname{span} \left\{C _ {n} ^ {(\lambda)} (x) \right\} _ {n = 0} ^ {N} = \operatorname{span} \left\{P _ {n} ^ {\left(\lambda - \frac {1}{2}, \lambda - \frac {1}{2}\right)} (x) \right\} _ {n = 0} ^ {N}.
$$

Corollary 4.4. For the Gegenbauer and Chebyshev expansions, it follows for $\mu \geq 0$ that 

(4.11) 

$$
\| f - \mathcal {P} _ {N} ^ {f, G e} \| _ {H ^ {m, \alpha , \beta} (\Omega)} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {m - \lambda - 2 \gamma - 1 / 2} \ln^ {\mu} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {m - s - 1 / 2} \ln^ {\mu} (N)), & f (x) = | x - z _ {0} | ^ {s} \ln^ {\mu} | x - z _ {0} | g (x), \\ \mathcal {O} (N ^ {m - \lambda - 2 \delta - 1 / 2} \ln^ {\mu} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \end{array} \right.
$$

(4.12) 

$$
\| f - \mathcal {P} _ {N} ^ {f, C h} \| _ {H ^ {m, \alpha , \beta} (\Omega)} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {m - 1 / 2 - 2 \gamma} \ln^ {\mu} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {m - s - 1 / 2} \ln^ {\mu} (N)), & f (x) = | x - z _ {0} | ^ {s} \ln^ {\mu} | x - z _ {0} | g (x), \\ \mathcal {O} (N ^ {m - 1 / 2 - 2 \delta} \ln^ {\mu} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \end{array} \right.
$$

where $g \in C ^ { \infty } [ - 1 , 1 ]$ , $z _ { 0 } ~ \in ~ ( - 1 , 1 )$ , min $\{ \lambda + \gamma , \lambda + 2 \gamma , \lambda + \delta , , \lambda + 2 \delta \} > m - \textstyle { \frac { 1 } { 2 } }$ , $\begin{array} { r } { s > m - \frac { 1 } { 2 } } \end{array}$ with $\lambda \geq 0$ . In particular, $i f \gamma , \delta$ are integers and $\mu \geq 1$ , then 

(4.13) 

$$
\| f - \mathcal {P} _ {N} ^ {f, G e} \| _ {H ^ {m, \alpha , \beta} (\Omega)} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {m - \lambda - 2 \gamma - 1 / 2} \ln^ {\mu - 1} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {m - \lambda - 2 \delta - 1 / 2} \ln^ {\mu - 1} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x), \end{array} \right.
$$

$$
\| f - \mathcal {P} _ {N} ^ {f, C h} \| _ {H ^ {m, \alpha , \beta} (\Omega)} = \left\{ \begin{array}{l l} \mathcal {O} (N ^ {m - 1 / 2 - 2 \gamma} \ln^ {\mu - 1} (N)), & f (x) = (1 - x) ^ {\gamma} \ln^ {\mu} (1 - x) g (x), \\ \mathcal {O} (N ^ {m - 1 / 2 - 2 \delta} \ln^ {\mu - 1} (N)), & f (x) = (1 + x) ^ {\delta} \ln^ {\mu} (1 + x) g (x). \end{array} \right.
$$

Numerical results for these estimates on the boundary or an interior point are illustrated in Figure 4.3 and Figure 4.4, which indicates the optimal orders of the estimates. 

Remark 8. In the case $\mu = 0$ and s even, all the three type projections $\parallel f -$ $\mathcal { P } _ { N } ^ { f } \| _ { H ^ { m , \alpha , \beta } ( \Omega ) }$ are exponentially decayed [5, 16, 18, 19, 22, 26, 27]. 

5. Final remarks. Applying the technique of the separation of singularities, the above results can be extended to the general functions with interior and boundary singularities for 

$$
f (x) = g (x) \prod_ {i = 1} ^ {m} | x - x _ {i} | ^ {\gamma_ {i}} \ln^ {\mu} | x - x _ {i} | = \sum_ {i = 1} ^ {m} | x - x _ {i} | ^ {\gamma_ {i}} \ln^ {\mu} | x - x _ {i} | g _ {i} (x), \tag {5.1}
$$

where $- 1 = x _ { m } < x _ { m - 1 } < \dots < x _ { 2 } < x _ { 1 } = 1 , g _ { i } \in C ^ { \infty } [ - 1 , 1 ] , \gamma _ { i } \geq 0$ for $i = 1 , m$ , $\gamma _ { i } > 0$ for $i = 2 , \ldots , m - 1$ and $\mu$ is a positive integer [17, pp. 219–220]. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/37df4ac86fc705ccb2f85e1ce2b7fb85bf2cb4dcf49c9371d4cce0ec262eafb7.jpg)



Fig. 4.3. The weighted norm errors of the truncated Chebyshev, Legendre, and Jacobi expansions kf − $\mathcal { P } _ { N } ^ { f } \Vert _ { H ^ { m , \alpha , \beta } ( \Omega ) }$ with m = 1 for $f ( x ) = | x - { \textstyle \frac { 1 } { 2 } } | ^ { 2 . 7 } \ln ^ { 2 } | x - { \textstyle \frac { 1 } { 2 } } |$ , respectively.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/7d27c324-eca9-4384-a14c-869a676d0714/03c0b0eda203e5c035f57bffc60b505d188358102425cc78d4f27c654f5794e5.jpg)



Fig. 4.4. The weighted norm errors of the truncated Chebyshev, Legendre, and Jacobi expansions $\| f - \mathcal { P } _ { N } ^ { f } \| _ { H ^ { m , \alpha , \beta } ( \Omega ) }$ with m = 1 for $f ( x ) = ( 1 - x ) ^ { 1 . 6 } \ln ^ { 2 } ( 1 - x )$ .


Corollary 5.1. Suppose f(x) is defined by (5.1); then the coefficients in the Jacobi series of f (x) satisfy 

$$
| a _ {n} (\alpha , \beta) | = \left\{ \begin{array}{l} \mathcal {O} \left(\max \left\{n ^ {- \min \left\{1 + \alpha + 2 \gamma_ {0}, 1 + \alpha + 2 \gamma_ {m} \right\}}, n ^ {- \min \left\{\gamma_ {2} + \frac {1}{2}, \ldots , \gamma_ {m - 1} + \frac {1}{2} \right\}} \ln (n) \right\} \ln^ {\mu - 1} (n)\right), \\ \qquad \qquad \qquad \qquad \gamma_ {1}, \gamma_ {m} \in \mathbf {Z} _ {+}, \\ \mathcal {O} \left(\max \left\{n ^ {- 1 - \alpha - 2 \gamma_ {1}}, n ^ {- \min \left\{1 + \alpha + 2 \gamma_ {m}, \gamma_ {2} + \frac {1}{2}, \ldots , \gamma_ {m - 1} + \frac {1}{2} \right\}} \ln (n) \right\} \ln^ {\mu - 1} (n)\right), \\ \qquad \qquad \qquad \qquad \gamma_ {1} \in \mathbf {Z} _ {+},   \gamma_ {m} \notin \mathbf {Z} _ {+}, \\ \mathcal {O} \left(\max \left\{\ln (n) n ^ {- 1 - \beta - 2 \gamma_ {m}}, n ^ {- \min \left\{1 + \alpha + 2 \gamma_ {1}, \gamma_ {2} + \frac {1}{2}, \ldots , \gamma_ {m - 1} + \frac {1}{2} \right\}} \ln (n) \right\} \ln^ {\mu - 1} (n)\right), \\ \qquad \qquad \qquad \qquad \gamma_ {1} \notin \mathbf {Z} _ {+},   \gamma_ {m} \in \mathbf {Z} _ {+}, \\ \mathcal {O} \left(n ^ {- \min \left\{1 + \alpha + 2 \gamma_ {1}, 1 + \beta + 2 \gamma_ {m}, \gamma_ {2} + \frac {1}{2}, \ldots , \gamma_ {m - 1} + \frac {1}{2} \right\}} \ln^ {\mu} (n)\right), \\ \qquad \qquad \qquad \qquad \gamma_ {1}, \gamma_ {m} \notin \mathbf {Z} _ {+}. \end{array} \right.
$$

From (2.8) and (2.10), one may obtain sharp bounds for the Gegenbauer coefficients (2.5) and Chebyshev coefficients (2.9), and optimal convergence rates for suitably differentiable functions in a similar way. 

Acknowledgments. The author is grateful for the referees’ helpful suggestions and insightful comments, which helped improve the manuscript significantly. The author thanks Prof. Li-Lian Wang, Dr. Guidong Liu, Desong Kong, Qingyang Zhang, and Jiangli Liang for many constructive discussion and comments that helped to improve the presentation of this paper. 

# REFERENCES



[1] M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions, National Bureau of Standards, Washington, D.C., 1964. 





[2] I. Babuska and B. Guo, Optimal estimates for lower and upper bounds of approximation errors in the p-version of the finite element method in two dimensions, Numer. Math., 85 (2000), pp. 219–255. 





[3] I. Babuska and B. Guo, Direct and inverse approximation theorems for the p-version of the finite element method in the framework of weighted Besov spaces. I. Approximability of functions in the weighted Besov spaces, SIAM J. Numer. Anal., 39 (2001/02), pp. 1512– 1538. 





[4] I. Babuska and B. Guo, Direct and inverse approximation theorems for the p-version of the finite element method in the framework of weighted Besov spaces. II. Optimal rate of convergence of the p-version finite element solutions, Math. Models Methods Appl. Sci., 12 (2002), pp. 689–719. 





[5] S. Bernstein, Sur l’order de la meilleure approximation des fonctions continues par des polynomes de degr´e donn´e, Memories de l’. Academie Royale de Belgique, 4 (1912), pp. 1–103. 





[6] P. Castillo, B. Cockburn, D. Schotzau, and C. Schwab ¨ , Optimal a priori error estimates for the hp-version of the local discontinuous Galerkin method for convection-diffusion problems, Math. Comp., 71 (2002), pp. 455–478. 





[7] G. Darboux, M´emoire sur l’approximation des fonctions de tr`es grands nombres, J. de. Math´ematiques, 4 (1878), pp. 5–56, 377–416. 





[8] B. Guo, J. Shen, and L. Wang, Optimal spectral-Galerkin methods using generalized Jacobi polynomials, J. Sci. Comput., 27 (2006), pp. 305–322. 





[9] B. Guo and L. Wang, Jacobi approximations in non-uniformly Jacobi-weighted Sobolev spaces, J. Appr. Theory, 128 (2004), pp. 1–41. 





[10] J. Hesthaven, S. Gottlieb, and D. Gottlieb, Spectral Methods for Time-Dependent Problems, Cambridge University Press, Cambridge, 2007. 





[11] W. Liu and L. Wang, Asymptotics of Generalized Gegenbauer functions of fractional degree, J. Appr. Theory, 253 (2020), 105378. 





[12] W. Liu, L. Wang, and H. Li, Optimal error estimates for Chebyshev approximations of functions with limited regularity in fractional Sobolev-type spaces, Math. Comp., 88 (2019), pp. 2857–2895. 





[13] V. Mackevicius ˆ , Integral and Measure: From Rather Simple to Rather Complex, John Wiley & Sons, New York, 2014. 





[14] J. Shen, T. Tang, and L. Wang, Spectral Methods: Algorithms, Analysis and Applications, Springer, Berlin, 2011. 





[15] G. Szego¨, Orthogonal Polynomial, Academic Mathematical Society, Providence, RI, 1939. 





[16] L. N. Trefethen, Approximation Theory and Approximation Practice, SIAM, Philadelphia, 2013. 





[17] P. D. Tuan and D. Elliott, Coefficients in series expansions for certain classes of functions, Math. Comp., 26 (1972), pp. 213–232. 





[18] H. Wang, On the optimal estimates and comparison of Gegenbauer expansion coefficients, SIAM, J. Numer. Anal., 34 (2016), pp. 1557–1580. 





[19] H. Wang and S. Xiang, On the convergence rates of Legendre approximation, Math. Comp., 81 (2012), pp. 861–877. 





[20] G. N. Watson, A Tretise on the Theory of Bessel Functions, Cambridge University Press, Cambridge, 1922. 





[21] S. Xiang, Numerical analysis of a fast integration methods for highly oscillatory functions, BIT Numer. Anal., 47 (2007), pp. 469–482. 





[22] S. Xiang, On error bounds for orthogonal polynomial expansions and Gauss-Type Quadrature, SIAM J. Numer. Anal., 50 (2012), pp. 1240–1263. 





[23] S. Xiang, On van der Corput-type lemmas for Bessel and Airy transforms and applications, J. Comput. Appl. Math., 351 (2019), pp. 179–185. 





[24] S. Xiang and H. Brunner, Efficient methods for Volterra integral equations with highly oscillatory Bessel kernels, BIT Numer. Anal., 53 (2013), pp. 4241–4263. 





[25] S. Xiang and G. Liu, Optimal decay rates on the asymptotics of orthogonal polynomial expansions for functions of limited regularities, Numer. Math., 145 (2020), pp.117–148. 





[26] Z. Xie, L. Wang, and X. Zhao, On exponential convergence of Gegenbauer interpolation and spectral differentiation, Math. Comp., 82 (2013), pp. 1017–1036. 





[27] X. Zhao, L. Wang, and Z. Xie, Sharp error bounds for Jacobi expansions and Gegenbauer-Gauss quadrature of analytic functions, SIAM J. Numer. Anal., 251 (2012), pp. 1443–1469. 

