# Optimal Error Estimates for Gegenbauer Approximations in Fractional Spaces

1 ·  1,2 ·  3,4 ·  1 

Received: 22 August 2024 / Revised: 10 November 2024 / Accepted: 3 December 2024 

© The Author(s), under exclusive licence to Springer Science+Business Media, LLC, part of Springer Nature 2024, corrected publication 2025 

# Abstract

Gegenbauer polynomial approximations play an important role in applied mathematics and computational physics. In this paper, we present sharp bounds for Gegenbauer expansion coefficients of functions belonging to fractional spaces and then derive some new and sharp error bounds for Gegenbauer approximations in weighted $L ^ { \infty }$ and $L ^ { 2 }$ norms. Ample numerical results are provided to demonstrate the sharpness of the estimates. 

Keywords Approximation by Gegenbauer polynomials Optimal estimates Singular functions · Fractional spaces 

Mathematics Subject Classification Primary 41A10 41A25 41A81 

# 1 Introduction

Polynomial approximation theory plays an important role in many areas of scientific computing, such as spectral and spectral-element methods, p- and hp-version of the finite element methods (see, e.g., [4–6, 14] and the references therein). Over the past decade, there has been a growing research interest in approximation theory of orthogonal polynomials, including Jacobi, Gegenbauer and their special cases like Legendre and Chebyshev polynomials (e.g., Wang et al. [34, 40, 42], Xiang et al. [35–38], Chen et al. [10–12], Wang [29–33] and Zhang et al. [41] etc.). One of the most remarkable advantages of orthogonal polynomials is that their approximation power depends solely on the regularity of the underlying function. When designing efficient numerical methods for differential equations, it is known that the convergence behavior of polynomial-based methods depends upon the regularity of the underlying solutions. In general, the presence of corner singularities (Fig. 1) caused by non-smoothness of the boundary will destroy the convergence rate of numerical methods (refer to [28] for more details). Analysis of corner singularities can be found in [3, 16, 20, 21] and many other references. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/dea76820-43c7-4583-82cf-5d601e885b86/273628666d9a2f0877ba8f4a7dfad6169e5bd16d22d5dc41e5d3b442cbc59d54.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/dea76820-43c7-4583-82cf-5d601e885b86/81b136d1d5c22cd93df99da57c9894a187022a4bc8dc127a23d8c1b7c71f306f.jpg)



Fig. 1 Left: Polygonal region - with boundary vertices $V _ { j } , j = 1 , 2 , \ldots , N$ , vertex subdomains $\omega _ { j }$ and interior vertex angles $\alpha _ { j }$ . Right: The neighbourhood $\omega _ { j }$ of the vertex $V _ { j }$


Here, we address the singularity in the neighbourhood $\omega _ { j }$ of a single corner $V _ { j }$ of E, with arms $E _ { j } , E _ { j - 1 }$ and interior angle $\alpha _ { j }$ as shown in Fig. 1, and we define neighbourhoods $\omega _ { j }$ of $V _ { j }$ for each j such that $\bar { \omega } _ { i } \cap \bar { \omega } _ { j } = \varnothing$ for $i \neq j$ . Also, we take local polar coordinates $( r , \theta )$ centered on the corner. It is interesting to consider the equation $\Delta u = 0$ with zero Dirichlet boundary conditions on $E _ { j }$ , and $E _ { j - 1 }$ . Writing $\Delta u = 0$ in these local polar coordinates we have that 

$$
\Delta u (r, \theta) \equiv \frac {\partial^ {2} u}{\partial r ^ {2}} + \frac {1}{r} \frac {\partial u}{\partial r} + \frac {1}{r ^ {2}} \frac {\partial^ {2} u}{\partial \theta^ {2}} = 0 \tag {1.1}
$$

with $u ( r , \theta ) = 0$ on $\theta \ : = \ : 0$ (corresponding to $E _ { j - 1 } ~ )$ and $\theta \ : = \ : \alpha \pi$ (corresponding to $E _ { j } )$ . For each $j ,$ each corner $V _ { j }$ has two adjoining edges, denoted here as $E _ { j - 1 }$ and $E _ { j }$ , with an interior vertex angle $\alpha _ { j }$ formed between these edges. From separation of variables the functions 

$$
r ^ {\beta_ {j}} \sin \beta_ {j} \theta \quad \mathrm{and} \quad r ^ {\beta_ {j}} \cos \beta_ {j} \theta ,
$$

with $\begin{array} { r } { \beta _ { j } = \frac { \pi } { \alpha _ { j } } } \end{array}$ satisfy (1.1) for any α. Now, $r ^ { \beta }$ sin $\beta \theta$ also satisfies the homogeneous Dirichlet boundary conditions on $E _ { j - 1 }$ and $E _ { j }$ , so that $r ^ { \beta _ { j } }$ sin $\beta _ { j } \theta$ is the leading term in the solution u. More specifically, the values of $\beta _ { j }$ for different boundary coefficient combinations are given in Table 1, and refer to [3, p. 44-46] for more details. 


Table 1 Singularity parameters $\beta _ { j }$ for different boundary conditions and interior vertex angles $\alpha _ { j }$ (see [3, p. 46])


<table><tr><td><eq>E_{j-1}</eq></td><td><eq>E_j</eq></td><td><eq>\beta_j</eq></td></tr><tr><td>Dirichlet</td><td>Dirichlet</td><td><eq>\pi/\alpha_j</eq></td></tr><tr><td>Neumann</td><td>Neumann</td><td><eq>\pi/\alpha_j</eq></td></tr><tr><td>Robin</td><td>Robin</td><td><eq>\pi/\alpha_j</eq></td></tr><tr><td>Dirichlet</td><td>Neumann</td><td><eq>\pi/2\alpha_j</eq></td></tr><tr><td>Dirichlet</td><td>Robin</td><td><eq>\pi/2\alpha_j</eq></td></tr><tr><td>Neumann</td><td>Robin</td><td><eq>\pi/\alpha_j</eq></td></tr></table>

Let $\mathcal { P } _ { N } ^ { f }$ denote an orthogonal projection or interpolation operator on the set of all polynomials of degree at most N . Generally, we have approximation results as $\| f - \mathcal { P } _ { N } ^ { f } \| _ { S } \leq$ $C N ^ { - \rho } | f | _ { B }$ , where $\rho \ge 0$ , and C is a positive constant independent of N and $f . S$ is a certain Sobolev space, B is a related Sobolev or Besov space, and $\rho$ depends on the regularity of $f .$ . However, for B is the standard weighted Sobolev space $H _ { \omega } ^ { k } ( \Omega )$ with integer $k \geq 0$ and certain weight function $\omega ( x )$ on $\Omega = ( - 1 , 1 )$ , it only has suboptimal order for functions with interior or endpoint singularities (e.g., [13, 27]). For more details, see Liu et al. [18, 19, 39]. To fill the above gap, Liu, Wang and Li introduced for the first time in [18] certain fractional Sobolev-type spaces and derived optimal Chebyshev polynomial approximation to functions with interior and endpoint singularities within this new framework. 

Recently, Xiang and Liu [38] derived optimal asymptotics on the decay of the coefficients for functions of limited regularity expanded in terms of Jacobi and Gegenbauer polynomial series. Subsequently, based on the Hilb type formula of Jacobi polynomials, Xiang [36] derived optimal decay rates on the Jacobi expansion coefficients by applying van der Corput type lemmas for functions of algebraic and logarithmatic singularities, which leads to the optimal convergence rates on the Jacobi, Gegenbauer, and Chebyshev orthogonal projections. More recently, Wang [31] proved that Chebyshev spectral differentiations converge faster than their best counterparts except in a neighborhood of the singularity, and in the particular case where the singularity is located in the interior of interval, they converge even faster than their best counterparts in the maximum norm. Motivated by comparing the convergence behavior of Gegenbauer projections and best approximations, Wang [33] analyzed optimal rates of convergence of Gegenbauer projections $S _ { n } ^ { \lambda } ( f )$ in the maximum norm. Xiang et al. [37] studied pointwise error estimates and local superconvergence of Jacobi expansions and gave some remarkable results. 

Inspired from Wang [33] and Liu et al. [17–19], we introduce fractional spaces to characterize functions with singularities and derive the optimal (weighted) $L ^ { \infty }$ -estimates and $L ^ { 2 }$ -estimates of the Gegenbauer polynomial approximations for functions in these spaces. We highlight the main contributions of this work as follows: 

(i) We refine the definition of fractional spaces and present sharp bounds for Gegenbauer expansion coefficients of functions in fractional spaces. 

(ii) We show the optimal (weighted) $L ^ { \infty }$ -estimates and $L ^ { 2 } \mathrm { - e s t i m a t e s }$ of the Gegenbauer polynomial approximations for functions with interior or endpoint singularities. 

The paper is organised as follows. In Sect. 2, we review the basics of Gegenbauer polynomials. In Sect. 3, we recall the basics of several spaces of functions and useful formulas that will be used to characterise and analyze the regularity of the class of functions. We derive the main results on Gegenbauer approximation of functions with interior or endpoint singularities in Sect. 4. We conclude this paper in Sect. 5. 

# 2 Preliminaries

In this section, we introduce some notations and review some relevant properties of Gegenbauer (ultraspherical) polynomials and the generalised Gegenbauer functions of fractional degree (GGF-Fs) that are necessary throughout the paper. 

# 2.1 Gegenbauer Polynomials

Let N and R be the set of positive integers and real numbers, respectively. Denote 

$$
\mathbb {N} _ {0} := \{0 \} \cup \mathbb {N}, \mathbb {R} ^ {+} := \bigl \{a \in \mathbb {R}: a > 0 \bigr \}, \mathbb {R} _ {0} ^ {+} := \{0 \} \cup \mathbb {R} ^ {+}.
$$

For $\alpha \in \mathbb { R }$ , the rising factorial in the Pochhammer symbol is related by 

$$
(\alpha) _ {n} = \left\{ \begin{array}{l l} 1, & n = 0, \\ \alpha (\alpha + 1) \dots (\alpha + n - 1), & n \geq 1, \end{array} \right.
$$

and has the property: 

$$
(\alpha) _ {n} = \frac {\Gamma (\alpha + n)}{\Gamma (\alpha)}, \quad \alpha \neq 0, - 1, - 2, \ldots ,
$$

where $\Gamma ( \cdot )$ is the usual Gamma function. Recall the hypergeometric function (cf. [2]) 

$$
{ } _ { 2 } F _ { 1 } \bigg [ \begin{array} { c } a , b \\ c \end{array} ; z \bigg ] = \sum _ { k = 0 } ^ { \infty } \frac { ( a ) _ { k } ( b ) _ { k } } { ( c ) _ { k } } \frac { z ^ { k } } { k ! } , | z | <   1 , a , b , c \in \mathbb { R } , - c \notin \mathbb { N } _ { 0 } .
$$

Definition 2.1 For $n \in  { \mathbb { N } } _ { 0 }$ and $\lambda > - 1 / 2$ , the Gegenbauer polynomial is defined by 

$$
G _ {n} ^ {(\lambda)} (x) := {} _ {2} F _ {1} \Big [ \begin{array}{c} - n, n + 2 \lambda \\ \lambda + \frac {1}{2} \end{array} ; \frac {1 - x}{2} \Big ], \quad x \in [ - 1, 1 ],
$$

which has a normalisation different from that in Szeg ˝o [26], i.e., $G _ { n } ^ { ( \lambda ) } ( 1 ) = 1$ . 

Notice that if $\lambda = 0 , G _ { n } ^ { ( \lambda ) } ( x )$ reduces to the Chebyshev polynomial 

$$
T _ {n} (x) = G _ {n} ^ {(0)} (x) = _ {2} F _ {1} \Big [ \begin{array}{c} - n, n \\ \frac {1}{2} \end{array} ; \frac {1 - x}{2} \Big ] = \cos (n \arccos (x)).
$$

If $\lambda = 1 / 2 , G _ { n } ^ { ( \lambda ) } ( x )$ reduces to the usual Legendre polynomial of degree n, i.e., $L _ { n } ( x ) =$ $G _ { n } ^ { ( 1 / 2 ) } ( x )$ . 

Gegenbauer polynomials satisfy the three-term recurrence relation 

$$
(n + 2 \lambda) G _ {n + 1} ^ {(\lambda)} (x) = 2 (n + \lambda) x G _ {n} ^ {(\lambda)} (x) - n G _ {n - 1} ^ {(\lambda)} (x), n \geq 1,
$$

where $G _ { 0 } ^ { ( \lambda ) } ( x ) = 1$ and $G _ { 1 } ^ { ( \lambda ) } ( x ) = x$ . Moreover, according to [24, Thm. 3.24], Gegenbauer polynomials satisfy the following inequality: 

$$
\left| G _ {n} ^ {(\lambda)} (x) \right| \leq G _ {n} ^ {(\lambda)} (1) = 1, \quad | x | \leq 1, \quad \lambda \geq 0. \tag {2.1}
$$

For a fixed λ, these Gegenbauer polynomials are orthogonal over the interval 1, 1 with respect to the weight function $\omega _ { \lambda } ( x ) : = \left( 1 - x ^ { 2 } \right) ^ { \lambda - \frac 1 2 }$ , namely, 

$$
\int_ {- 1} ^ {1} G _ {n} ^ {(\lambda)} (x) G _ {m} ^ {(\lambda)} (x) \omega_ {\lambda} (x) \mathrm{d} x = \gamma_ {n} ^ {(\lambda)} \delta_ {n m}, \tag {2.2}
$$

where $\delta _ { m n }$ is the Kronecker delta, and 

$$
\gamma_ {n} ^ {(\lambda)} := \frac {2 ^ {2 \lambda - 1} n ! \Gamma^ {2} (\lambda + 1 / 2)}{(n + \lambda) \Gamma (n + 2 \lambda)}, n \geq 1; \gamma_ {0} ^ {(\lambda)} := \frac {\sqrt {\pi} \Gamma (\lambda + 1 / 2)}{\Gamma (\lambda + 1)}.
$$

The Rodrigues formula of Gegenbauer polynomials reads 

$$
\omega_ {\lambda} (x) G _ {n} ^ {\lambda} (x) = - \frac {1}{2 \lambda + 1} \frac {\mathrm{d}}{\mathrm{d} x} \left\{\omega_ {\lambda + 1} (x) G _ {n - 1} ^ {\lambda + 1} (x) \right\}, \tag {2.3}
$$

which will be used in the following analysis of the Gegenbauer coefficients. 

# 2.2 Generalised Gegenbauer Functions of Fractional Degree

Throughout this paper, the GGF-Fs we use are defined by (cf. [18, Def. 2.1]) 

Definition 2.2 For real $\lambda > - 1 / 2$ and real $\mu \geq 0$ , the right GGF-F of degree $\mu$ is defined by 

$$
{ } ^ { r } G _ { \mu } ^ { ( \lambda ) } ( x ) : = _ { 2 } F _ { 1 } \Big [ \begin{array} { c } - \mu , \mu + 2 \lambda \\ \lambda + \frac { 1 } { 2 } \end{array} ; \frac { 1 - x } { 2 } \Big ] ,
$$

for $\begin{array} { r l r } { x } & { { } \in } & { ( - 1 , 1 ) } \end{array}$ ; while the left GGF-F of degree $\mu$ is defined by $^ l G _ { \mu } ^ { ( \lambda ) } ( x ) \quad : =$ $( - 1 ) ^ { \lfloor \mu \rfloor r _ { G _ { \mu } ^ { ( \lambda ) } } } ( - x )$ where $\lfloor \mu \rfloor$ is the largest integer $\leq \mu$ . 

Proposition 2.1 The GGF-Fs defined in Definition 2.2 satisfy 

$$
{ } ^ { r } G _ { n } ^ { ( \lambda ) } ( x ) = { } ^ { l } G _ { n } ^ { ( \lambda ) } ( x ) = G _ { n } ^ { ( \lambda ) } ( x ) , \quad n \in \mathbb { N } _ { 0 } ; \tag {2.4}
$$

$$
{ } ^ { r } G _ { \mu } ^ { ( \lambda ) } ( - x ) = ( - 1 ) ^ { \lfloor \mu \rfloor l } G _ { \mu } ^ { ( \lambda ) } ( x ) , \quad { } ^ { r } G _ { \mu } ^ { ( \lambda ) } ( 1 ) = 1 , \quad { } ^ { l } G _ { \mu } ^ { ( \lambda ) } ( - 1 ) = ( - 1 ) ^ { \lfloor \mu \rfloor } . \tag {2.5}
$$

The uniform upper bound of GGF-Fs is as follows. 

Lemma 2.1 (See [17, (4.16)] and [18, (4.30)]). For real $\lambda \geq 0$ and $\mu \geq 0$ , we have 

$$
\max _ {| x | \leq 1} \left\{\varpi_ {\lambda} (x) | ^ {r} G _ {\mu} ^ {(\lambda)} |, \varpi_ {\lambda} (x) | ^ {l} G _ {\mu} ^ {(\lambda)} | \right\} \leq \frac {\Gamma (\lambda + 1 / 2) \Gamma ((\mu + 1) / 2)}{\sqrt {\pi} \Gamma ((\mu + 1) / 2 + \lambda)}, \tag {2.6}
$$

where 

$$
\varpi_ {\lambda} (x) = \left(1 - x ^ {2}\right) ^ {\lambda / 2}, \quad 0 \leq \lambda <   1; \quad \varpi_ {\lambda} (x) = \left(1 - x ^ {2}\right) ^ {\lambda - \frac {1}{2}} = \omega_ {\lambda} (x), \quad \lambda \geq 1.
$$

# 3 Spaces of Functions

In this section, we recall the basics of several spaces of functions and useful formulas that will be used to characterise and analyze the regularity of the class of functions of interest. Then we give some definitions of the fractional spaces in order to deal with interior or endpoint singularities. 

# 3.1 Important Basics and Formulas

Denote by ${ \boldsymbol \Omega } = ( a , b )$ a finite open interval in R. Let AC(-) be the space of absolutely continuous functions on -, and BV(-) be the space of functions of bounded variation on -. Recall that (cf. [8, Chap. 2] and [7, p. 207]): a real function $f ( x ) \in \mathbf { A C } ( \Omega )$ , if for any $\varepsilon > 0$ , there exists a $\delta > 0$ , such that for every finite sequence of disjoint intervals $( \alpha _ { i } , \beta _ { i } )$ in $( a , b )$ , we have $\begin{array} { r } { \sum _ { i } \left| \beta _ { i } - \alpha _ { i } \right| < \delta } \end{array}$ , which implies $\textstyle \sum _ { i } | f ( \beta _ { i } ) - f ( \alpha _ { i } ) | < \varepsilon$ . Denote by $f ( y + )$ the limit of f for x approaching y from the right. Similarly, the left limit is denoted by $f ( y - )$ . Let $f ( a ) : = f ( a + )$ and $f ( b ) : = f ( b - )$ . According to [7, Chap. 8]: a function $f ( x ) \in \mathbf { B } \mathbf { V } ( \Omega )$ ), if there exists a positive constant M such that $\begin{array} { r } { V ( P , f ) : = \sum _ { i = 1 } ^ { k } | f ( x _ { i } ) - f ( x _ { i - 1 } ) | \leq M } \end{array}$ for each partition $P = \{ a = x _ { 0 } < x _ { 1 } < \cdot \cdot \cdot < x _ { k } = b \}$ . Hence the total variation of f on - as $V _ { \Omega } [ f ] : = \operatorname* { s u p } \{ V ( P , f ) \}$ , where the supreme is taken over all partitions. In general, $V _ { \Omega } [ f ]$ could be infinite. According to [25, Prop. 4.2], if $f \in \mathrm { A C } ( \Omega )$ , then $\begin{array} { r } { V _ { \Omega } [ f ] = \int _ { a } ^ { b } \left| f ^ { \prime } ( x ) \right| } \end{array}$ dx . Here, we review the definition of the Riemann-Stieltjes (RS) integral (cf. [15, Chap. X]). Given numbers $c _ { i }$ with $x _ { i - 1 } \leq c _ { i } \leq x _ { i }$ of the partition P, the Riemann-Stieltjes sum as: $\begin{array} { r } { S _ { c } ( P , f , g ) : = \sum _ { i = 1 } ^ { k } f ( c _ { i } ) \big ( g ( x _ { i } ) - g ( x _ { i - 1 } ) \big ) } \end{array}$ . Denote by $\delta ( P ) = \operatorname* { m a x } _ { 1 \leq i \leq k } \{ x _ { i } - x _ { i - 1 } \}$ 1 i k the maximum length of the partition P, if $\operatorname * { l i m } _ { \delta ( P )  0 } S _ { c } ( P , f , g )$ exists, denote by $\textstyle \int _ { a } ^ { b } f \mathrm { d } g$ , then f is ${ \mathrm { R S } } ( g )$ -integrable. The definition of the improper Riemann-Stieltjes integral is similar to the one discussed in Chapter 10 of [9], we omit it. Denote 

$$
\int_ {a} ^ {b} f | \mathrm{d} g | := \lim _ {\delta (P) \to 0} \sum_ {i = 1} ^ {k} f (c _ {i}) \bigl | g (x _ {i}) - g (x _ {i - 1}) \bigr |,
$$

we have $\begin{array} { r } { \int _ { a } ^ { b } | \mathrm { d } g | = V _ { \Omega } [ g ] } \end{array}$ 

Also, a function f is RS(g)-integrable, if $\textstyle \int _ { a } ^ { b } f \mathrm { d } g \ < \infty$ for $g ~ \in \mathrm { { \bf ~ B V } } ( \Omega )$ . According to [15, Prop. 1.3 on p. 282], we have the following important property: for all $g \in { \mathrm { B V } } ( \Omega )$ and f is RS(g)-integrable, then 

$$
\left| \int_ {a} ^ {b} f \mathrm{d} g \right| \leq \| f \| _ {\infty} V _ {\Omega} [ g ], \tag {3.1}
$$

where $\| f \| _ { \infty }$ is the $L ^ { \infty }$ -norm of f on -. Also, if f is continuous and $g \in \mathsf { A C } ( \Omega )$ , then $f \in \operatorname { R S } ( g )$ )-integrable, and (cf. [15, p. 284]) 

$$
\int_ {a} ^ {b} f (x) \mathrm{d} g (x) = \int_ {a} ^ {b} f (x) g ^ {\prime} (x) \mathrm{d} x. \tag {3.2}
$$

Assume that f is continuous and g of bounded variation on $( a , b )$ , then $f \in \operatorname { R S } ( g )$ (cf. [15, p. 283]). 

Denote the ordinary derivative by $\mathcal { D } ^ { k } = \mathrm { d } ^ { k } / \mathrm { d } x ^ { k }$ (with $k \in \mathbb { N } )$ , and let $f ^ { ( k ) } = { \mathcal { D } } ^ { k } f$ . We first recall the definition of the Riemann-Liouville fractional integrals (cf. [23, p. 33, p. 44]) as follows: For any $f \in L ^ { 1 } ( \Omega )$ , the left-sided and right-sided Riemann-Liouville fractional integrals of order $s \in \mathbb { R } ^ { + }$ are defined by 

$$
(I _ {a +} ^ {s} f) (x) = \frac {1}{\Gamma (s)} \int_ {a} ^ {x} \frac {f (y)}{(x - y) ^ {1 - s}} \mathrm{d} y, \quad (I _ {b -} ^ {s} f) (x) = \frac {1}{\Gamma (s)} \int_ {x} ^ {b} \frac {f (y)}{(y - x) ^ {1 - s}} \mathrm{d} y, \quad x \in \Omega . \tag {3.3}
$$

According to [23, p. 51], we have the following property: for $f ( x ) \in L ^ { 1 } ( \Omega )$ , then 

$$
\lim _ {s \to 0} (I _ {a +} ^ {s} f) (x) = f (x),
$$

for any Lebesgue point of a function $f ( x )$ and consequently almost everywhere on $\Omega$ . 

For $\sigma \in ( k - 1 , k ]$ with $k \in \mathbb N .$ , the left-sided and right-side Caputo fractional derivatives of order $\sigma$ are respectively defined by 

$$
\left(^ {C} D _ {a +} ^ {\sigma} f\right) (x) = \left(I _ {a +} ^ {k - \sigma} f ^ {(k)}\right) (x), \quad \left(^ {C} D _ {b -} ^ {\sigma} f\right) (x) = (- 1) ^ {k} \left(I _ {b -} ^ {k - \sigma} f ^ {(k)}\right) (x), \quad x \in \Omega , \tag {3.4}
$$

and the left-sided and right-side Riemann-Liouville fractional derivatives of order $\sigma$ are respectively defined by 

$$
(^ {R} D _ {a +} ^ {\sigma} f) (x) = \mathcal {D} ^ {k} \bigl \{I _ {a +} ^ {k - \sigma} f \bigr \} (x), \quad (^ {R} D _ {b -} ^ {\sigma} f) (x) = (- 1) ^ {k} \mathcal {D} ^ {k} \bigl \{I _ {b -} ^ {k - \sigma} f \bigr \} (x), \quad x \in \Omega . \tag {3.5}
$$

Recall the formula of fractional integration by parts (cf. [19]). 

Lemma 3.1 Let $\rho \geq 0 , f ( x ) \in L ^ { 1 } ( \Omega ) a n d g ( x ) \in \mathrm { A C } ( \Omega )$ 

(i) $I f I _ { b - } ^ { \rho } f ( x ) \in \mathbf { B } \mathbf { V } ( \Omega )$ , then 

$$
\int_ {a} ^ {b} f (x) I _ {a +} ^ {\rho} g ^ {\prime} (x) \mathrm{d} x = \left. \left\{g (x) I _ {b -} ^ {\rho} f (x) \right\} \right| _ {a +} ^ {b -} - \int_ {a} ^ {b} g (x) \mathrm{d} \left\{I _ {b -} ^ {\rho} f (x) \right\}. \tag {3.6}
$$

(ii) $I f I _ { a + } ^ { \rho } f ( x ) \in \mathbf { B } \mathbf { V } ( \Omega )$ , then 

$$
\int_ {a} ^ {b} f (x) I _ {b -} ^ {\rho} g ^ {\prime} (x) \mathrm{d} x = \left. \left\{g (x) I _ {a +} ^ {\rho} f (x) \right\} \right| _ {a +} ^ {b -} - \int_ {a} ^ {b} g (x) \mathrm{d} \left\{I _ {a +} ^ {\rho} f (x) \right\}. \tag {3.7}
$$

# 3.2 Fractional Spaces

Definition 3.1 We say u is of $\mathbf { A C } ^ { k } ( \Omega )$ if $u , u ^ { \prime } , \cdot \cdot \cdot , u ^ { ( k - 1 ) } \in \mathrm { A C } ( \Omega )$ with integer $k \geq 1$ Also, we denote $\mathrm { A C } ^ { 0 } ( \Omega ) : = L ^ { 1 } ( \Omega )$ . Equipped with the norm: 

$$
\| u \| _ {\mathrm{AC} ^ {k} (\Omega)} := \sum_ {i = 0} ^ {k} \| u ^ {(i)} \| _ {L ^ {1} (\Omega)}.
$$

Definition 3.2 We say u is of $\mathbb { W } ^ { k } ( \Omega )$ if $u \in \mathsf { A C } ^ { k - 1 } ( \Omega )$ and $u ^ { ( k - 1 ) } \in \mathrm { B V } ( \Omega )$ with integral $k \geq 1$ . Also, we denote $\mathbb { W } ^ { 0 } ( \Omega )$ is a Radon measure. Equipped with the norm: 

$$
\| u \| _ {\mathbb {W} ^ {k} (\Omega)} := \| u \| _ {\mathrm{AC} ^ {k - 1} (\Omega)} + V _ {\Omega} [ u ^ {(k - 1)} ], \quad k \geq 1.
$$

Definition 3.3 For a fixed $\theta \in \Omega = ( a , b )$ , we denote $\Omega _ { \theta } ^ { - } : = ( a , \theta )$ and $\Omega _ { \theta } ^ { + } : = ( \theta , b )$ . For $\sigma \in ( k - 1 , k )$ ) with integer $k \geq 2$ , we define the fractional space: 

$$
\mathbb {W} _ {\theta} ^ {\sigma} (\Omega) := \left\{u \in \mathrm{AC} ^ {k - 1} (\Omega), ^ {C} D _ {\theta -} ^ {\sigma - 1} u \in \mathrm{BV} (\Omega_ {\theta} ^ {-}) \text { and } ^ {C} D _ {\theta +} ^ {\sigma - 1} u \in \mathrm{BV} (\Omega_ {\theta} ^ {+}) \right\}, \tag {3.8}
$$

equipped with the norm: 

$$
\| u \| _ {\mathbb {W} _ {\theta} ^ {\sigma} (\Omega)} = \sum_ {i = 0} ^ {k} \| u ^ {(i)} \| _ {L ^ {1} (\Omega)} + U _ {\theta} ^ {\sigma}, \tag {3.9}
$$

where the semi-norm is defined by 

$$
U _ {\theta} ^ {\sigma} := V _ {\Omega_ {\theta} ^ {-}} \left[ ^ {C} D _ {\theta -} ^ {\sigma - 1} u \right] + V _ {\Omega_ {\theta} ^ {+}} \left[ ^ {C} D _ {\theta +} ^ {\sigma - 1} u \right] + \left| ^ {C} D _ {\theta -} ^ {\sigma - 1} u (\theta -) \right| + \left| ^ {C} D _ {\theta +} ^ {\sigma - 1} u (\theta +) \right|. \tag {3.10}
$$

For $\theta  a + ( \mathrm { o r } \theta  b - )$ the $\Omega _ { \theta } ^ { - }$ (or $\Omega _ { \theta } ^ { + } )$ and corresponding norm vanish. For $\sigma = k$ , the space $\mathbb { W } _ { \theta } ^ { \sigma } ( \Omega )$ defined by the same as $\mathbb { W } ^ { \sigma } ( \Omega )$ , if $\sigma \geq 1$ , denote $\| u \| _ { \mathbb { W } _ { \theta } ^ { \sigma } ( \Omega ) } : = \| u \| _ { \mathbb { W } ^ { \sigma } ( \Omega ) }$ and 

$U _ { \theta } ^ { \sigma } : = V _ { \Omega } [ u ^ { ( \sigma - 1 ) } ]$ . For $\sigma \in ( 0 , 1 )$ , the $^ C D _ { \theta - } ^ { \sigma - 1 }$ and $^ C D _ { \theta + } ^ { \sigma - 1 }$ are replace by $I _ { \theta - } ^ { 1 - \sigma }$ and $I _ { \theta + } ^ { 1 - \sigma }$ , respectively. 

Definition 3.4 For some m, $k \in  { \mathbb { N } } _ { 0 } , s \in ( 0 , 1 )$ and $\sigma : = m + s .$ , assume that $u \in \mathbb { W } _ { a + } ^ { \sigma } ( \Omega )$ (respectively, $u \in \mathbb { W } _ { b - } ^ { \sigma } ( \Omega ) ) , v \in \mathbb { W } ^ { k + 1 } ( \Omega )$ with $v ( x ) = I _ { a + } ^ { 1 - s } u ^ { ( m ) } ( x )$ (respectively, $v ( x ) =$ $I _ { b - } ^ { 1 - s } u ^ { ( m ) } ( x ) ;$ ). Also, we say u is of $\mathbb { W } _ { a + } ^ { \sigma , k } ( \Omega )$ (respectively, $\mathbb { W } _ { b - } ^ { \sigma , k } ( \Omega ) )$ ). Accordingly, we denote 

$$
U _ {a +} ^ {\sigma , k} := V _ {\Omega} \left[ v ^ {(k)} \right] + \sum_ {i = 0} ^ {k} \left| v ^ {(i)} (a +) \right|,
$$

and 

$$
U _ {b -} ^ {\sigma , k} := V _ {\Omega} \left[ v ^ {(k)} \right] + \sum_ {i = 0} ^ {k} \left| v ^ {(i)} (b -) \right|,
$$

are semi-norms of $\mathbb { W } _ { a + } ^ { \sigma , k } ( \Omega )$ ) and $\mathbb { W } _ { b - } ^ { \sigma , k } ( \Omega )$ , respectively. 

# 4 Gegenbauer Approximations of Functions in Fractional Spaces

In this section, we derive the main results on Gegenbauer approximations of functions with interior or endpoint singularities. 

# 4.1 Exact Formulas and Decay Rate of Gegenbauer Expansion Coefficients

In what follows, we focus on the estimation of the errors between u and its finite Gegenbauer expansion $\pi _ { N } ^ { ( \lambda ) }$ as follows 

$$
u (x) = \sum_ {n = 0} ^ {\infty} \hat {u} _ {n} ^ {(\lambda)} G _ {n} ^ {(\lambda)} (x), \quad \pi_ {N} ^ {(\lambda)} u (x) = \sum_ {n = 0} ^ {N} \hat {u} _ {n} ^ {(\lambda)} G _ {n} ^ {(\lambda)} (x), \tag {4.1}
$$

where the Gegenbauer coefficients are defined by 

$$
\hat {u} _ {n} ^ {(\lambda)} := \frac {1}{\gamma_ {n} ^ {(\lambda)}} \int_ {- 1} ^ {1} u (x) G _ {n} ^ {(\lambda)} (x) \omega_ {\lambda} (x) \mathrm{d} x. \tag {4.2}
$$

In this section, we are concerned with the Gegenbauer approximation of functions with interior singularities. For example, consider $u ( x ) = | x - \theta | ^ { \alpha } g ( x )$ , where $\alpha > - 1$ is a noninteger and $g ( x )$ being a smooth function. We remark that the results can be straightforwardly extended to multiple interior singularities. 

Theorem 4.1 Let $\Lambda = ( - 1 , 1 )$ in place of - in Definition 3.3. Assume that $u \in \mathbb { W } _ { \theta } ^ { \sigma } ( \Lambda )$ with $\sigma > 0$ and $\theta \in \Lambda$ , then for $n \geq \sigma$ , we have 

$$
\begin{array}{l} \hat {u} _ {n} ^ {(\lambda)} = \frac {2 ^ {1 - 2 \lambda - \sigma} (n + \lambda) \Gamma (n + 2 \lambda)}{n ! \Gamma (\lambda + 1 / 2) \Gamma (\lambda + \sigma + 1 / 2)} \left\{^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (\theta) v _ {+} (\theta +) + \int_ {\theta} ^ {1} ^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (x) d v _ {+} (x) \right. \tag {4.3} \\ \left. - ^ {l} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (\theta) v _ {-} (\theta -) + \int_ {- 1} ^ {\theta} ^ {l} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (x) \mathrm{d} v _ {-} (x) \right\}, \\ \end{array}
$$

where 

$$
{ } ^ { z } \mathcal { G } _ { \mu } ^ { ( \gamma , \kappa ) } ( x ) : = \omega _ { \gamma + \kappa } ( x ) { } ^ { z } G _ { \mu - \kappa } ^ { ( \gamma + \kappa ) } ( x ) , \quad z = r , l ; \quad v _ { \pm } ( x ) = { } ^ { C } D _ { \theta \pm } ^ { \sigma - 1 } u ( x ) .
$$

Moreover, we have the following bound: 

$$
\left| \hat {u} _ {n} ^ {(\lambda)} \right| \leq \frac {2 ^ {1 - 2 \lambda - \sigma}}{\sqrt {\pi} \Gamma (\lambda + 1 / 2)} \frac {(n + \lambda) \Gamma (n + 2 \lambda) \Gamma ((n - \sigma + 1) / 2)}{n ! \Gamma ((n + \sigma + 1) / 2 + \lambda)} \mathbb {U} _ {\theta} ^ {\sigma}, \tag {4.4}
$$

where 

(i) For $0 \leq \lambda + \sigma < 1$ , 

$$
\begin{array}{l} \mathbb {U} _ {\theta} ^ {\sigma} := \left| \left\{\omega_ {(\lambda + \sigma) / 2} v _ {-} \right\} (\theta -) \right| + \left| \left\{\omega_ {(\lambda + \sigma) / 2} v _ {+} \right\} (\theta +) \right| + \int_ {- 1} ^ {\theta} \omega_ {(\lambda + \sigma) / 2} | \mathrm{d} v _ {-} (x) | \\ + \int_ {\theta} ^ {1} \omega_ {(\lambda + \sigma) / 2} | \mathrm{d} v _ {+} (x) |. \\ \end{array}
$$

(ii) $F o r \lambda + \sigma \ge 1 , \mathbb { U } _ { \theta } ^ { \sigma } : = U _ { \theta } ^ { \sigma }$ . 

Proof We take two steps to carry out the proof. 

and (i) We split $u \in \mathrm { A C } ^ { m } ( \Omega )$ $\sigma = m + s$ = +    ∈ 0   ∈  . In view of (2.4), we simply denote with $m \in  { \mathbb { N } } _ { 0 }$ , and $s \in ( 0 , 1 ]$ $\mathcal { G } _ { n } ^ { ( \lambda , k ) } ( x ) = { } ^ { l } \mathcal { G } _ { n } ^ { ( \lambda , k ) } ( x ) = { } ^ { r } \mathcal { G } _ { n } ^ { ( \lambda , k ) } ( x )$ . In the first step, we consider $m \geq 1$ , if $\mu  n$ and $\kappa  k$ are integers. According to (2.3), we have 

$$
\mathcal {G} _ {n} ^ {(\lambda , k)} (x) = - \frac {1}{2 (k + \lambda) + 1} \left(\mathcal {G} _ {n} ^ {(\lambda , k + 1)} (x)\right) ^ {\prime}, \quad n \geq k + 1. \tag {4.5}
$$

Using (4.2) and (4.5) with $k = 0 , 1 , \cdots , m - 1$ in order, we obtain from integration by parts that for $n \geq m$ , 

$$
\begin{array}{l} \gamma_ {n} ^ {(\lambda)} \hat {u} _ {n} ^ {(\lambda)} = \int_ {- 1} ^ {1} u (x) \mathcal {G} _ {n} ^ {(\lambda , 0)} (x) d x = - \frac {1}{2 \lambda + 1} \int_ {- 1} ^ {1} u (x) \left(\mathcal {G} _ {n} ^ {(\lambda , 1)} (x)\right) ^ {\prime} d x \\ = \frac {1}{2 \lambda + 1} \int_ {- 1} ^ {1} u ^ {\prime} (x) \mathcal {G} _ {n} ^ {(\lambda , 1)} (x) d x = - \frac {1}{2 \lambda + 1} \frac {1}{2 \lambda + 3} \int_ {- 1} ^ {1} u ^ {\prime} (x) \left(\mathcal {G} _ {n} ^ {(\lambda , 2)} (x)\right) ^ {\prime} d x \\ = \frac {1}{2 \lambda + 1} \frac {1}{2 \lambda + 3} \int_ {- 1} ^ {1} u ^ {\prime \prime} (x) \mathcal {G} _ {n} ^ {(\lambda , 2)} (x) d x = \dots \\ = \frac {1}{(2 \lambda + 1) (2 \lambda + 3) \cdots (2 \lambda + 2 m - 1)} \int_ {- 1} ^ {1} u ^ {(m)} (x) \mathcal {G} _ {n} ^ {(\lambda , m)} (x) d x \\ = \frac {\Gamma (\lambda + 1 / 2)}{2 ^ {m} \Gamma (\lambda + m + 1 / 2)} \int_ {- 1} ^ {1} u ^ {(m)} (x) \mathcal {G} _ {n} ^ {(\lambda , m)} (x) \mathrm{d} x, \tag {4.6} \\ \end{array}
$$

which yields 

$$
\hat {u} _ {n} ^ {(\lambda)} = \frac {1}{\gamma_ {n} ^ {(\lambda)}} \frac {\Gamma (\lambda + 1 / 2)}{2 ^ {m} \Gamma (\lambda + m + 1 / 2)} \int_ {- 1} ^ {1} u ^ {(m)} (x) \mathcal {G} _ {n} ^ {(\lambda , m)} (x) \mathrm{d} x, \quad n \geq m \geq 1. \tag {4.7}
$$

For $n \geq m = 0$ , the (4.7) still holds. 

(ii) In the second step, we further assume that $I _ { \theta - } ^ { 1 - s } u ^ { ( m ) } \in \mathrm { B V } \left( \Omega _ { \theta } ^ { - } \right)$ and $I _ { \theta + } ^ { 1 - s } u ^ { ( m ) } \in$ $\mathrm { B V } \left( \Omega _ { \theta } ^ { + } \right)$ . Then we can conduct the fractional integration by parts using Lemma 3.1. For this purpose, we resort to the identities: for $\sigma > 0$ , and $n \geq \sigma$ , 

$$
\begin{array}{l} \mathcal {G} _ {n} ^ {(\lambda , m)} (x) = - \frac {\Gamma (\lambda + m + 1 / 2)}{2 ^ {s} \Gamma (\lambda + \sigma + 1 / 2)} I _ {1 -} ^ {1 - s} \left\{\binom {r} {\mathcal {G} _ {n} ^ {(\lambda , \sigma)} (x)} ^ {\prime} \right\} \\ = - \frac {\Gamma (\lambda + m + 1 / 2)}{2 ^ {s} \Gamma (\lambda + \sigma + 1 / 2)} I _ {- 1 +} ^ {1 - s} \left\{\binom {l} {\mathcal {G} _ {n} ^ {(\lambda , \sigma)}} (x) ^ {\prime} \right\}, \tag {4.8} \\ \end{array}
$$

which reduces to (4.5) when $s = 1$ . 

Recall [18, Thm. 3.1]: For real $\mu \geq \nu > 0$ and real $\gamma > - 1 / 2$ , the GGF-Fs on ( 1, 1) satisfy the Riemann-Liouville fractional integral formulas: 

$$
I _ {1 -} ^ {\nu} \left\{\omega_ {\gamma} (x) ^ {r} G _ {\mu} ^ {(\gamma)} (x) \right\} = h _ {\gamma} ^ {(- \nu)} \omega_ {\gamma + \nu} (x) ^ {r} G _ {\mu - \nu} ^ {(\gamma + \nu)} (x), \tag {4.9a}
$$

$$
I _ {- 1 +} ^ {\nu} \left\{\omega_ {\gamma} (x) ^ {l} G _ {\mu} ^ {(\gamma)} (x) \right\} = (- 1) ^ {\lfloor \mu \rfloor + \lfloor \mu - \nu \rfloor} h _ {\gamma} ^ {(- \nu)} \omega_ {\gamma + \nu} (x) ^ {l} G _ {\mu - \nu} ^ {(\gamma + \nu)} (x). \tag {4.9b}
$$

For real $\gamma > \nu - 1 / 2$ and real $\mu \geq 0$ , the GGF-Fs on (−1, 1) satisfy the Riemann-Liouville fractional derivative formulas: 

$$
{ } ^ { R } D _ { 1 - } ^ { \nu } \left\{ \omega _ { \gamma } ( x ) ^ { r } G _ { \mu } ^ { ( \gamma ) } ( x ) \right\} = h _ { \gamma } ^ { ( \nu ) } \omega _ { \gamma - \nu } ( x ) ^ { r } G _ { \mu + \nu } ^ { ( \gamma - \nu ) } ( x ) , \tag {4.10a}
$$

$$
{ } ^ { R } D _ { - 1 + } ^ { \nu } \left\{ \omega _ { \gamma } ( x ) ^ { l } G _ { \mu } ^ { ( \gamma ) } ( x ) \right\} = ( - 1 ) ^ { \lfloor \mu \rfloor + \lfloor \mu + \nu \rfloor } h _ { \gamma } ^ { ( \nu ) } \omega _ { \gamma - \nu } ( x ) ^ { l } G _ { \mu + \nu } ^ { ( \gamma - \nu ) } ( x ) , \tag {4.10b}
$$

where 

$$
\omega_ {\alpha} (x) = \left(1 - x ^ {2}\right) ^ {\alpha - \frac {1}{2}}, h _ {\gamma} ^ {(\beta)} = \frac {2 ^ {\beta} \Gamma (\gamma + 1 / 2)}{\Gamma (\gamma - \beta + 1 / 2)}.
$$

Taking $\nu = 1 - s , \gamma = \lambda + \sigma - 1$ and $\mu = n - \sigma + 1$ in (4.9a)-(4.9b) leads to 

$$
{ } ^ { r } \mathcal { G } _ { n } ^ { ( \lambda , m ) } ( x ) = \frac { 2 ^ { 1 - s } \Gamma ( \lambda + m + 1 / 2 ) } { \Gamma ( \lambda + \sigma - 1 / 2 ) } I _ { 1 - } ^ { 1 - s } \left\{ { } ^ { r } \mathcal { G } _ { n } ^ { ( \lambda , \sigma - 1 ) } ( x ) \right\} , \tag {4.11a}
$$

$$
{ } ^ { l } \mathcal { G } _ { n } ^ { ( \lambda , m ) } ( x ) = \frac { 2 ^ { 1 - s } \Gamma ( \lambda + m + 1 / 2 ) } { \Gamma ( \lambda + \sigma - 1 / 2 ) } I _ { - 1 + } ^ { 1 - s } \left\{ { } ^ { l } \mathcal { G } _ { n } ^ { ( \lambda , \sigma - 1 ) } ( x ) \right\} , \tag {4.11b}
$$

Taking $\nu = 1 - s , \gamma = \sigma + s + \lambda - 1$ and $\mu = n - \sigma - s + 1$ in (4.9a)-(4.9b), and then taking the derivative on both sides. By (3.5) and (4.10a)-(4.10b), we deduce that for $n \geq \sigma .$ , 

$$
{ } ^ { r } \mathcal { G } _ { n } ^ { ( \lambda , \sigma - 1 ) } ( x ) = - \frac { 1 } { 2 \lambda + 2 \sigma - 1 } \left\{ { } ^ { r } \mathcal { G } _ { n } ^ { ( \lambda , \sigma ) } ( x ) \right\} ^ { \prime } ,
$$

$$
{ } ^ { l } \mathcal { G } _ { n } ^ { ( \lambda , \sigma - 1 ) } ( x ) = - \frac { 1 } { 2 \lambda + 2 \sigma - 1 } \left\{ { } ^ { l } \mathcal { G } _ { n } ^ { ( \lambda , \sigma ) } ( x ) \right\} ^ { \prime } . \tag {4.12}
$$

Then substituting (4.12) into (4.11a)-(4.11b), we obtain (4.8) immediately. 

Using (4.8), we can rewrite (4.7) as 

$$
\begin{array}{l} \hat {u} _ {n} ^ {(\lambda)} = \frac {2 ^ {1 - 2 \lambda - m} (n + \lambda) \Gamma (n + 2 \lambda)}{n ! \Gamma (\lambda + 1 / 2) \Gamma (\lambda + m + 1 / 2)} \left\{\int_ {- 1} ^ {\theta} u ^ {(m)} (x) \mathcal {G} _ {n} ^ {(\lambda , m)} (x) \mathrm{d} x + \int_ {\theta} ^ {1} u ^ {(m)} (x) \mathcal {G} _ {n} ^ {(\lambda , m)} (x) \mathrm{d} x \right\} \\ = \frac {2 ^ {1 - 2 \lambda - \sigma} (n + \lambda) \Gamma (n + 2 \lambda)}{n ! \Gamma (\lambda + 1 / 2) \Gamma (\lambda + \sigma + 1 / 2)} \left\{\int_ {- 1} ^ {\theta} f (x) I _ {- 1 +} ^ {1 - s} g ^ {\prime} (x) d x + \int_ {\theta} ^ {1} f (x) I _ {1 -} ^ {1 - s} h ^ {\prime} (x) d x \right\}, \tag {4.13} \\ \end{array}
$$

where we denoted 

$$
f (x) = u ^ {(m)} (x), \quad g (x) = - ^ {l} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (x), \quad h (x) = - ^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (x).
$$

Observe from (2.5) and (4.12) that $g ^ { \prime } ( x ) \left( \mathrm { r e s p . ~ } h ^ { \prime } ( x ) \right)$ is continuous on $\left( - 1 , \theta \right] \left( \mathrm { r e s p . } \left[ \theta , 1 \right) \right)$ , and they are also integrable when $\sigma > 0$ . Given that $I _ { \theta - } ^ { 1 - s } f ( x ) \in \mathbf { B } \mathbf { V } \left( \Omega _ { \theta } ^ { - } \right)$ and $I _ { \theta + } ^ { 1 - s } f ( x ) \in$ $\mathrm { B V } \left( \Omega _ { \theta } ^ { + } \right)$ , we can apply fractional integration by parts to the two integrals in (4.13) using Lemma 3.1. Then, we obtain the identity (4.3) immediately. 

Using the identity (4.3), for $0 \leq \lambda + \sigma < 1$ and $\lambda + \sigma \geq 1$ , we obtain (4.4) from the bound (2.6) (with $\lambda \to \sigma + \lambda , \nu \to n - \sigma )$ and (3.1). This completes the proof. □ 

# 4.2 Approximation Results for Functions with Interior Singularities

With Theorem 4.1 at our disposal, we are now ready to estimate the Gegenbauer expansion errors in the weighted $L ^ { \infty }$ -norm and $L ^ { 2 }$ -norm for functions with interior singularities. 

Theorem 4.2 Assume that the conditions in Theorem 4.1 hold, we have that for $\lambda \geq 0$ and $1 < \sigma \le N + 1$ , 

$$
\left\| \varpi_ {\lambda} (x) \left(u (x) - \pi_ {N} ^ {(\lambda)} u (x)\right) \right\| _ {L ^ {\infty}} \leq \frac {U _ {\theta} ^ {\sigma}}{2 ^ {\sigma - 2} (\sigma - 1) \pi} \frac {\Gamma ((N - \sigma) / 2 + 1)}{\Gamma ((N + \sigma) / 2)}. \tag {4.14}
$$

Proof From Lemma 2.1, (4.1) and (4.4), we obtain 

$$
\begin{array}{l} \left\| \varpi_ {\lambda} (x) \left(u (x) - \pi_ {N} ^ {(\lambda)} u (x)\right) \right\| _ {L ^ {\infty}} = \left\| \sum_ {n = N + 1} ^ {\infty} \hat {u} _ {n} ^ {(\lambda)} \varpi_ {\lambda} (x) G _ {n} ^ {(\lambda)} (x) \right\| _ {L ^ {\infty}} \\ \leq \sum_ {n = N + 1} ^ {\infty} \left| \hat {u} _ {n} ^ {(\lambda)} \right| \frac {\Gamma (\lambda + 1 / 2) \Gamma ((n + 1) / 2)}{\sqrt {\pi} \Gamma ((n + 1) / 2 + \lambda)} \tag {4.15} \\ \leq \frac {U _ {\theta} ^ {\sigma}}{2 ^ {\sigma - 1} \pi} \sum_ {n = N + 1} ^ {\infty} \widehat {\Upsilon} _ {n} ^ {(\lambda , \sigma)} \frac {\Gamma ((n - \sigma + 1) / 2)}{\Gamma ((n + \sigma + 1) / 2)}, \\ \end{array}
$$

where 

$$
\widehat {\Upsilon} _ {n} ^ {(\lambda , \sigma)} = \frac {1}{2 ^ {2 \lambda}} \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)} \frac {\Gamma ((n + 1) / 2)}{\Gamma ((n + 1) / 2 + \lambda)} \frac {\Gamma ((n + \sigma + 1) / 2)}{\Gamma ((n + \sigma + 1) / 2 + \lambda)}.
$$

Thus, by $z \Gamma ( z ) = \Gamma ( z + 1 )$ and a simple calculation leads to 

$$
\frac {\widehat {\Upsilon} _ {n} ^ {(\lambda , \sigma)}}{\widehat {\Upsilon} _ {n + 2} ^ {(\lambda , \sigma)}} = \frac {(n + 2) (n + \lambda) (n + \sigma + 1 + 2 \lambda)}{(n + 2 \lambda) (n + 2 + \lambda) (n + \sigma + 1)} = 1 - \frac {2 \lambda ((n + \lambda) (\sigma - 1) + n + \sigma + 1)}{(n + 2 \lambda) (n + \lambda + 2) (n + \sigma + 1)} <   1, \tag {4.16}
$$

so the sequence $\big \{ \widehat { \Upsilon } _ { n } ^ { ( \lambda , \sigma ) } \big \} _ { n \geq N + 1 }$ is increasing for $\lambda \geq 0$ . Recall that (cf. [22, (5.11.13)]): for $a < b ,$ , 

$$
\frac {\Gamma (z + a)}{\Gamma (z + b)} = z ^ {a - b} + \frac {1}{2} (a - b) (a + b - 1) z ^ {a - b - 1} + O (z ^ {a - b - 2}), \quad z \gg 1. \tag {4.17}
$$

We can obtain 

$$
\lim _ {n \rightarrow \infty} \widehat {\Upsilon} _ {n} ^ {(\lambda , \sigma)} = 1, \quad \widehat {\Upsilon} _ {n} ^ {(\lambda , \sigma)} <   \widehat {\Upsilon} _ {n + 2} ^ {(\lambda , \sigma)} <   \dots \leq 1. \tag {4.18}
$$

Using (4.18), we can rewrite (4.15) as 

$$
\left\| \varpi_ {\lambda} (x) \left(u (x) - \pi_ {N} ^ {(\lambda)} u (x)\right) \right\| _ {L ^ {\infty}} \leq \frac {U _ {\theta} ^ {\sigma}}{2 ^ {\sigma - 1} \pi} \sum_ {n = N + 1} ^ {\infty} \frac {\Gamma ((n - \sigma + 1) / 2)}{\Gamma ((n + \sigma + 1) / 2)}.
$$

For simplicity, we denote 

$$
\mathcal {S} _ {n} ^ {\sigma} := \frac {\Gamma ((n - \sigma + 1) / 2)}{\Gamma ((n + \sigma + 1) / 2)}, \quad \mathcal {T} _ {n} ^ {\sigma} := \frac {\Gamma ((n - \sigma + 1) / 2)}{\Gamma ((n + \sigma - 1) / 2)}. \tag {4.19}
$$

A direct calculation leads to the identity: 

$$
\mathcal {T} _ {n} ^ {\sigma} - \mathcal {T} _ {n + 2} ^ {\sigma} = \frac {n + \sigma - 1}{2} \frac {\Gamma ((n - \sigma + 1) / 2)}{\Gamma ((n + \sigma + 1) / 2)} - \frac {n - \sigma + 1}{2} \frac {\Gamma ((n - \sigma + 1) / 2)}{\Gamma ((n + \sigma + 1) / 2)} \tag {4.20}
$$

$$
= (\sigma - 1) \mathcal {S} _ {n} ^ {\sigma}.
$$

We find from [1, (1.1) and Theorem 10] that for $0 \leq a \leq b$ , the ratio 

$$
R _ {b} ^ {a} (z) := \frac {\Gamma (z + a)}{\Gamma (z + b)}, \quad z \geq 0, \tag {4.21}
$$

is decreasing with respect to z. As sigma > 1,we have 

$$
\mathcal {T} _ {N + 2} ^ {\sigma} = R _ {\sigma - 1} ^ {0} ((N - \sigma + 1) / 2 + 1) \leq R _ {\sigma - 1} ^ {0} ((N - \sigma) / 2 + 1) = \mathcal {T} _ {N + 1} ^ {\sigma}. \tag {4.22}
$$

With (4.22), we obtain 

$$
\begin{array}{l} \left\| \varpi_ {\lambda} (x) \left(u (x) - \pi_ {N} ^ {(\lambda)} u (x)\right) \right\| _ {L ^ {\infty}} \leq \frac {U _ {\theta} ^ {\sigma}}{2 ^ {\sigma - 1} \pi} \sum_ {n = N + 1} ^ {\infty} \mathcal {S} _ {n} ^ {\sigma} = \frac {U _ {\theta} ^ {\sigma}}{2 ^ {\sigma - 1} \pi (\sigma - 1)} \sum_ {n = N + 1} ^ {\infty} \left\{\mathcal {T} _ {n} ^ {\sigma} - \mathcal {T} _ {n + 2} ^ {\sigma} \right\} \\ = \frac {U _ {\theta} ^ {\sigma}}{2 ^ {\sigma - 1} \pi (\sigma - 1)} \left\{\mathcal {T} _ {N + 1} ^ {\sigma} + \mathcal {T} _ {N + 2} ^ {\sigma} \right\} \\ \leq \frac {U _ {\theta} ^ {\sigma}}{2 ^ {\sigma - 2} \pi (\sigma - 1)} \frac {\Gamma ((N - \sigma) / 2 + 1)}{\Gamma ((N + \sigma) / 2)}. \\ \end{array}
$$

This completes the proof. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/dea76820-43c7-4583-82cf-5d601e885b86/84c4eaae9a87c834bd9ab97dc6b61ba25af7b1c2714a0133b48cafe8b9d8ba6a.jpg)


Using the Stirling’s formula for Gamma function, we can obtain the explicit convergence order in N . 

Corollary 4.1 Assume that the conditions in Theorem 4.1 hold. As $N  \infty ,$ , for $\lambda \geq 0$ , $\sigma > 1$ , we have 

$$
\left\| \varpi_ {\lambda} (x) \left(u (x) - \pi_ {N} ^ {(\lambda)} u (x)\right) \right\| _ {L ^ {\infty}} \leq C N ^ {1 - \sigma} U _ {\theta} ^ {\sigma},
$$

where C is a generic positive constant independent of N and u. 

Remark 4.1 Consider $u ( x ) = g ( x ) | x - \theta | ^ { \alpha }$ , where $g ( x ) \in C ^ { \infty } ( \Lambda )$ and $\theta \in \Lambda$ , if α is a non-integer greater than 0, then $u \in \mathbb { W } _ { \theta } ^ { \alpha + 1 } ( \Lambda )$ , we can directly apply Theorem 4.2 and (4.17) to get 

$$
\left\| \varpi_ {\lambda} (x) \left(u (x) - \pi_ {N} ^ {(\lambda)} u (x)\right) \right\| _ {L ^ {\infty}} \leq C N ^ {- \alpha} U _ {\theta} ^ {\alpha + 1}, \quad N \gg 1,
$$

where C is a positive constant independent of N and u. 

Remark 4.2 Here, we have presented the optimal (weighted) $L ^ { \infty }$ -estimates of the Gegenbauer polynomial approximations for functions with interior singularities. Wang [33] derived pointwise error estimates for the unweighted $L ^ { \infty }$ -norm of the model function. Theorem 6.5 in [33] provided the optimal unweighted L∞-norm estimate, which relied on the analysis of pointwise errors. For more general cases, the question of how to modify the fractional space to best characterise the optimal unweighted $L ^ { \infty }$ -estimate is left for future work. 

The weighted $L ^ { 2 }$ norm is defined by 

$$
\| u \| _ {L _ {\omega_ {\lambda}} ^ {2} (\Lambda)} := \left(\int_ {\Lambda} | u (x) | ^ {2} \omega_ {\lambda} (x) \mathrm{d} x\right) ^ {\frac {1}{2}}.
$$

Here, we present the optimal (weighted) $L ^ { 2 }$ -estimates of the Gegenbauer polynomial approximations. 

Theorem 4.3 Assume that the conditions in Theorem 4.1 hold. For $\begin{array} { r } { \frac { 1 } { 2 } < \sigma \leq N + 1 } \end{array}$ , then we have the following bounds. 


Table 2 Order of $\left\| \varpi _ { \lambda } ( x ) \left( u ( x ) - \pi _ { N } ^ { ( \lambda ) } u ( x ) \right) \right\| _ { L ^ { \infty } }$ with $\lambda = 0 . 8$ and $u ( x ) = | x | ^ { \alpha } e ^ { \sin { . } }$


<table><tr><td>N</td><td><eq>\alpha = 0.8</eq></td><td>Order</td><td><eq>\alpha = 1.2</eq></td><td>Order</td><td><eq>\alpha = 2.3</eq></td><td>Order</td><td><eq>\alpha = 3.1</eq></td><td>Order</td></tr><tr><td><eq>2^{2}</eq></td><td>1.69e-01</td><td>-</td><td>6.45e-02</td><td>-</td><td>8.59e-03</td><td>-</td><td>2.29e-02</td><td>-</td></tr><tr><td><eq>2^{3}</eq></td><td>1.13e-01</td><td>0.58</td><td>3.58e-02</td><td>0.85</td><td>1.80e-03</td><td>2.25</td><td>1.27e-03</td><td>4.17</td></tr><tr><td><eq>2^{4}</eq></td><td>7.04e-02</td><td>0.68</td><td>1.75e-02</td><td>1.03</td><td>4.51e-04</td><td>2.00</td><td>1.88e-04</td><td>2.76</td></tr><tr><td><eq>2^{5}</eq></td><td>4.22e-02</td><td>0.74</td><td>8.13e-03</td><td>1.11</td><td>1.03e-04</td><td>2.14</td><td>2.53e-05</td><td>2.89</td></tr><tr><td><eq>2^{6}</eq></td><td>2.48e-02</td><td>0.77</td><td>3.66e-03</td><td>1.15</td><td>2.22e-05</td><td>2.21</td><td>3.19e-06</td><td>2.99</td></tr><tr><td><eq>2^{7}</eq></td><td>1.45e-02</td><td>0.77</td><td>1.63e-03</td><td>1.17</td><td>4.65e-06</td><td>2.25</td><td>3.89e-07</td><td>3.04</td></tr><tr><td><eq>2^{8}</eq></td><td>8.57e-03</td><td>0.76</td><td>7.26e-04</td><td>1.16</td><td>9.64e-07</td><td>2.27</td><td>4.64e-08</td><td>3.07</td></tr><tr><td>Pred</td><td></td><td>0.80</td><td></td><td>1.20</td><td></td><td>2.30</td><td></td><td>3.10</td></tr></table>


Table 3 Order of $\left\| \varpi _ { \lambda } ( x ) \left( u ( x ) - \pi _ { N } ^ { ( \lambda ) } u ( x ) \right) \right\| _ { L ^ { \infty } }$ with $\theta = 0 . 2 , \lambda = 0 . 5 \mathrm { a n d } u ( x ) = | x - \theta | ^ { \alpha } e ^ { \sin x }$


<table><tr><td>N</td><td><eq>\alpha = 0.8</eq></td><td>Order</td><td><eq>\alpha = 1.2</eq></td><td>Order</td><td><eq>\alpha = 2.3</eq></td><td>Order</td><td><eq>\alpha = 3.1</eq></td><td>Order</td></tr><tr><td><eq>2^{2}</eq></td><td>2.49e-01</td><td>-</td><td>1.09e-01</td><td>-</td><td>1.48e-02</td><td>-</td><td>3.13e-02</td><td>-</td></tr><tr><td><eq>2^{3}</eq></td><td>1.50e-01</td><td>0.73</td><td>4.97e-02</td><td>1.13</td><td>2.87e-03</td><td>2.37</td><td>2.14e-03</td><td>3.87</td></tr><tr><td><eq>2^{4}</eq></td><td>8.52e-02</td><td>0.81</td><td>2.13e-02</td><td>1.22</td><td>5.54e-04</td><td>2.37</td><td>2.33e-04</td><td>3.20</td></tr><tr><td><eq>2^{5}</eq></td><td>5.07e-02</td><td>0.75</td><td>9.74e-03</td><td>1.13</td><td>1.22e-04</td><td>2.18</td><td>3.00e-05</td><td>2.96</td></tr><tr><td><eq>2^{6}</eq></td><td>2.97e-02</td><td>0.77</td><td>4.36e-03</td><td>1.16</td><td>2.61e-05</td><td>2.23</td><td>3.72e-06</td><td>3.01</td></tr><tr><td><eq>2^{7}</eq></td><td>1.73e-02</td><td>0.78</td><td>1.93e-03</td><td>1.18</td><td>5.44e-06</td><td>2.26</td><td>4.50e-07</td><td>3.05</td></tr><tr><td><eq>2^{8}</eq></td><td>1.02e-02</td><td>0.77</td><td>8.57e-04</td><td>1.17</td><td>1.12e-06</td><td>2.28</td><td>5.34e-08</td><td>3.08</td></tr><tr><td>Pred</td><td></td><td>0.80</td><td></td><td>1.20</td><td></td><td>2.30</td><td></td><td>3.10</td></tr></table>

(i) $H 0 \leq \lambda < 1 / 2$ , 

$$
\begin{array}{l} \left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Lambda)} \leq \left\{\frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{2 ^ {2 \lambda + 2 \sigma - 1} \pi \Gamma (N + 2)} \frac {\Gamma^ {2} ((N - \sigma) / 2 + 1)}{\Gamma^ {2} ((N + \sigma) / 2 + \lambda + 1)} \right. \\ + \frac {2}{(2 \sigma - 1) \pi} \widetilde {\Upsilon} _ {N + 2} ^ {(\lambda , \sigma)} \max \left\{1, \frac {N + \lambda + 1}{N + \sigma + 1} \right\} \frac {\Gamma (N - \sigma + 2)}{\Gamma (N + \sigma + 1)} \bigg \} ^ {1 / 2} \mathbb {U} _ {\theta} ^ {\sigma}, \tag {4.23} \\ \end{array}
$$

where 

$$
\widetilde {\Upsilon} _ {n} ^ {(\lambda , \sigma)} = \frac {\Gamma (n + 2 \lambda) \Gamma^ {2} ((n + \sigma) / 2 + 1)}{2 ^ {2 \lambda - 1} \Gamma (n + 1) \Gamma^ {2} ((n + \sigma + 1) / 2 + \lambda)}. \tag {4.24}
$$

(ii) $I f \lambda \ge 1 / 2$ , 

$$
\begin{array}{l} \left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Lambda)} \leq \left\{\frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{2 ^ {2 \lambda + 2 \sigma - 1} \pi \Gamma (N + 2)} \frac {\Gamma^ {2} ((N - \sigma) / 2 + 1)}{\Gamma^ {2} ((N + \sigma) / 2 + \lambda + 1)} \right. \\ + \frac {2}{(2 \sigma - 1) \pi} \max \left\{1, \frac {N + \lambda + 1}{N + \sigma + 1} \right\} \frac {\Gamma (N - \sigma + 2)}{\Gamma (N + \sigma + 1)} \bigg \} ^ {1 / 2} \mathbb {U} _ {\theta} ^ {\sigma}. \tag {4.25} \\ \end{array}
$$

Proof By the orthogonality of Gegenbauer polynomials, we derive from (4.1) and (4.4), 

$$
\begin{array}{l} \left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Lambda)} ^ {2} = \sum_ {n = N + 1} ^ {\infty} \frac {2 ^ {2 \lambda - 1} n ! \Gamma^ {2} (\lambda + 1 / 2)}{(n + \lambda) \Gamma (n + 2 \lambda)} \left| \hat {u} _ {n} ^ {(\lambda)} \right| ^ {2} \\ \leq \frac {\left(\mathbb {U} _ {\theta} ^ {\sigma}\right) ^ {2}}{2 ^ {2 \lambda + 2 \sigma - 1} \pi} \left\{\frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{\Gamma (N + 2)} \frac {\Gamma^ {2} ((N - \sigma) / 2 + 1)}{\Gamma^ {2} ((N + \sigma) / 2 + \lambda + 1)} \right. \\ \left. + \sum_ {n = N + 2} ^ {\infty} \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)} \frac {\Gamma^ {2} ((n - \sigma + 1) / 2)}{\Gamma^ {2} ((n + \sigma + 1) / 2 + \lambda)} \right\} \\ \leq Q _ {1} + Q _ {2}, \tag {4.26} \\ \end{array}
$$

where 

$$
Q _ {1} := \frac {\left(\mathbb {U} _ {\theta} ^ {\sigma}\right) ^ {2}}{2 ^ {2 \lambda + 2 \sigma - 1} \pi} \frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{\Gamma (N + 2)} \frac {\Gamma^ {2} ((N - \sigma) / 2 + 1)}{\Gamma^ {2} ((N + \sigma) / 2 + \lambda + 1)}, \tag {4.27}
$$

$$
Q _ {2} := \frac {\left(\mathbb {U} _ {\theta} ^ {\sigma}\right) ^ {2}}{2 ^ {2 \sigma} \pi} \sum_ {n = N + 2} ^ {\infty} \widetilde {\Upsilon} _ {n} ^ {(\lambda , \sigma)} \frac {(n + \lambda) \Gamma^ {2} ((n - \sigma + 1) / 2)}{\Gamma^ {2} ((n + \sigma) / 2 + 1)}. \tag {4.28}
$$

By simple calculations lead to 

$$
\frac {\widetilde {\Upsilon} _ {n} ^ {(\lambda , \sigma)}}{\widetilde {\Upsilon} _ {n + 2} ^ {(\lambda , \sigma)}} = \frac {(n + 1) (n + 2) (n + \sigma + 2 \lambda + 1) ^ {2}}{(n + 2 \lambda) (n + 2 \lambda + 1) (n + \sigma + 2) ^ {2}} \tag {4.29}
$$

$$
= 1 - \frac {(2 \lambda - 1) \mathcal {W}}{(n + 2 \lambda) (n + 2 \lambda + 1) (n + \sigma + 2) ^ {2}},
$$

where 

$$
\mathcal {W} = 2 \lambda \sigma^ {2} + (4 n + 8) \lambda \sigma + (2 n + 4) \lambda + (2 n + 2) \sigma^ {2} + (2 n ^ {2} + 6 n + 4) \sigma + n ^ {2} + 3 n + 2 > 0.
$$

We immediately obtain that the sequence $\left\{ \widetilde { \Upsilon } _ { n } ^ { ( \lambda , \sigma ) } \right\} _ { n > N + 1 }$ is decreasing for $- 1 / 2 < \lambda < 1 / 2$ , n N 1 and increasing for $\lambda \ge 1 / 2$ ≥. From (4.17), we have for $\dot { \lambda } \geq 1 / 2$ 

$$
\lim _ {n \rightarrow \infty} \widetilde {\Upsilon} _ {n} ^ {(\lambda , \sigma)} = 1, \quad \widetilde {\Upsilon} _ {n} ^ {(\lambda , \sigma)} <   \widetilde {\Upsilon} _ {n + 2} ^ {(\lambda , \sigma)} <   \dots \leq 1. \tag {4.30}
$$

By (4.21), we obtain 

$$
\frac {\Gamma ((n - \sigma + 1) / 2)}{\Gamma ((n + \sigma) / 2 + 1)} \leq \frac {\Gamma ((n - \sigma) / 2)}{\Gamma ((n + \sigma + 1) / 2)}.
$$

Using the identity: 

$$
\Gamma (2 z) = \pi^ {- 1 / 2} 2 ^ {2 z - 1} \Gamma (z) \Gamma (z + 1 / 2), \tag {4.31}
$$

we derive 

$$
\begin{array}{l} \frac {(n + \lambda) \Gamma^ {2} ((n - \sigma + 1) / 2)}{\Gamma^ {2} ((n + \sigma) / 2 + 1)} = \frac {2 (n + \lambda) \Gamma^ {2} ((n - \sigma + 1) / 2)}{(n + \sigma) \Gamma ((n + \sigma) / 2) \Gamma ((n + \sigma) / 2 + 1)} \\ \leq \frac {2 (n + \lambda)}{(n + \sigma)} \frac {\Gamma ((n - \sigma + 1) / 2) \Gamma ((n - \sigma) / 2)}{\Gamma ((n + \sigma) / 2) \Gamma ((n + \sigma + 1) / 2)} \\ = \frac {2 ^ {2 \sigma + 1} (n + \lambda)}{(n + \sigma)} \frac {\Gamma (n - \sigma)}{\Gamma (n + \sigma)} \\ = \frac {2 ^ {2 \sigma + 1} (n + \lambda)}{(2 \sigma - 1) (n + \sigma)} \left\{\frac {\Gamma (n - \sigma)}{\Gamma (n + \sigma - 1)} - \frac {\Gamma (n - \sigma + 1)}{\Gamma (n + \sigma)} \right\}, \tag {4.32} \\ \end{array}
$$

where we used the property: $z \Gamma ( z ) = \Gamma ( z + 1 )$ . Here for brevity, we denote 

$$
\widetilde {R} _ {n} ^ {(\lambda , \sigma)} := \frac {n + \lambda}{n + \sigma}, \tag {4.33}
$$

then for $- 1 / 2 < \lambda \leq \sigma , \widetilde { R } _ { n } ^ { ( \lambda , \sigma ) } \leq 1$ , for $\lambda > \sigma > 1 / 2 , \left\{ \widetilde { R } _ { n } ^ { ( \lambda , \sigma ) } \right\} _ { n > N + 1 }$ is decreasing. n N 1 Hence from (4.33) and (4.32), we can rewrite (4.26) for $- 1 / 2 < \lambda < 1 \bar { 7 } 2$ +  as follows: 

$$
Q _ {2} \leq \frac {\left(\mathbb {U} _ {\theta} ^ {\sigma}\right) ^ {2}}{2 ^ {2 \sigma} \pi} \widetilde {\Upsilon} _ {N + 2} ^ {(\lambda , \sigma)} \sum_ {n = N + 2} ^ {\infty} \frac {2 ^ {2 \sigma + 1} (n + \lambda)}{(2 \sigma - 1) (n + \sigma)} \left\{\frac {\Gamma (n - \sigma)}{\Gamma (n + \sigma - 1)} - \frac {\Gamma (n - \sigma + 1)}{\Gamma (n + \sigma)} \right\} \tag {4.34}
$$

$$
\leq \frac {2 \left(\mathbb {U} _ {\theta} ^ {\sigma}\right) ^ {2}}{(2 \sigma - 1) \pi} \widetilde {\Upsilon} _ {N + 2} ^ {(\lambda , \sigma)} \max \left\{1, \frac {N + \lambda + 1}{N + \sigma + 1} \right\} \frac {\Gamma (N - \sigma + 2)}{\Gamma (N + \sigma + 1)}.
$$

Similarly, for $\lambda \ge 1 / 2$ , we can rewrite (4.26) using (4.30) as follows: 

$$
Q _ {2} \leq \frac {2 \left(\mathbb {U} _ {\theta} ^ {\sigma}\right) ^ {2}}{(2 \sigma - 1) \pi} \max \left\{1, \frac {N + \lambda + 1}{N + \sigma + 1} \right\} \frac {\Gamma (N - \sigma + 2)}{\Gamma (N + \sigma + 1)}. \tag {4.35}
$$

We get (4.23) from (4.26), (4.27) and (4.34). Also, we get (4.25) from (4.26), (4.27) and (4.35). This completes the proof. □ 

Corollary 4.2 Assume that the conditions in Theorem 4.1 hold. As $N \to \infty , f o r \lambda \geq 0$ $\sigma > \frac { 1 } { 2 }$ , we have 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Lambda)} \leq C N ^ {1 / 2 - \sigma} \mathbb {U} _ {\theta} ^ {\sigma},
$$

where C is a generic positive constant independent of N and u. 

Remark 4.3 Consider $u ( x ) = g ( x ) | x - \theta | ^ { \alpha }$ , where $g ( x ) \in C ^ { \infty } ( \Lambda )$ and $\theta \in \Lambda$ , if α is a non-integer greater than 0, then $u \in \mathbb { W } _ { \theta } ^ { \alpha + 1 } ( \Lambda )$ , we can directly apply Theorem 4.3 and (4.17) to get 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Lambda)} \leq C N ^ {- \alpha - 1 / 2} U _ {\theta} ^ {\alpha + 1}, \quad N \gg 1,
$$

where C is a positive constant independent of N and u. 

# 4.3 Approximation Results for Functions with Endpoint Singularities

We now study the case of endpoint singularities. To fix the idea, we focus on the left endpoint singularityconditions, with $x = - 1$ . Fo ce of $x = 1$ , the following theorems hold under correspondingHere, we present the main results below. $\mathbb { U } _ { 1 - } ^ { \sigma , k }$ $\mathbb { U } _ { - 1 + } ^ { \sigma , k }$ 


Table 4 Order of u(x ) − π (λ)N u(x )  L2ω (- $\left. u ( x ) - \pi _ { N } ^ { ( \lambda ) } u ( x ) \right. _ { L _ { \omega _ { \lambda } } ^ { 2 } ( \Omega ) }$ with $\lambda = 0 . 3$ and $u ( x ) = | x | ^ { \alpha } e ^ { \sin { x } }$


<table><tr><td>N</td><td><eq>\alpha = 0.8</eq></td><td>Order</td><td><eq>\alpha = 1.2</eq></td><td>Order</td><td><eq>\alpha = 2.3</eq></td><td>Order</td><td><eq>\alpha = 3.1</eq></td><td>Order</td></tr><tr><td><eq>2^{2}</eq></td><td>7.40e-02</td><td>-</td><td>4.32e-02</td><td>-</td><td>9.43e-03</td><td>-</td><td>3.09e-02</td><td>-</td></tr><tr><td><eq>2^{3}</eq></td><td>3.41e-02</td><td>1.12</td><td>1.46e-02</td><td>1.57</td><td>1.22e-03</td><td>2.95</td><td>1.10e-03</td><td>4.82</td></tr><tr><td><eq>2^{4}</eq></td><td>1.52e-02</td><td>1.17</td><td>5.04e-03</td><td>1.53</td><td>2.01e-04</td><td>2.61</td><td>1.02e-04</td><td>3.43</td></tr><tr><td><eq>2^{5}</eq></td><td>6.48e-03</td><td>1.23</td><td>1.65e-03</td><td>1.61</td><td>3.18e-05</td><td>2.66</td><td>9.40e-06</td><td>3.44</td></tr><tr><td><eq>2^{6}</eq></td><td>2.68e-03</td><td>1.27</td><td>5.25e-04</td><td>1.65</td><td>4.81e-06</td><td>2.72</td><td>8.29e-07</td><td>3.50</td></tr><tr><td><eq>2^{7}</eq></td><td>1.08e-03</td><td>1.32</td><td>1.63e-04</td><td>1.69</td><td>7.10e-07</td><td>2.76</td><td>7.08e-08</td><td>3.55</td></tr><tr><td><eq>2^{8}</eq></td><td>4.03e-04</td><td>1.42</td><td>4.83e-05</td><td>1.75</td><td>1.03e-07</td><td>2.79</td><td>5.94e-09</td><td>3.58</td></tr><tr><td>Pred</td><td></td><td>1.30</td><td></td><td>1.70</td><td></td><td>2.80</td><td></td><td>3.60</td></tr></table>


Table 5 Order of u(x ) − π (λ)N u(x )  L2ω ( $\left. u ( x ) - \pi _ { N } ^ { ( \lambda ) } u ( x ) \right. _ { L _ { \omega _ { \lambda } } ^ { 2 } ( \Omega ) }$ with $\theta = 0 . 2 , \lambda = 0 . 5 \mathrm { a n d } u ( x ) = | x - \theta | ^ { \alpha } e ^ { \sin x }$


<table><tr><td>N</td><td><eq>\alpha = 0.8</eq></td><td>Order</td><td><eq>\alpha = 1.2</eq></td><td>Order</td><td><eq>\alpha = 2.3</eq></td><td>Order</td><td><eq>\alpha = 3.1</eq></td><td>Order</td></tr><tr><td><eq>2^{2}</eq></td><td>1.08e-01</td><td>-</td><td>6.42e-02</td><td>-</td><td>1.22e-02</td><td>-</td><td>3.20e-02</td><td>-</td></tr><tr><td><eq>2^{3}</eq></td><td>4.48e-02</td><td>1.27</td><td>1.95e-02</td><td>1.72</td><td>1.68e-03</td><td>2.87</td><td>1.46e-03</td><td>4.45</td></tr><tr><td><eq>2^{4}</eq></td><td>1.80e-02</td><td>1.31</td><td>5.93e-03</td><td>1.71</td><td>2.32e-04</td><td>2.85</td><td>1.17e-04</td><td>3.65</td></tr><tr><td><eq>2^{5}</eq></td><td>7.70e-03</td><td>1.23</td><td>1.95e-03</td><td>1.61</td><td>3.67e-05</td><td>2.66</td><td>1.07e-05</td><td>3.45</td></tr><tr><td><eq>2^{6}</eq></td><td>3.20e-03</td><td>1.27</td><td>6.21e-04</td><td>1.65</td><td>5.57e-06</td><td>2.72</td><td>9.46e-07</td><td>3.50</td></tr><tr><td><eq>2^{7}</eq></td><td>1.29e-03</td><td>1.31</td><td>1.93e-04</td><td>1.68</td><td>8.25e-07</td><td>2.75</td><td>8.11e-08</td><td>3.54</td></tr><tr><td><eq>2^{8}</eq></td><td>4.92e-04</td><td>1.39</td><td>5.80e-05</td><td>1.74</td><td>1.20e-07</td><td>2.79</td><td>6.81e-09</td><td>3.58</td></tr><tr><td>Pred</td><td></td><td>1.30</td><td></td><td>1.70</td><td></td><td>2.80</td><td></td><td>3.60</td></tr></table>

Theorem 4.4 For some m, $k \in  { \mathbb { N } } _ { 0 } , s \in ( 0 , 1 )$ and $\sigma : = m + s$ , assume that $u \in \mathbb { W } _ { - 1 + } ^ { \sigma , k } ( \Lambda )$ Then for $n > \sigma + k$ and $\sigma + \lambda > 1 / 2$ , we have 

$$
\begin{array}{l} \hat {u} _ {n} ^ {(\lambda)} = \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{n ! \Gamma (\lambda + 1 / 2)} \left\{\frac {2 ^ {1 - 2 \lambda - \sigma - k}}{\Gamma (\lambda + \sigma + k + 1 / 2)} \int_ {- 1} ^ {1} ^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + k)} (x)   \mathrm{d} v ^ {(k)} (x) \right. \\ \left. + \sum_ {j = 0} ^ {k} (- 1) ^ {n + j} \widehat {C} _ {n, \sigma + j} ^ {(\lambda)} \sin (\sigma \pi) v ^ {(j)} (- 1 +) \right\}, \tag {4.36} \\ \end{array}
$$

where $v ( x ) = I _ { - 1 + } ^ { 1 - s } u ^ { ( m ) } ( x )$ , and 

$$
\widehat {C} _ {n, \rho} ^ {(\lambda)} := \frac {2 ^ {\rho} \Gamma (\rho + \lambda - 1 / 2)}{\pi} \frac {\Gamma (n - \rho + 1)}{\Gamma (n + \rho + 2 \lambda)}. \tag {4.37}
$$

Moreover, we have the following bound: 

$$
\begin{array}{l} \left| \hat {u} _ {n} ^ {(\lambda)} \right| \leq \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1) \Gamma (\lambda + 1 / 2)} \left\{\frac {2 ^ {1 - 2 \lambda - \sigma - k} \Gamma ((n - \sigma - k + 1) / 2)}{\sqrt {\pi} \Gamma ((n + \sigma + k + 1) / 2 + \lambda)} \right. \\ \left. + \frac {2 ^ {\sigma} \Gamma (\sigma + \lambda - 1 / 2)}{\pi} \frac {n + \lambda}{n - \sigma - k + 1} \frac {\Gamma (n - \sigma + 1)}{\Gamma (n + \sigma + 2 \lambda)} \right\} \mathbb {U} _ {- 1 +} ^ {\sigma , k}, \tag {4.38} \\ \end{array}
$$

where 

(i) For $1 / 2 < \lambda + \sigma + k < 1$ , 

$$
\mathbb {U} _ {- 1 +} ^ {\sigma , k} := \int_ {- 1} ^ {1} \omega_ {(\lambda + \sigma + k) / 2} | \mathrm{d} v ^ {(k)} (x) | + | v (- 1 +) |.
$$

$F o r \lambda + \sigma + k \ge 1 , \mathbb { U } _ { - 1 + } ^ { \sigma , k } : = U _ { - 1 + } ^ { \sigma , k }$ . 

Proof As the conditions in Definition 3.4, we derive from [18, (6.13)] with (3.2) that 

$$
\hat {u} _ {n} ^ {(\lambda)} = \frac {2 ^ {1 - 2 \lambda - \sigma} (n + \lambda) \Gamma (n + 2 \lambda)}{n ! \Gamma (\lambda + 1 / 2) \Gamma (\lambda + \sigma + 1 / 2)} \left\{\int_ {- 1} ^ {1} ^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (x) v ^ {\prime} (x) d x + \lim _ {x \rightarrow - 1 +} \left\{^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (x) v (x) \right\}\right\}. \tag {4.39}
$$

Taking $\sigma  \sigma + j + 1$ in (4.12), we obtain for $n > \sigma + j + 1$ , 

$$
{ } ^ { r } \mathcal { G } _ { n } ^ { ( \lambda , \sigma + j ) } ( x ) = - \frac { 1 } { 2 ( \lambda + \sigma + j ) + 1 } \left\{ { } ^ { r } \mathcal { G } _ { n } ^ { ( \lambda , \sigma + j + 1 ) } ( x ) \right\} ^ { \prime } . \tag {4.40}
$$

Since $v \in \mathbb { W } ^ { k + 1 } ( \Lambda )$ , using (4.40) and Lemma 3.1, yields 

$$
\begin{array}{l} \int_ {- 1} ^ {1}   ^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma)} (x) v ^ {\prime} (x) \mathrm{d} x = - \frac {1}{2 (\lambda + \sigma) + 1} \int_ {- 1} ^ {1} \left\{^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + 1)} (x) \right\} ^ {\prime} v ^ {\prime} (x) \mathrm{d} x \\ = \frac {1}{2 (\lambda + \sigma) + 1} \left\{\int_ {- 1} ^ {1}   ^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + 1)} (x) v ^ {\prime \prime} (x) \mathrm{d} x + \left\{^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + 1)} (x) v ^ {\prime} (x) \right\} \big | _ {x = - 1 +} \right\} \\ = \dots = \frac {\Gamma (\lambda + \sigma + 1 / 2)}{2 ^ {k - 1} \Gamma (\lambda + \sigma + k - 1 / 2)} \int_ {- 1} ^ {1} r \mathcal {G} _ {n} ^ {(\lambda , \sigma + k - 1)} (x) v ^ {(k)} (x) d x \\ + \sum_ {j = 1} ^ {k - 1} \frac {\Gamma (\lambda + \sigma + 1 / 2)}{2 ^ {j} \Gamma (\lambda + \sigma + j + 1 / 2)} \left\{^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + j)} (x) v ^ {(j)} (x) \right\} | _ {x = - 1 +} \\ = \frac {\Gamma (\lambda + \sigma + 1 / 2)}{2 ^ {k} \Gamma (\lambda + \sigma + k + 1 / 2)} \int_ {- 1} ^ {1} r \mathcal {G} _ {n} ^ {(\lambda , \sigma + k)} (x)   \mathrm{d} v ^ {(k)} (x) \\ + \sum_ {j = 1} ^ {k} \frac {\Gamma (\lambda + \sigma + 1 / 2)}{2 ^ {j} \Gamma (\lambda + \sigma + j + 1 / 2)} \left\{^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + j)} (x) v ^ {(j)} (x) \right\} | _ {x = - 1 +}. \tag {4.41} \\ \end{array}
$$

Substituting (4.41) into (4.39), then we can get 

$$
\begin{array}{l} \hat {u} _ {n} ^ {(\lambda)} = \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{n ! \Gamma (\lambda + 1 / 2)} \left\{\frac {2 ^ {1 - 2 \lambda - \sigma}}{\Gamma (\lambda + \sigma + k + 1 / 2)} \int_ {- 1} ^ {1} ^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + k)} (x) d v ^ {(k)} (x) \right. \\ \left. + \sum_ {j = 0} ^ {k} \frac {2 ^ {1 - 2 \lambda - \sigma - j}}{\Gamma (\lambda + \sigma + j + 1 / 2)} \left\{^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + j)} (x) v ^ {(j)} (x) \right\} \right| _ {x = - 1 +} \Bigg \}. \tag {4.42} \\ \end{array}
$$

From [18, (6.15)], we readily can find 

$$
\lim _ {x \to - 1 +} \left\{^ {r} \mathcal {G} _ {n} ^ {(\lambda , \sigma + j)} (x) \right\} = (- 1) ^ {n + j} 2 ^ {2 (\sigma + j + \lambda) - 1} \frac {\sin (\sigma \pi)}{\pi} \Psi , \tag {4.43}
$$

where  is 

$$
\Psi = \frac {\Gamma (\sigma + j + \lambda - 1 / 2) \Gamma (\sigma + j + \lambda + 1 / 2) \Gamma (n - \sigma - j + 1)}{\Gamma (n + 2 \lambda + \sigma + j)}. \tag {4.44}
$$

Thus, using (4.42) and (4.43), we obtain (4.36). 

According to (4.37), from straightforward calculations, we find 

$$
\begin{array}{l} \sum_ {j = 0} ^ {k} \frac {\widehat {C} _ {n , \sigma + j} ^ {(\lambda)}}{\widehat {C} _ {n , \sigma} ^ {(\lambda)}} = \sum_ {j = 0} ^ {k} 2 ^ {j} \frac {\Gamma (\sigma + j + \lambda - 1 / 2)}{\Gamma (\sigma + \lambda - 1 / 2)} \frac {\Gamma (n - \sigma - j + 1)}{\Gamma (n - \sigma + 1)} \frac {\Gamma (n + 2 \lambda + \sigma)}{\Gamma (n + 2 \lambda + \sigma + j)} \\ = 1 + \sum_ {j = 1} ^ {k} \left(\prod_ {i = 1} ^ {j} \frac {2 (\sigma + \lambda + i - 3 / 2)}{(n - \sigma - i + 1) (n + 2 \lambda + \sigma + i - 1)}\right) \\ \leq 1 + \sum_ {j = 1} ^ {k} \left(\frac {\sigma + \lambda + k - 1}{n + \lambda}\right) ^ {j} = \sum_ {j = 0} ^ {k} \left(\frac {\sigma + \lambda + k - 1}{n + \lambda}\right) ^ {j} \\ = \frac {n + \lambda}{n - \sigma - k + 1}, \tag {4.45} \\ \end{array}
$$

where we used the fact: for $1 \leq i \leq j \leq k$ and $n > \sigma + k$ , 

$$
(n - \sigma - i + 1) (n + 2 \lambda + \sigma + i - 1) \geq (n - \sigma - k + 1) (n + 2 \lambda + \sigma + k - 1) \geq 2 (n + \lambda).
$$

It follows from Lemma 2.1 and (4.36), we can get the bound 

$$
\left| \hat {u} _ {n} ^ {(\lambda)} \right| \leq \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1) \Gamma (\lambda + 1 / 2)} \left\{\frac {2 ^ {1 - 2 \lambda - \sigma - k} \Gamma ((n - \sigma - k + 1) / 2)}{\sqrt {\pi} \Gamma ((n + \sigma + k + 1) / 2 + \lambda)} + \sum_ {j = 0} ^ {k} \widehat {C} _ {n, \sigma + j} ^ {(\lambda)} \right\} \mathbb {U} _ {- 1 +} ^ {\sigma , k}, \tag {4.46}
$$

then from (4.45) and (4.46), we obtain (4.38) immediately. 

Theorem 4.5 For som $k \in  { \mathbb { N } } _ { 0 } , s \in ( 0 , 1 )$ $\sigma : = m + s$ , assume that $u \in \mathbb { W } _ { - 1 + } ^ { \sigma , k } ( \Lambda )$ $0 \leq \lambda < \sigma + k \leq N + 1$ 

$$
\begin{array}{l} \left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L ^ {\infty}} \leq \left\{\frac {2 ^ {1 - \lambda} N ^ {\lambda + 1 - \sigma - k}}{(\sigma + k - \lambda - 1) \sqrt {\pi} \Gamma (\lambda + 1 / 2)} \right. \\ \left. + \frac {N + \lambda + 1}{N - \sigma - k + 2} \frac {2 ^ {\sigma - 1} \Gamma (\sigma + \lambda - 1 / 2)}{\pi (\sigma - 1) \Gamma (\lambda + 1 / 2)} \frac {\Gamma (N - \sigma + 2)}{\Gamma (N + \sigma)} \right\} \mathbb {U} _ {- 1 +} ^ {\sigma , k}. \tag {4.47} \\ \end{array}
$$

Proof By (2.1), (4.1) and (4.38), we obtain 

$$
\begin{array}{l} \left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L ^ {\infty}} \leq \sum_ {n = N + 1} ^ {\infty} \left| \hat {u} _ {n} ^ {(\lambda)} \right| \\ \leq \sum_ {n = N + 1} ^ {\infty} \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1) \Gamma (\lambda + 1 / 2)} \left\{\frac {2 ^ {1 - 2 \lambda - \sigma - k} \Gamma ((n - \sigma - k + 1) / 2)}{\sqrt {\pi} \Gamma ((n + \sigma + k + 1) / 2 + \lambda)} \right. \\ \left. + \frac {2 ^ {\sigma} \Gamma (\sigma + \lambda - 1 / 2)}{\pi} \frac {n + \lambda}{n - \sigma - k + 1} \frac {\Gamma (n - \sigma + 1)}{\Gamma (n + \sigma + 2 \lambda)} \right\} \mathbb {U} _ {- 1 +} ^ {\sigma , k} \\ = \widehat {Q} _ {1} + \widehat {Q} _ {2}, \tag {4.48} \\ \end{array}
$$

where 

$$
\widehat {Q _ {1}} := \frac {2 ^ {1 - \lambda}}{\sqrt {\pi} \Gamma (\lambda + 1 / 2)} \sum_ {n = N + 1} ^ {\infty} \left\{\frac {(n + \lambda) \Gamma (n + 2 \lambda)}{2 ^ {\sigma + k + \lambda} \Gamma (n + 1)} \frac {\Gamma ((n - \sigma - k + 1) / 2)}{\Gamma ((n + \sigma + k + 1) / 2 + \lambda)} \right\} \mathbb {U} _ {- 1 +} ^ {\sigma , k}, \tag {4.49}
$$

$$
\widehat {Q _ {2}} := \frac {2 ^ {\sigma} \Gamma (\sigma + \lambda - 1 / 2)}{\pi \Gamma (\lambda + 1 / 2)} \sum_ {n = N + 1} ^ {\infty} \left\{\frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)} \frac {n + \lambda}{n - \sigma - k + 1} \right. \tag {4.50}
$$

$$
\times \frac {\Gamma (n + \sigma)}{\Gamma (n + \sigma + 2 \lambda)} \frac {\Gamma (n - \sigma + 1)}{\Gamma (n + \sigma)} \Bigg \} \mathbb {U} _ {- 1 +} ^ {\sigma , k}.
$$

For brevity, we denote 

$$
M _ {n} ^ {(\lambda , \sigma , k)} := \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{2 ^ {\sigma + k + \lambda} \Gamma (n + 1)} \frac {\Gamma ((n - \sigma - k + 1) / 2)}{\Gamma ((n + \sigma + k + 1) / 2 + \lambda)}, \tag {4.51}
$$

$$
\widetilde {M} _ {n} ^ {(\lambda , \sigma)} := \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)} \frac {\Gamma (n + \sigma)}{\Gamma (n + \sigma + 2 \lambda)}, \tag {4.52}
$$

$$
\widehat {M} _ {n} ^ {(\lambda , \sigma , k)} := \frac {n + \lambda}{n - \sigma - k + 1}. \tag {4.53}
$$

We first consider (4.51). Using the Stirling’s formula yields 

$$
M _ {n} ^ {(\lambda , \sigma , k)} \cong \frac {1}{n ^ {\sigma + k - \lambda}}. \tag {4.54}
$$

Thus, we have 

$$
\sum_ {n = N + 1} ^ {\infty} \frac {1}{n ^ {\sigma + k - \lambda}} \leq \int_ {N} ^ {\infty} \frac {\mathrm{d} x}{x ^ {\sigma + k - \lambda}} = \frac {1}{(\sigma + k - \lambda - 1) N ^ {\sigma + k - \lambda - 1}}. \tag {4.55}
$$

Furthermore, we consider (4.52). By simple calculations lead to 

$$
\frac {\widetilde {M} _ {n} ^ {(\lambda , \sigma)}}{\widetilde {M} _ {n + 1} ^ {(\lambda , \sigma)}} = \frac {(n + 1) (n + \lambda) (n + \sigma + 2 \lambda)}{(n + 2 \lambda) (n + 1 + \lambda) (n + \sigma)} = 1 - \frac {\lambda (\sigma - n - 2 \lambda + 2 \lambda \sigma + 2 n \sigma)}{(n + 2 \lambda) (n + 1 + \lambda) (n + \sigma)} <   1,
$$

so the sequence $\left\{ { \widetilde M } _ { n } ^ { ( \lambda , \sigma ) } \right\} _ { n \ge N + 1 }$ is increasing for $\lambda \geq 0$ . From (4.17), we have 

$$
\lim _ {n \to \infty} \widetilde {M} _ {n} ^ {(\lambda , \sigma)} = 1. \tag {4.56}
$$

Finally, we consider (4.53). A direct calculation leads to 

$$
\frac {\check {D} _ {n} ^ {(\lambda , \sigma , k)}}{\check {D} _ {n + 1} ^ {(\lambda , \sigma , k)}} = \left(\frac {(n + \lambda) (n - \sigma - k + 2)}{(n + \lambda + 1) (n - \sigma - k + 1)}\right) ^ {2} > 1, \tag {4.57}
$$

so the sequence $\left\{ { \widehat M } _ { n } ^ { ( \lambda , \sigma , k ) } \right\} _ { n > N + 1 }$ is decreasing for $\lambda \geq 0$ . n N 1 

≥ +According to [39, Thm. 4.2], we can easily obtain 

$$
\sum_ {n = N + 1} ^ {\infty} \frac {\Gamma (n - \sigma + 1)}{\Gamma (n + \sigma)} = \frac {1}{2 (\sigma - 1)} \sum_ {n = N + 1} ^ {\infty} \left\{\frac {\Gamma (n - \sigma + 1)}{\Gamma (n + \sigma - 1)} - \frac {\Gamma (n - \sigma + 2)}{\Gamma (n + \sigma)} \right\} \tag {4.58}
$$

$$
= \frac {1}{2 (\sigma - 1)} \frac {\Gamma (N - \sigma + 2)}{\Gamma (N + \sigma)}.
$$

Using (4.54), and (4.55), we can rewrite (4.49) as 

$$
\widehat {Q _ {1}} \leq \frac {2 ^ {1 - \lambda} N ^ {\lambda + 1 - \sigma - k}}{(\sigma + k - \lambda - 1) \sqrt {\pi} \Gamma (\lambda + 1 / 2)} \mathbb {U} _ {- 1 +} ^ {\sigma , k}. \tag {4.59}
$$


Table 6 Order of $\left\| u ( x ) - \pi _ { N } ^ { ( \lambda ) } u ( x ) \right\| _ { L ^ { \infty } }$ with λ = 0.6 and $u = ( 1 + x ) ^ { \alpha } e ^ { \sin { \it { \cdot } } }$


<table><tr><td>N</td><td><eq>\alpha = 0.8</eq></td><td>Order</td><td><eq>\alpha = 1.2</eq></td><td>Order</td><td><eq>\alpha = 2.3</eq></td><td>Order</td><td><eq>\alpha = 3.1</eq></td><td>Order</td></tr><tr><td><eq>2^{2}</eq></td><td>1.87e-02</td><td>-</td><td>3.41e-02</td><td>-</td><td>8.43e-02</td><td>-</td><td>1.06e-01</td><td>-</td></tr><tr><td><eq>2^{3}</eq></td><td>4.81e-03</td><td>1.96</td><td>9.72e-04</td><td>5.13</td><td>3.43e-04</td><td>7.94</td><td>7.45e-04</td><td>7.15</td></tr><tr><td><eq>2^{4}</eq></td><td>1.75e-03</td><td>1.46</td><td>2.28e-04</td><td>2.09</td><td>4.25e-06</td><td>6.33</td><td>1.50e-07</td><td>12.28</td></tr><tr><td><eq>2^{5}</eq></td><td>6.11e-04</td><td>1.52</td><td>4.69e-05</td><td>2.28</td><td>2.04e-07</td><td>4.38</td><td>2.47e-09</td><td>5.93</td></tr><tr><td><eq>2^{6}</eq></td><td>2.07e-04</td><td>1.56</td><td>9.26e-06</td><td>2.34</td><td>9.08e-09</td><td>4.49</td><td>3.71e-11</td><td>6.05</td></tr><tr><td><eq>2^{7}</eq></td><td>6.93e-05</td><td>1.58</td><td>1.79e-06</td><td>2.37</td><td>3.89e-10</td><td>4.54</td><td>5.32e-13</td><td>6.12</td></tr><tr><td><eq>2^{8}</eq></td><td>2.30e-05</td><td>1.59</td><td>3.43e-07</td><td>2.38</td><td>1.64e-11</td><td>4.57</td><td>7.43e-15</td><td>6.16</td></tr><tr><td>Pred</td><td></td><td>1.60</td><td></td><td>2.40</td><td></td><td>4.60</td><td></td><td>6.20</td></tr></table>


Table 7 Order of $\left\| u ( x ) - \pi _ { N } ^ { ( \lambda ) } u ( x ) \right\| _ { L ^ { \infty } }$ with $\lambda = 1 . 2 \operatorname { a n d } u = ( 1 + x ) ^ { \alpha } e ^ { \sin x }$


<table><tr><td>N</td><td><eq>\alpha = 0.8</eq></td><td>Order</td><td><eq>\alpha = 1.2</eq></td><td>Order</td><td><eq>\alpha = 2.3</eq></td><td>Order</td><td><eq>\alpha = 3.1</eq></td><td>Order</td></tr><tr><td><eq>2^{2}</eq></td><td>2.95e-02</td><td>-</td><td>5.21e-02</td><td>-</td><td>1.29e-01</td><td>-</td><td>1.61e-01</td><td>-</td></tr><tr><td><eq>2^{3}</eq></td><td>6.37e-03</td><td>2.21</td><td>1.36e-03</td><td>5.26</td><td>6.52e-04</td><td>7.63</td><td>1.47e-03</td><td>6.77</td></tr><tr><td><eq>2^{4}</eq></td><td>2.41e-03</td><td>1.40</td><td>3.46e-04</td><td>1.97</td><td>7.76e-06</td><td>6.39</td><td>2.99e-07</td><td>12.27</td></tr><tr><td><eq>2^{5}</eq></td><td>8.60e-04</td><td>1.49</td><td>7.38e-05</td><td>2.23</td><td>3.97e-07</td><td>4.29</td><td>5.35e-09</td><td>5.80</td></tr><tr><td><eq>2^{6}</eq></td><td>2.96e-04</td><td>1.54</td><td>1.49e-05</td><td>2.31</td><td>1.84e-08</td><td>4.43</td><td>8.46e-11</td><td>5.98</td></tr><tr><td><eq>2^{7}</eq></td><td>9.95e-05</td><td>1.57</td><td>2.91e-06</td><td>2.35</td><td>8.04e-10</td><td>4.51</td><td>1.25e-12</td><td>6.09</td></tr><tr><td><eq>2^{8}</eq></td><td>3.32e-05</td><td>1.58</td><td>5.59e-07</td><td>2.38</td><td>3.41e-11</td><td>4.56</td><td>1.76e-14</td><td>6.14</td></tr><tr><td>Pred</td><td></td><td>1.60</td><td></td><td>2.40</td><td></td><td>4.60</td><td></td><td>6.20</td></tr></table>

From (4.56),(4.57) and (4.58), we can rewrite (4.50) as 

$$
\widehat {Q _ {2}} \leq \left\{\frac {2 ^ {\sigma - 1} \Gamma (\sigma + \lambda - 1 / 2)}{\pi (\sigma - 1) \Gamma (\lambda + 1 / 2)} \frac {N + \lambda + 1}{N - \sigma - k + 2} \frac {\Gamma (N - \sigma + 2)}{\Gamma (N + \sigma)} \right\} \mathbb {U} _ {- 1 +} ^ {\sigma , k}. \tag {4.60}
$$

Then from (4.48), (4.59) and (4.60), we obtain (4.47) immediately. 

Corollary 4.3 For some m, $k \in  { \mathbb { N } } _ { 0 } , s \in ( 0 , 1 )$ and $\sigma : = m + s$ , assume that $u \in \mathbb { W } _ { - 1 + } ^ { \sigma , k } ( \Lambda )$ As $N \to \infty$ , for $0 \leq \lambda < \sigma + k$ and fixed k, we have 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L ^ {\infty}} \leq C \bigl (N ^ {1 + \lambda - \sigma - k} + N ^ {2 (1 - \sigma)} \bigr) \mathbb {U} _ {- 1 +} ^ {\sigma , k},
$$

where C is a generic positive constant independent of N and u. 

Remark 4.4 Considgreater than 0, then $u ( x ) = ( 1 + x ) ^ { \alpha } g ( x )$ $g ( x ) \in C ^ { \infty } ( \Lambda )$ , if α is a non-integerirectly apply Theorem $u \in \mathbb { W } _ { - 1 + } ^ { \alpha + 1 , k } ( \Lambda )$ $k = \lfloor \alpha + \lambda + 1 \rfloor$ 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L ^ {\infty}} \leq C N ^ {- 2 \alpha} U _ {- 1 +} ^ {\alpha + 1, k}, \quad N \gg 1,
$$

where C is a positive constant independent of N and u. 

Here, we present the optimal (weighted) $L ^ { 2 }$ -estimates of the Gegenbauer polynomial approximations with endpoint singularities. 

Theorem 4.6 For some m, $k \in  { \mathbb { N } } _ { 0 } , s \in ( 0 , 1 )$ and $\sigma : = m + s$ , assume that $u \in \mathbb { W } _ { - 1 + } ^ { \sigma , k } ( \Lambda )$ Then for $\sigma + k \le N + 1$ , we have the following bounds. 

(i) $H 0 \leq \lambda < 1 / 2$ , 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Omega)}
$$

$$
\leq \left\{\frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{2 ^ {2 (\lambda + \sigma + k - 1)} \pi \Gamma (N + 2)} \frac {\Gamma^ {2} ((N - \sigma - k) / 2 + 1)}{\Gamma^ {2} ((N + \sigma + k) / 2 + \lambda + 1)} + \mathcal {Q} _ {1} + \mathcal {Q} _ {2} \right\} ^ {1 / 2} \mathbb {U} _ {- 1 +} ^ {\sigma , k}, \tag {4.61}
$$

where 

$$
\mathcal {Q} _ {1} = \frac {2 ^ {2 (1 - \lambda)} (N + \lambda + 2)}{(2 \sigma + 2 k - 1) (N + \sigma + k + 2) \pi} \frac {\Gamma (N + 2 \lambda + 2)}{\Gamma (N + 3)}
$$

$$
\frac {\Gamma^ {2} ((N + \sigma + k) / 2 + 2)}{\Gamma^ {2} ((N + \sigma + k + 1) / 2 + \lambda + 1)} \frac {\Gamma (N - \sigma - k + 2)}{\Gamma (N + \sigma + k + 1)},
$$

$$
\mathcal {Q} _ {2} = \frac {2 ^ {6 \sigma + 4 \lambda - 3} \Gamma^ {2} (\sigma + \lambda - 1 / 2)}{(2 \sigma + \lambda - 2) \pi^ {2}} \frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{(4 N + 6 \lambda + 3) \Gamma (N + 2)}
$$

$$
\frac {\Gamma^ {2} (N - \sigma + 2)}{\Gamma^ {2} (N - \sigma + \lambda + 2)} \left(\frac {N + \lambda + 1}{N - \sigma - k + 2}\right) ^ {2} \frac {\Gamma (2 N - 2 \sigma + 2 \lambda + 3)}{\Gamma (2 N + 2 \sigma + 4 \lambda - 1)}.
$$

(ii) $I f \lambda \geq 1 / 2$ 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Omega)}
$$

$$
\leq \left\{\frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{2 ^ {2 (\lambda + \sigma + k - 1)} \pi \Gamma (N + 2)} \frac {\Gamma^ {2} ((N - \sigma - k) / 2 + 1)}{\Gamma^ {2} ((N + \sigma + k) / 2 + \lambda + 1)} + \widetilde {\mathcal {Q}} _ {1} + \mathcal {Q} _ {2} \right\} ^ {1 / 2} \mathbb {U} _ {- 1 +} ^ {\sigma , k}, \tag {4.62}
$$

where 

$$
\widetilde {\mathcal {Q}} _ {1} = \frac {2 (N + \lambda + 2)}{(2 \sigma + 2 k - 1) (N + \sigma + k + 2) \pi} \frac {\Gamma (N - \sigma - k + 2)}{\Gamma (N + \sigma + k + 1)}.
$$

Proof By the orthogonality of Gegenbauer polynomials, we derive from (4.1), (4.38) and the Cauchy’s inequality, 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Omega)} ^ {2} = \sum_ {n = N + 1} ^ {\infty} \frac {2 ^ {2 \lambda - 1} n ! \Gamma^ {2} (\lambda + 1 / 2)}{(n + \lambda) \Gamma (n + 2 \lambda)} \left| \hat {u} _ {n} ^ {(\lambda)} \right| ^ {2}
$$

$$
\leq \frac {\left(\mathbb {U} _ {- 1 +} ^ {\sigma , k}\right) ^ {2}}{2 ^ {2 (\lambda + \sigma + k - 1)} \pi} \sum_ {n = N + 1} ^ {\infty} \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)}
$$

$$
\frac {\Gamma^ {2} ((n - \sigma - k + 1) / 2)}{\Gamma^ {2} ((n + \sigma + k + 1) / 2 + \lambda)} + \frac {2 ^ {2 (\sigma + \lambda)} \Gamma^ {2} (\sigma + \lambda - 1 / 2)}{\pi^ {2}} \bigl (\mathbb {U} _ {- 1 +} ^ {\sigma , k} \bigr) ^ {2}
$$

$$
\sum_ {n = N + 1} ^ {\infty} \left(\frac {n + \lambda}{n - \sigma - k + 1}\right) ^ {2} \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)} \frac {\Gamma^ {2} (n - \sigma + 1)}{\Gamma^ {2} (n + \sigma + 2 \lambda)}
$$

$$
= \frac {\left(\mathbb {U} _ {- 1 +} ^ {\sigma , k}\right) ^ {2}}{2 ^ {2 (\lambda + \sigma + k - 1)} \pi} \frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{\Gamma (N + 2)}
$$

$$
\frac {\Gamma^ {2} ((N - \sigma - k) / 2 + 1)}{\Gamma^ {2} ((N + \sigma + k) / 2 + \lambda + 1)} + \widetilde {Q} _ {1} + \widetilde {Q} _ {2}, \tag {4.63}
$$

where 

$$
\begin{array}{l} \widetilde {Q} _ {1} := \frac {\left(\mathbb {U} _ {- 1 +} ^ {\sigma , k}\right) ^ {2}}{2 ^ {2 (\lambda + \sigma + k - 1)} \pi} \sum_ {n = N + 2} ^ {\infty} \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)} \frac {\Gamma^ {2} ((n + \sigma + k) / 2 + 1)}{\Gamma^ {2} ((n + \sigma + k + 1) / 2 + \lambda)} \\ \times \frac {\Gamma^ {2} ((n - \sigma - k + 1) / 2)}{\Gamma^ {2} ((n + \sigma + k) / 2 + 1)}, (4.64) \\ \widetilde {Q} _ {2} := \frac {2 ^ {2 (\sigma + \lambda)} \Gamma^ {2} (\sigma + \lambda - 1 / 2)}{\pi^ {2}} \left(\mathbb {U} _ {- 1 +} ^ {\sigma , k}\right) ^ {2} \widehat {\Psi}, (4.65) \\ \end{array}
$$

$$
\widehat {\Psi} = \sum_ {n = N + 1} ^ {\infty} \left(\frac {n + \lambda}{n - \sigma - k + 1}\right) ^ {2} \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)} \frac {\Gamma^ {2} (n - \sigma + 1)}{\Gamma^ {2} (n - \sigma + \lambda + 1)} \frac {\Gamma^ {2} (n - \sigma + \lambda + 1)}{\Gamma^ {2} (n + \sigma + 2 \lambda)}.
$$

For brevity, we denote 

$$
\widetilde {D} _ {n} ^ {(\lambda , \sigma , k)} := \frac {\Gamma (n + 2 \lambda)}{2 ^ {2 \lambda - 1} \Gamma (n + 1)} \frac {\Gamma^ {2} ((n + \sigma + k) / 2 + 1)}{\Gamma^ {2} ((n + \sigma + k + 1) / 2 + \lambda)}, \tag {4.66}
$$

$$
\widehat {D} _ {n} ^ {(\lambda , \sigma)} := \frac {(n + \lambda) \Gamma (n + 2 \lambda)}{\Gamma (n + 1)} \frac {\Gamma^ {2} (n - \sigma + 1)}{\Gamma^ {2} (n - \sigma + \lambda + 1)}, \tag {4.67}
$$

$$
\check {D} _ {n} ^ {(\lambda , \sigma , k)} := \left(\frac {n + \lambda}{n - \sigma - k + 1}\right) ^ {2}. \tag {4.68}
$$

We first consider (4.66). Indeed, with direct calculations similar to those in (4.29), we have 

$$
\begin{array}{l} \frac {\widetilde {D} _ {n} ^ {(\lambda , \sigma , k)}}{\widetilde {D} _ {n + 2} ^ {(\lambda , \sigma , k)}} = \frac {(n + 1) (n + 2) (n + \sigma + k + 2 \lambda + 1) ^ {2}}{(n + 2 \lambda) (n + 2 \lambda + 1) (n + \sigma + k + 2) ^ {2}} \\ = 1 - \frac {(2 \lambda - 1) \widetilde {\mathcal {W}}}{(n + 2 \lambda) (n + 2 \lambda + 1) (n + \sigma + k + 2) ^ {2}}, \\ \end{array}
$$

where 

$$
\widetilde {\mathcal {W}} = 2 \lambda (\sigma + k) ^ {2} + (4 n + 8) (\sigma + k) \lambda + (2 n + 4) \lambda + (2 n + 2) (\sigma + k) ^ {2}
$$

$$
+ (2 n ^ {2} + 6 n + 4) (\sigma + k) + n ^ {2} + 3 n + 2 > 0.
$$

We can get that the sequence $\left\{ \widetilde D _ { n } ^ { ( \lambda , \sigma , k ) } \right\} _ { n \geq N + 2 }$ is decreasing for $- 1 / 2 < \lambda < 1 / 2$ , and n N 2 increasing for $\lambda \geq 1 / 2 .$ ≥ + . From (4.17), we have that 

$$
\lim _ {n \rightarrow \infty} \widetilde {D} _ {n} ^ {(\lambda , \sigma , k)} = 1, \quad \widetilde {D} _ {n} ^ {(\lambda , \sigma , k)} <   \widetilde {D} _ {n + 2} ^ {(\lambda , \sigma , k)} <   \dots \leq 1. \tag {4.69}
$$

Furthermore, we consider (4.67). By simple calculations lead to 

$$
\frac {\widehat {D} _ {n} ^ {(\lambda , \sigma)}}{\widehat {D} _ {n + 1} ^ {(\lambda , \sigma)}} = \frac {(n + 1) (n + \lambda) (n - \sigma + \lambda + 1) ^ {2}}{(n + 2 \lambda) (n + \lambda + 1) (n - \sigma + 1) ^ {2}} = 1 + \frac {\lambda \widehat {\mathcal {W}}}{(n + 2 \lambda) (n + \lambda + 1) (n - \sigma + 1) ^ {2}} > 1, \tag {4.70}
$$

where 

$$
\widehat {\mathcal {W}} = (n + 1) ^ {2} (2 \lambda - 1) + \lambda (n + 1) (\lambda + 2 \sigma + n) - \sigma^ {2} (2 \lambda + 2 n + 1) > 0.
$$

We can get that the sequence $\left\{ { \widehat D } _ { n } ^ { ( \lambda , \sigma ) } \right\} _ { n \geq N + 1 }$ is decreasing for $\lambda \geq 0$ 

≥ +Finally, we consider (4.68). From (4.57) we immediately obtain 

$$
\frac {\check {D} _ {n} ^ {(\lambda , \sigma , k)}}{\check {D} _ {n + 1} ^ {(\lambda , \sigma , k)}} = \left(\frac {(n + \lambda) (n - \sigma - k + 2)}{(n + \lambda + 1) (n - \sigma - k + 1)}\right) ^ {2} > 1, \tag {4.71}
$$

so the sequence $\left\{ \check { D } _ { n } ^ { ( \lambda , \sigma , k ) } \right\} _ { n > N + 1 }$ is decreasing for $\lambda \geq 0$ . n N 1 

According to (4.32), we can easily obtain 

$$
\begin{array}{l} \sum_ {n = N + 2} ^ {\infty} \frac {(n + \lambda) \Gamma^ {2} ((n - \sigma - k + 1) / 2)}{\Gamma^ {2} ((n + \sigma + k) / 2 + 1)} \leq \sum_ {n = N + 2} ^ {\infty} \frac {2 ^ {2 (\sigma + k + 1)} (n + \lambda)}{(n + \sigma + k)} \frac {\Gamma (n - \sigma - k)}{\Gamma (n + \sigma + k)} \\ = \sum_ {n = N + 2} ^ {\infty} \frac {2 ^ {2 (\sigma + k)} (n + \lambda)}{(2 \sigma + 2 k - 1) (n + \sigma + k)} \left\{\frac {\Gamma (n - \sigma - k)}{\Gamma (n + \sigma + k - 1)} - \frac {\Gamma (n - \sigma - k + 1)}{\Gamma (n + \sigma + k)} \right\} \\ = \frac {2 ^ {2 (\sigma + k)} (N + \lambda + 2)}{(2 \sigma + 2 k - 1) (N + \sigma + k + 2)} \frac {\Gamma (N - \sigma - k + 2)}{\Gamma (N + \sigma + k + 1)}. \tag {4.72} \\ \end{array}
$$

Through (4.21) and (4.31), we have 

$$
\begin{array}{l} \frac {\Gamma^ {2} (n - \sigma + \lambda + 1)}{\Gamma^ {2} (n + \sigma + 2 \lambda)} = \frac {1}{n + \sigma + 2 \lambda - 1} \frac {\Gamma^ {2} (n - \sigma + \lambda + 1)}{\Gamma (n + \sigma + 2 \lambda - 1) \Gamma (n + \sigma + 2 \lambda)} \\ \leq \frac {1}{n + \sigma + 2 \lambda - 1} \frac {\Gamma (n - \sigma + \lambda + 1 / 2) \Gamma (n - \sigma + \lambda + 1)}{\Gamma (n + \sigma + 2 \lambda - 1) \Gamma (n + \sigma + 2 \lambda - 1 / 2)} \\ = \frac {2 ^ {4 \sigma + 2 \lambda - 3}}{n + \sigma + 2 \lambda - 1} \frac {\Gamma (2 n - 2 \sigma + 2 \lambda + 1)}{\Gamma (2 n + 2 \sigma + 4 \lambda - 2)} \\ \leq \frac {2 ^ {4 \sigma + 2 \lambda - 3}}{(2 \sigma + \lambda - 2) (4 n + 6 \lambda - 1)} \left\{\frac {\Gamma (2 n - 2 \sigma + 2 \lambda + 1)}{\Gamma (2 n + 2 \sigma + 4 \lambda - 3)} - \frac {\Gamma (2 n - 2 \sigma + 2 \lambda + 3)}{\Gamma (2 n + 2 \sigma + 4 \lambda - 1)} \right\}. \tag {4.73} \\ \end{array}
$$

From (4.66) and (4.72), we can rewrite (4.64) for $- 1 / 2 < \lambda < 1 / 2$ as follows: 

$$
\widetilde {Q} _ {1} \leq \frac {\left(\mathbb {U} _ {- 1 +} ^ {\sigma , k}\right) ^ {2}}{2 ^ {2 (\lambda - 1)} \pi} \frac {(N + \lambda + 2)}{(2 \sigma + 2 k - 1) (N + \sigma + k + 2)} \frac {\Gamma (N + 2 \lambda + 2)}{\Gamma (N + 3)} \tag {4.74}
$$

$$
\times \frac {\Gamma^ {2} ((N + \sigma + k) / 2 + 2)}{\Gamma^ {2} ((N + \sigma + k + 1) / 2 + \lambda + 1)} \frac {\Gamma (N - \sigma - k + 2)}{\Gamma (N + \sigma + k + 1)}.
$$

Similarly, for $\lambda \geq 1 / 2$ , we can rewrite (4.64) using (4.69) and (4.72) as follows: 

$$
\widetilde {Q} _ {1} \leq \frac {2 (N + \lambda + 2)}{(2 \sigma + 2 k - 1) (N + \sigma + k + 2) \pi} \frac {\Gamma (N - \sigma - k + 2)}{\Gamma (N + \sigma + k + 1)} \left(\mathbb {U} _ {- 1 +} ^ {\sigma , k}\right) ^ {2}. \tag {4.75}
$$


Table 8 Order of u(x ) − π (λ)N u(x )  L2ω (- $\left\| u ( x ) - \pi _ { N } ^ { ( \lambda ) } u ( x ) \right\| _ { L _ { \omega _ { \lambda } } ^ { 2 } ( \Omega ) }$ with $\lambda = 0 . 3$ and $u = ( 1 + x ) ^ { \alpha } e ^ { \sin { y } }$


<table><tr><td>N</td><td><eq>\alpha = 0.8</eq></td><td>Order</td><td><eq>\alpha = 1.2</eq></td><td>Order</td><td><eq>\alpha = 2.3</eq></td><td>Order</td></tr><tr><td><eq>2^{2}</eq></td><td>7.33e-03</td><td>-</td><td>1.38e-02</td><td>-</td><td>2.84e-02</td><td>-</td></tr><tr><td><eq>2^{3}</eq></td><td>4.97e-04</td><td>3.88</td><td>1.17e-04</td><td>6.88</td><td>1.06e-04</td><td>8.06</td></tr><tr><td><eq>2^{4}</eq></td><td>1.06e-04</td><td>2.23</td><td>1.67e-05</td><td>2.81</td><td>4.08e-07</td><td>8.03</td></tr><tr><td><eq>2^{5}</eq></td><td>2.13e-05</td><td>2.31</td><td>1.98e-06</td><td>3.08</td><td>1.11e-08</td><td>5.20</td></tr><tr><td><eq>2^{6}</eq></td><td>4.17e-06</td><td>2.36</td><td>2.24e-07</td><td>3.14</td><td>2.82e-10</td><td>5.30</td></tr><tr><td><eq>2^{7}</eq></td><td>8.01e-07</td><td>2.38</td><td>2.49e-08</td><td>3.17</td><td>6.90e-12</td><td>5.35</td></tr><tr><td><eq>2^{8}</eq></td><td>1.51e-07</td><td>2.41</td><td>2.73e-09</td><td>3.19</td><td>1.66e-13</td><td>5.37</td></tr><tr><td>Pred</td><td></td><td>2.40</td><td></td><td>3.20</td><td></td><td>5.40</td></tr></table>

From (4.70), (4.71) and (4.73), we can rewrite (4.65) as 

$$
\begin{array}{l} \widetilde {Q} _ {2} \leq \frac {2 ^ {6 \sigma + 4 \lambda - 3} \Gamma^ {2} (\sigma + \lambda - 1 / 2)}{(2 \sigma + \lambda - 2) \pi^ {2}} \frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{\Gamma (N + 2)} \frac {\Gamma^ {2} (N - \sigma + 2)}{\Gamma^ {2} (N - \sigma + \lambda + 2)} \Big (\mathbb {U} _ {- 1 +} ^ {\sigma , k} \Big) ^ {2} \\ \times \left(\frac {N + \lambda + 1}{N - \sigma - k + 2}\right) ^ {2} \sum_ {n = N + 1} ^ {\infty} \frac {1}{(4 n + 6 \lambda - 1)} \left\{\frac {\Gamma (2 n - 2 \sigma + 2 \lambda + 1)}{\Gamma (2 n + 2 \sigma + 4 \lambda - 3)} - \frac {\Gamma (2 n - 2 \sigma + 2 \lambda + 3)}{\Gamma (2 n + 2 \sigma + 4 \lambda - 1)} \right\} \\ = \frac {2 ^ {6 \sigma + 4 \lambda - 3} \Gamma^ {2} (\sigma + \lambda - 1 / 2)}{(2 \sigma + \lambda - 2) \pi^ {2}} \frac {(N + \lambda + 1) \Gamma (N + 2 \lambda + 1)}{(4 N + 6 \lambda + 3) \Gamma (N + 2)} \frac {\Gamma^ {2} (N - \sigma + 2)}{\Gamma^ {2} (N - \sigma + \lambda + 2)} \\ \times \left(\frac {N + \lambda + 1}{N - \sigma - k + 2}\right) ^ {2} \frac {\Gamma (2 N - 2 \sigma + 2 \lambda + 3)}{\Gamma (2 N + 2 \sigma + 4 \lambda - 1)} \left(\mathbb {U} _ {- 1 +} ^ {\sigma , k}\right) ^ {2}. \tag {4.76} \\ \end{array}
$$

We get (4.61) from (4.63), (4.74) and (4.76). Also, we get (4.62) from (4.63), (4.75) and (4.76). This completes the proof. □ 

Co 4 Fo m , $k \in  { \mathbb { N } } _ { 0 } , s \in ( 0 , 1 )$ and $\sigma : = m + s$ , assume that $u \in \mathbb { W } _ { - 1 + } ^ { \sigma , k } ( \Lambda )$ $N \to \infty$ $\lambda \geq 0$ 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Omega)} \leq C \big (N ^ {1 / 2 - \sigma - k} + N ^ {3 / 2 - 2 \sigma - \lambda} \big) \mathbb {U} _ {- 1 +} ^ {\sigma , k},
$$

where c is a generic positive constant independent of N and u. 

Remark 4.5 Considgreater than 0, then $u ( x ) = ( 1 + x ) ^ { \alpha } g ( x )$ $g ( x ) \in C ^ { \infty } ( \Lambda )$ , if α is a non-integerirectly apply Theorem $u \in \mathbb { W } _ { - 1 + } ^ { \alpha + 1 , k } ( \Lambda )$ $k = \lfloor \alpha + \lambda + 1 \rfloor$ 4.6 and (4.17) to get 

$$
\left\| u (x) - \pi_ {N} ^ {(\lambda)} u (x) \right\| _ {L _ {\omega_ {\lambda}} ^ {2} (\Omega)} \leq C N ^ {- 2 \alpha - \lambda - 1 / 2} \mathbb {U} _ {- 1 +} ^ {\alpha + 1, k}, \quad N \gg 1,
$$

where C is a positive constant independent of N and u. 

Remark 4.6 Xiang [36] considered functions that belong to the fractional Sobolev space but exhibit logarithmic singularities, such as 

$$
f (x) = | x - \theta | ^ {\gamma} \ln^ {\mu} | x - \theta | g (x) \in \mathbb {W} _ {\theta} ^ {\gamma + 1 - \epsilon} (\Omega) \forall \epsilon \in (0, 1),
$$


Table 9 Order of u(x ) − π (λ)N u(x )  L2ω (- $\left\| u ( x ) - \pi _ { N } ^ { ( \lambda ) } u ( x ) \right\| _ { L _ { \omega _ { \lambda } } ^ { 2 } ( \Omega ) }$ with $\lambda = 0 . 6$ and $\begin{array} { r } { u = ( 1 + x ) ^ { \alpha } e ^ { \sin . } } \end{array}$ x


<table><tr><td>N</td><td><eq>\alpha = 0.8</eq></td><td>Order</td><td><eq>\alpha = 1.2</eq></td><td>Order</td><td><eq>\alpha = 2.3</eq></td><td>Order</td></tr><tr><td><eq>2^{2}</eq></td><td>6.15e-03</td><td>-</td><td>1.12e-02</td><td>-</td><td>2.29e-02</td><td>-</td></tr><tr><td><eq>2^{3}</eq></td><td>2.97e-04</td><td>4.37</td><td>7.27e-05</td><td>7.27</td><td>8.59e-05</td><td>8.06</td></tr><tr><td><eq>2^{4}</eq></td><td>5.41e-05</td><td>2.46</td><td>9.26e-06</td><td>2.97</td><td>2.59e-07</td><td>8.37</td></tr><tr><td><eq>2^{5}</eq></td><td>9.13e-06</td><td>2.57</td><td>9.23e-07</td><td>3.33</td><td>6.02e-09</td><td>5.43</td></tr><tr><td><eq>2^{6}</eq></td><td>1.47e-06</td><td>2.63</td><td>8.66e-08</td><td>3.41</td><td>1.28e-10</td><td>5.56</td></tr><tr><td><eq>2^{7}</eq></td><td>2.32e-07</td><td>2.67</td><td>7.89e-09</td><td>3.46</td><td>2.58e-12</td><td>5.63</td></tr><tr><td><eq>2^{8}</eq></td><td>3.58e-08</td><td>2.69</td><td>7.07e-10</td><td>3.48</td><td>5.08e-14</td><td>5.66</td></tr><tr><td>Pred</td><td></td><td>2.70</td><td></td><td>3.50</td><td></td><td>5.70</td></tr></table>

where $g \in C ^ { \infty } [ - 1 , 1 ] , \ \gamma > 0$ and $\theta \in ( - 1 , 1 )$ . However, the estimate is still suboptimal in general sense. How to modify the fractional space to best characterize this type singularity appears nontrivial and is still open [18, Remark 4.5]. We will report on the work related to optimal bounds for logarithmic singularities in future publications. 

# 5 Concluding Remarks

Along the way to Liu et al. [18, 19, 39], we have conducted a rigorous error analysis for functions in fractional spaces. We present sharp bounds for Gegenbauer expansion coefficients and derive optimal (weighted) $L ^ { \infty } -$ and $L ^ { 2 } .$ -estimates for the error of Gegenbauer approximations. Illustrative numerical results are provided to support the analysis. Here, we merely discussed the approximation results, but hopefully the analysis techniques can be extended to general Jacobi polynomial approximations and enrich the existing theory for $p$ and $h p$ methods for singular problems in multiple dimensions. 

Funding The research of the second author was supported by the National Natural Science Foundation of China (No. 12271128) and the Natural Science Foundation of Heilongjiang Province of China (No. YQ2023A002). The research of the third author was supported by the National Natural Science Foundation of China (No. 12371367) and by the Hubei Provincial Natural Science Foundation of China (No. 2023AFA083). The research of the fourth author was supported by the National Natural Science Foundation of China (No. 11971131) and the Natural Science Foundation of Heilongjiang Province of China (No. ZD2022A001). 

Data Availibility The code used in this work will be made available upon request to the authors. 

# Declarations

Conflict of interest The authors declare that they have no conflict of interest. 

# References



1. Alzer, H.: On some inequalities for the Gamma and Psi functions. Math. Comput. 66(217), 373–389 (1997) 





2. Andrews, G.E., Askey, R., Roy, R.: Special Functions, Encyclopedia of Mathematics and its Applications, vol. 71. Cambridge University Press, Cambridge (1999) 





3. Babuška, I., Whiteman, J., Strouboulis, T.: Finite Elements: An Introduction to the Method and Error Estimation. Oxford University Press, Oxford (2010) 





4. Babuška, I., Guo, B.Q.: Optimal estimates for lower and upper bounds of approximation errors in the p-version of the finite element method in two dimensions. Numer. Math. 85, 219–255 (2000) 





5. Babuška, I., Guo, B.Q.: Direct and inverse approximation theorems for the p-version of the finite element method in the framework of weighted Besov spaces . I. Approximability of functions in the weighted Besov spaces. SIAM J. Numer. Anal. 39(5), 1512–1538 (2001). https://doi.org/10.1137/S0036142901356551 





6. Babuška, I., Suri, M.: The optimal convergence rate of the p-version of the finite element method. SIAM J. Numer. Anal. 24(4), 750–776 (1987) 





7. Brezis, H.: Functional Analysis Sobolev Spaces and Partial Differential Equations. Springer, New York (2011) 





8. Buttazzo, G., Giaquinta, M., Hildebrandt, S.: One-Dimensional Variational Problems: An Introduction. Oxford University Press, New York (1998) 





9. Canuto, C., Tabacco, A.: Mathematical Analysis I, 2nd edn. Springer, Switzerland (2015) 





10. Chen, S., Shen, J.: Log orthogonal functions: approximation properties and applications. IMA J. Numer. Anal. 42(1), 712–743 (2022). https://doi.org/10.1093/imanum/draa087 





11. Chen, S., Shen, J.: Log orthogonal functions in semi-infinite intervals: approximation results and applications. SIAM J. Numer. Anal. 61(1), 110–134 (2023). https://doi.org/10.1137/21M1466840 





12. Chen, S., Shen, J., Wang, L.-L.: Generalized Jacobi functions and their applications to fractional differential equations. Math. Comp. 85(300), 1603–1638 (2016). https://doi.org/10.1090/mcom3035 





13. Guo, B.Y., Wang, L.-L.: Jacobi approximations in non-uniformly Jacobi-weighted Sobolev spaces. J. Approx. Theory 128(1), 1–41 (2004) 





14. Houston, P., Schwab, C., Süli, E.: Discontinuous hp-finite element methods for advection-diffusionreaction problems. SIAM J. Numer. Anal. 39(6), 2133–2163 (2002). https://doi.org/10.1137/ S0036142900374111 





15. Lang, S.: Real and Functional Analysis, 3rd edn. Springer, New York (1993) 





16. Li, H.G.: Elliptic Equations with Singularities: A Priori Analysis and Numerical Approaches. Ph.D. thesis (2008) 





17. Liu, W.J., Wang, L.-L.: Asymptotics of the generalized Gegenbauer functions of fractional degree. J. Approx. Theory 253, 105378 (2020) 





18. Liu, W.J., Wang, L.-L., Li, H.Y.: Optimal error estimates for Chebyshev approximations of functions with limited regularity in fractional Sobolev-type spaces. Math. Comp. 88(320), 2857–2895 (2019). https:// doi.org/10.1090/mcom/3456 





19. Liu, W.J., Wang, L.-L., Wu, B.Y.: Optimal error estimates for Legendre expansions of singular functions with fractional derivatives of bounded variation. Adv. Comput. Math. 47(79), 1–32 (2021). https://doi. org/10.1007/s10444-021-09905-3 





20. Melenk, J.: hp-Finite Element Methods for Singular Perturbations. Springer, Heidelberg (2002) 





21. Müller, F., Schötzau, D., Schwab, C.: Symmetric interior penalty discontinuous Galerkin methods for elliptic problems in polygons. SIAM J. Numer. Anal. 55(5), 2490–2521 (2017). https://doi.org/10.1137/ 17M1120634 





22. Olver, F.W.J., Lozier, D.W., Boisvert, R.F., Clark, C.W.: NIST Handbook of Mathematical Functions. Cambridge University Press, New York (2010) 





23. Samko, S.G., Kilbas, A.A., Marichev, O.I.: Fractional Integrals and Derivatives Theory and Applications. Gordan and Breach Science Publisher, New York (1993) 





24. Shen, J., Tang, T., Wang, L.-L.: Spectral Methods: Algorithms Analysis and Applications. Springer, New York (2011) 





25. Stein, E.M., Shakarchi, R.: Real Analysis: Measure Theory, Integration, and Hilbert Spaces, Princeton Lectures in Analysis, vol. III. Princeton University Press, New Jersey (2005) 





26. Szegö, G.: Orthogonal Polynomials, 4th edn. American Mathematical Society, Providence (1975) 





27. Trefethen, L.N.: Is Gauss quadrature better than Clenshaw-Curtis? SIAM Rev. 50(1), 67–87 (2008) 





28. Wahlbin, L.B.: On the sharpness of certain local estimates for $\mathring { H } ^ { 1 }$ projections into finite element spaces: influence of a re-entrant corner. Math. Comp. 42(165), 1–8 (1984). https://doi.org/10.2307/2007556 





29. Wang, H.Y.: On the optimal estimates and comparison of Gegenbauer expansion coefficients. SIAM J. Numer. Anal. 54(3), 1557–1581 (2016) 





30. Wang, H.Y.: How much faster does the best polynomial approximation converge than Legendre projection? Numer. Math. 147(2), 481–503 (2021). https://doi.org/10.1007/s00211-021-01173-z 





31. Wang, H.Y.: Analysis of error localization of Chebyshev spectral approximations. SIAM J. Numer. Anal. 61(2), 952–972 (2023) 





32. Wang, H.Y.: New error bounds for Legendre approximations of differentiable functions. J. Fourier Anal. Appl. (2023). https://doi.org/10.1007/s00041-023-10024-4 





33. Wang, H.Y.: Optimal rates of convergence and error localization of Gegenbauer projections. IMA J. Numer. Anal. 43(4), 2413–2444 (2023). https://doi.org/10.1093/imanum/drac047 





34. Wang, L.-L., Zhao, X.D., Zhang, Z.M.: Superconvergence of Jacobi-Gauss-type spectral interpolation. J. Sci. Comput. 59(3), 667–687 (2014). https://doi.org/10.1007/s10915-013-9777-x 





35. Xiang, S.H.: On interpolation approximation: convergence rates for polynomial interpolation for functions of limited regularity. SIAM J. Numer. Anal. 54(4), 2081–2113 (2016) 





36. Xiang, S.H.: Convergence rates on spectral orthogonal projection approximation for functions of algebraic and logarithmatic regularities. SIAM J. Numer. Anal. 59(3), 1374–1398 (2021) 





37. Xiang, S.H., Kong, D.S., Liu, G.D., Wang, L.-L.: Pointwise error estimates and local superconvergence of Jacobi expansions. Math. Comp. 92(342), 1747–1778 (2023). https://doi.org/10.1090/mcom/3835 





38. Xiang, S.H., Liu, G.D.: Optimal decay rates on the asymptotics of orthogonal polynomial expansions for functions of limited regularities. Numer. Math. 145(1), 117–148 (2020). https://doi.org/10.1007/s00211- 020-01113-3 





39. Xie, R.Y., Wu, B.Y., Liu, W.J.: Optimal error estimates for Chebyshev approximations of functions with endpoint singularities in fractional spaces. J. Sci. Comput. 96(3), 71 (2023) 





40. Xie, Z.Q., Wang, L.-L., Zhao, X.D.: On exponential convergence of Gegenbauer interpolation and spectral differentiation. Math. Comp. 82(282), 1017–1036 (2013) 





41. Zhang, X.L., Boyd, J.: Asymptotic coefficients and errors for Chebyshev polynomial approximations with weak endpoint singularities: effects of different bases. Sci. China Math. 65 (2022) 





42. Zhao, X.D., Wang, L.-L., Xie, Z.Q.: Sharp error bounds for Jacobi expansions and Gegenbauer-Gauss quadrature of analytic functions. SIAM J. Numer. Anal. 51(3), 1443–1469 (2013). https://doi.org/10. 1137/12089421X 



Publisher’s Note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations. 

Springer Nature or its licensor (e.g. a society or other partner) holds exclusive rights to this article under a publishing agreement with the author(s) or other rightsholder(s); author self-archiving of the accepted manuscript version of this article is solely governed by the terms of such publishing agreement and applicable law. 