# Log orthogonal functions: approximation properties and applications

Sheng Chen 

Beijing Computational Science Research Center, Beijing 100193, P.R. China; 

School of Mathematics and Statistics, Jiangsu Normal University, Xuzhou 221116, China 

and 

Jie Shen∗ 

Department of Mathematics, Purdue University, West Lafayette, IN 47907, USA 

∗Corresponding author: shen7@purdue.edu. 

[Received on 2 March 2020; revised on 20 July 2020] 

We present two new classes of orthogonal functions, log orthogonal functions and generalized log orthogonal functions, which are constructed by applying a log mapping to Laguerre polynomials. We develop basic approximation theory for these new orthogonal functions, and apply them to solve several typical fractional differential equations whose solutions exhibit weak singularities. Our error analysis and numerical results show that our methods based on the new orthogonal functions are particularly suitable for functions that have weak singularities at one endpoint and can lead to exponential convergence rate, as opposed to low algebraic rates if usual orthogonal polynomials are used. 

Keywords: log orthogonal functions; Laguerre functions; mapped spectral methods; fractional differential equations; weak singularity. 

# 1. Introduction

It is well known that classical spectral methods can provide high accuracy for problems with smooth solutions Gottlieb & Orszag (1977); Boyd (2001); Canuto et al. (2006); Shen et al. (2011), but may not have any advantage for problems with non-smooth solutions. On the other hand, solutions of many practical applications involve weakly singular solutions, such as in non-smooth domains, with non-matching boundary conditions, in integral equations with singular/weakly singular kernels, and in fractional differential equations (FDEs). One effective strategy in finite differences/finite elements is to employ a local adaptive procedure (Morin et al., 2002), but this strategy cannot be effectively extended to the global spectral method. Hence, in order to develop accurate spectral methods for problems with non-smooth solutions, one has to choose suitable basis functions that can effectively approximate the underlying non-smooth solutions. A popular strategy in this regard is to enrich the usual polynomial-based approximation space by adding special functions that capture the singular behavior of the underlying problem, for example, the so-called singular functions method (Strang & Fix, 1973), extended or generalized finite element method (GFEM/XFEM) (cf. Fries & Belytschko, 2010; Babuška & Banerjee, 2012, and the references therein), and in the context of spectral methods, the enriched spectral-tau method (Schumack et al., 1991) and the enriched spectral-Galerkin method (Chen & Shen, 2018, 2020; Chen, 2020). Another effective strategy in the context of spectral method is to construct special orthogonal functions that are suitable for a certain class of problems with singular behaviors. In addition to classical orthogonal polynomials one can use suitable mappings to classical orthogonal polynomials to construct orthogonal functions in weighted Sobolev spaces, leading to the so-called mapped spectral methods (cf. Boyd, 2001, for an extensive discussion). In Boyd (1987a,b); Shen & Wang (2004); Wang & Shen (2005) the mapped spectral methods have been successfully used in constructing efficient spectral methods for problems in unbounded domains, and in Müntz Galerkin method (Shen & Wang, 2016) for a special class of singular problems. On the other hand, Boyd (1986) briefly discussed several possible alternatives to deal with weak singularities at both endpoints through different mappings. 

In this paper we are concerned with problems that exhibit weakly singular behaviors at the initial time for initial value problems (IVPs), or at one endpoint for boundary value problems (BVPs). We construct special classes of orthogonal functions, through a suitable log mapping to Laguerre functions, which are capable of resolving weak singularities. We shall develop basic approximation results for two new classes of orthogonal functions, log orthogonal functions (LOFs) and generalized log orthogonal functions (GLOFs). In particular, these results indicate that approximation by the LOFs and GLOFs to functions behaving like $t ^ { r } ( - \log t ) ^ { k }$ near t = 0 will converge exponentially for any real $r \geq 0 , k \in  { \mathbb { N } } _ { 0 }$ . In fact, we believe that this is the first set of basis functions that can approximate regular polynomials $t ^ { k } ( k \in \mathbb { N } _ { 0 } )$ and weakly singular functions like $t ^ { r } ( - \log t ) ^ { k } ~ ( 0 ~ < ~ r ~ < ~ 1 , ~ k ~ \in ~ \mathbb { N } _ { 0 } )$ with exponential accuracy. Thus, LOFs and GLOFs are particularly suitable for problems whose solutions exhibit weak singularities behaving like $\textstyle { \sum _ { i } c _ { i } t ^ { r _ { i } } ( - \log t ) ^ { k _ { i } } }$ near t = 0 for small $r _ { i } > 0$ . In particular, solutions of time FDEs and BVPs with one-sided fractional derivatives fall into this category. Hence, the spectral methods using GLOFs that we propose in this paper can be used to deal with a large class of FDEs having weak singularities at the initial time or at one endpoint. 

Numerical solution of FDEs has been a subject of intensive investigation in recent years, cf. for instance Meerschaert & Tadjeran (2004); Sun & Wu (2006); Ervin et al. (2007); Ervin & Roop (2007); Jin et al. (2013) (and the references therein) for finite-difference and finite-element methods, and Li & Xu (2009, 2010); Li et al. (2012); Zayernouri & Karniadakis (2013); Chen et al. (2016) (and the references therein) for spectral methods. However, most of the error analysis are derived in the context of usual Sobolev spaces that are not quite suitable for FDEs. In our previous works (Chen et al., 2016; Mao et al., 2016) we developed an error analysis using the generalized Jacobi functions (GJFs) based on the non-uniformly weighted spaces that showed that, for some model FDEs whose solutions behave as $( t - a ) ^ { r } ( b - t ) ^ { s } h ( t )$ with known $r , s \ > \ 0$ and smooth h(t), the error may converge exponentially as long as the data function is smooth in the usual sense. However, for more general FDEs such as those with variable coefficients or nonlinearity, the singular behavior of their solutions is unknown a priori, so approximations by GJFs cannot achieve desired accuracy. However, GLOFs can handle functions with unknown endpoint singularity since they can approximate singular functions of the form $\textstyle { \sum _ { i } c _ { i } t ^ { r _ { i } } ( - \log t ) ^ { k _ { i } } }$ with exponential accuracy. In fact, we show in Corollary 3.1 that GLOFs can achieve exponential convergence for typical singular functions $t ^ { r } ( - \log t ) ^ { k }$ for any $r \geq 0 , k \in  { \mathbb { N } } _ { 0 }$ . To the best of our knowledge this is the first set of basis functions that can approximate the singular solutions of the form in (4.6) with exponential accuracy. 

The rest of the paper is organized as follows. In the next section we introduce the LOFs, derive optimal projection and interpolation errors in weighted pseudo-derivatives that are adapted to the involved mapping. In Section 3 we introduce the GLOFs that involve an additional parameter, so are more flexible than LOFs, and derive the corresponding optimal projection and interpolation errors. In Section 4 we apply GLOFs to solve several typical classes of FDEs, and derive optimal error estimates that indicate, in particular, that for solutions and data functions having weak singularities at t = 0 or one endpoint, errors of the proposed GLOF-Galerkin methods will converge exponentially. In each of the sections 2, 3 and 4 we also present numerical results to validate the theoretical estimates and to show the effectiveness of our new method. Some concluding remarks are given in the last section. 

# 2. Log orthogonal functions

In this section we introduce the LOFs, derive the corresponding approximation theory, and present numerical results to validate the theoretical estimates and to show their effectiveness. 

To fix the idea we consider the canonical time interval I = (0, 1). Throughout the paper we use the mapping 

$$
y (t) := - (\beta + 1) \log t, \quad t \in I, \tag {2.1}
$$

to map $I \tan \mathbb { R } ^ { + }$ . We shall make use of $\mathcal { L } _ { n } ^ { ( \alpha ) } ( y ) , \alpha > - 1$ , which is the Laguerre polynomial of $y \in \mathbb { R } ^ { + }$ , satisfying 

$$
\int_ {0} ^ {\infty} \mathscr {L} _ {n} ^ {(\alpha)} (y) \mathscr {L} _ {m} ^ {(\alpha)} (y) y ^ {\alpha} e ^ {- y} d y = \gamma_ {n} ^ {(\alpha)} \delta_ {m n}, \quad \gamma_ {n} ^ {(\alpha)} = \frac {\Gamma (n + \alpha + 1)}{\Gamma (n + 1)}. \tag {2.2}
$$

Some additional properties of Laguerre polynomials are listed in Appendix A, see also Szego (1975); Shen et al. (2011). 

# 2.1 Definition and properties

Definition 2.1 (LOFs). Let α, $\beta > - 1$ . We define the LOFs by 

$$
\mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) = \mathcal {L} _ {n} ^ {(\alpha)} (y (t)) = \mathcal {L} _ {n} ^ {(\alpha)} (- (\beta + 1) \log t), \quad n = 0, 1, \dots . \tag {2.3}
$$

From the properties of Laguerre polynomials listed in Appendix A and the following relations 

$$
y = - (\beta + 1) \log t, \quad \mathrm{d} y = - (\beta + 1) t ^ {- 1} \mathrm{d} t, \quad \partial_ {t} = - (\beta + 1) t ^ {- 1} \partial_ {y} \tag {2.4}
$$

we can easily derive the following lemma: 

Lemma 2.2 The LOFs satisfy the following properties: 

P1. Three-term recurrence relation: 

$$
\begin{array}{l} \mathcal {S} _ {0} ^ {(\alpha , \beta)} (t) = 1, \quad \mathcal {S} _ {1} ^ {(\alpha , \beta)} (t) = (\beta + 1) \log t + \alpha + 1, \\ \text {(2.5)} \end{array}
$$

$$
\mathcal {S} _ {n + 1} ^ {(\alpha , \beta)} (t) = \frac {2 n + \alpha + 1 + (\beta + 1) \log t}{n + 1} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) - \frac {n + \alpha}{n + 1} \mathcal {S} _ {n - 1} ^ {(\alpha , \beta)} (t).
$$

P2. Derivative relations: 

$$
(\beta + 1) ^ {- 1} t \partial_ {t} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) = \mathcal {S} _ {n - 1} ^ {(\alpha + 1, \beta)} (t) = \sum_ {l = 0} ^ {n - 1} \mathcal {S} _ {l} ^ {(\alpha , \beta)} (t), \quad n \geq 1. \tag {2.6}
$$

# P3. Orthogonality:

$$
\int_ {0} ^ {1} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) \mathcal {S} _ {m} ^ {(\alpha , \beta)} (t) (- \log t) ^ {\alpha} t ^ {\beta} \mathrm{d} t = \gamma_ {n} ^ {(\alpha , \beta)} \delta_ {m n}, \tag {2.7}
$$

where 

$$
\gamma_ {n} ^ {(\alpha , \beta)} := \frac {\Gamma (n + \alpha + 1)}{(\beta + 1) ^ {\alpha + 1} \Gamma (n + 1)}.
$$

# P4. Sturm–Liouville problem:

$$
(- \log t) ^ {- \alpha} t ^ {- \beta} \partial_ {t} \left((- \log t) ^ {\alpha + 1} t ^ {\beta + 2} \partial_ {t} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t)\right) + n (\beta + 1) \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) = 0. \tag {2.8}
$$

# P5. Gauss-LOFs quadrature: Let $\{ y _ { j } ^ { ( \alpha ) } , \omega _ { j } ^ { ( \alpha ) } \} _ { j = 0 } ^ { N }$ be the Gauss-node and -weight of $\mathcal { L } _ { n + 1 } ^ { ( \alpha ) } ( y )$ . Denote

$$
\left\{t _ {j} ^ {(\alpha , \beta)} := e ^ {- (\beta + 1) ^ {- 1} y _ {j} ^ {(\alpha)}}, \quad \chi_ {j} ^ {(\alpha , \beta)} = (\beta + 1) ^ {- \alpha - 1} \omega_ {j} ^ {(\alpha)} \right\} _ {j = 0} ^ {N}. \tag {2.9}
$$

Then, 

$$
\int_ {0} ^ {1} p (t) (- \log t) ^ {\alpha} t ^ {\beta} \mathrm{d} t = \sum_ {j = 0} ^ {N} p (t _ {j} ^ {(\alpha , \beta)})   \chi_ {j} ^ {(\alpha , \beta)}, \quad \forall   p \in \mathcal {P} _ {2 N + 1} ^ {\log t}, \tag {2.10}
$$

where 

$$
\mathcal {P} _ {K} ^ {\log t} := \text { span } \bigl \{1, \log t, (\log t) ^ {2}, \dots , (\log t) ^ {K} \bigr \}. \tag {2.11}
$$

Proof. The three-term recurrence relation (2.5) is a straightforward result from (A.1) with the variable transform (2.4). 

(2.6) can be obtained from the relations (A.5) and (2.4). Indeed, 

$$
(\beta + 1) ^ {- 1} t \partial_ {t} \{\mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) \} \stackrel {(2. 4)} {=} - \partial_ {y} \mathcal {L} _ {n} ^ {(\alpha)} (y) \stackrel {(L. 5)} {=} \mathcal {L} _ {n - 1} ^ {(\alpha + 1)} (y) = \mathcal {S} _ {n - 1} ^ {(\alpha + 1, \beta)} (t).
$$

We derive from $y = - ( \beta + 1 )$ log t that 

$$
\int_ {0} ^ {1} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) \mathcal {S} _ {m} ^ {(\alpha , \beta)} (t) (- \log t) ^ {\alpha} t ^ {\beta} d t = \frac {1}{(\beta + 1) ^ {\alpha + 1}} \int_ {0} ^ {\infty} \mathcal {L} _ {n} ^ {(\alpha)} (y) \mathcal {L} _ {m} ^ {(\alpha)} (y) y ^ {\alpha} e ^ {- y} d y.
$$

Hence, we have (2.7). 

(2.8) is valid since 

$$
y = - (\beta + 1) \log t, \quad \partial_ {y} = - (\beta + 1) ^ {- 1} t \partial_ {t}, \quad \mathcal {L} _ {n} ^ {(\alpha)} (y) = \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t)
$$

lead to 

$$
y ^ {- \alpha} e ^ {y} \partial_ {y} \left(y ^ {\alpha + 1} e ^ {- y} \partial_ {y} \mathcal {L} _ {n} ^ {(\alpha)} (y)\right) = \frac {(- \log t) ^ {- \alpha} t ^ {- \beta}}{\beta + 1} \partial_ {t} \left((- \log t) ^ {\alpha + 1} t ^ {\beta + 2} \partial_ {t} \mathcal {S} _ {n} ^ {(\alpha , \beta)}\right).
$$

Finally, setting $t = e ^ { - ( \beta + 1 ) ^ { - 1 } y }$ , we can obtain (2.9) from the Laguerre–Gauss quadrature: 

$$
\begin{array}{l} \int_ {0} ^ {1} p (t) (- \log t) ^ {\alpha} t ^ {\beta} \mathrm{d} t = \int_ {\mathbb {R} ^ {+}} p (x (y)) \frac {y ^ {\alpha} e ^ {- y}}{(\beta + 1) ^ {\alpha + 1}} \mathrm{d} y \\ = \sum_ {j = 0} ^ {N} p (x (y _ {j} ^ {(\alpha)})) \frac {\omega_ {j} ^ {(\alpha)}}{(\beta + 1) ^ {\alpha + 1}} = \sum_ {j = 0} ^ {N} p (t _ {j} ^ {(\alpha , \beta)}) \chi_ {j} ^ {(\alpha , \beta)}. \\ \end{array}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/1f1b1d0f23a5e67e33d98f462d13e14c132319237e6c7ca2b3cd9d293c452f65.jpg)


Remark 2.3 We used two parameters α and $\beta$ to provide better flexibility, $\mathrm { e . g . }$ , they allow us to effectively deal with problems with weight $( - \log t ) ^ { \alpha } t ^ { \beta }$ . In the applications considered in this paper the log term does not appear so we can take $\alpha = 0$ . On the other hand, taking $\beta = 0$ offers good approximation properties for problems with weight $( - \log t ) ^ { \alpha }$ , which cannot be well approximated by classical orthogonal polynomials. 

# 2.2 Projection estimate

Let α, $\beta > - 1$ , and $\chi ^ { \alpha , \beta } ( t ) : = ( - \log t ) ^ { \alpha } t ^ { \beta }$ . For any $u \in L _ { \chi ^ { \alpha , \beta } } ^ { 2 } ( I )$ we denote $\pi _ { N } ^ { \alpha , \beta } \iota$ u the projection from $L _ { \chi ^ { \alpha , \beta } } ^ { 2 }$ to $\mathcal { P } _ { N } ^ { \log t }$ by 

$$
(u - \pi_ {N} ^ {\alpha , \beta} u, v) _ {\chi^ {\alpha , \beta}} = \int_ {0} ^ {1} \left\{u - \pi_ {N} ^ {\alpha , \beta} u \right\} (t) v (t) \chi^ {\alpha , \beta} (t) d t = 0, \quad \forall v \in \mathcal {P} _ {N} ^ {\log t}. \tag {2.12}
$$

$\{ S _ { n } ^ { ( \alpha , \beta ) } \} _ { n = 0 } ^ { \infty }$ 

$$
\pi_ {N} ^ {\alpha , \beta} u = \sum_ {n = 0} ^ {N} \hat {u} _ {n} ^ {\alpha , \beta} \mathcal {S} _ {n} ^ {(\alpha , \beta)}, \quad \hat {u} _ {n} ^ {\alpha , \beta} = (\gamma_ {n} ^ {(\alpha , \beta)}) ^ {- 1} \int_ {0} ^ {1} u (t) \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) \chi^ {\alpha , \beta} (t) d t. \tag {2.13}
$$

To better describe the projection error π α, $\pi _ { N } ^ { \alpha , \beta }$ u we define a pseudo-derivative 

$$
\widehat {\partial} _ {t} u := t \partial_ {t} u \tag {2.14}
$$

and a non-uniformly weighted Sobolev space 

$$
A _ {\alpha , \beta} ^ {k} (I) := \{v \in L _ {\chi^ {\alpha , \beta}} ^ {2} (I): \widehat {\partial} _ {t} ^ {j} v \in L _ {\chi^ {\alpha + j, \beta}} ^ {2} (I), j = 1, 2, \dots , k \}, \quad k \in \mathbb {N}, \tag {2.15}
$$

equipped with semi-norm and norm 

$$
| v | _ {A _ {\alpha , \beta} ^ {m}} := \| \widehat {\partial} _ {t} ^ {m} v \| _ {\chi^ {\alpha + m, \beta}}, \quad \| v \| _ {A _ {\alpha , \beta} ^ {m}} := \left(\sum_ {k = 0} ^ {m} | v | _ {A _ {\alpha , \beta} ^ {k}} ^ {2}\right) ^ {1 / 2}.
$$

Theorem 2.4 Let m, N, $k \in \mathbb N$ and $\alpha , \beta > - 1$ . For any $u \in A _ { \alpha , \beta } ^ { m } ( I )$ and $0 \leq k \leq \widetilde { m } : = \operatorname* { m i n } \{ m , N + 1 \}$ , we have 

$$
\left\| \widehat {\partial} _ {t} ^ {k} \left(u - \pi_ {N} ^ {\alpha , \beta} u\right) \right\| _ {\chi^ {\alpha + k, \beta}} \leq \sqrt {(\beta + 1) ^ {k - \widetilde {m}} \frac {(N - \widetilde {m} + 1) !}{(N - k + 1) !}} \| \widehat {\partial} _ {t} ^ {\widetilde {m}} u \| _ {\chi^ {\alpha + \widetilde {m}, \beta}}. \tag {2.16}
$$

In particular, fixing $\alpha = \beta = k = 0$ and $m < N + 1$ , it holds that 

$$
\left\| u - \pi_ {N} u \right\| \leq c N ^ {- \frac {m}{2}} \left\| \widehat {\partial} _ {t} ^ {m} u \right\| _ {\chi^ {m}}, \tag {2.17}
$$

where $\pi _ { N } = \pi _ { N } ^ { 0 , 0 }$ and $\chi ^ { m } = \chi ^ { m , 0 } = ( - \log t ) ^ { m }$ for notational simplicity. 

Proof. For any $u \in A _ { \alpha , \beta } ^ { m } ( I )$ , via relations (2.6) and (2.14), we have 

$$
\widehat {\partial} _ {t} ^ {l} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) = (\beta + 1) ^ {l} \mathcal {S} _ {n - l} ^ {(\alpha + l, \beta)} (t), \quad l \leq n. \tag {2.18}
$$

Then, it can be easily detected from the orthogonality P3 that 

$$
u (t) = \sum_ {n = 0} ^ {\infty} \hat {u} _ {n} ^ {\alpha , \beta} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t), \quad \| \widehat {\partial} _ {t} ^ {l} u \| _ {\chi^ {\alpha + l, \beta}} ^ {2} = \sum_ {n = l} ^ {\infty} (\beta + 1) ^ {2 l} \gamma_ {n - l} ^ {(\alpha + l, \beta)} | \hat {u} _ {n} ^ {\alpha , \beta} | ^ {2}, \quad l \geq 1.
$$

Therefore, 

$$
\begin{array}{l} \| \widehat {\partial} _ {t} ^ {k} (u - \pi_ {N} ^ {\alpha , \beta} u) \| _ {\chi^ {\alpha + k, \beta}} ^ {2} = \sum_ {n = N + 1} ^ {\infty} (\beta + 1) ^ {2 k} \gamma_ {n - k} ^ {(\alpha + k, \beta)} | \hat {u} _ {n} ^ {\alpha , \beta} | ^ {2} \\ \leq \max \left\{\frac {\gamma_ {n - k} ^ {(\alpha + k , \beta)}}{\gamma_ {n - \widetilde {m}} ^ {(\alpha + \widetilde {m} , \beta)}} \right\} \sum_ {n = N + 1} ^ {\infty} (\beta + 1) ^ {2 k} \gamma_ {n - \widetilde {m}} ^ {(\alpha + \widetilde {m}, \beta)} | \hat {u} _ {n} ^ {\alpha , \beta} | ^ {2} \\ \leq (\beta + 1) ^ {2 (k - \widetilde {m})} \frac {\gamma_ {N + 1 - k} ^ {(\alpha + k , \beta)}}{\gamma_ {N + 1 - \widetilde {m}} ^ {(\alpha + \widetilde {m} , \beta)}} \| \widehat {\partial} _ {t} ^ {\widetilde {m}} u \| _ {\chi^ {\alpha + \widetilde {m}, \beta}} ^ {2} \\ \leq (\beta + 1) ^ {k - \widetilde {m}} \frac {(N - \widetilde {m} + 1) !}{(N - k + 1) !} \| \widehat {\partial_ {t}} ^ {\widetilde {m}} u \| _ {\chi^ {\alpha + \widetilde {m}, \beta}} ^ {2}. \\ \end{array}
$$

Finally, the proof of the special case can be proved by the following useful result: for any constant $a , b \in \mathbb { R } , n \in \mathbb { N } , n + a > 1$ and $n + b > 1$ (see Zhao et al., 2013, Lemma 2.1), 

$$
\frac {\Gamma (n + a)}{\Gamma (n + b)} \leq \nu_ {n} ^ {a, b} n ^ {a - b}, \tag {2.19}
$$

where 

$$
v _ {n} ^ {a, b} = \exp \left(\frac {a - b}{2 (n + b - 1)} + \frac {1}{1 2 (n + a - 1)} + \frac {(a - b) ^ {2}}{n}\right). \tag {2.20}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/b8f0ec71af9666c9cfc4bb854ff23362356500b1d25666f8c1833f0dc580efcf.jpg)


Remark 2.5 The essential difference between approximations by LOFs and traditional polynomials can be explained by the estimate (2.17). In fact, since $\widehat { \partial } _ { t } t ^ { r } = r t ^ { r }$ , it’s easy to check that $\| \widehat { \partial } _ { t } ^ { m } t ^ { r } \| _ { \chi ^ { m } } < \infty$ for all $r \geq 0$ and any positive integer m. So the LOFs can approximate a function whose singularity behave as $\sum _ { i } c _ { i } t ^ { r _ { i } }$ with exponential convergence. 

On the contrary, the polynomial approximation error depends on the regularity defined by the usual derivative. Specifically, for the classical polynomial projection $\Pi _ { N } \colon L ^ { 2 } \to P _ { N } ^ { t } : = s p a n \{ 1 , t , \ldots , t ^ { N } \}$ , it holds that 

$$
\| u - \Pi_ {N} u \| \leq c N ^ {- m} \| \partial_ {t} ^ {m} u \|.
$$

Hence, functions behaving as $\textstyle \sum _ { i } c _ { i } t ^ { r _ { i } }$ with several small $r _ { i } ~ > ~ 0$ cannot be well approximated by polynomials. 

# 2.3 Interpolation estimate

Let {t (α,β)j } Nj=0 $\{ t _ { j } ^ { ( \alpha , \beta ) } \} _ { j = 0 } ^ { N }$ be the mapped Gauss points defined in (2.9). We define the mapped Lagrange functions 

$$
l _ {j} (y (t)) = \frac {\prod_ {i \neq j} \left(y (t) - y (t _ {i} ^ {(\alpha , \beta)})\right)}{\prod_ {i \neq j} \left(y (t _ {j} ^ {(\alpha , \beta)}) - y (t _ {i} ^ {(\alpha , \beta)})\right)} = \frac {\prod_ {i \neq j} \log (t _ {i} ^ {(\alpha , \beta)} / t)}{\prod_ {i \neq j} \log (t _ {i} ^ {(\alpha , \beta)} / t _ {j} ^ {(\alpha , \beta)})}, \tag {2.21}
$$

and the interpolation operator $\mathcal { T } _ { N } ^ { \alpha , \beta } : C ( I )  P _ { N } ^ { \mathrm { l o g } t }$ by 

$$
\mathcal {I} _ {N} ^ {\alpha , \beta} v (t) = \sum_ {j = 0} ^ {N} v \left(t _ {j} ^ {(\alpha , \beta)}\right) l _ {j} (y (t)). \tag {2.22}
$$

Obviously, we have $\mathcal { Z } _ { N } ^ { \alpha , \beta } \nu ( t _ { i } ^ { ( \alpha , \beta ) } ) = \nu ( t _ { i } ^ { ( \alpha , \beta ) } ) , j = 0 , 1 , \cdots , N$ 

We first establish a stability result. 

Theorem 2.6 For any $\nu \in C ( I ) \cap A _ { \alpha , \beta } ^ { 1 } ( I )$ and $\widehat { \partial } _ { t } \nu \in L _ { \chi ^ { \alpha , \beta } } ^ { 2 } ( I )$ , we have 

$$
\left\| \mathcal {I} _ {N} ^ {\alpha , \beta} v \right\| _ {\chi^ {\alpha , \beta}} \leq c \sqrt {(\beta + 1) ^ {\alpha}} \left(c _ {1} ^ {\beta} N ^ {- \frac {1}{2}} \| \widehat {\partial} _ {t} v \| _ {\chi^ {\alpha , \beta}} + c _ {2} ^ {\beta} \sqrt {\log N} \| v \| _ {A _ {\alpha , \beta} ^ {1}}\right), \tag {2.23}
$$

where $c _ { 1 } ^ { \beta } = ( \beta + 1 ) ^ { - \frac { 1 } { 2 } } , ~ c _ { 2 } ^ { \beta } = 2 \sqrt { \operatorname* { m a x } \{ 1 , \beta + 1 \} } .$ 

Proof. Let $t ( y ) = e ^ { - ( \beta + 1 ) ^ { - 1 } y } \mathrm { ~ a n d ~ } \tilde { \nu } ( y ) = \nu ( t ( y ) )$ . Via relations (2.22) and (2.9), we have 

$$
\mathcal {I} _ {N} ^ {\alpha , \beta} v (t) = \mathcal {I} _ {N} ^ {\alpha} \tilde {v} (y) := \sum_ {j = 0} ^ {N} \tilde {v} (y _ {j} ^ {(\alpha)}) l _ {j} (y), \quad y \in \mathbb {R} ^ {+}.
$$

Thanks to (Guo et al., 2006, (3.12) with $\beta = 1 )$ , we have 

$$
\| \mathcal {I} _ {N} ^ {\alpha} \tilde {v} \| _ {y ^ {\alpha} e ^ {- x}} \leq c (N ^ {- \frac {1}{2}} \sqrt {M _ {1} ^ {\tilde {v}}} + 2 \sqrt {\log N} \sqrt {M _ {2} ^ {\tilde {v}}}),
$$

where 

$$
M _ {1} ^ {\tilde {\nu}} = \int_ {0} ^ {\infty} (\partial_ {y} \tilde {\nu} (y)) ^ {2} y ^ {\alpha} e ^ {- y} d y, \quad M _ {2} ^ {\tilde {\nu}} = \int_ {0} ^ {\infty} \left(\tilde {\nu} ^ {2} + y (\partial_ {y} \tilde {\nu}) ^ {2}\right) y ^ {\alpha} e ^ {- y} d y.
$$

On the other hand, we have 

$$
\begin{array}{l} \int_ {0} ^ {\infty} \left(\partial_ {y} \tilde {v} (y)\right) ^ {2} y ^ {\alpha} e ^ {- y} d y = \int_ {0} ^ {1} \left(\frac {t}{\beta + 1} \left(\partial_ {t} v (t)\right) ^ {2} (- (\beta + 1) \log t)\right) ^ {\alpha} t ^ {\beta + 1} \frac {\beta + 1}{t} d t \\ = (\beta + 1) ^ {\alpha - 1} \int_ {0} ^ {1} (\widehat {\partial} _ {t} v (t)) ^ {2} (- \log t) ^ {\alpha} t ^ {\beta} \mathrm{d} t = (\beta + 1) ^ {\alpha - 1} \| \widehat {\partial} _ {t} v \| _ {\chi^ {\alpha , \beta}} ^ {2} \\ \end{array}
$$

and 

$$
\begin{array}{l} \int_ {0} ^ {\infty} \left(\tilde {v} ^ {2} + y \left(\partial_ {y} \tilde {v}\right) ^ {2}\right) y ^ {\alpha} e ^ {- y} d y = (\beta + 1) ^ {\alpha + 1} \int_ {0} ^ {1} \left(v ^ {2} + \frac {(- \log t)}{\beta + 1} \left(\widehat {\partial} _ {t} v\right) ^ {2}\right) (- \log t) ^ {\alpha} t ^ {\beta} d t \\ \leq (\beta + 1) ^ {\alpha} \max \{1, (\beta + 1) \} \| v \| _ {A _ {\alpha , \beta} ^ {1}} ^ {2}. \\ \end{array}
$$

We can then derive (2.23) by combing the above relations. 

With the above stability result in hand we can now establish an estimate for the interpolation error. 

Theorem 2.7 Let m and N be positive integers, and $\alpha , \beta \ > \ - 1$ . For any $\nu \in C ( I ) \cap A _ { \alpha , \beta } ^ { m } ( I )$ and $\widehat { \partial } _ { t } \nu \in A _ { \alpha , \beta } ^ { m - 1 } ( I )$ , we have 

$$
\left\| \mathcal {I} _ {N} ^ {\alpha , \beta} v - v \right\| _ {\chi^ {\alpha , \beta}} \leq c \sqrt {\frac {(N + 1 - \widetilde {m}) !}{(\beta + 1) ^ {\widetilde {m}} N !}} \left\{c _ {1} ^ {\beta} \| \widehat {\partial} _ {t} ^ {\widetilde {m}} v \| _ {\chi^ {\alpha + m - 1, \beta}} + c _ {2} ^ {\beta} \sqrt {\log N} \| \widehat {\partial} _ {t} ^ {\widetilde {m}} v \| _ {\chi^ {\alpha + m, \beta}} \right\}, \tag {2.24}
$$

where $c _ { 1 } ^ { \beta } = ( \beta + 1 ) ^ { - \frac 1 2 } , \quad c _ { 2 } ^ { \beta } = 2 \sqrt { \operatorname* { m a x } \{ 1 , \beta + 1 \} } \mathrm { ~ a n d } \widetilde m = \operatorname* { m i n } \{ m , N + 1 \} .$ 

Proof. By the triangle inequality, we have 

$$
\left\| \mathcal {I} _ {N} ^ {\alpha , \beta} v - v \right\| _ {\chi^ {\alpha , \beta}} \leq \left\| \mathcal {I} _ {N} ^ {\alpha , \beta} v - \pi_ {N} ^ {\alpha , \beta} v \right\| _ {\chi^ {\alpha , \beta}} + \left\| \pi_ {N} ^ {\alpha , \beta} v - v \right\| _ {\chi^ {\alpha , \beta}}. \tag {2.25}
$$

We only need to estimate the first term since the estimate for the second term is already available in Theorem 2.4. Thanks to (2.23), 

$$
\begin{array}{l} \| \mathcal {I} _ {N} ^ {\alpha , \beta} v - \pi_ {N} ^ {\alpha , \beta} v \| _ {\chi^ {\alpha , \beta}} = \| \mathcal {I} _ {N} ^ {\alpha , \beta} (v - \pi_ {N} ^ {\alpha , \beta} v) \| _ {\chi^ {\alpha , \beta}} \\ \leq c \sqrt {(\beta + 1) ^ {\alpha}} \left(c _ {1} ^ {\beta} N ^ {- \frac {1}{2}} \| \widehat {\partial} _ {t} (v - \pi_ {N} ^ {\alpha , \beta} v) \| _ {\chi^ {\alpha , \beta}} + c _ {2} ^ {\beta} \sqrt {\log N} \| v - \pi_ {N} ^ {\alpha , \beta} v \| _ {A _ {\alpha , \beta} ^ {1}}\right). \tag {2.26} \\ \end{array}
$$

The term v − πα,βN $\Vert \nu - \pi _ { N } ^ { \alpha , \beta } \nu \Vert _ { A _ { \alpha , \beta } ^ { \mu } } \ ( \mu = 0 , 1 )$ vAμα,β can be estimated through Theorem 2.4. For the first term in the last inequality, we have 

$$
\| \widehat {\partial} _ {t} (v - \pi_ {N} ^ {\alpha , \beta} v) \| _ {\chi^ {\alpha , \beta}} \leq \| \widehat {\partial} _ {t} v - \pi_ {N} ^ {\alpha , \beta} \{\widehat {\partial} _ {t} v \} \| _ {\chi^ {\alpha , \beta}} + \| \pi_ {N} ^ {\alpha , \beta} \{\widehat {\partial} _ {t} v \} - \widehat {\partial} _ {t} \{\pi_ {N} ^ {\alpha , \beta} v \} \| _ {\chi^ {\alpha , \beta}}. \tag {2.27}
$$

We now follow a classical procedure as in Bernardi & Maday (1997); Guo (1998) to derive the desired estimate. The starting point is the relation 

$$
\widehat {\partial} _ {t} v = \sum_ {n = 0} ^ {\infty} \hat {v} _ {n} ^ {\alpha , \beta} \widehat {\partial} _ {t} \mathcal {S} _ {n} ^ {(\alpha , \beta)} \stackrel {(2. 2 6)} {=} \sum_ {n = 1} ^ {\infty} \hat {v} _ {n} ^ {\alpha , \beta} \left((\beta + 1) \sum_ {l = 0} ^ {n - 1} \mathcal {S} _ {l} ^ {(\alpha , \beta)}\right) = \sum_ {l = 0} ^ {\infty} \left((\beta + 1) \sum_ {n = l + 1} ^ {\infty} \hat {v} _ {n} ^ {\alpha , \beta}\right) \mathcal {S} _ {l} ^ {(\alpha , \beta)}.
$$

The above equation implies that 

$$
\pi_ {N} ^ {\alpha , \beta} \{\widehat {\partial} _ {t} v \} (t) = \sum_ {n = 0} ^ {N} \hat {v} _ {1, n} ^ {\alpha , \beta} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t), \quad \hat {v} _ {1, n} ^ {\alpha , \beta} := (\beta + 1) \sum_ {k = n + 1} ^ {\infty} \hat {v} _ {k} ^ {\alpha , \beta}.
$$

Similarly, we have 

$$
\widehat {\partial} _ {t} \{\pi_ {N} ^ {\alpha , \beta} v \} (t) = \sum_ {n = 0} ^ {N - 1} \left((\beta + 1) \sum_ {k = n + 1} ^ {N} \hat {v} _ {k} ^ {\alpha , \beta}\right) \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) = \sum_ {n = 0} ^ {N - 1} (\hat {v} _ {1, n} ^ {\alpha , \beta} - \hat {v} _ {1, N} ^ {\alpha , \beta}) \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t).
$$

Hence, 

$$
\begin{array}{l} \| \pi_ {N} ^ {\alpha , \beta} \{\widehat {\partial} _ {t} v \} - \widehat {\partial} _ {t} \{\pi_ {N} ^ {\alpha , \beta} v \} \| _ {\chi^ {\alpha , \beta}} ^ {2} \stackrel {(2. 2 7)} {=} \sum_ {n = 0} ^ {N} \gamma_ {n} ^ {(\alpha , \beta)} (\hat {v} _ {1, N} ^ {\alpha , \beta}) ^ {2} = \gamma_ {N} ^ {(\alpha , \beta)} (\hat {v} _ {1, N} ^ {\alpha , \beta}) ^ {2} \sum_ {n = 0} ^ {N} \gamma_ {n} ^ {(\alpha , \beta)} (\gamma_ {N} ^ {(\alpha , \beta)}) ^ {- 1} \\ \leq \| \widehat {\partial} _ {t} v - \pi_ {N - 1} ^ {\alpha , \beta} \{\widehat {\partial} _ {t} v \} \| _ {\chi^ {\alpha , \beta}} ^ {2} \sum_ {n = 0} ^ {N} \gamma_ {n} ^ {(\alpha , \beta)} (\gamma_ {N} ^ {(\alpha , \beta)}) ^ {- 1}. \tag {2.28} \\ \end{array}
$$

It remains to estit’s obvious that $\begin{array} { r } { s _ { N } : = \sum _ { n = 0 } ^ { N } \gamma _ { n } ^ { ( \alpha , \beta ) } ( \gamma _ { N } ^ { ( \alpha , \beta ) } ) ^ { - 1 } } \end{array}$ . For any se Stirlin $\alpha \geq 0 .$ , in view of the expression of ula to deduce that for a large $\gamma _ { n } ^ { ( \alpha , \beta ) }$ $s _ { N } \leq N + 1 . \mathrm { F o r } - 1 < \alpha < 0$ 

$$
M \leq k \leq N,
$$

$$
\frac {\gamma_ {k} ^ {(\alpha , \beta)}}{\gamma_ {N} ^ {(\alpha , \beta)}} = \frac {\Gamma (N + 1) \Gamma (k + \alpha + 1)}{\Gamma (N + \alpha + 1) \Gamma (k + 1)} \sim N ^ {- \alpha} k ^ {\alpha}.
$$

Therefore, there exists a constant $c _ { M }$ such that 

$$
s _ {N} = \sum_ {n = 0} ^ {N} \gamma_ {n} ^ {(\alpha , \beta)} (\gamma_ {N} ^ {(\alpha , \beta)}) ^ {- 1} \leq N ^ {- \alpha} \left(c _ {M} + c \sum_ {k = M} ^ {N} k ^ {\alpha}\right) \leq c N. \tag {2.29}
$$

Finally, combing (2.25)–(2.29) and Theorem 2.4 leads to (2.24). 

Remark 2.8 Let {t (α,β)j } Nj=0 $\{ t _ { j } ^ { ( \alpha , \beta ) } \} _ { j = 0 } ^ { N }$ and $\{ \omega _ { j } ^ { ( \alpha , \beta ) } \} _ { j = 0 } ^ { N }$ be the same as (2.9). Then, we have the following estimate for the quadrature error: 

$$
\left| \int_ {0} ^ {1} v (t) \chi^ {\alpha , \beta} (t) \mathrm{d} t - \sum_ {j = 0} ^ {N} v (t _ {j} ^ {(\alpha , \beta)})   \omega_ {j} ^ {(\alpha , \beta)} \right| \leq \sqrt {\frac {\Gamma (\alpha + 1)}{(\beta + 1) ^ {\alpha + 1}}} \| \mathcal {I} _ {N} ^ {\alpha , \beta} v - v \| _ {\chi^ {\alpha , \beta}}. \tag {2.30}
$$

Indeed, the above estimate can be derived from 

$$
\sum_ {j = 0} ^ {N} v (t _ {j} ^ {(\alpha , \beta)}) \omega_ {j} ^ {(\alpha , \beta)} = \int_ {0} ^ {1} \mathcal {I} _ {N} ^ {\alpha , \beta} v (t) \chi^ {\alpha , \beta} (t) \mathrm{d} t,
$$

and 

$$
\int_ {0} ^ {1} (- \log t) ^ {\alpha} t ^ {\beta} \mathrm{d} t = \frac {1}{(\beta + 1) ^ {\alpha + 1}} \int_ {0} ^ {\infty} y ^ {\alpha} e ^ {- y} \mathrm{d} y = \frac {\Gamma (\alpha + 1)}{(\beta + 1) ^ {\alpha + 1}}.
$$

To understand better why the singular function $t ^ { s } , s > 0$ can be well approximated by LOFs, we plot distribution of the Gauss-LOFs quadrature nodes $\{ t _ { j } ^ { ( \alpha , \beta ) } \} _ { j = 0 } ^ { N }$ with various N and α, $\beta$ in Figs 1 and 2. We observe from Fig. 1 that the nodes cluster near zero, with significant more points near zero than the usual Gauss–Radau points. Figure 2 exhibits the influence of the parameters $( \alpha , \beta )$ on the nodes distribution. In particular, as α increases with $\beta$ fixed, the quadrature nodes move towards zero; on the other hand, as $\beta$ increases with α fixed, the quadrature nodes move away from zero. 

# 2.4 Numerical examples

We first demonstrate the accuracy of Gauss-LOFs quadrature for computing 

$$
\int_ {0} ^ {1} f (t) (- \log t) ^ {\alpha} t ^ {\beta} \mathrm{d} t, \quad \alpha , \beta > - 1, \tag {2.31}
$$

with the following functions $f ( t ) = \sin t , \ e ^ { t } , \ t ^ { - 1 / 3 }$ and $t ^ { 1 / 1 0 }$ , respectively. The quadrature errors are shown to the left of Fig. 3. We observe that the errors decay exponentially in all cases. We note that $f ( t ) = t ^ { - 1 / 3 }$ is singular and cannot be computed efficiently by the classical Gauss quadrature. However, $t ^ { - 1 / 3 }$ is smooth with the norm defined through the pseudo-derivative (2.14), so we achieve exponential convergence for this case as well. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/40246e248abdef35eeaa988b999fb11e0e8872f6d663fc1af5caddd0fe43ca55.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/73f450e82eda6a29ea4245a88cd1b98106ded5035327a0b6687927abd83282bc.jpg)



Fig. 1. Nodes distribution of S(α,β) n (t): α = 0, β = 5 with different N.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/edffe82a8ab70575199c9ca304c67349b1f6b88d72ce91b540b52c4fefd0f210.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/6889efcb7bc7e916d9f99a8864b292606ee920eead2261b9e3c0504013c5293c.jpg)



Fig. 2. Nodes distribution of $S _ { n } ^ { ( \alpha , \beta ) } ( t ) \colon N = 8 0$ with different α, β.


Next, we compute the projection error for $f ( t ) = t ^ { 1 / 1 0 }$ , which is not smooth in the usual Sobolev space, but it is smooth with the norm defined through the pseudo-derivative. To the right of Fig. 3 we plot the projection errors by using the shifted Legendre polynomial $L _ { n } ( 2 t - 1 )$ ), $t \in I$ and LOFs for function $\dot { f } ( \dot { t } ) = t ^ { 1 / 1 0 }$ with the fixed degree of basis $N = 4 0$ . We observe that the projection error by using LOFs is uniformly small across the interval [0,1], while the error by using the shifted Legendre polynomial is very large. 

# 3. Generalized Log orthogonal functions

The LOFs introduced in the last section is capable of resolving certain type of singularities at $t = 0$ , but LOFs $S _ { n } ^ { ( \alpha , \beta ) } ( t )$ consist of $\{ ( - \log t ) ^ { k } \} _ { k = 0 } ^ { n }$ , so grow very fast near t = 0 (cf. Fig. 4). This behavior may severely affect the accuracy in many situations. In addition, derivatives of LOFs involve the singular term $t ^ { - 1 }$ , so they are not suitable to serve as basis functions to represent solutions of partial differential equations (PDEs) or fractional PDEs. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/efd73ebcfeda07b66b22b4dd62f16f4d25a7a352be11422b38da729277590c8a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/16c8c7222bad0de1bc4cc1c699c45fb95d393b8d7eaa8a4fa07a36dfc69b35b3.jpg)



Fig. 3. Left: Gauss-LOFs quadrature errors. Right: projection errors.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/d0086d6b9c7a9fe0b429d4b5a062ce3543c0400b96baaf021986ba134359f534.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/4315fc82cc9881376573cdcb60d4914a4c59aa60087bcec228d39aeb16d808f4.jpg)



Fig. 4. Left: graphs of $y = S _ { n } ^ { ( \alpha , \beta ) } ( t )$ . Right: graphs of $y = S _ { n } ^ { ( \alpha , \beta , { \lambda } ) } ( t )$


Therefore, we shall introduce below the GLOFs, which are more suitable for numerical approximations of functions with weak singularities at one endpoint. 

# 3.1 Definition and properties

Remark 3.1 As depicted to the left of Fig. 4 values of LOFs near t = 0 vary in a very large range. However, as shown to the right of Fig. 4, GLOFs are much better behaved. 

Definition 3.2 (GLOFs). Let $\alpha , \beta > - 1 , \lambda \in \mathbb { R }$ . We define the GLOFs by 

$$
\mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} (t) := t ^ {\frac {\beta - \lambda}{2}} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t), \quad n \geq 0. \tag {3.1}
$$

In particular, $S _ { n } ^ { ( \alpha , \beta , \beta ) } ( t ) = S _ { n } ^ { ( \alpha , \beta ) } ( t )$ . 

GLOFs enjoy similar properties as those listed in Lemma 2.2 for LOFs. Owing to the relation (2.7) it is obvious that 

$$
\int_ {0} ^ {1} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} (t) \mathcal {S} _ {m} ^ {(\alpha , \beta , \lambda)} (t) (- \log t) ^ {\alpha} t ^ {\lambda} d t = \gamma_ {n} ^ {(\alpha , \beta)} \delta_ {m n}, \tag {3.2}
$$

where γ (n $\gamma _ { n } ^ { ( \alpha , \beta ) }$ is the same as the definition in (2.7). 

The derivative relation can be derived from the relation (2.6) and Definition 3.2. Indeed, 

$$
\begin{array}{l} \partial_ {t} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} (t) = \frac {\beta - \lambda}{2} t ^ {\frac {\beta - \lambda}{2} - 1} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) + t ^ {\frac {\beta - \lambda}{2}} \partial_ {t} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) \\ = \frac {\beta - \lambda}{2} t ^ {\frac {\beta - \lambda - 2}{2}} \mathcal {S} _ {n} ^ {(\alpha , \beta)} (t) + (\beta + 1) t ^ {\frac {\beta - \lambda - 2}{2}} \mathcal {S} _ {n - 1} ^ {(\alpha + 1, \beta)} (t) \\ = \frac {\beta - \lambda}{2} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda + 2)} (t) + (\beta + 1) \mathcal {S} _ {n - 1} ^ {(\alpha + 1, \beta , \lambda + 2)} (t). \tag {3.3} \\ \end{array}
$$

The pseudo-derivative with respect to GLOFs should be defined as 

$$
\widehat {\partial} _ {\gamma , t} u = t ^ {1 + \gamma} \partial_ {t} \{t ^ {- \gamma} u \}. \tag {3.4}
$$

Then, thanks to the definition of $S _ { n } ^ { ( \alpha , \beta , { \lambda } ) }$ and (2.6), we have the following important derivative relation: 

$$
(\beta + 1) ^ {- 1} \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} (t) = \mathcal {S} _ {n - 1} ^ {(\alpha + 1, \beta , \lambda)} (t) = \sum_ {l = 0} ^ {n - 1} \mathcal {S} _ {l} ^ {(\alpha , \beta , \lambda)} (t), \quad n \geq 1. \tag {3.5}
$$

Let 

$$
t _ {j} ^ {(\alpha , \beta , \lambda)} := t _ {j} ^ {(\alpha , \beta)}, \chi_ {j} ^ {(\alpha , \beta , \lambda)} := (t _ {j} ^ {(\alpha , \beta)}) ^ {\lambda - \beta} \chi_ {j} ^ {(\alpha , \beta)}, j = 0, 1, \dots , N, \tag {3.6}
$$

and denote 

$$
\mathcal {P} _ {N} ^ {\gamma , \log t} := \{t ^ {\gamma} p (t): p \in \mathcal {P} _ {N} ^ {\log t} \}. \tag {3.7}
$$

Then, we have the following Gauss-GLOFs quadrature: 

$$
\int_ {0} ^ {1} f (t) (- \log t) ^ {\alpha} t ^ {\lambda} \mathrm{d} t = \sum_ {j = 0} ^ {N} f (t _ {j} ^ {(\alpha , \beta , \lambda)})   \chi_ {j} ^ {(\alpha , \beta , \lambda)}, \quad \forall f \in \mathcal {P} _ {2 N + 1} ^ {\beta - \lambda , \log t}. \tag {3.8}
$$

In addition, we derive from Definition 2.1 and the closed form of the Laguerre polynomial that 

$$
\mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} (t) = \sum_ {k = 0} ^ {n} \frac {(- 1) ^ {k}}{k !} \binom {n + \alpha} {n - k} t ^ {\frac {\beta - \lambda}{2}} [ - (\beta + 1) \log t ] ^ {k}, \quad t \in I := (0, 1).
$$

# 3.2 Projection estimate

Let α, β > −1 and λ ∈ R. We define the projection operator πN $\beta > - 1$ $\lambda \in \mathbb { R }$ $\pi _ { N } ^ { \alpha , \beta , \lambda } \colon L _ { \chi ^ { \alpha , \lambda } } ^ { 2 } \to \mathcal { P } _ { N } ^ { \frac { \beta - \lambda } { 2 } , \log t } \ u _ { \mathbf { b } \mathbf { y } }$ α,β ,λ : L2 χ α,λ β−λ ,log t by 

$$
(u - \pi_ {N} ^ {\alpha , \beta , \lambda} u, v) _ {\chi^ {\alpha , \lambda}} = 0, \quad \forall u \in L _ {\chi^ {\alpha , \lambda}} ^ {2}, v \in \mathcal {P} _ {N} ^ {\frac {\beta - \lambda}{2}, \log t}, \tag {3.9}
$$

where $\chi ^ { \alpha , \lambda } ( t ) : = ( - \log t ) ^ { \alpha } t ^ { \lambda }$ . 

Thanks to the orthogonality of the basis $\{ S _ { n } ^ { ( \alpha , \beta , \lambda ) } \} _ { n = 0 } ^ { \infty }$ ) }∞n=0, we have 

$$
\pi_ {N} ^ {\alpha , \beta , \lambda} u = \sum_ {n = 0} ^ {N} \hat {u} _ {n} ^ {\alpha , \beta , \lambda} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} \text {   with   } \hat {u} _ {n} ^ {\alpha , \beta , \lambda} = (\gamma_ {n} ^ {(\alpha , \beta)}) ^ {- 1} \int_ {0} ^ {1} u (t) \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} (t) \chi^ {\alpha , \lambda} (t) d t. \tag {3.10}
$$

To better describe the approximability of $\pi _ { N } ^ { \alpha , \beta , \lambda }$ we define non-uniformly weighted Sobolev spaces 

$$
A _ {\alpha , \beta , \lambda} ^ {k} (I) := \left\{v \in L _ {\chi^ {\alpha , \lambda}} ^ {2} (I): \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {j} v \in L _ {\chi^ {\alpha + j, \lambda}} ^ {2} (I), j = 1, 2, \dots , k \right\}, \quad k \in \mathbb {N}, \tag {3.11}
$$

with the corresponding semi-norm and norm defined by 

$$
| v | _ {A _ {\alpha , \beta , \lambda} ^ {m}} := \| \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {m} v \| _ {\chi^ {\alpha + m, \lambda}}, \quad \| v \| _ {A _ {\alpha , \beta , \lambda} ^ {m}} := \left(\sum_ {k = 0} ^ {m} | v | _ {A _ {\alpha , \beta , \lambda} ^ {k}} ^ {2}\right) ^ {\frac {1}{2}}.
$$

Theorem 3.3 Let m, N, $k \in \mathbb { N } , \lambda \in \mathbb { R }$ and $\alpha , \beta \ > \ - 1$ . For any $u \in A _ { \alpha , \beta , \lambda } ^ { m } ( I )$ and $0 \leq k \leq \widetilde { m } =$ min $\{ m , N + 1 \}$ , we have 

$$
\left\| \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {k} \left(u - \pi_ {N} ^ {\alpha , \beta , \lambda} u\right) \right\| _ {\chi^ {\alpha + k, \lambda}} \leq \sqrt {(\beta + 1) ^ {k - \widetilde {m}} \frac {(N - \widetilde {m} + 1) !}{(N - k + 1) !}} \| \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {\widetilde {m}} u \| _ {\chi^ {\alpha + \widetilde {m}, \lambda}}, \tag {3.12}
$$

where $\widehat { \partial } _ { \frac { \beta - \lambda } { 2 } , t }$ is the pseudo-derivative defined in (3.4). 

Proof. For any $u \in A _ { \alpha , \beta , \lambda } ^ { m } ( I )$ we can expand it as $\begin{array} { r } { u = \sum _ { n = 0 } ^ { \infty } \hat { u } _ { n } ^ { \alpha , \beta , \lambda } S _ { n } ^ { ( \alpha , \beta , \lambda ) } } \end{array}$ 

$$
\widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {l} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} (t) \stackrel {(3. 5)} {=} (\beta + 1) ^ {l} \mathcal {S} _ {n - l} ^ {(\alpha + l, \beta , \lambda)} (t), \quad l \leq n,
$$

and (3.2), we have 

$$
\| \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {l} u \| _ {\chi^ {\alpha + l, \lambda}} ^ {2} = \sum_ {n = l} ^ {\infty} (\beta + 1) ^ {2 l} \gamma_ {n - l} ^ {(\alpha + l, \beta)} | \hat {u} _ {n} ^ {\alpha , \beta , \lambda} | ^ {2}, \quad l \geq 1.
$$

Then, by following the same procedure as in the proof of Theorem 2.4, we can obtain the desired result (3.12). - 

# 3.3 Interpolation estimate

L et {t α ,β } N $\{ t _ { j } ^ { \alpha , \beta } \} _ { j = 0 } ^ { N }$ be the same set of collocation points as for the LOFs. We define the interpolation operator $\mathcal { T } _ { N } ^ { \alpha , \beta , \dot { \lambda } } : \dot { C } ( I )  \mathcal { P } _ { N } ^ { \frac { \beta - \lambda } { 2 } , \log t } \mathfrak { b y }$ λ : C(I) → P β−2 N 

$$
(\mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} v) (t _ {j} ^ {\alpha , \beta}) = v (t _ {j} ^ {\alpha , \beta}), j = 0, 1, \dots , N.
$$

It is easy to see that 

$$
\mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} v (t) = \sum_ {j = 0} ^ {N} v (t _ {j} ^ {(\alpha , \beta)}) l _ {j} ^ {\beta , \lambda} \big (y (t) \big), \quad y (t) = - (\beta + 1) \log t, \tag {3.13}
$$

where $\{ l _ { j } ^ { \beta , \lambda } \}$ are the Lagrange ‘polynomials’ defined by 

$$
l _ {j} ^ {\beta , \lambda} (y (t)) = \frac {t ^ {\frac {\beta - \lambda}{2}} \prod_ {i \neq j} \log \left(t _ {i} ^ {(\alpha , \beta)} / t\right)}{\left(t _ {j} ^ {(\alpha , \beta)}\right) ^ {\frac {\beta - \lambda}{2}} \prod_ {i \neq j} \log \left(t _ {i} ^ {(\alpha , \beta)} / t _ {j} ^ {(\alpha , \beta)}\right)}. \tag {3.14}
$$

In view of (2.22) and (2.21), we have 

$$
\mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} v (t) = t ^ {\frac {\beta - \lambda}{2}} \mathcal {I} _ {N} ^ {\alpha , \beta} \{t ^ {\frac {\lambda - \beta}{2}} v (t) \} \in \mathcal {P} _ {N} ^ {\frac {\beta - \lambda}{2}, \log t}.
$$

Hence, we can derive the following result from Theorem 2.7. 

Theorem 3.4 Let m and N be positive integers, $\alpha , \beta > - 1$ and $\lambda \in \mathbb { R }$ . For any $\nu \in C ( I ) \cap A _ { \alpha , \beta , \lambda } ^ { m } ( I )$ （2号and $\widehat { \partial } _ { \frac { \beta - \lambda } { 2 } , t } \nu \in A _ { \alpha , \beta , \lambda } ^ { m - 1 } ( I )$ , we have

$$
\| \mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} v - v \| _ {\chi^ {\alpha , \lambda}} \leq c \sqrt {\frac {(N + 1 - \widetilde {m}) !}{(\beta + 1) ^ {\widetilde {m} - \alpha} N !}} \left\{c _ {1} ^ {\beta} \| \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {\widetilde {m}} v \| _ {\chi^ {\alpha + m - 1, \lambda}} + c _ {2} ^ {\beta} \sqrt {\log N} \| \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {\widetilde {m}} v \| _ {\chi^ {\alpha + m, \lambda}} \right\},
$$

where $c _ { 1 } ^ { \beta } = ( \beta + 1 ) ^ { - \frac 1 2 } , \quad c _ { 2 } ^ { \beta } = 2 \sqrt { \operatorname* { m a x } \{ 1 , \beta + 1 \} } \mathrm { a n d } \widetilde m = \operatorname* { m i n } \{ m , N + 1 \}$ . 

Proof. Since 

$$
\| \mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} v - v \| _ {\chi^ {\alpha , \lambda}} = \| \mathcal {I} _ {N} ^ {\alpha , \beta} \{t ^ {\frac {\lambda - \beta}{2}} v \} - t ^ {\frac {\lambda - \beta}{2}} v \| _ {\chi^ {\alpha , \beta}},
$$

and 

$$
\widehat {\partial} _ {t} \{t ^ {\frac {\lambda - \beta}{2}} v \} = t ^ {\frac {\lambda - \beta}{2}} \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} v \implies \widehat {\partial} _ {t} ^ {\widetilde {m}} \{t ^ {\frac {\lambda - \beta}{2}} v \} = t ^ {\frac {\lambda - \beta}{2}} \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {\widetilde {m}} v.
$$

We can then derive the desired result from the above relation and Theorem 2.8. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/e6f13253feb2697137fc6eb163b9fd26f5abdb844f1dc4d4ae768cdcc9372aa6.jpg)


# 3.4 Explicit error estimate for a class of weakly singular functions

The result in Theorem 3.4 is not easy to interpret for general functions, so we consider the following typical weakly singular functions 

$$
f (t) = t ^ {r} (- \log t) ^ {k}, \quad r \geq 0, k \in \mathbb {N} _ {0}.
$$

We first present a very useful relation of the Laguerre polynomials $\mathcal { L } _ { n } ^ { ( \alpha ) } ( y ) , \alpha > - 1$ . 

Lemma 3.5 Let $s > 0 , \alpha > - 1$ . For $k , n \in  { \mathbb { N } } _ { 0 }$ and $n > k$ , there exists 

$$
\int_ {0} ^ {\infty} y ^ {\alpha + k} e ^ {- s y} \mathscr {L} _ {n} ^ {(\alpha)} (y) \mathrm{d} y = \left(\frac {s - 1}{s}\right) ^ {n} \frac {k !}{s ^ {\alpha + k + 1}} \sum_ {j = 0} ^ {k} \frac {\Gamma (n - j + k + \alpha + 1)}{(j !) ^ {2} \Gamma (n - j + 1)} \left(\frac {s}{1 - s}\right) ^ {j}. \tag {3.15}
$$

Proof. The case $k = 0$ is a direct result of Gradshteyn & Ryzhik (2007, 7.414: 8), i.e., 

$$
\int_ {0} ^ {\infty} e ^ {- s y} y ^ {\alpha} \mathcal {L} _ {n} ^ {(\alpha)} (y) \mathrm{d} y = \frac {\Gamma (\alpha + n + 1)}{\Gamma (n + 1)} (1 - 1 / s) ^ {n} s ^ {- \alpha - 1}.
$$

For $k > 0 .$ , owing to Rodrigues’ formula (see Szego, 1975, (5.1.5)), we have that 

$$
\begin{array}{l} \int_ {0} ^ {\infty} y ^ {k} e ^ {- s y} y ^ {\alpha} \mathscr {L} _ {n} ^ {(\alpha)} (y) \mathrm{d} y = \frac {1}{n !} \int_ {0} ^ {\infty} y ^ {k} e ^ {(1 - s) y} \partial_ {y} ^ {n} (y ^ {n + \alpha} e ^ {- y}) \mathrm{d} y \\ = \frac {(- 1) ^ {n}}{n !} \int_ {0} ^ {\infty} \partial_ {y} ^ {n} (y ^ {k} e ^ {(1 - s) y}) y ^ {n + \alpha} e ^ {- y} d y \\ = \frac {(- 1) ^ {n}}{n !} \sum_ {j = 0} ^ {k} \binom {n} {j} (1 - s) ^ {n - j} \int_ {0} ^ {\infty} \frac {k !}{j !} y ^ {k - j} y ^ {n + \alpha} e ^ {- s y} d y \\ = (- 1) ^ {n} \sum_ {j = 0} ^ {k} \frac {k !   \Gamma (n - j + k + \alpha + 1)}{(j !) ^ {2} \Gamma (n - j + 1)} \frac {(1 - s) ^ {n - j}}{s ^ {n - j + \alpha + k + 1}}. \tag {3.16} \\ \end{array}
$$

One can easily check the equivalence of the relations (3.15) and (3.16), which completes the proof. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/1793cc023c113527a6124b01c9d1dbc8bfa515060f80217fdd2af6248016ecf7.jpg)


With the above lemma in hand we have the following error estimate: 

Theorem 3.6 Given $f ( t ) = t ^ { r } ( - \log t ) ^ { k } , r \geq 0 , k \in \mathbb { N } _ { 0 } . \operatorname { L e t } \lambda > - 1 - 2 r , \alpha , \beta > - 1 \mathrm { ~ a n d ~ } \beta > \lambda$ . Then, we have 

$$
f \in L _ {\chi^ {\alpha , \lambda}} ^ {2} \text {and} R _ {\alpha , \beta , \lambda} = \left| \frac {2 r + \lambda - \beta}{2 r + 2 + \lambda + \beta} \right| <   1,
$$

and 

$$
\left\| f - \pi_ {N} ^ {\alpha , \beta , \lambda} f \right\| _ {\chi^ {\alpha , \lambda}} \leq c (k + 1)! N ^ {\frac {\alpha + 1}{2} + k} \left(R _ {\alpha , \beta , \lambda}\right) ^ {N} \text {   when   } N > - \frac {2 k + \alpha + 2}{2 \log \left(R _ {r , \beta , \lambda}\right)}, \tag {3.17}
$$

where 

$$
c \approx \sqrt {\frac {2 ^ {\alpha + 1 + k} (\beta + 1) ^ {2 \alpha + 2 - k}}{(\beta + \lambda + 2 r + 2) ^ {\alpha + 1 + k}}}.
$$

Proof. Since $\lambda > - 1 - 2 r$ it is easy to check that $f \in L _ { \chi ^ { \alpha , \lambda } } ^ { 2 }$ λ nd Rα,β,λ = a $\begin{array} { r } { R _ { \alpha , \beta , \lambda } = \left| \frac { 2 r + \lambda - \beta } { 2 r + 2 + \lambda + \beta } \right| < 1 } \end{array}$ . Thanks to the orthogonality of the basis {S(α,β,λ) n }∞n=0 $\{ S _ { n } ^ { ( \alpha , \beta , \lambda ) } \} _ { n = 0 } ^ { \infty }$ we can write 

$$
f = \sum_ {n = 0} ^ {\infty} \hat {f} _ {n} ^ {\alpha , \beta , \lambda} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)}, \quad \pi_ {N} ^ {\alpha , \beta , \lambda} f = \sum_ {n = 0} ^ {N} \hat {f} _ {n} ^ {\alpha , \beta , \lambda} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)}
$$

with coefficients 

$$
\hat {f} _ {n} ^ {\alpha , \beta , \lambda} = (\gamma_ {n} ^ {(\alpha , \beta)}) ^ {- 1} \int_ {0} ^ {1} f (t)   \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} (t)   \chi^ {\alpha , \lambda} (t) \mathrm{d} t, \quad n = 0, 1, \ldots , N.
$$

Let $y ( t ) = - ( \beta + 1 )$ log t and $\mathrm { d } y = - ( \beta + 1 ) t ^ { - 1 } \mathrm { d } t$ . It holds that 

$$
\begin{array}{l} \hat {f} _ {n} ^ {\alpha , \beta , \lambda} = \frac {(\beta + 1) ^ {\alpha + 1} \Gamma (n + 1)}{\Gamma (n + \alpha + 1)} \int_ {0} ^ {1} t ^ {r} (- \log t) ^ {k} t ^ {(\beta - \lambda) / 2} \mathcal {L} _ {n} ^ {(\alpha)} (y (t)) (- \log t) ^ {\alpha} t ^ {\lambda} \mathrm{d} t \\ = \frac {(\beta + 1) ^ {- k} \Gamma (n + 1)}{\Gamma (n + \alpha + 1)} \int_ {0} ^ {\infty} \exp \left(- \frac {\beta + \lambda + 2 r + 2}{2 (\beta + 1)} y\right) y ^ {\alpha + k} \mathcal {L} ^ {(\alpha)} (y) \mathrm{d} y. \\ \end{array}
$$

Taking $s = ( \beta + \lambda + 2 r + 2 ) / ( 2 \beta + 2 )$ into (3.15), we have 

$$
\hat {f} _ {n} ^ {\alpha , \beta , \lambda} = \left(\frac {s - 1}{s}\right) ^ {n} \frac {\Gamma (n + 1) k ! (\beta + 1) ^ {- k}}{\Gamma (n + \alpha + 1) s ^ {\alpha + k + 1}} \sum_ {j = 0} ^ {k} \frac {\Gamma (n - j + k + \alpha + 1)}{(j !) ^ {2} \Gamma (n - j + 1)} \left(\frac {s - 1}{s}\right) ^ {- j}. \tag {3.18}
$$

Owing to Zhao et al. (2013, Lemma 2.1), we have 

$$
\frac {\Gamma (n + a)}{\Gamma (n + b)} \leq \nu_ {n} ^ {a, b} n ^ {a - b}, \quad n + a > 1 \text {   and   } n + b > 1, \tag {3.19}
$$

where 

$$
\nu_ {n} ^ {a, b} = \exp \left(\frac {a - b}{2 (n + b - 1)} + \frac {1}{1 2 (n + a - 1)} + \frac {(a - b) ^ {2}}{n}\right). \tag {3.20}
$$

Combing (3.18)–(3.20) and the fact that $R _ { r , \beta , \lambda } < 1$ , we have 

$$
\begin{array}{l} \| f - \pi_ {N} ^ {\alpha , \beta , \lambda} f \| _ {\chi^ {\alpha , \lambda}} ^ {2} = \sum_ {N + 1} ^ {\infty} | \hat {f} _ {n} ^ {\alpha , \beta , \lambda} | ^ {2} \gamma_ {n} ^ {\alpha , \beta} \leq C _ {\alpha , \beta} ^ {\lambda , r} \sum_ {N + 1} ^ {\infty} \nu_ {n} ^ {\alpha , 0} n ^ {2 k + \alpha} (R _ {r, \beta , \lambda}) ^ {2 n - 2 k} \\ \leq C _ {\alpha , \beta} ^ {\lambda , r} \nu_ {N} ^ {\alpha , 0} \int_ {N} ^ {\infty} (R _ {r, \beta , \lambda}) ^ {2 x - 2 k} x ^ {2 k + \alpha} d x, \\ \end{array}
$$

where $\begin{array} { r } { C _ { \alpha , \beta } ^ { \lambda , r } = \frac { 2 ^ { \alpha + 1 + k } ( \beta + 1 ) ^ { 2 \alpha + 2 - k } } { ( \beta + \lambda + 2 r + 2 ) ^ { \alpha + 1 + k } } \big ( ( k + 1 ) ! \big ) ^ { 2 } } \end{array}$ (β+λ+2r+2)α+1+k . Finally, as $( R _ { r , \beta , \lambda } ) ^ { 2 x - 2 k } x ^ { 2 k + \alpha + 2 }$ is a decreasing function of x when $\begin{array} { r } { N > - \frac { 2 k + \alpha + 2 } { 2 \log ( R _ { r , \beta , \lambda } ) } } \end{array}$ , we conclude that 

$$
\| f - \pi_ {N} ^ {\alpha , \beta , \lambda} f \| _ {\chi^ {\alpha , \lambda}} \leq \sqrt {C _ {\alpha , \beta} ^ {\lambda , r}} N ^ {\frac {\alpha + 1}{2} + k} (R _ {r, \beta , \lambda}) ^ {N - k} \text {when} N > - \frac {2 k + \alpha + 2}{2 \log (R _ {r , \beta , \lambda})}.
$$

The proof is complete. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/059c54e95c9dd4b32f43a2c10cc38e6509668a05b622008228dd8753bd7c8f14.jpg)


The above theorem provides an accurate estimate for the GLOFs to a large class of singular functions. In particular, by setting $\alpha = \lambda = 0 .$ , we have an estimate in $L ^ { 2 } { \mathrm { - n o r m } } .$ 

Corollary 3.7 For $f ( t ) = t ^ { r } ( - \log t ) ^ { k } , \ r \geq 0 , \ k \in \mathbb { N }$ , it holds that 

$$
\left\| f - \pi_ {N} ^ {0, \beta , 0} f \right\| \leq \sqrt {2} ^ {k} (\beta + 1) ^ {- k} k! N ^ {k} \sqrt {2 (\beta + 1) N} \left| \frac {2 r - \beta}{2 r + \beta + 2} \right| ^ {N - k}. \tag {3.21}
$$

In particular, for $f = t ^ { r } , \ r \geq 0$ , we have 

$$
\left\| f - \pi_ {N} ^ {0, \beta , 0} f \right\| \leq \sqrt {2 (\beta + 1) N} \left| \frac {2 r - \beta}{2 r + \beta + 2} \right| ^ {N}. \tag {3.22}
$$

In order to verify the above theoretical results we plot the error curves for the GLOFs approximation $\mathbf { t o } f ( t ) = t ^ { r }$ with various r in Fig. 5, left with $r \in ( 0 , 1 )$ and right with r being integers. We observe exponential convergence for all $r \geq 0$ . We also plot the error curves for the GLOFs approximation to $f ( \bar { t } ) \ = \ t ( - \log t ) ^ { k }$ and $f ( t ) = t ^ { 2 } ( - \log t ) ^ { k }$ in Fig. 6. We also observe exponential convergence in all cases. All these numerical results are consistent with the approximation results in Corollary 3.7. 

# 4. Application to fractional differential equations

In this section we shall use GLOFs as the basis functions to solve some typical FDEs. 

We first review the definitions of Riemann–Liouville and Caputo fractional integrals and fractional derivatives (see e.g., Samko et al., 1993; Podlubny, 1999). 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/0a7676cf59434e6eca2ab3bbbba2af76bc87587023a3fffea8b4f90eea09c315.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/0b64df2d7da0539c02fa464c8d626c53184d1f3dbedafa929e666e23cc59dfed.jpg)



Fig. 5. GLOFs Approximation: $S _ { n } ^ { ( \alpha , \beta , { \lambda } ) }$ , $\alpha = \lambda = 0$ .


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/be65023a81134897a0fa2fd8b27721de5a595414ac2f5c9fbdeecfe2a007f06e.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/007c0123768025862585803f9ccac7d772700e2c44795960b9f2c83929957811.jpg)



Fig. 6. GLOFs Approximation: $S _ { n } ^ { ( \alpha , \beta , { \lambda } ) }$ , α = λ = 0.


Definition 4.1 (Fractional integrals and derivatives). For $t \in I = ( 0 , 1 )$ and $\rho \in \mathbb { R } ^ { + }$ the left and right fractional integrals are respectively defined as 

$$
{ } _ { 0 } \mathrm{I} _ { t } ^ { \rho } f ( t ) = \frac { 1 } { \Gamma ( \rho ) } \int _ { 0 } ^ { t } \frac { f ( \tau ) } { ( t - \tau ) ^ { 1 - \rho } } \mathrm{d} \tau , \quad { } _ { t } \mathrm{I} _ { 1 } ^ { \rho } f ( t ) = \frac { 1 } { \Gamma ( \rho ) } \int _ { t } ^ { 1 } \frac { f ( \tau ) } { ( \tau - t ) ^ { 1 - \rho } } \mathrm{d} \tau . \tag {4.1}
$$

For real $s \in [ k - 1 , k )$ with $k \in \mathbb N$ the Riemann–Liouville fractional derivatives are defined by 

$$
{ } _ { 0 } \mathrm{D} _ { t } ^ { s } f ( t ) = \frac { \mathrm{d} ^ { k } } { \mathrm{d} t ^ { k } } \left\{ _ { 0 } \mathrm{I} _ { t } ^ { k - s } f ( t ) \right\} , \quad { } _ { t } \mathrm{D} _ { 1 } ^ { s } f ( t ) = ( - 1 ) ^ { k } \frac { \mathrm{d} ^ { k } } { \mathrm{d} t ^ { k } } \left\{ _ { t } \mathrm{I} _ { 1 } ^ { k - s } f ( t ) \right\} . \tag {4.2}
$$

The Caputo fractional derivative of order s is defined by 

$$
{ } _ { 0 } ^ { C } \mathrm{D} _ { t } ^ { s } f ( t ) = { } _ { 0 } \mathrm{I} _ { t } ^ { k - s } \left\{ \frac { \mathrm{d} ^ { k } } { \mathrm{d} t ^ { k } } f ( t ) \right\} , \quad { } _ { t } ^ { C } \mathrm{D} _ { 1 } ^ { s } f ( t ) = ( - 1 ) ^ { k } { } _ { t } \mathrm{I} _ { 1 } ^ { k - s } \left\{ \frac { \mathrm{d} ^ { k } } { \mathrm{d} t ^ { k } } f ( t ) \right\} . \tag {4.3}
$$

# 4.1 An initial value problem

Given g, $q \in L ^ { 2 } ( I )$ we consider the following Caputo FDE of order $\nu \in ( 0 , 1 )$ : 

$$
{ } _ { 0 } ^ { C } \mathrm{D} _ { t } ^ { \nu } u ( t ) + q ( t ) u ( t ) = g ( t ) , \quad t \in I ; \qquad u ( 0 ) = u _ { 0 } . \tag {4.4}
$$

We shall first homogenize the initial condition. Setting $u = \nu + u _ { 0 }$ into the above equation we find that the problem (4.4) is equivalent to 

$$
{ } _ { 0 } ^ { C } \mathrm{D} _ { t } ^ { v } v ( t ) + q ( t ) v ( t ) = g ( t ) - u _ { 0 } q ( t ) , \quad v ( 0 ) = 0 . \tag {4.5}
$$

A main difficulty in obtaining accurate approximate solution of (4.5) is that the solution of this problem is weakly singular at $t = 0$ even if q and g are smooth. To design an effective approach to deal with this difficulty we need to understand the nature of this singularity. 

Applying $0 ^ { \mathrm { I } _ { t } ^ { \nu } }$ into both sides of (4.5) and using the fact that $_ 0 \mathbf { I } _ { t 0 } ^ { s } \mathbf { I } _ { t } ^ { r } = _ { 0 } \mathbf { I } _ { t } ^ { s + r }$ , we find 

$$
v (t) + \frac {1}{\Gamma (\nu)} \int_ {0} ^ {t} (t - \tau) ^ {\nu - 1} q (\tau) v (\tau) \mathrm{d} \tau = _ {0} \mathrm{I} _ {t} ^ {\nu} \{g - u _ {0} q \} (t).
$$

We then find from Cao et al. (2003, Theorem 2.1) that the solution near $t = 0$ behaves like 

$$
v (t) = \sum_ {i = 0} ^ {\infty} \sum_ {j = 1} ^ {\infty} \tilde {v} _ {i j} t ^ {i + j \nu}. \tag {4.6}
$$

This is why usual approximations based on global or piece-wise polynomials cannot approximate v(t) well. On the other hand, based on the analysis from the last section, the GLOFs are particularly suitable for this problem. 

Let us define $X _ { N } ^ { 0 } = \{ t ^ { \frac { \beta - \lambda } { 2 } } p : p \in \mathcal { P } _ { N } ^ { \log t } , \beta > \lambda \}$ . Then, the GLOF-Galerkin method for (4.5) is as follows: find $\nu _ { N } \in X _ { N } ^ { 0 }$ such that 

$$
\left(_ {0} ^ {C} \mathrm{D} _ {t} ^ {\nu} v _ {N}, w\right) + (q v _ {N}, w) = \left(\mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} \{g - u _ {0} q \}, w\right), \quad \forall w \in X _ {N} ^ {0}. \tag {4.7}
$$

Writing 

$$
v _ {N} = \sum_ {n = 0} ^ {N} \tilde {v} _ {n} ^ {\alpha , \beta , \lambda} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)}, \quad \bar {v} = \left(\tilde {v} _ {0} ^ {\alpha , \beta , \lambda}, \tilde {v} _ {1} ^ {\alpha , \beta , \lambda}, \dots , \tilde {v} _ {N} ^ {\alpha , \beta , \lambda}\right) ^ {t},
$$

and setting 

$$
S _ {k j} = \bigl (^ C _ {0} \mathrm{D} _ {t} ^ {\nu} \mathcal {S} _ {j} ^ {(\alpha , \beta , \lambda)}, \mathcal {S} _ {k} ^ {(\alpha , \beta , \lambda)} \bigr), \quad S = (S _ {k j}),
$$

$$
M _ {k j} = \big (q \mathcal {S} _ {j} ^ {(\alpha , \beta , \lambda)}, \mathcal {S} _ {k} ^ {(\alpha , \beta , \lambda)} \big), \quad M = (M _ {k j}), \tag {4.8}
$$

$$
f _ {j} = \Big (\mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} \{g - u _ {0} q \}, \mathcal {S} _ {j} ^ {(\alpha , \beta , \lambda)} \Big), \quad \bar {f} = (f _ {0}, f _ {1}, \ldots , f _ {N}) ^ {t},
$$

then (4.7) reduces to the following linear system 

$$
(S + M) \bar {v} = \bar {f}. \tag {4.9}
$$

The entries of M and $\bar { \boldsymbol { f } }$ can be computed accurately by using the Gauss-LOFs quadrature formula, but the computation of the stiffness matrix S needs special care. 

Indeed, for any v, $w \in X _ { N } ^ { 0 }$ 

$$
\left(_ {0} ^ {C} \mathrm{D} _ {t} ^ {\nu} v, w\right) = \int_ {0} ^ {1} \frac {1}{\Gamma (1 - \nu)} \int_ {0} ^ {t} \frac {v ^ {\prime} (s)}{(t - s) ^ {\nu}}   \mathrm{d} s   w (t) \mathrm{d} t \stackrel {{s = t \tau}} {{=}} \frac {1}{\Gamma (1 - \nu)} \int_ {0} ^ {1} \int_ {0} ^ {1} \frac {v ^ {\prime} (t \tau)}{(1 - \tau) ^ {\nu}}   \mathrm{d} \tau   w (t) t ^ {1 - \nu}   \mathrm{d} t.
$$

Note that the integrand in the above is weakly singular as $t  0$ and $\tau \ \to \ 1$ . In order to compute accurately the inner integral, we split it into two terms 

$$
\begin{array}{l} \int_ {0} ^ {1} \frac {v ^ {\prime} (t \tau)}{(1 - \tau) ^ {\nu}} d \tau = \int_ {0} ^ {\frac {1}{2}} v ^ {\prime} (t \tau) (1 - \tau) ^ {- \nu} d \tau + \int_ {\frac {1}{2}} ^ {1} v ^ {\prime} (t \tau) (1 - \tau) ^ {- \nu} d \tau \\ = \frac {1}{2} \int_ {0} ^ {1} v ^ {\prime} \left(\frac {t \tau}{2}\right) \left(1 - \frac {\tau}{2}\right) ^ {- \nu} d \tau + \frac {1}{4 ^ {1 - \nu}} \int_ {- 1} ^ {1} v ^ {\prime} \left(\frac {t (\xi + 3)}{4}\right) (1 - \xi) ^ {- \nu} d \xi . \\ \end{array}
$$

Hence, 

$$
\begin{array}{l} \left(_ {0} ^ {C} \mathrm{D} _ {t} ^ {\nu} v, w\right) = \frac {1}{2 \Gamma (1 - \nu)} \int_ {0} ^ {1} \int_ {0} ^ {1} v ^ {\prime} \left(\frac {t \tau}{2}\right) \left(1 - \frac {\tau}{2}\right) ^ {- \nu} \mathrm{d} \tau w (t) t ^ {1 - \nu} \mathrm{d} t \\ + \frac {1}{4 ^ {1 - \nu} \Gamma (1 - \nu)} \int_ {0} ^ {1} \int_ {- 1} ^ {1} v ^ {\prime} \left(\frac {t (\xi + 3)}{4}\right) (1 - \xi) ^ {- \nu} d \xi w (t) t ^ {1 - \nu} d t. \tag {4.10} \\ \end{array}
$$

The first term has weak singularity as $t , \tau \to 0$ while the second term has weak singularities as $t  0$ and $\xi \ \to \ 1$ . Therefore, the first term can be computed by using the tensor product of Gauss-GLOF quadratures (in t and τ ), which is effective with weak singularities as $t , \tau \to 0 .$ , and the second term can be computed by using the tensor product of Gauss-GLOF quadrature in t and of Gauss–Jacobi quadrature in $\xi$ with weight function $( 1 - \xi ) ^ { - \nu }$ . More precisely, 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/65697a6e510b28cb2b9351b8df10daf7281dacc69da1a6942b234dbd2fa531b6.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/1478ec1dc37b91312314b60148dd169c7602a41e78783eb177410a96de5b24f1.jpg)



Fig. 7. Left: (4.11) with $\kappa = 1 .$ . Right: $g = \cos t .$


$$
\begin{array}{l} (^ {C} _ {0} \mathrm{D} _ {t} ^ {\nu} v _ {N}, w) \approx \frac {1}{2 \Gamma (1 - \nu)} \sum_ {i = 0} ^ {N _ {I}} \sum_ {j = 0} ^ {N _ {I}} v _ {N} ^ {\prime} \left(\frac {t _ {i} t _ {j}}{2}\right) \left(1 - \frac {t _ {j}}{2}\right) ^ {- \nu} w (t _ {i}) t _ {i} ^ {1 - \nu} \chi_ {i} \chi_ {j} \\ + \frac {1}{4 ^ {1 - \nu} \Gamma (1 - \nu)} \sum_ {i = 0} ^ {N _ {I}} \sum_ {j = 0} ^ {N _ {I}} v _ {N} ^ {\prime} \left(\frac {t _ {i} (\xi_ {j} + 3)}{4}\right) w (t _ {i}) t _ {i} ^ {1 - \nu} \chi_ {i} \eta_ {j}, \\ \end{array}
$$

where $N _ { I } \geq N$ is a suitable number, $\{ t _ { i } , \chi _ { i } \} _ { i = 0 } ^ { N _ { I } }$ are the Gauss-LOFs nodes with weight function $\chi ^ { 0 , 0 } \equiv 1$ and $\{ \xi _ { i } , \eta _ { i } \} _ { i = 0 } ^ { N _ { I } }$ are the Gauss–Jacobi nodes with weight function $( 1 - \tau ) ^ { - \nu }$ . 

We present below some numerical results. We consider 

$$
{ } _ { 0 } ^ { C } \mathrm{D} _ { t } ^ { \nu } u ( t ) + \mathcal { K } u ( t ) = 0 , \quad u ( 0 ) = 1 , \tag {4.11}
$$

whose solution is Diethelm (2010, Theorem 4.3) $u ( t ) = E _ { \nu } ( - \mathcal { K } t ^ { \nu } )$ , where $E _ { \gamma } ( z )$ is the Mittag–Leffler function 

$$
E _ {\gamma} (z) = \sum_ {j = 0} ^ {\infty} \frac {z ^ {j}}{\Gamma (\gamma j + 1)}. \tag {4.12}
$$

We fix the parameters $\alpha = 0 , \beta = 5 , \lambda = 0$ , and plot to the left of Fig. 7 the convergence rates for various values of ν with $\kappa = 1$ . It is clear that the solution is not smooth in the classical Sobolev space, but it is smooth in the space defined through the pseudo-derivative, so we still obtain an exponential convergence rate. 

Next, we consider 

$$
{ } _ { 0 } ^ { C } \mathrm{D} _ { t } ^ { \nu } u ( t ) + ( 1 + \sin t ) u ( t ) = \cos t , \qquad u ( 0 ) = 1 ,
$$

for which the exact solution is unknown. Although the coefficients are smooth, the exact solution is expected to be weakly singular near zero, but smooth in the space defined through the pseudo-derivative. We fix the parameters $\alpha = 0 , \beta = 5 , \lambda = 0$ , and plot to the right of Fig. 7 the convergence rates for various values of ν. We obtain again exponential convergence rates. 

# 4.2 A boundary value problem

We consider 

$$
\left\{ \begin{array}{l l} - _ {0} \mathrm{D} _ {t} ^ {\mu} u (t) + q (t) u (t) = g (t), & t \in I, \\ u (0) = 0, & u (1) = 0, \end{array} \right. \tag {4.13}
$$

where $\mu \in ( 1 , 2 )$ and $g , q$ are given functions. 

Similar to the initial problem (4.4) the solution of the above problem is usually weakly singular, even with smooth g and $q .$ However, it can be approximated accurately by GLOFs since the solution is smooth in the space defined through the pseudo-derivative (3.4). 

Let us denote 

$$
X _ {N} ^ {0, 0} := \operatorname{span} \left\{\phi_ {n} = \frac {n}{n + \alpha} \mathcal {S} _ {n} ^ {(\alpha , \beta , \lambda)} - \mathcal {S} _ {n - 1} ^ {(\alpha , \beta , \lambda)}: n = 1, 2, \ldots , N, \beta > \lambda \right\}.
$$

Note that we have $\phi _ { n } ( 0 ) = \phi _ { n } ( 1 ) = 0$ for $n \geq 1$ . Our GLOF Galerkin method is as follows: find $u _ { N } \in X _ { N } ^ { 0 , 0 }$ s uch that 

$$
- (_ {0} \mathrm{D} _ {t} ^ {\mu} u _ {N}, w) + (q u _ {N}, w) = (\mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} g, w), \quad \forall w \in X _ {N} ^ {0, 0}. \tag {4.14}
$$

The stiffness and mass matrices of the above problem can be formulated as in the case of the IVP considered above. 

We now present some numerical results. We first take $q ( t ) ~ = ~ e ^ { t }$ and the exact solution to be $u ( t ) = t ^ { 3 / 2 } ( 1 - t )$ . The convergence rate is shown to the left of the Fig. 8. We then take $q ( t ) = e ^ { t }$ and $g ( t ) = t \sin t$ . In this case the exact solution is not known explicitly, so we used a very fine mesh to compute a reference solution. The convergence rate is shown to the right of the Fig. 8. We observe that the error converges exponentially in both cases despite the fact that the solutions are weakly singular near $t = 0$ . 

# 4.3 Error analysis

We carry out below error analysis for the GLOF Galerkin schemes for both the initial and boundary value problems. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/08fd66ee25f7ae906a3d697cc2d142e60614807041eb0789300c34841395df76.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/cfd8fec612e734a72fc909d24eb345e8535a2c8d77526a86c26bc5c68a699536.jpg)



Fig. 8. Left: $q ( t ) = e ^ { t }$ and $u = t ^ { 3 / 2 } ( 1$ − t). Right: $q ( t ) = e ^ { t }$ and $g = t \sin i$ .


We define 

• for $\ 0 \leq s < { \frac { 1 } { 2 } }$ , 

$$
H _ {0} ^ {s} (I) := \{f \in L ^ {2} (I): _ {0} \mathbf {D} _ {t} ^ {s} f, _ {t} \mathbf {D} _ {1} ^ {s} f \in L ^ {2} (I) \};
$$

• for $\begin{array} { r } { \frac { 1 } { 2 } < s \leq 1 } \end{array}$ 

$$
H _ {0} ^ {s} (I) := \{f \in L ^ {2} (I): _ {0} \mathrm{D} _ {t} ^ {s} f, _ {t} \mathrm{D} _ {1} ^ {s} f \in L ^ {2} (I), f (0) = f (1) = 0 \};
$$

equipped with norm and semi-norm 

$$
\| f \| _ {s} = \sqrt {\| f \| ^ {2} + | f | _ {s} ^ {2}} \quad \mathrm{with} \quad | f | _ {s} = \sqrt {\| _ {0} \mathrm{D} _ {t} ^ {s} f \| ^ {2} + \| _ {t} \mathrm{D} _ {1} ^ {s} f \| ^ {2}}.
$$

It can be shown that the above definitions coincide with the usual definitions by space interpolation. 

To avoid repetition we use the following weak formulation for both problems (4.5) and (4.13): find $P \in H _ { 0 } ^ { s } ( I )$ such that 

$$
a (P, w) := (- 1) ^ {\sigma (s)} (_ {0} \mathrm{D} _ {t} ^ {s} P, _ {t} \mathrm{D} _ {1} ^ {s} w) + (q P, w) = (Q, w), \quad \forall   w \in H _ {0} ^ {s} (I), \tag {4.15}
$$

where 

• for $( 4 . 5 ) \colon s = \nu / 2 , \sigma ( s ) = 0 , P ( t ) = \nu ( t ) , Q ( t ) = \{ g - u _ { 0 } q \} ( t ) ;$ 

• for (4.13): $s = \mu / 2 , \sigma ( s ) = 1 , P ( t ) = u ( t ) , Q ( t ) = g ( t ) .$ . 

The error analysis follows similar procedures used in Ervin & Roop (2006) and Li & Xu (2009, 2010). We first recall some useful results. 

Lemma 4.2 Let $s , r \in [ 0 , 1 ] / \{ \frac { 1 } { 2 } \}$ and $s \leq r .$ . For any f , $h \in H _ { 0 } ^ { r } ( I )$ , there exists 

• (Ervin & Roop, 2006, Corollary 2.15) 

$$
\| f \| \leq c _ {1} | f | _ {s} \leq c _ {2} | f | _ {r}. \tag {4.16}
$$

• (Li & Xu, 2010, Lemma 2.8) 

$$
(_ {0} \mathrm{D} _ {t} ^ {2 s} f, h) = (_ {0} \mathrm{D} _ {t} ^ {s} f, _ {t} \mathrm{D} _ {1} ^ {s} h). \tag {4.17}
$$

• (Li & Xu, 2010, Lemma 2.6) 

$$
c _ {1} | f | _ {s} ^ {2} \leq \frac {\left(_ {0} \mathrm{D} _ {t} ^ {s} f , {} _ {t} \mathrm{D} _ {1} ^ {s} f\right)}{\cos (s \pi)} \leq c _ {2} | f | _ {s} ^ {2}, \tag {4.18}
$$

where $c _ { 1 } , \ : c _ { 2 }$ are two positive constants independent of function f . 

Thanks to relation (4.17) and the identity below 

$$
(_ {0} \mathrm{D} _ {t} ^ {\mu} u, w) = (_ {0} \mathrm{D} _ {t} ^ {s} u, _ {t} \mathrm{D} _ {1} ^ {s} w), \quad \mu = 2 s, s \in (1 / 2, 1),
$$

we can rewrite (4.7) and (4.14) as follows: find $P _ { N } \in X _ { N }$ such that 

$$
a (P _ {N}, w) = (\mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} Q, w), \quad \forall   w \in X _ {N}, \tag {4.19}
$$

where $X _ { N } = X _ { N } ^ { 0 }$ for (4.7), and $X _ { N } = X _ { N } ^ { 0 , 0 }$ for (4.14). 

Lemma 4.3 If q(t) is bounded and min $q ( t ) \geq 0$ , then there exist $c ( s ) , c ( q ) > 0$ such that for any t∈[0,1] $P , \tilde { P } \in H _ { 0 } ^ { s } ( I )$ , we have 

$$
c (s) \| P \| _ {s} ^ {2} \leq a (P, P), \quad a (P, \tilde {P}) \leq c (q) \| P \| _ {s} \| \tilde {P} \| _ {s}. \tag {4.20}
$$

Proof. Due to the fact $( - 1 ) ^ { \sigma ( s ) } \cos ( s \pi ) > 0$ and the relation (4.18), it is easy to derive that 

$$
c _ {1} (- 1) ^ {\sigma (s)} \cos (s \pi) | P | _ {s} ^ {2} \leq (- 1) ^ {\sigma (s)} (_ {0} \mathrm{D} _ {t} ^ {s} P, _ {t} \mathrm{D} _ {1} ^ {s} P) \leq a (P, P),
$$

where $c _ { 1 }$ is a constant independent of the function P. Then, we can then derive the first inequality in (4.20) from the generalized Poincare inequality (4.16). The second inequality in (4.20) is a direct consequence of Cauchy–Schwarz inequality. - 

Thanks to (4.20) the existence uniqueness of the weak formulation (4.15) and the schemes (4.19) follows immediately from the Lax–Milgram Lemma. 

As for the error estimate, we have the following result. 

Theorem 4.4 $\mathrm { L e t } - 1 < \alpha \le 0 , \lambda \le 0$ and $\beta > 1 . \operatorname { L e t } P$ and $P _ { N }$ be respectively the solution of (4.15) and (4.19) with $\operatorname* { m i n } _ { t \in [ 0 , 1 ] } q ( t ) \geq 0$ . Then, we have 

$$
\begin{array}{l} \| P - P _ {N} \| _ {s} \leq c \sqrt {(\beta - 1) ^ {- \widetilde {m}} \frac {(N - \widetilde {m} + 1) !}{(N + 1) !}} \| \widehat {\partial} _ {\frac {\beta - \lambda - 2}{2}, t} ^ {\widetilde {m} + 1} P \| _ {\chi^ {\alpha + \widetilde {m}, \lambda}} \\ + c \sqrt {\frac {(N + 1 - \widetilde {m}) !}{(\beta + 1) ^ {\widetilde {m} - \alpha} N !}} \left\{c _ {1} ^ {\beta} \| \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {\widetilde {m}} Q \| _ {\chi^ {\alpha + m - 1, \lambda}} + c _ {2} ^ {\beta} \sqrt {\log N} \| \widehat {\partial} _ {\frac {\beta - \lambda}{2}, t} ^ {\widetilde {m}} Q \| _ {\chi^ {\alpha + m, \lambda}} \right\}, \tag {4.21} \\ \end{array}
$$

where $c _ { 1 } ^ { \beta } = ( \beta + 1 ) ^ { - \frac 1 2 } , \quad c _ { 2 } ^ { \beta } = 2 \sqrt { \operatorname* { m a x } \{ 1 , \beta + 1 \} } \mathrm { ~ a n d } \widetilde m = \operatorname* { m i n } \{ m , N + 1 \}$ 

Proof. For any $w _ { N } \in X _ { N }$ we derive from (4.15) and (4.19) that 

$$
a (P - P _ {N}, w _ {N}) = (Q - \mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} Q, w _ {N}). \tag {4.22}
$$

Let $e _ { N } = P _ { N } - w _ { N }$ , we have 

$$
c (s) \| e _ {N} \| _ {s} ^ {2} \leq a (e _ {N}, e _ {N}) = a (P - w _ {N}, e _ {N}) + a (P _ {N} - P, e _ {N}). \tag {4.23}
$$

Take $w _ { N } = e _ { N }$ in (4.22), we find 

$$
a (P _ {N} - P, e _ {N}) = (\mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} Q - Q, e _ {N}) \leq \| \mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} Q - Q \| \| e _ {N} \|. \tag {4.24}
$$

We then derive from (4.23) and (4.20) that 

$$
c (s) \| w _ {N} - P _ {N} \| _ {s} \leq c (\| P - w _ {N} \| _ {s} + \| \mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} Q - Q \|),
$$

which, along with (4.24), implies that 

$$
\| P - P _ {N} \| _ {s} \leq \| P - w _ {N} \| _ {s} + \| w _ {N} - P _ {N} \| _ {s} \leq c \inf _ {w _ {N} \in X _ {N}} \| P - w _ {N} \| _ {s} + c \| \mathcal {I} _ {N} ^ {\alpha , \beta , \lambda} Q - Q \|. \tag {4.25}
$$

Next, we set 

$$
w _ {N} (t) = - \int_ {t} ^ {1} \pi_ {N} ^ {\alpha , \beta - 2, \lambda} \{\partial_ {\tau} P \} (\tau) \mathrm{d} \tau .
$$

Obviously, $w _ { N } ( 1 ) = 0$ . Writing $\begin{array} { r } { \pi _ { N } ^ { \alpha , \beta - 2 , \lambda } \{ \partial _ { \tau } P \} ( \tau ) = \tau ^ { \frac { \beta - \lambda } { 2 } - 1 } \sum _ { k = 0 } ^ { N } \tilde { c } _ { k } ( \log \tau ) ^ { k } } \end{array}$ π α,β−2,λ {∂τ P}(τ ) = τ β , and integrating by parts, we find 

$$
\int_ {t} ^ {1} \pi_ {N} ^ {\alpha , \beta - 2, \lambda} \{\partial_ {\tau} P \} (\tau) \mathrm{d} \tau = \int_ {t} ^ {1} \tau^ {\frac {\beta - \lambda}{2} - 1} \sum_ {k = 0} ^ {N} \tilde {c} _ {k} (\log \tau) ^ {k} \mathrm{d} \tau = t ^ {\frac {\beta - \lambda}{2}} \sum_ {k = 0} ^ {N} \tilde {d} _ {k} (\log t) ^ {k},
$$

which implies $w _ { N } ( 0 ) = 0$ . Hence, $w _ { N } \in X _ { N }$ . On the other hand, we have 

$$
\partial_ {t} P - \partial_ {t} w _ {N} = (I - \pi_ {N} ^ {\alpha , \beta - 2, \lambda}) \partial_ {t} P.
$$

We can then derive from the above relations and Poincare inequality (4.16) that 

$$
\| P - w _ {N} \| _ {s} \leq c \| \partial_ {t} P - \partial_ {t} w _ {N} \| \leq c \| \partial_ {t} P - \partial_ {t} w _ {N} \| _ {\chi^ {\alpha , \lambda}} = c \| \partial_ {t} P - \pi_ {N} ^ {\alpha , \beta - 2, \lambda} \partial_ {t} P \| _ {\chi^ {\alpha , \lambda}}.
$$

Finally, combing the above and (4.25), and applying Theorems 3.3 and 3.4, we obtain the desired result. - 

Remark 4.5 As in Corollaries 2.1 and 3.1, we can show that for $P = t ^ { r }$ and $Q = t ^ { q }$ , the estimate (4.21) leads to exponential convergence rate if r and q are within certain range. 

# 4.4 Time-fractional diffusion equations

As a final application we consider the time-fractional diffusion equation 

$$
{ } _ { 0 } ^ { C } \mathrm{D} _ { t } ^ { v } u ( x , t ) - \Delta u ( x , t ) = f ( x , t ) , \quad x \in \Omega , t \in ( 0 , T ) , \tag {4.26}
$$

where $\varOmega$ be a bounded domain in $\mathbb { R } ^ { d } ( d = 1 , 2 , 3 )$ with suitable initial and boundary conditions. It is clear that the solution of the above equation will exhibit weak singularities at $t = 0$ so it is appropriate to use GLOFs for the time variable. As for the space variables, any consistent approximation can be used. The resulting linear system can be efficiently solved by using a matrix-diagonalization method Shen (1994); Shen et al. (2011). 

As a specific example we consider $\varOmega = ( - 1 , 1 ) ^ { 2 }$ with the following initial and boundary conditions: 

$$
u (x, t) = 0, \quad x \in \partial \Omega , t \in (0, T), \tag {4.27}
$$

$$
u (x, 0) = 0, \quad x \in \Omega , \tag {4.28}
$$

and we use a Legendre–Galerkin method (Shen, 1994) for the space variables. 

Let $N _ { t } , N _ { x }$ be respectively the degree of freedom of GLOFs in time and Legendre polynomials in each spatial direction. 

In the first test we choose the exact solution to be $u = ( t ^ { \mu } + t ^ { 2 \mu } )$ sin $\mathsf { \Omega } _ { 1 } ( \pi x _ { 1 } )$ sin(π x ), which is smooth in space, but has typical weak singularity in time. The errors in L2-norm with different ν are plotted in Fig. 9. We observe that the errors converge exponentially w.r.t. $N _ { t }$ and $N _ { x }$ . 

In the second test we take $\boldsymbol { f } = e ^ { x _ { 1 } x _ { 2 } t }$ . The explicit form of the exact solution is unknown, but is expected to be weakly singular at $\mathrm { \Delta t } = 0$ . We used a fine mesh to compute a reference solution and plotted the convergence rate in Fig. 10. Again, exponential convergence rates are observed for both $N _ { t }$ and $N _ { x }$ . 

Remark 4.6 In recent years numerous methods using the GJFs are developed for solving some model FDEs $\textstyle { \mathcal { L } } u = f , \quad t \in ( - 1 , 1 )$ , with $\mathcal { L } = \mathbf { \Phi } _ { - 1 } ^ { C } D _ { t } ^ { s } , \ \mathbf { o r } \mathbf { \Phi } _ { - 1 } D _ { t } ^ { s }$ , or $( - \varDelta ) ^ { s }$ and suitable boundary conditions, whose solutions behave like $u ( t ) = ( 1 + t ) ^ { r } ( 1 - t ) ^ { s } h ( t )$ with known $r , s \ge 0$ and smooth $h ( t )$ , see Zayernouri & Karniadakis (2013); Chen et al. (2016); Mao et al. (2016); Acosta et al. (2018); Hao & Zhang (2020) and the references therein. Specifically, one can define $\mathrm { G J F s } \omega P _ { n }$ such that 

$$
\mathcal {L} (\omega P _ {n}) = c _ {n} P _ {n}, \quad \mathcal {L} = \frac {C}{- 1} D _ {t} ^ {s}, \text {   or   } - 1 D _ {t} ^ {s}, \text {   or   } (- \Delta) ^ {s}, \tag {4.29}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/fa6561e66b508e00c240e23aa8203462245f6312ccbaeb5a221d9fa64cf39fa4.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/d5a61820d0a6ebb5c7f8e8cc62f8c57426288b8e077632ea5635af2abb486ca0.jpg)



Fig. 9. $u = ( t ^ { \mu } + t ^ { 2 \mu } )$ sin(π x1) sin(π x2), μ = 0.6, T = 1.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/07a590d95530356e866d7527bca9dd235016a0ead256c43f97595ca3091ab220.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-05-25/8aba2d3a-e5c4-458e-b43d-c17b10bd9c8a/9e294bb925f208312bb15cb1003b8cb26d606636f4d3b600932f852ffb9e2265.jpg)



Fig. 10. $f = e ^ { x _ { 1 } x _ { 2 } t } ;$ , T = 12 .


where Jacobi weight $\omega = ( 1 + t ) ^ { a } ( 1 - t ) ^ { b }$ and $P _ { n }$ is the related Jacobi polynomials. Obviously, the GJFs with $a = r , \ b = s$ can approximate the solution of the model problems $\mathcal { L } u = f$ very well. However, as is shown in (4.6), the singularity of the solution becomes much more complicated with just an additional zeroth-order term. In addition, the approximation results proved in Liu et al. (2019, 2020); Hao & Zhang (2020) show that the regularity is very limited when the weight function ω of the GJFs is different with the singularity of the solution. Therefore, methods using the GJFs can only lead to very limited accuracy for more general FDEs. On the other hand, GLOFs can handle functions with unknown endpoint singularity since they can approximate singular functions of the form $\textstyle { \sum _ { i } c _ { i } t ^ { r _ { i } } ( - \log t ) ^ { k _ { i } } }$ with exponential accuracy. This implies that the GLOFs can be used to obtain accurate solutions for a large class of fractional problems, such as equations with multiple different fractional derivatives. 

# 5. Concluding remarks

We constructed in this paper two new classes of orthogonal functions, the LOFs and the GLOFs, by applying a log mapping to the Laguerre functions. We developed basic approximation theory for these new orthogonal functions. The approximate results reveal that the new orthogonal functions are particularly suitable for functions that have weak singularities at one endpoint. In particular, for functions involving one or multiple terms of $t ^ { \alpha }$ with α in an adjustable range, its approximation by the new orthogonal functions will converge exponentially, as opposed to a low algebraic rate if usual orthogonal polynomials are used. 

As applications, we considered several typical FDEs whose solutions usually exhibit weak singularities at one endpoint. By using the GLOFs as basis functions, we constructed Galerkin methods for solving these FDEs and derived corresponding error analysis, which reveals the fact that exponential convergence rate can be achieved even if the solution is weakly singular at one endpoint. We provided ample numerical results to show that our methods based on GLOFs are very effective for problems with solutions having weak singularities at one endpoint, such as the cases in many FDEs. In particular, a special case of the GLOFs introduced in this paper has been used in Chen et al. (2020) to develop a very efficient and accurate spectral-Galerkin method (in the time direction) for solving the time-fractional subdiffusion equations. 

The methods presented in this paper are limited to problems with singularities at one endpoint. To deal with problems having singularities at both endpoints, one can use a two-domain approach with GLOFs on each subdomain, or to construct new classes of orthogonal functions that are suitable for problems having singularities at both endpoints. This topic will be the subject of a future study. 

# Funding

Postdoctoral Science Foundation of China (BX20180032, 2019M650459 to S.C.); Natural Science Foundation of China (11801235 to S.C.); Natural Science Foundation of the Jiangsu Higher Education Institutions of China (BK20181002 to S.C.); IN 47907-1957, USA; Natural Science Foundation (DMS-2012585 to J.S.); Air Force Office of Scientific Research (FA9550-20-1-0309 to J.S.). 

# Acknowledgements

S.C. would like to thank Professor Lilian Wang for many useful suggestions and enlightening discussions on the subject of 3.4 during his visit at Nanyang Technological University. 

# References



Acosta, G., Borthagaray, J. P., Bruno, O. & Maas, M. (2018) Regularity theory and high order numerical methods for the (1d)-fractional Laplacian. Math. Comput., 87, 1821–1857. 





Babuška, I. & Banerjee, U. (2012) Stable generalized finite element method (SGFEM). Comput. Methods Appl. Mech. Enger., 201, 91–111. 





Bernardi, C. & Maday, Y. (1997) Spectral methods. Handbook Numer. Anal., 5, 209–485. 





Boyd, J. P. (1986) Polynomial series versus sinc expansions for functions with corner or endpoint singularities. J. Comput. Phys., 64, 266–270. 





Boyd, J. P. (1987a) Exponentially convergent Fourier–Chebyshev quadrature schemes on bounded and infinite intervals. J. Sci. Comput., 2, 99–109. 





Boyd, J. P. (1987b) Spectral methods using rational basis functions on an infinite interval. J. Comput. Phys., 69, 112–142. 





Boyd, J. P. (2001) Chebyshev and Fourier Spectral Methods. New York: Dover. 





Canuto, C., Hussaini, M. Y., Quarteroni, A. & Zang, T. A. (2006) Spectral Methods: Fundamentals in Single Domains. Berlin-Heidelberg: Springer. 





Cao, Y., Herdman, T. & Xu, Y. (2003) A hybrid collocation method for Volterra integral equations with weakly singular kernels. SIAM J. Numer. Anal., 41, 364–381. 





Chen, S. (2020) High-accuracy numerical approximations to several singularly perturbed problems and singular integral equations by enriched spectral Galerkin methods. J. Math. Study., 53, 143–158. 





Chen, S. & Shen, J. (2018) Enriched spectral methods and applications to problems with weakly singular solutions. J. Sci. Comput., 77, 1468–1489. 





Chen, S. & Shen, J. (2020) An efficient and accurate numerical method for the spectral fractional Laplacian equation. J. Sci. Comput., 82, 1–25. 





Chen, S., Shen, J. & Wang, L. L. (2016) Generalized Jacobi functions and their applications to fractional differential equations. Math. Comput., 85, 1603–1638. 





Chen, S., Shen, J., Zhang, Z. M. & Zhou, Z. (2020) A spectrally accurate approximation to subdiffusion equations using the log orthogonal functions. SIAM J. Sci. Comput., 42, A849–A877. 





Diethelm, K. (2010) The Analysis of Fractional Differential Equations. Lecture Notes in Math., vol. 2004. Berlin: Springer. 





Ervin, V. J., Heuer, N. & Roop, J. P. (2007) Numerical approximation of a time dependent, nonlinear, spacefractional diffusion equation. SIAM J. Numer. Anal., 45, 572–591. 





Ervin, V. J. & Roop, J. P. (2006) Variational formulation for the stationary fractional advection dispersion equation. Numer. Meth. Part. Diff. Eqn., 22, 558–576. 





Ervin, V. J. & Roop, J. P. (2007) Variational solution of fractional advection dispersion equations on bounded domains in Rd. Numer. Meth. Part. Diff. Eqn., 23, 256–281. 





Fries, T. P. & Belytschko, T. (2010) The extended/generalized finite element method: an overview of the method and its applications. Int. J. Numer. Meth. Eng., 84, 253–304. 





Gottlieb, D. & Orszag, S. A. (1977) Numerical Analysis of Spectral Methods: Theory and Applications. Number 26 in CBMS-NSF Regional Conference Series in Applied Mathematics. Philadelphia: SIAM. 





Gradshteyn, I. S. & Ryzhik, I. M. (2007) In table of integrals, series, and products. Math. Comp., 20. 





Guo, B. Y. (1998) Spectral Methods and Their Applications. Singapore: World Scientific. 





Guo, B. Y., Wang, L. L. & Wang, Z. Q. (2006) Generalized Laguerre interpolation and pseudospectral method for unbounded domains. SIAM J. Numer. Anal., 43, 2567–2589. 





Hao, Z. P. & Zhang, Z. Q. (2020) Optimal regularity and error estimates of a spectral Galerkin method for fractional advection–diffusion–reaction equations. SIAM J. Numer. Anal., 58, 211–233. 





Jin, B. T., Lazarov, R. & Zhou, Z. (2013) Error estimates for a semidiscrete finite element method for fractional order parabolic equations. SIAM J. Numer. Anal., 51, 445–466. 





Li, C. P., Zeng, F. H. & Liu, F. W. (2012) Spectral approximations to the fractional integral and derivative. Fract. Calc. Appl. Anal., 15, 383–406. 





Li, X. & Xu, C. (2009) A space-time spectral method for the time fractional diffusion equation. SIAM J. Numer. Anal., 47, 2108–2131. 





Li, X. & Xu, C. (2010) Existence and uniqueness of the weak solution of the space-time fractional diffusion equation and a spectral method approximation. Commun. Comput. Phys., 8, 1016. 





Liu, W. J., Wang, L. L. & Li, H. Y. (2019) Optimal error estimates for Chebyshev approximations of functions with limited regularity in fractional Sobolev-type spaces. Math. Comput., 88, 2857–2895. 





Liu, W. J. , Wang, L. L. & Wu, B. Y. (2020) Optimal error estimates for Legendre approximation of singular functions with limited regularity. arXiv preprint arXiv:2006.00667. 





Mao, Z. P., Chen, S. & Shen, J. (2016) Efficient and accurate spectral method using generalized Jacobi functions for solving Riesz fractional differential equations. Appl. Numer. Math., 106, 165–181. 





Meerschaert, M. M. & Tadjeran, C. (2004) Finite difference approximations for fractional advection–dispersion flow equations. J. Comput. Appl. Math., 172, 65–77. 





Morin, P., Nochetto, R. & Siebert, K. G. (2002) Convergence of adaptive finite element methods. SIAM Rev., 44, 631–658. 





Podlubny, I. (1999) Fractional Differential Equations: An Introduction to Fractional Derivatives, Fractional Differential Equations to Methods of Their Solution and Some of Their Applications. USA: Academic Press. 





Samko, S. G., Kilbas, A. A. & Maricev ˇ , O. I. (1993) Fractional Integrals and Derivatives. Yverdon: Gordon and Breach Science Publ. 





Schumack, M. R., Schultz, W. W. & Boyd, J. P. (1991) Spectral method solution of the Stokes equations on nonstaggered grids. J. Comput. Phys, 94, 30–58. 





Shen, J. (1994) Efficient spectral-Galerkin method i. Direct solvers of second-and fourth-order equations using Legendre polynomials. SIAM J. Sci. Comput., 15, 1489–1505. 





Shen, J., Tang, T. & Wang, L. L. (2011) Spectral Methods: Algorithms, Analysis and Applications. vol. 41 of Series in Computational Mathematics. Berlin: Heidelberg. 





Shen, J. & Wang, L. L. (2004) Error analysis for mapped Legendre spectral and pseudospectral methods. SIAM J. Numer. Anal., 42, 326–349. 





Shen, J. & Wang, Y. (2016) Müntz–Galerkin methods and applications to mixed Dirichlet–Neumann boundary value problems. SIAM J. Sci. Comput., 38, A2357–A2381. 





Strang, G. & Fix, G. J. (1973) An Analysis of the Finite Element Method. Prentice-Hall Series in Automatic Computation. Englewood Cliffs, N. J: Prentice-Hall, Inc. 





Sun, Z. Z. & Wu, X. N. (2006) A fully discrete difference scheme for a diffusion-wave system. Appl. Numer. Math., 56, 193–209. 





Szego, G. (1975) Orthogonal Polynomials, Math. Soc. Colloq. Publ. Series. Math. Soc., vol. 23, 2nd edn., Providence, RI: Amer. Math. Soc. 





Wang, L. L. & Shen, J. (2005) Error analysis for mapped Jacobi spectral methods. J. Sci. Comput., 24, 183–218. 





Zayernouri, M. & Karniadakis, G. E. (2013) Fractional Sturm–Liouville eigen-problems: theory and numerical approximation. J. Comput. Phys., 252, 495–517. 





Zhao, X. D., Wang, L. L. & Xie, Z. Q. (2013) Sharp error bounds for Jacobi expansions and Gegenbauer–Gauss quadrature of analytic functions. SIAM J. Numer. Anal., 51, 1443–1469. 



# A. Some properties of Laguerre polynomials

The three-term recurrence 

$$
\mathcal {L} _ {0} ^ {(\alpha)} (y) = 1, \qquad \mathcal {L} _ {1} ^ {(\alpha)} = - y + \alpha + 1,
$$

$$
\mathcal {L} _ {n + 1} ^ {(\alpha)} (y) = \frac {2 n + \alpha + 1 - y}{n + 1} \mathcal {L} _ {n} ^ {(\alpha)} (y) - \frac {n + \alpha}{n + 1} \mathcal {L} _ {n - 1} ^ {(\alpha)} (y). \tag {A.1}
$$

Sturm–Liouville problem 

$$
y ^ {- \alpha} e ^ {y} \partial_ {y} \left(y ^ {\alpha + 1} e ^ {- y} \partial_ {y} \mathscr {L} _ {n} ^ {(\alpha)} (y)\right) + n \mathscr {L} _ {n} ^ {(\alpha)} (y) = 0. \tag {A.2}
$$

Derivative relations 

$$
\mathcal {L} _ {n} ^ {(\alpha)} (y) = \partial_ {y} \mathcal {L} _ {n} ^ {(\alpha)} (y) - \partial_ {y} \mathcal {L} _ {n + 1} ^ {(\alpha)} (y), \tag {A.3}
$$

$$
y \partial_ {y} \mathcal {L} _ {n} ^ {(\alpha)} (y) = n \mathcal {L} _ {n} ^ {(\alpha)} (y) - (n + \alpha) \mathcal {L} _ {n - 1} ^ {(\alpha)} (y), \tag {A.4}
$$

$$
\partial_ {y} \mathscr {L} _ {n} ^ {(\alpha)} (y) = - \mathscr {L} _ {n - 1} ^ {(\alpha + 1)} (y) = - \sum_ {k = 0} ^ {n - 1} \mathscr {L} _ {k} ^ {(\alpha)} (y). \tag {A.5}
$$

Laguerre–Gauss quadrature Let $\{ y _ { j } ^ { ( \alpha ) } \} _ { j = 0 } ^ { N }$ be the zeros of $\mathcal { L } _ { N + 1 } ^ { ( \alpha ) } ( y )$ , then the associated weights are 

$$
\omega_ {j} ^ {(\alpha)} = \frac {\Gamma (N + \alpha + 1)}{(N + \alpha + 1) N !} \frac {y _ {j} ^ {(\alpha)}}{[ \mathscr {L} _ {N} ^ {(\alpha)} (y _ {j} ^ {(\alpha)}) ] ^ {2}}, \quad 0 \leq j \leq N, \tag {A.6}
$$

the quadrature formula is 

$$
\int_ {\mathbb {R} ^ {+}} p (y) y ^ {\alpha} e ^ {- y} \mathrm{d} y = \sum_ {j = 0} ^ {N} p (y _ {j} ^ {(\alpha)}) \omega_ {j} ^ {(\alpha)}, \quad \forall p \in \mathcal {P} _ {2 N + 1} ^ {y}. \tag {A.7}
$$