# 端点对数奇异函数 $x^r(-\ln x)^k$ 的谱后处理与指数精度恢复研究报告

## 执行摘要

本文讨论的问题是：对区间 $[0,1]$ 上具有端点对数奇异性的函数
$$
f(x)=x^r(-\ln x)^k,\qquad r>0,\ k\in\mathbb Z_+,
$$
若只知道其傅里叶系数或谱配置法配点值，如何设计后处理使其重新获得“解析函数情形”下的指数精度。核心结论是：**纯粹照搬经典 Gegenbauer 重构或 Chen–Shu 的代数型映射，并不能直接把 $x^r(-\ln x)^k$ 在整个区间上变成解析函数，因此一般不能单独保证全区间指数精度；要实现真正的指数恢复，必须把“已知的对数奇异因子”显式剥离，或者改用对数坐标 $t=-\ln x$ 的半无限区间重构，再对解析残差做谱后处理。** 这一点与 Gottlieb–Shu 框架中的“**存储在一个基底中的高阶信息，需要在另一个更适配的基底中读取出来**”的思想完全一致，但对数奇异要求“读取基底”进一步与奇异结构精确匹配。citeturn19view1turn19view2turn3view1turn12view0turn17search2

从理论上看，$x^r(-\ln x)^k$ 在 $x=0$ 属于典型的**代数–对数端点奇异**。在 Chebyshev/移位 Chebyshev 展开中，其系数只按代数–对数律衰减；若 $r\notin\mathbb N_0$，典型量级是 $n^{-(2r+1)}(\ln n)^k$，若 $r\in\mathbb N_0$，由于 $\sin(\pi r)$ 的零点抵消，最高对数幂会下降一阶，变为 $n^{-(2r+1)}(\ln n)^{k-1}$。这类结论与 Wang 关于单对数奇异的精确渐近公式、以及 Xiang–Liu 关于 $\ln^\mu$ 端点奇异的最优衰减率是一致的。对 Fourier 系数，端点 $x=0$ 的代数–对数奇异贡献可由 Sidi 的一般端点振荡积分渐近理论给出，但在 $[0,1]$ 的普通 Fourier 系数情形下，还必须同时考虑 $x=1$ 处解析端的周期拼接效应；后者常常会与 $x=0$ 的奇异贡献竞争，因而在 Fourier 情形下必须额外做“右端点中和”或基底变换。citeturn10view0turn15view0turn16view1turn12view0turn7search0

因此，本文给出的首选方案不是“纯 Gegenbauer”，而是**奇异因子剥离 + 解析残差谱重构**。若已知奇异结构恰为 $x^r(-\ln x)^k$，则写成
$$
f(x)=w(x)g(x),\qquad w(x)=x^r(-\ln x)^k,
$$
对当前目标函数甚至有 $g\equiv 1$。此时只要从配点值中先除去 $w(x)$，或从 Fourier 系数中通过一个正则化连接矩阵把数据投影到 $w(x)$ 加权基底，剩余的 $g$ 就是解析函数；对 $g$ 做普通 Chebyshev、Gegenbauer 或 Jacobi 重构即可恢复指数精度，最后再乘回 $w(x)$。若不想显式除以 $w$，另一条同样强有力的路线是做对数坐标变换 $t=-\ln x$，则
$$
f(e^{-t})=e^{-rt}t^k,
$$
奇异性被“平直化”为半无限区间上的指数衰减解析函数，此时配合 Laguerre / rational-Chebyshev / Adcock–Richardson 型变换后处理，同样可获得指数或近指数恢复。citeturn17search0turn17search2turn18search4turn20view1turn20view2

就文献脉络而言，Gottlieb–Shu 1995 及其 1997 SIAM Review 综述奠定了“Gegenbauer 重构 + $\lambda,m\sim N$”的基本框架，证明了从 Fourier 系数、Galerkin 系数或 collocation 点值中恢复分段解析函数局部指数精度的机制。Chen–Shu 2014 将此法扩展到配点值已知、端点存在**有理代数型**奇异 $a(x)+b(x)(1+x)^{p/q}$ 的情形；Chen–Shu 2015 又把它扩展到 Fourier 系数情形，并展示了对噪声的相对鲁棒性。但这些工作都依赖于一个关键步骤：通过**构造性代数映射**把端点奇异转成 $y$-变量中的解析函数；这一条件对纯对数奇异并不满足，因此本报告建议的新研究方向，是把 Chen–Shu 的“构造性映射思想”与“奇异因子剥离”或“对数坐标重构”结合起来，建立针对 $x^r(-\ln x)^k$ 的严格理论。citeturn3view1turn19view1turn19view2 fileciteturn0file0 fileciteturn0file1

## 背景与动机

Gottlieb–Shu 对 Gibbs 现象的核心认识，不是“截断谱和本身坏了”，而是“**错误的读取方式**导致了坏的点值重构”。他们在综述中明确指出：Fourier 系数本身包含足够的信息，问题出在用原始 Fourier 偏和去读取一个非周期或分段解析函数；若改用 Gegenbauer 基底，并令参数 $\lambda$ 与保留项数 $m$ 随 $N$ 线性增长，就可以把局部解析区域中的点值恢复到与解析情形同量级的指数精度。这个“storage vs. retrieval”的思路，就是所谓“谱系数与配点值中隐含高阶谱信息”的精确数学表述。citeturn19view1turn19view2turn19view3

在 collocation 情形，Gottlieb–Shu 1995 证明：若已知 Gauss 或 Gauss–Lobatto 点上的配点值，那么可以从插值多项式出发，重构出在解析子区间上一致指数收敛的 Gegenbauer 近似。这个结论后来被 1997 年 SIAM Review 系统化为两步：第一步从给定系数/配点值近似出前 $m$ 个 Gegenbauer 系数，第二步再由这些系数做 Gegenbauer 正则化；两步误差分别称为 truncation error 与 regularization error，二者都可在适当参数标定下指数衰减。citeturn3view1turn19view0turn19view2

Chen–Shu 2014/2015 的贡献，在于把这个框架从“分段解析但边界跳跃”推进到“开区间解析、端点有奇异”。2014 年论文处理的是配点值已知的情形，目标函数型如
$$
f(x)=a(x)+b(x)(1+x)^s,\qquad 0<s=p/q<1,
$$
它通过代数映射
$$
(2^{q-1}(1+x))^{1/q}=1+y
$$
把 $x$-域中的分数次幂奇异变成 $y$-域中的解析函数，然后再做 Gegenbauer 重构，并证明只要 $\lambda=\alpha N$、$m=\gamma\lambda$，则 truncation error 与 regularization error 都指数衰减。2015 年论文把这个想法进一步推广到了前 $2N+1$ 个 Fourier 系数已知的情形，并展示了对含噪数据的数值鲁棒性。fileciteturn0file0 fileciteturn0file1

但是，Chen–Shu 的构造性映射有一个非常强的适用条件：必须能把端点奇异通过**代数型变量变换**推成新变量中的解析对象。对 $(1+x)^{p/q}$ 这成立，因为分数幂在适当开 $q$ 次根后会变成整数幂；而对 $(-\ln x)^k$ 或 $\ln(1+x)$ 这不成立，因为任何有限次代数代换都只会把对数奇异变成“另一个对数奇异”，不会凭空消失。因此，对 $x^r(-\ln x)^k$ 而言，**纯 Chen–Shu 映射并不足以直接恢复全区间指数精度**。要么先剥离对数奇异因子，要么改用 $t=-\ln x$ 这样的非代数变换。fileciteturn0file0 fileciteturn0file1 citeturn17search0turn17search2

## 目标函数特性与谱系数衰减

函数
$$
f(x)=x^r(-\ln x)^k
$$
在每个 $x_0\in(0,1)$ 的邻域内都是实解析的；其复解析延拓可在去掉负实轴的割线平面上选定分支，因此从复变观点看，$x=0$ 是一个典型的**代数–对数分支点**。若 $r\notin\mathbb N_0$，奇异性同时来自 $x^r$ 的代数分支与 $\ln x$ 的对数分支；若 $r=m\in\mathbb N_0$，函数本身在 $x=0$ 仍趋于 $0$，但第 $m$ 阶导数会出现对数奇异，因此更准确地说是“**第 $m$ 阶导数的对数奇异**”。这正是 Sidi 2009 关于端点 algebraic-logarithmic singularities 的框架所覆盖的类型。citeturn7search0turn8view1

把区间线性移到 $[-1,1]$，令 $t=2x-1$，并记移位 Chebyshev 系数为
$$
f(x)=\sum_{n\ge0} c_n^*\,T_n^*(x),\qquad T_n^*(x)=T_n(2x-1).
$$
则
$$
f\!\left(\frac{1+t}{2}\right)=2^{-r}(1+t)^r\Bigl(\ln\frac{2}{1+t}\Bigr)^k.
$$
于是 $t=-1$ 是唯一的端点奇异点，属于 Xiang–Liu 所研究的 $(1+t)^\delta \ln^\mu(1+t)\,g(t)$ 型。由其最优衰减率理论可直接推出：若 $r\notin\mathbb N_0$，则
$$
|c_n^*|=O\!\bigl((\ln n)^k\,n^{-(2r+1)}\bigr),
$$
而若 $r\in\mathbb N_0$，则最高对数幂下降一阶，
$$
|c_n^*|=O\!\bigl((\ln n)^{k-1}n^{-(2r+1)}\bigr).
$$
这两条就是本问题最基本的“不可直接指数收敛”的谱障碍。citeturn15view0turn16view1turn14view2

进一步地，Wang 对单对数端点奇异给出了 Chebyshev 系数的显式渐近展开。将其代数情形的主项对参数 $r$ 连续求导 $k$ 次，可得到本问题的一个强形式推论：对非整数 $r$，移位 Chebyshev 系数的主导项满足
$$
c_n^*
\sim
(-1)^{n+1}\frac{2^{k+1-2r}\Gamma(2r+1)\sin(\pi r)}{\pi}\,
n^{-(2r+1)}(\ln n)^k,
\qquad n\to\infty,
$$
其后续修正项是同阶但对数次数更低的多项式；若 $r=m\in\mathbb N_0$，则因为 $\sin(\pi r)$ 在整数点具有简单零点，主导项中的 $(\ln n)^k$ 被消掉，首项降为 $n^{-(2m+1)}(\ln n)^{k-1}$。这一步是由 Wang 的单对数精确展开与参数微分得到的推论，因此应看作“由一般定理推出的针对本函数族的特化公式”。citeturn10view0turn15view0turn16view1

对 Fourier 系数也有类似结论。令
$$
\widehat f_n=\int_0^1 x^r(-\ln x)^k e^{-2\pi i n x}\,dx.
$$
Sidi 2010 给出了有限区间上带端点 algebraic-logarithmic singularities 的 Fourier 振荡积分完整渐近框架。应用到 $x=0$ 的奇异端点，其首要贡献可写成
$$
\widehat f_n^{(0)}
\sim
(-1)^k\partial_r^k\!\Bigl(
e^{-i\pi(r+1)/2}\Gamma(r+1)(2\pi n)^{-(r+1)}
\Bigr),
$$
因此在量级上有
$$
|\widehat f_n^{(0)}|
\asymp
n^{-(r+1)}(\ln n)^k
$$
（整数 $r$ 时同样会出现对数次数下降一阶的现象）。但这里有一个 Fourier 情形独有的麻烦：$[0,1]$ 上的普通 Fourier 系数同时还会收到 $x=1$ 解析端的贡献；由于 $f(x)$ 在 $x=1$ 处是解析的并且零点阶数恰为 $k$，该端点会产生 $O(n^{-(k+1)})$ 的常规边界项。因此，**未经中和的 Fourier 系数，其主导衰减未必由 $x=0$ 的奇异端控制，而常常由 $x=1$ 的周期拼接效应控制**。这正是 Fourier 后处理比 collocation 后处理更难的核心原因。citeturn12view0

这也解释了为什么：如果目标是“恢复全区间指数精度”，仅靠改写原 Fourier 偏和并不足够；必须先在模型层面识别并剥离 $x=0$ 的已知奇异因子，同时对 $x=1$ 的光滑端做周期化中和，或者干脆改用在 $[0,1]$ 更自然的 Chebyshev/shifted Jacobi 数据结构。换句话说，对 $x^r(-\ln x)^k$ 来说，**正确的后处理对象不是原函数 $f$，而是“去奇异后的解析残差” $g$**。citeturn19view3turn12view0turn16view1

## 关键文献综述

下表压缩列出与本问题最相关的原始论文与综述。表中“适用性”一栏专门针对 $x^r(-\ln x)^k$ 这一类端点对数奇异进行评价。

| 文献 | 方法核心 | 主要结论 | 对 $x^r(-\ln x)^k$ 的适用性 |
|---|---|---|---|
| Gottlieb–Shu 1995, *On the Gibbs Phenomenon V* | 从配点值近似 Gegenbauer 系数，再做 Gegenbauer 正则化 | 对分段解析函数，可由 Gauss/Gauss–Lobatto 配点值恢复解析子区间上一致指数精度。 citeturn3view1 | 奠定框架，但原文不含对数型端点奇异。 |
| Gottlieb–Shu 1997, *SIAM Review* | 系统总结 Fourier / Gegenbauer / collocation / Galerkin 的 Gibbs 消除 | 明确提出“存储在一个基底、在另一个基底读取”的机制；$\lambda,m$ 须与 $N$ 同阶增长。 citeturn19view1turn19view2 | 是方法论基础，但不直接覆盖对数奇异。 |
| Chen–Shu 2014 | 构造性代数映射 + collocation Gegenbauer 重构 | 端点代数奇异 $a+b(1+x)^{p/q}$ 的配点值可恢复指数精度；误差分解为 truncation 与 regularization 两部分，并都可指数衰减。 fileciteturn0file0 | 对有理代数幂非常合适；对纯对数奇异不直接适用。 |
| Chen–Shu 2015 | 从 Fourier 系数出发的映射型 Gegenbauer 重构 | 对端点导数无界的代数奇异函数，可从 Fourier 系数恢复指数精度；数值上对噪声有一定鲁棒性。 fileciteturn0file1 | 对 Fourier 情形最接近本问题，但仍要求代数型可映射奇异。 |
| Adcock–Richardson 2014 | 端点奇异的变量变换近似 | 证明传统映射分辨率次优，并提出新变换以改进端点奇异函数的近似与分辨率；可获得根指数或用户可控精度。 citeturn17search0turn17search2 | 对“对数坐标/非代数映射”路线很关键。 |
| Wang 2014 | Chebyshev 系数与 Clenshaw–Curtis 对端点 algebraic / algebraic-log singularity 的渐近分析 | 给出单对数端点奇异的 Chebyshev 系数渐近展开，并由此得到相应积分公式的最优收敛率。 citeturn10view0turn9view3 | 可直接为 $k=1$ 提供精确系数渐近；对一般 $k$ 可做参数微分推广。 |
| Sidi 2009 / 2010 | Legendre 系数与 Fourier 振荡积分的端点 algebraic-log 渐近分析 | 建立了端点 algebraic-log singularities 下 Legendre 系数与 Fourier 积分的完整渐近理论。 citeturn8view1turn12view0 | 为本问题的 Fourier 渐近、以及“为什么普通 Fourier 不干净”提供理论支撑。 |
| Xiang–Liu 2020 | Jacobi / Gegenbauer / Chebyshev 展开的最优衰减率 | 对有限正则与端点对数奇异给出最优系数衰减率与投影误差阶；指出 Jacobi/Gegenbauer 可比 Chebyshev 更快。 citeturn13search1turn15view0turn15view2 | 直接给出 $(1\pm x)^\delta\ln^\mu(1\pm x)$ 的最优衰减率，是本问题最关键的现代参考。 |
| Zhang–Boyd 2023 | 弱端点奇异下 Chebyshev / difference / quadratic-factor 基底比较 | 说明基底选择会显著改变系数衰减、误差分布及边界层结构。 citeturn4search0turn3view4 | 虽主攻弱代数奇异，但对“先抽出奇异因子再逼近解析残差”的思路很有启发。 |
| Guo–Shen–Wang 2006；Chen–Shen–Wang 2016；Mao–Chen–Shen 2016 | generalized Jacobi / generalized Jacobi functions | 针对端点奇异或分数阶 PDE，构造与边界/奇异度匹配的自然基底，可获得稀疏或良态线性系统，并在光滑数据下达到真正的谱精度。 citeturn20view1turn18search4turn20view2 | 为“加权基底/奇异富集基”路线提供强有力的替代框架。 |

从这些文献的共同结论可以看出，**真正能恢复指数精度的机制，从来不是“简单重和原级数”，而是“构造一个让剩余对象变解析的表示”**。对代数奇异，这个表示可以是 Chen–Shu 的代数映射；对对数奇异，更自然的表示则是“已知奇异因子的显式剥离”或“对数坐标下的半无限区间表示”。这一判断同时符合 Gottlieb–Shu 的信息重读思想、Adcock–Richardson 的非代数变换思路、以及 generalized Jacobi functions 的奇异匹配基底理念。citeturn19view3turn17search0turn18search4turn20view2

## 技术问题与挑战

第一个挑战是**纯对数奇异无法被有限次代数映射消除**。Chen–Shu 2014/2015 之所以能证明指数恢复，是因为 $x^s$ 在 $x\mapsto y$ 的适当有理次幂变换后可变成解析项；但 $(-\ln x)^k$ 无论经过多少次代数代换，仍然会保留 $\ln(\cdot)$ 结构。这意味着：若坚持在有限区间上只用“多项式基 + 代数映射”，则全区间解析性通常不会出现， ফলে classical Gegenbauer regularization 的理论闭环断掉。fileciteturn0file0 fileciteturn0file1

第二个挑战是**Fourier 数据的双端点耦合**。对 $[0,1]$ 上的函数，普通 Fourier 系数并不只反映 $x=0$ 的奇异，还会受到 $x=1$ 处周期拼接边界项的影响。对于 $x^r(-\ln x)^k$，$x=1$ 虽然不是奇异点，但因为函数在该端恰有 $k$ 阶零点，周期扩展后往往只剩有限阶光滑性；因此若直接在 Fourier 系数上做 Gegenbauer 或最小二乘拟合，拟合对象会被两种不同机制共同驱动，参数识别与理论估计都明显复杂化。citeturn12view0

第三个挑战是**参数选择与条件数**。Gottlieb–Shu 与 Chen–Shu 理论都要求 $\lambda\sim N$、$m\sim N$，但这同时意味着 Gegenbauer 基底的条件数会迅速增大。Chen–Shu 2015 的数值结果已经表明：按“最优”参数选取可以把误差压得更小，但对噪声也更敏感；“线性规则”参数虽略保守，却往往更稳健。对本问题若再叠加奇异因子剥离、矩阵连接或对数映射，这种条件数–精度–噪声三角关系会更尖锐。fileciteturn0file1

第四个挑战是**数值下溢与局部放大**。在接近 $x=0$ 的区域，$w(x)=x^r(-\ln x)^k$ 既很小又变化很快，直接做 $g(x)=f(x)/w(x)$ 会放大测量误差与舍入误差。若原始数据来自含噪 Fourier 系数，或某种谱 PDE 数值解而非精确投影，这个问题会尤为突出。因此任何“奇异因子剥离”方法都必须伴随：对数尺度计算、局部加权最小二乘、端点点值舍弃/约束、以及正则化参数选择。citeturn12view0turn20view2

第五个挑战是**理论证明尚不完整**。现有可直接引用的严格定理，一类针对代数型端点奇异的 Gegenbauer/映射恢复，另一类针对对数奇异的 Jacobi/Chebyshev 系数衰减与投影误差；但把两者合并成“从 Fourier 系数或标准配点值出发、通过显式对数奇异剥离后，恢复全区间指数精度”的完整定理，目前公开文献中并没有直接现成版本。这正是本研究最值得推进的地方。citeturn15view0turn15view2turn12view0turn19view0

## 可行后处理策略与算法

从现有理论看，针对 $x^r(-\ln x)^k$ 最有希望形成完整理论与高质量数值方法的，不是单一方案，而是三条互补路线：**奇异因子剥离**、**对数坐标重构**、以及**加权/富集基底重构**。下面先给出比较结论，再给出推荐算法。

### 方法比较

| 方法 | 基本思想 | 预期收敛 | 稳定性 | 实现复杂度 |
|---|---|---:|---:|---:|
| 纯 Gegenbauer 重构 | 直接从原 Fourier/配点数据重构 | 一般仅代数–对数改进，难证全区间指数 | 对 $\lambda,m$ 敏感 | 中 |
| Chen–Shu 代数映射 | 先用代数映射消奇异，再做 Gegenbauer | 对纯对数奇异通常**不能**直接指数恢复 | 中等 | 中 |
| 奇异因子剥离 + 解析残差谱重构 | 先写 $f=w\,g$，对解析 $g$ 做谱重构，再乘回 $w$ | **若 $w$ 已知且 $g$ 解析，可达指数精度** | 端点需正则化，但总体最好 | 低到中 |
| 对数坐标 $t=-\ln x$ + Laguerre/rational-Chebyshev | 在半无限区间把奇异点平直化 | 指数或近指数，取决于基底与截断 | 好，但映射参数需调 | 中到高 |
| 加权 Jacobi / generalized Jacobi / singular basis | 直接在匹配奇异的基底中展开 | 若基底精确匹配，可达谱精度 | 好于纯 Gegenbauer | 中到高 |
| 局部减奇 + 分段谱/最小二乘 | 拟合若干奇异首项后重构其余部分 | 若只减有限项，多为改进的代数阶；减全奇异因子则可指数 | 依赖拟合质量 | 中 |
| Fourier 系数上的正则化连接矩阵 | 直接把 Fourier 数据反投影到 $w\times$ 多项式基 | 理论上可指数，实践上受病态与噪声限制 | 需 Tikhonov 或 TSVD | 高 |

综合来看，**对当前这个已知奇异结构的函数族，首选是“奇异因子剥离 + 解析残差谱重构”**。原因非常简单：一旦把 $w(x)=x^r(-\ln x)^k$ 剥掉，残差 $g$ 就是解析对象；此时恢复指数精度已经不是“困难的 Gibbs 问题”，而退化成“普通解析函数的谱逼近问题”。这比试图让 Gegenbauer 基底自己去“吸收”对数奇异要直接得多，也更容易得到证明。citeturn19view3turn16view1turn20view2

### 推荐方案

我建议把研究重点放在下面这个**混合型后处理框架**上：

```mermaid
flowchart TD
A[输入数据: Fourier系数或谱配点值] --> B{数据类型}
B -->|配点值| C[构造奇异因子 w(x)=x^r(-ln x)^k]
B -->|Fourier系数| D[右端点周期化中和 或 直接构造 Fourier-加权基连接矩阵]
C --> E[对 x>0 的样本做 g=f/w]
D --> E
E --> F[对解析残差 g 做 Chebyshev/Gegenbauer/Jacobi 重构]
F --> G[必要时在 x=0 附近做局部正则化或分区拼接]
G --> H[输出 \u007Ef(x)=w(x)\u007Eg(x)]
H --> I[参数扫描与误差评估]
```

这一框架有两个版本。

#### 配点值版本

若输入是 Chebyshev、Legendre 或其他正交 Gauss 型配点值 $f(x_j)$，则最简单的处理是直接做
$$
g(x_j)=\frac{f(x_j)}{x_j^r(-\ln x_j)^k},\qquad x_j>0.
$$
如果节点包含 $x=0$，则该点不参与除法，而作为硬约束 $f(0)=0$ 保留。对 $g$ 而言，只要原问题确实属于 $w(x)\times$ 解析函数类，则 $g$ 在 $[0,1]$ 上解析；于是普通移位 Chebyshev 插值、Gegenbauer 投影，甚至不需要 $\lambda\sim N$ 的重型参数，只用标准解析函数谱逼近即可达到指数收敛。此时再乘回 $w(x)$ 即得重构。由于真正困难的奇异性已被显式处理，这个版本最容易证明，并且能直接吸收 generalized Jacobi / weighted basis 的现有理论。citeturn20view1turn18search4turn20view2

一个干净的实现是：在解析残差层面直接使用移位 Chebyshev 基 $T_n^*(x)$ 或 Gegenbauer 基 $C_n^\lambda(2x-1)$。若原数据就是 Chebyshev 类节点，残差系数可用 DCT 级别的快速变换获得；若是一般 Gauss 节点，则可做稳定的加权最小二乘。数值上应在对数尺度计算权函数，
$$
\log w(x)=r\log x + k\log(-\log x),
$$
避免靠近 $x=0$ 时的下溢。citeturn18search4turn20view2

#### Fourier 系数版本

若输入是普通 Fourier 系数，则建议分两步。第一步做**端点中和**：在系数层面减去 $x=1$ 端解析部分的若干 Taylor 项对 Fourier 系数的贡献，或者先用一个温和滤波器把 Fourier 数据转成过采样物理点值。第二步再把剩余数据投影到
$$
\psi_m(x)=x^r(-\ln x)^k p_m(x)
$$
这样的加权谱基上，其中 $p_m$ 可取移位 Chebyshev / Gegenbauer / Jacobi 多项式。若选择“直接系数反投影”，则需求解
$$
\widehat f_n \approx \sum_{m=0}^{M} a_m \int_0^1 x^r(-\ln x)^k p_m(x)e^{-2\pi i n x}\,dx,
$$
这是一个稠密但结构化的连接矩阵问题；实际实现中需要 Tikhonov 正则化或截断奇异值分解。若选择“先转物理点值”，则可直接回到上一节的配点值算法。二者相比，前者理论上更纯粹，后者工程上更容易落地。citeturn12view0turn19view3turn16view1

### 伪代码

下面给出我认为最值得首先实现并写理论的两个伪代码。

**算法一：配点值后处理**

```text
Input:
    节点 x_j ∈ [0,1], 配点值 f_j = f(x_j), 已知奇异参数 r,k, 重构阶 M
Choose:
    解析残差基 p_m(x) = T_m^*(x) 或 C_m^λ(2x-1)

Step A:
    对所有 x_j > ε，计算
        g_j = f_j / exp(r log x_j + k log(-log x_j))
    若节点包含 x=0，则不做除法，仅保留约束 f(0)=0

Step B:
    用 DCT 或加权最小二乘求解析残差系数 a_m:
        g_M(x) = Σ_{m=0}^M a_m p_m(x)

Step C:
    输出
        f_M(x) = x^r(-log x)^k g_M(x)

Optional:
    在 x < x_c 的极小邻域用局部展开或高精度算术评估
```

**算法二：Fourier 系数后处理**

```text
Input:
    Fourier 系数 \hat f_n, |n|≤N; 已知 r,k; 重构阶 M

Step A:
    右端点中和:
        从 \hat f_n 中减去 x=1 处前 J 个 Taylor 边界项的 Fourier 贡献
    或者
        用逆 FFT 得到过采样物理点值

Step B:
    构造基函数 ψ_m(x)=x^r(-log x)^k p_m(x)

Step C:
    若直接系数法:
        形成连接矩阵 B_{n,m}=∫_0^1 ψ_m(x)e^{-2πinx}dx
        求解正则化问题
            min_a ||Ba-\hat f||_2^2 + τ^2 ||Da||_2^2
    若物理点值法:
        转入算法一

Step D:
    输出
        f_M(x)=Σ_{m=0}^M a_m ψ_m(x)
```

### 对其他策略的判断

若坚持使用“纯 Gegenbauer 重构”，那么最保守、也最符合现有理论的判断是：**它可以改善点值质量，但对 $x^r(-\ln x)^k$ 通常不能单独恢复严格意义的全区间指数精度**。理由并不神秘：Chen–Shu 理论中的决定性步骤是“变换后在 $y$-域内解析”，而对数项恰恰破坏了这一点。与之相比，局部减奇、加权正交、对数映射、nonstandard singular basis 等方法，虽然形式不同，但本质上都在做同一件事：把“非解析部分”显式剥离出去。fileciteturn0file0 fileciteturn0file1 citeturn18search16turn20view1turn20view2

## 数值实验设计建议

数值实验不应只测试 $f(x)=x^r(-\ln x)^k$ 本身，而应测试整个“已知奇异因子 $\times$ 未知解析残差”的函数类。建议最少包含三组基准函数：
$$
f_1(x)=x^r(-\ln x)^k,
\qquad
f_2(x)=x^r(-\ln x)^k e^x,
\qquad
f_3(x)=x^r(-\ln x)^k\frac{1+x}{1+x/2}.
$$
若想评估 Fourier 情形中特有的右端点中和问题，建议再加入
$$
f_4(x)=x^r(-\ln x)^k + (1-x)^k,
$$
因为它能明确分离“左端对数奇异”与“右端解析零点”的竞争。参数上建议 $r$ 取 $\{0.2,\,0.5,\,1,\,1.5,\,2,\,2.5\}$，$k$ 取 $\{1,2,3\}$，从而同时覆盖非整数 $r$、整数 $r$、弱奇异与高阶导数奇异四种情况。citeturn15view0turn16view1turn12view0

误差度量不应只看一个 $L^\infty$ 数。建议至少同时报告：全区间 $L^\infty$ 误差、去掉极小端点邻域后的 interior $L^\infty$ 误差、$L^2$ 误差、谱系数误差，以及若干对数聚点
$$
x=10^{-1},10^{-2},10^{-4},10^{-6}
$$
处的点值相对误差。这样可以区分“全局 max norm 被端点层控制”与“主体区域已指数恢复”这两类截然不同的现象。Zhang–Boyd 的研究特别强调，对端点奇异函数，误差常高度局部化在边界层中，因此只看全局 $L^\infty$ 容易掩盖方法真正的优劣。citeturn4search0turn3view4

参数扫描方面，若测试纯 Gegenbauer 或 Chen–Shu 型方法，建议取
$$
\lambda=\alpha N,\qquad m=\gamma\lambda,
$$
其中 $\alpha\in\{0.05,0.1,0.2,0.3\}$，$\gamma\in\{0.3,0.5,0.8,1.0,1.5\}$。若测试“奇异因子剥离 + 解析残差重构”，则主参数反而应转向残差多项式次数 $M$、截断阈值 $\varepsilon$、局部端点舍弃半径 $x_c$ 与正则化参数 $\tau$。如果采用 Fourier 连接矩阵求解，可把 $M/N$ 从 $0.2$ 扫到 $0.8$，并配合 $\tau\in[10^{-16},10^{-4}]$ 的对数网格搜索。Chen–Shu 2015 的算例值得借鉴：他们的“最优参数”通常给出更小的无噪误差，但对噪声更敏感，因此实验中务必做不同噪声水平下的稳定性图。fileciteturn0file1

建议至少绘制四类图。第一类是误差随 $N$ 变化的半对数图；若出现接近直线的下降，就意味着指数恢复已经实现。第二类是谱系数衰减图，直接检查标准 Chebyshev/Fourier 系数与“加权残差系数”的衰减差异。第三类是点态误差分布图，用来观察边界层是否被成功压缩。第四类是参数敏感性热图，例如 $(\alpha,\gamma)$ 或 $(M,\tau)$ 平面上的 $L^\infty$ 误差分布。citeturn19view0turn15view2turn3view4

最后，建议把实验结果至少按下表组织，以便直观看出“纯方法”与“剥离型方法”的差别。

| 方法 | 输入类型 | 先验信息 | 主要参数 | 预期现象 |
|---|---|---|---|---|
| 标准 Fourier / Chebyshev 偏和 | 原始系数/配点值 | 无 | $N$ | 仅代数–对数收敛，端点附近误差大 |
| 纯 Gegenbauer 重构 | 原始系数/配点值 | 奇异点位置 | $\alpha,\gamma$ | 有改善，但通常难达全区间指数 |
| Chen–Shu 代数映射 | 原始系数/配点值 | 代数指数 | $q,\alpha,\gamma$ | 对纯代数奇异有效；对数奇异仍残留 |
| 奇异因子剥离 + 残差 Chebyshev | 配点值 | 已知 $r,k$ | $M,\varepsilon$ | 最有希望出现稳定指数收敛 |
| 奇异因子剥离 + Fourier 连接矩阵 | Fourier 系数 | 已知 $r,k$ | $M,\tau,J$ | 无噪时近指数，含噪时受正则化控制 |
| 对数坐标 + Laguerre / rational-Chebyshev | 配点值/采样值 | 已知 $r,k$ | 截断长度、映射参数 | 近端点非常强，适合进一步理论推广 |

## 开放问题与研究方向

最重要的理论空白，是证明如下命题：**若 $f(x)=x^r(-\ln x)^k g(x)$，其中 $g$ 在 $[0,1]$ 上解析，而 $r,k$ 已知，那么由有限 Fourier 系数或标准谱配点值出发，经“奇异因子剥离 + 残差谱重构”后，是否能在全区间 $L^\infty$ 范数下得到明确的指数误差界？** 现有理论分别覆盖了“代数奇异的恢复”和“对数奇异的系数衰减”，但这二者之间还缺少一条完整、可发表的连接定理。citeturn19view0turn12view0turn15view0

第二个方向是**参数自适应**。在实践中，$r,k$ 可能并非完全已知，或者仅知道候选模型。这时可以先利用高频 Chebyshev/Fourier 系数的渐近斜率与对数修正项，反推出 $r$ 与 $k$ 的估计，再做奇异因子剥离；这实际上把“后处理”与“奇异性诊断”结合到了一起。Boyd 关于由谱系数反推奇异类型的讨论、以及现代关于弱端点奇异系数渐近的工作，都说明这一路线是可行的。citeturn11search2turn4search0turn13search1

第三个方向是**稳定性改进**。对于 Fourier 连接矩阵方法，当前最需要的是可证明的条件数估计、适合的基函数归一化、以及自动选取正则化参数 $\tau$ 的规则。Chen–Shu 2015 已经表明，在 endpoint singularity 的 Fourier 恢复中，参数越激进并不总是越好，特别是在有噪场景下；这提示未来理论应把噪声与离散病态一起纳入误差分析，而不只是做无噪极限。fileciteturn0file1

第四个方向是**多维推广与非均匀网格**。若二维或三维问题在某个边、角、界面附近具有类似 $d(\mathbf x,\Gamma)^r(\ln d)^k$ 的局部奇异结构，那么“局部距离函数奇异因子剥离 + 解析残差谱重构”有望自然推广；而 rational / lightning approximation 最近在 singularity resolution 上取得的进展，也提示对含对数奇异的多维问题，纯多项式或许并非最优表示。对非均匀网格，generalized Jacobi / weighted basis 的思想也比纯 Gegenbauer 更容易与局部网格适配。citeturn17search5turn17search19turn20view1turn20view2

归结起来，若把研究目标压缩成一句话，那么最值得推进的主线是：

**把 Chen–Shu“映射后解析”的思想，改写成“剥离已知对数奇异因子后解析”的思想；然后对解析残差建立从 Fourier 系数或标准配点值出发的严格指数重构理论。**

这条路线既忠实继承了 Gottlieb–Shu 的信息重读哲学，也真正针对了 $x^r(-\ln x)^k$ 这一类函数的结构本质。对当前问题而言，它比试图强行把纯对数奇异塞进原始 Gegenbauer 机械框架，更自然，也更有可能得到完整而漂亮的理论结果。citeturn19view1turn19view3turn12view0turn16view1

## 参考文献选编

Gottlieb, D., Shu, C.-W. “On the Gibbs phenomenon V: recovering exponential accuracy from collocation point values of a piecewise analytic function.” *Numerische Mathematik* 71 (1995). citeturn3view1

Gottlieb, D., Shu, C.-W. “On the Gibbs phenomenon and its resolution.” *SIAM Review* 39 (1997). citeturn3view0turn19view1turn19view2

Chen, Z., Shu, C.-W. “Recovering exponential accuracy from collocation point values of smooth functions with end-point singularities.” *Journal of Computational and Applied Mathematics* 265 (2014). fileciteturn0file0

Chen, Z., Shu, C.-W. “Recovering Exponential Accuracy in Fourier Spectral Methods Involving Piecewise Smooth Functions with Unbounded Derivative Singularities.” *Journal of Scientific Computing* 65 (2015). fileciteturn0file1

Adcock, B., Richardson, M. “New exponential variable transform methods for functions with endpoint singularities.” *SIAM Journal on Numerical Analysis* 52 (2014). citeturn17search0turn17search2

Wang, H. “Convergence rate and acceleration of Clenshaw–Curtis quadrature for functions with endpoint singularities.” arXiv and related journal version (2014). citeturn8view2turn10view0

Sidi, A. “Asymptotic expansions of Legendre series coefficients for functions with endpoint singularities.” *Journal of Applied Analysis* / SAGE entry (2009). citeturn7search0turn8view1

Sidi, A. “A simple approach to asymptotic expansions for Fourier integrals of singular functions.” *Applied Mathematics and Computation* 216 (2010). citeturn12view0

Xiang, S., Liu, G. “Optimal decay rates on the asymptotics of orthogonal polynomial expansions for functions of limited regularities.” *Numerische Mathematik* 145 (2020). citeturn13search1turn15view0turn15view2

Zhang, X., Boyd, J. P. “Asymptotic coefficients and errors for Chebyshev polynomial approximations with weak endpoint singularities: Effects of different bases.” *Science China Mathematics* 66 (2023). citeturn4search0turn3view4

Guo, B.-Y., Shen, J., Wang, L.-L. “Optimal Spectral-Galerkin Methods Using Generalized Jacobi Polynomials.” *Journal of Scientific Computing* 27 (2006). citeturn20view1

Chen, S., Shen, J., Wang, L.-L. “Generalized Jacobi Functions and Their Applications to Fractional Differential Equations.” *Mathematics of Computation* 85 (2016) / arXiv 2014. citeturn18search4

Mao, Z., Chen, S., Shen, J. “Efficient and accurate spectral method using generalized Jacobi functions for solving Riesz fractional differential equations.” *Applied Numerical Mathematics* 106 (2016). citeturn20view2