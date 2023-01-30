# 二次规划 {#sec:7.2}

具有二次目标函数和线性约束的约束优化问题被称为二次规划 (quadratic
programming, QP) 问题. 这是一类非常常见的约束优化问题,
且它们也常常作为一般的约束优化问题的子问题出现. 二次规划问题的一般形式为
$$\label{eq:quadratic-programming-1}
\begin{array}{cl}
\text{minimize} & q({\bm{x}}) := \frac{1}{2} {\bm{x}}^T G {\bm{x}} + {\bm{d}}^T {\bm{x}}, \\
\text{subject to} & {\bm{a}}_i^T {\bm{x}} = b_i, ~~ i \in \mathcal{E} = \{1, \ldots, m_1\}, \\
& {\bm{a}}_i^T {\bm{x}} \leqslant b_i, ~~ i \in \mathcal{I} = \{m_1 + 1, \ldots, m\},
\end{array}$$ 其中$G$是$n$阶对称方阵, ${\bm{d}}, {\bm{a}}_i,$
是$n$维列向量, $b_i$是常数, $i \in \mathcal{E} \cup \mathcal{I},$
约束条件有等式约束以及不等式约束.
以上形式的二次规划问题可以将约束条件以矩阵的形式表达为更加紧凑的形式：
$$\label{eq:quadratic-programming-2}
\begin{array}{cl}
\text{minimize} & q({\bm{x}}) := \frac{1}{2} {\bm{x}}^T G {\bm{x}} + {\bm{d}}^T {\bm{x}}, \\
\text{subject to} & A^T_{\mathcal{E}} {\bm{x}} = {\bm{b}}_{\mathcal{E}}, ~~ A^T_{\mathcal{I}} x \preccurlyeq {\bm{b}}_{\mathcal{I}},
\end{array}$$
其中$A_{\mathcal{E}} = ({\bm{a}}_1, \ldots, {\bm{a}}_{m_1}), {\bm{b}}_{\mathcal{E}} = (b_1, \ldots, b_{m_1})^T;$
$A_{\mathcal{I}} = ({\bm{a}}_{m_1+1}, \ldots, {\bm{a}}_m), {\bm{b}}_{\mathcal{I}} = (b_{m_1+1}1, \ldots, b_m)^T.$
$A^T_{\mathcal{I}} {\bm{x}} \preccurlyeq {\bm{b}}_{\mathcal{I}}$表示向量$A^T_{\mathcal{I}} {\bm{b}}_{\mathcal{I}}$的每个元素都小于等于向量${\bm{b}}_{\mathcal{I}}$相应位置元素.
问题[\[eq:quadratic-programming-2\]](#eq:quadratic-programming-2){reference-type="eqref"
reference="eq:quadratic-programming-2"}的KKT条件为
$$\label{eq:quadratic-programming-kkt}
\begin{aligned}
& \nabla q({\bm{x}}^*) + A {\bm{\lambda}}^* = {\bm{0}}, \\
& A^T_{\mathcal{E}} {\bm{x}}^* - {\bm{b}}_{\mathcal{E}} = {\bm{0}}, \\
& {\bm{\lambda}}^*_{\mathcal{I}} \succcurlyeq {\bm{0}}, ~ A_{\mathcal{I}}^T {\bm{x}}^* - {\bm{b}}_{\mathcal{I}} \preccurlyeq {\bm{0}}, ~ ({\bm{\lambda}}^*_{\mathcal{I}})^T (A_{\mathcal{I}}^T {\bm{x}}^* - {\bm{b}}_{\mathcal{I}}) = {\bm{0}},
\end{aligned}$$ 其中$A = ({\bm{a}}_1, \ldots, {\bm{a}}_m);$
${\bm{\lambda}}^* = (\lambda_1^*, \ldots, \lambda_m^*)^T$为拉格朗日乘子,
${\bm{\lambda}}^*_{\mathcal{I}} = (\lambda_{m_1+1}^*, \ldots, \lambda_m^*)^T;$
${\bm{x}}^*$为KKT点. 值得注意的是, 由于二次规划问题线性约束规范总成立,
二次规划问题的任一解都满足KKT条件[\[eq:quadratic-programming-kkt\]](#eq:quadratic-programming-kkt){reference-type="eqref"
reference="eq:quadratic-programming-kkt"}.

与线性规划问题类似, 一个二次规划问题总是能在有限步内求解,
或者证明其是不可行的或者无(下）界的. 当二次规划问题有解时,
若目标函数$q({\bm{x}}) = \frac{1}{2} {\bm{x}}^T G {\bm{x}} + {\bm{d}}^T {\bm{x}}$的Hessian阵$G$是半正定的,
那么$q({\bm{x}})$是凸函数, 则该二次规划问题是凸规划问题,
其KKT点即为全局极小值点; 如果$G$是正定的,
则解是唯一的(唯一的全局极小值点). 当$G$不定时,
目标函数$q({\bm{x}})$非凸, 可能有多个鞍点以及局部极小值点,
这些局部极小值点可能不是全局最优解. 此时,
求解全局最优解是一个NP-难的问题[@Murty_1987].

## 等式约束的二次规划问题 {#subsec:7.2.1}

我们首先考虑最简单的等式约束二次规划问题, 问题形式如下
$$\label{eq:quadratic-programming-eq-constrained}
\begin{array}{cl}
\text{minimize} & q({\bm{x}}) := \frac{1}{2} {\bm{x}}^T G {\bm{x}} + {\bm{d}}^T {\bm{x}}, \\
\text{subject to} & A^T {\bm{x}} = {\bm{b}},
\end{array}$$
其中$A = ({\bm{a}}_1, \ldots, {\bm{a}}_m)^T$是一个$n \times m$的矩阵,
${\bm{a}}_i \in \mathcal{R}^n,$ ${\bm{b}} \in \mathcal{R}^m.$
我们假设$A$列满秩, 即$A$的秩为$m \leqslant n.$ 如若不然,
只需要去掉线性相关的约束即可.

求解等式约束的二次规划问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained"}的一个最直接的想法是通过等约束条件消去部分变量,
将问题化为无约束极小化问题进行求解. 我们将矩阵$A, ~ G,$ 向量${\bm{d}}$,
以及变量${\bm{x}}$进行分块 $$\label{eq:quadratic-programming-eq-blocks}
A = \begin{bmatrix} A_1 \\ A_2 \end{bmatrix}, ~ G = \begin{bmatrix} G_{11} & G_{12} \\ G_{21} & G_{22} \end{bmatrix}, ~ {\bm{d}} = \begin{bmatrix} {\bm{d}}_1 \\ {\bm{d}}_2 \end{bmatrix}, ~ {\bm{x}} = \begin{bmatrix} {\bm{x}}_1 \\ {\bm{x}}_2 \end{bmatrix},$$
其中$A_1$是可逆的$m$阶方阵, ${\bm{x}}_1 \in \mathcal{R}^{m},$
${\bm{x}}_2 \in \mathcal{R}^{n-m},$
矩阵$G$的分块以及向量${\bm{d}}$的分块与${\bm{x}}$的分块相容.
那么等式约束的二次规划问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained"}中的约束方程组$A^T {\bm{x}} = {\bm{b}}$就可写作$A_1^T {\bm{x}}_1 + A_2^T {\bm{x}}_2 = {\bm{b}}.$
由于$A_1$可逆,
可以利用高斯消元法求解关于${\bm{x}}_1$的线性方程组$A_1^T {\bm{x}}_1 = {\bm{b}} - A_2^T {\bm{x}}_2$得唯一解
$$\label{eq:quadratic-programming-eq-constrained-subs}
{\bm{x}}_1 = A_1^{-T} ({\bm{b}} - A_2^T {\bm{x}}_2).$$
将上式代入目标函数$q({\bm{x}}),$ 可得关于${\bm{x}}_2$的二次函数
$$\begin{aligned}
\psi ({\bm{x}}_2) := & \frac{1}{2} {\bm{x}}_2^T \left( G_{22} - G_{21} A_1^{-T} A_2^T - A_2 A_1^{-1} G_{12} + A_2 A_1^{-1} G_{12} A_1^{-T} A_2^T \right) {\bm{x}}_2 \\
& + {\bm{x}}_2^T \left( G_{21} - A_2 A_1^{-1} G_{11} \right) A_1^{-T} {\bm{b}} + \frac{1}{2} {\bm{b}}^T A_1^{-1} G_{11} A_1^{-T} {\bm{b}} \\
& + {\bm{x}}_2^T \left( {\bm{d}}_2 - A_2 A_1^{-1} {\bm{d}}_1 \right) + {\bm{d}}_1^T A_1^{-T} {\bm{b}}.
\end{aligned}$$ 于是,
等式约束的二次规划问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained"}转化为了无约束极小化问题
$$\text{minimize} ~~ \psi({\bm{x}}_2)$$ 如果$\nabla^2 \psi$正定,
解(一阶条件)线性方程组$\nabla \psi ({\bm{x}}_2) = 0$可得唯一解${\bm{x}}_2^*,$
回代入[\[eq:quadratic-programming-eq-constrained-subs\]](#eq:quadratic-programming-eq-constrained-subs){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained-subs"}可得${\bm{x}}_1^* = A_1^{-1} ( {\bm{b}} - A_2 {\bm{x}}_2^* ).$
拉格朗日乘子${\bm{\lambda}}^*$由KKT条件[\[eq:quadratic-programming-kkt\]](#eq:quadratic-programming-kkt){reference-type="eqref"
reference="eq:quadratic-programming-kkt"}给出, 即由
$$\label{eq:quadratic-programming-eq-lambda}
\begin{bmatrix} A_1 \\ A_2 \end{bmatrix} {\bm{\lambda}}^* = A {\bm{\lambda}}^* = -\nabla q({\bm{x}}^*) = - G {\bm{x}}^* - {\bm{d}} = - \begin{bmatrix} G_{11} & G_{12} \\ G_{21} & G_{22} \end{bmatrix} \begin{bmatrix} {\bm{x}}_1^* \\ {\bm{x}}_2^* \end{bmatrix} - \begin{bmatrix} {\bm{d}}_1 \\ {\bm{d}}_2 \end{bmatrix}$$
确定. 解得拉格朗日乘子的显式表达
$${\bm{\lambda}}^* = -A_1^{-1} \left( G_{11} {\bm{x}}_1^* + G_{12} {\bm{x}}_2^* + {\bm{d}}_1 \right).$$

::: {#eg:7.2.1 .eg}
**例 1**. *我们用直接消元法求解下列等式约束的二次规划问题
$$\begin{array}{cl}
\text{minimize} & q({\bm{x}}) = x_1^2 + x_2^2 + 2x_3^2 + 5x_1x_2 - 3x_2 - 7x_3, \\
\text{subject to} & x_1 + x_2 + x_3 = 1, \\
& x_1 - 2x_2 - 3x_3 = -2,
\end{array}$$ 那么相应的矩阵及其划分为
$$A = \left[\begin{array}{@{}cc@{}} 1 & 1 \\ 1 & -2 \\ \hdashline 1 & -3 \end{array}\right], ~~ G = \left[\begin{array}{@{}cc:c@{}} 2 & 5 & 0 \\ 5 & 2 & 0 \\ \hdashline 0 & 0 & 4 \end{array}\right], ~~ {\bm{d}} = \left[\begin{array}{@{}c@{}} 0 \\ -3 \\ \hdashline -7 \end{array}\right].$$*

*我们将约束条件写为 $$\begin{aligned}
x_1 + x_2 & = 1 - x_3, \\
x_1 - 2x_2 & = -2 + 3x_3.
\end{aligned}$$ 利用高斯消元法,
解得$x_1 = \frac{1}{3} x_3, ~ x_2 = 1 - \frac{4}{3} x_3.$
代入原目标函数得 $$\psi(x_3) = \frac{5 x_{3}^{2}}{3} - 4 x_{3} - 2.$$
$\psi$的Hessian阵$\nabla^2 \psi = \psi'' = \frac{10}{3} > 0,$
于是有唯一解.
由$0 = \nabla \psi (x_3) = \frac{10 x_{3}}{3} - 4$解得$x_3^* = \frac{6}{5}.$
回代得$x_1^* = \frac{2}{5},$ $x_2^* = -\frac{3}{5},$
即${\bm{x}}^* = \left( \frac{2}{5}, -\frac{3}{5}, \frac{6}{5} \right)^T.$
最后,
解方程组[\[eq:quadratic-programming-eq-lambda\]](#eq:quadratic-programming-eq-lambda){reference-type="eqref"
reference="eq:quadratic-programming-eq-lambda"}, 即
$$\begin{bmatrix} 1 & 1 \\ 1 & -2 \\ 1 & -3 \end{bmatrix} ~ \begin{bmatrix} \lambda_1^* \\ \lambda_2^* \end{bmatrix} = - \begin{bmatrix} 2 & 5 & 0 \\ 5 & 2 & 0 \\ 0 & 0 & 4 \end{bmatrix} \begin{bmatrix} \frac{2}{5} \\ -\frac{3}{5} \\ \frac{6}{5} \end{bmatrix} - \begin{bmatrix} 0 \\ -3 \\ -7 \end{bmatrix},$$
算得拉格朗日乘子$\lambda_1^* = \frac{11}{5}, \lambda_2^* = 0.$*
:::

直接消元法是一种朴素初等的方法, 非常直观,
但不是求解等式约束二次规划问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained"}的最优方法.
我们还可以采用广义消元法(generalized elimination)来求解.
这个方法的名字来源于矩阵的广义逆,
利用矩阵的广义逆来表达问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained"}等式约束(线性方程组)$A^T {\bm{x}} = {\bm{b}}$的通解,
回代入目标函数$q({\bm{x}}),$ 得到一个无约束的极小化问题.
广义消元法的本质是对变量做了线性变换,
是直接消元法[\[eq:quadratic-programming-eq-constrained-subs\]](#eq:quadratic-programming-eq-constrained-subs){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained-subs"}的推广.

我们知道, 对于非齐次线性方程组$A^T {\bm{x}} = {\bm{b}},$ 它如果有解,
那么它的通解可以表示为 $$\label{eq:general-elim-1}
{\bm{x}} = Y {\bm{b}} + {\bm{s}}$$
其中$Y \in \mathcal{R}^{n\times m}$是系数矩阵$A^T$的某个广义逆,
${\bm{s}}$是齐次线性方程组$A^T {\bm{x}} = 0$的解,
或者说属于$A$的列零空间(null column space).
${\bm{s}}$也被称作可行点${\bm{x}}$处的可行增量. 我们回忆一下,
$m\times n$的矩阵$A^T$的广义逆, 指的是一个$n\times m$的矩阵$Y,$
满足$A^T Y A^T = A^T.$ 同时,
非齐次线性方程组$A^T {\bm{x}} = {\bm{b}}$有解当且仅当$A^T Y {\bm{b}} = {\bm{b}}.$
假设我们取定了$n \times m$的矩阵$Y$以及$n \times (n-m)$的矩阵$Z,$ 满足
$$\label{eq:eq:general-elim-req}
\begin{bmatrix} Y & Z\end{bmatrix} \text{ 非奇异, 且~} A^T Y = I_m, ~ A^T Z = 0.$$
也就是说$Y$是$A^T$的一个右广义逆(比广义逆的定义$A^T Y A^T = A^T$更进一步要求$A^T Y = I_m$).
而矩阵$Z$的列${\bm{z}}_1, \ldots, {\bm{z}}_{n-m}$构成了$A$的列零空间的一组基,
也被称作既约(reduced)坐标方向, 矩阵$Z$也被称作零空间矩阵(null-space
matrix).
矩阵$Z^T G Z$被称作二次规划问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained"}的既约Hessian阵.
那么,
可行增量${\bm{s}}$可以表示为${\bm{z}}_1, \ldots, {\bm{z}}_{n-m}$的线性组合
$$\label{eq:general-elim-2}
{\bm{s}} = Z {\bm{y}} = \sum\limits_{i=1}^{n-m} y_i {\bm{z}}_i,$$
其中$y_1, \ldots, y_{n-m}$是每个既约坐标方向的分量,
我们称${\bm{y}} = (y_1, \ldots, y_{n-m})^T$为既约变量.
将式[\[eq:general-elim-2\]](#eq:general-elim-2){reference-type="eqref"
reference="eq:general-elim-2"}代入式[\[eq:general-elim-1\]](#eq:general-elim-1){reference-type="eqref"
reference="eq:general-elim-1"}, 于是每一个可行点${\bm{x}}$都可以表示为
$$\label{eq:general-elim-3}
{\bm{x}} = Y {\bm{b}} + Z {\bm{y}}.$$ 这样,
我们就可以用既约变量${\bm{y}}$代换原变量${\bm{x}},$ 将等式约束消去,
将原问题转化为了无约束的既约二次函数的极小化问题
$$\label{eq:general-elim-4}
\text{minimize} ~~ \psi({\bm{y}}) := \frac{1}{2} {\bm{y}}^T Z^T G Z {\bm{y}} + \left( {\bm{d}} + G Y {\bm{b}} \right)^T Z {\bm{y}} + \frac{1}{2} \left( 2{\bm{d}} + G Y {\bm{b}} \right)^T Y {\bm{b}}.$$

下面我们来解无约束的极小化问题[\[eq:general-elim-4\]](#eq:general-elim-4){reference-type="eqref"
reference="eq:general-elim-4"}. 考虑一阶条件 $$\label{eq:general-elim-5}
0 = \nabla \psi ({\bm{y}}) = Z^T G Z {\bm{y}} + Z^T \left( {\bm{d}} + G Y {\bm{b}} \right),$$
如果既约Hessian阵$Z^T G Z$是正定的, 那么有唯一解.
解关于${\bm{y}}$的线性方程组 $$\label{eq:general-elim-y}
Z^T G Z {\bm{y}} = - Z^T \left( {\bm{d}} + G Y {\bm{b}} \right)$$
得唯一解${\bm{y}}^*,$
再代入变量替换式[\[eq:quadratic-programming-eq-constrained-subs\]](#eq:quadratic-programming-eq-constrained-subs){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained-subs"}得到原问题的唯一解${\bm{x}}^* = Y {\bm{b}} + Z {\bm{y}}^*.$
由KKT条件[\[eq:quadratic-programming-kkt\]](#eq:quadratic-programming-kkt){reference-type="eqref"
reference="eq:quadratic-programming-kkt"},
拉格朗日乘子${\bm{\lambda}}^*$满足 $$\label{eq:general-elim-lagrange}
A {\bm{\lambda}}^* = -\left( G {\bm{x}}^* + {\bm{d}} \right) =: - {\bm{g}}^*,$$
两边同时左乘$Y^T,$ 有 $$\label{eq:general-elim-lambda}
{\bm{\lambda}}^* = \left( A^T Y \right)^T {\bm{\lambda}}^* = -Y^T {\bm{g}}^*.$$

到目前为止, 我们还没有介绍如何选取, 或者说构造矩阵$Y$和$Z.$ 矩阵$Y,$
$Z$的构造方法有很多种, 通常首选的是正交分解法(orthogonal factorization
method). 我们对矩阵$A$进行QR分解, 即
$$\label{eq:quadratic-programming-qr-decomp-1}
A = Q \begin{bmatrix} R \\ 0 \end{bmatrix} = \begin{bmatrix} Q_1 & Q_2 \end{bmatrix} \begin{bmatrix} R \\ 0 \end{bmatrix} = Q_1 R,$$
其中$Q$是$n \times n$的正交阵, $R$是$m \times m$的可逆上三角阵.
$Q_1, Q_2$分别是$n \times m$和$n \times (n - m)$的矩阵. 然后选取
$$\label{eq:quadratic-programming-qr-decomp-2}
Y = Q_1 R^{-T}, ~~ Z = Q_2$$
即可满足[\[eq:eq:general-elim-req\]](#eq:eq:general-elim-req){reference-type="eqref"
reference="eq:eq:general-elim-req"}中关于矩阵$Y,$ $Z$的要求.
通过解线性方程组$R^T {\bm{v}} = {\bm{b}}$得到一个特解${\bm{v}},$
进而得到[\[eq:general-elim-3\]](#eq:general-elim-3){reference-type="eqref"
reference="eq:general-elim-3"}中的特解$Y {\bm{b}} = Y R^T {\bm{v}} = Q_1 {\bm{v}}.$
接下来我们只要按照广义消元法的程序,
解关于既约变量${\bm{y}}$的方程组[\[eq:general-elim-y\]](#eq:general-elim-y){reference-type="eqref"
reference="eq:general-elim-y"}得到唯一解${\bm{y}}^*,$
回代得到原问题最优解${\bm{x}}^* = Y {\bm{b}} + Z {\bm{y}}^*,$
以及用式[\[eq:general-elim-lambda\]](#eq:general-elim-lambda){reference-type="eqref"
reference="eq:general-elim-lambda"}计算得到相应的拉格朗日乘子.

矩阵$Y,$ $Z$的通用计算格式如下:
任选一个$n \times (n - m)$的矩阵$V$使得$\begin{bmatrix} A & V \end{bmatrix}$是非奇异的,
且逆可以分块表示为 $$\label{eq:quadratic-programming-qr-decomp-general}
\begin{bmatrix} A & V \end{bmatrix}^{-1} = \begin{bmatrix} Y^T \\ Z^T \end{bmatrix},$$
其中$Y$和$Z$分别是$n \times m$的矩阵和$n \times (n - m)$的矩阵.
依照矩阵逆的定义, 有
$$I_n = \begin{bmatrix} Y^T \\ Z^T \end{bmatrix} \begin{bmatrix} A & V \end{bmatrix} = \begin{bmatrix} Y^T A & Y^T V \\ Z^T A & Z^T V \end{bmatrix},$$
即有$Y^T A = I_m, Z^T A = 0,$
满足要求[\[eq:eq:general-elim-req\]](#eq:eq:general-elim-req){reference-type="eqref"
reference="eq:eq:general-elim-req"}, 适用于广义消元法.
这也可以解释之前介绍的直接消元法,
即如果取$V = \begin{bmatrix} 0 \\ I_{n - m} \end{bmatrix},$
同时假设矩阵$A$有形如[\[eq:quadratic-programming-eq-blocks\]](#eq:quadratic-programming-eq-blocks){reference-type="eqref"
reference="eq:quadratic-programming-eq-blocks"}中的分块,
那么矩阵$\begin{bmatrix} A & V \end{bmatrix}$的逆可以显式地表示为
$$\begin{bmatrix} Y^T \\ Z^T \end{bmatrix} = \begin{bmatrix} A_1 & 0 \\ A_2 & I_{n - m} \end{bmatrix}^{-1} = \begin{bmatrix} A_1^{-1} & 0 \\ -A_2A_1^{-1} & I_{n - m} \end{bmatrix}.$$
这样我们就得到了之前的直接消元法. 此外, 如果取$V = Q_2,$
其中$Q_2$由矩阵$A$的QR分解式[\[eq:quadratic-programming-qr-decomp-1\]](#eq:quadratic-programming-qr-decomp-1){reference-type="eqref"
reference="eq:quadratic-programming-qr-decomp-1"}给出, 那么由
$$\begin{bmatrix} Y^T \\ Z^T \end{bmatrix} = \begin{bmatrix} A & V \end{bmatrix}^{-1} = \begin{bmatrix} Q_1R & Q_2 \end{bmatrix}^{-1} = \begin{bmatrix} R^{-1} Q_1^T \\ Q_2^T \end{bmatrix},$$
我们可得正交分解法. 尽管关于矩阵$Y$和$Z$有各种各样的选取方式,
然而正交分解法通常是首选的. 原因之一是计算矩阵$Z$时涉及正交矩阵的操作,
而正交矩阵具有良好的数值稳定性；
其次是选取$Z = Q_2$可以给出条件数$\kappa(Z^T G Z)$的最优上界, 即
$$\kappa(Z^T G Z) \leqslant \kappa(G).$$

::: {#eg:7.2.2 .eg}
**例 2**.
*我们用广义消元法(正交分解法)来求解例[\[eg:7.2.1\]](#eg:7.2.1){reference-type="eqref"
reference="eg:7.2.1"}中的二次规划问题.
我们首先得到约束条件系数矩阵的转置矩阵$A$的QR分解
$$A = \begin{bmatrix} 1 & 1 \\ 1 & -2 \\ 1 & -3 \end{bmatrix} = \left[\begin{array}{@{}cc:c@{}} \frac{\sqrt{3}}{3} & \frac{7 \sqrt{78}}{78} & \frac{\sqrt{26}}{26} \\ \frac{\sqrt{3}}{3} & - \frac{\sqrt{78}}{39} & - \frac{2 \sqrt{26}}{13} \\ \frac{\sqrt{3}}{3} & - \frac{5 \sqrt{78}}{78} & \frac{3 \sqrt{26}}{26} \end{array}\right] ~ \left[\begin{array}{@{}cc@{}} \sqrt{3} & - \frac{4 \sqrt{3}}{3} \\ 0 & \frac{\sqrt{78}}{3} \\ \hdashline 0 & 0 \end{array}\right]$$
依照式[\[eq:quadratic-programming-qr-decomp-2\]](#eq:quadratic-programming-qr-decomp-2){reference-type="eqref"
reference="eq:quadratic-programming-qr-decomp-2"}, 令
$$Y = \left[\begin{matrix} \frac{\sqrt{3}}{3} & \frac{7 \sqrt{78}}{78} \\ \frac{\sqrt{3}}{3} & - \frac{\sqrt{78}}{39} \\ \frac{\sqrt{3}}{3} & - \frac{5 \sqrt{78}}{78} \end{matrix}\right] ~ \left[\begin{matrix} \sqrt{3} & - \frac{4 \sqrt{3}}{3} \\ 0 & \frac{\sqrt{78}}{3} \end{matrix}\right]^{-T} = \frac{1}{26} \left[\begin{matrix} 18 & 7 \\ 6 & -2 \\ 2 & -5 \end{matrix}\right], ~~ Z = \frac{1}{\sqrt{26}} \left[\begin{matrix} 1 \\ -4 \\ 3 \end{matrix}\right].$$
此时, $Z^T G Z = \frac{15}{13} > 0,$ 有唯一解.
代入关于既约变量$y$的线性方程组[\[eq:general-elim-y\]](#eq:general-elim-y){reference-type="eqref"
reference="eq:general-elim-y"}有
$$\frac{15}{13} y = - Z^T \left( d + G Y b \right) = \frac{48 \sqrt{26}}{169},$$
解得$y^* = \frac{16 \sqrt{26}}{65}.$ 回代得到原问题的最优解
$${\bm{x}}^* = Y {\bm{b}} + Z y^* = \left[\begin{matrix} \frac{2}{5} \\ - \frac{3}{5} \\ \frac{6}{5} \end{matrix}\right],$$
同时,
利用式[\[eq:general-elim-lambda\]](#eq:general-elim-lambda){reference-type="eqref"
reference="eq:general-elim-lambda"}得到拉格朗日乘子
$${\bm{\lambda}}^* = - Y^T \left( G {\bm{x}}^* + {\bm{d}} \right) = \left[\begin{matrix} \frac{11}{5} \\ 0 \end{matrix}\right].$$
这与我们在例[\[eg:7.2.1\]](#eg:7.2.1){reference-type="eqref"
reference="eg:7.2.1"}中用直接消元法得到的结果是一样的.*
:::

我们还可以直接从等式约束的二次规划问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained"}的KKT条件入手,
来寻求该问题的解.
等式约束的二次规划问题的KKT条件可以从式[\[eq:quadratic-programming-kkt\]](#eq:quadratic-programming-kkt){reference-type="eqref"
reference="eq:quadratic-programming-kkt"}简化为
$$\label{eq:quadratic-programming-kkt-2}
\begin{aligned}
& G {\bm{x}}^* + {\bm{d}} + A {\bm{\lambda}}^* = {\bm{0}}, \\
& A^T {\bm{x}}^* - {\bm{b}} = {\bm{0}}.
\end{aligned}$$ 写为矩阵形式, 即为
$$\label{eq:quadratic-programming-kkt-2-mat}
\begin{bmatrix} G & A \\ A^T & 0 \end{bmatrix} \begin{bmatrix} {\bm{x}}^* \\ {\bm{\lambda}}^* \end{bmatrix} = \begin{bmatrix} -{\bm{d}} \\ {\bm{b}} \end{bmatrix}.$$
上式中的系数矩阵$K := \begin{bmatrix} G & A \\ A^T & 0 \end{bmatrix}$被称作拉格朗日矩阵,
也被称作KKT矩阵, 很容易看出KKT矩阵是对称不定的.
那么等式约束的二次规划问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
reference="eq:quadratic-programming-eq-constrained"}的求解就转化为了线性方程组[\[eq:quadratic-programming-kkt-2-mat\]](#eq:quadratic-programming-kkt-2-mat){reference-type="eqref"
reference="eq:quadratic-programming-kkt-2-mat"}的求解.

如果KKT矩阵可逆, 我们将逆矩阵按相同的方式进行分块, 即
$$\label{eq:quadratic-programming-kkt-inverse}
\begin{bmatrix} G & A \\ A^T & 0 \end{bmatrix}^{-1} = \begin{bmatrix} C & E \\ E^T & F \end{bmatrix}$$
那么KKT系统[\[eq:quadratic-programming-kkt-2-mat\]](#eq:quadratic-programming-kkt-2-mat){reference-type="eqref"
reference="eq:quadratic-programming-kkt-2-mat"}的解就可写为
$$\label{eq:quadratic-programming-kkt-sol-1}
\begin{aligned}
{\bm{x}}^* & = - C {\bm{d}} + E {\bm{b}}, \\
{\bm{\lambda}}^* & = - E^T {\bm{d}} + F {\bm{b}}. \\
\end{aligned}$$ 当$G$可逆时, 我们可以利用分块矩阵的消元法算得
$$\label{eq:quadratic-programming-kkt-inv-1}
\begin{aligned}
C & = G^{-1} - G^{-1} A \left( A^TG^{-1}A \right)^{-1} A^T G^{-1},\\
E & = G^{-1} A \left( A^TG^{-1}A \right)^{-1},\\
F & = - \left( A^TG^{-1}A \right)^{-1}.\\
\end{aligned}$$ 相关的方法被称作Schur补法(Schur-complement method),
这个名字来源于, 执行分块消元操作的时候,
我们会得到矩阵$G$的Schur补$-A^TG^{-1}A.$ 要注意的是,
Schur补法要求$G$可逆, 而这并不是KKT矩阵非奇异的必要条件;
同时基于数值稳定性等方面的考虑,
只有当$G$正定以及$A^TG^{-1}A$的计算开销足够小等一系列条件得到满足,
Schur补法才有应用的价值.

我们可以利用式[\[eq:quadratic-programming-qr-decomp-general\]](#eq:quadratic-programming-qr-decomp-general){reference-type="eqref"
reference="eq:quadratic-programming-qr-decomp-general"}定义的矩阵$Y,$
$Z$, 来求KKT矩阵的逆矩阵, 有
$$\label{eq:quadratic-programming-kkt-inv-2}
\begin{aligned}
C & = Z \left( Z^T G Z \right)^{-1} Z^T,\\
E & = Y - Z \left( Z^T G Z \right)^{-1} Z^T G Y,\\
F & = Y^T G Z \left( Z^T G Z \right)^{-1} Z^T G Y - Y^T G Y.\\
\end{aligned}$$ 这样, 我们可以利用矩阵$Y,$
$Z$和$Z^T G Z$的Cholesky分解$L L^T$等计算方法来分解KKT矩阵,
相关的方法被称作零空间法(null-space method).
这一名称来源于矩阵$Z$的选取,
即$Z$的列构成约束条件系数矩阵的转置矩阵$A$的列零空间.
对应于不同的$Y$和$Z$的选取, 有不同的零空间法, 如上面介绍的广义消元法.

## 求解二次规划问题的积极集法 {#subsec:7.2.2}

对于求解一般的带不等式约束的二次规划问题[\[eq:quadratic-programming-1\]](#eq:quadratic-programming-1){reference-type="eqref"
reference="eq:quadratic-programming-1"},
我们有一系列的实用算法来求解这些问题. 经典的积极集法(active-set
methods)自从20世纪70年代起被广泛应用于求解二次规划问题.
积极集法适用于求解中小规模(成百上千个变量)的凸和非凸的二次规划问题.
梯度投影法(gradient-projection methods)是属于积极集法的一种特殊的算法,
是经典积极集法的推广, 能够非常高效地求解简单约束的二次规划问题,
例如对每个变量的约束都是区间约束的二次规划问题($a_i \leqslant x_i \leqslant b_i,$
我们称这样的问题为BoxQP). 还有一类方法是内点法(interior-point methods),
广泛应用于求解二次规划问题的时间比经典积极集法稍晚,
大概始于20世纪90年代. 内点法适用于求解大规模的凸二次规划问题.
本节主要介绍如何利用积极集法将等式约束问题的求解方法推广以求解带不等式约束的二次规划问题.
为描述简单, 我们假设原问题是凸二次规划问题.

我们回忆一下,
对于二次规划问题[\[eq:quadratic-programming-1\]](#eq:quadratic-programming-1){reference-type="eqref"
reference="eq:quadratic-programming-1"}, 积极集$\mathcal{A}(x)$的定义为
$$\label{eq:qp-active-set}
\mathcal{A} = \mathcal{A}({\bm{x}}) = \left\{ i : ~ {\bm{a}}_i^T {\bm{x}} = b_i, ~ i \in \mathcal{E} \cup \mathcal{I} \right\},$$
即在点${\bm{x}}$处, 等式成立的约束条件的指标(index)组成的集合.
我们将一般的二次规划的KKT条件[\[eq:quadratic-programming-kkt\]](#eq:quadratic-programming-kkt){reference-type="eqref"
reference="eq:quadratic-programming-kkt"}根据积极集$\mathcal{A}({\bm{x}}^*)$改写一下,
即有 $$\label{eq:qp-active-set-kkt}
\begin{aligned}
& G {\bm{x}}^* + d + \sum\limits_{i \in \mathcal{A}({\bm{x}}^*)} \lambda_i^* {\bm{a}}_i = {\bm{0}}, \\
& {\bm{a}}_i {\bm{x}}^* = b_i, ~~ \forall i \in \mathcal{A}({\bm{x}}^*), \\
& {\bm{a}}_i {\bm{x}}^* \leqslant b_i, ~~ \forall i \in \mathcal{I} \setminus \mathcal{A}({\bm{x}}^*), \\
& \lambda_i^* \geqslant 0, ~~ \forall i \in \mathcal{I} \cap \mathcal{A}({\bm{x}}^*).
\end{aligned}$$ 不难发现, ${\bm{x}}^*$也是下面等式问题的KKT点
$$\label{eq:qp-active-set-1}
\begin{array}{cl}
\text{minimize} & \frac{1}{2} {\bm{x}}^T G {\bm{x}} + {\bm{d}}^T {\bm{x}}, \\
\text{subject to} & {\bm{a}}_i^T {\bm{x}} = b_i, ~ i \in \mathcal{A}({\bm{x}}^*).
\end{array}$$ 这说明,
求解含不等式约束的二次规划问题几乎(注意KKT条件[\[eq:qp-active-set-kkt\]](#eq:qp-active-set-kkt){reference-type="eqref"
reference="eq:qp-active-set-kkt"}中的最后两个条件)等价于求解一个等式约束的二次规划问题,
如果我们事先知道$\mathcal{A}({\bm{x}}^*).$ 但通常这是不可能的,
因此不能通过求解等式问题[\[eq:qp-active-set-1\]](#eq:qp-active-set-1){reference-type="eqref"
reference="eq:qp-active-set-1"}来求解原二次规划问题[\[eq:quadratic-programming-1\]](#eq:quadratic-programming-1){reference-type="eqref"
reference="eq:quadratic-programming-1"}.

在积极集法中, 依据以上的观察,
我们将积极集$\mathcal{A}$确定的约束看作等式约束, 而暂时忽略其余约束条件,
并通过某种迭代的方式不断修正调整这个集合,
直到识别出原问题[\[eq:quadratic-programming-1\]](#eq:quadratic-programming-1){reference-type="eqref"
reference="eq:quadratic-programming-1"}的解处的正确的积极约束.

更具体来说, 在第$k$次迭代, 我们从可行点${\bm{x}}^{(k)}$出发,
积极集为$\mathcal{A} = \mathcal{A}({\bm{x}}^{(k)}).$ 在这一步迭代中,
我们求解等式问题[\[eq:qp-active-set-1\]](#eq:qp-active-set-1){reference-type="eqref"
reference="eq:qp-active-set-1"}.
更方便的做法是将原点平移到${\bm{x}}^{(k)},$
令${\bm{s}} = {\bm{x}} - {\bm{x}}^{(k)},$ 求解问题
$$\label{eq:qp-active-set-2}
\begin{array}{cl}
\text{minimize} & \frac{1}{2} {\bm{s}}^T G {\bm{s}} + \left( {\bm{g}}^{(k)} \right)^T {\bm{s}}, \\
\text{subject to} & {\bm{a}}_i^T {\bm{s}} = 0, ~ i \in \mathcal{A},
\end{array}$$
其中${\bm{g}}^{(k)} = \nabla q({\bm{x}}^{(k)}) = G {\bm{x}}^{(k)} + {\bm{d}}$是原二次规划问题[\[eq:quadratic-programming-1\]](#eq:quadratic-programming-1){reference-type="eqref"
reference="eq:quadratic-programming-1"}的目标函数$q({\bm{x}})$在点${\bm{x}}^{(k)}$处的梯度向量.
这个问题是一个等式约束的二次规划问题,
可以用上一小节§[1.1](#subsec:7.2.1){reference-type="ref"
reference="subsec:7.2.1"} 中介绍的等式约束二次规划问题的求解方法进行求解.

记问题[\[eq:qp-active-set-2\]](#eq:qp-active-set-2){reference-type="eqref"
reference="eq:qp-active-set-2"}的解为${\bm{s}}^{(k)},$
我们需要对各种可能的情况进行分类讨论.

如果${\bm{s}}^{(k)} = {\bm{0}},$
即${\bm{x}}^{(k)}$是当前等式约束问题[\[eq:qp-active-set-1\]](#eq:qp-active-set-1){reference-type="eqref"
reference="eq:qp-active-set-1"}的解,
那么可以根据式[\[eq:general-elim-lagrange\]](#eq:general-elim-lagrange){reference-type="eqref"
reference="eq:general-elim-lagrange"}计算积极约束的拉格朗日乘子,
记为${\bm{\lambda}}^{(k)},$ 即有 $$\label{eq:qp-active-set-lambda}
{\bm{g}}^{(k)} + \sum\limits_{i \in \mathcal{A}} \lambda_i^{(k)} \alpha_i = 0.$$
此时, 除对偶可行性条件$\lambda_i \geqslant 0, ~ i \in \mathcal{I}$以外,
其余KKT条件[\[eq:qp-active-set-kkt\]](#eq:qp-active-set-kkt){reference-type="eqref"
reference="eq:qp-active-set-kkt"}均满足.
当与原问题不等式约束对应的拉格朗日乘子均非负, 即
$$\lambda_i^{(k)} \geqslant 0, ~ \forall i \in \mathcal{I} \cap \mathcal{A},$$
则$x^{(k)}$是原问题的KKT点, 迭代结束, 求解完毕. 若不然, 则设
$$\label{eq:qp-active-set-inactive-index}
q = \mathop{\mathrm{arg\,min}}_{i \in \mathcal{I} \cap \mathcal{A}} \lambda_i^{(k)},$$
有$\lambda_q^{(k)} < 0.$
令$\mathcal{A} = \mathcal{A} \setminus \{ q \},$
代入问题[\[eq:qp-active-set-2\]](#eq:qp-active-set-2){reference-type="eqref"
reference="eq:qp-active-set-2"}, 并求解此新的子问题.

如果${\bm{s}}^{(k)} \neq {\bm{0}},$
即${\bm{x}}^{(k)}$不是当前等式约束问题[\[eq:qp-active-set-1\]](#eq:qp-active-set-1){reference-type="eqref"
reference="eq:qp-active-set-1"}的解, 那么我们进一步检验试探点
$$\label{eq:qp-active-set-test-point}
\bar{{\bm{x}}}^{(k)} = {\bm{x}}^{(k)} + {\bm{s}}^{(k)}$$
是否满足其他不在积极集$\mathcal{A}$中的不等式约束条件. 如果都满足的话,
我们令 $$\label{eq:qp-active-set-next-step-1}
{\bm{x}}^{(k+1)} = \bar{{\bm{x}}}^{(k)} = {\bm{x}}^{(k)} + {\bm{s}}^{(k)},$$
积极集保持不变, 进入下一步迭代搜索. 如果不然,
即存在指标$i \not\in \mathcal{A},$ 使得
$${\bm{a}}_i^T {\bm{x}}^{(k)} - {\bm{b}} + {\bm{s}}^{(k)} > 0$$
此时试探点$\bar{{\bm{x}}}^{(k)}$不是原问题的可行点,
需要将其投影到原问题的可行域.
我们沿着方向${\bm{p}}^{(k)} = {\bm{s}}^{(k)}$进行线搜索,
选一个小于$1$但尽可能大的步长 $$\label{eq:qp-active-set-step-len-1}
\begin{aligned}
\bar{\alpha}_k & = \max \left\{ \alpha : ~ \alpha > 0, ~ {\bm{a}}_i^T {\bm{x}}^{(k)} - b_i + \alpha {\bm{a}}_i^T {\bm{p}}^{(k)} \leqslant 0, ~ {\bm{a}}_i^T {\bm{p}}^{(k)} > 0, ~ \forall i \not\in \mathcal{A} \right\}, \\
& = \min_{\substack{i: i \not\in \mathcal{A} \\ {\bm{a}}_i^T {\bm{p}}^{(k)} > 0}} \frac{b_i - {\bm{a}}_i^T {\bm{x}}^{(k)}}{{\bm{a}}_i^T {\bm{p}}^{(k)}}.
\end{aligned}$$ 注意$\bar{\alpha}_k$是严格小于$1$的,
因为此时的试探点${\bm{x}}^{(k)} + {\bm{s}}^{(k)}$不可行.
取指标$j$使得第$j$个约束取得上式中$\frac{b_i - {\bm{a}}_i^T {\bm{x}}^{(k)}}{{\bm{a}}_i^T {\bm{p}}^{(k)}}$的最大值,
即 $$\label{eq:qp-active-set-step-len-2}
\bar{\alpha}_k = \frac{b_j - {\bm{a}}_j^T {\bm{x}}^{(k)}}{{\bm{a}}_j^T {\bm{p}}^{(k)}}, ~~ j = \mathop{\mathrm{arg\,min}}_{\substack{i: i \not\in \mathcal{A} \\ {\bm{a}}_i^T {\bm{p}}^{(k)} > 0}} \frac{b_i - {\bm{a}}_i^T {\bm{x}}^{(k)}}{{\bm{a}}_i^T {\bm{p}}^{(k)}},$$
并称指标$j$对应的约束为阻滞(blocking)约束. 我们取
$$\label{eq:qp-active-set-next-step-2}
{\bm{x}}^{(k+1)} = {\bm{x}}^{(k)} + \bar{\alpha}_k {\bm{s}}^{(k)}, ~~ \mathcal{A} \gets \mathcal{A} \cup \{ j \}.$$
这里我们把指标$j$添加到了积极集, 是因为对于${\bm{x}}^{(k+1)},$
指标$j$对应的约束条件等式成立, 非积极约束$j$变成积极的.
迭代格式[\[eq:qp-active-set-next-step-1\]](#eq:qp-active-set-next-step-1){reference-type="eqref"
reference="eq:qp-active-set-next-step-1"}与[\[eq:qp-active-set-next-step-2\]](#eq:qp-active-set-next-step-2){reference-type="eqref"
reference="eq:qp-active-set-next-step-2"}可以统一表述为:
以等式问题[\[eq:qp-active-set-2\]](#eq:qp-active-set-2){reference-type="eqref"
reference="eq:qp-active-set-2"}的解${\bm{s}}^{(k)}$为搜索方向${\bm{p}}^{(k)},$
以 $$\label{eq:qp-active-set-step-len-uniform}
\alpha_k = \min (1, \bar{\alpha}_k)$$ 为迭代步长, 得下一步的迭代点
$$\label{eq:qp-active-set-next-step-uniform}
{\bm{x}}^{(k+1)} = {\bm{x}}^{(k)} + \alpha_k {\bm{p}}^{(k)},$$
同时根据是否有阻滞约束确定是否更新积极集: 若$\alpha_k < 1,$
有约束阻滞$j,$ 则将$j$添加到积极集$\mathcal{A}.$

至此, 我们可以以伪代码的形式,
将经典的积极集法表示为算法[\[algo:active-set\]](#algo:active-set){reference-type="ref"
reference="algo:active-set"}.

::: algorithm
随机选取一个可行点${\bm{x}}^{(0)},$ 确定相应的积极集$\mathcal{A},$
$k \gets 0.$
:::

::: {#eg:qp-active-set-algo .eg}
**例 3**. *我们来举例说明用积极集法求解含不等式约束的二次规划问题的步骤.
考虑如下的二次规划问题 $$\begin{array}{cl}
\text{minimize} & q({\bm{x}}) = (x_1 - 1)^2 + (x_2 - 2.5)^2, \\
\text{subject to} & -x_1 + 2x_2 - 2 \leqslant 0, \\
& x_1 + 2x_2 - 6 \leqslant 0, \\
& x_1 - 2x_2 - 2 \leqslant 0, \\
& -x_1 \leqslant 0, \\
& -x_2 \leqslant 0,
\end{array}$$ 我们用积极集法进行求解.
该问题的可行域可见图[\[fig:active-set-eg\]](#fig:active-set-eg){reference-type="ref"
reference="fig:active-set-eg"}, 由其中实线以及坐标轴围成的阴影区域.*

*我们选取初始点${\bm{x}}^{(0)} = (2, 0)^T,$
用$1$至$5$依次作为约束条件的指标. 在初始点${\bm{x}}^{(0)}$处,
约束$3$和$5$满足等式关系, 是积极约束,
所以初始积极集$\mathcal{A} = \{ 3, 5 \}.$
当前需要求解的等式问题[\[eq:qp-active-set-1\]](#eq:qp-active-set-1){reference-type="eqref"
reference="eq:qp-active-set-1"}为 $$\begin{array}{cl}
\text{minimize} & q({\bm{x}}) = (x_1 - 1)^2 + (x_2 - 2.5)^2, \\
\text{subject to} & x_1 - 2x_2 - 2 = 0, \\
& x_2 = 0,
\end{array}$$
或者经过平移的问题[\[eq:qp-active-set-2\]](#eq:qp-active-set-2){reference-type="eqref"
reference="eq:qp-active-set-2"} $$\begin{array}{cl}
\text{minimize} & q({\bm{s}}) = (s_1 + 1)^2 + (s_2 - 2.5)^2, \\
\text{subject to} & s_1 - 2s_2 = 0, \\
& s_2 = 0.
\end{array}$$
很容易看到${\bm{x}}^{(0)}$(即${\bm{s}}^{(0)} = {\bm{0}}$)是该问题的解.
由式[\[eq:qp-active-set-lambda\]](#eq:qp-active-set-lambda){reference-type="eqref"
reference="eq:qp-active-set-lambda"}求解积极约束的拉格朗日乘子,
即求解方程组
$$\begin{bmatrix} 1 \\ -2 \end{bmatrix} \lambda_3^{(0)} + \begin{bmatrix} 0 \\ -1 \end{bmatrix} \lambda_5^{(0)} = \begin{bmatrix} -2 \\ 5 \end{bmatrix},$$
得$\lambda_3^{(0)} = -2, \lambda_5^{(0)} = -1.$
指标$q = \mathop{\mathrm{arg\,min}}\limits_{i \in \{ 3, 5 \}} \lambda_i^{(0)} = 3.$
由于$\lambda_q^{(0)} = \lambda_3^{(0)} = -2 < 0,$
因此我们置${\bm{x}}^{(1)} = {\bm{x}}^{(0)},$
同时将指标$q = 3$从积极集中删去, 进入下一步迭代.
接下来我们需要求解等式问题 $$\begin{array}{cl}
\text{minimize} & q({\bm{s}}) = (s_1 + 1)^2 + (s_2 - 2.5)^2, \\
\text{subject to} & s_2 = 0.
\end{array}$$ 容易解得${\bm{s}}^{(1)} = (-1, 0)^T.$ 此时,
试探点${\bm{x}}^{(1)} + {\bm{s}}^{(1)} = (1, 0)^T$是可行点,
置${\bm{x}}^{(2)} = {\bm{x}}^{(1)} + {\bm{s}}^{(1)} = (1, 0)^T,$
同时积极集$\mathcal{A} = \{ 5 \}$保持不变, 进入下一步迭代.
容易验证${\bm{x}}^{(2)}$是这一步要解的等式问题的可行解,
进而可计算得相应的积极约束的拉格朗日乘子$\lambda_5^{(2)} = -5.$ 此时,
约束5变成非积极的, 积极集$\mathcal{A} = \emptyset$变为空集.
置${\bm{x}}^{(3)} = {\bm{x}}^{(2)} = (1, 0)^T$进入下一步迭代.
再次求解当前的等式约束问题(实际上已成为无约束问题)
$$\text{minimize} ~~ q({\bm{s}}) = s_1^2 + (s_2 - 2.5)^2,$$
得解${\bm{s}}^{(3)} = (0, 2.5)^T.$
试探点${\bm{x}}^{(3)} + {\bm{s}}^{(3)} = (1, 2.5)^T$不是可行点,
因此需要以${\bm{p}}^{(3)} = {\bm{s}}^{(3)}$为方向进行线搜索${\bm{x}}^{(3)} + \alpha_3 {\bm{p}}^{(3)},$
并由式[\[eq:qp-active-set-step-len-2\]](#eq:qp-active-set-step-len-2){reference-type="eqref"
reference="eq:qp-active-set-step-len-2"}以及式[\[eq:qp-active-set-step-len-uniform\]](#eq:qp-active-set-step-len-uniform){reference-type="eqref"
reference="eq:qp-active-set-step-len-uniform"}算得最优步长$\alpha_3 = 0.6,$
以及相应阻滞约束的指标$j = 1.$
置${\bm{x}}^{(4)} = {\bm{x}}^{(3)} + \alpha_3 {\bm{p}}^{(3)} = (1, 1.5)^T,$
并将阻滞约束的指标$j = 1$添加到积极集得$\mathcal{A} = \{ 1 \},$
进入下一步迭代. 再次求解当前的等式问题,
得${\bm{s}}^{(4)} = (0.4, 0.2)^T.$
相应的试探点${\bm{x}}^{(4)} + {\bm{s}}^{(4)} = (1.4, 1.7)^T$可行,
于是得新的迭代点${\bm{x}}^{(5)} = (1.4, 1.7)^T.$
由于${\bm{x}}^{(5)}$是当前等式问题的可行点,
且解得积极约束的拉格朗日乘子$\lambda_1^{(5)} = 0.8 > 0,$
算法终止条件达成, 得原问题的最优解
$${\bm{x}}^* = {\bm{x}}^{(5)} = (1.4, 1.7)^T.$$*

*我们将每一步的数值结果总结在表[1](#tab:active-set-eg){reference-type="ref"
reference="tab:active-set-eg"} 中.
表中的"$\backslash$"表示当前迭代步不需要计算相应的量.*

::: {#tab:active-set-eg}
   *$k$*   *${\bm{x}}^{(k)}$*   *$\mathcal{A}$*   *${\bm{s}}^{(k)}$*   *$\bar{{\bm{x}}}^{(k)}$可行*   *${\bm{\lambda}}^{(k)}$*        *$q$*         *$\alpha_k$*        *$j$*        *$q({\bm{x}}^{(k)})$*
  ------- -------------------- ----------------- -------------------- ------------------------------ --------------------------- ---------------- ---------------- ---------------- -----------------------
   *$0$*      *$(2, 0)^T$*      *$\{ 3, 5 \}$*       *$(0, 0)^T$*             *$\checkmark$*          *$\lambda_3^{(0)} = -2$*        *$3$*            *$0$*        *$\backslash$*         *$7.25$*
                                                                                                      *$\lambda_5^{(0)} = -1$*                                                      
   *$1$*      *$(2, 0)^T$*        *$\{ 5 \}$*       *$(-1, 0)^T$*             *$\checkmark$*               *$\backslash$*         *$\backslash$*       *$1$*        *$\backslash$*         *$7.25$*
   *$2$*      *$(1, 0)^T$*        *$\{ 5 \}$*        *$(0, 0)^T$*             *$\checkmark$*          *$\lambda_5^{(2)} = -5$*        *$5$*            *$0$*        *$\backslash$*         *$6.25$*
   *$3$*      *$(1, 0)^T$*       *$\emptyset$*      *$(0, 2.5)^T$*              *$\times$*                 *$\backslash$*         *$\backslash$*      *$0.6$*           *$1$*              *$6.25$*
   *$4$*     *$(1, 1.5)^T$*       *$\{ 1 \}$*      *$(0.4, 0.2)^T$*           *$\checkmark$*               *$\backslash$*         *$\backslash$*       *$1$*        *$\backslash$*           *$1$*
   *$5$*    *$(1.4, 1.7)^T$*      *$\{ 1 \}$*        *$(0, 0)^T$*             *$\checkmark$*          *$\lambda_1^{(5)} = 0.8$*   *$\backslash$*   *$\backslash$*   *$\backslash$*          *$0.8$*

  : *积极集法求解例[例 3](#eg:qp-active-set-algo){reference-type="ref"
  reference="eg:qp-active-set-algo"} 数值结果*
:::
:::

我们接下来对积极集法进行进一步的说明, 以阐明一些尚未讨论的问题.
首先是关于初始可行点${\bm{x}}^{(0)}$以及初始积极集$\mathcal{A}$的选取.
我们可以利用人工变量法, 求解辅助问题得${\bm{x}}^{(0)}.$
$\mathcal{A}$可以选取为$\mathcal{A}({\bm{x}}^{(0)}),$
也可以选取$\mathcal{A}({\bm{x}}^{(0)})$的任何一个子集,
即$\mathcal{A} \subseteq \mathcal{A}({\bm{x}}^{(0)}).$
不同的选取方式会导致不同的迭代过程, 但算得的最终结果是一样的,
这里不作过多讨论.

其次是积极集法的有限终止性. 事实上,
我们已经从表[1](#tab:active-set-eg){reference-type="ref"
reference="tab:active-set-eg"} 中观察到了, 如果$\alpha_k \neq 0,$
那么目标函数值总是下降的. 更具体来说,
我们假设初始积极集$\mathcal{A}$中约束的梯度是线性无关的,
那么对积极集的修正策略可以确保梯度线性无关这一性质对后面所有的积极集都成立.
这样一来,
等式问题[\[eq:qp-active-set-2\]](#eq:qp-active-set-2){reference-type="eqref"
reference="eq:qp-active-set-2"}总是适定的. 一个重要的观察是,
如果某一步(设为第$k$步)的等式问题[\[eq:qp-active-set-2\]](#eq:qp-active-set-2){reference-type="eqref"
reference="eq:qp-active-set-2"}解为${\bm{s}}^{(k)} = {\bm{0}},$
那么当前点${\bm{x}}^{(k)}$是等式问题[\[eq:qp-active-set-1\]](#eq:qp-active-set-1){reference-type="eqref"
reference="eq:qp-active-set-1"}的全局极小值点. 如果算法在这一步未终止,
即当前点${\bm{x}}^{(k)}$不是原问题[\[eq:quadratic-programming-1\]](#eq:quadratic-programming-1){reference-type="eqref"
reference="eq:quadratic-programming-1"}的解,
此时会有某个积极约束变成非积极的,
下一步的搜索方向${\bm{p}}^{(k+1)}$会是目标函数$q({\bm{x}})$的严格下降方向.
那么此后迭代步中的目标函数值都会严格小于$q({\bm{x}}^{(k)}),$
进而知此后迭代步中的积极集不可能与当前第$k$步的积极集相同,
否则与${\bm{x}}^{(k)}$是等式问题[\[eq:qp-active-set-1\]](#eq:qp-active-set-1){reference-type="eqref"
reference="eq:qp-active-set-1"}的全局极小值点这一事实产生矛盾.
另一个重要观察是, 设${\bm{s}}^{(k)} = {\bm{0}},$
那么至多$n$步迭代之后的第$k'$步,
我们又会遇到${\bm{s}}^{(k')} = {\bm{0}}$是等式问题[\[eq:qp-active-set-2\]](#eq:qp-active-set-2){reference-type="eqref"
reference="eq:qp-active-set-2"}的解这种情况. 证明我们留给读者完成.
那么综上所述, 积极集法的每$n$次迭代中, 我们至少会遇到其中某一步,
使得当前点是等式问题[\[eq:qp-active-set-1\]](#eq:qp-active-set-1){reference-type="eqref"
reference="eq:qp-active-set-1"}的全局极小值, 并且在往后的迭代中,
这一步的积极集不会再次出现. 但是, 所有可能的积极集的数目是有限的,
也就是说, 积极集法必然在有限步之内以找到原问题的最优解而结束.

我们最后要强调的是, 积极集法的重要特征是: 积极集改变后,
要有效求解等式约束的二次规划问题[\[eq:qp-active-set-2\]](#eq:qp-active-set-2){reference-type="eqref"
reference="eq:qp-active-set-2"}. 在每次迭代中,
我们并不重新分解KKT矩阵($O(n3)$的运算复杂度),
而是根据积极集$\mathcal{A}$的变化更新算法所需要的因子($O(n2)$的运算复杂度).
不同的方法采用不同的矩阵分解以便这种更新是可行且有效的.

## 习题 {#习题 .unnumbered}

1.  写出二次规划问题 $$\begin{array}{cl}
    \text{maximize} & q({\bm{x}}) = 6x_1 + 4x_2 - 13 - x_1^2 - x_2^2, \\
    \text{subject to} & x_1 + x_2 \leqslant 3, \\
    & x_1, ~ x_2 \geqslant 0
    \end{array}$$ 的KKT条件. 写出上述规划问题的对偶问题,
    以及对偶问题的对偶问题.

2.  分别利用直接消元法以及正交分解法(广义消元法)求解如下的等式约束的二次规划问题,
    并检查它们得到的解是否相同. $$\begin{array}{cl}
    \text{minimize} & q(x) = x_1^2 + x_2^2 + x_3^2 + 2x_3, \\
    \text{subject to} & x_1 + 2x_2 - x_3 = 4, \\
    & x_1 - x_2 + x_3 = -2.
    \end{array}$$

3.  设${\bm{x}}^*$是二次规划问题[\[eq:quadratic-programming-1\]](#eq:quadratic-programming-1){reference-type="eqref"
    reference="eq:quadratic-programming-1"}的解,
    证明${\bm{x}}^*$也是下列线性规划问题的解 $$\begin{array}{cl}
    \text{minimize} & {\bm{x}}^T \left( G {\bm{x}}^* + {\bm{d}} \right), \\
    \text{subject to} & {\bm{a}}_i^T {\bm{x}} = b_i, ~~ i \in \mathcal{E}, \\
    & {\bm{a}}_i^T {\bm{x}} \leqslant b_i, ~~ i \in \mathcal{I},
    \end{array}$$

4.  令$A$是一个$n \times m$的列满秩的矩阵,
    求$n$维空间中一点${\bm{x}}_0$到$\{ x : ~ A^T {\bm{x}} = {\bm{b}} \}$最短距离的问题可以表述为如下的等式约束二次规划问题
    $$\begin{array}{cl}
    \text{minimize} & \frac{1}{2} \left( {\bm{x}} - {\bm{x}}_0 \right)^T \left( {\bm{x}} - {\bm{x}}_0 \right), \\
    \text{subject to} & A^T x = b,
    \end{array}$$ 证明该问题的解以及拉格朗日乘子为 $$\begin{aligned}
    {\bm{x}}^* & = {\bm{x}}_0 + A \left( A^T A \right)^{-1} \left( {\bm{b}} - A^T {\bm{x}}_0 \right), \\
    {\bm{\lambda}}^* & = \left( A^T A \right)^{-1} \left( {\bm{b}} - A^T {\bm{x}}_0 \right),
    \end{aligned}$$ 并进一步证明当$A = {\bm{a}}$是一个列向量时,
    从点${\bm{x}}_0$到$\mathcal{R}^n$中的超平面$\{ x : ~ {\bm{a}}^T {\bm{x}} = b \}$最短距离等于$\dfrac{\lvert b - {\bm{a}}^T {\bm{x}} \rvert}{\lVert {\bm{a}} \rVert_2}.$

5.  设$G$是$n$阶可逆方阵, $A$是$n\times m$的矩阵且$A$列满秩.
    令$K = \begin{bmatrix} G & A \\ A^T & 0 \end{bmatrix}.$ 证明

    -   方阵$K$可逆;

    -   方阵$K$的逆矩阵可以由式[\[eq:quadratic-programming-kkt-inv-1\]](#eq:quadratic-programming-kkt-inv-1){reference-type="eqref"
        reference="eq:quadratic-programming-kkt-inv-1"}给出.

6.  考虑等式约束的二次规划问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
    reference="eq:quadratic-programming-eq-constrained"}. 设$G$半正定,
    $A \neq 0$列满秩, $Z$为矩阵$A$的列零空间的一组基构成的矩阵,
    $K = \begin{bmatrix} G & A \\ A^T & 0 \end{bmatrix}$是KKT矩阵. 证明

    -   KKT矩阵$K$总是不定的;

    -   若进一步假设既约Hessian阵$Z^T G Z$正定, 则KKT矩阵$K$非奇异,
        并且其逆可由式[\[eq:quadratic-programming-kkt-inv-2\]](#eq:quadratic-programming-kkt-inv-2){reference-type="eqref"
        reference="eq:quadratic-programming-kkt-inv-2"}给出;

    -   若既约Hessian阵$Z^T G Z$正定, 该二次规划问题有唯一极小点,
        该点就是唯一的全局最优解;

    -   假设方程$A^⊤ {\bm{x}} = {\bm{b}}, ~ G {\bm{x}} + A {\bm{\lambda}} = − {\bm{d}}$有解.
        如果既约Hessian阵$Z^T G Z$半正定且奇异,
        那么问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
        reference="eq:quadratic-programming-eq-constrained"}有无穷多解;

    -   如果既约Hessian阵$Z^T G Z$不定,
        或者方程$G {\bm{x}} + A {\bm{\lambda}} = − {\bm{d}}$无解,
        则问题[\[eq:quadratic-programming-eq-constrained\]](#eq:quadratic-programming-eq-constrained){reference-type="eqref"
        reference="eq:quadratic-programming-eq-constrained"}的目标函数无界.

7.  证明二次规划问题[\[eq:quadratic-programming-1\]](#eq:quadratic-programming-1){reference-type="eqref"
    reference="eq:quadratic-programming-1"}的KKT条件[\[eq:quadratic-programming-kkt\]](#eq:quadratic-programming-kkt){reference-type="eqref"
    reference="eq:quadratic-programming-kkt"}与利用积极集表达的KKT条件[\[eq:qp-active-set-kkt\]](#eq:qp-active-set-kkt){reference-type="eqref"
    reference="eq:qp-active-set-kkt"}等价.

8.  分别以$\{3\},$ $\{5\},$ 以及$\emptyset$为初始积极集$\mathcal{A},$
    用积极集法求解例[例 3](#eg:qp-active-set-algo){reference-type="ref"
    reference="eg:qp-active-set-algo"} 中的不等式约束的二次规划问题.

9.  假设积极集法选取的初始积极集$\mathcal{A}$中约束的梯度是线性无关的,
    即$\{ {\bm{a}}_i: ~ i \in \mathcal{A} \}$是线性无关组,
    并且在随后某一步,
    根据式[\[eq:qp-active-set-step-len-2\]](#eq:qp-active-set-step-len-2){reference-type="eqref"
    reference="eq:qp-active-set-step-len-2"}添加了一条阻滞约束$j$到积极集$\mathcal{A}$中.
    证明阻滞约束的梯度${\bm{a}}_j$与该集合中其余向量线性无关.

    另一方面,
    当我们从积极集中按式[\[eq:qp-active-set-inactive-index\]](#eq:qp-active-set-inactive-index){reference-type="eqref"
    reference="eq:qp-active-set-inactive-index"}删去一条非积极约束,
    显然不会改变当前积极集中梯度向量组成线性无关组这一性质. 由此,
    由归纳法可知, 积极集法的迭代步中关于积极集的修正策略,
    可以移植保持积极集对应的约束的梯度向量组线性无关这一性质.

10. 假设在积极集法的第$k$步,
    我们按式[\[eq:qp-active-set-inactive-index\]](#eq:qp-active-set-inactive-index){reference-type="eqref"
    reference="eq:qp-active-set-inactive-index"}删去一条非积极约束$q.$
    证明下一步的搜索方向${\bm{p}}^{(k+1)}$是目标函数$q({\bm{x}})$的严格下降方向(即$\left( {\bm{p}}^{(k+1)} \right)^T {\bm{g}}^{(k+1)} < 0$),
    且该方向是严格可行的(即$\left( {\bm{p}}^{(k+1)} \right)^T {\bm{a}}_q > 0$).

11. 利用矩阵增加或者减少一行时的QR分解的更新方法编写积极集法[\[algo:active-set\]](#algo:active-set){reference-type="ref"
    reference="algo:active-set"} 的程序, 并用其求解二次规划问题
    $$\begin{array}{cl}
    \text{minimize} & x_1^2 + 2x_2^2 -2x_1x_2 - 2x_1 - 6x_2, \\
    \text{subject to} & -x_1 - x_2 \leqslant -2, \\
    & -x_1 + 2x_2 \leqslant 2, \\
    & -x_1 \leqslant 0,\\
    & -x_2 \leqslant 0.
    \end{array}$$ 选取三个初始点:
    一个在可行域的内部(例如$(3, 1)^T$)、一个在顶点(例如$(\frac{2}{3}, \frac{4}{3})^T$)、一个在可行域的边界上,
    但不是顶点(例如$(4, 0)^T$),

    提示：可以利用Matlab的`qr.m`, `qrinsert.m`和`qrdelete.m`函数,
    或者利用Python `scipy`软件包的`scipy.linalg.qr`,
    `scipy.linalg.qr_insert`和`scipy.linalg.qr_delete`函数.
