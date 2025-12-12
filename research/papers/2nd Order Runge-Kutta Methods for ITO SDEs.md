# Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations

Andreas Rößler  

---

## 1. Introduction

We consider a probability space $(\Omega, \mathcal{F}, \mathbb{P})$ with a filtration $(\mathcal{F}_t)_{t \ge 0}$ satisfying the usual conditions, and a time interval
$$
I = [t_0, T], \quad 0 \le t_0 < T < \infty.
$$

Let $X = (X_t)_{t \in I}$ be the $d$-dimensional solution of the Itô SDE
$$
X_t = X_{t_0}
+ \int_{t_0}^t a(s, X_s)\,ds
+ \sum_{j=1}^m \int_{t_0}^t b_j(s, X_s)\,dW_s^j,
$$
where $W = (W_t)_{t \ge 0}$ is an $m$-dimensional Wiener process, and
$$
a,\, b_j : I \times \mathbb{R}^d \to \mathbb{R}^d, \quad j = 1,\dots,m,
$$
are at least Lipschitz with linear growth in the state variable. Under these assumptions there exists a unique strong solution.

We are interested in **weak approximation**: given a sufficiently smooth test function $f$, approximate $\mathbb{E}[f(X_T)]$ using a discrete-time process $Y$.

Let the grid be
$$
I_h = \{t_0, t_1, \dots, t_N\}, \quad t_{n+1} - t_n = h_n, \quad h = \max_n h_n.
$$

**Definition (Weak order).**  
A numerical method with approximation $Y(T)$ converges weakly with order $p$ if, for all $f \in C_P^{2(p+1)}(\mathbb{R}^d, \mathbb{R})$, there exist $C_f > 0$ and $\delta_0 > 0$ such that
$$
\bigl|\mathbb{E}[f(X_T)] - \mathbb{E}[f(Y(T))]\bigr| \le C_f h^p \quad\text{for } 0 < h < \delta_0.
$$

The goal of the paper is to construct **derivative–free stochastic Runge–Kutta (SRK) methods** of weak order $2$ for Itô SDEs with **multidimensional** noise, where:

- The number of stages is **independent** of the noise dimension $m$.
- Only $2m - 1$ random variables are needed per step (vs. $m(m+1)/2$ in standard schemes).
- The methods are efficient and practical for high-dimensional problems.

---

## 2. General Stochastic Runge–Kutta Methods

A very general class of $s$–stage SRK methods for weak approximation is defined as follows.

Let $M$ be a finite index set with $\kappa = |M|$ and random variables $\theta_\nu(h)$, $\nu \in M$, satisfying the **moment condition**
$$
\mathbb{E}\!\left[
\theta_{\nu_1}^{p_1}(h)\cdots \theta_{\nu_\kappa}^{p_\kappa}(h)
\right]
= O\left(h^{(p_1 + \cdots + p_\kappa)/2}\right)
$$
for all nonnegative integers $p_1,\dots,p_\kappa$. This ensures that $\theta_\nu(h)$ behave like Brownian increments of order $\sqrt{h}$.

Let $Y_n \approx X_{t_n}$. The method is
$$
Y_{n+1} = Y_n
+ \sum_{i=1}^s z_i^{(0,0)}\, a\bigl(t_n + c_i^{(0,0)} h_n,\, H_i^{(0,0)}\bigr)
+ \sum_{i=1}^s \sum_{k=1}^m \sum_{\nu\in M}
z_i^{(k,\nu)}\, b_k\bigl(t_n + c_i^{(k,\nu)} h_n,\, H_i^{(k,\nu)}\bigr),
$$
where the **stage values** satisfy
$$
\begin{aligned}
H_i^{(k,\nu)}
&= Y_n
  + \sum_{j=1}^s Z_{ij}^{(k,\nu),(0,0)} 
      a\bigl(t_n + c_j^{(0,0)} h_n,\, H_j^{(0,0)}\bigr) \\
&\quad + \sum_{j=1}^s \sum_{r=1}^m \sum_{\mu\in M}
    Z_{ij}^{(k,\nu),(r,\mu)}\,
    b_r\bigl(t_n + c_j^{(r,\mu)} h_n,\, H_j^{(r,\mu)}\bigr),
\end{aligned}
$$
with $k = 0,1,\dots,m$ and $\nu \in M \cup \{0\}$.

The coefficients are linear in the random variables:
$$
z_i^{(0,0)} = \alpha_i h_n,\quad
z_i^{(k,\nu)} = \sum_{\iota\in M} \gamma_i^{(\iota)(k,\nu)}\,\theta_\iota(h_n),
$$
$$
Z_{ij}^{(k,\nu),(0,0)} = A_{ij}^{(k,\nu),(0,0)} h_n,\quad
Z_{ij}^{(k,\nu),(r,\mu)} = \sum_{\iota\in M}
  B_{ij}^{(\iota)(k,\nu),(r,\mu)}\,\theta_\iota(h_n).
$$

Collecting the coefficients into vectors and matrices, we define the **stage abscissae**
$$
c^{(k,\nu)} = A^{(k,\nu),(0,0)} e,\quad e = (1,\dots,1)^T.
$$

The method is:
- **Explicit** if $A_{ij}^{(k,\nu),(0,0)} = B_{ij}^{(\iota)(k,\nu),(r,\mu)} = 0$ for all $j \ge i$.
- **Implicit** otherwise.

In the special case $b_j \equiv 0$ for all $j$, the scheme reduces to a standard deterministic Runge–Kutta method.

---

## 3. Colored Rooted Tree Analysis

To characterize weak order, we use a **colored rooted tree** formalism.

### 3.1 Tree Types and Order

We consider trees in a set $T_S(\Delta)$ whose nodes are of three types:

- Root: $\gamma$
- Deterministic node: $\tau$
- Stochastic node: $\sigma_{j_k}$, with color $j_k \in \{1,\dots,m\}$

For a tree $t$, let:
- $l(t)$: number of nodes,
- $d(t)$: number of deterministic nodes,
- $s(t)$: number of stochastic nodes.

The **order** of a tree is
$$
\rho(t) = d(t) + \frac{1}{2} s(t), \quad \rho(\gamma) = 0.
$$

Examples of trees of order $2$ (notation simplified):

- $t_{2,4} = (\sigma_{j_1}, \sigma_{j_2})$,
- $t_{2,8} = (\{\tau\}_{j_1}, \sigma_{j_2})$,
- $t_{2,13} = (\sigma_{j_1}, \{\sigma_{j_3},\sigma_{j_4}\}_{j_2})$.

Trees can be constructed recursively:

- $(t_1,\dots,t_k)$: attach $t_1,\dots,t_k$ to the root $\gamma$,
- $[t_1,\dots,t_k]$: attach to a deterministic node $\tau$,
- $\{t_1,\dots,t_k\}_j$: attach to a stochastic node $\sigma_j$.

### 3.2 Elementary Differentials

Each tree $t$ corresponds to an **elementary differential** $F(t)(x)$, defined recursively:

- $F(\gamma)(x) = f(x)$,
- $F(\tau)(x) = a(x)$,
- $F(\sigma_j)(x) = b_j(x)$.

For composite trees:
$$
F\big((t_1,\dots,t_k)\big)(x)
= f^{(k)}(x)\big(F(t_1)(x),\dots,F(t_k)(x)\big),
$$
$$
F\big([t_1,\dots,t_k]\big)(x)
= a^{(k)}(x)\big(F(t_1)(x),\dots,F(t_k)(x)\big),
$$
$$
F\big(\{t_1,\dots,t_k\}_j\big)(x)
= b_j^{(k)}(x)\big(F(t_1)(x),\dots,F(t_k)(x)\big).
$$

### 3.3 Labelings and Coefficients

For each tree $t$, we consider:

- $LTS(\Delta)$: set of **monotonically labeled** trees (labels increase along paths from the root).
- $\alpha_\Delta(t)$: number of such labelings.

For Itô SDE expectations, only a subset of trees appear, denoted $T_S(I) \subset T_S(\Delta)$, where:

- stochastic nodes appear in **pairs** with the same color,
- the paired nodes are not directly connected.

The number of valid Itô labelings for a tree $t$ is denoted $\alpha_I(t)$. For trees not in $T_S(I)$, $\alpha_I(t) = 0$.

We also define a combinatorial factor $\beta(t)$ counting how many trees in $T_S(\Delta)$ are equivalent to a given Itô tree in $T_S(I)$, once noise index equalities/inequalities are fixed.

### 3.4 Itô–Taylor Expansion of the Exact Solution

Assume $f, a_i, b_{i,j}$ have sufficient smoothness and polynomial growth. Then for small $h$,
$$
\mathbb{E}_{t_0,x_0}[f(X_{t_0+h})]
= \sum_{\substack{t \in T_S(I)\\ \rho(t) \le p}}
  \sum_{j_1,\dots,j_{s(t)/2}=1}^m
  \frac{\alpha_I(t)}{2^{s(t)/2}\,\rho(t)!}
  F(t)(x_0)\, h^{\rho(t)}
+ O(h^{p+1}).
$$

### 3.5 Expansion of the Numerical Method

For the SRK method, the expectation expands as
$$
\mathbb{E}_{t_0,x_0}[f(Y_{t_0+h})]
= \sum_{\substack{t \in T_S(\Delta)\\ \rho(t) \le p+\frac{1}{2}}}
  \sum_{j_1,\dots,j_{s(t)}=1}^m
  \frac{\alpha_\Delta(t)\,\gamma(t)}{(l(t)-1)!}
  F(t)(x_0)\, \mathbb{E}[\Phi_S(t)]
+ O(h^{p+1}),
$$
where:
- $\gamma(t)$ is a density factor defined recursively from the tree shape,
- $\Phi_S(t)$ is an **elementary weight** determined by the SRK tableau and the random variables.

Matching coefficients between the exact and numerical expansions yields the **order conditions**.

---

## 4. Order Conditions for SRK Methods

Let $p \in \mathbb{N}$ and assume sufficient smoothness of $a_i$ and $b_{i,j}$. The SRK method has weak order $p$ if:

For every tree $t \in T_S(\Delta)$ with $\rho(t) \le p + \frac{1}{2}$, and for every pattern of equalities/inequalities between the noise indices $j_1,\dots,j_{s(t)}$, the following condition holds:
$$
\mathbb{E}[\Phi_S(t)]
= \frac{\alpha_I(t)\,(l(t)-1)!}{\alpha_\Delta(t)\,\beta(t)\,\gamma(t)\,2^{s(t)/2}\,\rho(t)!}\, h^{\rho(t)}.
$$

Additionally, to control moments of the numerical solution uniformly in the number of steps, it is sufficient to use **bounded** random variables satisfying the moment condition, and to enforce that the total weight of each stochastic increment is centered, e.g.
$$
\mathbb{E}\bigl[z^{(k,\nu)T} e\bigr] = 0
\quad\text{for all } 1 \le k \le m,\; \nu\in M.
$$

These conditions specialize to a large but finite system of algebraic equations for the coefficients $\alpha$, $A^{(i)}$, $B^{(i)}$, and the distributions of the random variables.

---

## 5. A New Class of Efficient Weak Order Two SRK Methods

### 5.1 Motivation: Inefficiency of Standard Weak Order 2 Schemes

A classical derivative–free weak order two SRK scheme for autonomous SDEs is PL1WM (Platen’s method). For step $n$ with $Y_n = Y(t_n)$, its update can be written in terms of drift $a$, diffusion components $b_k$, stage values $H^{(0)}, H^{(k)}_\pm$, and discrete random variables $\hat{I}^{(k)}$, $\hat{J}^{(k,\ell)}$.

This scheme has two drawbacks:

- Per step, it needs
  - $2$ evaluations of $a$,
  - $2m+1$ evaluations of each $b_k$,
  - $m(m+1)/2$ random variables,
- So the cost grows **at least quadratically** with $m$.

The new class of methods aims to replace certain derivative approximations in the Itô–Taylor expansion by more efficient **stochastic difference approximations**, so that:

- The number of stages $s$ is independent of $m$,
- Only $2m-1$ random variables are needed,
- Each $b_k$ is evaluated only a fixed number of times per step (independent of $m$).

### 5.2 Key Approximation Idea

A central term in the weak Taylor expansion is of the form
$$
\sum_{k,\ell=1}^m \frac{\partial b_k}{\partial x}(x_0)\, b_\ell(x_0)\, I^{(\ell,k)},
$$
where $I^{(\ell,k)}$ denotes iterated Itô integrals.

A naive finite-difference approximation,
$$
\frac{b_k(x_0 + b_\ell(x_0)\sqrt{h}) - b_k(x_0 - b_\ell(x_0)\sqrt{h})}{2\sqrt{h}}\,I^{(\ell,k)},
$$
requires $2m+1$ evaluations of each $b_k$ for each pair $(k,\ell)$, which is costly.

Instead, the paper proposes a more efficient *combined* stochastic approximation:
$$
\sum_{k=1}^m 
\left(
b_k\Bigl(x_0 + \sum_{\ell=1}^m b_\ell(x_0)\,\frac{I^{(\ell,k)}}{\sqrt{h}}\Bigr)
-
b_k\Bigl(x_0 - \sum_{\ell=1}^m b_\ell(x_0)\,\frac{I^{(\ell,k)}}{\sqrt{h}}\Bigr)
\right) \frac{\sqrt{h}}{2},
$$
which can be shown (by Taylor expansion) to match the required term while using **only 3 evaluations** of each $b_k$ (for fixed $k$), regardless of $m$.

This leads to the new SRK family.

### 5.3 Definition of the New SRK Family

We define an $s$–stage SRK method with
- Stage abscissae $c^{(0)}, c^{(1)}, c^{(2)}$,
- Coefficient matrices $A^{(0)}, A^{(1)}, A^{(2)}$,
- Diffusion matrices $B^{(0)}, B^{(1)}, B^{(2)}$,
- Weight vectors $\alpha$, $\beta^{(1)}, \beta^{(2)}, \beta^{(3)}, \beta^{(4)}$.

The update reads:
$$
\begin{aligned}
Y_{n+1}
&= Y_n
 + \sum_{i=1}^s \alpha_i\, a\bigl(t_n + c_i^{(0)}h_n, H_i^{(0)}\bigr)\,h_n \\
&\quad + \sum_{i=1}^s \sum_{k=1}^m \beta_i^{(1)}\, b_k\bigl(t_n + c_i^{(1)}h_n, H_i^{(k)}\bigr)\,\hat{I}^{(k)} \\
&\quad + \sum_{i=1}^s \sum_{k=1}^m \beta_i^{(2)}\, b_k\bigl(t_n + c_i^{(1)}h_n, H_i^{(k)}\bigr)\,
  \frac{\hat{I}^{(k,k)}}{\sqrt{h_n}} \\
&\quad + \sum_{i=1}^s \sum_{k=1}^m \beta_i^{(3)}\, b_k\bigl(t_n + c_i^{(2)}h_n, \hat{H}_i^{(k)}\bigr)\,\hat{I}^{(k)} \\
&\quad + \sum_{i=1}^s \sum_{k=1}^m \beta_i^{(4)}\, b_k\bigl(t_n + c_i^{(2)}h_n, \hat{H}_i^{(k)}\bigr)\,\sqrt{h_n}.
\end{aligned}
$$

The stage values are
$$
\begin{aligned}
H_i^{(0)} &= Y_n
+ \sum_{j=1}^s A_{ij}^{(0)}\, a\bigl(t_n + c_j^{(0)}h_n, H_j^{(0)}\bigr)\,h_n
+ \sum_{j=1}^s \sum_{\ell=1}^m B_{ij}^{(0)}\, b_\ell\bigl(t_n + c_j^{(1)}h_n, H_j^{(\ell)}\bigr)\,\hat{I}^{(\ell)},\$$4pt]
H_i^{(k)} &= Y_n
+ \sum_{j=1}^s A_{ij}^{(1)}\, a\bigl(t_n + c_j^{(0)}h_n, H_j^{(0)}\bigr)\,h_n
+ \sum_{j=1}^s B_{ij}^{(1)}\, b_k\bigl(t_n + c_j^{(1)}h_n, H_j^{(k)}\bigr)\,\sqrt{h_n},\$$4pt]
\hat{H}_i^{(k)} &= Y_n
+ \sum_{j=1}^s A_{ij}^{(2)}\, a\bigl(t_n + c_j^{(0)}h_n, H_j^{(0)}\bigr)\,h_n
+ \sum_{j=1}^s \sum_{\substack{\ell=1\\ \ell\neq k}}^m B_{ij}^{(2)}\, b_\ell\bigl(t_n + c_j^{(1)}h_n, H_j^{(\ell)}\bigr)\,
  \frac{\hat{I}^{(k,\ell)}}{\sqrt{h_n}}.
\end{aligned}
$$

#### Discrete Random Variables

We use $2m - 1$ independent random variables per step:

- $\hat{I}^{(k)}$, $k = 1,\dots,m$, with moments
  $$
  \mathbb{E}[\hat{I}^{(k)}] = 0,\quad
  \mathbb{E}\bigl[(\hat{I}^{(k)})^2\bigr] = h,\quad
  \mathbb{E}\bigl[(\hat{I}^{(k)})^3\bigr] = 0,\quad
  \mathbb{E}\bigl[(\hat{I}^{(k)})^4\bigr] = 3h^2,
  $$
  and higher moments of order $O(h^{q/2})$.

  A simple choice is a three–point distribution:
  $$
  \mathbb{P}\bigl(\hat{I}^{(k)} = \pm\sqrt{3h}\bigr) = \tfrac{1}{6},
  \quad
  \mathbb{P}\bigl(\hat{I}^{(k)} = 0\bigr) = \tfrac{2}{3}.
  $$

- $\tilde{I}^{(k)}$, $k = 1,\dots,m-1$, with
  $$
  \mathbb{E}[\tilde{I}^{(k)}] = 0,\quad
  \mathbb{E}\bigl[(\tilde{I}^{(k)})^2\bigr] = h,
  $$
  e.g. a symmetric two–point distribution $\pm\sqrt{h}$ with probability $1/2$.

The mixed random variables $\hat{I}^{(k,\ell)}$ are defined by
$$
\hat{I}^{(k,\ell)} =
\begin{cases}
\frac{1}{2}\bigl(\hat{I}^{(k)}\hat{I}^{(\ell)} - \sqrt{h}\,\tilde{I}^{(k)}\bigr), & k < \ell,\$$4pt]
\frac{1}{2}\bigl(\hat{I}^{(k)}\hat{I}^{(\ell)} + \sqrt{h}\,\tilde{I}^{(\ell)}\bigr), & \ell < k,\$$4pt]
\frac{1}{2}\bigl((\hat{I}^{(k)})^2 - h\bigr), & k = \ell,
\end{cases}
$$
which reproduce the required Itô integral moments up to the needed order.

---

## 6. Weak Order Two Conditions for the New SRK Family

Using the colored rooted tree analysis described earlier, one can derive all weak order two conditions for the coefficients of the new SRK family.

Let
- $c^{(i)} = A^{(i)} e$, $i = 0,1,2$,
- $\cdot$ between vectors denote componentwise (Hadamard) products.

### 6.1 Order 1 Conditions

To achieve **weak order 1**, the coefficients must satisfy:

1. $\alpha^T e = 1$.
2. $\beta^{(4)T} e = 0$.
3. $\beta^{(3)T} e = 0$.
4. $(\beta^{(1)T} e)^2 = 1$.
5. $\beta^{(2)T} e = 0$.
6. $\beta^{(1)T} B^{(1)} e = 0$.
7. $\beta^{(4)T} A^{(2)} e = 0$.
8. $\beta^{(3)T} B^{(2)} e = 0$.
9. $\beta^{(4)T} \bigl(B^{(2)} e \cdot B^{(2)} e\bigr) = 0$.

### 6.2 Additional Conditions for Weak Order 2

Assuming more smoothness, **weak order 2** requires, in addition:

10. $\alpha^T A^{(0)} e = \tfrac{1}{2}$.
11. $\alpha^T \bigl(B^{(0)} e \cdot B^{(0)} e\bigr) = \tfrac{1}{2}$.
12. $(\beta^{(1)T} e)\bigl(\alpha^T B^{(0)} e\bigr) = \tfrac{1}{2}$.
13. $(\beta^{(1)T} e)\bigl(\beta^{(1)T} A^{(1)} e\bigr) = \tfrac{1}{2}$.
14. $\beta^{(3)T} A^{(2)} e = 0$.
15. $\beta^{(2)T} B^{(1)} e = 1$.
16. $\beta^{(4)T} B^{(2)} e = 1$.

Plus a number of higher-order algebraic relations (quadratic and cubic in the $B^{(i)}$-rows), ensuring consistency of all trees up to order $2.5$. In total, for non-commutative noise and $m>1$, this yields **59 conditions**.

For scalar noise ($m=1$) and $A^{(2)} = 0$, the system reduces to **28 conditions**, and explicit SRK methods with $s \ge 3$ stages can satisfy them.

---

## 7. Concrete Schemes and Efficiency

By solving the order conditions, one obtains explicit 3–stage SRK schemes of weak order 2:

- Some have deterministic order $(p_D, p_S) = (3,2)$, i.e. order 3 for ODEs and order 2 for SDEs.
- Another scheme, denoted RI6, has $(p_D, p_S) = (2,2)$.

All these schemes need:

- $3$ drift evaluations (RI1, RI3, RI5) or $2$ drift evaluations (RI6) per step.
- Exactly $5$ evaluations of each diffusion function $b_k$ per step, **independent of $m$**.
- Only $2m-1$ random variables per time step.

This is significantly cheaper than classical order 2 methods whose stage count and random-variable count scale with $m$.

---

## 8. Numerical Examples

Several numerical experiments validate the theoretical order and compare efficiency against:

- Euler–Maruyama (EM) (order 1),
- Extrapolated Euler (ExEu) (weak order 2 via Richardson extrapolation),
- Platen’s weak order 2 SRK scheme (PL1WM),
- The new schemes (RI5, RI6).

### 8.1 Linear 2D System with Commutative Noise

For a linear system
$$
dX_t
= A X_t\,dt + B X_t\,dW_t,
$$
with explicitly known moments $\mathbb{E}[X_t]$ and $\mathbb{E}(X_t X_t^T)$, the error
$$
\mu_h = \mathbb{E}[f(X_T)] - \frac{1}{M}\sum_{k=1}^M f(Y_T^{(k)})
$$
is estimated over many Monte Carlo trajectories and plotted versus stepsize $h$.

The new schemes RI5 and RI6 show the expected second-order slope in log–log graphs, with smaller error constants than extrapolated Euler and PL1WM.

### 8.2 Nonlinear 2D System with Non-Commutative Noise

For a nonlinear system with cross-coupled drift and diffusion, exact moments can still be computed. The tests confirm:

- Weak order 2 for RI5 and RI6.
- Clear efficiency gains in terms of computational effort vs. error.

### 8.3 4D Nonlinear System with Up to 6 Noise Dimensions

A more complex system with
- dimension $d=4$,
- up to $m=6$ independent Wiener processes,

is used to demonstrate the **scaling with noise dimension**. As $m$ increases from 2 to 6:

- The cost of PL1WM grows roughly like $m^2$,
- The cost of RI5/RI6 grows only linearly in $m$ (through the $2m-1$ random variables),
- For a fixed error tolerance, the new schemes are substantially faster.

---

## 9. Conclusion

The paper introduces a new class of explicit stochastic Runge–Kutta methods of weak order 2 for Itô SDEs with multidimensional noise, characterized by:

- **Stage count independent of the noise dimension $m$**,
- **Only $2m-1$ random variables per step**,
- **Fixed number of diffusion evaluations per step**, independent of $m$,
- The possibility of **higher deterministic order** (e.g. order 3 for ODEs),
- Good practical performance in numerical experiments.

These schemes are particularly attractive for high-dimensional SDEs, such as those in mathematical finance or physics, where traditional weak order 2 schemes become prohibitively expensive.

Future work mentioned in the paper includes:

- Development of higher-order SRK methods in this framework,
- Stability analysis,
- Implicit versions,
- Embedded error estimators and adaptive step-size control.

---

## References

(Labels [1]–[24] correspond to the original article’s bibliography, including works by Burrage, Kloeden–Platen, Milstein–Tretyakov, Talay–Tubaro, Tocino–Vigo-Aguiar, etc.)

