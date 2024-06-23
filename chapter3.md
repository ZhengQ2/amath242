## Chapter 3: Numerical linear algebra
### §3.1 Introduction

<u>Problem</u>: Solve for $\vec{x}$, given the linear system $A\vec{x} = \vec{b}$, with $A \in \mathbb{R}^{n \times n}$, $\vec{x}, \vec{b} \in \mathbb{R}^n$.
$$\begin{pmatrix}a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n}\\\vdots & \vdots & \ddots & \vdots \\ a_{n1} & a_{n2} & \cdots & a_{nn} \end{pmatrix} \begin{pmatrix}x_1\\x_2\\\vdots\\x_n\end{pmatrix}=\begin{pmatrix}b_1\\b_2\\\vdots\\b_n\end{pmatrix}$$

> Definition 3.1: (Determinant)
> The determinant of a  matrix $A \in \mathbb{R}^{n \times n}$ is given by
> $$\det(A) = \sum_{j=1}^n (-1)^{i+j}a_{ij}\det(A_{ij})$$
> - $i$ can be any number from 1 to $n$, then $a_{ij} (j = 1, 2, \ldots, n)$ comes from the $i$-th row of $A$.
> - $A_{ij}$ is the $(n-1) \times (n-1)$ matrix obtained by removing the $i$-th row and $j$-th column from the original matrix $A$.

> Example 3.1: 
> $$\begin{align*}&\det\begin{pmatrix}1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9\end{pmatrix} \\ =& (-1)^{1+1} a_{11} \det\begin{pmatrix}5 & 6 \\ 8 & 9\end{pmatrix} \\&+ (-1) ^{1+2}a_{12} \det\begin{pmatrix}4 & 6 \\ 7 & 9\end{pmatrix} \\&+ (-1)^{1+3}a_{13} \det\begin{pmatrix}4 & 5 \\ 7 & 8\end{pmatrix} \\ =& (-1)^{1+1} \times 1 \times (5 \times 9 - 6 \times 8) \\&+ (-1)^{1+2} \times 2 \times (4 \times 9 - 6 \times 7) \\&+ (-1)^{1+3} \times 3 \times (4 \times 8 - 5 \times 7) \\ =& 0\end{align*}$$

> Theorem 3.1: (Existance and uniqueness of solution of $A\vec{x} = \vec{b}$)
> - Case 1: If $\det(A) \neq 0$, then $\vec{x} = A^{-1}\vec{b}$ is the unique solution of $A\vec{x} = \vec{b}$.
> - Case 2: If $\det(A) = 0$, then
>    - Case 2a: If $\vec{b} \in \text{range}(A)$, then $A\vec{x} = \vec{b}$ has infinitely many solutions.
>    - Case 2b: If $\vec{b} \notin \text{range}(A)$, then $A\vec{x} = \vec{b}$ has no solution.

### §3.2 Gaussian elimination
#### §3.2.1 LU factorization
> Algorithm 3.1: (Gaussian elimination)
> - Phase 1: Reduce the matrix A to upper triangular form.
> - Phase 2: Solve the reduced system.

> Definition 3.2: A matrix $A \in \mathbb{R}^{n \times n}$ with components $a_{ij}$ is said to be
> - upper triangular if $a_{ij} = 0$ for $i > j$,
> - lower triangular if $a_{ij} = 0$ for $i < j$,
>
> $A \in \mathbb{R}^{n \times n}$ is said to be a triangular matrix if it is either upper or lower triangular.
> Solving these triangular systems are relatively straightforward.

> Algorithm 3.2: (Forward and backward substitution)
> $A\vec{x} = \vec{b}, A \in \mathbb{R}^{n \times n}, \vec{b} \in \mathbb{R}^{n\times 1}$n
> - A: upper triangular $\Rightarrow$ backward substitution)
> ![image3.1](images/image3.1.png)
> $$\begin{align*}x_n &= a_{n,n}^{-1} b_n \\ x_{n-1} &= a_{n-1,n-1}^{-1} (b_{n-1} - a_{n-1,n}x_n) \\ x_{n-2} &= a_{n-2,n-2}^{-1} (b_{n-2} - a_{n-2,n}x_n - a_{n-2,n-1}x_{n-1}) \\ &\vdots \\ x_i &= a_{ii}^{-1} (b_i - \sum_{k=i+1}^n a_{ik}x_k) \end{align*}$$
>
> - A: lower triangular $\Rightarrow$ forward substitution)
> ![image3.2](images/image3.2.png)
> $$\begin{align*}x_1 &= a_{11}^{-1} b_1 \\ x_2 &= a_{22}^{-1} (b_2 - a_{21}x_1) \\ x_3 &= a_{33}^{-1} (b_3 - a_{31}x_1 - a_{32}x_2) \\ &\vdots \\ x_i &= a_{ii}^{-1} (b_i - \sum_{k=1}^{i-1} a_{ik}x_k) \end{align*}$$
> Now we see how a dense matrix $A$ is reduced to an (upper-)triangular form: (with an example first)

> Example 3.2: Consider the system $A\vec{x} = \vec{b}$ with $$A = \begin{pmatrix}1 & 2 & 3 \\ \textcolor{red}{4} & 5 & 6 \\ \textcolor{red}{7} & \textcolor{blue}{8} & 1\end{pmatrix}, \vec{b} = \begin{pmatrix}1 \\ 2 \\ 3\end{pmatrix}$$
> <span style="color:red">1st</span>, <span style="color:blue">2nd</span> columns are to be eliminated.
> - Step 1: Denote $A^{(1)} = \begin{pmatrix}1 & 2 & 3 \\ \textcolor{red}{4} & 5 & 6 \\ \textcolor{red}{7} & 8 & 1\end{pmatrix}$
><span style="color:red"> 1st: eliminate this sub-column $\Rightarrow A^{(2)}$</span>.
> We use $a_{11}^{(1)}$ as the point element to compute $A^{(2)}$ by eliminating $a_{21}^{(1)}$ and $a_{31}^{(1)}$.
> $$A^{(1)} = \begin{pmatrix}1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 1\end{pmatrix} \begin{matrix} \\ R2 + (-4)R1 \rightarrow R2 \\ R3 + (-7)R1 \rightarrow R3 \end{matrix} $$
> >"To apply elementary row operations to $A$, we multiply $A$ by a matrix on the left: $EA$; and $E$ for any row operation is obtained by applying the same row operation to the identity matrix."
>
> $\begin{pmatrix}1 & 0 & 0 \\ -4 & 1 & 0 \\ -7 & 0 & 1\end{pmatrix} \cdot \begin{pmatrix}1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 1\end{pmatrix} = \begin{pmatrix}1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & -6 & -20\end{pmatrix}$
> - In general for $\mathbb{R}^{3 \times 3}$, we write $M_1 = \begin{pmatrix}1 & 0 & 0 \\ -a_{21}^{(1)}/a_{11}^{(1)} & 1 & 0 \\ -a_{31}^{(1)}/a_{11}^{(1)} & 0 & 1\end{pmatrix}$, then $A^{(2)} = \begin{pmatrix}a_{11}^{(2)} & a_{12}^{(2)} & a_{13}^{(2)} \\ 0 & a_{22}^{(2)} & a_{23}^{(2)} \\ 0 & a_{32}^{(2)} & a_{33}^{(2)}\end{pmatrix}$
> $\implies M_1A^{(1)} = A^{(2)}$
>
> Step 2: We have
> $$A^{(2)} = \begin{pmatrix}1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & -6 & -20\end{pmatrix} \begin{matrix} \\ \\ R3 + (-2)R2 \rightarrow R3 \end{matrix} = \begin{pmatrix}1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & 0 & -8\end{pmatrix}$$
> - In general for $\mathbb{R}^{3 \times 3}$, we write $M_2 = \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -a_{32}^{(2)}/a_{22}^{(2)} & 1\end{pmatrix}$, then $A^{(3)} =  \begin{pmatrix}a_{11}^{(2)} & a_{12}^{(2)} & a_{13}^{(2)} \\ 0 & a_{22}^{(2)} & a_{23}^{(2)} \\ 0 & 0 & a_{33}^{(2)}\end{pmatrix}$
> $\implies M_2A^{(2)} = A^{(3)}$
> - We also denote the final upper triangular matrix as $U = A^{(3)}$.
> 
> Summary of what happened so far:
> $$\begin{align*}A &= \begin{pmatrix}1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 1\end{pmatrix} \\ U &= \begin{pmatrix}1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & 0 & -8\end{pmatrix} \\ M_1 &= \begin{pmatrix}1 & 0 & 0 \\ -4 & 1 & 0 \\ -7 & 0 & 1\end{pmatrix} \\ M_2 &= \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -2 & 1\end{pmatrix} \\ M_2 & \cdot (M_1\cdot A) = U \\ (M_2 &\cdot M_1) \cdot A = U \\ A &= (M_2 \cdot M_1)^{-1} \cdot U  \\ A &= M_1^{-1} \cdot M_2^{-1} \cdot U \end{align*}$$
We denote $L_i = M_i^{-1}$, $L = L_1 \cdot L_2 = M_1^{-1} \cdot M_2^{-1}$, then $A = LU$.

Q: We know how to compute $M_1, M_2, U$, to compute $L$, we need:
- inverse of $M_1, M_2$;
- multiplication of $M_1^{-1}, M_2^{-1}$;

both are not trivial generally.

A: They turns out to be trivial due to the following "two strokes of luck":

> Proposition 3.1: (Inversion property)
> $L_i = M_i^{-1}$ can be obtained from $M_i$ by swapping the signs of the off-diagonal elements.
> Take $R^{3 \times 3}$ as an example:
> $$\begin{pmatrix}1 & 0 & 0 \\ c_1 & 1 & 0 \\ c_2 & 0 & 1\end{pmatrix} \cdot \begin{pmatrix}1 & 0 & 0 \\ -c_1 & 1 & 0 \\ -c_2 & 0 & 1\end{pmatrix} = \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1\end{pmatrix}$$
> $$\begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & c_3 & 1\end{pmatrix} \cdot \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -c_3 & 1\end{pmatrix} = \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1\end{pmatrix}$$

> Proposition 3.2: (Combination property)
> $L = L_1 \cdot L_2 \cdot L_3 \cdot \cdots \cdot L_{n-1} = M_1^{-1} \cdot M_2^{-1} \cdot M_3^{-1} \cdot \cdots \cdot M_{n-1}^{-1}$
> $\Rightarrow$ $L$ can be obtained by replacing all off-diagonal elements of $L_i$ with the corresponding positions in $L$.
> Take $R^{3 \times 3}$ as an example:
> $$\begin{pmatrix}1 & 0 & 0 \\ c_1 & 1 & 0 \\ c_2 & 0 & 1\end{pmatrix} \cdot \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & c_3 & 1\end{pmatrix} = \begin{pmatrix}1 & 0 & 0 \\ c_1 & 1 & 0 \\ c_2 & c_3 & 1\end{pmatrix}$$

> Example 3.2 (cont'd): 
> $$A = \begin{pmatrix}1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 1\end{pmatrix} U = \begin{pmatrix}1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & 0 & -8\end{pmatrix}$$
> $$M_1 = \begin{pmatrix}1 & 0 & 0 \\ -4 & 1 & 0 \\ -7 & 0 & 1\end{pmatrix} \implies L_1 = \begin{pmatrix}1 & 0 & 0 \\ 4 & 1 & 0 \\ 7 & 0 & 1\end{pmatrix}$$
> $$M_2 = \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -2 & 1\end{pmatrix} \implies L_2 = \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 2 & 1\end{pmatrix}$$
> $$L = L_1 \cdot L_2 = \begin{pmatrix}1 & 0 & 0 \\ 4 & 1 & 0 \\ 7 & 0 & 1\end{pmatrix} \cdot \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 2 & 1\end{pmatrix} = \begin{pmatrix}1 & 0 & 0 \\ 4 & 1 & 0 \\ 7 & 2 & 1\end{pmatrix}$$
> Therefore, we have LU-factorization of A:
> $$A = \begin{pmatrix}1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 1\end{pmatrix} = \begin{pmatrix}1 & 0 & 0 \\ 4 & 1 & 0 \\ 7 & 2 & 1\end{pmatrix} \begin{pmatrix}1 & 2 & 3 \\ 0 & -3 & -6 \\ 0 & 0 & -8\end{pmatrix} = LU$$

> Definition 3.3: (LU-factorization)
> For $A \in \mathbb{R}^{n \times n}$, LU-factorization may be computed as follows:
> $$\begin{align*} A^{(1)} &= A \\ A^{(2)} &= M_1A^{(1)} \\ A^{(3)} &= M_2A^{(2)} = M_2M_1A^{(1)} \\ &\vdots \\ A^{(n)} &= M_{n-1}A^{(n-1)} = M_{n-1}M_{n-2} \cdots M_1A^{(1)} = U \end{align*}$$
> $M_j$ looks like following:
> $$M_j = \begin{pmatrix}1 &&&&\\ & \ddots &&&\\ && 1 &&\\ &&&\ddots&\\ &&c_{ij} && 1\end{pmatrix}, c_{ij} = -\frac{a_{ij}^{(j)}}{a_{jj}^{(j)}}$$ for $j+1 \leq i \leq n$.
$$\begin{align*}M_{n-1}&M_{n-2} \cdots M_1A = U \\ \implies A &= M_1^{-1}M_2^{-1} \cdots M_{n-1}^{-1}U\\&= L_1L_2 \cdots L_{n-1}U \\&= L\cdot U\end{align*}$$
- $L$'s diagonal elements are all ones, which means $L$ is a unit lower triangular matrix.

> Algorithm 3.3: (Gaussian elimination (full version))
> - Phase 1: Factorize $A = L\cdot U \implies LU \vec{x} = \vec{b}$.
> - Phase 2: Solve $L\vec{y} = \vec{b}$ for $\vec{y}$ by forward substitution.
> - Phase 3: Solve $U\vec{x} = \vec{y}$ for $\vec{x}$ by backward substitution.

#### §3.2.2 Computational cost of Gaussian elimination

> Algorithm 3.3 Phase 1: Factorize $A = L\cdot U$
> ```python
> L = diag(1) # indentity matrix
> U = A
> for p = 1 : n-1
>     for r = p+1 : n
>         t = - U(r, p) / U(p, p)
>         U(r, p) = 0
>         for c = p+1 : n
>             U(r, c) = U(r, c) + t * U(p, c)
>         end for
>         L(r, p) = -t
>     end for
> end for
> ```
> - `p`: row number of the pivot element
> - `r`: number of the current row
> - `c`: number of the current column

<u>Computational cost</u>:
We count the number of:
- A: addition or subtraction;
- M: multiplication or division.

Useful identities:
- $\sum_{k=1}^{n-1} 1 = n - 1$
- $\sum_{k=1}^{n-1} k = \frac{1}{2}n(n-1)$
- $\sum_{k=1}^{n-1} k^2 = \frac{1}{6}n(n-1)(2n-1)$

- Count A:
```python
for p = 1 : n-1
    for r = p+1 : n
        for c = p+1 : n
            U(r, c) = U(r, c) + t * U(p, c)
```

$$\begin{aligned}A &= \sum_{p=1}^{n-1} \sum_{r=p+1}^{n} \sum_{c=p+1}^{n} 1\cdot A\\ &= \sum_{p=1}^{n-1} \sum_{r=p+1}^{n} (n - p) \cdot A\\&= \sum_{p=1}^{n-1} (n - p)^2 \cdot A \\ &= \sum_{p=1}^{n-1} \left(n^2 - 2np + p^2\right) \cdot A \\ &= \sum_{p=1}^{n-1} n^2 \cdot A - 2n \sum_{p=1}^{n-1} p \cdot A + \sum_{p=1}^{n-1} p^2 \cdot A \\ &= n^2 \sum_{p=1}^{n-1} 1 \cdot A - 2n \sum_{p=1}^{n-1} p \cdot A + \sum_{p=1}^{n-1} p^2 \cdot A \\ &= n^2 (n - 1) \cdot A - 2n \frac{1}{2}n(n-1) \cdot A + \frac{1}{6}n(n-1)(2n-1) \cdot A \\ &= \left[n^2(n-1) - 2n\frac{1}{2}n(n-1) + \frac{1}{6}n(n-1)(2n-1)\right] \cdot A \\ &= \left(\frac{1}{3}n^3 - \frac{1}{2}n^2 + \frac{1}{6}n\right) \cdot A \\ &= \left(\frac{1}{3}n^3 + \mathcal{O}\left(n^2\right)\right) \cdot A \end{aligned}$$

- Count M: left for exercise.
$M = \left(\frac{1}{3}n^3 + \mathcal{O}\left(n^2\right)\right) \cdot M$

- Total number of floating point operations (flops)
$W = \frac{2}{3}n^3 + \mathcal{O}\left(n^2\right)$ flops.

Algorithm 3.3 Phase 2: $L \vec{y} = \vec{b}$ (forward substitution)
```python
y = b
for r = 2 : n
    for c = 1 : r-1
        y(r) = y(r) - L(r, c) * y(c)
    end for
end for
```
<u>Computational cost</u>:
- Total number of floating point operations (flops)
$W = n^2 + \mathcal{O}\left(n\right)$ flops.

Algorithm 3.3 Phase 3: $U \vec{x} = \vec{y}$ (backward substitution)
```python
x = y
for r = n : -1 : 1
    for c = r-1 : -1 : 1
        x(r) = x(r) - U(r, c) * x(c)
    end for
    x(r) = x(r) / U(r, r)
end for
```

<u>Computational cost</u>:
- Total number of floating point operations (flops)
$W = n^2 + \mathcal{O}\left(n\right)$ flops.

<u>Total computational cost of Algorithm 3.3 (Gaussian elimination)</u>:
$$\begin{aligned}W &= \frac{2}{3}n^3 + \mathcal{O}\left(n^2\right) + 2 \cdot \left( n^2 + \mathcal{O}\left(n\right) \right) \text{ flops} \\ & = \frac{2}{3}n^3 + \mathcal{O}\left(n^2\right) \text{ flops} \end{aligned}$$

### §3.2.3 Pivoting

<u> Observation</u>: LU-factorization breaks down at step `i` if the pivot element $a_{ii}^{(i)} = 0$.

> Example 3.3: Consider the system $A\vec{x} = \vec{b}$ with $$\begin{pmatrix}0 & 1 \\ 2 & 1\end{pmatrix} \begin{pmatrix}x_1 \\ x_2\end{pmatrix} = \begin{pmatrix}1 \\ 3\end{pmatrix}$$
> - LU breaks at first step since $a_{11}^{(1)} = 0$.
> - We thus swap 1st and 2nd rows first
> $$\begin{pmatrix} 2 & 1 \\ 0 & 1\end{pmatrix} \begin{pmatrix} x_1 \\ x_2\end{pmatrix} = \begin{pmatrix} 3 \\ 1\end{pmatrix}$$
> and find the LU-factorization of the swapped matrix instead.
> - In matrix form:
> $\mathcal{P}A\vec{x} = \mathcal{P}\vec{b}$, with $\mathcal{P} = \begin{pmatrix}0 & 1 \\ 1 & 0\end{pmatrix}$.

> Definition 3.4: $\mathcal{P}$ is a permutation matrix if and only if $\mathcal{P}$ is obtained from the identity matrix by swapping any number of rows.

> Theorem 3.2: For all $A \in \mathbb{R}^{n \times n}$, there exists a permutation matrix $\mathcal{P}$, a unit lower triangular matrix $L$, and an upper triangular matrix $U$ such that $\mathcal{P}A = LU$.

> Corollary 3.3: If $A$ is nonsingular, then $A\vec{x} = \vec{b}$ can be solved by LU-factorization applied to $\mathcal{P}A$.

### §3.3 Conditioning of $A\vec{x} = \vec{b}$ and stability of Gaussian Elimination

- Recall: condition is a property of the exact mathematical problem, while stability is a property of the algorithm used to solve the problem.

- Condition of $A\vec{x} = \vec{b}$:
$$\vec{x} = f_P(A, \vec{b}) = A^{-1}\vec{b}$$
  - Absolute condition number: $\kappa_A = \left\|\Delta \vec{x}\right\| / \left\|\Delta(A, \vec{b})\right\|$
  - Relative condition number: $\kappa_R = \frac{\left\|\Delta \vec{x}\right\| / \left\|\vec{x}\right\|}{\left\|\Delta(A, \vec{b})\right\| / \left\|(A, \vec{b})\right\|}$

- Stability of Gaussian elimination (LU-factorization) can be still unable even $A\vec{x} = \vec{b}$ is well-conditioned.

#### §3.3.1 Matrix norm
> Definition 3.5: The natual matrix $p$-norm:
> $$\left\|A\right\|_p = \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|A\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p}$$
> - "Natural", since its induced by the vector $p$-norm.

> Theorem 3.4: $\left\|A\right\|_p$ is a norm, i.e., for all $A\in \mathbb{R}^{m \times n}$:
> 1. $\left\|A\right\|_p \geq 0$, $\left\|A\right\|_p = 0$ iff $A = 0$;
> 2. $\left\|\alpha A\right\|_p = \left|\alpha\right| \left\|A\right\|_p$ for all $\alpha \in \mathbb{R}$;
> 3. $\left\|A + B\right\|_p \leq \left\|A\right\|_p + \left\|B\right\|_p$ for all $A, B \in \mathbb{R}^{m \times n}$.
>
> Proof:
> 1. $\left\|A\right\|_p := \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|A\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} \geq 0$ is trivial.
> ($\Leftarrow$) If $A = 0$, then $\left\|A\vec{x}\right\|_p = \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|0\right\|_p}{\left\|\vec{x}\right\|_p} = 0$ since $\left\|A\vec{x}\right\|_p = 0$.
> ($\Rightarrow$) If $\left\|A\right\|_p = 0$, then $\left\|A\vec{x}\right\|_p = 0$ for all $\vec{x} \in \mathbb{R}^n$. We pick $\vec{x} = \vec{e}_i \forall i = 1, 2, \ldots, n$, then we have $\left\|\vec{A_i} \right\|_p = 0 \forall i = 1, 2, \ldots, n$, where $\vec{A_i}$ is the $i$-th column of $A$. This implies $A = 0$.
> 2. $\left\|\alpha A\right\|_p = \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|\alpha A\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} = \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left|\alpha\right| \left\|A\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p}  - |\alpha| \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|A\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} = |\alpha| \left\|A\right\|_p$.
> 3. 
> $$\begin{aligned}\left\|(A+B) \vec{x}\right\|_p &= \left\|A\vec{x} + B\vec{x}\right\|_p \\ &\leq \left\|A\vec{x}\right\|_p + \left\|B\vec{x}\right\|_p \\ &\leq \left\|A\right\|_p \left\|\vec{x}\right\|_p + \left\|B\right\|_p \left\|\vec{x}\right\|_p \\ \implies \frac{\left\|(A+B)\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} &\leq \left\|A\right\|_p + \left\|B\right\|_p \forall \vec{x} \neq 0 \\ \implies \left\|A+B\right\|_p & = \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|(A+B)\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} \leq \left\|A\right\|_p + \left\|B\right\|_p \end{aligned}$$

> Proposition 3.3:
> 1. $\left\| A\vec{x} \right\|_p \leq \left\|A\right\|_p \left\|\vec{x}\right\|_p$
> 2. $\left\|AB\right\|_p \leq \left\|A\right\|_p \left\|B\right\|_p$
>
> Proof:
> 1. For any $\tilde{x} \in \mathbb{R}^n$, we have $$\begin{aligned} \left\|A\right\|_p &= \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|A\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} \geq \frac{\left\|A\tilde{x}\right\|_p}{\left\|\tilde{x}\right\|_p} \\ \text{i.e. } \left\|A\right\|_p &\geq \frac{\left\|A\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} \\ \iff \left\|A\tilde{x}\right\|_p &\leq \left\|A\right\|_p \left\|\tilde{x}\right\|_p \forall \tilde{x} \in \mathbb{R}^n \end{aligned}$$
> 2. For any $\vec{x} \in \mathbb{R}^n$, we have $$\begin{aligned} \left\|AB\vec{x}\right\|_p = \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|AB\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} &\leq \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|A\right\|_p \left\|B\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} \\ &\leq \left\|A\right\|_p \sup_{\left\|\vec{x}\right\|_p \neq 0} \frac{\left\|B\vec{x}\right\|_p}{\left\|\vec{x}\right\|_p} \\ &\leq \left\|A\right\|_p \left\|B\right\|_p \end{aligned}$$

> Proposition 3.4:
> 1. $\left\|A\right\|_1 = \max_{1 \leq j \leq n} \sum_{i=1}^n \left|a_{ij}\right|$
>> Ex: A = $\begin{pmatrix}1 & 2 & 3 \\ -4 & -5 & -6 \\ 7 & 8 & 9\end{pmatrix}$
>> $\left\|A\right\|_1 = \max\left\{|1| + |-4| + |7|, |2| + |-5| + |8|, |3| + |-6| + |9|\right\} = \max\left\{12, 15, 18\right\} = 18$
> 2. $\left\|A\right\|_{\infty} = \max_{1 \leq i \leq n} \sum_{j=1}^n \left|a_{ij}\right|$
>> Ex: A = $\begin{pmatrix}1 & 2 & 3 \\ -4 & -5 & -6 \\ 7 & 8 & 9\end{pmatrix}$
>> $\left\|A\right\|_{\infty} = \max\left\{|1| + |2| + |3|, |-4| + |-5| + |-6|, |7| + |8| + |9|\right\} = \max\left\{6, 15, 24\right\} = 24$
> 3. $\left\|A\right\|_2 = \max_{1 \leq i \leq n} \sigma_i$, where $\sigma_i$'s are the singular values of $A$.

___
Intermission: Singular value decomposition (SVD)

- Given $A \in \mathbb{R}^{n \times n}$, a singular value decomposition (SVD) of $A$ is a factorization of the form:
  $$A = U\Sigma V^T$$
  where:
  - $U \in \mathbb{R}^{n \times n}$ is orthogonal (i.e., $U^TU = I$);
  - $V \in \mathbb{R}^{n \times n}$ is orthogonal;
  - $\Sigma \in \mathbb{R}^{n \times n}$ is diagonal
    $\Sigma = \begin{pmatrix}\sigma_1 & 0 & \cdots & 0 \\ 0 & \sigma_2 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \sigma_n\end{pmatrix}$
    The singular values $\{\sigma_j\}_{j=1}^n$ are uniquely determined.

- Geometric interpretation of SVD:
  $A = U\Sigma V^T \implies AV = U\Sigma(V^TV) = U\Sigma$
    - ![image3.3](images/image3.3.png)
    - $A \vec{v}_j = \sigma_j \vec{u}_j \forall i\leq j \leq n$
    - $\{\vec{u}_j\}_{j=1}^n$ is a set of unit orthogonal basis vectors of $\mathbb{R}^n$.
    - $\{\vec{v}_j\}_{j=1}^n$ is a set of unit orthogonal basis vectors of $\mathbb{R}^n$.
    - $A$ maps a given $\vec{v}_j$ to $\sigma_j \vec{u}_j$ and stretches it by $\sigma_j$.

> Example:
> A = $\begin{pmatrix}1.25 & 1.75 \\ 1.75 & 1.25\end{pmatrix}$
> SVD of A:
> ![image3.4](images/image3.4.png)

In general, SVD describes that:
> The image of the unit sphere under $A \in \mathbb{R}^{n \times n}$ is a hyperellipse with lengths of principal semiaxes determined by the singular values $\{\sigma_j\}_{j=1}^n$.

- Relation to eigenvalue decomposition:
  - SVD applies to any $m \times n$ matrix; eigen-decomposition applies to only square diagonalizable matrices.
  - $$\begin{aligned}A = U\Sigma V^T &\implies AA^T = U\Sigma V^T(V\Sigma^T U^T) = U(\Sigma\Sigma^T)U^T \\ \implies AA^T &\text{ has eigenvalues } \{\sigma_j^2\}_{j=1}^n \end{aligned}$$
  ___

#### §3.3.2 Conditioning of $A\vec{x} = \vec{b}$

<u>Problem:</u> Find $\vec{x}$ such that $A\vec{x} = \vec{b}$, i.e. $\vec{x} = A^{-1}\vec{b} (= f(A, \vec{b}))$.

To find the condition number of $f_\mathcal{P}(A, \vec{b})$, we consider perturbations on $A$ and $\vec{b}$ and the resulting perturbation of $\vec{x}$.

The relative condition number:
$$\kappa_R = \frac{\left\|\Delta \vec{x}\right\| / \left\|\vec{x}\right\|}{\left\|\Delta(A, \vec{b})\right\| / \left\|(A, \vec{b})\right\|}$$

To simplify the derivation,

Case 1 (No perturbation in A)

$$
\begin{aligned}
\vec{x} + \Delta \vec{x} &= A^{-1} (\vec{b} + \Delta \vec{b}) \\
\vec{x} + \Delta \vec{x} &= A^{-1} \vec{b} + A^{-1} \Delta \vec{b} \\
\Delta \vec{x} &= A^{-1} \Delta \vec{b} \\
\|\Delta \vec{x}\| &= \|A^{-1} \Delta \vec{b}\| \leq \|A^{-1}\| \|\Delta \vec{b}\| \quad \ldots (1)
\end{aligned}
$$

Moreover:

$$
\begin{aligned}
A \vec{x} &= \vec{b} \\
\|A \vec{x}\| &= \|\vec{b}\| \\
\|A\| \|\vec{x}\| &\geq \|\vec{b}\| \\
\|\vec{x}\| &\geq \|A\|^{-1} \|\vec{b}\| \\
\|\vec{x}\|^{-1} &\leq \|A\| \|\vec{b}\|^{-1} \quad \ldots (2)
\end{aligned}
$$

Combining (1) & (2):

$$
\frac{\|\Delta \vec{x}\|}{\|\vec{x}\|} \leq \|A^{-1}\| \|\Delta \vec{b}\| \|A\| \|\vec{b}\|^{-1} = \|A^{-1}\| \|A\| \frac{\|\Delta \vec{b}\|}{\|\vec{b}\|}
$$

$$
\kappa_R = \frac{\|\Delta \vec{x}\| / \|\vec{x}\|}{\|\Delta \vec{b}\| / \|\vec{b}\|} \leq \|A^{-1}\| \|A\|
$$

*See [LT] Lecture 12 for special cases where equality is attained*

Thus, $\|A^{-1}\|\|A\|$ is an upper bound for $\kappa_R$.

Case 2 (No perturbation in $\vec{b}$)

$$
\kappa_R = \frac{\|\Delta \vec{x}\| / \|\vec{x}\|}{\|\Delta A\| / \|A\|}
$$

$$
\begin{aligned}
\vec{x} + \Delta \vec{x} &= (A + \Delta A)^{-1} \vec{b} \\
(A + \Delta A)(\vec{x} + \Delta \vec{x}) &= \vec{b} \\
A \vec{x} + A \Delta \vec{x} + \Delta A \vec{x} + \Delta A \Delta \vec{x} &= \vec{b} \\
A \Delta \vec{x} + \Delta A \vec{x} + \Delta A \Delta \vec{x} &= 0 \\
A \Delta \vec{x} &= -\Delta A (\vec{x} + \Delta \vec{x}) \\
\|\Delta \vec{x}\| &\leq \|A^{-1}\Delta A (\vec{x} + \Delta \vec{x})\| \leq \|A^{-1}\| \|\Delta A\| \|\vec{x} + \Delta \vec{x}\|
\end{aligned}
$$

$$
\Rightarrow \frac{\|\Delta \vec{x}\|}{\|\vec{x} + \Delta \vec{x}\|} \leq \|A^{-1}\| \|\Delta A\|
$$

Consider $\Delta \vec{x}$ to be small compared to $\vec{x}$ and thus $\|\vec{x} + \Delta \vec{x}\| \approx \|\vec{x}\|$, which yields,

$$
\frac{\|\Delta \vec{x}\|}{\|\vec{x}\|} \leq \|A^{-1}\| \|\Delta A\|
$$

$$
\kappa_R = \frac{\|\Delta \vec{x}\| / \|\vec{x}\|}{\|\Delta A\| / \|A\|} \leq \|A^{-1}\| \|A\|
$$

*See [LT] Lecture 12 for special cases where equality is attained*

Case 3 (Perturbations in both $A$ and $\vec{b}$)

$$
\kappa_R \leq \|A^{-1}\| \|A\| \quad (\text{derivations omitted})
$$

*Since $\|A\|\|A^{-1}\|$ pops up a lot, especially for conditioning of linear algebra, we give it a name.*

> Definition 3.6 (Condition number of a matrix A):
> The condition number of a matrix \( A \) (invertible) is
> 
> $$
> \kappa(A) = \|A\| \|A^{-1}\|
> $$
> 
> $$
> \kappa_p(A) = \|A\|_p \|A^{-1}\|_p
> $$

>Proposition 3.5:
> For \( A \in \mathbb{R}^{n \times n} \), \(\det(A) \neq 0\), and \( p = 2 \),
> 
> $$
> \kappa_2(A) = \|A\|_2 \|A^{-1}\|_2 = \frac{\sigma_{\text{max}}(A)}{\sigma_{\text{min}}(A)}
> $$
> Proof:
> SVD of \( A = U \Sigma V^T \)
> 
> then \( A^{-1} = (V^T)^{-1} \Sigma^{-1} U^{-1} \)
> 
> $$
> A^{-1} = V \begin{pmatrix}
> \sigma_1^{-1} & 0 \\
> 0 & \sigma_n^{-1}
> \end{pmatrix} U^T
> $$
> 
> Therefore: singular values of \( A^{-1} \): \( \{\sigma_i^{-1}\} \)
> 
> Due to \( \|A\|_2 = \max_i \sigma_i \) (Proposition 3.4(iii))
> 
> $$
> \|A\|_2 \|A^{-1}\|_2 = \max_i \sigma_i \cdot \max_i (\sigma_i^{-1}) = \frac{\sigma_{\text{max}}(A)}{\sigma_{\text{min}}(A)}
> $$

> Back to Intuition (SVD) Example:
> $$
> A = \begin{pmatrix}
> 1.25 & 1.75 \\
> 1.75 & 1.25
> \end{pmatrix}
> \quad \text{SVD} \quad \sigma_1 = 3, \sigma_2 = 0.5
> $$
> 
> $$
> \kappa_2(A) = \|A\|_2 \|A^{-1}\|_2 = \frac{\sigma_1}{\sigma_2} = 6
> $$
> 
> - Matrix \( A \) is well-conditioned
> - Linear algebra problems like solving \( A \vec{x} = \vec{b} \) or matrix-vector multiplication \( A \vec{x} \) are well-conditioned

> Example 3.4 (An ill-conditioned matrix):
> $$
> A = \begin{pmatrix}
> 1000 & 2000 \\
> 499 & 1001
> \end{pmatrix}
> \implies \sigma_1 \approx 2500, \sigma_2 \approx 1.2
> $$
> 
> $$
> \kappa_2(A) = \frac{\sigma_1}{\sigma_2} \approx 2084 \implies \text{ill-conditioned}
> $$
> 
> Let's illustrate the ill-conditionedness:
> 
> Consider the linear system:
> 
> $$
> A \vec{x} = \vec{b}
> $$
> 
> $$
> \begin{pmatrix}
> 1000 & 2000 \\
> 499 & 1001
> \end{pmatrix}
> \vec{x} = \begin{pmatrix}
> 3000 \\
> 1500
> \end{pmatrix}
> \quad \text{exact solution} \quad \vec{x} = \begin{pmatrix}
> 1 \\
> 1
> \end{pmatrix}
> $$
> 
> If we perturb \( A \) a bit, (the last element \( 1001 \rightarrow 1000 \)):
> 
> $$
> \begin{pmatrix}
> 1000 & 2000 \\
> 499 & 1000
> \end{pmatrix}
> \vec{x} = \begin{pmatrix}
> 3000 \\
> 1500
> \end{pmatrix}
> \quad \text{exact solution} \quad \vec{x} = \begin{pmatrix}
> 0 \\
> 1.5
> \end{pmatrix}
> $$

*Small change in \( A \) \(\Rightarrow\) large change in \( \vec{x} \)*

*In practice, when having ill-conditioned matrices, before solving for anything, we need to perform **preconditioning** to bring the condition number down! (Beyond the scope of this course)*

### §3.3.3 Stability of the Gaussian Elimination algorithm

*Gaussian Elimination can be unstable even if \( A \) (Alg 3.3) is well-conditioned.*

> Example 3.5 (Unstable G.E. on a well-conditioned matrix):
> 
> $$
> A = \begin{pmatrix}
> \epsilon & 1 \\
> 1 & 1
> \end{pmatrix}
> \implies \begin{pmatrix}
> \epsilon & 1 \\
> 1 & 1
> \end{pmatrix}
> \begin{pmatrix}
> x_1 \\
> x_2
> \end{pmatrix}
> = \begin{pmatrix}
> 1 \\
> 2
> \end{pmatrix}
> \quad \text{with} \quad 0 < \epsilon \ll 1
> $$
> 
> - \( A \) is well-conditioned
> 
> To find singular values of \( A \), we (equivalently) find the eigenvalues of \( A^T A \):
> 
> $$
> A = U \Sigma V^T \implies A^T A = V \Sigma^T U^T U \Sigma V^T = V \Sigma^T \Sigma V^T
> $$
> 
> - \( V \Sigma^T \Sigma V^T \) is the eigen-decomposition of \( A^T A \)
> - eigenvalues of \( A^T A \) are squares of singular values of \( A \)
> 
> $$
> A^T A = \begin{pmatrix}
> \epsilon & 1 \\
> 1 & 1
> \end{pmatrix}
> \begin{pmatrix}
> \epsilon & 1 \\
> 1 & 1
> \end{pmatrix}
> = \begin{pmatrix}
> \epsilon^2 + 1 & \epsilon + 1 \\
> \epsilon + 1 & 2
> \end{pmatrix}
> $$
> 
> $$
> \lambda (A^T A) = \frac{1}{2} \left[ (\epsilon^2 + 3) \pm \sqrt{(\epsilon^2 + 3)^2 - 4(\epsilon^2 - 1)^2} \right]
> $$
>
> $$
> \epsilon^2 + 3 \pm \sqrt{(\epsilon^2 + 3)^2 - 4(\epsilon - 1)^2}
> $$
> 
> Observations:
> - \( \sigma_1 = (\epsilon^2 + 3) \pm \sqrt{(\epsilon^2 + 3)^2 - 4(\epsilon - 1)^2} \geq 0 \quad \forall \epsilon \in \mathbb{R} \)
> - \( \sigma_2 \approx 5 \) when \( 0 < \epsilon \ll 1 \)
> 
> Thus: \( \kappa_2(A) = \frac{\sigma_{\text{max}}}{\sigma_{\text{min}}} = \sqrt{\frac{\lambda_{\text{max}} (A^T A)}{\lambda_{\text{min}} (A^T A)}} \approx \frac{\sqrt{3 + \epsilon}}{\sqrt{3 - \epsilon}} \approx 2.618 \)
> 
> Therefore \( A \) is well-conditioned when \( 0 < \epsilon \ll 1 \).
> 
> - Gaussian Elimination is unstable on \( A \)
> 
> Phase 1:
> 
> $$
> LU \text{ on } A \implies L = \begin{pmatrix}
> 1 & 0 \\
> \frac{1}{\epsilon} & 1
> \end{pmatrix}, \quad U = \begin{pmatrix}
> \epsilon & 1 \\
> 0 & 1 - \frac{1}{\epsilon}
> \end{pmatrix}
> $$
> 
> Phase 2:
> 
> $$
> Ly = \vec{b}
> $$
> 
> $$
> \begin{pmatrix}
> 1 & 0 \\
> \frac{1}{\epsilon} & 1
> \end{pmatrix}
> \begin{pmatrix}
> y_1 \\
> y_2
> \end{pmatrix}
> = \begin{pmatrix}
> 1 \\
> 2
> \end{pmatrix}
> \implies y_1 = 1, \quad y_2 = 2 - \frac{1}{\epsilon}
> $$
> Phase 3:
> 
> $$
> U \vec{x} = \vec{y}
> $$
> 
> $$
> \begin{pmatrix}
> \epsilon & 1 \\
> 0 & 1 - \frac{1}{\epsilon}
> \end{pmatrix}
> \begin{pmatrix}
> x_1 \\
> x_2
> \end{pmatrix}
> = \begin{pmatrix}
> 1 \\
> 2 - \frac{1}{\epsilon}
> \end{pmatrix}
> \implies
> \begin{cases}
> x_1 = \frac{1}{\epsilon} \left( 1 - \frac{2 - \frac{1}{\epsilon}}{1 - \frac{1}{\epsilon}} \right) \\
> x_2 = \frac{2 - \frac{1}{\epsilon}}{1 - \frac{1}{\epsilon}}
> \end{cases}
> $$
> 
> - Consider using the floating point system \( F \langle 10, 4, -3, 5 \rangle \) with \( \epsilon = 10^{-5} \). We would have:
> 
> $$
> \hat{x}_2 = \text{fl} \left( \frac{\text{fl}(2 - \epsilon)}{\text{fl}(1 - \epsilon)} \right)
> $$
> 
> $$
> = \text{fl} \left( \frac{\text{fl}(2 - 10^{-5})}{\text{fl}(1 - 10^{-5})} \right)
> $$
> 
> $$
> = \text{fl} \left( \frac{\text{fl}(1.99998)}{\text{fl}(0.99999)} \right)
> $$
> 
> $$
> = \text{fl} \left( \frac{-1.000 \times 10^5}{-1.000 \times 10^5} \right) = \text{fl}(1) = 1
> $$
> 
> Exact \( x_2 = \frac{2 - \epsilon}{1 - \epsilon} \approx 1 \) when \( 0 < \epsilon \ll 1 \).
> 
> $$
> \hat{x}_1 = \text{fl} \left( \frac{1}{\epsilon} \left( 1 - \hat{x}_2 \right) \right) = 0
> $$
> 
> Exact \( x_1 = \frac{1}{\epsilon} \left( 1 - \frac{2 - \epsilon}{1 - \epsilon} \right) = \frac{1}{\epsilon} \left( \frac{1 - \epsilon}{1 - \epsilon} - \frac{2 - \epsilon}{1 - \epsilon} \right) = 1 \).
> 
> Hence, large error in \( \hat{x}_1 \).

> Algorithm 3.4 (Gaussian elimination with partial pivoting):
> Essentially, during each step of LU factorization, we rearrange the rows such that we get the largest pivoting element (in its absolute value).

> Example 3.5 (With partial pivoting):
> 
> Partial pivoting gives:
> 
> $$
> \begin{pmatrix}
> 1 & 1 \\
> \epsilon & 1
> \end{pmatrix}
> \begin{pmatrix}
> x_1 \\
> x_2
> \end{pmatrix}
> = \begin{pmatrix}
> 2 \\
> 1
> \end{pmatrix}
> $$
> Phase 1:
> $$
> L = \begin{pmatrix}
> 1 & 0 \\
> \epsilon & 1
> \end{pmatrix}, \quad U = \begin{pmatrix}
> 1 & 1 \\
> 0 & 1 - \epsilon
> \end{pmatrix}
> $$
> Phase 2:
> $$
> \begin{pmatrix}
> 1 & 0 \\
> \epsilon & 1
> \end{pmatrix}
> \begin{pmatrix}
> y_1 \\
> y_2
> \end{pmatrix}
> = \begin{pmatrix}
> 2 \\
> 1
> \end{pmatrix}
> \implies
> \begin{cases}
> y_1 = 2 \\
> y_2 = 1 - 2 \epsilon
> \end{cases}
> $$
> Phase 3:
> $$
> \begin{pmatrix}
> 1 & 1 \\
> 0 & 1 - \epsilon
> \end{pmatrix}
> \begin{pmatrix}
> x_1 \\
> x_2
> \end{pmatrix}
> = \begin{pmatrix}
> 2 \\
> 1 - 2 \epsilon
> \end{pmatrix}
> \implies
> \begin{cases}
> x_1 = 2 - \frac{1 - 2 \epsilon}{1 - \epsilon} \\
> x_2 = \frac{1 - 2 \epsilon}{1 - \epsilon}
> \end{cases}
> $$
> 
> With \( F [b=10, m=4, e = -5 ] \) and \( \epsilon = 10^{-5} \):
> 
> $$
> \hat{x}_2 = \text{fl} \left( \frac{1 - 2 \epsilon}{1 - \epsilon} \right) = \text{fl} \left( \frac{\text{fl}(0.99998)}{\text{fl}(0.99999)} \right) = \frac{0.1000 \times 10^1}{0.1000 \times 10^1} = 1
> $$
> 
> $$
> \hat{x}_1 = \text{fl}(2-\hat{x_2}) = 1$$
