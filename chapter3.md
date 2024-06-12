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
> Case 1: If $\det(A) \neq 0$, then $\vec{x} = A^{-1}\vec{b}$ is the unique solution of $A\vec{x} = \vec{b}$.
> Case 2: If $\det(A) = 0, >    Case 2a: If $\vec{b} \in \text{range}(A)$, then $A\vec{x} = \vec{b}$ has infinitely many solutions.
>    Case 2b: If $\vec{b} \notin \text{range}(A)$, then $A\vec{x} = \vec{b}$ has no solution.

### §3.2 Gaussian elimination
#### §3.2.1 LU factorization
> Algorithm 3.1: (Gaussian elimination)
> - Phase 1: Reduce the matrix A to upper triangular form.
> - Phase 2: Solve the reduced system.

> Definition 3.2: A matrix $A \in \mathbb{R}^{n \times n}$ with components $a_{ij}$ is said to be
> - upper triangular if $a_{ij} = 0$ for $i > j$,
>> - lower triangular if $a_{ij} = 0$ for $i < j$,>
> $A \in \mathbb{R}^{n \times n}$ is said to be a triangular matrix if it is either upper or lower triangular.
.> Solving these triangular systems are relatively straightforward.
)> Algorithm 3.2: (Forward and backward substitution)$> $A\vec{x} = \vec{b}, A \in \mathbb{R}^{n \times n}, \vec{b} \in \mathbb{R}^{n\times 1}$n> - A: upper triangular $\Rightarrow$ backward substitution)> ![image3.1](images/image3.1.png)$> $$\begin{align*}x_n &= a_{n,n}^{-1} b_n \\ x_{n-1} &= a_{n-1,n-1}^{-1} (b_{n-1} - a_{n-1,n}x_n) \\ x_{n-2} &= a_{n-2,n-2}^{-1} (b_{n-2} - a_{n-2,n}x_n - a_{n-2,n-1}x_{n-1}) \\ &\vdots \\ x_i &= a_{ii}^{-1} (b_i - \sum_{k=i+1}^n a_{ik}x_k) \end{align*}$$
>
> - A: lower triangular $\Rightarrow$ forward substitution)
> ![image3.2](images/image3.2.png)
> $$\begin{align*}x_1 &= a_{11}^{-1} b_1 \\ x_2 &= a_{22}^{-1} (b_2 - a_{21}x_1) \\ x_3 &= a_{33}^{-1} (b_3 - a_{31}x_1 - a_{32}x_2) \\ &\vdots \\ x_i &= a_{ii}^{-1} (b_i - \sum_{k=1}^{i-1} a_{ik}x_k) \end{align*}$$
> Now we see how a dense matrix $A$ is reduced to an (upper-)triangular form: (with an example first)
> Example 3.2: Consider the system $A\vec{x} = \vec{b}$ with $$A = \begin{pmatrix}1 & 2 & 3 \\ \textcolor{red}{4} & 5 & 6 \\ \textcolor{red}{7} & \textcolor{blue}{8} & 1\end{pmatrix}, \vec{b} = \begin{pmatrix}1 \\ 2 \\ 3\end{pmatrix}$$
> <span style="color:red">1st</span>, <span style="color:blue">2nd</span> columns are to be eliminated.
> - Step 1: Denote $A^{(1)} = \begin{pmatrix}1 & 2 & 3 \\ \textcolor{red}{4} & 5 & 6 \\ \textcolor{red}{7} & 8 & 1\end{pmatrix}$><span style="color:red"> 1st: eliminate this sub-column $\Rightarrow A^{(2)}$</span>.
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