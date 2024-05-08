# AMATH242/CS371: Introduction to Computational Mathematics

## Chapter 1: Floating Point Systems
- Numarical algorithms and computers operate on <u>finite percision arithmetic</u>.
- We do not have the totality of $\mathbb{R}$. In fact, we only have a "tiny" portion of $\mathbb{R}$.

> Definition 1.1: Let $\hat{x}$ be an approximation of a real number $x$.
> - Absolute error: $\Delta x = x - \hat{x}$
> - Relative error: $\delta x = \frac{x - \hat{x}}{x}$

### §1.1 Floating Point System
#### §1.1.1 Intro
> Definition 1.2: A floating point system $F \subseteq \mathbb{R}$ is a subset of the real numbers whose elements have the form $z = \pm (0.x_1x_2\dots x_m)_b \times b^{\pm(y_1y_2\dots y_e)_b}$, where $b$ is the base of the system, $x_1x_2\dots x_m$ is the mantissa, and $y_1y_2\dots y_e$ is the exponent
> which is categorized by:
> - Base $b_f$
> - Mantissa $m_f$
> - Exponent $e_f$
> noted as $F[b = b_f, m = m_f, e = e_f]$
> and $0 \leq x_i < b-1$, $0 \leq y_i < b-1$, and $1 \leq i \leq m$ and $1 \leq j \leq e$.

> Example 1.1: $F[b = 10, m = 3, e = 2]$
> - $z = \pm (0.x_1x_2x_3) \times 10^{\pm(y_1y_2)}$
> An example: $z = 0.127 \times 10^{19}$, in which $x_1 = 1, x_2 = 2, x_3 = 7, y_1 = 1, y_2 = 9$.

> Definition 1.3 (Normalization): A floating point number $z \in F \subseteq \mathbb{R}$ $\left(z = \pm (0.x_1x_2\dots x_m)_b \times b^{\pm(y_1y_2\dots y_e)_b}\right)$ is normalized if $x_1 \geq 1$.

> Example 1.1 (cont'd): $F[b = 10, m = 3, e = 2]$
> - $z_1 = 0.127 \times 10^{19}$ is normalized.
> - $z_2 = 0.034 \times 10^{-5}$ is not normalized, but can be normalized to $0.340 \times 10^{-6}$.

___
<u>Intermission</u>: (Converting binary to decimal)

$$\begin{align*}&(a_n a_{n-1} \dots a_1 a_0.a_{-1} a_{-2} \dots a_{-m})_2\\=& a_n \times 2^n + \dots + a_1 \times 2^1 + a_0 \times 2^0 + a_{-1} \times 2^{-1} + \dots + a_{-m} \times 2^{-m}\end{align*}$$
> Example:
> - $(101)_2 = 1 \times 2^2 + 0 \times 2^1 + 1 \times 2^0 = 5$.
> - $(1.11)_2 = 1 \times 2^0 + 1 \times 2^{-1} + 1 \times 2^{-2} = 1.75$.

We can further generalize this:

$$\begin{align*}&(a_n a_{n-1} \dots a_1 a_0.a_{-1} a_{-2} \dots a_{-m})_b\\=& a_n \times b^n + \dots + a_1 \times b^1 + a_0 \times b^0 + a_{-1} \times b^{-1} + \dots + a_{-m} \times b^{-m}\end{align*}$$
___

> Example 1.2: $F[b = 2, m = 2, e = 2]$
> Let's consider all normalized positive numbers in $F$.
> - $z = (0.x_1x_2)_2 \times 2^{\pm(y_1y_2)_2}$
> - $z_1 = (0.11)_2 \times 2^{(11)_2} = 0.75 \times 2^3 = 6$
> - $z_2 = (0.10)_2 \times 2^{(11)_2} = 0.5 \times 2^3 = 4$
> - $z_3 = (0.11)_2 \times 2^{(10)_2} = 0.75 \times 2^2 = 3$
> - $z_4 = (0.10)_2 \times 2^{(10)_2} = 0.5 \times 2^2 = 2$
> - $z_5 = (0.11)_2 \times 2^{(01)_2} = 0.75 \times 2^1 = 1.5$
> - $z_6 = (0.10)_2 \times 2^{(01)_2} = 0.5 \times 2^1 = 1$
> - $z_7 = (0.11)_2 \times 2^{(00)_2} = 0.75 \times 2^0 = 0.75$
> - $z_8 = (0.10)_2 \times 2^{(00)_2} = 0.5 \times 2^0 = 0.5$
> - $z_9 = (0.11)_2 \times 2^{(-01)_2} = 0.75 \times 2^{-1} = 0.375$
> - $z_{10} = (0.10)_2 \times 2^{(-01)_2} = 0.5 \times 2^{-1} = 0.25$
> - $z_{11} = (0.11)_2 \times 2^{(-10)_2} = 0.75 \times 2^{-2} = 0.1875$
> - $z_{12} = (0.10)_2 \times 2^{(-10)_2} = 0.5 \times 2^{-2} = 0.125$
> - $z_{13} = (0.11)_2 \times 2^{(-11)_2} = 0.75 \times 2^{-3} = 0.09375$
> - $z_{14} = (0.10)_2 \times 2^{(-11)_2} = 0.5 \times 2^{-3} = 0.0625$
> ![Image 1.1](images/image1.1.png/)
> ![Image 1.2](images/image1.2.png)

Observations:
1. Floating point numbers are not eqally spaced. The spacing "jumps" by a factor of 2 at each power of 2.
2. There is an awkward gap between 0 and the smallest normalized number.
3. 0 is unrepresentable in this system.

> Definition 1.4: The distant from 1 to the next larger normalized floating point number is called the machine epsilon, denoted as $\epsilon_{mach}$.
> We have the following,
> - $1 = (0.10\dots00)_b \times b^{(0\dots01)_b}$
> - $\text{next} = (0.10\dots01)_b \times b^{(0\dots01)_b}$
> - $\epsilon_{mach} = (0.00\dots01)_b \times b^{(0\dots01)_b} = b^{1-m}$
>
> Therefore we also found following properties:
> - number $m$ (# of digits in mantissa) is called <u>precision</u>.
> - $\epsilon_{mach}$ is also called <u>machine precision</u>.
> - $\epsilon_{mach} = b^{1-m}$.
> - IMPORTANT: the formula $\epsilon_{mach} = b^{1-m}$ is subject to slight change in single and double precision formats.

> Definition 1.5: The system $F$ can be extended by including <b>subnormal numbers</b> which are implemented by: $\pm (0.0x_2x_3\cdots x_m)_b \times b^{-(b-1,b-1,\cdots,b-1)_b}$, where $0 \leq x_2, x_3,\cdots,x_m \leq b-1$ and $(0.0x_2x_3\cdots x_m)_b \neq 0$.
> * Recall: Closest to zero normalized number: $\pm (0.10\cdots0)_b \times b^{-(b-1,b-1,\cdots,b-1)_b}$.

> Example 1.2 (cont'd): $F[b = 2, m = 2, e = 2]$
> $\pm (0.01)_2 \times 2^{-(11)_2} = \pm 0.03125$ is the only subnormal number in this system.
> ![Image 1.3](images/image1.3.png)

* If we denote the smallest normalized positive number as $\lambda$, the subnormal numbers fill the gap between 0 and $\lambda$ with the same spacing between $\lambda$ and $b\lambda$.
* Let's see another exmaple: $F[b = 2, m = 3, e = 2]$
 ![Image 1.4](images/image1.4.png)


#### §1.1.2 Rounding, overflow, and underflow
> Definition 1.6: Let $G \subseteq \mathbb{R}$ denote all real numbers that have the form $z = \pm (0.x_1x_2\dots x_m)_b \times b^y, y\in Z$, i.e., we life the lower and upper limits for the exponent.
> For $\forall x \in \mathbb{R}$, then $fl(x)$ denotes the nearest number to $x$ in $G$ and the operation $x \mapsto fl(x)$ is called <u>rounding</u>.

> Example 1.2 (cont'd): $F[b = 2, m = 2, e = 2]$
> ![Image 1.5](images/image1.5.png)
> Here, $z_1 = 8 \in G \not\in F$.
> $x_1 = 2.8 \rightarrow fl(x_1) = 3$
> $x_2 = 4.5 \rightarrow fl(x_2) = 4$
> $x_3 = 1.25 \rightarrow$ Tie!
> * Two common tie breakers:
>   1. round away from zero: $x_3 = 1.25 \rightarrow fl(x_3) = (0.11)_2 \times 2^{(01)_2} = 1.5$
>   2. round to the one with an even last digit: $x_3 = 1.25 \rightarrow fl(x_3) = (0.10)_2 \times 2^{(01)_2} = 1$
>
> $x_4 = 7.7 \rightarrow fl(x_4) = 8 \not\in F$

> Definition 1.7: We say $fl(x)$ overflows if $|fl(x)| > \max \{|z|: z \in F\}$, and underflows if $0 < |fl(x)| < \min \{|z|: z \in F\}$.

> Example 1.2 (cont'd): $F[b = 2, m = 2, e = 2]$
> $x_5 = 6.1 \rightarrow fl(x_5) = 6$.
> Thus, the following statement is false:
> "Overflow occurs when $x$ is bigger than the biggest normalized number in $F$."

> Theorem 1.1 (Unit roundoff): Each real number $x$ such that $fl(x)$ is a normalized number in $F$ has a relative error no larger than $u = \frac{1}{2}\epsilon_{mach}$, which is called unit roundoff.
> If $x\in \mathbb{R}$ such that $fl(x)$ is normalized in $F$, then $|\delta| = \left|\frac{x - fl(x)}{x}\right| \leq u = \frac{1}{2}\epsilon_{mach}$.
> * subnormal numbe correspond to bigger relative error.

#### §1.1.3 Standard floating point systems

- Single precision format (32-bit)
  ![Image 1.6](images/image1.6.png)
    - Sign bit: $s=1$ for negative, $s=0$ for positive.
    - Exponent:
        - We have $2^8 = 256$ exponents, therefore face value $[0,255]$.
        - We want a range of signed exponents. The convention is face value subtracted by a bias of 127 $\rightarrow [-127,128]$.
        - When $e = \begin{cases} (00000000)_2 \rightarrow -127 \\ (11111111)_2 \rightarrow 128 \end{cases}$ reserved for subnormals and special numbers.
        - exponent | mantissa (all zero) | mantissa (not all zero)
          |-----|-----|-----|
          | $(00\cdots 00)_2$ |  $\pm 0$ | subnormals |
          | $(11\cdots 11)_2$ | $\pm \infty$ | NaN (Not a number) |
- Double precision format (64-bit)