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