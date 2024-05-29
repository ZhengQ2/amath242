## Chapter 2: Root finding algorithms
### §2.1 Introduction
> <u>Problem</u>: Given any function $f(x)$, find $x^*$ such that $f(x^*) = 0$. Then $x^*$ is called a root of $f(x)=0$.
> - Computationally, we can only aim for finding $x^*$ such that $|f(x^*)| < \tau$ for some error tolerance $\tau$.
> - There is no method that guarantees finding the root(s) of any function.

> Theorem 2.1: (Intermediate Value Theorem)
> If $f(x)$ is continuous on a closed interval $[a, b]$ and $c \in [f(a), f(b)]$, then there exists $x^* \in [a, b]$ such that $f(x^*) = c$.

> Corollary 2.2: If $f(a)f(b) < 0$ for a closed interval $[a, b]$, then there exists at least one root $x^*$ as long as $f(x)$ is continuous$.

### §2.2 Four algorithms for root finding
#### §2.2.1 Bisection method
> Theorem 2.3: If $f(x)$ is continuous on $[a_0, b_0]$ such that $f(a_0)f(b_0) \leq 0$, then $[a_k, b_k]$, defined by $[a_k, b_k] := \begin{cases} [a_{k-1}, \frac{a_{k-1}+b_{k-1}}{2}] & \text{if } f(a_{k-1})f(\frac{a_{k-1}+b_{k-1}}{2}) \leq 0 \\ [\frac{a_{k-1}+b_{k-1}}{2}, b_{k-1}] & \text{if } f(a_{k-1})f(\frac{a_{k-1}+b_{k-1}}{2}) > 0 \end{cases}$ holds that $f(a_k)f(b_k) \leq 0$ for any $k \geq 1$.
![image2.1](images/image2.1.png)

> Algorithm 2.1: Bisection algorithm
> ```python
> # Input: f(x), a, b, tolerance tau.
> # Output: x, an approximant of x^* (f(x^*) = 0).
> while |b-a| > tau: # alternatively, |f((a+b)/x)| > tau
>     c = (a+b)/2
>     if f(a)*f(c) <= 0:
>         b = c
>     else:
>         a = c
>     end if
> end while
> x = (a+b)/2
> ```

#### §2.2.2 Fixed-point iteration
> Definition 2.1: We say that $x^*$ is a fixed point of $g(x)$ if $g(x^*) = x^*$.

> <u>Problem</u>: Root finding of f(x) is equivalent to finding the fixed point of $g(x) = x - f(x)$.
> $f(x^*) = 0 \Leftrightarrow x^* - f(x^*) = x^* \Leftrightarrow g(x^*) = x^*$.

- The fixed point interation: $x_{n+1} = g(x_n), n = 0, 1, 2, \cdots$.

> Algorithm 2.2: Fixed-point iteration
> ```python
> # Input: g(x), x_0, tolerance tau.
> # Output: x, an approximant of x^* (f(x^*) = 0).
> i = 0
> repeat
>    i = i + 1
>    x[i] = g(x[i-1])
> until |x[i] - x[i-1]| < tau
> x = x[i]
> ```

Convergence? If $|g'(x)| < 1$ and $x_0$ is close enough to $x^*$, then algorithm 2.2 will converge to $x^*$. More on this in §2.4.

#### §2.2.3 Newton's method
Consider the Taylor series of $f(x^*)$ about an initial estimate $x_0$:
$$\begin{aligned} 0 = f(x^*) &= f(x_0) + f'(x_0)(x^*-x_0) + \mathcal{O}((x^*-x_0)^2) \\ f(x^*) &\approx f(x_0) + f'(x_0)(x^*-x_0) \\ 0 &= f(x_0) + f'(x_0)(x_1-x_0) \\ x_1 &= x_0 - \frac{f(x_0)}{f'(x_0)} \end{aligned}$$
![image2.2](images/image2.2.png)
Thus we have iteration:
$ 0 = f(x_k) \approx f(x_{k-1}) + f'(x_{k-1})(x_k-x_{k-1}) \Rightarrow x_k = x_{k-1} - \frac{f(x_{k-1})}{f'(x_{k-1})}$.

> Algorithm 2.3: Newton's method
> ```python
> # Input: f(x), f'(x), x_0, tolerance tau.
> # Output: x, an approximant of x^* (f(x^*) = 0).
> i = 0
> x[0] = x_0
> repeat
>    i = i + 1
>    if f'(x[i-1]) == 0 stop
>    x[i] = x[i-1] - f(x[i-1])/f'(x[i-1])
> until |x[i] - x[i-1]| < tau # or |f(x[i])| < tau
> x = x[i]
> ```
Convergence? Under some conditions on $f(x)$, $f'(x)$ and that $x_0$ is close enough to $x^*$, Newton's method converges and converges "quadratically"! More on this in §2.4.