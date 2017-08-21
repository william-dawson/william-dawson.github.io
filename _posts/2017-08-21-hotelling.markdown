---
layout: post
title:  "Hotelling's Method"
date:   2017-08-21 11:00:00 +0900
categories: blog method
---

I would like to use this first blog post to introduce Hotelling's method
for computing the inverse of a matrix. We'll focus on the symmetric case.
First let's begin by asking why
we would want to directly compute the inverse of a matrix. After all,
one of the first things you'll learn in any numerical linear algebra class
is to always avoid explicitly inverting a matrix.

However, in quantum chemistry it's actually fairly normal to compute
a matrix inverse. In particular, when solving the generalized eigenvalue
problem:

$$\begin{equation}H\phi = \lambda S \phi.\end{equation}$$

We can reduce this generalized eigenvalue problem to the standard eigenvalue
problem using the inverse square root of the overlap matrix:

$$\begin{equation}
S^{-\frac{1}{2}}HS^{-\frac{1}{2}}\phi = \lambda \phi.
\end{equation}$$

Performing this calculation is acceptable because we can reuse the inverse
over many scf loops, and because the overlap matrix is usually well conditioned.
In the canonical density matrix purification method of Palser[1] (which
we will probably discuss in more detail later), one also finds the matrix
inverse:

$$\begin{equation}
p_0 = \frac{\lambda}{2}(\mu S^{-1} - S^{-1}HS^{-1}) + \frac{1}{2}S^{-1}
\end{equation}$$

where $$p_0$$ is an initial guess at the density matrix.

Hotelling's method is very simple to implement. It works through the
following iteration:

$$\begin{equation}
X_{n+1} = 2X_{n} - X_{n}SX_{n}
\end{equation}$$

where $$\lim_{n \to \infty} X_n = S^{-1}$$. Palser cites the famous Numerical
Recipes book for this method. Disappointingly, it seems that this method
is called "Hotelling's Method" because it was invented by statistician Harold
Hotelling, and not because it leaves chocolates on your pillows. I've seen some
people mention reference [2] as the original paper. Hotelling himself notes that
it was "noticed" in reference [3].

The condition for convergence is (from [2]) $$ \|1 - SX_0 \| < 1$$. So
how should we pick an initial $$X_0$$? One simple way is to just scale the
initial matrix. Consider the eigendecomposition of our initial matrix:

$$\begin{equation}
S = UDU^T.
\end{equation}$$

Now let's plug the decomposition into the convergence condition, with a
scaling value $$\alpha$$:

$$\begin{equation}
\|1 - UDU^T\alpha UDU^T \| < 1.
\end{equation}$$

$$\begin{equation}
\|1 - \alpha DD \| < 1.
\end{equation}$$

The effect of multiplying the diagonal matrix of eigenvalues is to square all
the eigenvalues. Hence they are all positive. Thus is $$\alpha$$ is equal
to the inverse of the largest eigenvalue squared, we satisfy the equation.
The largest eigenvalue can be cheaply computed using the power method.

This is just a simple starting guess that I've introduced, there are better ones
out there in the literature, such as the guess in reference [4] which is
specific to overlap matrices (homework question: how can we improve our
$$\alpha$$ value for overlap matrices?).

Let's take a look at a simple implementation:
{% highlight python %}
# Libraries
from scipy.sparse.linalg import eigsh,norm
from scipy.sparse import identity, rand
from sys import argv

# Input Parameters
dimension = int(argv[1])
sparsity = float(argv[2])

# Initial Matrix
matrix = rand(dimension, dimension, density=sparsity, format="csr")
matrix = matrix + matrix.T

# Initial Guess
largest_eigenvalue = eigsh(matrix, k=1, which='LM',
  return_eigenvectors=False)
inverse_mat = matrix * 1.0/(largest_eigenvalue[0]**2)

# Iterate
for i in range(0,100):
  inverse_mat = 2*inverse_mat - \
    inverse_mat.dot(matrix.dot(inverse_mat))
  norm_value = norm(inverse_mat.dot(matrix) - identity(dimension))
  if norm_value < 1e-8:
    break

print norm(inverse_mat.dot(matrix) - identity(dimension))

{% endhighlight %}

The implementation highlights two further features of this algorithm. First,
the main computational kernel is matrix multiplication, which is great for
high performance computing. Second, it's trivial to extend to the sparse case
by simply replacing the dense multiplies with sparse multiplies.

> [1] Palser, Adam HR, and David E. Manolopoulos. "Canonical purification of the
> density matrix in electronic-structure theory." Physical Review B 58, no. 19
> (1998): 12704.

> [2] Hotelling, Harold. "Some new methods in matrix calculation." The Annals
> of Mathematical Statistics 14, no. 1 (1943): 1-34.

> [3] Frazer, Robert Alexander, William Jolly Duncan, Arthur Roderich Collar,
> and A. A. Mullin. "Elementary matrices and some applications to dynamics and
> differential equations." American Journal of Physics 29, no. 8 (1961): 555-556.

> [4] Ozaki, T. "Efficient recursion method for inverting an overlap matrix."
> Physical Review B 64, no. 19 (2001): 195110.
