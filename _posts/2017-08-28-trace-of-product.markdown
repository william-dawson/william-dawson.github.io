---
layout: post
title:  "Computing The Trace of a Product of Matrices"
date:   2017-08-30 11:00:00 +0900
categories: blog method
---

I'm going to cover a very simple concept here, but it's one that I wish I had
learned earlier. In a number of algorithms, we're confronted with the problem
of computing the trace of a product of matrices:

$$\begin{equation}
\DeclareMathOperator{\Tr}{Tr}
\Tr(AB).
\end{equation}$$

One example of this is in quantum chemistry, where we're interested in computing
the energy of the system:

$$\begin{equation}
\DeclareMathOperator{\Tr}{Tr}
\text{energy} = \Tr(\rho H)
\end{equation}$$

where $$\rho$$ is the density matrix, and $$H$$ the hamiltonian.

What's the big deal? Well, let's say we're running an algorithm that iteratively
computes $$\rho$$. One way to check if we're converged would be to compute the
energy at each update. What is the cost of doing this? Naively, we would
need to multiply two matrices, and compute the trace. Unfortunately, those
multiplications add up, and can significantly increase the run time.
Fortunately, there's a very simple way to avoid this.

|A|A|A|A|A|.........|B|B|B|B|B|.........|C|C|C|C|C|
|A|A|A|A|A|...X.X..|B|B|B|B|B|..===..|C|C|C|C|C|
|A|A|A|A|A|....X....|B|B|B|B|B|.........|C|C|C|C|C|
|A|A|A|A|A|...X.X..|B|B|B|B|B|..===..|C|C|C|C|C|
|A|A|A|A|A|.........|B|B|B|B|B|.........|C|C|C|C|C|

To compute the trace, we only need to compute the diagonal elements of the
matrix. How do we compute a given diagonal element?

|O|O|O|O|O|.........|O|B|O|O|O|.........|O|O|O|O|O|
|A|A|A|A|A|...X.X..|O|B|O|O|O|..===..|O|C|O|O|O|
|O|O|O|O|O|....X....|O|B|O|O|O|.........|O|O|O|O|O|
|O|O|O|O|O|...X.X..|O|B|O|O|O|..===..|O|O|O|O|O|
|O|O|O|O|O|.........|O|B|O|O|O|.........|O|O|O|O|O|

And in the case of symmetric matrices:

|O|O|O|O|O|.........|O|O|O|O|O|.........|O|O|O|O|O|
|A|A|A|A|A|...X.X..|B|B|B|B|B|..===..|O|C|O|O|O|
|O|O|O|O|O|....X....|O|O|O|O|O|.........|O|O|O|O|O|
|O|O|O|O|O|...X.X..|O|O|O|O|O|..===..|O|O|O|O|O|
|O|O|O|O|O|.........|O|O|O|O|O|.........|O|O|O|O|O|

As we can see then, each diagonal element can be computed by just taking the
dot product of each row (element-wise multiplication). And for the full trace,
we just sum up those dot products. This means not only can we avoid performing
the matrix multiplication, but we can even do this perfectly in parallel, and
avoid any communication besides the final all-reduce of a single value.

{% highlight python %}
# Libraries
from sys import argv
from numpy import multiply, trace
from numpy import sum as np_sum
from numpy.random import rand

# Make Matrices
matrix_dimension = int(argv[1])
MatrixA = rand(matrix_dimension, matrix_dimension)
MatrixA = MatrixA + MatrixA.T
MatrixB = rand(matrix_dimension, matrix_dimension)
MatrixB = MatrixB + MatrixB.T

# Compute Traces
check_value = trace(MatrixA.dot(MatrixB))
new_value = np_sum(multiply(MatrixA, MatrixB))

# Check Answer
print(check_value)
print(new_value)
{% endhighlight %}
