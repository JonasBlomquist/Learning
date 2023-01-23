# Logistic Regression

This project is an implementation of Logistic Regression in C++.

Logistic regression assigns the label $y$ to a point $x$ using a previously chosen linear function $f(x) = w^T x$. To find the label, we use $y = sign(f(x))$ (i.e. if $f(x)$ is positive, the label is $1$, else the label is $-1$). We find $w$ by minimizing the function:

$$
    L(w) = \frac{1}{N} \sum_{i=1}^{N} \log (1 + \exp(-y_i w^T x_i))
$$

For this project we use gradient descent to find the optimal $w$. This algorithm uses the gradient of $L(w)$

$$
    \frac{\partial}{\partial w} L(w) = - \frac{1}{N} \sum_{i=1}^{N} y_i \frac{1}{1 + \exp(y_i w^T x_i} x_i
$$

That means for an initial guess of $w$ (e.g. $w=0$), we can find a better estimate using

$$
  w \leftarrow w - \alpha \cdot \frac{\partial}{\partial w} L(w)
$$

Where $\alpha \in [0,1]$ is a chosen learning rate. We repeat this until our solution is close enough to the optimal solution. We define our algorithm to be converged when $\lVert \frac{\partial}{partial w} L(w) \rVert < \epsilon$ for a chosen tolerance of $\epsilon$.

### C++ implementation

See [code](LogisticRegression.cpp).

The code will:

* Read "dataX.dat", "dataXtest.dat" and "dataY.dat" from disk. 
* Find the optimal $w$ using the rules above using $\epsilon = 10e-7$ 
* For each point in "dataXtest.dat", assign the apropriate label using the assignment rule $y = sign(f(x))$. 
* Write a new file "LogReg.dat" containing the assigned labels in the same format as "dataY.dat". 


### How to run:

I have used clang++ as compiler:

  `clang++ LogisticRegression.cpp *.o -o LogisticRegression`


### Dependencies

This project makes use of [Armadillo](https://arma.sourceforge.net/) for linear algebra implementation.
