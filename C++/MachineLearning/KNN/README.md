# K-Nearest Neighbours

This project is an implementation of k-nearest neighbours in C++.

The nearest neighbour algorithm is a simple algorithm that assigns a label using a majority voting scheme. A previously unseen point is compared to a set of points with known labels. The algorithm picks the $k$ points with closest distance and assigns the label of the majority of points. More formally, given a set of points ${(x_1, y_i), \dots, (x_N, y_N)}$ and a number of neighbours, $k$ we assign a label $y$ to a new point $x$ as:

* Compute distances $d_i = \lVert x_i - x \rVert$
* Create a set of label distance pairs $P$ with $P_i = (y_i, d_i), i = 1 \dots N$
* Sort $P$ ascending by $d_i$ (i.e. afterwards $P_1$ is the pair with smallest $d_i$)
* Assign to $y$ the label the label which appears most often among $P_1 \dots P_k$

### C++ implementation

See [code](NearestNeighbours.cpp).

The code will:

* Read "dataX.dat", "dataXtest.dat" and "dataY.dat" from disk. "dataX.dat" contains the points $x_i$ in its rows. They are stored as numbers separated by a white space (i.e. tab or space). "dataY.dat" contains the label $y_i$, one in each row. We assume that the only labels appearing are 1 or -1. "dataXtest.dat" contains test points. 
* For each point in "dataXtest.dat", assign the apropriate label using the nearest neighbour algorithm with $k = 5$. 
* Write a new file "NN.dat" containing the assigned labels in the same format as "dataY.dat".



### How to run:

I have used clang++ as compiler:

  `clang++ NearestNeighbours.cpp -o NearestNeihbours`

### Dependencies

This project makes use of [Armadillo](https://arma.sourceforge.net/) for linear algebra implementation.


