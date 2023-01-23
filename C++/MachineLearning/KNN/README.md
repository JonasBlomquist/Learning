# K-Nearest Neighbours

This project is an implementation of k-nearest neighbours in C++.

### Algorithm

* Compute distances $d_i = \lvert x_i - x \rvert$
* Create a set of label distance pairs $P$ with $P_i = (y_i, d_i), i = 1 \hdots N$



### How to run:

**Input files:** Training set (commaseparated file; see example_training.csv), Test set (commaseparated file of same format as training set, but without class)

**Output:** Prediction file - Predicted class for all points in test set file.


### C++ implementation

See [code](NearestNeighbours.cpp)
