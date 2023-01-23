#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include "armadillo.hpp"

arma::mat ReadX(std::string const input_x){
	
	std::ifstream file_x(input_x);
	assert(file_x.is_open());
	arma::mat points_matrix;

	int i=0;
	for(std::string line; getline(file_x, line);){
		std::stringstream ss(line);
		arma::rowvec point(34);
		int j=0;
		for ( double d; ss >> d; ){
			point(j) = d;
			j++;
		}
		points_matrix.insert_rows(i, point);
		i++;
	}
	file_x.close();
	return points_matrix;
}


arma::rowvec ReadY(std::string const input_y, int n){

	std::ifstream file_y(input_y);
	assert(file_y.is_open());
	arma::rowvec labels(n);

	int i=0;
	for (std::string point; getline(file_y, point);){
		labels(i) = atof(point.c_str());
		i++;
	}
	file_y.close();
	return labels;
}



arma::rowvec gradientDescent(arma::rowvec y, arma::mat x, arma::rowvec w, double alpha, double epsilon){
	int n = y.size();
	int len = w.size();
	double diff_norm = 10;
	double normalize_const = -1.0/(1.0*n);
	int i = 0;
	do {
		arma::rowvec differential(len, arma::fill::zeros);
		for (int j=0; j<n; j++) {
			arma::rowvec xi = x.row(j);
			double yi = y(j);
			differential = differential + normalize_const * ( yi * 1/(1 + std::exp(yi * dot(w, xi.t()))) * xi );
		}
		w = w - alpha * differential;
		diff_norm = norm(differential, 2);
		i++;
	} while ( epsilon < diff_norm );

	std::cout << "gradient descent done in: " << i << " iterations \n";
	
	return w;
};


int assignLabel(arma::rowvec x, arma::rowvec w){
	double prod = dot(w, x.t());
	if (prod < 0){ return -1; }
	else { return 1; }
};



int main()
{
	std::string input_x = "dataX.dat";
	std::string input_y = "dataY.dat";
	std::string input_test = "dataXtest.dat";

	double alpha = 0.9;
	double epsilon = 1e-7;
	
	arma::mat X = ReadX(input_x);
	int n = X.n_rows;
	arma::rowvec Y(n);
	Y = ReadY(input_y, n);
	
	arma::mat X_test = ReadX(input_test);
	int n_test = X_test.n_rows;

	arma::rowvec weight(34);
	
	weight = gradientDescent(Y, X, weight, alpha, epsilon);
	
	std::ofstream write_out("LogReg.dat");
	assert(write_out.is_open());

	for (int i=0; i<n_test; i++) {
		arma::rowvec test_row(34); test_row = X_test.row(i);
		int point_label = assignLabel(test_row, weight);
		write_out << point_label << "\n";
	}


	return 0;
}
