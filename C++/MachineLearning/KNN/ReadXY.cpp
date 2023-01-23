#include <fstream>
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>
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



