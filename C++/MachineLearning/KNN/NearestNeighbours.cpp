#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <vector>
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
};


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
};


struct points {
	double distance;
	int label;
};


bool comparePoints(points p1, points p2){
	return (p1.distance < p2.distance);
};

int assignLabel(std::vector<points> p, int const k){
	sort(p.begin(), p.end(), comparePoints);
	std::vector<points>::iterator start = p.begin();
	std::vector<points>::iterator end = p.begin() + k;
	std::vector<points> sub(k); 
	std::copy(start, end, sub.begin());
	int label_sum = 0;
	for (auto x : sub)
		label_sum += x.label;
	if (label_sum < 0){ return -1; }
	else { return 1; }
};


int main()
{
	int k=5;

	std::string input_x = "dataX.dat";
	arma::mat X = ReadX(input_x);
	
	int n = X.n_rows;

	std::string input_y = "dataY.dat";
	arma::rowvec Y(n);
	Y = ReadY(input_y, n);
	
	std::string input_test = "dataXtest.dat";
	arma::mat X_test = ReadX(input_test);
	int n_test = X_test.n_rows;
	
	std::ofstream write_out("NN.dat");
	assert(write_out.is_open());

	for(int i=0; i<n_test; i++){
		arma::rowvec test_row(34); test_row = X_test.row(i);
		arma::mat diff(n, 34); diff = X.each_row() - test_row;
		
		std::vector<points> labelling_vec;

		for(int j=0; j<n; j++){
			points point;
			point.distance = norm(diff.row(j), 2);
			point.label = Y(j);
			labelling_vec.push_back(point);
		}
		int point_label = assignLabel(labelling_vec, k);
		write_out << point_label << "\n";

	}
		

	return 0;
}
