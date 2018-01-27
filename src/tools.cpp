#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  cout << "rmse" << endl;
	VectorXd rmse = VectorXd(4);
	rmse << 0, 0, 0, 0;

	if (estimations.size() != ground_truth.size())
	{
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}
	if (estimations.size() == 0)
	{
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}



	//accumulate squared residuals

	VectorXd sum_residuals_squared(4);
	sum_residuals_squared << 0, 0, 0, 0;

	for (int i = 0; i < estimations.size(); ++i) {
		// ... your code here

		VectorXd diff = estimations[i] - ground_truth[i];


		VectorXd residuals_squared = diff.array() * diff.array();

		sum_residuals_squared += residuals_squared;

	}

	//calculate the mean
	// ... your code here
	VectorXd mean = sum_residuals_squared / estimations.size();

	//calculate the squared root
	// ... your code here
	rmse = mean.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3, 4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float px_2 = px * px;
	float py_2 = py * py;
	float px_py_2_3 = sqrt((px_2 + py_2))*(px_2 + py_2);



	//TODO: YOUR CODE HERE 

	//check division by zero
	if (fabs(px_2 + py_2) < 0.0001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix

	Hj << px / sqrt(px_2 + py_2), py / sqrt(px_2 + py_2), 0, 0,
		-py / (px_2 + py_2), px / (px_2 + py_2), 0, 0,
		py*(vx*py - vy * px) / px_py_2_3, px*(vy*px - vx * py) / px_py_2_3, px / sqrt(px_2 + py_2), py / sqrt(px_2 + py_2);

	return Hj;
}
