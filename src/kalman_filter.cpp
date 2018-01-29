#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {


  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

	x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

	VectorXd y = z - H_ * x_;
	MatrixXd S = H_ * P_*H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();
	P_ = (I - K * H_) * P_;

	x_ = x_ + K * y;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

	VectorXd h_x = RadarMeasEstimation(x_);

	double dPhi = z(1) - h_x(1);

	// deal with the case where prediction and measured value straddle the pi/-pi transition
	if (dPhi > M_PI)
	{
		h_x(1) = h_x(1) + 2 * M_PI;
	}
	if (dPhi < -M_PI)
	{
		h_x(1) = h_x(1) - 2 * M_PI;
	}
	VectorXd y = z - h_x;
	MatrixXd S = H_ * P_*H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();

	x_ = x_ + K * y;
	P_ = (I - K * H_) * P_;

}

VectorXd KalmanFilter::RadarMeasEstimation(VectorXd &x)
{
	VectorXd h = VectorXd(3);

	double px = x[0];
	double py = x[1];
	double vx = x[2];
	double vy = x[3];


	double rho = 0;
	double phi = 0;
	double rhoDot = 0;

	rho = sqrt(px*px + py * py);
	phi = atan2(py, px);
	rhoDot = (px * vx + py * vy) / rho;

	h << rho, phi, rhoDot;

	return h;

}

void KalmanFilter::UpdateF(double dt)
{
	MatrixXd F = MatrixXd(4,4);
	F << 1,0,dt,0,
		 0,1,0 ,dt,
		 0,0,1 ,0,
		 0,0,0 ,1;
	F_ = F;

}

void KalmanFilter::UpdateQ(double dt, double noise_ax, double noise_ay)
{
  MatrixXd Q = MatrixXd(4,4);

  double dt_2 = dt*dt;
  double dt_3 = dt_2*dt;
  double dt_4 = dt_3*dt;




  Q << dt_4/4*noise_ax, 0              , dt_3/2*noise_ax, 0,
        0              , dt_4/4*noise_ay, 0              , dt_3/2*noise_ay,
        dt_3/2*noise_ax, 0              , dt_2*noise_ax  , 0,
        0              , dt_3/2*noise_ay, 0              , dt_2*noise_ay;
  Q_ = Q;
}

void KalmanFilter::UpdateHj(void)
{
  MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    //check division by zero
    if(fabs(c1) < 0.0001){
      c1 = 0.0001;
    }

    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    H_ =  Hj;
}
