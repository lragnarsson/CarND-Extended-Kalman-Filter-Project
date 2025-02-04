#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // Verify that the estimation vector size is not zero
  // Verify that the estimation vector size is equal to the ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  //cout << "RMSE: " << rmse << endl;
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  // Based on EKF lesson
  MatrixXd Hj(3,4);
  Hj.fill(0.0);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float c1 = px * px + py * py;

  // check division by zero
  if (fabs(c1) < 0.0001) {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    //return Hj;
    c1 += c1 > 0 ? 0.0001 : -0.0001;
  }

  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  // compute the Jacobian matrix
  Hj <<  (px/c2),               (py/c2),               0,     0,
        -(py/c1),               (px/c1),               0,     0,
         py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
