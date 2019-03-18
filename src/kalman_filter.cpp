#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {debug_log_ = true;}

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
  MatrixXd Ft = F_.transpose();
  x_ = F_ * x_;
  P_ = F_ * P_ * Ft + Q_;

  if (debug_log_) {
    cout << "Q: " << Q_ << endl;
    cout << "------------------" << endl;
  }
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_hat = H_ * x_;
  VectorXd y = z - z_hat;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  // Calculate optimal Kalman Gain:
  MatrixXd K = PHt * Si;

  if (debug_log_) {
    cout << "x: " << x_ << endl;
    cout << "z_hat: " << z_hat << endl;
    cout << "y: " << y << endl;
    cout << "R: " << R_ << endl;
    cout << "K: " << K << endl;
    cout << "P: " << P_ << endl;
    cout << "------------------" << endl;
  }

  // Fused state estimate
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  // New state covariance:
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  // Predict next state z_hat in polar coordinates
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float rho, rho_dot, phi;
  rho = sqrt(px * px + py * py);
  if(rho < 0.00001) {
    rho = 0.00001;
    if (debug_log_) cout << "WARNING: Small Rho" << endl;
  }

  rho_dot = (px * vx + py * vy) / rho;
  phi = atan2(py, px);

  VectorXd z_hat  = VectorXd(3);
  z_hat << rho, phi, rho_dot;

  // Perform Kalman update:
  VectorXd y = z - z_hat;
  // Set  -pi/2 < Phi pi/2
  y(1) = BoundedAngle(y(1));

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  // Calculate optimal Kalman Gain:
  MatrixXd K = PHt * Si;

  if (debug_log_) {
    cout << "x: " << x_ << endl;
    cout << "z_hat: " << z_hat << endl;
    cout << "y: " << y << endl;
    cout << "R: " << R_ << endl;
    cout << "K: " << K << endl;
    cout << "P: " << P_ << endl;
    cout << "------------------" << endl;
  }

  // Fused state estimate
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  // New state covariance:
  P_ = (I - K * H_) * P_;

}

float KalmanFilter::BoundedAngle(float angle) {
  while (angle > M_PI) {
    angle -= 2 * M_PI;
  }
  while (angle < -M_PI) {
    angle += 2 * M_PI;
  }
  return angle;
}