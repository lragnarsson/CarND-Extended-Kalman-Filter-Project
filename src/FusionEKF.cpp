#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);


  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  // set initial uncertainty quite high:
  ekf_.P_ << 1, 0, 0,   0,
             0, 1, 0,   0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  debug_log_ = true;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */    
    cout << "EKF: " << endl;

    //ekf_.x_ = VectorXd(4);
    // If the first measurement is from LASER, then the velocity will be initialized to 0
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0; 

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // px = rho = cos(phi)
      ekf_.x_(0) = measurement_pack.raw_measurements_(0)     * cos(measurement_pack.raw_measurements_(1));
      // py = rho = sin(phi)
      ekf_.x_(1) = measurement_pack.raw_measurements_(0)     * sin(measurement_pack.raw_measurements_(1));
      // vx = rho_dot = cos(phi)
      ekf_.x_(2) = measurement_pack.raw_measurements_(2)     * cos(measurement_pack.raw_measurements_(1));
      // vy = rho_dot = cos(phi)
      ekf_.x_(3) = measurement_pack.raw_measurements_(2)     * sin(measurement_pack.raw_measurements_(1));
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
      ekf_.x_(2) = 0.0;
      ekf_.x_(3) = 0.0;
    }

    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  float Ts = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // Sample time for this iteration
  Ts = Ts > 1 ? 1 : Ts; // Prevent the Q-matrix from blowing up on start
  float Ts2 = Ts * Ts;
  float Ts3 = Ts * Ts * Ts / 2.0;
  float Ts4 = Ts * Ts * Ts * Ts / 4.0;

  previous_timestamp_ = measurement_pack.timestamp_;

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   // Update transition matrix with new sampel time:
   ekf_.F_(0, 2) = Ts;
   ekf_.F_(1, 3) = Ts;

  // Set process noise covariance: 
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << Ts4 * noise_ax, 0,              Ts3 * noise_ax, 0,
             0,              Ts4 * noise_ay, 0,              Ts3 * noise_ay,
             Ts3 * noise_ax, 0,              Ts2 * noise_ax, 0,
             0,              Ts3 * noise_ay, 0,              Ts2 * noise_ay;

  ekf_.Predict();
  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    if (debug_log_) cout << "Radar" << endl;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // TODO: Laser updates
    if (debug_log_) cout << "Laser" << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  if (debug_log_) {
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
    cout << "----------------------------------" << endl;
  }
}
