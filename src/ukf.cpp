#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1; // 1

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1; // 1

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // intially set to false
  is_initialized_ = false;

  // time tag in us
  time_us_ = 0;

    // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // predicted sigma points matrix
  Xsig_pred_= MatrixXd(n_x_, 2*n_aug_+1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_; // 3 is good for Gaussian

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill( 0.5/(n_aug_ + lambda_));
  weights_(0) = lambda_/(lambda_ + n_aug_);

  //some coefficients to avoid repeated computation
  sqrt_lambda_n_aug_ = sqrt(lambda_ + n_aug_);

  H_laser_ = MatrixXd(2, 5);
  H_laser_ << 1, 0, 0, 0, 0,
			        0, 1, 0, 0, 0;

  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;
  // NIS
  NIS_Lidar_ = 0;
  NIS_Radar_ = 0;



}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "UKF: " << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double r = meas_package.raw_measurements_[0];
      double cos_phi = cos(meas_package.raw_measurements_[1]);
      double sin_phi  = sin(meas_package.raw_measurements_[1]);
      x_ << r * cos_phi, r * sin_phi, 0, 0, 0;
      // Approximate P using Jacobian (first-order linearization)
      MatrixXd J_ = MatrixXd(2,2);
      J_ << cos_phi, -r * sin_phi, sin_phi, r * cos_phi;
      MatrixXd R_ = MatrixXd(2,2);
      R_ << std_radr_*std_radr_, 0, 0, std_radphi_ * std_radphi_;
      //update the up-left corner block
      P_.topLeftCorner(2,2) << J_ * R_ * J_.transpose();
      cout << "Radar Initialization" << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      //Initialize state from Laser.
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      P_(0,0) = std_laspx_ * std_laspx_;
      P_(1,1) = std_laspy_ * std_laspy_;
      P_(2,2) = 1;
      P_(3,3) = 1;
      P_(4,4) = 1;
      cout << "Laser Initialization" << endl;
    }
    cout << "x_init = " << x_ << endl;
    cout << "P_init = " << P_ << endl;

    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  /*****************************************************************************
  * Prediction & Update
  ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {

    double dt = (meas_package.timestamp_ - time_us_)/1000000.0;	//dt in seconds
    time_us_ = meas_package.timestamp_ ;
    //Preparing sigma points
    Prediction(dt);
    // Radar updates
    UpdateRadar(meas_package);
    time_us_ = meas_package.timestamp_;
    cout << "----------Radar Update----------" << endl;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    double dt = (meas_package.timestamp_ - time_us_)/1000000.0;	//dt in seconds
    time_us_ = meas_package.timestamp_ ;
    //Preparing sigma points
    Prediction(dt);
    // Laser updates
    UpdateLidar(meas_package);
    time_us_ = meas_package.timestamp_;
    cout << "----------Laser Update-----------" << endl;
  }
  x_(3) = Angle_Normalization(x_(3));

  // print the output
  cout << "x_posteri = " << x_ << endl;
  cout << "P_posteri = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt_lambda_n_aug_ * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_lambda_n_aug_ * L.col(i);
  }

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //some precomupted coefficients
    sin_yaw_ = sin(yaw);
    cos_yaw_ = cos(yaw);
    double delta_t2 = delta_t * delta_t;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin_yaw_);
        py_p = p_y + v/yawd * ( cos_yaw_ - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos_yaw_;
        py_p = p_y + v*delta_t*sin_yaw_;
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t2 * cos_yaw_;
    py_p = py_p + 0.5*nu_a*delta_t2 * sin_yaw_;
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t2;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  //compute projected mean
  x_ = Xsig_pred_ * weights_;
  x_(3) = Angle_Normalization(x_(3));

  //compute projected covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = Angle_Normalization(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  cout << "x_prior = " << x_ << endl;
  cout << "P_prior = " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  MatrixXd H = H_laser_;
  MatrixXd R = R_laser_;
  VectorXd z_pred = H * x_;
  VectorXd z = meas_package.raw_measurements_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(5, 5);
  P_ = (I - K * H) * P_;

  NIS_Lidar_ = y.transpose() * Si * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos_yaw_*v;
    double v2 = sin_yaw_*v;

    double eps = 1e-6;
    double r = max(sqrt(p_x*p_x + p_y*p_y), eps);

    // measurement model
    Zsig(0,i) = r;                                              //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / r;                         //r_dot
  }
  //compute projected mean
  VectorXd z_pred = Zsig * weights_;
  z_pred(1) = Angle_Normalization(z_pred(1));

  //measurement covariance matrix S and Cross Covariance matrix Tc
  MatrixXd S = MatrixXd(3,3);
  S.fill(0.0);

  MatrixXd Tc = MatrixXd(5, 3);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = Angle_Normalization(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = Angle_Normalization(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_radar_;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //innovation
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  //angle normalization
  y(1) = Angle_Normalization(y(1));

  //update state mean and covariance matrix
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();

}

double UKF::Angle_Normalization(double angle) {
  while (angle >  M_PI) angle -= 2.*M_PI;
  while (angle < -M_PI) angle += 2.*M_PI;
  return angle;
}
