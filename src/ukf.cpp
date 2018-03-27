#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"
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

    time_us_ = 0;

    //set state dimension
    n_x_ = 5;

    //set augmented dimension
    n_aug_ = 7;

    //define spreading parameter
    lambda_ = 3 - n_aug_;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(n_x_);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2.90;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = .75;

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

    //create augmented mean vector
    x_aug_ = VectorXd(n_aug_);

    //create augmented state covariance
    P_aug_ = MatrixXd(n_aug_, n_aug_);

    //create sigma point matrix
    Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    //create matrix with predicted sigma points as columns
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    //create vector for weights
    weights_ = VectorXd(2 * n_aug_ + 1);

    //create vector for predicted state
    x_predict_ = VectorXd(n_x_);

    //create covariance matrix for prediction
    P_predict_ = MatrixXd(n_x_, n_x_);

    //set measurement dimension, radar can measure r, phi, and r_dot
    n_z_radar_ = 3;

    //set measurement dimension, lidar can measure px and py
    n_z_lidar_ = 2;

    // set weights
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    weights_(0) = weight_0;
    for (int i = 1; i < 2 * n_aug_ + 1; i++) { //2n+1 weights
        double weight = 0.5 / (n_aug_ + lambda_);
        weights_(i) = weight;
    }
}

UKF::~UKF() {
}

Eigen::VectorXd h_Radar(const Eigen::VectorXd x_) {
    float px = x_(0); // position in x
    float py = x_(1); // position in y
    float vx = x_(2); // velocity in x
    float vy = x_(3); // velocity in y
    // Coordinates conversion from cartesian to polar
    float rho = sqrt(pow(px, 2) + pow(py, 2));
    if (rho < 0.0001) {
        px += 0.0001;
        py += 0.0001;
        rho = sqrt(px * px + py * py);
    }
    float phi = atan2(py, px);
    float dot = (px * vx + py * vy) / rho;

    VectorXd h(3);
    h << rho, phi, dot;
    return h;
}

Eigen::VectorXd h_inv_Radar(const Eigen::VectorXd x_) {
    double rho = x_[0]; // range
    double phi = x_[1]; // bearing
    // Coordinates conversion from polar to cartesian
    double x = rho * cos(phi);
    if (x < 0.0001) {
        x = 0.0001;
    }
    double y = rho * sin(phi);
    if (y < 0.0001) {
        y = 0.0001;
    }

    VectorXd h_inv(5);
    h_inv << x, y, 0, 0, 0;
    return h_inv;
}

void UKF::estimate(MeasurementPackage& measurement_pack) {
    AugmentedSigmaPoints();
    long long dt = measurement_pack.timestamp_ - time_us_;
    SigmaPointPrediction(dt / 1000000.0);
    time_us_ = measurement_pack.timestamp_;
    PredictMeanAndCovariance();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2], 0, 0;
            x_ = h_inv_Radar(x_);
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
        }
        // done initializing, no need to predict or update
        time_us_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }
    /*****************************************************************************
     *  Prediction & Update
     ****************************************************************************/
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        estimate(measurement_pack);
        PredictMeasurement(measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2]);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        estimate(measurement_pack);
        PredictMeasurement(measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1]);
    }
}

void UKF::AugmentedSigmaPoints() {
    //create augmented mean state
    x_aug_.head(5) = x_;
    x_aug_(5) = 0;
    x_aug_(6) = 0;

    //create augmented covariance matrix
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(5, 5) = P_;
    P_aug_(5, 5) = std_a_*std_a_;
    P_aug_(6, 6) = std_yawdd_*std_yawdd_;

    //create square root matrix
    MatrixXd L = P_aug_.llt().matrixL();

    //create augmented sigma points
    Xsig_aug_.col(0) = x_aug_;
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
    }
}

void UKF::SigmaPointPrediction(double delta_t) {
    //predict sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //extract values for better readability
        double p_x = Xsig_aug_(0, i);
        double p_y = Xsig_aug_(1, i);
        double v = Xsig_aug_(2, i);
        double yaw = Xsig_aug_(3, i);
        double yawd = Xsig_aug_(4, i);
        double nu_a = Xsig_aug_(5, i);
        double nu_yawdd = Xsig_aug_(6, i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        } else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;
    }
}

void UKF::PredictMeanAndCovariance() {
    //predicted state mean
    x_predict_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //iterate over sigma points
        x_predict_ = x_predict_ + weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P_predict_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_predict_;
        //angle normalization
        x_diff(3) = Tools::normalizeAngle(x_diff(3));

        P_predict_ = P_predict_ + weights_(i) * x_diff * x_diff.transpose();
    }
}

void UKF::PredictMeasurement(double rho, double phi, double rho_dot) {
    //create matrix for sigma points in measurement space
    Zsig_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

    //mean predicted measurement
    z_pred_ = VectorXd(n_z_radar_);

    //innovation covariance matrix S
    S_ = MatrixXd(n_z_radar_, n_z_radar_);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        // measurement model
        Zsig_(0, i) = sqrt(p_x * p_x + p_y * p_y); //r
        Zsig_(1, i) = atan2(p_y, p_x); //phi
        Zsig_(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y); //r_dot
    }


    z_pred_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
    }


    S_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig_.col(i) - z_pred_;

        //angle normalization
        z_diff(1) = Tools::normalizeAngle(z_diff(1));

        S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z_radar_, n_z_radar_);
    R << std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_*std_radrd_;
    S_ = S_ + R;
    UpdateState(rho, phi, rho_dot);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateState(double rho, double phi, double rho_dot) {

    //create example vector for incoming radar measurement
    VectorXd z_ = VectorXd(n_z_radar_);
    z_ <<
            rho,
            phi,
            rho_dot;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

    /*******************************************************************************
     * Student part begin
     ******************************************************************************/

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //2n+1 simga points

        //residual
        VectorXd z_diff = Zsig_.col(i) - z_pred_;
        //angle normalization
        //        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        //        while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;
        z_diff(1) = Tools::normalizeAngle(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_predict_;
        //angle normalization
        //        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        //        while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;
        x_diff(3) = Tools::normalizeAngle(x_diff(3));

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S_.inverse();

    //residual
    VectorXd z_diff = z_ - z_pred_;

    //angle normalization
    //    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    //    while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;
    z_diff(1) = Tools::normalizeAngle(z_diff(1));

    //update state mean and covariance matrix
    x_ = x_predict_ + K * z_diff;
    P_ = P_predict_ - K * S_ * K.transpose();

    /*******************************************************************************
     * Student part end
     ******************************************************************************/

}

void UKF::PredictMeasurement(double px, double py) {


    /*******************************************************************************
     * Student part begin
     ******************************************************************************/
    //create matrix for sigma points in measurement space
    Zsig_ = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);

    //mean predicted measurement
    z_pred_ = VectorXd(n_z_lidar_);

    //innovation covariance matrix S
    S_ = MatrixXd(n_z_lidar_, n_z_lidar_);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);

        // measurement model
        Zsig_(0, i) = p_x; //r
        Zsig_(1, i) = p_y; //phi
    }


    z_pred_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
    }


    S_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig_.col(i) - z_pred_;

        //angle normalization
        //        while (z_diff(1) > M_PI) z_diff(1) -= 2.f * M_PI;
        //        while (z_diff(1)<-M_PI) z_diff(1) += 2.f * M_PI;
        z_diff(1) = Tools::normalizeAngle(z_diff(1));

        S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z_lidar_, n_z_lidar_);
    R << std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;
    S_ = S_ + R;


    UpdateState(px, py);
}

void UKF::UpdateState(double px, double py) {

    //create example vector for incoming radar measurement
    VectorXd z_ = VectorXd(n_z_lidar_);
    z_ <<
            px,
            py;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z_lidar_);

    /*******************************************************************************
     * Student part begin
     ******************************************************************************/

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //2n+1 simga points

        //residual
        VectorXd z_diff = Zsig_.col(i) - z_pred_;
        //angle normalization
        //        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        //        while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;
        z_diff(1) = Tools::normalizeAngle(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_predict_;
        //angle normalization
        //        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        //        while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;
        x_diff(3) = Tools::normalizeAngle(x_diff(3));

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S_.inverse();

    //residual
    VectorXd z_diff = z_ - z_pred_;

    //angle normalization
    //    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    //    while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;
    z_diff(1) = Tools::normalizeAngle(z_diff(1));

    //update state mean and covariance matrix
    x_ = x_predict_ + K * z_diff;
    P_ = P_predict_ - K * S_ * K.transpose();

    /*******************************************************************************
     * Student part end
     ******************************************************************************/

}
