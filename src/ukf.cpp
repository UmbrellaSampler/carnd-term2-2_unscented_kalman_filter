#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

    // set initially to false
    is_initialized_ = false;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);
    // initial covariance matrix
    P_ = MatrixXd(5, 5);
    P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1000, 0, 0,
            0, 0, 0, 100, 0,
            0, 0, 0, 0, 1;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

    /**
     * DO NOT MODIFY measurement noise values below.
     * These are provided by the sensor manufacturer.
     */

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

    /**
     * End DO NOT MODIFY section for measurement noise values
     */

    // State dimension
    n_x_ = 5;

    // Augmented state dimension
    n_aug_ = 7;

    // Sigma point spreading parameter
    lambda_ = n_aug_ - 3;

    // augmented sigma points matrix
    Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    // predicted sigma points as columns
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    // Measurement state dimension for radar
    n_z_radar_ = 3;

    // Measurement state dimension for laser
    n_z_laser_ = 2;

    R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
    R_radar_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;

    R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
    R_laser_ << std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;

    // init weights_
    weights_ = VectorXd(2 * n_aug_ + 1);   // 15 sigma points
    weights_(0) = double(lambda_ / (lambda_ + n_aug_));
    for (int i = 1; i < 2 * n_aug_ + 1; i++) {
        weights_(i) = double(0.5 / (lambda_ + n_aug_));
    }

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
     * Initialization
     */

    if (!is_initialized_) {

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            const auto &meas = meas_package.raw_measurements_;
            auto x_cart = tools_.PolarToCartesian(meas);

            auto px = x_cart(0);
            auto py = x_cart(1);
            auto vx = x_cart(2);
            auto vy = x_cart(3);

            x_ << px, py, sqrt(vx * vx + vy * vy), 0.0, 0.0;
        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            const auto &meas = meas_package.raw_measurements_;

            auto p_x = meas(0);
            auto p_y = meas(1);
            x_ << p_x, p_y, 0.0, 0.0, 0.0;
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        // set previous timestamp
        prev_timestamp_ = meas_package.timestamp_;
        return;
    }

    // set time difference in seconds
    dt_ = (meas_package.timestamp_ - prev_timestamp_) / 1000000.0;
    // set previous timestamp
    prev_timestamp_ = meas_package.timestamp_;

    /***********
     * Predict
     ***********/
    Prediction();

    /***********
     * Update
     ***********/
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    }
}

void UKF::Prediction() {
    AugmentedSigmaPoints();
    SigmaPointPrediction();
    PredictMeanAndCovariance();
}


void UKF::UpdateLidar(MeasurementPackage meas_package) {
    MatrixXd Zsig = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        Zsig(0, i) = Xsig_pred_(0, i); // px
        Zsig(1, i) = Xsig_pred_(1, i); // py
    }

    Update(meas_package.raw_measurements_, Zsig, R_laser_);
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
    MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
    // transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // extract values for better readability
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        // measurement model
        Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y); // r
        Zsig(1, i) = atan2(p_y, p_x); // phi
        Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y); // r_dot
    }
    Update(meas_package.raw_measurements_, Zsig, R_radar_);
}

void UKF::AugmentedSigmaPoints() {
    // create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    // create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;

    // create square root matrix
    const MatrixXd L = P_aug.llt().matrixL();

    // create augmented sigma points
    Xsig_aug_.col(0) = x_aug;
    for (int i = 0; i < n_aug_; ++i) {
        Xsig_aug_.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug_.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }
}

void UKF::SigmaPointPrediction() {

    // predict sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        // extract values for better readability
        double p_x = Xsig_aug_(0, i);
        double p_y = Xsig_aug_(1, i);
        double v = Xsig_aug_(2, i);
        double yaw = Xsig_aug_(3, i);
        double yawd = Xsig_aug_(4, i);
        double nu_a = Xsig_aug_(5, i);
        double nu_yawdd = Xsig_aug_(6, i);

        // predicted state values
        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v / yawd * (sin(yaw + yawd * dt_) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * dt_));
        } else {
            px_p = p_x + v * dt_ * cos(yaw);
            py_p = p_y + v * dt_ * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd * dt_;
        double yawd_p = yawd;

        // add noise
        px_p = px_p + 0.5 * nu_a * dt_ * dt_ * cos(yaw);
        py_p = py_p + 0.5 * nu_a * dt_ * dt_ * sin(yaw);
        v_p = v_p + nu_a * dt_;

        yaw_p = yaw_p + 0.5 * nu_yawdd * dt_ * dt_;
        yawd_p = yawd_p + nu_yawdd * dt_;

        // write predicted sigma point into right column
        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;
    }
}

void UKF::PredictMeanAndCovariance() {
    // predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    // predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
}

void UKF::Update(const VectorXd &z_meas, const MatrixXd &Zsig, const MatrixXd &R) {

    int n_z = z_meas.size();
    VectorXd z_pred = VectorXd(n_z);
    // mean predicted measurement
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);

    // cross correlation matrix Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);

    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    S = S + R;

    // Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    // residual
    VectorXd z_diff = z_meas - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    // update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
}
