#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if (estimations.size() != ground_truth.size()
        || estimations.size() == 0) {
        std::cout << "Invalid estimation or ground_truth data" << std::endl;
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
    return rmse;
}

VectorXd Tools::PolarToCartesian(const VectorXd &x_polar) const {
    const auto rho = x_polar(0);
    const auto theta = x_polar(1);
    const auto rho_dot = x_polar(2);
    const auto px = rho * cos(theta);
    const auto py = rho * sin(theta);
    const auto vx = rho_dot * cos(theta);
    const auto vy = rho_dot * sin(theta);

    VectorXd x_cart(4);
    x_cart << px, py, vx, vy;

    return x_cart;
}