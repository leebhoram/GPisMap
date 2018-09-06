//////////////////////////////////////////
// covFnc.h
// This header includes definitions of covariance functions for Gaussian Processes.
//
// by Bhoram Lee
#ifndef __COVFNC_H__
#define __COVFNC_H__

#include <vector>
#include <Eigen/Dense>

#include "strct.h"

typedef Eigen::MatrixXf EMatrixX;
typedef Eigen::VectorXf EVectorX;


//////////////////////////////////////////
// Convariance matrix computation using the Ornstein-Uhlenbeck cov function.
// Note:
//       - Dim(x) MUST BE 2 or 3 (See covFnc.cpp)
//       - Used for GP Implicit Surface (GPIS) jointly with Gradient
//       - gradflag indicates whether the gradient is available or not for each point
//////////////////////////////////////////

// covariances for x1 (input points) with different noise params for inputs
EMatrixX matern32_sparse_deriv1(EMatrixX const& x1, std::vector<float> gradflag,
                                           float scale_param, EVectorX const& sigx, EVectorX const& siggrad);

// covariances for x1 (input points) and x2 (test points)
EMatrixX matern32_sparse_deriv1(EMatrixX const& x1, std::vector<float> gradflag,
                                          EMatrixX const& x2, float scale_param);

//////////////////////////////////////////
// Convariance matrix computation using the Ornstein-Uhlenbeck cov function.
// Note:
//       - Dim(x) > 0
//       - Used for Observation Regression
//////////////////////////////////////////

// covariances for x1 (input points) with a constanct noise param
EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, float scale_param, float sigx);

// covariances for x1 (input points) with different noise params for inputs
EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, float scale_param, EVectorX const& sigx);

// covariances for x1 (input points) and x2 (test points)
EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, EMatrixX const& x2, float scale_param);

#endif
