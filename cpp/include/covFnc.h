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

#ifdef __USE_DOUBLE_PRECISION__     // This is in use
typedef Eigen::MatrixXd EMatrixX;
typedef Eigen::VectorXd EVectorX;
#else                               // For potential transition to float
typedef Eigen::MatrixXf EMatrixX;
typedef Eigen::VectorXf EVectorX;
#endif

//////////////////////////////////////////
// Convariance matrix computation using the Ornstein-Uhlenbeck cov function.
// Note:
//       - Dim(x) MUST BE 2 or 3 (See covFnc.cpp)
//       - Used for GP Implicit Surface (GPIS) jointly with Gradient
//       - gradflag indicates whether the gradient is available or not for each point
//////////////////////////////////////////

// covariances for x1 (input points) with different noise params for inputs
EMatrixX matern32_sparse_deriv1(EMatrixX const& x1, std::vector<FLOAT> gradflag,
                                           FLOAT scale_param, EVectorX const& sigx, EVectorX const& siggrad);

// covariances for x1 (input points) and x2 (test points)
EMatrixX matern32_sparse_deriv1(EMatrixX const& x1, std::vector<FLOAT> gradflag,
                                          EMatrixX const& x2, FLOAT scale_param);


//////////////////////////////////////////
// Convariance matrix computation using the Ornstein-Uhlenbeck cov function.
// Note:
//       - Dim(x) > 0
//       - Used for Observation Regression
//////////////////////////////////////////

// covariances for x1 (input points) with a constanct noise param
EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, FLOAT scale_param, FLOAT sigx);

// covariances for x1 (input points) with different noise params for inputs
EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, FLOAT scale_param, EVectorX const& sigx);

// covariances for x1 (input points) and x2 (test points)
EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, EMatrixX const& x2, FLOAT scale_param);

#endif
