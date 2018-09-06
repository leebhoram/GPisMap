/*
 * GPisMap - Online Continuous Mapping using Gaussian Process Implicit Surfaces
 * https://github.com/leebhoram/GPisMap
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License v3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of any FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License v3 for more details.
 *
 * You should have received a copy of the GNU General Public License v3
 * along with this program; if not, you can access it online at
 * http://www.gnu.org/licenses/gpl-3.0.html.
 *
 * covFnc.cpp
 * This file implements covariance functions for Gaussian Processes.
 *
 * Authors: Bhoram Lee <bhoram.lee@gmail.com>
 */
#include "covFnc.h"
#include <cmath>
#include <iostream>

using namespace Eigen;

// inline functions for matern32_sparse_deriv1
inline float kf(float r, float a){return (1.0+a*r)*exp(-a*r);}
inline float kf1(float r, float dx,float a){return a*a*dx*exp(-a*r);}
inline float kf2(float r, float dx1, float dx2, float delta, float a){
    return a*a*(delta-a*dx1*dx2/r)*exp(-a*r);
}

// Dimension-specific implementations are required for covariances with derivatives.
// 2D
static EMatrixX matern32_sparse_deriv1_2D(EMatrixX const& x1, std::vector<float> gradflag,
                                           float scale_param, EVectorX const& sigx, EVectorX const& siggrad);
static EMatrixX matern32_sparse_deriv1_2D(EMatrixX const& x1, std::vector<float> gradflag,
                                          EMatrixX const& x2, float scale_param);
// 3D
static EMatrixX matern32_sparse_deriv1_3D(EMatrixX const& x1, std::vector<float> gradflag,
                                           float scale_param, EVectorX const& sigx, EVectorX const& siggrad);
static EMatrixX matern32_sparse_deriv1_3D(EMatrixX const& x1, std::vector<float> gradflag,
                                          EMatrixX const& x2, float scale_param);

EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, float scale_param, float sigx)
{
     int dim = x1.rows();
    int n = x1.cols();
    float a = 1/scale_param;
    EMatrixX K = EMatrixX::Zero(n,n);

     for (int k=0;k<n;k++){
        for (int j=k;j<n;j++){
            if (k==j){
                K(k,k) = 1.0+sigx;
            }
            else{
                float r = (x1.col(k)-x1.col(j)).norm();
                K(k,j) = exp(-a*r);
                K(j,k) = K(k,j);
            }
        }
     }

    return K;
}

EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, float scale_param, EVectorX const& sigx)
{
    int dim = x1.rows();
    int n = x1.cols();
    float a = 1/scale_param;
    EMatrixX K = EMatrixX::Zero(n,n);

     for (int k=0;k<n;k++){
        for (int j=k;j<n;j++){
            if (k==j){
                K(k,k) = 1.0+sigx(k);
            }
            else{
                float r = (x1.col(k)-x1.col(j)).norm();
                K(k,j) = exp(-a*r);
                K(j,k) = K(k,j);
            }
        }
     }

    return K;
}

EMatrixX ornstein_uhlenbeck(EMatrixX const& x1, EMatrixX const& x2, float scale_param)
{
    int dim = x1.rows();
    int n = x1.cols();
    int m = x2.cols();
    float a = 1/scale_param;
    EMatrixX K = EMatrixX::Zero(n,m);

     for (int k=0;k<n;k++){
        for (int j=0;j<m;j++){
            float r = (x1.col(k)-x2.col(j)).norm();
            K(k,j) = exp(-a*r);
        }
     }

    return K;
}

EMatrixX matern32_sparse_deriv1(EMatrixX const& x1, std::vector<float> gradflag,
                                           float scale_param, EVectorX const& sigx, EVectorX const& siggrad)
{
    int dim = x1.rows();

    EMatrixX K;

    if (dim==2)
        K = matern32_sparse_deriv1_2D(x1,gradflag,scale_param,sigx,siggrad);
    else if (dim==3)
        K = matern32_sparse_deriv1_3D(x1,gradflag,scale_param,sigx,siggrad);

    return K;
}

EMatrixX matern32_sparse_deriv1(EMatrixX const& x1, std::vector<float> gradflag,
                                          EMatrixX const& x2, float scale_param)
{
    int dim = x1.rows();

    EMatrixX K;

    if (dim==2)
        K = matern32_sparse_deriv1_2D(x1,gradflag,x2,scale_param);
    else if (dim==3)
        K = matern32_sparse_deriv1_3D(x1,gradflag,x2,scale_param);

    return K;
}

// 3D
EMatrixX matern32_sparse_deriv1_3D(EMatrixX const& x1, std::vector<float> gradflag,
                                           float scale_param, EVectorX const& sigx, EVectorX const& siggrad)
{
     int dim = x1.rows();
    int n = x1.cols();
    float sqr3L = sqrt(3)/scale_param;
    float sqr3L2 = sqr3L*sqr3L;
    EMatrixX K;

    int ng = 0;
    for (auto it = gradflag.begin();it!=gradflag.end();it++){
        if ((*it) > 0.5){
            (*it) = ng;
            ng++;
        }
        else
        {
            (*it) = -1.0;
        }
    }

    K = EMatrixX::Zero(n+ng*dim,n+ng*dim);

    for (int k=0;k<n;k++){
        int kind1=gradflag[k]+n;
        int kind2 = kind1+ng;
        int kind3 = kind2+ng;

        for (int j=k;j<n;j++){
            if (k==j){
                K(k,k) = 1.0+sigx(k);
                if (gradflag[k] > -0.5){
                    K(k,kind1) = 0.0;
                    K(kind1,k) = 0.0;
                    K(k,kind2) = 0.0;
                    K(kind2,k) = 0.0;
                    K(k,kind3) = 0.0;
                    K(kind3,k) = 0.0;

                    K(kind1,kind1) = sqr3L2+siggrad(k);;
                    K(kind1,kind2) = 0.0;
                    K(kind1,kind3) = 0.0;
                    K(kind2,kind1) = 0.0;
                    K(kind2,kind2) = sqr3L2+siggrad(k);
                    K(kind2,kind3) = 0.0;
                    K(kind3,kind1) = 0.0;
                    K(kind3,kind2) = 0.0;
                    K(kind3,kind3) = sqr3L2+siggrad(k);
                }
            }
            else{
                float r = (x1.col(k)-x1.col(j)).norm();
                K(k,j) = kf(r,sqr3L);
                K(j,k) = K(k,j);
                if (gradflag[k] > -1){

                    K(kind1,j) = -kf1(r,x1(0,k)-x1(0,j),sqr3L);
                    K(j,kind1) = K(kind1,j);
                    K(kind2,j) = -kf1(r,x1(1,k)-x1(1,j),sqr3L);
                    K(j,kind2) = K(kind2,j);
                    K(kind3,j) = -kf1(r,x1(2,k)-x1(2,j),sqr3L);
                    K(j,kind3) = K(kind3,j);

                    if (gradflag[j] > -1){

                        int jind1=gradflag[j]+n;
                        int jind2 = jind1+ng;
                        int jind3 = jind2+ng;
                        K(k,jind1) = -K(j,kind1);
                        K(jind1,k) =  K(k,jind1);
                        K(k,jind2) = -K(j,kind2);
                        K(jind2,k) =  K(k,jind2);
                        K(k,jind3) = -K(j,kind3);
                        K(jind3,k) =  K(k,jind3);

                        K(kind1,jind1) = kf2(r,x1(0,k)-x1(0,j),x1(0,k)-x1(0,j),1.0,sqr3L);
                        K(jind1,kind1) = K(kind1,jind1);
                        K(kind1,jind2) = kf2(r,x1(0,k)-x1(0,j),x1(1,k)-x1(1,j),0.0,sqr3L);
                        K(jind1,kind2) = K(kind1,jind2);
                        K(kind1,jind3) = kf2(r,x1(0,k)-x1(0,j),x1(2,k)-x1(2,j),0.0,sqr3L);
                        K(jind1,kind3) = K(kind1,jind3);

                        K(kind2,jind1) = K(kind1,jind2);
                        K(jind2,kind1) = K(kind1,jind2);
                        K(kind2,jind2) = kf2(r,x1(1,k)-x1(1,j),x1(1,k)-x1(1,j),1.0,sqr3L);
                        K(jind2,kind2) = K(kind2,jind2);
                        K(kind2,jind3) = kf2(r,x1(1,k)-x1(1,j),x1(2,k)-x1(2,j),0.0,sqr3L);
                        K(jind2,kind3) = K(kind2,jind3);

                        K(kind3,jind1) = K(kind1,jind3);
                        K(jind3,kind1) = K(kind1,jind3);
                        K(kind3,jind2) = K(kind2,jind3);
                        K(jind3,kind2) = K(kind2,jind3);
                        K(kind3,jind3) = kf2(r,x1(2,k)-x1(2,j),x1(2,k)-x1(2,j),1.0,sqr3L);
                        K(jind3,kind3) = K(kind3,jind3);
                    }
                }
                else if (gradflag[j] > -1){

                    int jind1=gradflag[j]+n;
                    int jind2 = jind1+ng;
                    int jind3 = jind2+ng;
                    K(k,jind1) = kf1(r,x1(0,k)-x1(0,j),sqr3L);
                    K(jind1,k) = K(k,jind1);
                    K(k,jind2) = kf1(r,x1(1,k)-x1(1,j),sqr3L);
                    K(jind2,k) = K(k,jind2);
                    K(k,jind3) = kf1(r,x1(2,k)-x1(2,j),sqr3L);
                    K(jind3,k) = K(k,jind3);
                }
            }
        }
    }

    return K;
}

EMatrixX matern32_sparse_deriv1_3D(EMatrixX const& x1, std::vector<float> gradflag,
                                          EMatrixX const& x2, float scale_param)
{
    int dim = x1.rows();
    int n = x1.cols();
    float sqr3L = sqrt(3)/scale_param;
    float sqr3L2 = sqr3L*sqr3L;
    EMatrixX K;

    int ng = 0;
    for (auto it = gradflag.begin();it!=gradflag.end();it++){
        if ((*it) > 0.5){
            (*it) = ng;
            ng++;
        }
        else
        {
            (*it) = -1.0;
        }
    }

    int m = x2.cols();
    int m2 = m+m;
    int m3 = m2+m;

    K = EMatrixX::Zero(n+ng*dim,m*(1+dim));

    for (int k=0;k<n;k++){
        int kind1=gradflag[k]+n;
        int kind2 = kind1+ng;
        int kind3 = kind2+ng;
        for (int j=0;j<m;j++){
            float r = (x1.col(k)-x2.col(j)).norm();

            K(k,j) = kf(r,sqr3L);
            K(k,j+m) = kf1(r,x1(0,k)-x2(0,j),sqr3L);
            K(k,j+m2) = kf1(r,x1(1,k)-x2(1,j),sqr3L);
            K(k,j+m3) = kf1(r,x1(2,k)-x2(2,j),sqr3L);
            if (gradflag[k] > -0.5){
                K(kind1,j) = -K(k,j+m);
                K(kind2,j) = -K(k,j+m2);
                K(kind3,j) = -K(k,j+m3);
                K(kind1,j+m) = kf2(r,x1(0,k)-x2(0,j),x1(0,k)-x2(0,j),1.0,sqr3L);
                K(kind1,j+m2) =  kf2(r,x1(0,k)-x2(0,j),x1(1,k)-x2(1,j),0.0,sqr3L);
                K(kind1,j+m3) =  kf2(r,x1(0,k)-x2(0,j),x1(2,k)-x2(2,j),0.0,sqr3L);
                K(kind2,j+m) = K(kind1,j+m2);
                K(kind2,j+m2) = kf2(r,x1(1,k)-x2(1,j),x1(1,k)-x2(1,j),1.0,sqr3L);
                K(kind2,j+m3) = kf2(r,x1(1,k)-x2(1,j),x1(2,k)-x2(2,j),0.0,sqr3L);
                K(kind3,j+m) = K(kind1,j+m3);
                K(kind3,j+m2) = K(kind2,j+m3);
                K(kind3,j+m3) = kf2(r,x1(2,k)-x2(2,j),x1(2,k)-x2(2,j),1.0,sqr3L);
            }
        }
    }

    return K;
}

// 2D
EMatrixX matern32_sparse_deriv1_2D(EMatrixX const& x1,std::vector<float> gradflag, float scale_param,
                                EVectorX const& sigx,EVectorX const& siggrad)
{
    int dim = x1.rows();
    int n = x1.cols();
    float sqr3L = sqrt(3)/scale_param;
    float sqr3L2 = sqr3L*sqr3L;
    EMatrixX K;

    int ng = 0;
    for (auto it = gradflag.begin();it!=gradflag.end();it++){
        if ((*it) > 0.5){
            (*it) = ng;
            ng++;
        }
        else
        {
            (*it) = -1.0;
        }
    }

    K = EMatrixX::Zero(n+ng*dim,n+ng*dim);

    for (int k=0;k<n;k++){
        int kind1=gradflag[k]+n;
        int kind2 = kind1+ng;

        for (int j=k;j<n;j++){
            if (k==j){
                K(k,k) = 1.0+sigx(k);
                if (gradflag[k] > -0.5){
                    K(k,kind1) = 0.0;
                    K(kind1,k) = 0.0;
                    K(k,kind2) = 0.0;
                    K(kind2,k) = 0.0;
                    K(kind1,kind1) = sqr3L2+sqrt(sigx(k)*siggrad(k));
                    K(kind1,kind2) = 0.0;
                    K(kind2,kind1) = 0.0;
                    K(kind2,kind2) = sqr3L2+siggrad(k);
                }
            }
            else{
                float r = (x1.col(k)-x1.col(j)).norm();
                K(k,j) = kf(r,sqr3L);
                K(j,k) = K(k,j);
                if (gradflag[k] > -1){

                    K(kind1,j) = -kf1(r,x1(0,k)-x1(0,j),sqr3L);
                    K(j,kind1) = K(kind1,j);
                    K(kind2,j) = -kf1(r,x1(1,k)-x1(1,j),sqr3L);
                    K(j,kind2) = K(kind2,j);

                    if (gradflag[j] > -1){

                        int jind1=gradflag[j]+n;
                        int jind2 = jind1+ng;
                        K(k,jind1) = -K(j,kind1);
                        K(jind1,k) =  K(k,jind1);
                        K(k,jind2) = -K(j,kind2);
                        K(jind2,k) =  K(k,jind2);

                        K(kind1,jind1) = kf2(r,x1(0,k)-x1(0,j),x1(0,k)-x1(0,j),1.0,sqr3L);
                        K(jind1,kind1) = K(kind1,jind1);
                        K(kind1,jind2) = kf2(r,x1(0,k)-x1(0,j),x1(1,k)-x1(1,j),0.0,sqr3L);
                        K(jind1,kind2) = K(kind1,jind2);
                        K(kind2,jind1) = K(kind1,jind2);
                        K(jind2,kind1) = K(kind1,jind2);
                        K(kind2,jind2) = kf2(r,x1(1,k)-x1(1,j),x1(1,k)-x1(1,j),1.0,sqr3L);
                        K(jind2,kind2) = K(kind2,jind2);
                    }
                }
                else if (gradflag[j] > -1){

                    int jind1=gradflag[j]+n;
                    int jind2 = jind1+ng;
                    K(k,jind1) = kf1(r,x1(0,k)-x1(0,j),sqr3L);
                    K(jind1,k) = K(k,jind1);
                    K(k,jind2) = kf1(r,x1(1,k)-x1(1,j),sqr3L);
                    K(jind2,k) = K(k,jind2);
                }
            }
        }
    }

    return K;
}

EMatrixX matern32_sparse_deriv1_2D(EMatrixX const& x1,std::vector<float> gradflag,
                                EMatrixX const& x2, float scale_param)
{
    int dim = x1.rows();
    int n = x1.cols();
    float sqr3L = sqrt(3)/scale_param;
    float sqr3L2 = sqr3L*sqr3L;
    EMatrixX K;

    int ng = 0;
    for (auto it = gradflag.begin();it!=gradflag.end();it++){
        if ((*it) > 0.5){
            (*it) = ng;
            ng++;
        }
        else
        {
            (*it) = -1.0;
        }
    }

    int m = x2.cols();
    int m2 = m*2;

    K = EMatrixX::Zero(n+ng*dim,m*(1+dim));
    for (int k=0;k<n;k++){
        int kind1=gradflag[k]+n;
        int kind2 = kind1+ng;
        for (int j=0;j<m;j++){
            float r = (x1.col(k)-x2.col(j)).norm();

            K(k,j) = kf(r,sqr3L);
            K(k,j+m) = kf1(r,x1(0,k)-x2(0,j),sqr3L);
            K(k,j+m2) = kf1(r,x1(1,k)-x2(1,j),sqr3L);
            if (gradflag[k] > -0.5){
                K(kind1,j) = -K(k,j+m);
                K(kind2,j) = -K(k,j+m2);
                K(kind1,j+m) = kf2(r,x1(0,k)-x2(0,j),x1(0,k)-x2(0,j),1.0,sqr3L);
                K(kind1,j+m2) = kf2(r,x1(0,k)-x2(0,j),x1(1,k)-x2(1,j),0.0,sqr3L);
                K(kind2,j+m) =  K(kind1,j+m2);
                K(kind2,j+m2) = kf2(r,x1(1,k)-x2(1,j),x1(1,k)-x2(1,j),1.0,sqr3L);
            }
        }
    }

    return K;
}
