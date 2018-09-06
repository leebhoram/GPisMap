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
 * Authors: Bhoram Lee <bhoram.lee@gmail.com>
 *          Huang Zonghao<ac@hzh.io>
 */

#include "OnGPIS.h"
#include "covFnc.h"
#include <Eigen/Cholesky>

#define SQRT_3  1.732051

void OnGPIS::reset(){
    nSamples = 0;
    trained = false;

    return;
}

void OnGPIS::train(const vecNode& samples){
    reset();

    int N = samples.size();
    int dim = 2;

    if (N > 0){
        nSamples = N;
        x = EMatrixX::Zero(dim,N);
        EMatrixX grad= EMatrixX::Zero(dim,N);
        EVectorX f = EVectorX::Zero(N);
        EVectorX sigx = EVectorX::Zero(N);
        EVectorX siggrad = EVectorX::Zero(N);

        gradflag.clear();
        gradflag.resize(N,0.0);

        EMatrixX grad_valid(N,2);

        int k=0;
        int count = 0;
        for (auto it = samples.begin(); it!=samples.end(); it++, k++){
            x(0,k) = (*it)->getPosX();
            x(1,k) = (*it)->getPosY();
            grad(0,k) = (*it)->getGradX();
            grad(1,k) = (*it)->getGradY();
            f(k) = (*it)->getVal();
            sigx(k) = (*it)->getPosNoise();
            siggrad(k) = (*it)->getGradNoise();
            if (siggrad(k) > 0.1001 || (fabs(grad(0,k)) < 1e-6 && fabs(grad(1,k)) < 1e-6)){
                gradflag[k] = 0.0;
                sigx(k) = 2.0;
            }
            else{
                gradflag[k] = 1.0;
                grad_valid(count,0) = grad(0,k);
                grad_valid(count,1) = grad(1,k) ;
                count++;
            }
        }
        grad_valid.conservativeResize(count,2);
        EVectorX y(N+2*count);
        y << f,grad_valid.col(0),grad_valid.col(1);
        EMatrixX K = matern32_sparse_deriv1(x, gradflag, param.scale, sigx, siggrad);

        L = K.llt().matrixL();

        alpha = y;
        L.template triangularView<Eigen::Lower>().solveInPlace(alpha);
        L.transpose().template triangularView<Eigen::Upper>().solveInPlace(alpha);

        trained = true;

    }
    return;
}

void OnGPIS::train(const vecNode3& samples){
    reset();

    int N = samples.size();
    int dim = 3;

    if (N > 0){
        nSamples = N;
        x = EMatrixX::Zero(dim,N);
        EMatrixX grad= EMatrixX::Zero(dim,N);
        EVectorX f = EVectorX::Zero(N);
        EVectorX sigx = EVectorX::Zero(N);
        EVectorX siggrad = EVectorX::Zero(N);

        gradflag.clear();
        gradflag.resize(N,0.0);

        EMatrixX grad_valid(N,dim);

        int k=0;
        int count = 0;
        for (auto it = samples.begin(); it!=samples.end(); it++, k++){
            x(0,k) = (*it)->getPosX();
            x(1,k) = (*it)->getPosY();
            x(2,k) = (*it)->getPosZ();
            grad(0,k) = (*it)->getGradX();
            grad(1,k) = (*it)->getGradY();
            grad(2,k) = (*it)->getGradZ();
            f(k) = (*it)->getVal();
            sigx(k) = (*it)->getPosNoise();
            siggrad(k) = (*it)->getGradNoise();
            if (siggrad(k) > 0.1001 || (fabs(grad(0,k)) < 1e-6 && fabs(grad(1,k)) < 1e-6 && fabs(grad(2,k)) < 1e-6)){
                gradflag[k] = 0.0;
                sigx(k) = 2.0;
            }
            else{
                gradflag[k] = 1.0;
                grad_valid(count,0) = grad(0,k);
                grad_valid(count,1) = grad(1,k) ;
                grad_valid(count,2) = grad(2,k) ;
                count++;
            }
        }
        grad_valid.conservativeResize(count,3);
        EVectorX y(N+dim*count);
        y << f,grad_valid.col(0),grad_valid.col(1),grad_valid.col(2);
        EMatrixX K = matern32_sparse_deriv1(x, gradflag, param.scale, sigx, siggrad);

        L = K.llt().matrixL();

        alpha = y;
        L.template triangularView<Eigen::Lower>().solveInPlace(alpha);
        L.transpose().template triangularView<Eigen::Upper>().solveInPlace(alpha);

        trained = true;

    }
    return;
}


void OnGPIS::test(const EMatrixX& xt,EVectorX& val, EMatrixX& gradval, EVectorX& var){

    if (!isTrained()){
        return;
    }

    EMatrixX K = matern32_sparse_deriv1(x, gradflag, xt, param.scale);

    EVectorX res = K.transpose()*alpha;
    val = res.head(xt.cols());
    gradval.row(0) =  res.segment(xt.cols(),xt.cols()).transpose();
    gradval.row(1) =  res.segment(2*xt.cols(),xt.cols()).transpose();

    L.template triangularView<Eigen::Lower>().solveInPlace(K);

    K = K.array().pow(2);
    EVectorX v = K.colwise().sum();

    var.row(0) = 1.01-v.head(xt.cols()).array();
    var.row(1) = three_over_scale + 0.1 - v.segment(xt.cols(),xt.cols()).array();
    var.row(2) = three_over_scale + 0.1 - v.segment(2*xt.cols(),xt.cols()).array();

    return;
}

void OnGPIS::testSinglePoint(const EVectorX& xt, float& val, float grad[],float var[])
{
    if (!isTrained())
        return;

    if (x.rows() != xt.size())
        return;

    EMatrixX K = matern32_sparse_deriv1(x, gradflag, xt, param.scale);

    EVectorX res = K.transpose()*alpha;
    val = res(0);
    if (res.size() == 3){
        grad[0] = res(1);
        grad[1] = res(2);
    }
    else if (res.size() == 4){
        grad[0] = res(1);
        grad[1] = res(2);
        grad[2] = res(3);
    }

    L.template triangularView<Eigen::Lower>().solveInPlace(K);
    K = K.array().pow(2);
    EVectorX v = K.colwise().sum();

    if (v.size() == 3){
        var[0] = 1.01-v(0);
        var[1] = three_over_scale + 0.1 - v(1);
        var[2] = three_over_scale + 0.1 - v(2);
    }
    else if (v.size() == 4){ // Noise param!
        var[0] = 1.001-v(0);
        var[1] = three_over_scale + 0.001 - v(1);
        var[2] = three_over_scale + 0.001 - v(2);
        var[3] = three_over_scale + 0.001 - v(3);
    }

    return;
}

void OnGPIS::test2Dpoint(const EVectorX& xt,float& val, float& gradx, float& grady, float& varval ,float& vargradx, float &vargrady)
{
    if (!isTrained()){
        return;
    }

    EMatrixX K = matern32_sparse_deriv1(x, gradflag, xt, param.scale);

    EVectorX res = K.transpose()*alpha;
    val = res(0);
    gradx = res(1);
    grady = res(2);

    L.template triangularView<Eigen::Lower>().solveInPlace(K);
    K = K.array().pow(2);
    EVectorX v = K.colwise().sum();

    varval = 1.01-v(0);
    vargradx = three_over_scale + 0.1 - v(1);
    vargrady = three_over_scale + 0.1 - v(2);

}

void OnGPIS::test2Dpoint(float px, float py,float& val, float& gradx, float& grady, float& varval ,float& vargradx, float &vargrady)
{
    if (!isTrained()){
        return;
    }

    EVectorX xt(2);
    xt << px,py;
    EMatrixX K = matern32_sparse_deriv1(x, gradflag, xt, param.scale);

    EVectorX res = K.transpose()*alpha;
    val = res(0);
    gradx = res(1);
    grady = res(2);

    L.template triangularView<Eigen::Lower>().solveInPlace(K);
    K = K.array().pow(2);
    EVectorX v = K.colwise().sum();

    varval = 1.01-v(0);
    vargradx = three_over_scale + 0.1 - v(1);
    vargrady = three_over_scale + 0.1 - v(2);
}
