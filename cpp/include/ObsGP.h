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

#ifndef __ObsGP_H__
#define __ObsGP_H__

#include <vector>
#include <memory>
#include <iostream>
#include <cstdint>
#include <Eigen/Dense>
#include "strct.h"
#include "params.h"

typedef Eigen::MatrixXf EMatrixX;
typedef Eigen::VectorXf EVectorX;
typedef Eigen::RowVectorXf ERowVectorX;

// This class builds a GP regressor using the Ornstein-Uuhlenbeck covariance function.
// NOTE: See covFnc.h)
class GPou{
    EMatrixX x;
    EMatrixX L;
    EVectorX alpha;

    int dim; // need this?
    const float scale = DEFAULT_OBSGP_SCALE_PARAM;
    const float noise = DEFAULT_OBSGP_NOISE_PARAM;
    bool trained = false;

public:

    GPou(){}

    void reset(){trained = false;}

    bool isTrained(){return trained;}

    int  getNumSamples(){return x.cols();}
    void train(const EMatrixX& xt,const EVectorX& f);
    void test(const EMatrixX& xt, EVectorX& f, EVectorX& var);

};

// This is a base class to build a partitioned GP regressor, holding multiple local GPs using GPou.
class ObsGP{
protected:

    bool trained;

    std::vector<std::shared_ptr<GPou> > gps;          // pointer to the local GPs

public:
    ObsGP(){}

    bool isTrained(){return trained;}
    void reset();

    virtual void train( float xt[],  float f[], int N[]) = 0;
    virtual void test(const EMatrixX& xt,EVectorX& val, EVectorX& var) = 0;
};

// This class implements ObsGP for 1D input.
class ObsGP1D : public ObsGP{
    int nGroup;         // number of local GPs
    int nSamples;       // number of total input points.

    std::vector<float> range;   // partitioned range for test

    const obsGPparam param = {DEFAULT_OBSGP_SCALE_PARAM,
                              DEFAULT_OBSGP_NOISE_PARAM,
                              DEFAULT_OBSGP_MARGIN,
                              DEFAULT_OBSGP_OVERLAP_SZ,
                              DEFAULT_OBSGP_GROUP_SZ};

public:
    ObsGP1D():nSamples(0){}

    void reset();

    // NOTE: In 1D, it must be f > 0.
    void train( float xt[],  float f[], int N[]) override;
    void test(const EMatrixX& xt,EVectorX& val, EVectorX& var) override;

};

// This class implements ObsGP for regular 2D input.
class ObsGP2D : public ObsGP{
    int    nGroup[2];       // dimension of local GPs
    int    szSamples[2];    // dimension of input data
    bool   repartition;

    // pre-computed partition indices
    std::vector<int>  Ind_i0;
    std::vector<int>  Ind_i1;
    std::vector<int>  Ind_j0;
    std::vector<int>  Ind_j1;

    // pre-computed partition values
    std::vector<float>  Val_i;
    std::vector<float>  Val_j;

    void clearGPs();
    void computePartition(float val[], int ni, int nj);
    void trainValidPoints(float xt[], float f[]);

    const obsGPparam param = {DEFAULT_OBSGP_SCALE_PARAM,
                              DEFAULT_OBSGP_NOISE_PARAM,
                              DEFAULT_OBSGP_MARGIN2,
                              DEFAULT_OBSGP_OVERLAP_SZ2,
                              DEFAULT_OBSGP_GROUP_SZ2};
public:

    void reset();

    void getNumValidPoints(std::vector<int> &nPts);

    // NOTE: In 2D, the input xt must be a regular 2D array of size N[0] x N[1].
    //       If not f > 0, the point is considered invalid.
    void train( float xt[],  float f[], int N[]) override;
    void train( float xt[],  float f[], int N[], std::vector<int> &numSamples);
    void test(const EMatrixX& xt,EVectorX& val, EVectorX& var) override;

private:
    void test_kernel(int thread_idx,
                     int start_idx,
                     int end_idx,
                     const EMatrixX &xt,
                     EVectorX &val,
                     EVectorX &var);
};

#endif
