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
 */
////////////////////////////////////////////////////
// Observation GP Module
//
// Bhoram Lee (bhoram.lee@gmail.com)

#ifndef __ObsGP_H__
#define __ObsGP_H__

#include <vector>
#include <memory>
#include <iostream>
#include <cstdint>
#include <Eigen/Dense>
#include "strct.h"    // FLOAT and __USE_DOUBLE_PRECISION__ is defined in this file
                      // Currently, "typedef double FLOAT"

#define DEFAULT_SCALE_PARAM 0.5
#define DEFAULT_NOISE_PARAM 0.01

// 1D
#define DEFAULT_OVERLAP_SZ  6
#define DEFAULT_GROUP_SZ    20
#define DEFAULT_MARGIN      0.0175

// 2D
#define DEFAULT_MARGIN2      0.005
#define DEFAULT_OVERLAP_SZ2  3
#define DEFAULT_GROUP_SZ2    5

#ifdef __USE_DOUBLE_PRECISION__     // This is in use
typedef Eigen::MatrixXd EMatrixX;
typedef Eigen::VectorXd EVectorX;
typedef Eigen::RowVectorXd ERowVectorX;
#else                               // For potential transition to float
typedef Eigen::MatrixXf EMatrixX;
typedef Eigen::VectorXf EVectorX;
typedef Eigen::RowVectorXf ERowVectorX;
#endif

// This class builds a GP regressor using the Ornstein-Uuhlenbeck covariance function.
// NOTE: See covFnc.h)
class GPou{
    EMatrixX x;
    EMatrixX L;
    EVectorX alpha;

    int dim; // need this?
    const FLOAT scale = DEFAULT_SCALE_PARAM;
    const FLOAT noise = DEFAULT_NOISE_PARAM;
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

  /* :REMARKS:Tue May 29 15:02:07 EDT 2018:huangzonghao:
   *  Moving the param def to the derived class, since we want the params to be
   *  determined at compile time and they each has different params
   */
    // obsGPparam param;   // defined in strct.h

    bool trained;

    std::vector<std::shared_ptr<GPou> > gps;          // pointer to the local GPs

public:
    ObsGP(){}

    bool isTrained(){return trained;}
    void reset();

    virtual void train( FLOAT xt[],  FLOAT f[], int N[]) = 0;
    virtual void test(const EMatrixX& xt,EVectorX& val, EVectorX& var) = 0;
    virtual void test_mt(const EMatrixX& xt,EVectorX& val, EVectorX& var) = 0;
};

// This class implements ObsGP for 1D input.
class ObsGP1D : public ObsGP{
    int nGroup;         // number of local GPs
    int nSamples;       // number of total input points.

    std::vector<FLOAT> range;   // partitioned range for test

    const obsGPparam param = {DEFAULT_SCALE_PARAM,
                              DEFAULT_NOISE_PARAM,
                              DEFAULT_MARGIN,
                              DEFAULT_OVERLAP_SZ,
                              DEFAULT_GROUP_SZ};

public:
    ObsGP1D():nSamples(0){}

    void reset();

    // NOTE: In 1D, it must be f > 0.
    void train( FLOAT xt[],  FLOAT f[], int N[]) override;
    void test(const EMatrixX& xt,EVectorX& val, EVectorX& var) override;
    void test_mt(const EMatrixX& xt,EVectorX& val, EVectorX& var) override;

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
    std::vector<FLOAT>  Val_i;
    std::vector<FLOAT>  Val_j;

    void clearGPs();
    void computePartition(FLOAT val[], int ni, int nj);
    void trainValidPoints(FLOAT xt[], FLOAT f[]);

    const obsGPparam param = {DEFAULT_SCALE_PARAM,
                              DEFAULT_NOISE_PARAM,
                              DEFAULT_MARGIN2,
                              DEFAULT_OVERLAP_SZ2,
                              DEFAULT_GROUP_SZ2};
public:

    void reset();

    void getNumValidPoints(std::vector<int> &nPts);

    // NOTE: In 2D, the input xt must be a regular 2D array of size N[0] x N[1].
    //       If not f > 0, the point is considered invalid.
    void train( FLOAT xt[],  FLOAT f[], int N[]) override;
    void train( FLOAT xt[],  FLOAT f[], int N[], std::vector<int> &numSamples);
    void test(const EMatrixX& xt,EVectorX& val, EVectorX& var) override;
    void test_mt(const EMatrixX& xt,EVectorX& val, EVectorX& var) override;

private:
    void test_kernel(int thread_idx,
                     int start_idx,
                     int end_idx,
                     const EMatrixX &xt,
                     EVectorX &val,
                     EVectorX &var);
};

#endif
