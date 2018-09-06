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
#ifndef __GPIS_MAP__H__
#define __GPIS_MAP__H__

#include "ObsGP.h"
#include "OnGPIS.h"
#include "quadtree.h"
#include <chrono>

typedef struct GPisMapParam_{
    FLOAT delx;         // numerical step delta (e.g. surface normal sampling)
    FLOAT fbias;        // constant map bias values (mean of GP)
    FLOAT sensor_offset[2];
    FLOAT angle_obs_limit[2];
    FLOAT obs_var_thre; // threshold for variance of ObsGP
                        //  - If var(prediction) > v_thre, then don't rely on the prediction.
    FLOAT min_position_noise;
    FLOAT min_grad_noise;

    FLOAT map_scale_param;
    FLOAT map_noise_param;

    int useCuda;

    GPisMapParam_(){
        delx = 1e-2;
        fbias = 0.2;
        obs_var_thre = 0.1;
        sensor_offset[0] = 0.08;
        sensor_offset[1] = 0.0;
        angle_obs_limit[0] = -135.0*M_PI/180.0;
        angle_obs_limit[1] = 135.0*M_PI/180.0;
        min_position_noise = 1e-2;
        min_grad_noise = 1e-2;
        map_scale_param = 1.2;
        map_noise_param = 1e-2;
        useCuda = 0;
    }

    GPisMapParam_( GPisMapParam_& par){
        delx = par.delx;
        fbias = par.fbias;
        obs_var_thre = par.obs_var_thre;
        sensor_offset[0] = par.sensor_offset[0];
        sensor_offset[1] = par.sensor_offset[1];
        min_position_noise = par.min_position_noise;
        min_grad_noise = par.min_grad_noise;
        map_scale_param = par.map_scale_param;
        map_noise_param = par.map_noise_param;
        useCuda = par.useCuda;
    }
}GPisMapParam;

class GPisMap{
protected:
    GPisMapParam setting;

    QuadTree* t;
    std::unordered_set<QuadTree*>  activeSet;
    const int mapDimension = 2;

    void init();
    bool preproData( FLOAT * datax,  FLOAT * dataf, int N, std::vector<FLOAT> & pose);
    bool regressObs();
    void updateMapPoints();
    void reEvalPoints(std::vector<std::shared_ptr<Node> >& nodes);
    void evalPoints();
    void addNewMeas();
    void updateGPs();
    void updateGPs_mt();

    ObsGP* gpo;
    std::vector<FLOAT> obs_theta;
    std::vector<FLOAT> obs_range;
    std::vector<FLOAT> obs_f;
    std::vector<FLOAT> obs_xylocal;
    std::vector<FLOAT> obs_xyglobal;
    std::vector<FLOAT> pose_tr;
    std::vector<FLOAT> pose_R;
    int obs_numdata;
    FLOAT range_obs_max;

    double runtime[4];

public:
    GPisMap();
    GPisMap(GPisMapParam par);
    ~GPisMap();
    void reset();

    void update(  FLOAT * datax,  FLOAT * dataf, int N, std::vector<FLOAT> & pose);
    void update_mt(  FLOAT * datax,  FLOAT * dataf, int N, std::vector<FLOAT> & pose);
    bool test( FLOAT* x, int dim, int leng, FLOAT * res);
    bool test_mt( FLOAT* x, int dim, int leng, FLOAT * res);

    int getMapDimension(){return mapDimension;}
    void getAllPoints(std::vector<FLOAT> & pos);
    // void getAllPoints(std::vector<FLOAT> & pos, std::vector<FLOAT> &var);
    void getAllPoints(std::vector<FLOAT> & pos, std::vector<FLOAT> &var, std::vector<FLOAT> &grad,  std::vector<FLOAT> &grad_var);

    double getRuntime0(){return runtime[0];}
    double getRuntime1(){return runtime[1];}
    double getRuntime2(){return runtime[2];}
    double getRuntime3(){return runtime[3];}

private:
    void test_kernel(int thread_idx,
                     int start_idx,
                     int end_idx,
                     FLOAT *x,
                     FLOAT *res);

    void updateGPs_kernel(int thread_idx,
                          int start_idx,
                          int end_idx,
                          QuadTree **nodes_to_update);
};

#endif
