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
#ifndef __GPIS_MAP3__H__
#define __GPIS_MAP3__H__

#include "ObsGP.h"
#include "OnGPIS.h"
#include "octree.h"

typedef struct camParam_{
    FLOAT fx;
    FLOAT fy;
    FLOAT cx;
    FLOAT cy;
    int width;
    int height;

    camParam_(){
        width = 640;
        height = 480;
        fx = 568.0; // 570.9361;
        fy = 568.0; // 570.9361;
        cx = 310;// 307;
        cy = 224; //240;
    }
    camParam_(FLOAT fx_, FLOAT fy_, FLOAT cx_, FLOAT cy_, FLOAT w_, FLOAT h_): fx(fx_),fy(fy_), cx(cx_), cy(cy_), width(w_), height(h_){}
} camParam;

typedef struct GPisMap3Param_{
    FLOAT delx;         // numerical step delta (e.g. surface normal sampling)
    FLOAT fbias;        // constant map bias values (mean of GP)
//    FLOAT sensor_offset[3];
    FLOAT obs_var_thre; // threshold for variance of ObsGP
                        //  - If var(prediction) > v_thre, then don't rely on the prediction.
    int obs_skip;     // use every 'skip'-th pixel
    FLOAT min_position_noise;
    FLOAT min_grad_noise;

    FLOAT map_scale_param;
    FLOAT map_noise_param;

    GPisMap3Param_(){
        delx = 1e-3;
        fbias = 0.2;
        obs_skip = 2;
        obs_var_thre = 0.04;
//         sensor_offset[0] = 0.0;
//         sensor_offset[1] = 0.0;
//         sensor_offset[2] = 0.0;
        min_position_noise = 1e-3;
        min_grad_noise = 1e-2;
        map_scale_param = 0.04; // 0.1;
        map_noise_param = 5e-3;
    }

    GPisMap3Param_(GPisMap3Param_& par){
        delx = par.delx;
        fbias = par.fbias;
        obs_skip = par.obs_skip;
        obs_var_thre = par.obs_var_thre;
//         sensor_offset[0] = par.sensor_offset[0];
//         sensor_offset[1] = par.sensor_offset[1];
//         sensor_offset[2] = par.sensor_offset[2];
        min_position_noise = par.min_position_noise;
        min_grad_noise = par.min_grad_noise;
        map_scale_param = par.map_scale_param;
        map_noise_param = par.map_noise_param;
        //useCuda = par.useCuda;
    }
}GPisMap3Param;

class GPisMap3{
protected:
    GPisMap3Param setting;
    camParam cam;

    FLOAT u_obs_limit[2];
    FLOAT v_obs_limit[2];

    std::vector<FLOAT> vu_grid;

    OcTree* t;
    std::unordered_set<OcTree*> activeSet;
    const int mapDimension = 3;
    FLOAT runtime[4];

    void init();
    bool preprocData( FLOAT * dataz, int N, std::vector<FLOAT> & pose);
    bool regressObs();
    void updateMapPoints();
    void reEvalPoints(std::vector<std::shared_ptr<Node3> >& nodes);
    void evalPoints();
    void addNewMeas();
    void updateGPs();
    void updateGPs_mt();

    ObsGP* gpo;
    std::vector<FLOAT> obs_valid_u;
    std::vector<FLOAT> obs_valid_v;
    std::vector<FLOAT> obs_zinv;
    std::vector<FLOAT> obs_valid_xyzlocal;
    std::vector<FLOAT> obs_valid_xyzglobal;
    std::vector<FLOAT> pose_tr;
    std::vector<FLOAT> pose_R;
    int obs_numdata;
    FLOAT range_obs_max;

public:
    GPisMap3();
    GPisMap3(GPisMap3Param par);
    GPisMap3(GPisMap3Param par, camParam c);
    ~GPisMap3();
    void reset();

    void update( FLOAT * dataz, int N, std::vector<FLOAT> & pose);
    void update_mt( FLOAT * dataz, int N, std::vector<FLOAT> & pose);
    bool test( FLOAT* x, int dim, int leng, FLOAT * res);
    bool test_mt( FLOAT* x, int dim, int leng, FLOAT * res);
    void resetCam(camParam c);

    void getAllPoints(std::vector<FLOAT> & pos);
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
                          OcTree **nodes_to_update);
};
#endif
