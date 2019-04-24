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

#ifndef __GPIS_MAP3__H__
#define __GPIS_MAP3__H__

#include "ObsGP.h"
#include "OnGPIS.h"
#include "octree.h"
#include "params.h"

typedef struct camParam_{
    float fx;
    float fy;
    float cx;
    float cy;
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
    camParam_(float fx_, float fy_, float cx_, float cy_, float w_, float h_): fx(fx_),fy(fy_), cx(cx_), cy(cy_), width(w_), height(h_){}
} camParam;

typedef struct GPisMap3Param_{
    float delx;         // numerical step delta (e.g. surface normal sampling)
    float fbias;        // constant map bias values (mean of GP)
    float obs_var_thre; // threshold for variance of ObsGP
                        //  - If var(prediction) > v_thre, then don't rely on the prediction.
    int   obs_skip;     // use every 'skip'-th pixel
    float min_position_noise;
    float min_grad_noise;

    float map_scale_param;
    float map_noise_param;

    GPisMap3Param_(){
        delx = GPISMAP3_DELX;
        fbias = GPISMAP3_FBIAS;
        obs_skip = GPISMAP3_OBS_SKIP;
        obs_var_thre = GPISMAP3_OBS_VAR_THRE;
        min_position_noise = GPISMAP3_MIN_POS_NOISE;
        min_grad_noise = GPISMAP3_MIN_GRAD_NOISE;
        map_scale_param = GPISMAP3_MAP_SCALE;
        map_noise_param = GPISMAP3_MAP_NOISE;
    }

    GPisMap3Param_(GPisMap3Param_& par){
        delx = par.delx;
        fbias = par.fbias;
        obs_skip = par.obs_skip;
        obs_var_thre = par.obs_var_thre;
        min_position_noise = par.min_position_noise;
        min_grad_noise = par.min_grad_noise;
        map_scale_param = par.map_scale_param;
        map_noise_param = par.map_noise_param;
    }
}GPisMap3Param;

class GPisMap3{
protected:
    GPisMap3Param setting;
    camParam cam;

    float u_obs_limit[2];
    float v_obs_limit[2];

    std::vector<float> vu_grid;

    OcTree* t;
    std::unordered_set<OcTree*> activeSet;
    const int mapDimension = 3;
 
    void init();
    bool preprocData( float * dataz, int N, std::vector<float> & pose);
    bool regressObs();
    void updateMapPoints();
    void reEvalPoints(std::vector<std::shared_ptr<Node3> >& nodes);
    void evalPoints();
    void addNewMeas();
    void updateGPs();

    ObsGP* gpo;
    std::vector<float> obs_valid_u;
    std::vector<float> obs_valid_v;
    std::vector<float> obs_zinv;
    std::vector<float> obs_valid_xyzlocal;
    std::vector<float> obs_valid_xyzglobal;
    std::vector<float> pose_tr;
    std::vector<float> pose_R;
    int obs_numdata;
    float range_obs_max;

public:
    GPisMap3();
    GPisMap3(GPisMap3Param par);
    GPisMap3(GPisMap3Param par, camParam c);
    ~GPisMap3();
    void reset();

    void getAllPoints(std::vector<float> & pos);
    void update( float * dataz, int N, std::vector<float> & pose);
    bool test( float* x, int dim, int leng, float * res);
    void resetCam(camParam c);

private:
    void test_kernel(int thread_idx,
                     int start_idx,
                     int end_idx,
                     float *x,
                     float *res);

    void updateGPs_kernel(int thread_idx,
                          int start_idx,
                          int end_idx,
                          OcTree **nodes_to_update);
};

#endif
