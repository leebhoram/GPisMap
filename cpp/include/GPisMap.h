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

#ifndef __GPIS_MAP__H__
#define __GPIS_MAP__H__

#include "ObsGP.h"
#include "OnGPIS.h"
#include "quadtree.h"
#include "params.h"

typedef struct GPisMapParam_{
    float delx;         // numerical step delta (e.g. surface normal sampling)
    float fbias;        // constant map bias values (mean of GP)
    float sensor_offset[2];
    float angle_obs_limit[2];
    float obs_var_thre; // threshold for variance of ObsGP
                        //  - If var(prediction) > v_thre, then don't rely on the prediction.
    float min_position_noise;
    float min_grad_noise;

    float map_scale_param;
    float map_noise_param;

    GPisMapParam_(){
        delx = GPISMAP_DELX;
        fbias = GPISMAP_FBIAS;
        obs_var_thre = GPISMAP_OBS_VAR_THRE;
        sensor_offset[0] = GPISMAP_SENSOR_OFFSET_0;
        sensor_offset[1] = GPISMAP_SENSOR_OFFSET_1;
        angle_obs_limit[0] = GPISMAP_ANGLE_OBS_LIMIT_0;
        angle_obs_limit[1] = GPISMAP_ANGLE_OBS_LIMIT_1;
        min_position_noise = GPISMAP_MIN_POS_NOISE;
        min_grad_noise = GPISMAP_MIN_GRAD_NOISE;
        map_scale_param = GPISMAP_MAP_SCALE;
        map_noise_param = GPISMAP_MAP_NOISE;
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
    }
}GPisMapParam;

class GPisMap{
protected:
    GPisMapParam setting;

    QuadTree* t;
    std::unordered_set<QuadTree*>  activeSet;
    const int mapDimension = 2;

    void init();
    bool preproData( float * datax,  float * dataf, int N, std::vector<float> & pose);
    bool regressObs();
    void updateMapPoints();
    void reEvalPoints(std::vector<std::shared_ptr<Node> >& nodes);
    void evalPoints();
    void addNewMeas();
    void updateGPs();

    ObsGP* gpo;
    std::vector<float> obs_theta;
    std::vector<float> obs_range;
    std::vector<float> obs_f;
    std::vector<float> obs_xylocal;
    std::vector<float> obs_xyglobal;
    std::vector<float> pose_tr;
    std::vector<float> pose_R;
    int obs_numdata;
    float range_obs_max;

public:
    GPisMap();
    GPisMap(GPisMapParam par);
    ~GPisMap();
    void reset();

    void update(  float * datax,  float * dataf, int N, std::vector<float> & pose);
    bool test( float* x, int dim, int leng, float * res);

    int getMapDimension(){return mapDimension;}

private:
    void test_kernel(int thread_idx,
                     int start_idx,
                     int end_idx,
                     float *x,
                     float *res);

    void updateGPs_kernel(int thread_idx,
                          int start_idx,
                          int end_idx,
                          QuadTree **nodes_to_update);
};

#endif
