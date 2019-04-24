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

#ifndef PARAMS_H_
#define PARAMS_H_

/*-----------------------------------------------------------------------------
 *  Tree Parameters
 *-----------------------------------------------------------------------------*/
#define DEFAULT_TREE_INIT_ROOT_HALFLENGTH 12.8  // 0.05*2^7
#define DEFAULT_TREE_MIN_HALFLENGTH       0.05  // 0.05*2^0
#define DEFAULT_TREE_MAX_HALFLENGTH       102.4 // 0.05*2^10
#define DEFAULT_TREE_CLUSTER_HALFLENGTH   0.1   // 0.05*2^2

/* GPisMap */
// Note: the parameters must be (some common base)^n (n: integer)
#define GPISMAP_TREE_CLUSTER_HALF_LENGTH   (0.8)
#define GPISMAP_TREE_MIN_HALF_LENGTH       (0.2)
#define GPISMAP_TREE_MAX_HALF_LENGTH       (102.4)
#define GPISMAP_TREE_INIT_ROOT_HALF_LENGTH (12.8)

/* GPisMap3 */
#define GPISMAP3_RTIMES                     2.0
#define GPISMAP3_TREE_CLUSTER_HALF_LENGTH   (0.025)
#define GPISMAP3_TREE_MIN_HALF_LENGTH       (0.0125/2.0)
#define GPISMAP3_TREE_MAX_HALF_LENGTH       (1.6)
#define GPISMAP3_TREE_INIT_ROOT_HALF_LENGTH (0.4)
// Note: 1.6 = 0.0125*(2.0^7)
//       0.4 = 0.0125*(2.0^5)

/*-----------------------------------------------------------------------------
 *  GP Map Params
 *-----------------------------------------------------------------------------*/

/* Mapping GP */
#define DEFAULT_MAP_SCALE_PARAM 1 //0.5
#define DEFAULT_MAP_NOISE_PARAM 1e-2

/* GPisMap */

/* delx: numerical step delta (e.g. surface normal sampling)
 * fbias: constant map bias values (mean of GP)
 * obs_var_thre: threshold for variance of ObsGP
 *                - If var(prediction) > v_thre, then don't rely on the prediction.
 * obs_skip: use every 'skip'-th pixel
 */
#define GPISMAP_DELX              1e-2
#define GPISMAP_FBIAS             0.2
#define GPISMAP_OBS_VAR_THRE      0.1
#define GPISMAP_SENSOR_OFFSET_0   0.08  // the hokuyo sensor position (0.08, 0) on a turtlebot for simulation
#define GPISMAP_SENSOR_OFFSET_1   0.0   // 
#define GPISMAP_ANGLE_OBS_LIMIT_0 (-135.0*M_PI/180.0)
#define GPISMAP_ANGLE_OBS_LIMIT_1 (135.0*M_PI/180.0)
#define GPISMAP_MIN_POS_NOISE     1e-2
#define GPISMAP_MIN_GRAD_NOISE    1e-2
#define GPISMAP_MAP_SCALE         1.2
#define GPISMAP_MAP_NOISE         1e-2

/* GPisMap3 */
#define GPISMAP3_MAX_RANGE   4e0
#define GPISMAP3_MIN_RANGE   4e-1

/* delx: numerical step delta (e.g. surface normal sampling)
 * fbias: constant map bias values (mean of GP)
 * obs_var_thre: threshold for variance of ObsGP
 *                - If var(prediction) > v_thre, then don't rely on the prediction.
 * obs_skip: use every 'skip'-th pixel
 */
#define GPISMAP3_DELX           1e-3
#define GPISMAP3_FBIAS          0.2
#define GPISMAP3_OBS_SKIP       2
#define GPISMAP3_OBS_VAR_THRE   0.04
#define GPISMAP3_MIN_POS_NOISE  1e-3
#define GPISMAP3_MIN_GRAD_NOISE 1e-2
#define GPISMAP3_MAP_SCALE      0.04
#define GPISMAP3_MAP_NOISE      5e-3

/*-----------------------------------------------------------------------------
 *  ObsGP Params
 *-----------------------------------------------------------------------------*/

#define DEFAULT_OBSGP_SCALE_PARAM 0.5
#define DEFAULT_OBSGP_NOISE_PARAM 0.01

/* 1D */
#define DEFAULT_OBSGP_OVERLAP_SZ  6
#define DEFAULT_OBSGP_GROUP_SZ    20
#define DEFAULT_OBSGP_MARGIN      0.0175

/* 2D */
#define DEFAULT_OBSGP_MARGIN2      0.005
#define DEFAULT_OBSGP_OVERLAP_SZ2  3
#define DEFAULT_OBSGP_GROUP_SZ2    5

#endif   /* ----- #ifndef PARAMS_H_  ----- */
