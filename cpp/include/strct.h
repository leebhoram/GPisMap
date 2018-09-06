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

#ifndef __STRCT_H_
#define __STRCT_H_

#include <vector>
#include <memory>
#include "params.h"

enum NODE_TYPE {NONE=0, HIT=1, FREE, CLUSTER};

template <typename T>
struct Point
{
    T x;
    T y;

    Point(T _x, T _y)
    {
        x = _x;
        y = _y;
    }
    Point()
    {
        x = 0;
        y = 0;
    }
};

template <typename T>
struct Point3
{
    T x;
    T y;
    T z;

    Point3(T _x, T _y, T _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    Point3()
    {
        x = 0;
        y = 0;
        z = 0;
    }
};

class Node
{
    Point<float> pos;
    Point<float> grad;
    float val;
    float pose_sig;
    float grad_sig;
    NODE_TYPE nt;

public:

    Node(Point<float> _pos, float _val, float _pose_sig, Point<float> _grad, float _grad_sig, NODE_TYPE n);
    Node(Point<float> _pos, NODE_TYPE _nt = NODE_TYPE::NONE);
    Node();

    void updateData(float _val, float _pose_sig, Point<float> _grad, float _grad_sig, NODE_TYPE n);
    void updateNoise( float _pose_sig, float _grad_sig);

    const Point<float>& getPos(){return pos;}
    const Point<float>& getGrad(){return grad;}
    float getPosX(){return pos.x;}
    float getPosY(){return pos.y;}
    float getGradX(){return grad.x;}
    float getGradY(){return grad.y;}
    float getVal(){return val;}
    float getPosNoise(){return pose_sig;};
    float getGradNoise(){return grad_sig;};
    NODE_TYPE getType(){return nt;}
};

class Node3
{
    Point3<float> pos;
    Point3<float> grad;
    float val;
    float pose_sig;
    float grad_sig;
    NODE_TYPE nt;

public:

    Node3(Point3<float> _pos, float _val, float _pose_sig, Point3<float> _grad, float _grad_sig, NODE_TYPE n = NODE_TYPE::NONE);
    Node3(Point3<float> _pos, NODE_TYPE _nt = NODE_TYPE::NONE);
    Node3();

    void updateData(float _val, float _pose_sig, Point3<float> _grad, float _grad_sig, NODE_TYPE n = NODE_TYPE::NONE);
    void updateNoise( float _pose_sig, float _grad_sig);

    const Point3<float>& getPos(){return pos;}
    const Point3<float>& getGrad(){return grad;}
    float getPosX(){return pos.x;}
    float getPosY(){return pos.y;}
    float getPosZ(){return pos.z;}
    float getGradX(){return grad.x;}
    float getGradY(){return grad.y;}
    float getGradZ(){return grad.z;}
    float getVal(){return val;}
    float getPosNoise(){return pose_sig;};
    float getGradNoise(){return grad_sig;};
    NODE_TYPE getType(){return nt;}
};

//////////////////////////////////////////////////////////////////////////
// Parameters

// Observation GP
typedef struct obsGPparam_{
    // Npte:
    // ObsGP is implemented to use the Ornstein-Uhlenbeck covariance function,
    // which has a form of k(r)=exp(-r/l) (See covFnc.h)
    float scale;            // the scale parameter l
    float noise;            // the noise parameter of the measurement
                            // currently use a constant value
                            // could be potentially modified to have heteroscedastic noise
    // Note:
    // ObsGP is implemented to have overlapping partitioned GPs.
    float margin;           // used to decide if valid range
                            // (don't use if too close to boundary
                            //  because the derivates are hard to sample)
    int   overlap;          // the overlapping parameters: number of samples to overlap
    int   group_size;       // the number of samples to group together
                            // (the actual group size will be (group_size+overlap)
    obsGPparam_(){}
    obsGPparam_(float s, float n, float m, int ov, int gsz):
                scale(s),
                noise(n),
                margin(m),
                overlap(ov),
                group_size(gsz){ }
} obsGPparam;

// GPIS (SDF)
typedef struct onGPISparam_{
    // Note:
    // OnlineGPIS is implemented to use the Matern class covariance function with (nu=2/3),
    // which has a form of k(r)=(1+sqrt(3)*r/l)exp(-sqrt(3)*r/l) (See covFnc.h)
    float scale;            // the scale parameter l
    float noise;            // the default noise parameter of the measurement
                            // currently use heteroscedastic noise acoording to a noise model   
    float noise_deriv;      // the default noise parameter of the derivative measurement
                            // currently use a noise model by numerical computation.
    onGPISparam_(){}
    onGPISparam_(float s, float n, float nd): scale(s), noise(n), noise_deriv(nd){}
} onGPISparam;

// QuadTree (2D) and OcTree (3D)
typedef struct tree_param_{
    float initroot_halfleng;
    float min_halfleng;         // minimum (leaf) resolution of tree
    float min_halfleng_sqr;
    float max_halfleng;         // maximum (root) resolution of tree
    float max_halfleng_sqr;
    float cluster_halfleng;    // the resolution of GP clusters
    float cluster_halfleng_sqr;
public:
    tree_param_():min_halfleng(DEFAULT_TREE_MIN_HALFLENGTH),
                min_halfleng_sqr(DEFAULT_TREE_MIN_HALFLENGTH*DEFAULT_TREE_MIN_HALFLENGTH),
                max_halfleng(DEFAULT_TREE_MAX_HALFLENGTH),
                max_halfleng_sqr(DEFAULT_TREE_MAX_HALFLENGTH*DEFAULT_TREE_MAX_HALFLENGTH),
                initroot_halfleng(DEFAULT_TREE_INIT_ROOT_HALFLENGTH),
                cluster_halfleng(DEFAULT_TREE_CLUSTER_HALFLENGTH),
                cluster_halfleng_sqr(DEFAULT_TREE_CLUSTER_HALFLENGTH*DEFAULT_TREE_CLUSTER_HALFLENGTH){}
    tree_param_(float mi, float ma, float ini, float c):
                min_halfleng(mi),
                min_halfleng_sqr(mi*mi),
                max_halfleng(ma),
                max_halfleng_sqr(ma*ma),
                initroot_halfleng(ini),
                cluster_halfleng(c),
                cluster_halfleng_sqr(c*c){}
} tree_param;

#endif
