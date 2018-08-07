#ifndef __STRCT_H_
#define __STRCT_H_

// strct.h
#include <vector>
#include <memory>

#define DEFAULT_TREE_INIT_ROOT_HALFLENGTH       12.8    // 0.05*2^7
#define DEFAULT_TREE_MIN_HALFLENGTH             0.05     // 0.05*2^0
#define DEFAULT_TREE_MAX_HALFLENGTH             102.4   // 0.05*2^10
#define DEFAULT_TREE_CLUSTER_HALFLENGTH         0.1     // 0.05*2^2

////////////////////////////////////////////////
// #define __USE_DOUBLE_PRECISION__  // Comment this out if using single precision

#ifdef __USE_DOUBLE_PRECISION__
typedef double FLOAT;
#else
typedef float FLOAT;
#endif
////////////////////////////////////////////////

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
    Point<FLOAT> pos;
    Point<FLOAT> grad;
    FLOAT val;
    FLOAT pose_sig;
    FLOAT grad_sig;
    //int level;
    NODE_TYPE nt;

public:

    Node(Point<FLOAT> _pos, FLOAT _val, FLOAT _pose_sig, Point<FLOAT> _grad, FLOAT _grad_sig, NODE_TYPE n);
    Node(Point<FLOAT> _pos, NODE_TYPE _nt = NODE_TYPE::NONE);
    Node();


    void updateData(FLOAT _val, FLOAT _pose_sig, Point<FLOAT> _grad, FLOAT _grad_sig, NODE_TYPE n);
    void updateNoise( FLOAT _pose_sig, FLOAT _grad_sig);

    const Point<FLOAT>& getPos(){return pos;}
    const Point<FLOAT>& getGrad(){return grad;}
    FLOAT getPosX(){return pos.x;}
    FLOAT getPosY(){return pos.y;}
    FLOAT getGradX(){return grad.x;}
    FLOAT getGradY(){return grad.y;}
    FLOAT getVal(){return val;}
    FLOAT getPosNoise(){return pose_sig;};
    FLOAT getGradNoise(){return grad_sig;};
    NODE_TYPE getType(){return nt;}
};

class Node3
{
    Point3<FLOAT> pos;
    Point3<FLOAT> grad;
    FLOAT val;
    FLOAT pose_sig;
    FLOAT grad_sig;
    //int level;
    NODE_TYPE nt;

public:

    Node3(Point3<FLOAT> _pos, FLOAT _val, FLOAT _pose_sig, Point3<FLOAT> _grad, FLOAT _grad_sig, NODE_TYPE n = NODE_TYPE::NONE);
    Node3(Point3<FLOAT> _pos, NODE_TYPE _nt = NODE_TYPE::NONE);
    Node3();

    void updateData(FLOAT _val, FLOAT _pose_sig, Point3<FLOAT> _grad, FLOAT _grad_sig, NODE_TYPE n = NODE_TYPE::NONE);
    void updateNoise( FLOAT _pose_sig, FLOAT _grad_sig);

    const Point3<FLOAT>& getPos(){return pos;}
    const Point3<FLOAT>& getGrad(){return grad;}
    FLOAT getPosX(){return pos.x;}
    FLOAT getPosY(){return pos.y;}
    FLOAT getPosZ(){return pos.z;}
    FLOAT getGradX(){return grad.x;}
    FLOAT getGradY(){return grad.y;}
    FLOAT getGradZ(){return grad.z;}
    FLOAT getVal(){return val;}
    FLOAT getPosNoise(){return pose_sig;};
    FLOAT getGradNoise(){return grad_sig;};
    NODE_TYPE getType(){return nt;}
};

//////////////////////////////////////////////////////////////////////////
// Parameters

// Observation GP
typedef struct obsGPparam_{
    // Npte:
    // ObsGP is implemented to use the Ornstein-Uhlenbeck covariance function,
    // which has a form of k(r)=exp(-r/l) (See covFnc.h)
    FLOAT scale;            // the scale parameter l
    FLOAT noise;            // the noise parameter of the measurement
                            // currently use a constant value
                            // could be potentially modified to have heteroscedastic noise
    // Note:
    // ObsGP is implemented to have overlapping partitioned GPs.
    FLOAT margin;           // used to decide if valid range
                            // (don't use if too close to boundary
                            //  because the derivates are hard to sample)
    int   overlap;          // the overlapping parameters: number of samples to overlap
    int   group_size;       // the number of samples to group together
                            // (the actual group size will be (group_size+overlap)
    obsGPparam_(){}
    obsGPparam_(FLOAT s, FLOAT n, FLOAT m, int ov, int gsz):
                scale(s),
                noise(n),
                margin(m),
                overlap(ov),
                group_size(gsz){ }
} obsGPparam;

// GPIS (SDF)
typedef struct onGPISparam_{
    // Npte:
    // OnlineGPIS is implemented to use the Matern class covariance function with (nu=2/3),
    // which has a form of k(r)=(1+sqrt(3)*r/l)exp(-sqrt(3)*r/l) (See covFnc.h)
    FLOAT scale;            // the scale parameter l
    FLOAT noise;            // the default noise parameter of the measurement
                            // currently use heteroscedastic noise acoording to a noise model
    // TO-DO: add actually parameters (weight, threshold, etc. if any) of the noise model
    FLOAT noise_deriv;      // the default noise parameter of the derivative measurement
                            // currently use a noise model by numerical computation.
    // TO-DO: add actually parameters (weight, threshold, etc. if any) of the noise model
    onGPISparam_(){}
    onGPISparam_(FLOAT s, FLOAT n, FLOAT nd): scale(s), noise(n), noise_deriv(nd){}
} onGPISparam;


// QuadTree (2D) and OcTree (3D)
typedef struct tree_param_{
    FLOAT initroot_halfleng;
    FLOAT min_halfleng;         // minimum (leaf) resolution of tree
    FLOAT min_halfleng_sqr;
    FLOAT max_halfleng;         // maximum (root) resolution of tree
    FLOAT max_halfleng_sqr;
    FLOAT cluster_halfleng;    // the resolution of GP clusters
    FLOAT cluster_halfleng_sqr;
public:
    tree_param_():min_halfleng(DEFAULT_TREE_MIN_HALFLENGTH),
                min_halfleng_sqr(DEFAULT_TREE_MIN_HALFLENGTH*DEFAULT_TREE_MIN_HALFLENGTH),
                max_halfleng(DEFAULT_TREE_MAX_HALFLENGTH),
                max_halfleng_sqr(DEFAULT_TREE_MAX_HALFLENGTH*DEFAULT_TREE_MAX_HALFLENGTH),
                initroot_halfleng(DEFAULT_TREE_INIT_ROOT_HALFLENGTH),
                cluster_halfleng(DEFAULT_TREE_CLUSTER_HALFLENGTH),
                cluster_halfleng_sqr(DEFAULT_TREE_CLUSTER_HALFLENGTH*DEFAULT_TREE_CLUSTER_HALFLENGTH){}
    tree_param_(FLOAT mi, FLOAT ma, FLOAT ini, FLOAT c):
                min_halfleng(mi),
                min_halfleng_sqr(mi*mi),
                max_halfleng(ma),
                max_halfleng_sqr(ma*ma),
                initroot_halfleng(ini),
                cluster_halfleng(c),
                cluster_halfleng_sqr(c*c){}
} tree_param;

#endif
