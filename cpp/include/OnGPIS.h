#ifndef __ONGPIS_H__
#define __ONGPIS_H__

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

typedef std::vector<std::shared_ptr<Node> > vecNode;
typedef std::vector<std::shared_ptr<Node3> > vecNode3;

class OnGPIS{
      EMatrixX x;
      EMatrixX L;
      EVectorX alpha;
      std::vector<float> gradflag;

      onGPISparam param;     // defined in strct.h
                             // currently noise param is not effective
      float three_over_scale;
      bool trained;
      int nSamples;

public:
    OnGPIS():param(DEFAULT_MAP_SCALE_PARAM,DEFAULT_MAP_NOISE_PARAM,DEFAULT_MAP_NOISE_PARAM),
             three_over_scale(3.0/(DEFAULT_MAP_SCALE_PARAM*DEFAULT_MAP_SCALE_PARAM)),
             trained(false),
             nSamples(0){ }

    OnGPIS(float s,float n):param(s,n,n),
                            three_over_scale(3.0/(s*s)),
                            trained(false),
                            nSamples(0){}

    void reset();
    bool isTrained(){return trained;}
    void setGPScaleParam(float l){param.scale = l;}

    void train(const vecNode& samples);
    void train(const vecNode3& samples);

    void test(const EMatrixX& xt,EVectorX& val, EMatrixX& gradval, EVectorX& var);
    void testSinglePoint(const EVectorX& xt, float& val, float grad[],float var[]);
    void test2Dpoint(const EVectorX& xt, float& val, float& gradx, float& grady, float& varval ,float& vargradx, float& vargrady);
    void test2Dpoint(float x, float y, float& val, float& gradx, float& grady, float& varval ,float& vargradx, float& vargrady);
};

#endif
