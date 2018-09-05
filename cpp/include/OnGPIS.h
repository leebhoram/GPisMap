#ifndef __ONGPIS_H__
#define __ONGPIS_H__

#include <vector>
#include <memory>
#include <iostream>
#include <cstdint>
#include <Eigen/Dense>
#include "strct.h"


#ifdef __USE_DOUBLE_PRECISION__     // This is in use
typedef Eigen::MatrixXd EMatrixX;
typedef Eigen::VectorXd EVectorX;
typedef Eigen::RowVectorXd ERowVectorX;
#else                               // For potential transition to float
typedef Eigen::MatrixXf EMatrixX;
typedef Eigen::VectorXf EVectorX;
typedef Eigen::RowVectorXf ERowVectorX;
#endif


#define DEFAULT_MAP_SCALE_PARAM 1 //0.5
#define DEFAULT_MAP_NOISE_PARAM 1e-2

typedef std::vector<std::shared_ptr<Node> > vecNode;
typedef std::vector<std::shared_ptr<Node3> > vecNode3;

class OnGPIS{
      EMatrixX x;
      EMatrixX L;
      EVectorX alpha;
      std::vector<FLOAT> gradflag;

      onGPISparam param;     // defined in strct.h
                             // currently noise param is not effective
      FLOAT three_over_scale;
      bool trained;
      int nSamples;

public:
    OnGPIS():param(DEFAULT_MAP_SCALE_PARAM,DEFAULT_MAP_NOISE_PARAM,DEFAULT_MAP_NOISE_PARAM),
             three_over_scale(3.0/(DEFAULT_MAP_SCALE_PARAM*DEFAULT_MAP_SCALE_PARAM)),
             trained(false),
             nSamples(0){ }

    OnGPIS(FLOAT s,FLOAT n):param(s,n,n),
                            three_over_scale(3.0/(s*s)),
                            trained(false),
                            nSamples(0){}

    void reset();
    bool isTrained(){return trained;}
    void setGPScaleParam(FLOAT l){param.scale = l;}

    void train(const vecNode& samples);
    void train(const vecNode3& samples);

    void test(const EMatrixX& xt,EVectorX& val, EMatrixX& gradval, EVectorX& var);
    void testSinglePoint(const EVectorX& xt, FLOAT& val, FLOAT grad[],FLOAT var[]);
    void test2Dpoint(const EVectorX& xt,FLOAT& val, FLOAT& gradx, FLOAT& grady, FLOAT& varval ,FLOAT& vargradx, FLOAT &vargrady);
    void test2Dpoint(FLOAT x, FLOAT y,FLOAT& val, FLOAT& gradx, FLOAT& grady, FLOAT& varval ,FLOAT& vargradx, FLOAT &vargrady);

   // void trainTemp(const vecNode& samples);
   // void test2DpointTemp(const EVectorX& xt,FLOAT& val, FLOAT& gradx, FLOAT& grady, FLOAT& varval ,FLOAT& vargradx, FLOAT &vargrady);
};
#endif
