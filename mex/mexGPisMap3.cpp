#include "mex.h"
#include <vector>
#include <memory>
#include <iostream>
#include <chrono>
#include "GPisMap3.h"

static GPisMap3* gpm = 0;

static double test_time = 0.0;
static std::chrono::high_resolution_clock::time_point t1;
static std::chrono::high_resolution_clock::time_point t2;
static const double bigbirdCams_fx[5] = {570.9361, 572.3318, 568.9403 , 567.9881, 572.7638};
static const double bigbirdCams_fy[5] = {570.9376, 572.3316, 568.9419 , 567.9995, 572.7567};
static const double bigbirdCams_cx[5] = {306.8789, 309.9968, 308.4583, 310.5243, 310.4192};
static const double bigbirdCams_cy[5] = {238.8476, 230.6296, 225.8232, 223.9443, 214.8762};
static const int bigbirdCams_width = 640;
static const int bigbirdCams_height = 480;
static const double YCBCams_fx[5] = {570.2590, 571.8461, 568.4464 , 566.9790, 574.0641};
static const double YCBCams_fy[5] = {570.2636, 571.8428, 568.4494 , 566.9812, 574.0598};
static const double YCBCams_cx[5] = {313.7235, 314.9134, 310.3805, 314.3801, 314.6690};
static const double YCBCams_cy[5] = {236.0783, 229.4538, 224.6232, 223.9443, 220.7985};
static const int YCBCams_width = 640;
static const int YCBCams_height = 480;

void mexFunction (int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {

    char command[128];
    mxGetString(prhs[0],command,128);
    std::string commstr(command);

//     if (commstr.compare("init")==0) {
//
//         if (gpm != 0){
//             delete gpm;
//             gpm = 0;
//         }
//
//         if (nrhs < 2){
//             gpm = new GPisMap3();
//         }
//         else{
//             size_t numel = mxGetNumberOfElements(prhs[1]);
//
//             if (numel == 8){
//                 if (sizeof(FLOAT) == sizeof(double)){ // double-precision
//                     if (mxDOUBLE_CLASS != mxGetClassID(prhs[1]) ){
//                         std::cout << "The input data must be double type. @ train()" << std::endl;
//                         return;
//                     }
//                 }
//                 else if (sizeof(FLOAT) == sizeof(float)) { // single-precision
//                     if (mxSINGLE_CLASS != mxGetClassID(prhs[1]) ){
//                         std::cout << "The input data must be float (single) type.  @ train()" << std::endl;
//                         return;
//                     }
//                 }
//                 else {
//                     std::cout << "The input data must be float (single or double) type. @ train()" << std::endl;
//                     return;
//                 }
//
//                 GPisMap3Param p;
//                 camParam c;
//                 // TO_DO: take struct as inpur and set the params;
//                 //FLOAT *px = (FLOAT *)mxGetData(prhs[1]);
//                 //GPisMapParam p(px[0], px[1], px[2], (int)px[3],(px+4), px[6], px[7]);
//                 gpm = new GPisMap3(p,c);
//             }
//         }
//         return;
//     }
//     else
    if (commstr.compare("update")==0) {

        if (gpm == 0){
            gpm = new GPisMap3();
        }
        mwSize numDim = mxGetNumberOfDimensions(prhs[1]);
        const mwSize *dims = mxGetDimensions(prhs[1]);

        if ( numDim == 2 && (nrhs == 3))
        {
            if (sizeof(FLOAT) == sizeof(double)){ // double-precision
                if (mxDOUBLE_CLASS != mxGetClassID(prhs[1]) ||
                    mxDOUBLE_CLASS != mxGetClassID(prhs[2]) ){
                    std::cout << "The input data must be double type. @ train()" << std::endl;
                    return;
                }
            }
            else if (sizeof(FLOAT) == sizeof(float)) { // single-precision
                if (mxSINGLE_CLASS != mxGetClassID(prhs[1]) ||
                    mxSINGLE_CLASS != mxGetClassID(prhs[2]) ){
                    std::cout << "The input data must be float (single) type.  @ train()" << std::endl;
                    return;
                }
            }
            else {
                std::cout << "The input data must be float (single or double) type. @ train()" << std::endl;
                return;
            }

            FLOAT * pz = (FLOAT *)mxGetData(prhs[1]);
            FLOAT * ppose = (FLOAT *)mxGetData(prhs[2]);

            std::vector<FLOAT> pose(ppose, ppose + mxGetNumberOfElements(prhs[2]) );

            if (nlhs > 0){
                plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            }

            int numel =  mxGetNumberOfElements(prhs[1]);
            t1 = std::chrono::high_resolution_clock::now();
            gpm->update(pz, numel, pose);
            t2= std::chrono::high_resolution_clock::now();

            // test time
            if (nlhs > 0){
                std::chrono::duration<double> time_collapsed = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1); // reset
                test_time = time_collapsed.count();
                double *ptime   = (double*)mxGetPr(plhs[0]);
                ptime[0] = test_time;
            }

        }
        else
            std::cout << "Error: Check the input dimension." << std::endl;
        return;
    }
    else if (commstr.compare("update-mt")==0) {

        if (gpm == 0){
            gpm = new GPisMap3();
        }
        mwSize numDim = mxGetNumberOfDimensions(prhs[1]);
        const mwSize *dims = mxGetDimensions(prhs[1]);

        if ( numDim == 2 && (nrhs == 3))
        {
            if (sizeof(FLOAT) == sizeof(double)){ // double-precision
                if (mxDOUBLE_CLASS != mxGetClassID(prhs[1]) ||
                    mxDOUBLE_CLASS != mxGetClassID(prhs[2]) ){
                    std::cout << "The input data must be double type. @ train()" << std::endl;
                    return;
                }
            }
            else if (sizeof(FLOAT) == sizeof(float)) { // single-precision
                if (mxSINGLE_CLASS != mxGetClassID(prhs[1]) ||
                    mxSINGLE_CLASS != mxGetClassID(prhs[2]) ){
                    std::cout << "The input data must be float (single) type.  @ train()" << std::endl;
                    return;
                }
            }
            else {
                std::cout << "The input data must be float (single or double) type. @ train()" << std::endl;
                return;
            }

            FLOAT * pz = (FLOAT *)mxGetData(prhs[1]);
            FLOAT * ppose = (FLOAT *)mxGetData(prhs[2]);

            std::vector<FLOAT> pose(ppose, ppose + mxGetNumberOfElements(prhs[2]) );

            if (nlhs > 0){
                plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            }

            int numel =  mxGetNumberOfElements(prhs[1]);
            t1 = std::chrono::high_resolution_clock::now();
            gpm->update_mt(pz, numel, pose);
            t2= std::chrono::high_resolution_clock::now();

            // test time
            if (nlhs > 0){
                std::chrono::duration<double> time_collapsed = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1); // reset
                test_time = time_collapsed.count();
                double *ptime   = (double*)mxGetPr(plhs[0]);
                ptime[0] = test_time;
            }

        }
        else
            std::cout << "Error: Check the input dimension." << std::endl;
        return;
    }
    else if (commstr.compare("test")==0) {
        const mwSize *dims = mxGetDimensions(prhs[1]);
        int dim = dims[0];
        int N = dims[1];

        if (dim !=3)
            return;

        if (sizeof(FLOAT) == sizeof(double)){ // double-precision
            if (mxDOUBLE_CLASS != mxGetClassID(prhs[1]) ){
                std::cout << "The input data must be double type.  @ test()" << std::endl;
                return;
            }
        }
        else if (sizeof(FLOAT) == sizeof(float))  { // single-precision
            if (mxSINGLE_CLASS != mxGetClassID(prhs[1])){
                std::cout << "The input data must be float (single) type. @ test()" << std::endl;
                return;
            }
        }
        else {
            std::cout << "The input data must be float (single or double) type. @ test()" << std::endl;
            return;
        }

        if ((gpm != 0) && (N > 0)){
            FLOAT* px = (FLOAT *)mxGetData(prhs[1]);

            if (sizeof(FLOAT) == sizeof(double)){
                plhs[0] = mxCreateDoubleMatrix(2*(1+dim),N,mxREAL);   // f-value[0],  grad-value[1,2], variances[3-5]
            }
            else if (sizeof(FLOAT) == sizeof(float)){
                plhs[0] = mxCreateNumericMatrix(2*(1+dim),N,mxSINGLE_CLASS, mxREAL);
            }

            if (nlhs == 2){
                plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
                // test time
                t1 = std::chrono::high_resolution_clock::now();
            }

            //std::cout << "Number of testing point:" <<  N << std::endl;

            FLOAT *pRes   = (FLOAT*)mxGetPr(plhs[0]);
            gpm->test(px, dim, N, pRes);

            // test time
            if (nlhs == 2){
                t2= std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_collapsed = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1); // reset
                test_time = time_collapsed.count();
                double *ptime   = (double*)mxGetPr(plhs[1]);
                ptime[0] = test_time;
            }

        }
        else if (gpm ==0){
            std::cout << "Error: the map is not initialized." << std::endl;
        }

    }
    else if (commstr.compare("test-mt")==0) {
        const mwSize *dims = mxGetDimensions(prhs[1]);
        int dim = dims[0];
        int N = dims[1];

        if (dim !=3)
            return;

        if (sizeof(FLOAT) == sizeof(double)){ // double-precision
            if (mxDOUBLE_CLASS != mxGetClassID(prhs[1]) ){
                std::cout << "The input data must be double type.  @ test()" << std::endl;
                return;
            }
        }
        else if (sizeof(FLOAT) == sizeof(float))  { // single-precision
            if (mxSINGLE_CLASS != mxGetClassID(prhs[1])){
                std::cout << "The input data must be float (single) type. @ test()" << std::endl;
                return;
            }
        }
        else {
            std::cout << "The input data must be float (single or double) type. @ test()" << std::endl;
            return;
        }

        if ((gpm != 0) && (N > 0)){
            FLOAT* px = (FLOAT *)mxGetData(prhs[1]);

            if (sizeof(FLOAT) == sizeof(double)){
                plhs[0] = mxCreateDoubleMatrix(2*(1+dim),N,mxREAL);   // f-value[0],  grad-value[1,2], variances[3-5]
            }
            else if (sizeof(FLOAT) == sizeof(float)){
                plhs[0] = mxCreateNumericMatrix(2*(1+dim),N,mxSINGLE_CLASS, mxREAL);
            }

            if (nlhs == 2){
                plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
                // test time
                t1 = std::chrono::high_resolution_clock::now();
            }

            //std::cout << "Number of testing point:" <<  N << std::endl;

            FLOAT *pRes   = (FLOAT*)mxGetPr(plhs[0]);
            gpm->test_mt(px, dim, N, pRes);

            // test time
            if (nlhs == 2){
                t2= std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_collapsed = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1); // reset
                test_time = time_collapsed.count();
                double *ptime   = (double*)mxGetPr(plhs[1]);
                ptime[0] = test_time;
            }

        }
        else if (gpm ==0){
            std::cout << "Error: the map is not initialized." << std::endl;
        }

    }
    else if (commstr.compare("getAllPoints")==0){
        if (gpm != 0){
            std::vector<FLOAT> pos;
            gpm->getAllPoints(pos);

            int N = pos.size()/3; // 2D
            if (N > 0){
                if (sizeof(FLOAT) == sizeof(double)){
                    plhs[0] = mxCreateDoubleMatrix(3,N,mxREAL);
                }
                else if (sizeof(FLOAT) == sizeof(float)){
                    plhs[0] = mxCreateNumericMatrix(3,N,mxSINGLE_CLASS, mxREAL);
                }
                FLOAT *pRes   = (FLOAT*)mxGetPr(plhs[0]);
                memcpy(pRes,pos.data(),sizeof(FLOAT)*3*N);
            }
        }
    }
    else if (commstr.compare("setCamera")==0) {

        double *camID = (double *)mxGetData(prhs[1]);
        int n = (int)(*camID)-1;

        char cam_name[12];
        mxGetString(prhs[2],cam_name,12);
        std::string camnamestr(cam_name);
        if (gpm != 0){
            if (n>=0 && n<5){
                if (camnamestr.compare("bigbird")==0){
                    camParam c(bigbirdCams_fx[n],bigbirdCams_fy[n],bigbirdCams_cx[n],bigbirdCams_cy[n], bigbirdCams_width, bigbirdCams_height);
                    gpm->resetCam(c);
                }
                else if (camnamestr.compare("ycb")==0){
                    camParam c(YCBCams_fx[n],YCBCams_fy[n],YCBCams_cx[n],YCBCams_cy[n],YCBCams_width,YCBCams_height);
                    gpm->resetCam(c);
                }
            }
        }
        else
        {
            GPisMap3Param p;
            if (camnamestr.compare("bigbird")==0){
                camParam c(bigbirdCams_fx[n],bigbirdCams_fy[n],bigbirdCams_cx[n],bigbirdCams_cy[n], bigbirdCams_width, bigbirdCams_height);
                gpm = new GPisMap3(p,c);
            }
            else if (camnamestr.compare("ycb")==0){
                camParam c(YCBCams_fx[n],YCBCams_fy[n],YCBCams_cx[n],YCBCams_cy[n],YCBCams_width,YCBCams_height);
                gpm = new GPisMap3(p,c);
            }
        }
        return;
    }
    /*
    else if (commstr.compare("getAllPoints2")==0){
        if (gpm != 0){
            std::vector<FLOAT> pos;
            std::vector<FLOAT> var;
            std::vector<FLOAT> grad;
            std::vector<FLOAT> grad_var;
            gpm->getAllPoints(pos,var,grad,grad_var);

            int N = pos.size()/2; // 2D
            if (N > 0){
                if (sizeof(FLOAT) == sizeof(double)){
                    plhs[0] = mxCreateDoubleMatrix(2,N,mxREAL);
                    plhs[1] = mxCreateDoubleMatrix(1,N,mxREAL);
                    plhs[2] = mxCreateDoubleMatrix(2,N,mxREAL);
                    plhs[3] = mxCreateDoubleMatrix(1,N,mxREAL);
                }
                else if (sizeof(FLOAT) == sizeof(float)){
                    plhs[0] = mxCreateNumericMatrix(2,N,mxSINGLE_CLASS, mxREAL);
                    plhs[1] = mxCreateNumericMatrix(1,N,mxSINGLE_CLASS, mxREAL);
                    plhs[2] = mxCreateNumericMatrix(2,N,mxSINGLE_CLASS, mxREAL);
                    plhs[3] = mxCreateNumericMatrix(1,N,mxSINGLE_CLASS, mxREAL);
                }
                FLOAT *ppos = (FLOAT*)mxGetPr(plhs[0]);
                FLOAT *pvar = (FLOAT*)mxGetPr(plhs[1]);
                FLOAT *pgrad = (FLOAT*)mxGetPr(plhs[2]);
                FLOAT *pgradvar = (FLOAT*)mxGetPr(plhs[3]);
                memcpy(ppos,pos.data(),sizeof(FLOAT)*2*N);
                memcpy(pvar,var.data(),sizeof(FLOAT)*N);
                memcpy(pgrad,grad.data(),sizeof(FLOAT)*2*N);
                memcpy(pgradvar,grad_var.data(),sizeof(FLOAT)*N);
            }
        }
    }  */
    else if (commstr.compare("getRuntimes")==0) {
        if (gpm != 0){
            plhs[0] = mxCreateDoubleMatrix(1,4,mxREAL);
            double *ptime = (double*)mxGetPr(plhs[0]);
            ptime[0] = gpm->getRuntime0();
            ptime[1] = gpm->getRuntime1();
            ptime[2] = gpm->getRuntime2();
            ptime[3] = gpm->getRuntime3();
        }
        else if (gpm ==0){
            std::cout << "Error: the map is not initialized." << std::endl;
        }
        return;
    }
    else if (commstr.compare("reset")==0) {
        if (gpm != 0){
            gpm->reset();
            delete gpm;
            gpm = 0;
            std::cout << "Map cleared\n" << std::endl;
        }
        return;
    }

    return;
}
