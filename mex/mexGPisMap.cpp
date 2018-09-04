#include "mex.h"
#include <vector>
#include <memory>
#include <iostream>
#include <chrono>
#include "GPisMap.h"

static GPisMap* gpm = 0;

static double test_time = 0.0;
static std::chrono::high_resolution_clock::time_point t1;
static std::chrono::high_resolution_clock::time_point t2;

void mexFunction (int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {

    char command[128];
    mxGetString(prhs[0],command,128);
    std::string commstr(command);

    if (commstr.compare("init")==0) {

        if (gpm != 0){
            delete gpm;
            gpm = 0;
        }

        if (nrhs < 2){
            gpm = new GPisMap();
        }
        else{
            size_t numel = mxGetNumberOfElements(prhs[1]);

            if (numel == 8){
                if (sizeof(FLOAT) == sizeof(double)){ // double-precision
                    if (mxDOUBLE_CLASS != mxGetClassID(prhs[1]) ){
                        std::cout << "The input data must be double type. @ train()" << std::endl;
                        return;
                    }
                }
                else if (sizeof(FLOAT) == sizeof(float)) { // single-precision
                    if (mxSINGLE_CLASS != mxGetClassID(prhs[1]) ){
                        std::cout << "The input data must be float (single) type.  @ train()" << std::endl;
                        return;
                    }
                }
                else {
                    std::cout << "The input data must be float (single or double) type. @ train()" << std::endl;
                    return;
                }

                GPisMapParam p;
                // TO_DO: take struct as inpur and set the params;
                //FLOAT *px = (FLOAT *)mxGetData(prhs[1]);
                //GPisMapParam p(px[0], px[1], px[2], (int)px[3],(px+4), px[6], px[7]);
                gpm = new GPisMap(p);
            }
        }
        return;
    }
    else if (commstr.compare("update")==0) {

        if (gpm == 0){
            gpm = new GPisMap();
        }
        mwSize numDim = mxGetNumberOfDimensions(prhs[1]);
        const mwSize *dims = mxGetDimensions(prhs[1]);

        if ( ((dims[0]+1) == gpm->getMapDimension()) && (nrhs == 4))
        {
            if (sizeof(FLOAT) == sizeof(double)){ // double-precision
                if (mxDOUBLE_CLASS != mxGetClassID(prhs[1]) ||
                    mxDOUBLE_CLASS != mxGetClassID(prhs[2]) ||
                    mxDOUBLE_CLASS != mxGetClassID(prhs[3]) ){
                    std::cout << "The input data must be double type. @ train()" << std::endl;
                    return;
                }
            }
            else if (sizeof(FLOAT) == sizeof(float)) { // single-precision
                if (mxSINGLE_CLASS != mxGetClassID(prhs[1]) ||
                    mxSINGLE_CLASS != mxGetClassID(prhs[2]) ||
                    mxSINGLE_CLASS != mxGetClassID(prhs[3]) ){
                    std::cout << "The input data must be float (single) type.  @ train()" << std::endl;
                    return;
                }
            }
            else {
                std::cout << "The input data must be float (single or double) type. @ train()" << std::endl;
                return;
            }

            FLOAT * pa = (FLOAT *)mxGetData(prhs[1]);
            FLOAT * pr = (FLOAT *)mxGetData(prhs[2]);
            FLOAT * ppose = (FLOAT *)mxGetData(prhs[3]);

            std::vector<FLOAT> pose(ppose, ppose + mxGetNumberOfElements(prhs[3]) );

            if (nlhs > 0){
                plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

            }

            int numel =  mxGetNumberOfElements(prhs[2]);
            t1 = std::chrono::high_resolution_clock::now();
            gpm->update(pa, pr, numel, pose);
            t2= std::chrono::high_resolution_clock::now();

            // test time
            if (nlhs > 0){
                std::chrono::duration<double> time_collapsed = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1); // reset
                test_time = time_collapsed.count();
                double *ptime   = (double*)mxGetPr(plhs[0]);
                ptime[0] = test_time;
            }

        }
    }
    else if (commstr.compare("update-mt")==0) {

        if (gpm == 0){
            gpm = new GPisMap();
        }
        mwSize numDim = mxGetNumberOfDimensions(prhs[1]);
        const mwSize *dims = mxGetDimensions(prhs[1]);

        if ( ((dims[0]+1) == gpm->getMapDimension()) && (nrhs == 4))
        {
            if (sizeof(FLOAT) == sizeof(double)){ // double-precision
                if (mxDOUBLE_CLASS != mxGetClassID(prhs[1]) ||
                    mxDOUBLE_CLASS != mxGetClassID(prhs[2]) ||
                    mxDOUBLE_CLASS != mxGetClassID(prhs[3]) ){
                    std::cout << "The input data must be double type. @ train()" << std::endl;
                    return;
                }
            }
            else if (sizeof(FLOAT) == sizeof(float)) { // single-precision
                if (mxSINGLE_CLASS != mxGetClassID(prhs[1]) ||
                    mxSINGLE_CLASS != mxGetClassID(prhs[2]) ||
                    mxSINGLE_CLASS != mxGetClassID(prhs[3]) ){
                    std::cout << "The input data must be float (single) type.  @ train()" << std::endl;
                    return;
                }
            }
            else {
                std::cout << "The input data must be float (single or double) type. @ train()" << std::endl;
                return;
            }

            FLOAT * pa = (FLOAT *)mxGetData(prhs[1]);
            FLOAT * pr = (FLOAT *)mxGetData(prhs[2]);
            FLOAT * ppose = (FLOAT *)mxGetData(prhs[3]);

            std::vector<FLOAT> pose(ppose, ppose + mxGetNumberOfElements(prhs[3]) );

            if (nlhs > 0){
                plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

            }

            int numel =  mxGetNumberOfElements(prhs[2]);
            t1 = std::chrono::high_resolution_clock::now();
            gpm->update_mt(pa, pr, numel, pose);
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
            std::cout << "Error: Check the inpur dimension." << std::endl;
        return;
    }
    else if (commstr.compare("test")==0) {
        const mwSize *dims = mxGetDimensions(prhs[1]);
        int dim = dims[0];
        int N = dims[1];

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

            int N = pos.size()/2; // 2D
            if (N > 0){
                if (sizeof(FLOAT) == sizeof(double)){
                    plhs[0] = mxCreateDoubleMatrix(2,N,mxREAL);
                }
                else if (sizeof(FLOAT) == sizeof(float)){
                    plhs[0] = mxCreateNumericMatrix(2,N,mxSINGLE_CLASS, mxREAL);
                }
                FLOAT *pRes   = (FLOAT*)mxGetPr(plhs[0]);
                memcpy(pRes,pos.data(),sizeof(FLOAT)*2*N);
            }
        }
    }
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
    }
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
