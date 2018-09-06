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

  if (commstr.compare("update")==0) {

        if (gpm == 0){
            gpm = new GPisMap();
        }
        mwSize numDim = mxGetNumberOfDimensions(prhs[1]);
        const mwSize *dims = mxGetDimensions(prhs[1]);

        if ( ((dims[0]+1) == gpm->getMapDimension()) && (nrhs == 4))
        {
            if (mxSINGLE_CLASS != mxGetClassID(prhs[1]) ||
                mxSINGLE_CLASS != mxGetClassID(prhs[2]) ||
                mxSINGLE_CLASS != mxGetClassID(prhs[3]) ){
                std::cout << "The input data must be float (single) type.  @ train()" << std::endl;
                return;
            }

            float * pa = (float *)mxGetData(prhs[1]);
            float * pr = (float *)mxGetData(prhs[2]);
            float * ppose = (float *)mxGetData(prhs[3]);

            std::vector<float> pose(ppose, ppose + mxGetNumberOfElements(prhs[3]) );

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
        else
            std::cout << "Error: Check the inpur dimension." << std::endl;
        return;
    }
    else if (commstr.compare("test")==0) {
        const mwSize *dims = mxGetDimensions(prhs[1]);
        int dim = dims[0];
        int N = dims[1];

        if (mxSINGLE_CLASS != mxGetClassID(prhs[1])){
            std::cout << "The input data must be float (single) type. @ test()" << std::endl;
            return;
        }

        if ((gpm != 0) && (N > 0)){
            float* px = (float *)mxGetData(prhs[1]);

            plhs[0] = mxCreateNumericMatrix(2*(1+dim),N,mxSINGLE_CLASS, mxREAL);

            if (nlhs == 2){
                plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
                // test time
                t1 = std::chrono::high_resolution_clock::now();
            }

            float *pRes   = (float*)mxGetPr(plhs[0]);
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
