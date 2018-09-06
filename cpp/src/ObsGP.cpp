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

#include "ObsGP.h"
#include "covFnc.h"
#include <Eigen/Cholesky>
#include <thread>

using namespace Eigen;

///////////////////////////////////////////////////////////
// GPou
///////////////////////////////////////////////////////////

void GPou::train(const EMatrixX& xt,const EVectorX& f)
{
    int dim = xt.rows();
    int N = xt.cols();

    if (N > 0){
        x = xt;
        EMatrixX K = ornstein_uhlenbeck(xt, scale, noise);

        L = K.llt().matrixL();
        alpha = f;
        L.template triangularView<Lower>().solveInPlace(alpha);
        L.transpose().template triangularView<Upper>().solveInPlace(alpha);

        trained = true;
    }
}

void GPou::test(const EMatrixX& xt,EVectorX& f, EVectorX& var) // Is is different?
{

    EMatrixX K = ornstein_uhlenbeck(x, xt, scale);
    f = K.transpose()*alpha;

    L.template triangularView<Lower>().solveInPlace(K);

    K = K.array().pow(2);
    EVectorX v = K.colwise().sum();

    var = 1 + noise -v.head(xt.cols()).array();
}

///////////////////////////////////////////////////////////
// ObsGP
///////////////////////////////////////////////////////////

void ObsGP::reset(){
    trained = false;
    gps.clear();
    return;
}

///////////////////////////////////////////////////////////
// ObsGP 1D
///////////////////////////////////////////////////////////

void ObsGP1D::reset(){
    ObsGP::reset();
    range.clear();
    nSamples = 0;
    return;
}

void ObsGP1D::train( float xt[],  float f[], int N[])
{
    reset();

    if ((N[0] > 0) && (xt !=0)){
        nSamples = N[0];
        nGroup = nSamples/(param.group_size) + 1;

        range.push_back(xt[0]);
        for (int n=0;n<(nGroup-1);n++){
            // Make sure there are enough overlap

            if (n < nGroup-2)
            {
                int i1 = n*param.group_size;
                int i2 = i1 + param.group_size + param.overlap;

                range.push_back(xt[i2-param.overlap/2]);

                Map<ERowVectorX> x_(xt+i1,param.group_size + param.overlap);
                Map<EVectorX> f_(f+i1,param.group_size + param.overlap);
                // train each gp group
                std::shared_ptr<GPou> g(new GPou());
                g->train(x_,f_);

                gps.push_back(std::move(g));

            }
            else{ // the last two groups split in half
                // the second to last
                int i1 = n*param.group_size;
                int i2 = i1 + (nSamples - i1)/2 + param.overlap;
                range.push_back(xt[i2-param.overlap/2]);

                Map<ERowVectorX> x_(xt+i1,i2-i1+1);
                Map<EVectorX> f_(f+i1,i2-i1+1);
                std::shared_ptr<GPou> g(new GPou());
                g->train(x_,f_);
                gps.push_back(std::move(g));
                n++;

                // the last one
                i1 = i1 + (nSamples - i1)/2;
                i2 = nSamples-1;
                range.push_back(xt[i2]);
                new (&x_) Map<ERowVectorX>(xt+i1,i2-i1+1);
                new (&f_) Map<EVectorX>(f+i1,i2-i1+1);

                std::shared_ptr<GPou> glast(new GPou());
                glast->train(x_,f_);
                gps.push_back(std::move(glast));
            }
        }

        trained = true;
    }

    return;
}

void ObsGP1D::test(const EMatrixX& xt,EVectorX& val, EVectorX& var){

    if (!isTrained()){
        return;
    }

    int dim = xt.rows();
    int N = xt.cols();

    if (dim ==1){
        float liml = (*(range.begin())+param.margin);
        float limr = (*(range.end()-1)-param.margin);
        for (int k=0;k<N;k++){

            EVectorX f = val.segment(k,1);
            EVectorX v = var.segment(k,1);
            var(k) = 1e6;
            // find the corresponding group
            if (xt(0,k) < liml){ // boundary 1
                ;
            }
            else if (xt(0,k) > limr){ // boundary 2
                ;
            }
            else{ // in-between
                int j = 0;
                for (auto it = (range.begin()+1); it!=range.end() ; it++,j++){
                    if (xt(0,k) > *(it-1) && xt(0,k) < *it){
                        // and test
                        if (gps[j]->isTrained()){
                            gps[j]->test(xt.block(0,k,1,1), f,v);
                            val(k) = f(0);
                            var(k) = v(0);
                        }
                        break;
                    }
                }
            }
        }
    }

    return;
}

///////////////////////////////////////////////////////////
// ObsGP 2D
///////////////////////////////////////////////////////////

void ObsGP2D::clearGPs(){
    ObsGP::reset();
}

void ObsGP2D::reset(){
    clearGPs();

    repartition = true;
    return;
}

void ObsGP2D::computePartition(float val[], int ni, int nj)
{
    // number of data grid
    szSamples[0] = ni;
    szSamples[1] = nj;

    // number group gp grid
    nGroup[0] = (szSamples[0]-param.overlap)/param.group_size +1;
    nGroup[1] = (szSamples[1]-param.overlap)/param.group_size + 1;

    Ind_i0.clear();
    Ind_i1.clear();
    Ind_j0.clear();
    Ind_j1.clear();
    Val_i.clear();
    Val_j.clear();

    // [0]-Range for each group
    Val_i.push_back(val[0]);
    for (int n=0;n<nGroup[0];n++){

        int i0 = n*param.group_size;
        int i1 = i0 + param.group_size + param.overlap -1;
        if (n < nGroup[0]-1)
        {
            Val_i.push_back(val[2*(i1-param.overlap/2)]);
        }
        else if (n == nGroup[0]-1)
        {
            i1 = szSamples[0]-1;
            Val_i.push_back(val[2*i1]);
        }
        Ind_i0.push_back(i0);
        Ind_i1.push_back(i1);

    }

    // [1]-Range for each group
    Val_j.push_back(val[1]);
    for (int m=0;m<nGroup[1];m++){

        int j0 = m*param.group_size;
        int j1 = j0 + param.group_size + param.overlap-1;
        if (m < nGroup[1]-1)
        {
            Val_j.push_back(val[2*(j1-param.overlap/2)*szSamples[0]+1]);
        }
        else {
            // the last one
            j1 = szSamples[1]-1;
            Val_j.push_back(val[2*j1*szSamples[0]+1]);
        }

        Ind_j0.push_back(j0);
        Ind_j1.push_back(j1);
    }

    if ( (Ind_i0.size() > 0) && (Ind_i1.size() > 0) && (Ind_j0.size() > 0) && (Ind_j1.size() > 0))
        repartition = false;

    return;
}

void ObsGP2D::getNumValidPoints(std::vector<int> &nPts)
{
    nPts.clear();
    for (auto it = gps.begin(); it!=gps.end(); it++){
        if (*it != nullptr)
            nPts.push_back((*it)->getNumSamples());
        else
            nPts.push_back(0);
    }

    return;
}

void ObsGP2D::trainValidPoints(float xt[], float f[])
{
    if (repartition)
        return;

     clearGPs();
     gps.resize(nGroup[0]*nGroup[1],nullptr);

     auto itj0 = Ind_j0.begin();
     auto itj1 = Ind_j1.begin();
     int m=0;
     for (;(itj0 != Ind_j0.end() && itj1 != Ind_j1.end() && m < nGroup[1]);itj0++,itj1++,m++){

        auto iti0 = Ind_i0.begin();
        auto iti1 = Ind_i1.begin();
        int n=0;
        for (;(iti0 != Ind_i0.end() && iti1 != Ind_i1.end() && n < nGroup[0]);iti0++,iti1++,n++){
            // Dynamic array for valid inputs
            std::vector<float> x_valid; // 2D array
            std::vector<float> f_valid;

            for (int j=*itj0; j<=*itj1; j++){
                for (int i=*iti0; i<=*iti1; i++){
                    int ind = j*szSamples[0]+i;
                    if (f[ind] > 0){
                        x_valid.push_back(xt[ind*2]);
                        x_valid.push_back(xt[ind*2+1]);
                        f_valid.push_back(f[ind]);
                    }
                }
            }

            // If not empty
            if (x_valid.size() > 1){
                // matrix/vector map from vector
                Map<EMatrixX> x_(x_valid.data(),2,f_valid.size());
                Map<EVectorX> f_(f_valid.data(),f_valid.size());

                // train each gp group
                std::shared_ptr<GPou> g(new GPou());
                g->train(x_,f_);

                gps[m*nGroup[0]+n] = std::move(g);
            }
        }
    }

    trained = true;
    return;
}

void ObsGP2D::train( float xt[], float f[], int N[])
{
    if ((N[0] > 0) && (N[1] > 0) && (xt !=0)){

        if ((szSamples[0] != N[0]) || (szSamples[1] != N[1]) || repartition){
            computePartition(xt,N[0],N[1]);
        }
        trainValidPoints(xt, f);
    }

    return;
}

void ObsGP2D::train( float xt[],  float f[], int N[], std::vector<int>& numSamples)
{
    train(xt, f, N);
    getNumValidPoints(numSamples);

    return;
}

void ObsGP2D::test_kernel(int thread_idx,
                          int start_idx,
                          int end_idx,
                          const EMatrixX &xt,
                          EVectorX &val,
                          EVectorX &var){

    for (int k = start_idx; k < end_idx; ++k){

        EVectorX f = val.segment(k,1);
        EVectorX v = var.segment(k,1);
        var(k) = 1e6;
        // find the corresponding group
        if (xt(0,k) < *(Val_i.begin())+param.margin ){ // boundary 1
            ;

        }
        else if (xt(0,k) > *(Val_i.end()-1)-param.margin){ // boundary 2
            ;
        }
        else if (xt(1,k) < *(Val_j.begin())+param.margin ){ // boundary 1
            ;
        }
        else if (xt(1,k) > *(Val_j.end()-1)-param.margin){ // boundary 2
            ;
        }
        else{ // in-between

            int n = 0;
            for (auto it = (Val_i.begin()+1); it!=Val_i.end() ; it++,n++){
                if (xt(0,k) < *it){
                    break;
                }
            }

            int m = 0;
            for (auto it = (Val_j.begin()+1); it!=Val_j.end() ; it++,m++){
                if (xt(1,k) < *it){
                    break;
                }
            }

            if (1){
                int gp_ind = m*nGroup[0]+n;
                if (gp_ind < gps.size() && gps[gp_ind] != nullptr){
                    if (gps[gp_ind]->isTrained()){
                        // and test
                        gps[gp_ind]->test(xt.block(0,k,2,1), f,v);
                        val(k) = f(0);
                        var(k) = v(0);
                    }
                }
            }
        }
    }
    return;
}

void ObsGP2D::test(const EMatrixX& xt,EVectorX& val, EVectorX& var){

    if (!isTrained() || xt.rows() != 2){
        return;
    }

    int N = xt.cols();

    int num_threads = std::thread::hardware_concurrency();
    int num_threads_to_use = num_threads;
    if (N < num_threads){
        num_threads_to_use = N;
    }
    else{
        num_threads_to_use = num_threads;
    }
    std::thread *threads = new std::thread[num_threads_to_use];

    int num_leftovers = N % num_threads_to_use;
    int batch_size = N / num_threads_to_use;
    int element_cursor = 0;

    for(int i = 0; i < num_leftovers; ++i){
        threads[i] = std::thread(&ObsGP2D::test_kernel,
                                 this,
                                 i,
                                 element_cursor,
                                 element_cursor + batch_size + 1,
                                 std::ref(xt),
                                 std::ref(val),
                                 std::ref(var));
        element_cursor += batch_size + 1;

    }
    for (int i = num_leftovers; i < num_threads_to_use; ++i){
        threads[i] = std::thread(&ObsGP2D::test_kernel,
                                 this,
                                 i,
                                 element_cursor,
                                 element_cursor + batch_size,
                                 std::ref(xt),
                                 std::ref(val),
                                 std::ref(var));
        element_cursor += batch_size;
    }

    for (int i = 0; i < num_threads_to_use; ++i){
        threads[i].join();
    }

    delete [] threads;

    return;
}
