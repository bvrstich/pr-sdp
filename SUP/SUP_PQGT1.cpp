#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "../headers/SUP/SUP_PQGT1.h"
#include "../headers/EIG/EIG_PQGT1.h"

/**
 * Standard constructor, allocates an extra DPM matrix in the pointer SZ_dp on top of calling the constructor of SUP_PQG. 
 * @param M dim of sp space
 * @param N nr of particles
 */
SUP_PQGT1::SUP_PQGT1(int M,int N) : SUP_PQG(M,N){

   this->n_dp = M*(M - 1)*(M - 2)/6;

   dim += n_dp;

   SZ_dp = new DPM(M,N);

}

/**
 * Standard constructor, allocates an extra DPM matrix in the pointer SZ_dp on top of calling the copy constructor of SUP_PQG with input SZ_c.
 * SZ_c.dpm() is copied into SZ_dp
 * @param SZ_c input SUP_PQGT1
 */
SUP_PQGT1::SUP_PQGT1(SUP_PQGT1 &SZ_c) : SUP_PQG(SZ_c){

   this->n_dp = SZ_c.n_dp;

   dim += n_dp;

   SZ_dp = new DPM(M,N);

   (*SZ_dp) = (*SZ_c.SZ_dp);

}

/**
 * destructor, deallocates the SZ_dp memory
 */
SUP_PQGT1::~SUP_PQGT1(){

   delete SZ_dp;

}

/**
 * += operator overloaded
 * @param SZ_pl SUP_PQGT1 to add to this
 */
SUP_PQGT1 &SUP_PQGT1::operator+=(SUP_PQGT1 &SZ_pl){

   this->SUP_PQG::operator+=(SZ_pl);

   (*SZ_dp) += (*SZ_pl.SZ_dp);

   return *this;

}

/**
 * -= operator overloaded
 * @param SZ_pl SUP_PQGT1 to be subtracted from this
 */
SUP_PQGT1 &SUP_PQGT1::operator-=(SUP_PQGT1 &SZ_pl){

   this->SUP_PQG::operator-=(SZ_pl);

   (*SZ_dp) -= (*SZ_pl.SZ_dp);

   return *this;

}

/**
 * Overload equality operator, copy SZ_c into this
 * @param SZ_c SUP_PQGT1 to be copied into this
 */
SUP_PQGT1 &SUP_PQGT1::operator=(SUP_PQGT1 &SZ_c){

   this->SUP_PQG::operator=(SZ_c);

   (*SZ_dp) = (*SZ_c.SZ_dp);

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP_PQGT1 &SUP_PQGT1::operator=(double &a){

   this->SUP_PQG::operator=(a);

   (*SZ_dp) = a;

   return *this;

}

ostream &operator<<(ostream &output,SUP_PQGT1 &SZ_p){

   for(int i = 0;i < 2;++i)
      output << SZ_p.tpm(i) << std::endl;

   output << SZ_p.phm() << std::endl;

   output << SZ_p.dpm() << std::endl;

   return output;

}

/**
 * @return dimension of dp space
 */
int SUP_PQGT1::gn_dp(){

   return n_dp;

}

/**
 * @return pointer to DPM object SZ_dp
 */
DPM &SUP_PQGT1::dpm(){

   return *SZ_dp;

}

/**
 * @param SZ_i input SUP_PQGT1 SZ_i
 * @return inproduct between this and input matrix SZ_i, defined as Tr(this SZ_i)
 */
double SUP_PQGT1::ddot(SUP_PQGT1 &SZ_i){

   double ward = this->SUP_PQG::ddot(SZ_i);

   ward += SZ_dp->ddot(SZ_i.dpm());

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 */
void SUP_PQGT1::invert(){

   this->SUP_PQG::invert();

   SZ_dp->invert();

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP_PQGT1::dscal(double alpha){

   this->SUP_PQG::dscal(alpha);

   SZ_dp->dscal(alpha);

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP_PQGT1::sqrt(int option){

   this->SUP_PQG::sqrt(option);

   SZ_dp->sqrt(option);

}

/**
 * Multiply symmetric SUP_PQGT1 blockmatrix object left en right with symmetric SUP_PQGT1 blockmatrix map to 
 * form another symmetric SUP_PQGT1 matrix and put it in (*this): this = map*object*map
 * @param map SUP_PQGT1 that will be multiplied to the left en to the right of matrix object
 * @param object central SUP_PQGT1
 */
void SUP_PQGT1::L_map(SUP_PQGT1 &map,SUP_PQGT1 &object){

   this->SUP_PQG::L_map(map,object);

   SZ_dp->L_map(*map.SZ_dp,*object.SZ_dp);

}

/**
 * add the SUP_PQGT1 SZ_p times the constant alpha to this
 * @param alpha the constant to multiply the SZ_p with
 * @param SZ_p the SUP_PQGT1 to be multiplied by alpha and added to (*this)
 */
void SUP_PQGT1::daxpy(double alpha,SUP_PQGT1 &SZ_p){

   this->SUP_PQG::daxpy(alpha,SZ_p);

   SZ_dp->daxpy(alpha,*SZ_p.SZ_dp);

}

/**
 * @return trace of SUP_PQGT1, which is defined as the trace of SUP_PQG plus the trace of the DPM SZ_dp
 */
double SUP_PQGT1::trace(){

   double ward = this->SUP_PQG::trace();

   ward += SZ_dp->trace();

   return ward;

}

/**
 * Fill the SUP_PQGT1 with TPM tpm, is defined as the SUP_PQG::fill function and adds:\n\n
 * SZ_dp = DPM::T (tpm)
 * @param tpm input TPM
 */
void SUP_PQGT1::fill(TPM &tpm){

   this->SUP_PQG::fill(tpm);

   SZ_dp->T(tpm);

}
