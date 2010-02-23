#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "../headers/SUP/SUP_PQG.h"
#include "../headers/EIG.h"

/**
 * Standard constructor, allocates an extra PHM matrix in the pointer SZ_ph on top of constructor of SUP_PQ. 
 * @param M dim of sp space
 * @param N nr of particles
 */
SUP_PQG::SUP_PQG(int M,int N) : SUP_PQ(M,N){

   this->n_ph = M*M;

   dim += n_ph;

   SZ_ph = new PHM(M,N);

}

/**
 * Copy constructor, allocates an extra PHM matrix in the pointer SZ_ph on top of the copy constructor of SUP_PQ with input SZ_c.
 * SZ_c.phm() is copied into SZ_ph
 * @param SZ_c input SUP_PQG
 */
SUP_PQG::SUP_PQG(SUP_PQG &SZ_c) : SUP_PQ(SZ_c){

   this->n_ph = SZ_c.n_ph;

   dim += n_ph;

   SZ_ph = new PHM(M,N);

   (*SZ_ph) = (*SZ_c.SZ_ph);

}

/**
 * destructor, deallocates the SZ_ph memory
 */
SUP_PQG::~SUP_PQG(){

   delete SZ_ph;

}

/**
 * += operator overloaded
 * @param SZ_pl SUP_PQG to add to this
 */
SUP_PQG &SUP_PQG::operator+=(SUP_PQG &SZ_pl){

   this->SUP_PQ::operator+=(SZ_pl);

   (*SZ_ph) += (*SZ_pl.SZ_ph);

   return *this;

}

/**
 * -= operator overloaded
 * @param SZ_pl SUP_PQG to be subtracted from this
 */
SUP_PQG &SUP_PQG::operator-=(SUP_PQG &SZ_pl){

   this->SUP_PQ::operator-=(SZ_pl);

   (*SZ_ph) -= (*SZ_pl.SZ_ph);

   return *this;

}

/**
 * Overload equality operator, copy SZ_c into this
 * @param SZ_c SUP_PQG to be copied into this
 */
SUP_PQG &SUP_PQG::operator=(SUP_PQG &SZ_c){

   this->SUP_PQ::operator=(SZ_c);

   (*SZ_ph) = (*SZ_c.SZ_ph);

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP_PQG &SUP_PQG::operator=(double &a){

   this->SUP_PQ::operator=(a);

   (*SZ_ph) = a;

   return *this;

}

ostream &operator<<(ostream &output,SUP_PQG &SZ_p){

   for(int i = 0;i < 2;++i)
      output << SZ_p.tpm(i) << std::endl;

   output << (*SZ_p.SZ_ph) << std::endl;

   return output;

}

/**
 * @return dimension of particle hole space
 */
int SUP_PQG::gn_ph(){

   return n_ph;

}

/**
 * @return pointer to PHM object SZ_ph
 */
PHM &SUP_PQG::phm(){

   return *SZ_ph;

}

/**
 * @param SZ_i input SUP_PQG SZ_i
 * @return inproduct between this and input matrix SZ_i, defined as Tr(this SZ_i)
 */
double SUP_PQG::ddot(SUP_PQG &SZ_i){

   double ward = this->SUP_PQ::ddot(SZ_i);

   ward += SZ_ph->ddot(*SZ_i.SZ_ph);

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 */
void SUP_PQG::invert(){

   this->SUP_PQ::invert();

   SZ_ph->invert();

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP_PQG::dscal(double alpha){

   this->SUP_PQ::dscal(alpha);

   SZ_ph->dscal(alpha);

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP_PQG::sqrt(int option){

   this->SUP_PQ::sqrt(option);

   SZ_ph->sqrt(option);

}

/**
 * Multiply symmetric SUP_PQG blockmatrix object left en right with symmetric SUP_PQG blockmatrix map to 
 * form another symmetric SUP_PQG matrix and put it in (*this): this = map*object*map
 * @param map SUP_PQG that will be multiplied to the left en to the right of matrix object
 * @param object central SUP_PQG
 */
void SUP_PQG::L_map(SUP_PQG &map,SUP_PQG &object){

   this->SUP_PQ::L_map(map,object);

   SZ_ph->L_map(*map.SZ_ph,*object.SZ_ph);

}

/**
 * add the SUP_PQG SZ_p times the constant alpha to this
 * @param alpha the constant to multiply the SZ_p with
 * @param SZ_p the SUP_PQG to be multiplied by alpha and added to (*this)
 */
void SUP_PQG::daxpy(double alpha,SUP_PQG &SZ_p){

   this->SUP_PQ::daxpy(alpha,SZ_p);

   SZ_ph->daxpy(alpha,*SZ_p.SZ_ph);

}

/**
 * @return trace of SUP_PQG, which is defined as the trace of SUP_PQ plus the trace of the PHM SZ_ph
 */
double SUP_PQG::trace(){

   double ward = this->SUP_PQ::trace();

   ward += SZ_ph->trace();

   return ward;

}

/**
 * Fill the SUP_PQG with TPM tpm, is defined as the SUP_PQ::fill function and adds:\n\n
 * SZ_ph = PHM::G (tpm)
 * @param tpm input TPM
 */
void SUP_PQG::fill(TPM &tpm){

   this->SUP_PQ::fill(tpm);

   SZ_ph->G(tpm);

}
