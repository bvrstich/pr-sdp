#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "../headers/SUP/SUP_PQ.h"
#include "../headers/EIG/EIG_PQ.h"

/**
 * Standard constructor, allocates two TPM matrices in the pointer SZ_tp.
 * @param M dim of sp space
 * @param N nr of particles
 */
SUP_PQ::SUP_PQ(int M,int N){

   this->M = M;
   this->N = N;
   this->n_tp = M*(M - 1)/2;

   dim = 2*n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

}

/**
 * Copy constructor, allocates two TPM matrices in the pointer SZ_tp, and copies the content of SZ_c into it.
 * @param SZ_c input SUP_PQ
 */
SUP_PQ::SUP_PQ(SUP_PQ &SZ_c){

   this->M = SZ_c.M;
   this->N = SZ_c.N;
   this->n_tp = SZ_c.n_tp;

   dim = 2*n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) = (*SZ_c.SZ_tp[i]);

}

/**
 * Destructor, deallocates memory in SZ_tp
 */
SUP_PQ::~SUP_PQ(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

}

/**
 * += operator overloaded
 * @param SZ_pl SUP_PQ to add to this
 */
SUP_PQ &SUP_PQ::operator+=(SUP_PQ &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

   return *this;

}

/**
 * -= operator overloaded
 * @param SZ_pl SUP_PQ to be subtracted from this
 */
SUP_PQ &SUP_PQ::operator-=(SUP_PQ &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

   return *this;

}

/**
 * Overload equality operator, copy SZ_c into this
 * @param SZ_c SUP_PQ to be copied into this
 */
SUP_PQ &SUP_PQ::operator=(SUP_PQ &SZ_c){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) = (*SZ_c.SZ_tp[i]);

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP_PQ &SUP_PQ::operator=(double &a){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) = a;

   return *this;

}

ostream &operator<<(ostream &output,SUP_PQ &SZ_p){

   for(int i = 0;i < 2;++i)
      output << (*SZ_p.SZ_tp[i]) << std::endl;

   return output;

}

/**
 * @return particle number
 */
int SUP_PQ::gN(){

   return N;

}

/**
 * @return dim of sp space
 */
int SUP_PQ::gM(){

   return M;

}

/**
 * @return dim of tp space
 */
int SUP_PQ::gn_tp(){

   return n_tp;

}

/**
 * @param i which block you want to have the pointer to.
 * @return pointer to the individual TPM blocks: SZ_tp[i]
 */
TPM &SUP_PQ::tpm(int i){

   return *SZ_tp[i];

}

/**
 * @param SZ_i input SUP_PQ SZ_i
 * @return inproduct between this and input matrix SZ_i, defined as Tr(this SZ_i)
 */
double SUP_PQ::ddot(SUP_PQ &SZ_i){

   double ward = 0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 */
void SUP_PQ::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP_PQ::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP_PQ::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

}

/**
 * Multiply symmetric SUP_PQ blockmatrix object left en right with symmetric SUP_PQ blockmatrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map SUP_PQ that will be multiplied to the left en to the right of matrix object
 * @param object central SUP_PQ
 */
void SUP_PQ::L_map(SUP_PQ &map,SUP_PQ &object){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->L_map(*map.SZ_tp[i],*object.SZ_tp[i]);

}

/**
 * add the SUP_PQ SZ_p times the constant alpha to this
 * @param alpha the constant to multiply the SZ_p with
 * @param SZ_p the SUP_PQ to be multiplied by alpha and added to (*this)
 */
void SUP_PQ::daxpy(double alpha,SUP_PQ &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,*SZ_p.SZ_tp[i]);

}

/**
 * @return trace of SUP_PQ, which is defined as the sum of the traces of the individual TPM's
 */
double SUP_PQ::trace(){

   double ward = 0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->trace();

   return ward;

}

/**
 * Fill the SUP_PQ with TPM tpm, which means put:\n\n
 * SZ_tp[0] = tpm\n
 * SZ_tp[1] = TPM::Q (tpm)\n\n
 * @param tpm input TPM
 */
void SUP_PQ::fill(TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(tpm);

}
