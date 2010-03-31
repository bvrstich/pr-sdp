#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "../headers/SUP/SUP_PQGT2.h"
#include "../headers/EIG/EIG_PQGT2.h"

/**
 * Standard constructor, allocates an extra PPHM matrix in the pointer SZ_pph on top of calling the constructor of SUP_PQG.
 * @param M dim of sp space
 * @param N nr of particles
 */
SUP_PQGT2::SUP_PQGT2(int M,int N) : SUP_PQG(M,N)
{
   this->n_pph = M*M*(M - 1)/2;

   dim += n_pph;

   SZ_pph = new PPHM(M,N);
}

/**
 * Standard constructor, allocates an extra PPHM matrix in the pointer SZ_pph on top of calling the copy constructor of SUP_PQG with input SZ_c.
 * SZ_c.PPHM() is copied into SZ_pph
 * @param SZ_c input SUP_PQGT2
 */
SUP_PQGT2::SUP_PQGT2(SUP_PQGT2 &SZ_c) : SUP_PQG(SZ_c)
{
   this->n_pph = SZ_c.n_pph;

   dim += n_pph;

   SZ_pph = new PPHM(M,N);

   (*SZ_pph) = (*SZ_c.SZ_pph);
}

/**
 * destructor, deallocates the SZ_pph memory
 */
SUP_PQGT2::~SUP_PQGT2(){

   delete SZ_pph;

}

/**
 * += operator overloaded
 * @param SZ_pl SUP_PQGT2 to add to this
 */
SUP_PQGT2 &SUP_PQGT2::operator+=(SUP_PQGT2 &SZ_pl){

   this->SUP_PQG::operator+=(SZ_pl);

   (*SZ_pph) += (*SZ_pl.SZ_pph);

   return *this;

}

/**
 * -= operator overloaded
 * @param SZ_pl SUP_PQGT2 to be subtracted from this
 */
SUP_PQGT2 &SUP_PQGT2::operator-=(SUP_PQGT2 &SZ_pl){

   this->SUP_PQG::operator-=(SZ_pl);

   (*SZ_pph) -= (*SZ_pl.SZ_pph);

   return *this;

}

/**
 * Overload equality operator, copy SZ_c into this
 * @param SZ_c SUP_PQGT2 to be copied into this
 */
SUP_PQGT2 &SUP_PQGT2::operator=(SUP_PQGT2 &SZ_c){

   this->SUP_PQG::operator=(SZ_c);

   (*SZ_pph) = (*SZ_c.SZ_pph);

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP_PQGT2 &SUP_PQGT2::operator=(double &a){

   this->SUP_PQG::operator=(a);

   (*SZ_pph) = a;

   return *this;

}

ostream &operator<<(ostream &output,SUP_PQGT2 &SZ_p){

   for(int i = 0;i < 2;++i)
      output << SZ_p.tpm(i) << std::endl;

   output << SZ_p.phm() << std::endl;

   output << SZ_p.pphm() << std::endl;

   return output;

}

/**
 * @return dimension of dp space
 */
int SUP_PQGT2::gn_pph(){

   return n_pph;

}

/**
 * @return pointer to PPHM object SZ_pph
 */
PPHM &SUP_PQGT2::pphm(){

   return *SZ_pph;

}

/**
 * @param SZ_i input SUP_PQGT2 SZ_i
 * @return inproduct between this and input matrix SZ_i, defined as Tr(this SZ_i)
 */
double SUP_PQGT2::ddot(SUP_PQGT2 &SZ_i){

   double ward = this->SUP_PQG::ddot(SZ_i);

   ward += SZ_pph->ddot(SZ_i.pphm());

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 */
void SUP_PQGT2::invert(){

   this->SUP_PQG::invert();

   SZ_pph->invert();

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP_PQGT2::dscal(double alpha){

   this->SUP_PQG::dscal(alpha);

   SZ_pph->dscal(alpha);

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP_PQGT2::sqrt(int option){

   this->SUP_PQG::sqrt(option);

   SZ_pph->sqrt(option);

}

/**
 * Multiply symmetric SUP_PQGT2 blockmatrix object left en right with symmetric SUP_PQGT2 blockmatrix map to 
 * form another symmetric SUP_PQGT2 matrix and put it in (*this): this = map*object*map
 * @param map SUP_PQGT2 that will be multiplied to the left en to the right of matrix object
 * @param object central SUP_PQGT2
 */
void SUP_PQGT2::L_map(SUP_PQGT2 &map,SUP_PQGT2 &object){

   this->SUP_PQG::L_map(map,object);

   SZ_pph->L_map(*map.SZ_pph,*object.SZ_pph);

}

/**
 * add the SUP_PQGT2 SZ_p times the constant alpha to this
 * @param alpha the constant to multiply the SZ_p with
 * @param SZ_p the SUP_PQGT2 to be multiplied by alpha and added to (*this)
 */
void SUP_PQGT2::daxpy(double alpha,SUP_PQGT2 &SZ_p){

   this->SUP_PQG::daxpy(alpha,SZ_p);

   SZ_pph->daxpy(alpha,*SZ_p.SZ_pph);

}

/**
 * @return trace of SUP_PQGT2, which is defined as the trace of SUP_PQG plus the trace of the PPHM SZ_pph
 */
double SUP_PQGT2::trace(){

   double ward = this->SUP_PQG::trace();

   ward += SZ_pph->trace();

   return ward;

}

/**
 * Fill the SUP_PQGT2 with TPM tpm, is defined as the SUP_PQG::fill function and adds:\n\n
 * SZ_pph = PPHM::T (tpm)
 * @param tpm input TPM
 */
void SUP_PQGT2::fill(TPM &tpm){

   this->SUP_PQG::fill(tpm);

   SZ_pph->T(tpm);

}
