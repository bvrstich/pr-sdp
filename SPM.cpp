#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "headers/include.h"

/**
 * constructor, makes matrix of dimension M
 * @param M dimension of single particle space and dimension of the Matrix
 * @param N Nr of particles
 */
SPM::SPM(int M,int N) : Matrix(M) {

   this->M = M;
   this->N = N;

}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(SPM &spm_copy) : Matrix(spm_copy) {

   this->M = spm_copy.gM();
   this->N = spm_copy.gN();

}

/**
 * destructor
 */
SPM::~SPM(){

}

/**
 * @return nr of particles
 */
int SPM::gN(){

   return N;

}

/**
 * @return dimension of sp space and of matrix
 */
int SPM::gM(){

   return M;

}

ostream &operator<<(ostream &output,SPM &spm_p){

   for(int i = 0;i < spm_p.M;++i)
      for(int j = 0;j < spm_p.M;++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}

/**
 * This function is a template specialization of bar() for the PPHM matrix.
 * It's the so called A dubble bar
 * @param MT A PPHM matrix to contract
 */
template<> void SPM::bar(PPHM &MT)
{
    for(int a = 0;a < M;a++)
	for(int b = a;b < M;b++)
	{
	    (*this)(a,b) = 0.0;

	    for(int l=0;l<M;l++)
		for(int k=0;k<M;k++)
		    (*this)(a,b) += MT(l,k,a,l,k,b);

	    (*this)(b,a) = (*this)(a,b);
	}
}

/* vim: set ts=3 sw=3 tw=3 expandtab :*/
