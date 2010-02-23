#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "../headers/EIG/EIG_PQ.h"
#include "../headers/SUP/SUP_PQ.h"
#include "../headers/lapack.h"

/**
 * standard constructor, allocates the memory for the eigenvalues of an SUP_PQ object and makes the pointers point
 * to the right place in the array.
 * @param M dimension of sp space
 * @param N nr of particles
 */
EIG_PQ::EIG_PQ(int M,int N){
   
   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

   this->dim = 2*n_tp;

   eig = new double * [2];

   eig[0] = new double [dim];
   eig[1] = eig[0] + n_tp;

}

/**
 * copy constructor, allocates the memory for the eigenvalues of an SUP_PQ object and makes the pointers point
 * to the right place in the array. Copies the content of eig_c into this.
 * @param eig_c EIG_PQ object that will be copied into this
 */
EIG_PQ::EIG_PQ(EIG_PQ &eig_c){

   this->N = eig_c.N;
   this->M = eig_c.M;
   this->n_tp = eig_c.n_tp;

   this->dim = 2*n_tp;

   eig = new double * [2];

   eig[0] = new double [dim];
   eig[1] = eig[0] + n_tp;

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

}

/**
 * standard constructor with initialization on the eigenvalues of a SUP_PQ object.
 * @param SZ input SUP_PQ object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP_PQ matrix.
 */
EIG_PQ::EIG_PQ(SUP_PQ &SZ){

   this->N = SZ.gN();
   this->M = SZ.gM();
   this->n_tp = SZ.gn_tp();

   this->dim = 2*n_tp;

   eig = new double * [2];

   eig[0] = new double [dim];
   eig[1] = eig[0] + n_tp;

   //dit vernietigd de originele matrix!
   for(int i = 0;i < 2;++i)
      SZ.tpm(i).diagonalize(eig[i]);

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG_PQ &EIG_PQ::operator=(EIG_PQ &eig_c){

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

   return *this;

}

/** 
 * @param i == 0, the eigenvalues of the P block will be returned, i == 1, the eigenvalues of the Q block will be returned.
 * @return array of eigenvalues
 */
double *EIG_PQ::operator[](int i){

   return eig[i];

}

/** 
 * @return nr of particles
 */
int EIG_PQ::gN(){

   return N;

}

/** 
 * @return dimension of sp space
 */
int EIG_PQ::gM(){

   return M;

}

/** 
 * @return dimension of tp space
 */
int EIG_PQ::gn_tp(){

   return n_tp;

}

/**
 * Destructor, deallocates the memory.
 */
EIG_PQ::~EIG_PQ(){

   delete [] eig[0];
   delete [] eig;

}

ostream &operator<<(ostream &output,EIG_PQ &eig_p){

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[0][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.n_tp;++i)
      output << i << "\t" << eig_p.eig[1][i] << std::endl;

   return output;

}

/**
 * acces to the numbers
 * @param block == 0, get element "index" from P block, == 1 get element index from Q block 
 * @param index which element in block you want
 * @return eig[block][index]
 */
double EIG_PQ::operator()(int block,int index){

   return eig[block][index];

}

/**
 * @param alpha step length along the Newton direction
 * @return The line search function, gradient of the potential in the Newton direction as a function of the step length alpha
 */
double EIG_PQ::lsfunc(double alpha){

   double ward = 0.0;

   for(int i = 0;i < n_tp;++i)
      ward += eig[0][i]/(1.0 + alpha*eig[0][i]);

   for(int i = 0;i < n_tp;++i)
      ward += eig[1][i]/(1.0 + alpha*eig[1][i]);

   return ward;

}

/**
 * @return the minimal element present in this EIG_PQ object.
 */
double EIG_PQ::min(){

   double ward;

   if(eig[0][0] < eig[1][0])
      ward = eig[0][0];
   else
      ward = eig[1][0];

   return ward;

}
