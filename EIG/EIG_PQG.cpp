#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "../headers/EIG/EIG_PQG.h"
#include "../headers/SUP/SUP_PQG.h"
#include "../headers/lapack.h"

/**
 * standard constructor, calls the constructor of EIG_PQ and expands it by allocating the memory for the 
 * G-space vector eig_ph.
 * @param M dimension of sp space
 * @param N nr of particles
 */
EIG_PQG::EIG_PQG(int M,int N) : EIG_PQ(M,N){
   
   this->n_ph = M*M;

   dim += n_ph;

   eig_ph = new double [n_ph];

}

/**
 * standard constructor, calls the constructor of EIG_PQ and expands it by allocating the memory for the 
 * G-space vector eig_ph. copies the content of eig_c into this
 * @param eig_c input EIG_PQG
 */
EIG_PQG::EIG_PQG(EIG_PQG &eig_c) : EIG_PQ(eig_c){

   this->n_ph = eig_c.n_ph;

   dim += n_ph;

   eig_ph = new double [n_ph];

   int inc = 1;

   dcopy_(&n_ph,eig_c.eig_ph,&inc,eig_ph,&inc);

}

/**
 * standard constructor, with initialization on the eigenvalues of a SUP_PQG object.
 * @param SZ input SUP_PQG object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP_PQG matrix.
 */
EIG_PQG::EIG_PQG(SUP_PQG &SZ) : EIG_PQ(SZ){

   this->n_ph = SZ.gn_ph();
   dim += n_ph;

   eig_ph = new double [n_ph];

   SZ.phm().diagonalize(eig_ph);

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG_PQG &EIG_PQG::operator=(EIG_PQG &eig_c){

   this->EIG_PQ::operator=(eig_c);

   int inc = 1;

   dcopy_(&n_ph,eig_c.eig_ph,&inc,eig_ph,&inc);

   return *this;

}

/** 
 * @param i < 2 calls EIG_PQ::operator[], == 2 returns eig_ph
 * @return array of eigenvalues
 */
double *EIG_PQG::operator[](int i){

   if(i < 2)
      return this->EIG_PQ::operator[](i);
   else
      return eig_ph;

}

/**
 * @return dimension of particle hole space
 */
int EIG_PQG::gn_ph(){

   return n_ph;

}
 
/**
 * destructor: deallocates the memory of eig_ph
 */
EIG_PQG::~EIG_PQG(){

   delete [] eig_ph;

}

ostream &operator<<(ostream &output,EIG_PQG &eig_p){

   for(int i = 0;i < eig_p.gn_tp();++i)
      output << i << "\t" << eig_p[0][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_tp();++i)
      output << i << "\t" << eig_p[1][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_ph();++i)
      output << i << "\t" << eig_p[2][i] << std::endl;

   return output;

}

/**
 * acces to the numbers
 * @param block < 2, call EIG_PQ::operator() , == 2 get element "index" from G block
 * @param index which element in block you want
 * @return the element on place index in block "block"
 */
double EIG_PQG::operator()(int block,int index){

   if(block < 2)
      return this->EIG_PQ::operator()(block,index);
   else
      return eig_ph[index];

}

/**
 * expands the line search function of EIG_PQ by adding the G constraint
 * @param alpha step length along the Newton direction
 * @return The line search function, gradient of the potential in the Newton direction as a function of the step length alpha
 */
double EIG_PQG::lsfunc(double alpha){

   double ward = this->EIG_PQ::lsfunc(alpha);

   for(int i = 0;i < n_ph;++i)
      ward += eig_ph[i]/(1.0 + alpha*eig_ph[i]);

   return ward;

}

/**
 * @return the minimal element present in this EIG_PQG object.
 */
double EIG_PQG::min(){

   double ward = this->EIG_PQ::min();

   if(ward > eig_ph[0])
      ward = eig_ph[0];

   return ward;

}
