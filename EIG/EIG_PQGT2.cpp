#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "../headers/EIG/EIG_PQGT2.h"
#include "../headers/SUP/SUP_PQGT2.h"
#include "../headers/lapack.h"

/**
 * standard constructor, calls the constructor of EIG_PQG and expands it by allocating the memory for the 
 * T1-space vector eig_pph.
 * @param M dimension of sp space
 * @param N nr of particles
 */
EIG_PQGT2::EIG_PQGT2(int M,int N) : EIG_PQG(M,N){

   this->n_pph = M*M*(M - 1)/2;

   dim += n_pph;

   eig_pph = new double [n_pph];

}

/**
 * standard constructor, calls the constructor of EIG_PQG and expands it by allocating the memory for the 
 * T1-space vector eig_pph. copies the content of eig_c into this
 * @param eig_c input EIG_PQGT2
 */
EIG_PQGT2::EIG_PQGT2(EIG_PQGT2 &eig_c) : EIG_PQG(eig_c){

   this->n_pph = eig_c.n_pph;

   dim += n_pph;

   eig_pph = new double [n_pph];

   int inc = 1;

   dcopy_(&n_pph,eig_c.eig_pph,&inc,eig_pph,&inc);

}

/**
 * standard constructor, with initialization on the eigenvalues of a SUP_PQGT2 object.
 * @param SZ input SUP_PQGT2 object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP_PQGT2 matrix.
 */
EIG_PQGT2::EIG_PQGT2(SUP_PQGT2 &SZ) : EIG_PQG(SZ){

   this->n_pph = SZ.gn_pph();
   dim += n_pph;

   eig_pph = new double [n_pph];

   SZ.pphm().diagonalize(eig_pph);

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG_PQGT2 &EIG_PQGT2::operator=(EIG_PQGT2 &eig_c){

   this->EIG_PQG::operator=(eig_c);

   int inc = 1;

   dcopy_(&n_pph,eig_c.eig_pph,&inc,eig_pph,&inc);

   return *this;

}

/** 
 * @param i < 3 calls EIG_PQG::operator[], == 3 returns eig_pph
 * @return array of eigenvalues
 */
double *EIG_PQGT2::operator[](int i){

   if(i < 3)
      return this->EIG_PQG::operator[](i);
   else
      return eig_pph;

}

/**
 * @return the dimension of three particle space
 */
int EIG_PQGT2::gn_pph(){

   return n_pph;

}
 
/**
 * destructor: deallocates the memory of eig_pph
 */
EIG_PQGT2::~EIG_PQGT2(){

   delete [] eig_pph;

}

ostream &operator<<(ostream &output,EIG_PQGT2 &eig_p){

   for(int i = 0;i < eig_p.gn_tp();++i)
      output << i << "\t" << eig_p[0][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_tp();++i)
      output << i << "\t" << eig_p[1][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_ph();++i)
      output << i << "\t" << eig_p[2][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_pph();++i)
      output << i << "\t" << eig_p[3][i] << std::endl;

   return output;

}

/**
 * acces to the numbers
 * @param block < 3, call EIG_PQG::operator() , == 2 get element "index" from T1 block
 * @param index which element in block you want
 * @return the element on place index in block "block"
 */
double EIG_PQGT2::operator()(int block,int index){

   if(block < 3)
      return this->EIG_PQG::operator()(block,index);
   else
      return eig_pph[index];

}

/**
 * expands the line search function of EIG_PQG by adding the T1 constraint
 * @param alpha step length along the Newton direction
 * @return The line search function, gradient of the potential in the Newton direction as a function of the step length alpha
 */
double EIG_PQGT2::lsfunc(double alpha){

   double ward = this->EIG_PQG::lsfunc(alpha);

   for(int i = 0;i < n_pph;++i)
      ward += eig_pph[i]/(1.0 + alpha*eig_pph[i]);

   return ward;

}

/**
 * @return the minimal element present in this EIG_PQGT2 object.
 */
double EIG_PQGT2::min(){

   double ward = this->EIG_PQG::min();

   if(ward > eig_pph[0])
      ward = eig_pph[0];

   return ward;

}
