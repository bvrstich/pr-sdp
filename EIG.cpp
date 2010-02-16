#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "headers/EIG.h"
#include "headers/SUP.h"
#include "headers/lapack.h"

//constructor:
EIG::EIG(int M,int N){
   
   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

#ifdef PQ

   this->n = 2*n_tp;

   eig = new double * [2];

   eig[0] = new double [n];
   eig[1] = eig[0] + n_tp;

#else

   this->n_ph = M*M;

   n = 2*n_tp + n_ph;

   eig = new double * [3];

   eig[0] = new double [n];

   for(int i = 1;i < 3;++i)
      eig[i] = eig[i - 1] + n_tp;

#endif

}

//copy constructor:
EIG::EIG(EIG &eig_c){

   this->N = eig_c.N;
   this->M = eig_c.M;
   this->n_tp = eig_c.n_tp;

#ifdef PQ

   this->n = 2*n_tp;

   eig = new double * [2];

   eig[0] = new double [n];
   eig[1] = eig[0] + n_tp;

#else

   this->n_ph = eig_c.n_ph;

   n = 2*n_tp + n_ph;

   eig = new double * [3];

   eig[0] = new double [n];

   for(int i = 1;i < 3;++i)
      eig[i] = eig[i - 1] + n_tp;

#endif

   int inc = 1;

   dcopy_(&n,eig_c.eig[0],&inc,eig[0],&inc);

}

//constructor met initialisatie door SUP matrix:
EIG::EIG(SUP &SZ){

   this->N = SZ.gN();
   this->M = SZ.gM();
   this->n_tp = SZ.gn_tp();

#ifdef PQ

   this->n = 2*n_tp;

   eig = new double * [2];

   eig[0] = new double [n];
   eig[1] = eig[0] + n_tp;

#else

   this->n_ph = SZ.gn_ph();
   this->n = 2*n_tp + n_ph;

   eig = new double * [3];

   eig[0] = new double [n];

   eig[1] = eig[0] + n_tp;
   eig[2] = eig[1] + n_tp;

#endif

   //dit vernietigd de originele matrix!
   for(int i = 0;i < 2;++i)
      SZ.tpm(i).diagonalize(eig[i]);

#ifndef PQ

   SZ.phm().diagonalize(eig[2]);

#endif

}

EIG &EIG::operator=(EIG &eig_c){

   int inc = 1;

   dcopy_(&n,eig_c.eig[0],&inc,eig[0],&inc);

   return *this;

}

double *EIG::operator[](int i){

   return eig[i];

}

int EIG::gN(){

   return N;

}

int EIG::gM(){

   return M;

}

int EIG::gn_tp(){

   return n_tp;

}

#ifndef PQ

int EIG::gn_ph(){

   return n_ph;

}
 
#endif

//destructor
EIG::~EIG(){

   delete [] eig[0];
   delete [] eig;

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,const EIG &eig_p){

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[0][i] << std::endl;

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[1][i] << std::endl;

#ifndef PQ

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_ph;++i)
      std::cout << i << "\t" << eig_p.eig[2][i] << std::endl;

#endif

   return output;

}

double &EIG::operator()(int block,int index){

   return eig[block][index];

}

double EIG::lsfunc(double alpha){

   double ward = 0.0;

   for(int i = 0;i < n_tp;++i)
      ward += eig[0][i]/(1.0 + alpha*eig[0][i]);

   for(int i = 0;i < n_tp;++i)
      ward += eig[1][i]/(1.0 + alpha*eig[1][i]);

#ifndef PQ

   for(int i = 0;i < n_ph;++i)
      ward += eig[2][i]/(1.0 + alpha*eig[2][i]);

#endif

   return ward;

}

double EIG::min(){

   double ward;

   if(eig[0][0] < eig[1][0])
      ward = eig[0][0];
   else
      ward = eig[1][0];

#ifndef PQ

   if(ward > eig[2][0])
      ward = eig[2][0];

#endif

   return ward;

}
