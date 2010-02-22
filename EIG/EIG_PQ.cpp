#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "../headers/EIG/EIG_PQ.h"
#include "../headers/SUP/SUP_PQ.h"
#include "../headers/lapack.h"

//constructor:
EIG_PQ::EIG_PQ(int M,int N){
   
   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

   this->dim = 2*n_tp;

   eig = new double * [2];

   eig[0] = new double [dim];
   eig[1] = eig[0] + n_tp;

}

//copy constructor:
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

//constructor met initialisatie door SUP matrix:
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

EIG_PQ &EIG_PQ::operator=(EIG_PQ &eig_c){

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

   return *this;

}

double *EIG_PQ::operator[](int i){

   return eig[i];

}

int EIG_PQ::gN(){

   return N;

}

int EIG_PQ::gM(){

   return M;

}

int EIG_PQ::gn_tp(){

   return n_tp;

}

//destructor
EIG_PQ::~EIG_PQ(){

   delete [] eig[0];
   delete [] eig;

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,EIG_PQ &eig_p){

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[0][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.n_tp;++i)
      output << i << "\t" << eig_p.eig[1][i] << std::endl;

   return output;

}

double EIG_PQ::operator()(int block,int index){

   return eig[block][index];

}

double EIG_PQ::lsfunc(double alpha){

   double ward = 0.0;

   for(int i = 0;i < n_tp;++i)
      ward += eig[0][i]/(1.0 + alpha*eig[0][i]);

   for(int i = 0;i < n_tp;++i)
      ward += eig[1][i]/(1.0 + alpha*eig[1][i]);

   return ward;

}

double EIG_PQ::min(){

   double ward;

   if(eig[0][0] < eig[1][0])
      ward = eig[0][0];
   else
      ward = eig[1][0];

   return ward;

}
