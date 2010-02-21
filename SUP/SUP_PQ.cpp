#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "../headers/SUP/SUP_PQ.h"
#include "../headers/EIG.h"

//constructor
SUP_PQ::SUP_PQ(int M,int N){

   this->M = M;
   this->N = N;
   this->n_tp = M*(M - 1)/2;

   dim = 2*n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

}

//copy constructor
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

//destructor
SUP_PQ::~SUP_PQ(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

}

SUP_PQ &SUP_PQ::operator+=(SUP_PQ &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

   return *this;

}

SUP_PQ &SUP_PQ::operator-=(SUP_PQ &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

   return *this;

}

//overload equality operator
SUP_PQ &SUP_PQ::operator=(SUP_PQ &SZ_c){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) = (*SZ_c.SZ_tp[i]);

   return *this;

}

SUP_PQ &SUP_PQ::operator=(double &a){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) = a;

   return *this;

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,const SUP_PQ &SZ_p){

   for(int i = 0;i < 2;++i)
      output << (*SZ_p.SZ_tp[i]) << std::endl;

   return output;

}

int SUP_PQ::gN(){

   return N;

}

int SUP_PQ::gM(){

   return M;

}

int SUP_PQ::gn_tp(){

   return n_tp;

}

TPM &SUP_PQ::tpm(int i){

   return *SZ_tp[i];

}

double SUP_PQ::ddot(SUP_PQ &SZ_i){

   double ward = 0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

   return ward;

}

void SUP_PQ::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

}

void SUP_PQ::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

}

void SUP_PQ::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

}

void SUP_PQ::L_map(SUP_PQ &map,SUP_PQ &object){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->L_map(*map.SZ_tp[i],*object.SZ_tp[i]);

}

void SUP_PQ::daxpy(double alpha,SUP_PQ &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,*SZ_p.SZ_tp[i]);

}

double SUP_PQ::trace(){

   double ward = 0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->trace();

   return ward;

}

void SUP_PQ::fill(TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(tpm);

}
