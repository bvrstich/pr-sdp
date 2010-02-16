#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "SUP.h"
#include "EIG.h"

//constructor
SUP::SUP(int M,int N){

   this->M = M;
   this->N = N;
   this->n_tp = M*(M - 1)/2;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

#ifndef PQ

   this->n_ph = M*M;

   SZ_ph = new PHM(M,N);

#endif

#ifdef T_1
   
   this->n_dp = M*(M - 1)*(M - 2)/6;

   SZ_dp = new DPM(M,N);

#endif

}

//copy constructor
SUP::SUP(SUP &SZ_c){

   this->M = SZ_c.M;
   this->N = SZ_c.N;
   this->n_tp = SZ_c.n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) = (*SZ_c.SZ_tp[i]);

#ifndef PQ

   this->n_ph = SZ_c.n_ph;

   SZ_ph = new PHM(M,N);

   (*SZ_ph) = (*SZ_c.SZ_ph);

#endif

#ifdef T_1

   this->n_dp = SZ_c.n_ph;

   SZ_dp = new DPM(M,N);

   (*SZ_dp) = (*SZ_c.SZ_dp);

#endif

}

//destructor
SUP::~SUP(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

#ifndef PQ

   delete SZ_ph;

#endif

#ifdef T_1

   delete SZ_dp;

#endif

}

SUP &SUP::operator+=(SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

#ifndef PQ

   (*SZ_ph) += (*SZ_pl.SZ_ph);

#endif

#ifdef T_1

   (*SZ_dp) += (*SZ_pl.SZ_dp);

#endif

   return *this;

}

SUP &SUP::operator-=(SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

#ifndef PQ

   (*SZ_ph) -= (*SZ_pl.SZ_ph);

#endif

#ifdef T_1

   (*SZ_dp) -= (*SZ_pl.SZ_dp);

#endif

   return *this;

}

//overload equality operator
SUP &SUP::operator=(SUP &SZ_c){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) = (*SZ_c.SZ_tp[i]);

#ifndef PQ

   (*SZ_ph) = (*SZ_c.SZ_ph);

#endif

#ifdef T_1

   (*SZ_dp) = (*SZ_pl.SZ_dp);

#endif

   return *this;

}

SUP &SUP::operator=(double &a){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) = a;

#ifndef PQ

   (*SZ_ph) = a;

#endif

#ifdef T_1
   
   (*SZ_dp) = a;

#endif

   return *this;

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,const SUP &SZ_p){

   for(int i = 0;i < 2;++i)
      output << (*SZ_p.SZ_tp[i]) << std::endl;

#ifndef PQ

   output << (*SZ_p.SZ_ph) << std::endl;

#endif

#ifdef T_1
   
   output << (*SZ_p.SZ_dp) << std::endl;

#endif

   return output;

}

int SUP::gN(){

   return N;

}

int SUP::gM(){

   return M;

}

int SUP::gn_tp(){

   return n_tp;

}

#ifndef PQ

int SUP::gn_ph(){

   return n_ph;

}

PHM &SUP::phm(){

   return *SZ_ph;

}

#endif

#ifdef T_1

int SUP::gn_dp(){

   return n_dp;

}

DPM &SUP::dpm(){

   return *SZ_dp;

}

#endif

TPM &SUP::tpm(int i){

   return *SZ_tp[i];

}

double SUP::ddot(SUP &SZ_i){

   double ward = 0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

#ifndef PQ

   ward += SZ_ph->ddot(*SZ_i.SZ_ph);

#endif

   return ward;

}

void SUP::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

#ifndef PQ

   SZ_ph->invert();

#endif

}

void SUP::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

#ifndef PQ

   SZ_ph->dscal(alpha);

#endif

}

void SUP::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

#ifndef PQ

   SZ_ph->sqrt(option);

#endif

}

void SUP::L_map(SUP &map,SUP &object){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->L_map(*map.SZ_tp[i],*object.SZ_tp[i]);

#ifndef PQ

   SZ_ph->L_map(*map.SZ_ph,*object.SZ_ph);

#endif

}

void SUP::daxpy(double alpha,SUP &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,*SZ_p.SZ_tp[i]);

#ifndef PQ

   SZ_ph->daxpy(alpha,*SZ_p.SZ_ph);

#endif

}

double SUP::trace(){

   double ward = 0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->trace();

#ifndef PQ

   ward += SZ_ph->trace();

#endif

   return ward;

}

void SUP::fill(TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(tpm);

#ifndef PQ

   SZ_ph->G(tpm);

#endif

}
