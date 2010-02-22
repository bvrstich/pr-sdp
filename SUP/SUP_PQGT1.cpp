#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "../headers/SUP/SUP_PQGT1.h"
#include "../headers/EIG/EIG_PQGT1.h"

//constructor
SUP_PQGT1::SUP_PQGT1(int M,int N) : SUP_PQG(M,N){

   this->n_dp = M*(M - 1)*(M - 2)/6;

   dim += n_dp;

   SZ_dp = new DPM(M,N);

}

//copy constructor
SUP_PQGT1::SUP_PQGT1(SUP_PQGT1 &SZ_c) : SUP_PQG(SZ_c){

   this->n_dp = SZ_c.n_dp;

   dim += n_dp;

   SZ_dp = new DPM(M,N);

   (*SZ_dp) = (*SZ_c.SZ_dp);

}

//destructor
SUP_PQGT1::~SUP_PQGT1(){

   delete SZ_dp;

}

SUP_PQGT1 &SUP_PQGT1::operator+=(SUP_PQGT1 &SZ_pl){

   this->SUP_PQG::operator+=(SZ_pl);

   (*SZ_dp) += (*SZ_pl.SZ_dp);

   return *this;

}

SUP_PQGT1 &SUP_PQGT1::operator-=(SUP_PQGT1 &SZ_pl){

   this->SUP_PQG::operator-=(SZ_pl);

   (*SZ_dp) -= (*SZ_pl.SZ_dp);

   return *this;

}

//overload equality operator
SUP_PQGT1 &SUP_PQGT1::operator=(SUP_PQGT1 &SZ_c){

   this->SUP_PQG::operator=(SZ_c);

   (*SZ_dp) = (*SZ_c.SZ_dp);

   return *this;

}

SUP_PQGT1 &SUP_PQGT1::operator=(double &a){

   this->SUP_PQG::operator=(a);

   (*SZ_dp) = a;

   return *this;

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,SUP_PQGT1 &SZ_p){

   for(int i = 0;i < 2;++i)
      output << SZ_p.tpm(i) << std::endl;

   output << SZ_p.phm() << std::endl;

   output << SZ_p.dpm() << std::endl;

   return output;

}

int SUP_PQGT1::gn_dp(){

   return n_dp;

}

DPM &SUP_PQGT1::dpm(){

   return *SZ_dp;

}

double SUP_PQGT1::ddot(SUP_PQGT1 &SZ_i){

   double ward = this->SUP_PQG::ddot(SZ_i);

   ward += SZ_dp->ddot(SZ_i.dpm());

   return ward;

}

void SUP_PQGT1::invert(){

   this->SUP_PQG::invert();

   SZ_dp->invert();

}

void SUP_PQGT1::dscal(double alpha){

   this->SUP_PQG::dscal(alpha);

   SZ_dp->dscal(alpha);

}

void SUP_PQGT1::sqrt(int option){

   this->SUP_PQG::sqrt(option);

   SZ_dp->sqrt(option);

}

void SUP_PQGT1::L_map(SUP_PQGT1 &map,SUP_PQGT1 &object){

   this->SUP_PQG::L_map(map,object);

   SZ_dp->L_map(*map.SZ_dp,*object.SZ_dp);

}

void SUP_PQGT1::daxpy(double alpha,SUP_PQGT1 &SZ_p){

   this->SUP_PQG::daxpy(alpha,SZ_p);

   SZ_dp->daxpy(alpha,*SZ_p.SZ_dp);

}

double SUP_PQGT1::trace(){

   double ward = this->SUP_PQG::trace();

   ward += SZ_dp->trace();

   return ward;

}

void SUP_PQGT1::fill(TPM &tpm){

   this->SUP_PQG::fill(tpm);

   SZ_dp->T(tpm);

}
