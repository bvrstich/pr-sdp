#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "headers/SUP_PQG.h"
#include "headers/EIG.h"

//constructor
SUP_PQG::SUP_PQG(int M,int N) : SUP_PQ(M,N){

   this->n_ph = M*M;

   dim += n_ph;

   SZ_ph = new PHM(M,N);

}

//copy constructor
SUP_PQG::SUP_PQG(SUP_PQG &SZ_c) : SUP_PQ(SZ_c){

   this->n_ph = SZ_c.n_ph;

   dim += n_ph;

   SZ_ph = new PHM(M,N);

   (*SZ_ph) = (*SZ_c.SZ_ph);

}

//destructor
SUP_PQG::~SUP_PQG(){

   delete SZ_ph;

}

SUP_PQG &SUP_PQG::operator+=(SUP_PQG &SZ_pl){

   this->SUP_PQ::operator+=(SZ_pl);

   (*SZ_ph) += (*SZ_pl.SZ_ph);

   return *this;

}

SUP_PQG &SUP_PQG::operator-=(SUP_PQG &SZ_pl){

   this->SUP_PQ::operator-=(SZ_pl);

   (*SZ_ph) -= (*SZ_pl.SZ_ph);

   return *this;

}

//overload equality operator
SUP_PQG &SUP_PQG::operator=(SUP_PQG &SZ_c){

   this->SUP_PQ::operator=(SZ_c);

   (*SZ_ph) = (*SZ_c.SZ_ph);

   return *this;

}

SUP_PQG &SUP_PQG::operator=(double &a){

   this->SUP_PQ::operator=(a);

   (*SZ_ph) = a;

   return *this;

}

int SUP_PQG::gn_ph(){

   return n_ph;

}

PHM &SUP_PQG::phm(){

   return *SZ_ph;

}

double SUP_PQG::ddot(SUP_PQG &SZ_i){

   double ward = this->SUP_PQ::ddot(SZ_i);

   ward += SZ_ph->ddot(*SZ_i.SZ_ph);

   return ward;

}

void SUP_PQG::invert(){

   this->SUP_PQ::invert();

   SZ_ph->invert();

}

void SUP_PQG::dscal(double alpha){

   this->SUP_PQ::dscal(alpha);

   SZ_ph->dscal(alpha);

}

void SUP_PQG::sqrt(int option){

   this->SUP_PQ::sqrt(option);

   SZ_ph->sqrt(option);

}

void SUP_PQG::L_map(SUP_PQG &map,SUP_PQG &object){

   this->SUP_PQ::L_map(map,object);

   SZ_ph->L_map(*map.SZ_ph,*object.SZ_ph);

}

void SUP_PQG::daxpy(double alpha,SUP_PQG &SZ_p){

   this->SUP_PQ::daxpy(alpha,SZ_p);

   SZ_ph->daxpy(alpha,*SZ_p.SZ_ph);

}

double SUP_PQG::trace(){

   double ward = this->SUP_PQ::trace();

   ward += SZ_ph->trace();

   return ward;

}

void SUP_PQG::fill(TPM &tpm){

   this->SUP_PQ::fill(tpm);

   SZ_ph->G(tpm);

}
