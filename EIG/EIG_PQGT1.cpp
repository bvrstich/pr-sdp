#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "../headers/EIG/EIG_PQGT1.h"
#include "../headers/SUP/SUP_PQGT1.h"
#include "../headers/lapack.h"

//constructor:
EIG_PQGT1::EIG_PQGT1(int M,int N) : EIG_PQG(M,N){
   
   this->n_dp = M*(M - 1)*(M - 2)/6;

   dim += n_dp;

   eig_dp = new double [n_dp];

}

//copy constructor:
EIG_PQGT1::EIG_PQGT1(EIG_PQGT1 &eig_c) : EIG_PQG(eig_c){

   this->n_dp = eig_c.n_dp;

   dim += n_dp;

   eig_dp = new double [n_dp];

   int inc = 1;

   dcopy_(&n_dp,eig_c.eig_dp,&inc,eig_dp,&inc);

}

//constructor met initialisatie door SUP matrix:
EIG_PQGT1::EIG_PQGT1(SUP_PQGT1 &SZ) : EIG_PQG(SZ){

   this->n_dp = SZ.gn_dp();
   dim += n_dp;

   eig_dp = new double [n_dp];

   SZ.dpm().diagonalize(eig_dp);

}

EIG_PQGT1 &EIG_PQGT1::operator=(EIG_PQGT1 &eig_c){

   this->EIG_PQG::operator=(eig_c);

   int inc = 1;

   dcopy_(&n_dp,eig_c.eig_dp,&inc,eig_dp,&inc);

   return *this;

}

double *EIG_PQGT1::operator[](int i){

   if(i < 3)
      return this->EIG_PQG::operator[](i);
   else
      return eig_dp;

}

int EIG_PQGT1::gn_dp(){

   return n_dp;

}
 
//destructor
EIG_PQGT1::~EIG_PQGT1(){

   delete [] eig_dp;

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,EIG_PQGT1 &eig_p){

   for(int i = 0;i < eig_p.gn_tp();++i)
      output << i << "\t" << eig_p[0][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_tp();++i)
      output << i << "\t" << eig_p[1][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_ph();++i)
      output << i << "\t" << eig_p[2][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_dp();++i)
      output << i << "\t" << eig_p[3][i] << std::endl;

   return output;

}

double EIG_PQGT1::operator()(int block,int index){

   if(block < 3)
      return this->EIG_PQG::operator()(block,index);
   else
      return eig_dp[index];

}

double EIG_PQGT1::lsfunc(double alpha){

   double ward = this->EIG_PQG::lsfunc(alpha);

   for(int i = 0;i < n_dp;++i)
      ward += eig_dp[i]/(1.0 + alpha*eig_dp[i]);

   return ward;

}

double EIG_PQGT1::min(){

   double ward = this->EIG_PQG::min();

   if(ward > eig_dp[0])
      ward = eig_dp[0];

   return ward;

}
