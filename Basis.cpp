#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * constructor 
 */
Basis::Basis(int M,int N){

   this->M = M;
   this->N = N;

   basis = new TPM * [TPTPM::gn()];

   int I,J;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      basis[i] = new TPM(M,N);

      *basis[i] = 0.0;

      if(I == J)
         (*basis[i])(I,I) = 1.0;
      else
         (*basis[i])(I,J) = (*basis[i])(J,I) = 1.0/std::sqrt(2.0);

   }

}

/**
 * copy constructor 
 */
Basis::Basis(const Basis &basis_c){

   this->M = basis_c.M;
   this->N = basis_c.N;

   basis = new TPM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i)
      basis[i] = new TPM(basis_c[i]);

}

/**
 * Destructor
 */
Basis::~Basis(){

   for(int i = 0;i < TPTPM::gn();++i)
      delete basis[i];

   delete [] basis;

}

TPM &Basis::operator[](int i){

   return *basis[i];

}

const TPM &Basis::operator[](int i) const{

   return *basis[i];

}

ostream &operator<<(ostream &output,const Basis &basis_p){

   int I,J;

   for(int i = 0;i < TPTPM::gn();++i){
      output << "FUCK YOU" << endl;

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      output << endl;
      output << "basismatrix\t" << i << "\t|\t" << I << "\t" << J << endl;
      output << endl;
      output << basis_p[i] << endl;

   }

   return output;

}
