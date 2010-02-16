#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::endl;

#include "headers/PHM.h"
#include "headers/SPM.h"
#include "headers/lapack.h"

int PHM::counter = 0;

int **PHM::ph2s;
int **PHM::s2ph;

//constructor:
PHM::PHM(int M,int N) : Matrix(M*M) {
   
   this->N = N;
   this->M = M;
   this->n = M*M;

   if(counter == 0){

      //allocatie van s2ph
      s2ph = new int * [M];
      s2ph[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2ph[i] = s2ph[i - 1] + M;

      //allocatie van ph2s
      ph2s = new int * [n];

      for(int i = 0;i < n;++i)
         ph2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            s2ph[a][b] = teller;

            ph2s[teller][0] = a;
            ph2s[teller][1] = b;

            ++teller;

         }

   }

   ++counter;

}

//copy constructor
PHM::PHM(PHM &phm_c) : Matrix(phm_c){

   this->N = phm_c.N;
   this->M = phm_c.M;
   this->n = M*M;

   if(counter == 0){

      //allocatie van sp2tp
      s2ph = new int * [M];
      s2ph[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2ph[i] = s2ph[i - 1] + M;

      //allocatie van tp2sp
      ph2s = new int * [n];

      for(int i = 0;i < n;++i)
         ph2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            s2ph[a][b] = teller;

            ph2s[teller][0] = a;
            ph2s[teller][1] = b;

            ++teller;

         }

   }

   ++counter;

}

//destructor
PHM::~PHM(){

   if(counter == 1){

      delete [] s2ph[0];
      delete [] s2ph;

      for(int i = 0;i < n;++i)
         delete [] ph2s[i];

      delete [] ph2s;

   }

   --counter;

}

//access the numbers: sp indices
double PHM::operator()(int a,int b,int c,int d) const{

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(i,j);

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,const PHM &phm_p){

   for(int i = 0;i < phm_p.n;++i)
      for(int j = 0;j < phm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << phm_p.ph2s[i][0] << "\t" << phm_p.ph2s[i][1]

            << "\t" << phm_p.ph2s[j][0] << "\t" << phm_p.ph2s[j][1] << "\t" << phm_p(i,j) << endl;

      }

   return output;

}

int PHM::gN(){

   return N;

}

int PHM::gM(){

   return M;

}

int PHM::gn(){

   return n;

}

void PHM::G(TPM &tpm){

   SPM spm(tpm);

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < n;++j){

         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = -tpm(a,d,c,b);

         if(b == d)
            (*this)(i,j) += spm(a,c);

      }
   }

   this->symmetrize();

}
