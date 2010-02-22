#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::endl;

#include "headers/include.h"

int DPM::counter = 0;

int **DPM::dp2s;
int ***DPM::s2dp;

//constructor:
DPM::DPM(int M,int N) : Matrix(M*(M - 1)*(M - 2)/6) {
   
   this->N = N;
   this->M = M;
   this->n = M*(M - 1)*(M - 2)/6;

   if(counter == 0){

      //allocatie van s2dp
      s2dp = new int ** [M];

      for(int i = 0;i < M;++i){

         s2dp[i] = new int * [M];

         for(int j = 0;j < M;++j)
            s2dp[i][j] = new int [M];

      }

      //allocatie van dp2s
      dp2s = new int * [n];

      for(int i = 0;i < n;++i)
         dp2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = b + 1;c < M;++c){

               s2dp[a][b][c] = teller;

               dp2s[teller][0] = a;
               dp2s[teller][1] = b;
               dp2s[teller][2] = c;

               ++teller;

            }

   }

   ++counter;

}

//copy constructor
DPM::DPM(DPM &dpm_c) : Matrix(dpm_c){

   this->N = dpm_c.N;
   this->M = dpm_c.M;
   this->n = M*(M - 1)*(M - 2)/6;

   if(counter == 0){

      //allocatie van s2dp
      s2dp = new int ** [M];

      for(int i = 0;i < M;++i){

         s2dp[i] = new int * [M];

         for(int j = 0;j < M;++j)
            s2dp[i][j] = new int [M];

      }

      //allocatie van dp2s
      dp2s = new int * [n];

      for(int i = 0;i < n;++i)
         dp2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = b + 1;c < M;++c){

               s2dp[a][b][c] = teller;

               dp2s[teller][0] = a;
               dp2s[teller][1] = b;
               dp2s[teller][2] = c;

               ++teller;

            }

   }

   ++counter;

}

//destructor
DPM::~DPM(){

   if(counter == 1){

      for(int i = 0;i < M;++i){

         for(int j = 0;j < M;++j)
            delete [] s2dp[i][j];

         delete [] s2dp[i];

      }

      delete [] s2dp;

      for(int i = 0;i < n;++i)
         delete [] dp2s[i];

      delete [] dp2s;

   }

   --counter;

}

//access the numbers: sp indices
double DPM::operator()(int a,int b,int c,int d,int e,int z) const{

   //eerst kijken of er geen indices gelijk zijn:
   if(a == b || a == c || b == c)
      return 0;

   if(d == e || d == z || e == z)
      return 0;

   //dan kijken wel dp index met welke fase moet genomen worden:
   //eerst voor de i
   int i;

   int phase = 1;

   if(a < b){

      if(b < c)
         i = s2dp[a][b][c];
      else if(c < a)
         i = s2dp[c][a][b];
      else{

         i = s2dp[a][c][b];
         phase *= -1;

      }

   }
   else{

      if(a < c){

         i = s2dp[b][a][c];
         phase *= -1;

      }
      else if(c < b){

         i = s2dp[c][b][a];
         phase *= -1;

      }
      else
         i = s2dp[b][c][a];

   }

   //idem voor j maar met d e z
   int j;

   if(d < e){

      if(e < z)
         j = s2dp[d][e][z];
      else if(z < d)
         j = s2dp[z][d][e];
      else{

         j = s2dp[d][z][e];
         phase *= -1;

      }

   }
   else{

      if(d < z){

         j = s2dp[e][d][z];
         phase *= -1;

      }
      else if(z < e){

         j = s2dp[z][e][d];
         phase *= -1;

      }
      else
         j = s2dp[e][z][d];

   }

   return phase*(*this)(i,j);

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,const DPM &dpm_p){

   for(int i = 0;i < dpm_p.n;++i)
      for(int j = 0;j < dpm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << dpm_p.dp2s[i][0] << "\t" << dpm_p.dp2s[i][1] << "\t" << dpm_p.dp2s[i][2]

            << "\t" << dpm_p.dp2s[j][0] << "\t" << dpm_p.dp2s[j][1] << "\t" << dpm_p.dp2s[j][2] << "\t" << dpm_p(i,j) << endl;

      }

   return output;

}

int DPM::gN(){

   return N;

}

int DPM::gM(){

   return M;

}

int DPM::gn(){

   return n;

}

void DPM::T(TPM &tpm){

   SPM spm(tpm);

   double ward = 2.0*tpm.trace()/(N*(N - 1.0));

   int a,b,c,d,e,z;

   for(int i = 0;i < n;++i){

      a = dp2s[i][0];
      b = dp2s[i][1];
      c = dp2s[i][2];

      for(int j = i;j < n;++j){

         d = dp2s[j][0];
         e = dp2s[j][1];
         z = dp2s[j][2];

         (*this)(i,j) = 0;

         //no particle stuk
         if(i == j)
            (*this)(i,i) = ward;

         if(a == d){

            //tp stuk
            (*this)(i,j) += tpm(b,c,e,z);

            //4 sp stukken
            if(c == z)
               (*this)(i,j) -= spm(e,b);

            if(b == z)
               (*this)(i,j) += spm(c,e);

            if(c == e)
               (*this)(i,j) += spm(b,z);

            if(b == e)
               (*this)(i,j) -= spm(c,z);

         }

         if(b == d){

            //tp stuk
            (*this)(i,j) -= tpm(a,c,e,z);

            //2 sp stukken
            if(c == z)
               (*this)(i,j) += spm(a,e);

            if(c == e)
               (*this)(i,j) -= spm(a,z);

         }

         if(b == e){

            //tp stuk
            (*this)(i,j) += tpm(a,c,d,z);

            //sp stuk
            if(c == z)
               (*this)(i,j) -= spm(a,d);

         }

         //nu enkel nog tp stukken
         if(c == z)
            (*this)(i,j) += tpm(a,b,d,e);

         if(b == z)
            (*this)(i,j) -= tpm(a,c,d,e);

         if(c == e)
            (*this)(i,j) -= tpm(a,b,d,z);

         if(c == d)
            (*this)(i,j) += tpm(a,b,e,z);

      }
   }

   //niet vergeten!
   this->symmetrize();

}
