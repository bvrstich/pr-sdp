#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int PPHM::counter = 0;

int **PPHM::pph2s;
int ***PPHM::s2pph;

/**
 * standard constructor: constructs Matrix object of dimension M*M*(M - 1)/2 and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and pph basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
PPHM::PPHM(int M,int N) : Matrix(M*M*(M - 1)/2) {

   this->N = N;
   this->M = M;
   this->n = M*M*(M - 1)/2;

   if(counter == 0){

      //allocatie van s2pph
      s2pph = new int ** [M];

      for(int i = 0;i < M;++i){

         s2pph[i] = new int * [M];

         for(int j = 0;j < M;++j)
            s2pph[i][j] = new int [M];

      }

      //allocatie van pph2s
      pph2s = new int * [n];

      for(int i = 0;i < n;++i)
         pph2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = 0;c < M;++c){

               s2pph[a][b][c] = teller;

               pph2s[teller][0] = a;
               pph2s[teller][1] = b;
               pph2s[teller][2] = c;

               ++teller;

            }

   }

   ++counter;

}

/**
 * copy constructor: constructs Matrix object of dimension M*M*(M - 1)/2 and copies the content of pphm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and pph basis.
 * @param pphm_c input PPHM to be copied
 */
PPHM::PPHM(const PPHM &pphm_c) : Matrix(pphm_c){

   this->N = pphm_c.N;
   this->M = pphm_c.M;
   this->n = pphm_c.n;

   if(counter == 0){

      //allocatie van s2pph
      s2pph = new int ** [M];

      for(int i = 0;i < M;++i){

         s2pph[i] = new int * [M];

         for(int j = 0;j < M;++j)
            s2pph[i][j] = new int [M];

      }

      //allocatie van pph2s
      pph2s = new int * [n];

      for(int i = 0;i < n;++i)
         pph2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = 0;c < M;++c){

               s2pph[a][b][c] = teller;

               pph2s[teller][0] = a;
               pph2s[teller][1] = b;
               pph2s[teller][2] = c;

               ++teller;

            }

   }

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists pph2s en s2pph will be deleted.
 */
PPHM::~PPHM(){

   if(counter == 1){

      for(int i = 0;i < M;++i){

         for(int j = 0;j < M;++j)
            delete [] s2pph[i][j];

         delete [] s2pph[i];

      }

      delete [] s2pph;

      for(int i = 0;i < n;++i)
         delete [] pph2s[i];

      delete [] pph2s;

   }

   --counter;

}

/**
 * access the elements of the matrix in sp mode, antisymmetry in the first two indices is automatically accounted for:\n\n
 * PPHM(a,b,c,d,e,f) = -PPHM(b,a,c,d,e,f)  but PPHM(a,b,c,d,e,f) != - PPHM(a,c,b,d,e,f) \n\n
 * PPHM(a,a,c,d,e,f) = 0 but PPHM(a,b,b,d,e,f) != 0 \n\n
 * @param a first sp index that forms the pph row index i together with b and c
 * @param b second sp index that forms the pph row index i together with a and c
 * @param c third sp index that forms the pph row index i together with a and b
 * @param d first sp index that forms the pph column index j together with e and z
 * @param e second sp index that forms the pph column index j together with d and z
 * @param z third sp index that forms the pph column index j together with d and e
 * @return the number on place PPHM(i,j) with the right phase.
 */
double PPHM::operator()(int a,int b,int c,int d,int e,int z) const{

   //eerst kijken of de eerste twee indices gelijk zijn:
   if(a == b || d == e)
      return 0;

   //dan kijken wel pph index met welke fase moet genomen worden:
   //eerst voor de i
   int i;

   int phase = 1;

   if(a > b){

      i = s2pph[b][a][c];
      phase *= -1;

   }
   else
      i = s2pph[a][b][c];

   int j;

   if(d > e){

      j = s2pph[e][d][z];
      phase *= -1;

   }
   else
      j = s2pph[d][e][z];

   return phase*(*this)(i,j);

}

ostream &operator<<(ostream &output,const PPHM &pphm_p){

   for(int i = 0;i < pphm_p.n;++i)
      for(int j = 0;j < pphm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << pphm_p.pph2s[i][0] << "\t" << pphm_p.pph2s[i][1] << "\t" << pphm_p.pph2s[i][2]

            << "\t" << pphm_p.pph2s[j][0] << "\t" << pphm_p.pph2s[j][1] << "\t" << pphm_p.pph2s[j][2] << "\t" << pphm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return nr of particles
 */
int PPHM::gN() const{

   return N;

}

/**
 * @return dimension of sp space
 */
int PPHM::gM() const{

   return M;

}

/**
 * The T2-map: maps a TPM object (tpm) on a PPHM object (*this). see primal_dual.pdf for more information
 * @param option == 0, regular T2, == 1, special (incorrect) T2, keep for test in program with regular T2
 * @param tpm input TPM
 */
void PPHM::T(int option,const TPM &tpm){

   if(option == 0){

      //construct the spm
      SPM spm(1.0/(N - 1.0),tpm);

      int a,b,c,d,e,z;

      for(int i = 0;i < n;++i){

         a = pph2s[i][0];
         b = pph2s[i][1];
         c = pph2s[i][2];

         for(int j = i;j < n;++j){

            d = pph2s[j][0];
            e = pph2s[j][1];
            z = pph2s[j][2];

            //initialize
            (*this)(i,j) = 0.0;

            if(a == d){

               //sp part
               if(b == e)
                  (*this)(i,j) += spm(c,z);

               //tp part
               (*this)(i,j) -= tpm(c,e,z,b);

            }

            //now only tp parts left:
            if(c == z)
               (*this)(i,j) += tpm(a,b,d,e);

            if(b == d)
               (*this)(i,j) += tpm(c,e,z,a);

            if(b == e)
               (*this)(i,j) -= tpm(c,d,z,a);

         }
      }
   }
   else if(option == 1){

      TPM Q(M,N);
      Q.Q(1,tpm);

      DPM T(M,N);
      T.T(1,tpm);

      int a,b,c,d,e,z;

      for(int i = 0;i < n;++i){

         a = pph2s[i][0];
         b = pph2s[i][1];
         c = pph2s[i][2];

         for(int j = i;j < n;++j){

            d = pph2s[j][0];
            e = pph2s[j][1];
            z = pph2s[j][2];

            //initialize
            (*this)(i,j) = 0.0;

            if(c == z)
               (*this)(i,j) += tpm(a,b,d,e) + Q(a,b,d,e);

            (*this)(i,j) -= T(a,b,z,d,e,c);

         }
      }

   }
   else{

      
      //construct the spm
      SPM spm(1.0/(N - 1.0),tpm);

      int a,b,c,d,e,z;

      for(int i = 0;i < n;++i){

         a = pph2s[i][0];
         b = pph2s[i][1];
         c = pph2s[i][2];

         for(int j = i;j < n;++j){

            d = pph2s[j][0];
            e = pph2s[j][1];
            z = pph2s[j][2];

            //initialize
            (*this)(i,j) = 0.0;

            if(a == d){

               //sp part
               if(b == e)
                  (*this)(i,j) += spm(c,z);

               //tp part
               (*this)(i,j) -= tpm(c,e,z,b);

            }

            //now only tp parts left:
            if(c == z)
               (*this)(i,j) += tpm(a,b,d,e);

            if(b == d)
               (*this)(i,j) += tpm(c,e,z,a);

            if(b == e)
               (*this)(i,j) -= tpm(c,d,z,a);

         }
      }

   }

   //and symmetrize
   this->symmetrize();

}

/** 
 * yet another way to avoid the very slow access operator with the ordering
 */
void PPHM::convert(double *ppharray) const{

   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b)
         for(int c = 0;c < M;++c)
            for(int d = 0;d < M;++d)
               for(int e = 0;e < M;++e)
                  for(int z = 0;z < M;++z)
                     ppharray[a*M5 + b*M4 + c*M3 + d*M2 + e*M + z] = (*this)(a,b,c,d,e,z);


}
