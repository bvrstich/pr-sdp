#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ifstream;
using std::endl;

#include "include.h"

int cPHM::counter = 0;

int **cPHM::ph2s;
int **cPHM::s2ph;

/**
 * standard constructor: constructs Matrix object of dimension D*D and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 */
cPHM::cPHM(int D) : Matrix(D*D) {
   
   this->D = D;

   if(counter == 0){

      //allocatie van s2ph
      s2ph = new int * [D];
      s2ph[0] = new int [D*D];

      for(int i = 1;i < D;++i)
         s2ph[i] = s2ph[i - 1] + D;

      //allocatie van ph2s
      ph2s = new int * [D*D];

      for(int i = 0;i < D*D;++i)
         ph2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < D;++a)
         for(int b = 0;b < D;++b){

            s2ph[a][b] = teller;

            ph2s[teller][0] = a;
            ph2s[teller][1] = b;

            ++teller;

         }

   }

   ++counter;

}

/**
 * copy constructor: constructs Matrix object of dimension M*M and copies the content of cphm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 * @param cphm_c cPHM to be copied into (*this)
 */
cPHM::cPHM(const cPHM &cphm_c) : Matrix(cphm_c){

   this->D = cphm_c.gD();

   if(counter == 0){

      //allocatie van sp2tp
      s2ph = new int * [D];
      s2ph[0] = new int [D*D];

      for(int i = 1;i < D;++i)
         s2ph[i] = s2ph[i - 1] + D;

      //allocatie van tp2sp
      ph2s = new int * [D*D];

      for(int i = 0;i < D*D;++i)
         ph2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < D;++a)
         for(int b = 0;b < D;++b){

            s2ph[a][b] = teller;

            ph2s[teller][0] = a;
            ph2s[teller][1] = b;

            ++teller;

         }

   }

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists ph2s en s2ph will be deleted.
 */
cPHM::~cPHM(){

   if(counter == 1){

      delete [] s2ph[0];
      delete [] s2ph;

      for(int i = 0;i < D*D;++i)
         delete [] ph2s[i];

      delete [] ph2s;

   }

   --counter;

}

/**
 * access the elements of the matrix in sp mode, 
 * @param a first sp index that forms the ph row index i together with b
 * @param b second sp index that forms the ph row index i together with a
 * @param c first sp index that forms the ph column index j together with d
 * @param d second sp index that forms the ph column index j together with c
 * @return the number on place cPHM(i,j)
 */
double cPHM::operator()(int a,int b,int c,int d) const{

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,const cPHM &cphm_p){

   for(int i = 0;i < cphm_p.gn();++i)
      for(int j = 0;j < cphm_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << cphm_p.ph2s[i][0] << "\t" << cphm_p.ph2s[i][1]

            << "\t" << cphm_p.ph2s[j][0] << "\t" << cphm_p.ph2s[j][1] << "\t" << cphm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return number of single particle oribals
 */
int cPHM::gD() const{

   return D;

}

void cPHM::fill(const PHM &phm){

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < gn();++j){

         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = phm(a,b,c,d);

      }
   }

   this->symmetrize();

}


