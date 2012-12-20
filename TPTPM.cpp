#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

int **TPTPM::t2tpmm;
vector< vector<int> > TPTPM::tpmm2t;

int TPTPM::M;
int TPTPM::N;

/**
 * initialize the static lists
 */
void TPTPM::init(int M_in,int N_in){

   M = M_in;
   N = N_in;

   t2tpmm = new int * [TPM::gn()];

   for(int i = 0;i < TPM::gn();++i)
      t2tpmm[i] = new int [TPM::gn()];

   vector<int> v(2);

   int tpmm = 0;

   for(int i = 0;i < TPM::gn();++i)
      for(int j = i;j < TPM::gn();++j){

         v[0] = i;
         v[1] = j;

         tpmm2t.push_back(v);

         t2tpmm[i][j] = tpmm;
         t2tpmm[j][i] = tpmm;

         ++tpmm;

      }

}

/**
 * deallocate the static lists
 */
void TPTPM::clear(){

   for(int i = 0;i < TPM::gn();++i)
      delete [] t2tpmm[i];

   delete [] t2tpmm;

}

/**
 * standard constructor:
 */
TPTPM::TPTPM() : Matrix(tpmm2t.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param tpmm_c object that will be copied into this.
 */
TPTPM::TPTPM(const TPTPM &tpmm_c) : Matrix(tpmm_c){ }

/**
 * destructor
 */
TPTPM::~TPTPM(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param d second sp index that forms the tp column index j together with c
 * @param e first sp index that forms the tp row index i together with b
 * @param z second sp index that forms the tp row index i together with a
 * @param t first sp index that forms the tp column index j together with d
 * @param h second sp index that forms the tp column index j together with c
 * @return the number on place TPTPM(i,j) with the right phase.
 */
double TPTPM::operator()(int a,int b,int c,int d,int e,int z,int t,int h) const{

   if( (a == b) || (c == d) || (e == z) || (t == h))
      return 0;
   else{

      int I = TPM::gs2t(a,b);
      int J = TPM::gs2t(c,d);
      int K = TPM::gs2t(e,z);
      int L = TPM::gs2t(t,h);

      int phase = 1;

      if(a > b)
         phase *= -1;

      if(c > d)
         phase *= -1;

      if(e > z)
         phase *= -1;

      if(t > h)
         phase *= -1;

      int i = t2tpmm[I][J];
      int j = t2tpmm[K][L];

      return phase * (*this)(i,j);

   }

}

/**
 * access the elements of the matrix in tp mode
 * @param I first tp index that forms the tpmm row index i together with J
 * @param J second tp index that forms the tpmm row index i together with I
 * @param K first tp index that forms the tpmm column index j together with L
 * @param L second tp index that forms the tpmm column index j together with K
 * @return the number on place TPTPM(i,j)
 */
double TPTPM::operator()(int I,int J,int K,int L) const{

   int i = t2tpmm[I][J];
   int j = t2tpmm[K][L];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,const TPTPM &tpmm_p){

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = tpmm_p.tpmm2t[i][0];
      J = tpmm_p.tpmm2t[i][1];

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = tpmm_p.tpmm2t[j][0];
         L = tpmm_p.tpmm2t[j][1];

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         output << i << "\t" << j << "\t|\t" << I << "\t" << J << "\t" << K << "\t" << L << "\t|\t" << 
         
            "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << tpmm_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * @return the dimension of a TPTPM matrix
 */
int TPTPM::gn(){

   return tpmm2t.size();

}

/**
 * access to the lists from outside the class
 */
int TPTPM::gt2tpmm(int i,int j){

   return t2tpmm[i][j];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return a, == 1 return b
 */
int TPTPM::gtpmm2t(int i,int option){

   return tpmm2t[i][option];

}

 /**
 * construct the antisymmetrized "symmetric direct product" of two PHM matrix
 */
void TPTPM::dp(const PHM &phm){ 

   int I,J,K,L;

   int a,b,c,d,e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         //first 16 direct product terms of the PHM object
         (*this)(i,j) = phm(a,d,e,h) * phm(c,b,t,z) + phm(a,d,t,z) * phm(c,b,e,h) - phm(a,d,z,h) * phm(c,b,t,e) - phm(a,d,t,e) * phm(c,b,z,h)

               - phm(a,d,e,t) * phm(c,b,h,z) - phm(a,d,h,z) * phm(c,b,e,t) + phm(a,d,z,t) * phm(c,b,h,e) + phm(a,d,h,e) * phm(c,b,z,t)

               - phm(b,d,e,h) * phm(c,a,t,z) - phm(b,d,t,z) * phm(c,a,e,h) + phm(b,d,z,h) * phm(c,a,t,e) + phm(b,d,t,e) * phm(c,a,z,h)

               + phm(b,d,e,t) * phm(c,a,h,z) + phm(b,d,h,z) * phm(c,a,e,t) - phm(b,d,z,t) * phm(c,a,h,e) - phm(b,d,h,e) * phm(c,a,z,t)

               - phm(a,c,e,h) * phm(d,b,t,z) - phm(a,c,t,z) * phm(d,b,e,h) + phm(a,c,z,h) * phm(d,b,t,e) + phm(a,c,t,e) * phm(d,b,z,h)

               + phm(a,c,e,t) * phm(d,b,h,z) + phm(a,c,h,z) * phm(d,b,e,t) - phm(a,c,z,t) * phm(d,b,h,e) - phm(a,c,h,e) * phm(d,b,z,t)

               + phm(b,c,e,h) * phm(d,a,t,z) + phm(b,c,t,z) * phm(d,a,e,h) - phm(b,c,z,h) * phm(d,a,t,e) - phm(b,c,t,e) * phm(d,a,z,h)

               - phm(b,c,e,t) * phm(d,a,h,z) - phm(b,c,h,z) * phm(d,a,e,t) + phm(b,c,z,t) * phm(d,a,h,e) + phm(b,c,h,e) * phm(d,a,z,t) ;

      }

   }

   this->symmetrize();

}

/**
 * yet another double trace of a direct product of DPM's
 */
void TPTPM::dpt2(const DPM &dpm){

   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;

   double *dparray = new double [M5*M];

   dpm.convert(dparray);

   int I,J,K,L;

   int a,b,c,d,e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         (*this)(i,j) = 0.0;

         for(int k = 0;k < M;++k)
            for(int l = 0;l < M;++l){

               (*this)(i,j) += dparray[a*M5 + b*M4 + k*M3 + e*M2 + z*M + l] * dparray[c*M5 + d*M4 + k*M3 + t*M2 + h*M + l]

                  + dparray[a*M5 + b*M4 + k*M3 + t*M2 + h*M + l] * dparray[c*M5 + d*M4 + k*M3 + e*M2 + z*M + l]; 

            }

      }
   }

   this->symmetrize();

   delete [] dparray;

}

/**
 * construct the antisymmetrized, double 'tilde' of a symmetric direct product of two PPHM matrices
 * @param pphm input PPHM
 */
void TPTPM::dpw2(const PPHM &pphm){ 

   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;

   double *ppharray = new double [M5*M];

   pphm.convert(ppharray);

   int I,J,K,L;

   int a,b,c,d,e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         (*this)(i,j) = 0.0;

         for(int k = 0;k < M;++k)
            for(int l = 0;l < M;++l){

               (*this)(i,j) += ppharray[k*M5 + d*M4 + a*M3 + l*M2 + h*M + e] * ppharray[k*M5 + b*M4 + c*M3 + l*M2 + z*M + t]

                  + ppharray[k*M5 + d*M4 + a*M3 + l*M2 + z*M + t] * ppharray[k*M5 + b*M4 + c*M3 + l*M2 + h*M + e]

                  - ppharray[k*M5 + d*M4 + b*M3 + l*M2 + h*M + e] * ppharray[k*M5 + a*M4 + c*M3 + l*M2 + z*M + t]

                  - ppharray[k*M5 + d*M4 + b*M3 + l*M2 + z*M + t] * ppharray[k*M5 + a*M4 + c*M3 + l*M2 + h*M + e]

                  - ppharray[k*M5 + c*M4 + a*M3 + l*M2 + h*M + e] * ppharray[k*M5 + b*M4 + d*M3 + l*M2 + z*M + t]

                  - ppharray[k*M5 + c*M4 + a*M3 + l*M2 + z*M + t] * ppharray[k*M5 + b*M4 + d*M3 + l*M2 + h*M + e]

                  + ppharray[k*M5 + c*M4 + b*M3 + l*M2 + h*M + e] * ppharray[k*M5 + a*M4 + d*M3 + l*M2 + z*M + t]

                  + ppharray[k*M5 + c*M4 + b*M3 + l*M2 + z*M + t] * ppharray[k*M5 + a*M4 + d*M3 + l*M2 + h*M + e]

                  - ppharray[k*M5 + d*M4 + a*M3 + l*M2 + h*M + z] * ppharray[k*M5 + b*M4 + c*M3 + l*M2 + e*M + t]

                  - ppharray[k*M5 + d*M4 + a*M3 + l*M2 + e*M + t] * ppharray[k*M5 + b*M4 + c*M3 + l*M2 + h*M + z]

                  + ppharray[k*M5 + d*M4 + b*M3 + l*M2 + h*M + z] * ppharray[k*M5 + a*M4 + c*M3 + l*M2 + e*M + t]

                  + ppharray[k*M5 + d*M4 + b*M3 + l*M2 + e*M + t] * ppharray[k*M5 + a*M4 + c*M3 + l*M2 + h*M + z]

                  + ppharray[k*M5 + c*M4 + a*M3 + l*M2 + h*M + z] * ppharray[k*M5 + b*M4 + d*M3 + l*M2 + e*M + t]

                  + ppharray[k*M5 + c*M4 + a*M3 + l*M2 + e*M + t] * ppharray[k*M5 + b*M4 + d*M3 + l*M2 + h*M + z]

                  - ppharray[k*M5 + c*M4 + b*M3 + l*M2 + h*M + z] * ppharray[k*M5 + a*M4 + d*M3 + l*M2 + e*M + t]

                  - ppharray[k*M5 + c*M4 + b*M3 + l*M2 + e*M + t] * ppharray[k*M5 + a*M4 + d*M3 + l*M2 + h*M + z]

                  - ppharray[k*M5 + d*M4 + a*M3 + l*M2 + t*M + e] * ppharray[k*M5 + b*M4 + c*M3 + l*M2 + z*M + h]

                  - ppharray[k*M5 + d*M4 + a*M3 + l*M2 + z*M + h] * ppharray[k*M5 + b*M4 + c*M3 + l*M2 + t*M + e]

                  + ppharray[k*M5 + d*M4 + b*M3 + l*M2 + t*M + e] * ppharray[k*M5 + a*M4 + c*M3 + l*M2 + z*M + h]

                  + ppharray[k*M5 + d*M4 + b*M3 + l*M2 + z*M + h] * ppharray[k*M5 + a*M4 + c*M3 + l*M2 + t*M + e]

                  + ppharray[k*M5 + c*M4 + a*M3 + l*M2 + t*M + e] * ppharray[k*M5 + b*M4 + d*M3 + l*M2 + z*M + h]

                  + ppharray[k*M5 + c*M4 + a*M3 + l*M2 + z*M + h] * ppharray[k*M5 + b*M4 + d*M3 + l*M2 + t*M + e]

                  - ppharray[k*M5 + c*M4 + b*M3 + l*M2 + t*M + e] * ppharray[k*M5 + a*M4 + d*M3 + l*M2 + z*M + h]

                  - ppharray[k*M5 + c*M4 + b*M3 + l*M2 + z*M + h] * ppharray[k*M5 + a*M4 + d*M3 + l*M2 + t*M + e]

                  + ppharray[k*M5 + d*M4 + a*M3 + l*M2 + t*M + z] * ppharray[k*M5 + b*M4 + c*M3 + l*M2 + e*M + h]

                  + ppharray[k*M5 + d*M4 + a*M3 + l*M2 + e*M + h] * ppharray[k*M5 + b*M4 + c*M3 + l*M2 + t*M + z]

                  - ppharray[k*M5 + d*M4 + b*M3 + l*M2 + t*M + z] * ppharray[k*M5 + a*M4 + c*M3 + l*M2 + e*M + h]

                  - ppharray[k*M5 + d*M4 + b*M3 + l*M2 + e*M + h] * ppharray[k*M5 + a*M4 + c*M3 + l*M2 + t*M + z]

                  - ppharray[k*M5 + c*M4 + a*M3 + l*M2 + t*M + z] * ppharray[k*M5 + b*M4 + d*M3 + l*M2 + e*M + h]

                  - ppharray[k*M5 + c*M4 + a*M3 + l*M2 + e*M + h] * ppharray[k*M5 + b*M4 + d*M3 + l*M2 + t*M + z]

                  + ppharray[k*M5 + c*M4 + b*M3 + l*M2 + t*M + z] * ppharray[k*M5 + a*M4 + d*M3 + l*M2 + e*M + h]

                  + ppharray[k*M5 + c*M4 + b*M3 + l*M2 + e*M + h] * ppharray[k*M5 + a*M4 + d*M3 + l*M2 + t*M + z];

            }

      }
   }

   delete [] ppharray;

}

/**
 * construct the "symmetric direct product" of two PPHM matrices, once traced and once 'tilded'.
 * BE CAREFULL! NOT SYMMETRIC!
 */
void TPTPM::dptw(const PPHM &pphm){ 

   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;

   double *ppharray = new double [M5*M];

   pphm.convert(ppharray);

   int I,J,K,L;

   int a,b,c,d,e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = 0;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         (*this)(i,j) =  0.0;

         for(int k = 0;k < M;++k)
            for(int l = 0;l < M;++l){

               (*this)(i,j) += ppharray[a*M5 + b*M4 + k*M3 + l*M2 + h*M + e] * ppharray[c*M5 + d*M4 + k*M3 + l*M2 + z*M + t]   

                  + ppharray[a*M5 + b*M4 + k*M3 + l*M2 + z*M + t] * ppharray[c*M5 + d*M4 + k*M3 + l*M2 + h*M + e]

                  - ppharray[a*M5 + b*M4 + k*M3 + l*M2 + h*M + z] * ppharray[c*M5 + d*M4 + k*M3 + l*M2 + e*M + t]   

                  - ppharray[a*M5 + b*M4 + k*M3 + l*M2 + e*M + t] * ppharray[c*M5 + d*M4 + k*M3 + l*M2 + h*M + z]

                  - ppharray[a*M5 + b*M4 + k*M3 + l*M2 + t*M + e] * ppharray[c*M5 + d*M4 + k*M3 + l*M2 + z*M + h]   

                  - ppharray[a*M5 + b*M4 + k*M3 + l*M2 + z*M + h] * ppharray[c*M5 + d*M4 + k*M3 + l*M2 + t*M + e]

                  + ppharray[a*M5 + b*M4 + k*M3 + l*M2 + t*M + z] * ppharray[c*M5 + d*M4 + k*M3 + l*M2 + e*M + h]   

                  + ppharray[a*M5 + b*M4 + k*M3 + l*M2 + e*M + h] * ppharray[c*M5 + d*M4 + k*M3 + l*M2 + t*M + z];

            }

      }

   }

   delete [] ppharray;

}

/**
 * construct the "symmetric direct product" of two PPHM matrices, doubly traced
 */
void TPTPM::dpt2(const PPHM &pphm){ 

   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;

   double *ppharray = new double [M5*M];

   pphm.convert(ppharray);

   int I,J,K,L;

   int a,b,c,d,e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      a = TPM::gt2s(I,0);
      b = TPM::gt2s(I,1);

      c = TPM::gt2s(J,0);
      d = TPM::gt2s(J,1);

      for(int j = i;j < TPTPM::gn();++j){

         K = TPTPM::gtpmm2t(j,0);
         L = TPTPM::gtpmm2t(j,1);

         e = TPM::gt2s(K,0);
         z = TPM::gt2s(K,1);

         t = TPM::gt2s(L,0);
         h = TPM::gt2s(L,1);

         (*this)(i,j) =  0.0;

         for(int k = 0;k < M;++k)
            for(int l = 0;l < M;++l){

               (*this)(i,j) += ppharray[a*M5 + b*M4 + k*M3 + e*M2 + z*M + l] * ppharray[c*M5 + d*M4 + k*M3 + t*M2 + h*M + l]

                  + ppharray[a*M5 + b*M4 + k*M3 + t*M2 + h*M + l] * ppharray[c*M5 + d*M4 + k*M3 + e*M2 + z*M + l];

            }

      }

   }

   delete [] ppharray;

   this->symmetrize();

}

void TPTPM::I(const TPM &tpm ){

   Basis basis(M,N);

   TPM **lmap = new TPM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      lmap[i] = new TPM(M,N);

      lmap[i]->L_map(tpm,basis[i]);

   }

   for(int i = 0;i < TPTPM::gn();++i){

      for(int j = i;j < TPTPM::gn();++j){

         (*this)(i,j) =  lmap[i]->ddot(basis[j]);

      }
   }

   for(int i = 0;i < TPTPM::gn();++i)
      delete lmap[i];

   delete [] lmap;

   this->symmetrize();

}

void TPTPM::Q(const TPM &tpm ){

   Basis basis(M,N);

   TPM **qbasis = new TPM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      qbasis[i] = new TPM(M,N);

      qbasis[i]->Q(1,basis[i]);

   }

   TPM **lmap = new TPM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      lmap[i] = new TPM(M,N);

      lmap[i]->L_map(tpm,*qbasis[i]);

   }

   for(int i = 0;i < TPTPM::gn();++i){

      for(int j = i;j < TPTPM::gn();++j){

         (*this)(i,j) =  lmap[i]->ddot(*qbasis[j]);

      }
   }

   for(int i = 0;i < TPTPM::gn();++i){

      delete lmap[i];
      delete qbasis[i];

   }

   delete [] lmap;
   delete [] qbasis;

   this->symmetrize();

}


void TPTPM::G(const PHM &phm){

   Basis basis(M,N);

   PHM **gbasis = new PHM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      gbasis[i] = new PHM(M,N);

      gbasis[i]->G(1,basis[i]);

   }

   PHM **lmap = new PHM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      lmap[i] = new PHM(M,N);

      lmap[i]->L_map(phm,*gbasis[i]);

   }

   for(int i = 0;i < TPTPM::gn();++i){

      for(int j = i;j < TPTPM::gn();++j){

         (*this)(i,j) =  lmap[i]->ddot(*gbasis[j]);

      }
   }

   for(int i = 0;i < TPTPM::gn();++i){

      delete lmap[i];
      delete gbasis[i];

   }

   delete [] lmap;
   delete [] gbasis;

   this->symmetrize();

}


void TPTPM::T(const DPM &dpm){

   Basis basis(M,N);

   DPM **tbasis = new DPM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      tbasis[i] = new DPM(M,N);

      tbasis[i]->T(1,basis[i]);

   }

   DPM **lmap = new DPM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      lmap[i] = new DPM(M,N);

      lmap[i]->L_map(dpm,*tbasis[i]);

   }

   for(int i = 0;i < TPTPM::gn();++i){

      for(int j = i;j < TPTPM::gn();++j){

         (*this)(i,j) =  lmap[i]->ddot(*tbasis[j]);

      }
   }

   for(int i = 0;i < TPTPM::gn();++i){

      delete lmap[i];
      delete tbasis[i];

   }

   delete [] lmap;
   delete [] tbasis;

   this->symmetrize();

}

void TPTPM::T(const PPHM &pphm){

   Basis basis(M,N);

   PPHM **tbasis1 = new PPHM * [TPTPM::gn()];
   PPHM **tbasis2 = new PPHM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      tbasis1[i] = new PPHM(M,N);

      tbasis1[i]->T(0,basis[i]);

   }

   for(int i = 0;i < TPTPM::gn();++i){

      tbasis2[i] = new PPHM(M,N);

      tbasis2[i]->T(2,basis[i]);

   }

   PPHM **lmap = new PPHM * [TPTPM::gn()];

   for(int i = 0;i < TPTPM::gn();++i){

      lmap[i] = new PPHM(M,N);

      lmap[i]->L_map(pphm,*tbasis1[i]);

   }

   for(int i = 0;i < TPTPM::gn();++i){

      for(int j = 0;j < TPTPM::gn();++j){

         (*this)(i,j) =  lmap[i]->ddot(*tbasis2[j]);

      }
   }

   for(int i = 0;i < TPTPM::gn();++i){

      delete lmap[i];
      delete tbasis2[i];
      delete tbasis1[i];

   }

   delete [] lmap;
   delete [] tbasis1;
   delete [] tbasis2;

   //this->symmetrize();

}
