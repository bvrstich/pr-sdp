#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;

#include "headers/TPM.h"
#include "headers/SPM.h"
#include "headers/SUP.h"
#include "headers/EIG.h"
#include "headers/lapack.h"

#ifdef PQG

#include "headers/PHM.h"

#endif

#ifdef PQGT1

#include "headers/DPM.h"

#endif

int TPM::counter = 0;

int **TPM::t2s;
int **TPM::s2t;

//constructor:
TPM::TPM(int M,int N) : Matrix(M*(M - 1)/2){

   this->N = N;
   this->M = M;
   this->n = M*(M - 1)/2;

   if(counter == 0){

      //allocatie van sp2tp
      s2t = new int * [M];
      s2t[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2t[i] = s2t[i - 1] + M;

      //allocatie van t2s
      t2s = new int * [n];

      for(int i = 0;i < n;++i)
         t2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j){

            s2t[i][j] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j)
            s2t[j][i] = s2t[i][j];

   }

   ++counter;

}

//copy constructor
TPM::TPM(TPM &tpm_c) : Matrix(tpm_c){

   this->N = tpm_c.N;
   this->M = tpm_c.M;
   this->n = M*(M - 1)/2;

   if(counter == 0){

      //allocatie van sp2tp
      s2t = new int * [M];
      s2t[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2t[i] = s2t[i - 1] + M;

      //allocatie van tp2sp
      t2s = new int * [n];

      for(int i = 0;i < n;++i)
         t2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j){

            s2t[i][j] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j)
            s2t[j][i] = s2t[i][j];

   }

   ++counter;

}

//destructor
TPM::~TPM(){

   if(counter == 1){

      delete [] s2t[0];
      delete [] s2t;

      for(int i = 0;i < n;++i)
         delete [] t2s[i];

      delete [] t2s;

   }

   --counter;

}

//access the numbers: sp indices
double TPM::operator()(int a,int b,int c,int d) const{

   if( (a == b) || (c == d) )
      return 0;
   else{

      int i = s2t[a][b];
      int j = s2t[c][d];

      int phase = 1;

      if(a > b)
         phase *= -1;
      if(c > d)
         phase *= -1;

      return phase*(*this)(i,j);

   }

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,const TPM &tpm_p){

   for(int i = 0;i < tpm_p.n;++i)
      for(int j = 0;j < tpm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << tpm_p.t2s[i][0] << "\t" << tpm_p.t2s[i][1]

            << "\t" << tpm_p.t2s[j][0] << "\t" << tpm_p.t2s[j][1] << "\t" << tpm_p(i,j) << endl;

      }

   return output;

}

int TPM::gN(){

   return N;

}

int TPM::gM(){

   return M;

}

int TPM::gn(){

   return n;

}

void TPM::hubbard(double U){

   int a,b,c,d;//sp orbitals

   double ward = 1.0/(N - 1.0);

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0;

         //eerst hopping
         if( (a == c) && ( ( (b + 2)%M == d ) || ( b == (d + 2)%M ) ) )
            (*this)(i,j) -= ward;

         if( (b == c) && ( ( (a + 2)%M == d ) || ( a == (d + 2)%M ) ) )
            (*this)(i,j) += ward;

         if( (b == d) && ( ( (a + 2)%M == c ) || ( a == (c + 2)%M ) ) )
            (*this)(i,j) -= ward;

         //on site interaction
         if( (a % 2) == 0 && (c % 2) == 0 )
            if(a == (b - 1) && c == (d - 1) && a == c)
               (*this)(i,j) += U;

      }

   }
   
   this->symmetrize();

}

void TPM::Q(TPM &tpm_d){

   SPM spm(tpm_d);

   double ward = tpm_d.trace()/(N*(N - 1.0)/2.0);

   for(int i = 0;i < n;++i){

      int a = t2s[i][0];
      int b = t2s[i][1];

      for(int j = i;j < n;++j){

         int c = t2s[j][0];
         int d = t2s[j][1];

         (*this)(i,j) = tpm_d(i,j);

         if(i == j)
            (*this)(i,i) += ward;

         if(a == c)
            (*this)(i,j) -= spm(b,d);

         if(b == c)
            (*this)(i,j) += spm(a,d);

         if(b == d)
            (*this)(i,j) -= spm(a,c);

      }
   }

   this->symmetrize();

}

#ifndef PQ

void TPM::G(PHM &phm){

   SPM spm(phm);

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = phm(b,d,c,a) - phm(a,d,c,b) - phm(b,c,d,a) + phm(a,c,d,b);

         if(b == d)
            (*this)(i,j) += spm(a,c);

         if(b == c)
            (*this)(i,j) -= spm(a,d);

         if(a == c)
            (*this)(i,j) += spm(b,d);

      }

   }

   this->symmetrize();

}

#endif

void TPM::init(){

   double ward = N*(N - 1.0)/(2.0*n);

   for(int i = 0;i < n;++i){

      (*this)(i,i) = ward;

      for(int j = i + 1;j < n;++j)
         (*this)(j,i) = (*this)(i,j) = 0.0;

   }

}

void TPM::proj_Tr(){

   double ward = (this->trace())/(double)n;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= ward;

}

//maak de gradient van de barrierefunctie
void TPM::constr_grad(double t,TPM &ham,SUP &P){

   //eerst P conditie 
   *this = P.tpm(0);

   //dan de Q conditie
   TPM hulp(M,N);

   hulp.Q(P.tpm(1));

   *this += hulp;

#ifndef PQ

   //tenslotte de G conditie
   hulp.G(P.phm());

   *this += hulp;

#endif

   this->dscal(t);

   *this -= ham;

   this->proj_Tr();

}

//los het Newton stelsel H delta = -G op:
int TPM::solve(double t,SUP &P,TPM &b){

   int iter = 0;

   //delta = 0
   *this = 0;

   //residu:
   TPM r(b);

   //norm van het residu
   double rr = r.ddot(r);

   //enkele variabelen
   double rr_old,ward;

   TPM Hb(M,N);

   while(rr > 1.0e-7){ 

      Hb.H(t,b,P);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe variabelen berekenen en oude overdragen
      rr_old = rr;
      rr = r.ddot(r);

      //nieuwe b nog:
      b.dscal(rr/rr_old);

      b += r;

      ++iter;

   }

   return iter;

}

double TPM::line_search(double t,SUP &P,TPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta(M,N);

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp(M,N);

   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

void TPM::H(double t,TPM &b,SUP &P){

   //hulpje
   TPM hulp(M,N);

   //maak Q(b)
   TPM Q_b(M,N);
   Q_b.Q(b);

   //stop Q(rdm)^{-1}Q(b)Q(rdm)^{-1} in hulp
   hulp.L_map(P.tpm(1),Q_b);

   //maak Q(hulp) en stop in this 
   this->Q(hulp);

   //stop rdm^{-1} b rdm^{-1} in hulp
   hulp.L_map(P.tpm(0),b);

   //en tel op bij this
   *this += hulp;

#ifndef PQ

   //nu de G conditie:
   PHM hulp_ph(M,N);
   PHM G_b(M,N);

   //stop G(b) in G_b
   G_b.G(b);

   //bereken G(rdm)^{-1}G(b)G(rdm)^{-1} en stop in hulp_ph
   hulp_ph.L_map(P.phm(),G_b);

   //tenslotte nog de antisymmetrische G hierop:
   hulp.G(hulp_ph);

   //en optellen bij this
   *this += hulp;

#endif

   //nog schalen met t:
   this->dscal(t);

   //en projecteren op spoorloze ruimte
   this->proj_Tr();

}

double TPM::line_search(double t,TPM &rdm,TPM &ham){

   SUP P(M,N);

   P.fill(rdm);

   P.invert();

   return this->line_search(t,P,ham);

}

#ifdef T_1

void TPM::bar(DPM &dpm){

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < M;++l)
            (*this)(i,j) += dpm(a,b,l,c,d,l);

      }
   }

   this->symmetrize();

}

//de T_1 down eigenlijk een Q-like afbeelding
void TPM::T(DPM &dpm){

   double ward = 2.0*dpm.trace()/(N*(N - 1.0));

   int a,b,c,d;

   TPM tpm(M,N);
   tpm.bar(dpm);

   SPM spm(tpm);

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = tpm(i,j);

         if(i == j)
            (*this)(i,j) += ward;

         if(b == d)
            (*this)(i,j) -= 0.5*spm(a,c);

         if(b == c)
            (*this)(i,j) += 0.5*spm(a,d);

         if(a == c)
            (*this)(i,j) -= 0.5*spm(b,d);

      }

   }

   this->symmetrize();

}

#endif
