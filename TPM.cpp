#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;

//include headers files and define some defines for the conditions
#include "headers/include.h"

int TPM::counter = 0;

int **TPM::t2s;
int **TPM::s2t;

/**
 * standard constructor: constructs Matrix object of dimension M*(M - 1)/2 and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and tp basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
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
            s2t[j][i] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

   }

   ++counter;

}

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpm_c
 * if counter == 0, the lists containing the relationship between sp and tp basis.
 * @param tpm_c object that will be copied into this.
 */
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
            s2t[j][i] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

   }

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists t2s en s2t will be deleted.
 * 
 */
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

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * TPM(a,b,c,d) = -TPM(b,a,c,d) = -TPM(a,b,d,c) = TPM(b,a,c,d)
 * @param a first sp index that forms the tp row index i together with b
 * @param b second sp index that forms the tp row index i together with a
 * @param c first sp index that forms the tp column index j together with d
 * @param d second sp index that forms the tp column index j together with c
 * @return the number on place TPM(i,j) with the right phase.
 */
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

ostream &operator<<(ostream &output,const TPM &tpm_p){

   for(int i = 0;i < tpm_p.n;++i)
      for(int j = 0;j < tpm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << tpm_p.t2s[i][0] << "\t" << tpm_p.t2s[i][1]

            << "\t" << tpm_p.t2s[j][0] << "\t" << tpm_p.t2s[j][1] << "\t" << tpm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return the number of particles
 */
int TPM::gN(){

   return N;

}

/**
 * @return the dimension of sp space
 */
int TPM::gM(){

   return M;

}

/**
 * @return the dimension of tp space and of the matrix
 */
int TPM::gn(){

   return n;

}

/**
 * construct the hubbard hamiltonian with on site repulsion U
 * @param U onsite repulsion term
 */
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

         (*this)(j,i) = (*this)(i,j);
      }

   }

}

/**
 * The Q-map, maps a TPM (tpm_d) on a TPM (*this), definition in notes:
 * @param tpm_d The input TPM
 */
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

         (*this)(j,i) = (*this)(i,j);
      }
   }

}

#ifdef __G_CON

/**
 * The G down map, maps a PHM (phm) on a TPM (*this)
 * @param phm input PHM
 */
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

         (*this)(j,i) = (*this)(i,j);
      }

   }

}

#endif

/**
 * Initializes the TPM on the a unit matrix with trace N*(N - 1)/2, i.e. the number of tp pairs
 */
void TPM::init(){

   double ward = N*(N - 1.0)/(2.0*n);

   for(int i = 0;i < n;++i){

      (*this)(i,i) = ward;

      for(int j = i + 1;j < n;++j)
         (*this)(j,i) = (*this)(i,j) = 0.0;

   }

}

/**
 * Project the TPM on a traceless matrix:\n\n
 * this = this - Tr(this)/n 1
 */
void TPM::proj_Tr(){

   double ward = (this->trace())/n;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= ward;

}

/**
 * Construct the right hand side of the Newton equation for the determination of the search direction, 
 * the gradient of the potential:
 * @param t scaling factor of the potential
 * @param ham Hamiltonian of the current problem
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 */
void TPM::constr_grad(double t,TPM &ham,SUP &P){

   //eerst P conditie 
   *this = P.tpm(0);

   //de Q conditie toevoegen

#ifdef __Q_CON

   TPM hulp(M,N);

   hulp.Q(P.tpm(1));

   *this += hulp;

#endif

   //de G conditie indien nodig

#ifdef __G_CON

   hulp.G(P.phm());

   *this += hulp;

#endif

   //de T_1 conditie toevoegen indien nodig

#ifdef __T1_CON

   hulp.T(P.dpm());

   *this += hulp;

#endif

#ifdef __T2_CON

   hulp.T(P.pphm());

   *this +=hulp;
#endif

#ifdef __T2P_CON

   hulp.T(P.t2pm());

   *this +=hulp;
#endif

   this->dscal(t);

   *this -= ham;

   this->proj_Tr();

}

/**
 * solve the Newton equations for the determination of the search direction,
 * @param t scaling factor of the potential
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 * @param b right hand side (the gradient constructed int TPM::constr_grad)
 * @return nr of iterations needed to converge to the desired accuracy
 */
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

/**
 * perform a line search what step size in along the Newton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
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

/**
 * The hessian-map of the Newton system:
 * @param t potential scaling factor
 * @param b the TPM on which the hamiltonian will work, the image will be put in (*this)
 * @param P the SUP matrix containing the constraints, (can be seen as the metric).
 */
void TPM::H(double t,TPM &b,SUP &P){

   //eerst de P conditie:

   this->L_map(P.tpm(0),b);

   //de Q conditie toevoegen:
#ifdef __Q_CON

   //hulpje
   TPM hulp(M,N);

   //maak Q(b)
   TPM Q_b(M,N);
   Q_b.Q(b);

   //stop Q(rdm)^{-1}Q(b)Q(rdm)^{-1} in hulp
   hulp.L_map(P.tpm(1),Q_b);

   //maak Q(hulp) en stop in Q_b
   Q_b.Q(hulp);

   //en tel op bij this
   *this += Q_b;

#endif

   //de G conditie toevoegen:
   
#ifdef __G_CON

   //hulpje voor het PHM stuk
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
   
   //de T_1 conditie toevoegen

#ifdef __T1_CON

   //hulpjes voor het DPM stuk
   DPM hulp_dp(M,N);
   DPM T1_b(M,N);

   //stop T1(b) in T1_b
   T1_b.T(b);

   hulp_dp.L_map(P.dpm(),T1_b);

   hulp.T(hulp_dp);

   *this += hulp;

#endif

#ifdef __T2_CON

   PPHM hulp_pph(M,N);
   PPHM T2_b(M,N);

   T2_b.T(b);

   hulp_pph.L_map(P.pphm(),T2_b);

   hulp.T(hulp_pph);

   *this+=hulp;

#endif

#ifdef __T2P_CON

   T2PM hulp_t2p(M,N);
   T2PM T2P_b(M,N);

   T2P_b.T(b);

   hulp_t2p.L_map(P.t2pm(),T2P_b);

   hulp.T(hulp_t2p);

   *this+=hulp;

#endif

   //nog schalen met t:
   this->dscal(t);

   //en projecteren op spoorloze ruimte
   this->proj_Tr();

}

/**
 * perform a line search what step size in along the Newton direction is ideal, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param rdm TPM containing the current approximation of the rdm.
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,TPM &rdm,TPM &ham){

   SUP P(M,N);

   P.fill(rdm);

   P.invert();

   return this->line_search(t,P,ham);

}

/**
 * calculate the trace of one pair of sp indices of a DPM an put in (*this):\n\n
 * TPM(a,b,d,e) = \\sum_{c} DPM(a,b,c,d,e,c)
 * @param dpm input DPM
 */
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

	 (*this)(j,i) = (*this)(i,j);
      }
   }
}

/**
 * map a DPM (dpm) on a TPM (*this) with T1_down map:
 * @param dpm The input DPM
 */
void TPM::T(DPM &dpm){

   double ward = 2.0*dpm.trace()/(N*(N - 1.0));

   int a,b,c,d;

   TPM tpm(M,N);
   tpm.bar(dpm);

   SPM spm(tpm);

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

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

         (*this)(j,i) = (*this)(i,j);
      }

   }
}

void TPM::bar(PPHM &pphm)
{
    int a,b,c,d;

    for(int i=0;i<n;i++)
    {
	a = t2s[i][0];
	b = t2s[i][1];

	for(int j=i;j<n;j++)
	{
	    c = t2s[j][0];
	    d = t2s[j][1];

	    (*this)(i,j) = 0.0;

	    for(int l=0;l<M;l++)
		(*this)(i,j) += pphm(a,b,l,c,d,l);

	    (*this)(j,i) = (*this)(i,j);
	}
    }
}

void TPM::T(PPHM &pphm)
{
    double brecht= 1.0/(2*(N-1));
    int a,b,c,d;

    // A bar
    TPM tpm(M,N);
    tpm.bar(pphm);

    // A dubble bar
    SPM spm(M,N);
    spm.bar(pphm);

    // A tilde bar
    PHM phm(M,N);
    phm.bar(pphm);

    for(int i=0;i<n;i++)
    {
	a = t2s[i][0];
	b = t2s[i][1];

	for(int j=i;j<n;j++)
	{
	    c = t2s[j][0];
	    d = t2s[j][1];

	    (*this)(i,j) = 0;

	    if(b==d)
		(*this)(i,j) += spm(a,c);

	    if(a==d)
		(*this)(i,j) -= spm(b,c);

	    if(b==c)
		(*this)(i,j) -= spm(a,d);

	    if(a==c)
		(*this)(i,j) += spm(b,d);

	    (*this)(i,j) *= brecht;

	    (*this)(i,j) += tpm(i,j);

	    (*this)(i,j) -= phm(d,a,b,c)-phm(d,b,a,c)-phm(c,a,b,d)+phm(c,b,a,d);

	    (*this)(j,i) = (*this)(i,j);
	}
    }
}

void TPM::bar(T2PM &t2pm)
{
   int a,b,c,d;

   for(int i=0;i<n;i++)
   {
      a = t2s[i][0];
      b = t2s[i][1];

      for(int j=i;j<n;j++)
      {
         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l=0;l<M;l++)
            (*this)(i,j) += t2pm(a,b,l,c,d,l);

         (*this)(j,i) = (*this)(i,j);
      }
   }
}

void TPM::T(T2PM &t2pm)
{
   int a,b,c,d;
   int n_pph = M*M*(M-1)/2;

   double brecht = 1.0/(2*(N-1));
   double matthias = 2*brecht;

   // A bar
   TPM tpm(M,N);
   tpm.bar(t2pm);

   // A dubble bar
   SPM spm(M,N);
   spm.bar(t2pm);

   // A tilde bar
   PHM phm(M,N);
   phm.bar(t2pm);

   for(int i=0;i<n;i++)
   {
      a = t2s[i][0];
      b = t2s[i][1];

      for(int j=i;j<n;j++)
      {
         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0;

         // T2 part
         if(b==d)
            (*this)(i,j) += spm(a,c);

         if(a==d)
            (*this)(i,j) -= spm(b,c);

         if(b==c)
            (*this)(i,j) -= spm(a,d);

         if(a==c)
            (*this)(i,j) += spm(b,d);

         (*this)(i,j) *= brecht;

         (*this)(i,j) += tpm(i,j);

         (*this)(i,j) -= phm(d,a,b,c)-phm(d,b,a,c)-phm(c,a,b,d)+phm(c,b,a,d);


         // rho part
         if( b == d )
            (*this)(i,j) += matthias*t2pm(n_pph+c,n_pph+a);

         if( a == d )
            (*this)(i,j) -= matthias*t2pm(n_pph+c,n_pph+b);

         if( b == c )
            (*this)(i,j) -= matthias*t2pm(n_pph+d,n_pph+a);

         if( a == c )
            (*this)(i,j) += matthias*t2pm(n_pph+d,n_pph+b);


         // up diagonal part
         (*this)(i,j) += t2pm(a,b,d,c)+t2pm(d,c,a,b);
         (*this)(i,j) -= t2pm(a,b,c,d)+t2pm(d,c,b,a);

         (*this)(j,i) = (*this)(i,j);
      }
   }
}

/* vim: set ts=3 sw=3 expandtab :*/
