#include "headers/SUP.h"

using std::endl;

ostream &operator<<(ostream &output,SUP &SZ)
{
    output << (*SZ.SZ_tp[0]) << endl;
    output << (*SZ.SZ_tp[1]) << endl;

#ifdef __G_CON
    output << (*SZ.SZ_ph) << endl;
#endif

#ifdef __T1_CON
    output << (*SZ.SZ_dp) << endl;
#endif

#ifdef __T2_CON
    output << (*SZ.SZ_pph) << endl;
#endif

    return output;
}

//constructor
SUP::SUP(int M,int N)
{
    this->M=M;
    this->N=N;
    this->n_tp=M*(M-1)/2;

    dim=2*n_tp;

    SZ_tp = new TPM *[2];
    SZ_tp[0]=new TPM(M,N);
    SZ_tp[1]=new TPM(M,N);

#ifdef __G_CON
      n_ph = M*M;
      dim+=n_ph;
      //! double pointer to PHM
      SZ_ph = new PHM(M,N);
#endif

#ifdef __T1_CON
      n_dp=M*(M-1)*(M-2)/6;
      dim+=n_dp;
      SZ_dp=new DPM(M,N);
#endif

#ifdef __T2_CON
      n_pph=M*M*(M-1)/2;
      dim+=n_pph;
      SZ_pph=new PPHM(M,N);
#endif
}


//copy constructor
SUP::SUP(SUP &SZ)
{
    M=SZ.M;
    N=SZ.N;
    n_tp=SZ.n_tp;
    dim=SZ.dim;

    SZ_tp = new TPM *[2];
    SZ_tp[0]=new TPM(M,N);
    SZ_tp[1]=new TPM(M,N);

    (*SZ_tp[0])=(*SZ.SZ_tp[0]);
    (*SZ_tp[1])=(*SZ.SZ_tp[1]);

#ifdef __G_CON
      n_ph = SZ.n_ph;
      SZ_ph = new PHM(M,N);
      (*SZ_ph) = (*SZ.SZ_ph);
#endif

#ifdef __T1_CON
      n_dp=SZ.n_dp;
      SZ_dp=new DPM(M,N);
      (*SZ_dp)=(*SZ.SZ_dp);
#endif

#ifdef __T2_CON
      n_pph=SZ.n_pph;
      SZ_pph=new PPHM(M,N);
      (*SZ_pph)=(*SZ.SZ_pph);
#endif
}

//destructor
SUP::~SUP()
{
    delete SZ_tp[0];
    delete SZ_tp[1];
    delete [] SZ_tp;

#ifdef __G_CON
      delete SZ_ph;
#endif

#ifdef __T1_CON
      delete SZ_dp;
#endif

#ifdef __T2_CON
      delete SZ_pph;
#endif
}

//overload += operator
SUP &SUP::operator+=(SUP &SZ)
{
    (*SZ_tp[0])+=(*SZ.SZ_tp[0]);
    (*SZ_tp[1])+=(*SZ.SZ_tp[1]);

#ifdef __G_CON
      (*SZ_ph)+=(*SZ.SZ_ph);
#endif

#ifdef __T1_CON
      (*SZ_dp)+=(*SZ.SZ_dp);
#endif

#ifdef __T2_CON
      (*SZ_pph)+=(*SZ.SZ_pph);
#endif

    return *this;
}

//overload -= operator
SUP &SUP::operator-=(SUP &SZ)
{
    (*SZ_tp[0])-=(*SZ.SZ_tp[0]);
    (*SZ_tp[1])-=(*SZ.SZ_tp[1]);

#ifdef __G_CON
      (*SZ_ph)-=(*SZ.SZ_ph);
#endif

#ifdef __T1_CON
      (*SZ_dp)-=(*SZ.SZ_dp);
#endif

#ifdef __T2_CON
      (*SZ_pph)-=(*SZ.SZ_pph);
#endif

    return *this;
}

//overload equality operator
SUP &SUP::operator=(SUP &SZ)
{
    (*SZ_tp[0])=(*SZ.SZ_tp[0]);
    (*SZ_tp[1])=(*SZ.SZ_tp[1]);

#ifdef __G_CON
    (*SZ_ph)=(*SZ.SZ_ph);
#endif

#ifdef __T1_CON
    (*SZ_dp)=(*SZ.SZ_dp);
#endif

#ifdef __T2_CON
    (*SZ_pph)=(*SZ.SZ_pph);
#endif

    return *this;
}

//overload equality operator
SUP &SUP::operator=(double &a)
{
    (*SZ_tp[0])=a;
    (*SZ_tp[1])=a;

#ifdef __G_CON
    (*SZ_ph)=a;
#endif

#ifdef __T1_CON
    (*SZ_dp)=a;
#endif

#ifdef __T2_CON
    (*SZ_pph)=a;
#endif

    return *this;
}

int SUP::gN()
{
    return N;
}

int SUP::gM()
{
    return M;
}

double SUP::ddot(SUP &SZ)
{
   double brecht;

   brecht = SZ_tp[0]->ddot(*SZ.SZ_tp[0]);
   brecht += SZ_tp[1]->ddot(*SZ.SZ_tp[1]);

#ifdef __G_CON
   brecht += SZ_ph->ddot(*SZ.SZ_ph);
#endif

#ifdef __T1_CON
   brecht += SZ_dp->ddot(*SZ.SZ_dp);
#endif

#ifdef __T2_CON
   brecht += SZ_pph->ddot(*SZ.SZ_pph);
#endif

   return brecht;
}

void SUP::invert()
{
    SZ_tp[0]->invert();
    SZ_tp[1]->invert();

#ifdef __G_CON
    SZ_ph->invert();
#endif

#ifdef __T1_CON
    SZ_dp->invert();
#endif

#ifdef __T2_CON
    SZ_pph->invert();
#endif
}

void SUP::dscal(double alpha)
{
    SZ_tp[0]->dscal(alpha);
    SZ_tp[1]->dscal(alpha);

#ifdef __G_CON
    SZ_ph->dscal(alpha);
#endif

#ifdef __T1_CON
    SZ_dp->dscal(alpha);
#endif

#ifdef __T2_CON
    SZ_pph->dscal(alpha);
#endif
}

//positieve of negatieve vierkantswortel uit een supermatrix
void SUP::sqrt(int option)
{
    SZ_tp[0]->sqrt(option);
    SZ_tp[1]->sqrt(option);

#ifdef __G_CON
    SZ_ph->sqrt(option);
#endif

#ifdef __T1_CON
    SZ_dp->sqrt(option);
#endif

#ifdef __T2_CON
    SZ_pph->sqrt(option);
#endif
}

void SUP::L_map(SUP &map,SUP &obj)
{
    SZ_tp[0]->L_map(*map.SZ_tp[0],*obj.SZ_tp[0]);
    SZ_tp[1]->L_map(*map.SZ_tp[1],*obj.SZ_tp[1]);

#ifdef __G_CON
    SZ_ph->L_map(*map.SZ_ph,*obj.SZ_ph);
#endif

#ifdef __T1_CON
    SZ_dp->L_map(*map.SZ_dp,*obj.SZ_dp);
#endif

#ifdef __T2_CON
    SZ_pph->L_map(*map.SZ_pph,*obj.SZ_pph);
#endif
}

void SUP::daxpy(double alpha,SUP &SZ)
{
    SZ_tp[0]->daxpy(alpha,*SZ.SZ_tp[0]);
    SZ_tp[1]->daxpy(alpha,*SZ.SZ_tp[1]);

#ifdef __G_CON
    SZ_ph->daxpy(alpha,*SZ.SZ_ph);
#endif

#ifdef __T1_CON
    SZ_dp->daxpy(alpha,*SZ.SZ_dp);
#endif

#ifdef __T2_CON
    SZ_pph->daxpy(alpha,*SZ.SZ_pph);
#endif
}

double SUP::trace()
{
   double brecht;

   brecht = SZ_tp[0]->trace();
   brecht += SZ_tp[1]->trace();

#ifdef __G_CON
   brecht += SZ_ph->trace();
#endif

#ifdef __T1_CON
   brecht += SZ_dp->trace();
#endif

#ifdef __T2_CON
   brecht += SZ_pph->trace();
#endif

   return brecht;
}

void SUP::fill(TPM &tpm)
{
    *SZ_tp[0] = tpm;
    SZ_tp[1]->Q(tpm);

#ifdef __G_CON
    SZ_ph->G(tpm);
#endif

#ifdef __T1_CON
    SZ_dp->T(tpm);
#endif

#ifdef __T2_CON
    SZ_pph->T(tpm);
#endif
}

TPM &SUP::tpm(int i)
{
    return *SZ_tp[i];
}

int SUP::gn_tp()
{
    return n_tp;
}

#ifdef __G_CON
PHM &SUP::phm()
{
    return *SZ_ph;
}
int SUP::gn_ph()
{
    return n_ph;
}
#endif

#ifdef __T1_CON
DPM &SUP::dpm()
{
    return *SZ_dp;
}
int SUP::gn_dp()
{
    return n_dp;
}
#endif

#ifdef __T2_CON
PPHM &SUP::pphm()
{
    return *SZ_pph;
}
int SUP::gn_pph()
{
    return n_pph;
}
#endif

