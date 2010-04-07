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

/**
 * Standard constructor, allocates the matrices.
 * @param M dim of sp space
 * @param N nr of particles
 */

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


/**
 * Copy constructor, allocates the matrices and copies the content of SZ into it.
 * @param SZ input SUP
 */
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

/**
 * Destructor, deallocates memory in matrices.
 */
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

/**
 * += operator overloaded
 * @param SZ SUP to add to this
 */
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

/**
 * -= operator overloaded
 * @param SZ SUP to be subtracted from this
 */
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

/**
 * Overload equality operator, copy SZ into this
 * @param SZ SUP to be copied into this
 */
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

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. SZ = 0 makes all the Matrix elements zero.
 * @param a the number
 */
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

/**
 * @return particle number
 */
int SUP::gN()
{
    return N;
}

/**
 * @return dim of sp space
 */
int SUP::gM()
{
    return M;
}

/**
 * @param SZ input SUP SZ
 * @return inproduct between this and input matrix SZ, defined as Tr(this SZ)
 */
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

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 */
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

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
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

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
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

/**
 * Multiply symmetric SUP blockmatrix object left en right with symmetric SUP blockmatrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map SUP that will be multiplied to the left en to the right of matrix object
 * @param obj central SUP
 */
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

/**
 * add the SUP SZ times the constant alpha to this
 * @param alpha the constant to multiply the SZ with
 * @param SZ the SUP to be multiplied by alpha and added to (*this)
 */
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

/**
 * @return trace of SUP, which is defined as the sum of the traces of the individual matrices.
 */
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

/**
 * Fill the SUP with TPM tpm, which means put:\n\n
 * SZ_tp[0] = tpm\n
 * SZ_tp[1] = TPM::Q (tpm)\n
 * SZ_ph = TPM::G (tpm)\n
 * SZ_dp = DPM::T (tpm)\n
 * SZ_pph = PPHM::T (tpm)\n
 * @param tpm input TPM
 */
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

/**
 * @param i which block you want to have the pointer to.
 * @return pointer to the individual TPM blocks: SZ_tp[i]
 */
TPM &SUP::tpm(int i)
{
    return *SZ_tp[i];
}

/**
 * @return dim of tp space
 */
int SUP::gn_tp()
{
    return n_tp;
}

#ifdef __G_CON
/**
 * @return pointer to the PHM block
 */
PHM &SUP::phm()
{
    return *SZ_ph;
}

/**
 * @return dim of ph space
 */
int SUP::gn_ph()
{
    return n_ph;
}
#endif

#ifdef __T1_CON
/**
 * @return pointer to the DPM block
 */
DPM &SUP::dpm()
{
    return *SZ_dp;
}

/**
 * @return dim of dp space
 */
int SUP::gn_dp()
{
    return n_dp;
}
#endif

#ifdef __T2_CON
/**
 * @return pointer to the PPHM block
 */
PPHM &SUP::pphm()
{
    return *SZ_pph;
}

/**
 * @return dim of pph space
 */
int SUP::gn_pph()
{
    return n_pph;
}
#endif

