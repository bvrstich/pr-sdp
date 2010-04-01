
#include "headers/EIG.h"

ostream &operator<<(ostream &output,EIG &eig)
{
    for(int i = 0;i < eig.n_tp;i++)
	std::cout << i << "\t" << eig.eig[0][i] << std::endl;

    output << std::endl;

    for(int i = 0;i < eig.n_tp;i++)
	output << i << "\t" << eig.eig[1][i] << std::endl;

#ifdef __G_CON
    output << std::endl;

    for(int i = 0;i < eig.n_ph;i++)
	output << i << "\t" << eig.eig_ph[i] << std::endl;
#endif

#ifdef __T1_CON
    output << std::endl;

    for(int i = 0;i < eig.n_dp;i++)
	output << i << "\t" << eig.eig_dp[i] << std::endl;
#endif

#ifdef __T2_CON
    output << std::endl;

    for(int i = 0;i < eig.n_pph;i++)
	output << i << "\t" << eig.eig_pph[i] << std::endl;
#endif

    return output;
}

//constructor
EIG::EIG(int M,int N)
{
    this->N=N;
    this->M=M;

    n_tp=M*(M-1)/2;
    dim=2*n_tp;

    eig = new double *[2];

    eig[0] = new double [dim];
    eig[1] = eig[0] + n_tp;

#ifdef __G_CON
    n_ph=M*M;
    dim+=n_ph;
    eig_ph=new double[n_ph];
#endif

#ifdef __T1_CON
    n_dp=M*(M-1)*(M-2)/6;
    dim+=n_dp;
    eig_dp=new double[n_dp];
#endif

#ifdef __T2_CON
    n_pph=M*M*(M-1)/2;
    dim+=n_pph;
    eig_pph=new double[n_pph];
#endif
}

//copy constructor
EIG::EIG(EIG &eig_c)
{
    N=eig_c.N;
    M=eig_c.M;

    n_tp=eig_c.n_tp;
    dim=2*n_tp;

    eig = new double *[2];

    eig[0] = new double [dim];
    eig[1] = eig[0] + n_tp;

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

#ifdef __G_CON
    n_ph=eig_c.n_ph;
    dim+=n_ph;
    eig_ph=new double[n_ph];

    dcopy_(&n_ph,eig_c.eig_ph,&inc,eig_ph,&inc);
#endif

#ifdef __T1_CON
    n_dp=eig_c.n_dp;
    dim+=n_dp;
    eig_dp=new double[n_dp];

    dcopy_(&n_dp,eig_c.eig_dp,&inc,eig_dp,&inc);
#endif

#ifdef __T2_CON
    n_pph=eig_c.n_pph;
    dim+=n_pph;
    eig_pph=new double[n_pph];

    dcopy_(&n_pph,eig_c.eig_pph,&inc,eig_pph,&inc);
#endif
}

//constructor met initialisatie op 
EIG::EIG(SUP &SZ)
{
    N=SZ.gN();
    M=SZ.gM();

    n_tp=SZ.gn_tp();
    dim=2*n_tp;

    eig = new double *[2];

    eig[0] = new double [dim];
    eig[1] = eig[0] + n_tp;

    SZ.tpm(0).diagonalize(eig[0]);
    SZ.tpm(1).diagonalize(eig[1]);

#ifdef __G_CON
    n_ph=SZ.gn_ph();
    dim+=n_ph;
    eig_ph=new double[n_ph];
    SZ.phm().diagonalize(eig_ph);
#endif

#ifdef __T1_CON
    n_dp=SZ.gn_dp();
    dim+=n_dp;
    eig_dp=new double[n_dp];
    SZ.dpm().diagonalize(eig_dp);
#endif

#ifdef __T2_CON
    n_pph=SZ.gn_pph();
    dim+=n_pph;
    eig_pph=new double[n_pph];
    SZ.pphm().diagonalize(eig_pph);
#endif
}

//destructor
EIG::~EIG()
{
    delete [] eig[0];
    delete [] eig;

#ifdef __G_CON
    delete [] eig_ph;
#endif

#ifdef __T1_CON
    delete [] eig_dp;
#endif

#ifdef __T2_CON
    delete [] eig_pph;
#endif
}

int EIG::gN()
{
    return N;
}

int EIG::gM()
{
    return M;
}

//overload equality operator
EIG &EIG::operator=(EIG &eig)
{
   int inc = 1;

   dcopy_(&dim,eig.eig[0],&inc,eig[0],&inc);

   return *this;
}

double *EIG::operator[](int i)
{
    switch(i)
    {
	case 0:
	    return eig[0];
	case 1:
	    return eig[1];
#ifdef __G_CON
	case 2:
	    return eig_ph;
#endif
#ifdef __T1_CON
	case 3:
	    return eig_dp;
#endif
#ifdef __T2_CON
	case 4:
	    return eig_pph;
#endif
	default:
	    std::cout << "Error: invalid block" << std::endl;
    }
    return 0;
}

//acces to the numbers
double EIG::operator()(int block,int index)
{
    //return eig[block][index];
    switch(block)
    {
	case 0:
	    return eig[0][index];
	case 1:
	    return eig[1][index];
#ifdef __G_CON
	case 2:
	    return eig_ph[index];
#endif
#ifdef __T1_CON
	case 3:
	    return eig_dp[index];
#endif
#ifdef __T2_CON
	case 4:
	    return eig_pph[index];
#endif
	default:
	    std::cout << "Error: invalid block" << std::endl;
    }
    return -1.0;
}

double EIG::lsfunc(double alpha)
{
   double brecht = 0.0;

   for(int i = 0;i < n_tp;i++)
      brecht += eig[0][i]/(1.0 + alpha*eig[0][i]);

   for(int i = 0;i < n_tp;i++)
      brecht += eig[1][i]/(1.0 + alpha*eig[1][i]);

#ifdef __G_CON
   for(int i = 0;i < n_ph;i++)
      brecht += eig_ph[i]/(1.0 + alpha*eig_ph[i]);
#endif

#ifdef __T1_CON
   for(int i = 0;i < n_dp;i++)
      brecht += eig_dp[i]/(1.0 + alpha*eig_dp[i]);
#endif

#ifdef __T2_CON
   for(int i = 0;i < n_pph;i++)
      brecht += eig_pph[i]/(1.0 + alpha*eig_pph[i]);
#endif

   return brecht;
}

double EIG::min()
{
    double brecht=eig[0][0];

    if(brecht>eig[1][0])
	brecht=eig[1][0];

#ifdef __G_CON
    if(brecht>eig_ph[0])
	brecht=eig_ph[0];
#endif

#ifdef __T1_CON
    if(brecht>eig_dp[0])
	brecht=eig_dp[0];
#endif

#ifdef __T2_CON
    if(brecht>eig_pph[0])
	brecht=eig_pph[0];
#endif

    return brecht;
}

int EIG::gn_tp()
{
    return n_tp;
}

#ifdef __G_CON
int EIG::gn_ph()
{
    return n_ph;
}
#endif

#ifdef __T1_CON
int EIG::gn_dp()
{
    return n_dp;
}
#endif

#ifdef __T2_CON
int EIG::gn_pph()
{
    return n_pph;
}
#endif
