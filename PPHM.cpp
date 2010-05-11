#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;
    using std::endl;

#include "headers/include.h"

int PPHM::counter = 0;

int **PPHM::pph2s;
int ***PPHM::s2pph;

/**
 * standard constructor: constructs Matrix object of dimension M^2*(M - 1)/2 and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
PPHM::PPHM(int M,int N) : Matrix(M*M*(M - 1)/2)
{
    this->N = N;
    this->M = M;
    this->n = M*M*(M - 1)/2;

    if(counter == 0)
    {
	//allocatie van s2pph
	s2pph = new int ** [M];

	for(int i = 0;i < M;i++)
	{
	    s2pph[i] = new int * [M];

	    for(int j = 0;j < M;j++)
		s2pph[i][j] = new int [M];
	}

	//allocatie van pph2s
	pph2s = new int * [n];

	for(int i = 0;i < n;i++)
	    pph2s[i] = new int [3];

	//initialisatie van de twee arrays
	int teller = 0;

	for(int a = 0;a < M;a++)
	    for(int b = a + 1;b < M;b++)
		for(int c = 0;c < M;c++)
		{
		    s2pph[a][b][c] = teller;

		    pph2s[teller][0] = a;
		    pph2s[teller][1] = b;
		    pph2s[teller++][2] = c;
		}
    }
    counter++;
}

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)*(M - 2)/6 and copies the content of PPHM_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param PPHM_c input PPHM to be copied
 */
PPHM::PPHM(PPHM &PPHM_c) : Matrix(PPHM_c)
{
    this->N = PPHM_c.N;
    this->M = PPHM_c.M;
    this->n = M*M*(M - 1)/2;

    if(counter == 0)
    {
	//allocatie van s2pph
	s2pph = new int ** [M];

	for(int i = 0;i < M;i++)
	{
	    s2pph[i] = new int * [M];

	    for(int j = 0;j < M;j++)
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
		for(int c = 0;c < M;++c)
		{
		    s2pph[a][b][c] = teller;

		    pph2s[teller][0] = a;
		    pph2s[teller][1] = b;
		    pph2s[teller++][2] = c;
		}
    }
    counter++;
}

/**
 * destructor: if counter == 1 the memory for the static lists pph2s en s2pph twill be deleted.
 */
PPHM::~PPHM()
{
    if(counter == 1)
    {
	for(int i = 0;i < M;++i)
	{
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
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * PPHM(a,b,c,d,e,f) = -PPHM(b,a,c,d,e,f) = ...\n\n
 * PPHM(a,a,c,d,e,f) = 0\n\n
 * @param a first sp index that forms the pph row index i together with b and c
 * @param b second sp index that forms the pph row index i together with a and c
 * @param c third sp index that forms the pph row index i together with a and b
 * @param d first sp index that forms the pph column index j together with e and z
 * @param e second sp index that forms the pph column index j together with d and z
 * @param z third sp index that forms the pph column index j together with d and e
 * @return the number on place PPHM(i,j) with the right phase.
 */
double PPHM::operator()(int a,int b,int c,int d,int e,int z) const
{
    //eerst kijken of er geen indices gelijk zijn:
    if(a == b || d == e)
	return 0;

    //dan kijken wel pph index met welke fase moet genomen worden:
    int i,j;

    int phase = 1;

    if(a < b)
	i = s2pph[a][b][c];
    else
    {
	i = s2pph[b][a][c];
	phase *= -1;
    }

    if(d < e)
	j = s2pph[d][e][z];
    else
    {
	j = s2pph[e][d][z];
	phase *= -1;
    }

    return phase*(*this)(i,j);
}

ostream &operator<<(ostream &output,const PPHM &PPHM_p)
{
    for(int i = 0;i < PPHM_p.n;++i)
	for(int j = 0;j < PPHM_p.n;++j)
	{
	    output << i << "\t" << j << "\t|\t" << PPHM_p.pph2s[i][0] << "\t" << PPHM_p.pph2s[i][1] << "\t" << PPHM_p.pph2s[i][2]
		<< "\t" << PPHM_p.pph2s[j][0] << "\t" << PPHM_p.pph2s[j][1] << "\t" << PPHM_p.pph2s[j][2] 
		<< "\t" << PPHM_p(i,j) << endl;
	}

    return output;
}

/**
 * @return nr of particles
 */
int PPHM::gN()
{
    return N;
}

/**
 * @return dimension of sp space
 */
int PPHM::gM()
{
    return M;
}

/**
 * @return dimension of dp space and of Matrix
 */
int PPHM::gn()
{
    return n;
}

/**
 * The T2-map: maps a TPM object (tpm) on a PPHM object (*this)
 * @param tpm input TPM
 */
void PPHM::T(TPM &tpm)
{
    SPM spm(tpm);

    int a,b,c,d,e,z;

    for(int i = 0;i < n;i++)
    {
	a = pph2s[i][0];
	b = pph2s[i][1];
	c = pph2s[i][2];

	for(int j = i;j < n;j++)
	{
	    d = pph2s[j][0];
	    e = pph2s[j][1];
	    z = pph2s[j][2];

	    (*this)(i,j) = 0;

	    // sp terms
	    if( a == d && b == e )
		(*this)(i,j) += spm(c,z);

	    if( a == e && b == d )
		(*this)(i,j) -= spm(c,z);

	    // tp terms
	    if( a == d )
		(*this)(i,j) -= tpm(z,b,c,e);

	    if( b == d )
		(*this)(i,j) += tpm(z,a,c,e);

	    if( c == z )
		(*this)(i,j) += tpm(a,b,d,e);

	    if( b == e )
		(*this)(i,j) -= tpm(z,a,c,d);

	    if( a == e )
		(*this)(i,j) += tpm(z,b,c,d);

	    (*this)(j,i) = (*this)(i,j);
	}
    }
}

/* vim: set ts=3 sw=3 tw=3 expandtab :*/
