#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SUP.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, EIG is a "block"-vector over the carrierspace's of the P,Q,G,T1,T2 
 * conditions (depending on the macro's that are set). It contains room to store the eigenvalues
 * of all conditions, and special member function that work with these eigenvalues.
 * This class should only be used when a SUP matrix has been diagonalized, some functions 
 * could give strange results when the EIG object is filled with random numbers.\n\n
 */
class EIG
{
    /**
     * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
     * ifstream object and type:\n\n
     * object << eig_p << endl;\n\n
     * For output onto the screen type: \n\n
     * cout << eig_p << endl;\n\n
     * @param output The stream to which you are writing (e.g. cout)
     * @param eig_p the EIG you want to print
     */
    friend ostream &operator<<(ostream &output,EIG &eig_p);

    public:
    //constructor
    EIG(int M,int N);

    //copy constructor
    EIG(EIG &);

    //constructor met initialisatie op 
    EIG(SUP &);

    //destructor
    ~EIG();

    int gN();

    int gM();

    //overload equality operator
    EIG &operator=(EIG &);

    double *operator[](int);

    //acces to the numbers
    double operator()(int block,int index);

    double lsfunc(double );

    double min();

    int gn_tp();

#ifdef __G_CON
    int gn_ph();
#endif

#ifdef __T1_CON
    int gn_dp();
#endif

#ifdef __T2_CON
    int gn_pph();
#endif

    private:

    //!double pointer that will store the eigenvalues of the P (eig[0][0 -> n_tp - 1]) and Q (eig[1][0 -> n_tp - 1]) blocks of the SUP_PQ matrix
    double **eig;

    //!nr of particles
    int N;

    //!dim sp space
    int M;

    //!total dimension
    int dim;

    //!dim tp space
    int n_tp;

#ifdef __G_CON
    //! dimension of the ph space
    int n_ph;
    //! double pointer to PHM
    double *eig_ph;
#endif

#ifdef __T1_CON
    //! dimension of the dp space
    int n_dp;
    double *eig_dp;
#endif

#ifdef __T2_CON
    //! dimension of the pph space
    int n_pph;
    double *eig_pph;
#endif
};

#endif /* EIG_H */

/* vim: set ts=3 sw=3 tw=3 expandtab :*/
