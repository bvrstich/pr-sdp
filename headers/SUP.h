#ifndef SUP_H
#define SUP_H

#include "include.h"

class EIG;

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, SUP is a blockmatrix over the carrierspace's of the P and Q conditions. This is the
 * motherclass of all the SUP* classes because I assume that the P and Q condtional will always be used.
 * This class contains two TPM objects, that are independent of each other (by which I mean that TPM::Q(SUP::tpm (0))
 * is not neccesarily equal to SUP::tpm (1)).
 */
class SUP {

    /**
     * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
     * ifstream object and type:\n\n
     * object << sup_p << endl;\n\n
     * For output onto the screen type: \n\n
     * cout << sup_p << endl;\n\n
     * @param output The stream to which you are writing (e.g. cout)
     * @param sup_p the SUP you want to print
     */
    friend ostream &operator<<(ostream &output,SUP &sup_p);

    public:

    //constructor
    SUP(int M,int N);

    //copy constructor
    SUP(SUP &);

    //destructor
    ~SUP();

    //overload += operator
    SUP &operator+=(SUP &);

    //overload -= operator
    SUP &operator-=(SUP &);

    //overload equality operator
    SUP &operator=(SUP &);

    //overload equality operator
    SUP &operator=(double &);

    SUP operator*(SUP &);

    int gN();

    int gM();

    double ddot(SUP &);

    void invert();

    void dscal(double alpha);

    //positieve of negatieve vierkantswortel uit een supermatrix
    void sqrt(int option);

    void L_map(SUP &,SUP &);

    void daxpy(double alpha,SUP &);

    double trace();

    void fill(TPM &);

    TPM &tpm(int i);

    int gn_tp();

#ifdef __G_CON
    PHM &phm();
    int gn_ph();
#endif

#ifdef __T1_CON
    DPM &dpm();
    int gn_dp();
#endif

#ifdef __T2_CON
    PPHM &pphm();
    int gn_pph();
#endif

#ifdef __T2P_CON
    T2PM &t2pm();
    int gn_t2p();
#endif

    private:

    //!double pointer to TPM's. will containt the P space matrix in SZ_tp[0] and the Q space matrix in SZ_tp[1]
    TPM **SZ_tp;

    //!nr of sp orbitals
    int M;

    //!nr of particles
    int N;

    //!total dimension of the block matrix
    int dim;

    //!dimension of tp space
    int n_tp;

#ifdef __G_CON
    //! dimension of the ph space
    int n_ph;
    //! double pointer to PHM
    PHM *SZ_ph;
#endif

#ifdef __T1_CON
    //! dimension of the dp space
    int n_dp;
    //! double pointer to DPM
    DPM *SZ_dp;
#endif

#ifdef __T2_CON
    //! dimension of the pph space
    int n_pph;
    //! double pointer to PPHM
    PPHM *SZ_pph;
#endif

#ifdef __T2P_CON
    //! dimension of the T2' space
    int n_t2p;
    //! double pointer to T2PM
    T2PM *SZ_t2p;
#endif
};

#endif /* SUP_H */

/* vim: set ts=3 sw=3 tw=3 expandtab :*/
