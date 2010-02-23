#ifndef SUP_PQ_H
#define SUP_PQ_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "../TPM.h"

class EIG;

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, SUP_PQ is a blockmatrix over the carrierspace's of the P and Q conditions. This is the
 * motherclass of all the SUP_PQ* classes because I assume that the P and Q condtional will always be used.
 * This class contains two TPM objects, that are independent of each other (by which I mean that TPM::Q(SUP_PQ::tpm (0))
 * is not neccesarily equal to SUP_PQ::tpm (1)).
 */
class SUP_PQ{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param sup_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,SUP_PQ &sup_p);

   public:

      //constructor
      SUP_PQ(int M,int N);

      //copy constructor
      SUP_PQ(SUP_PQ &);

      //destructor
      ~SUP_PQ();

      //overload += operator
      SUP_PQ &operator+=(SUP_PQ &);

      //overload -= operator
      SUP_PQ &operator-=(SUP_PQ &);

      //overload equality operator
      SUP_PQ &operator=(SUP_PQ &);

      //overload equality operator
      SUP_PQ &operator=(double &);

      SUP_PQ operator*(SUP_PQ &);

      TPM &tpm(int i);

      int gN();

      int gM();

      int gn_tp();

      double ddot(SUP_PQ &);

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP_PQ &,SUP_PQ &);

      void daxpy(double alpha,SUP_PQ &);

      double trace();

      void fill(TPM &);

   protected:

      //!double pointer to TPM's. will containt the P space matrix in SZ_tp[0] and the Q space matrix in SZ_tp[1]
      TPM **SZ_tp;

      //!nr of sp orbitals
      int M;

      //!nr of particles
      int N;

      //!dimension of tp space
      int n_tp;

      //!total dimension of the block matrix
      int dim;

};

#endif
