#ifndef EIG_PQ_H
#define EIG_PQ_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "../SUP/SUP_PQ.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, EIG_PQ is a "block"-vector over the carrierspace's of the P and Q conditions. It contains room
 * to store the eigenvalues of a P and a Q block, and special member function that work with these eigenvalues.
 * This class should only be used when a SUP_PQ matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 * This function is the mother class of all EIG_PQ* classes.
 */
class EIG_PQ{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << eig_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << eig_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,EIG_PQ &eig_p);

   public:
      
      //constructor
      EIG_PQ(int M,int N);

      //copy constructor
      EIG_PQ(EIG_PQ &);

      //constructor met initialisatie op 
      EIG_PQ(SUP_PQ &);

      //destructor
      ~EIG_PQ();

      int gN();

      int gM();

      int gn_tp();

      //overload equality operator
      EIG_PQ &operator=(EIG_PQ &);

      double *operator[](int);

      //acces to the numbers
      double operator()(int block,int index);

      double lsfunc(double );

      double min();

   protected:

      //!double pointer that will store the eigenvalues of the P (eig[0][0 -> n_tp - 1]) and Q (eig[1][0 -> n_tp - 1]) blocks of the SUP_PQ matrix
      double **eig;

      //!nr of particles
      int N;

      //!dim sp space
      int M;

      //!dim tp space
      int n_tp;

      //!total dimension
      int dim;

};

#endif
