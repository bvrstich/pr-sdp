#ifndef EIG_PQGT2_H
#define EIG_PQGT2_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "../SUP/SUP_PQGT2.h"
#include "EIG_PQG.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, EIG_PQGT2 is a "block"-vector over the carrierspace's of the P, Q, G and T1 conditions. It inherits from EIG_PQG and expands it with
 * a vector over the T1-space (of dimension n_dp). Some functions of EIG_PQG are reimplemented to include for the T1-condition.
 * This class should only be used when a SUP_PQGT2 matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 */
class EIG_PQGT2 : public EIG_PQG {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << eig_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << eig_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,EIG_PQGT2 &eig_p);

   public:
      
      //constructor
      EIG_PQGT2(int M,int N);

      //copy constructor
      EIG_PQGT2(EIG_PQGT2 &);

      //constructor met initialisatie op 
      EIG_PQGT2(SUP_PQGT2 &);

      //destructor
      ~EIG_PQGT2();

      int gn_pph();

      //overload equality operator
      EIG_PQGT2 &operator=(EIG_PQGT2 &);

      double *operator[](int);

      //acces to the numbers
      double operator()(int block,int index);

      double lsfunc(double );

      double min();

   protected:

      //!pointer to the eigenvalues of the T1 part of SUP_PQGT2 object, containt the number in: eig_dp[0 -> n_dp - 1]
      double *eig_pph;

      //!dimension of three particle space
      int n_pph;

};

#endif
