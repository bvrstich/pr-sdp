#ifndef EIG_PQG_H
#define EIG_PQG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "../SUP/SUP_PQG.h"
#include "EIG_PQ.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, EIG_PQG is a "block"-vector over the carrierspace's of the P, Q and G conditions. It inherits from EIG_PQ and expands it with
 * a vector over the G-space (of dimension n_ph). Some functions of EIG_PQ are reimplemented to include the G-condition.
 * This class should only be used when a SUP_PQG matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 */
class EIG_PQG : public EIG_PQ {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << eig_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << eig_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,EIG_PQG &eig_p);

   public:
      
      //constructor
      EIG_PQG(int M,int N);

      //copy constructor
      EIG_PQG(EIG_PQG &);

      //constructor met initialisatie op 
      EIG_PQG(SUP_PQG &);

      //destructor
      ~EIG_PQG();

      int gn_ph();

      //overload equality operator
      EIG_PQG &operator=(EIG_PQG &);

      double *operator[](int);

      //acces to the numbers
      double operator()(int block,int index);

      double lsfunc(double );

      double min();

   protected:

      //!pointer to double, will contain the eigenvector of the G part of a SUP_PQG matrix eig[0 -> n_ph - 1]
      double *eig_ph;

      //!dimension of particle hole space
      int n_ph;

};

#endif
