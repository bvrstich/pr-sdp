#ifndef EIG_PQGT1_H
#define EIG_PQGT1_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "../SUP/SUP_PQGT1.h"
#include "EIG_PQG.h"

class EIG_PQGT1 : public EIG_PQG {

   friend ostream &operator<<(ostream &,EIG_PQGT1 &);

   public:
      
      //constructor
      EIG_PQGT1(int M,int N);

      //copy constructor
      EIG_PQGT1(EIG_PQGT1 &);

      //constructor met initialisatie op 
      EIG_PQGT1(SUP_PQGT1 &);

      //destructor
      ~EIG_PQGT1();

      int gn_dp();

      //overload equality operator
      EIG_PQGT1 &operator=(EIG_PQGT1 &);

      double *operator[](int);

      //acces to the numbers
      double operator()(int block,int index);

      double lsfunc(double );

      double min();

   protected:

      double *eig_dp;

      int n_dp;//dim three particle space

};

#endif
