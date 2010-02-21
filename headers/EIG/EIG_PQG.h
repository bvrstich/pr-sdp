#ifndef EIG_PQG_H
#define EIG_PQG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "../SUP.h"
#include "EIG_PQ.h"

class EIG_PQG : public EIG_PQ {

   public:
      
      //constructor
      EIG_PQG(int M,int N);

      //copy constructor
      EIG_PQG(EIG_PQG &);

      //constructor met initialisatie op 
      EIG_PQG(SUP &);

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

      double *eig_ph;

      int n_ph;//dim ph space

};

#endif
