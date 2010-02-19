#ifndef EIG_PQ_H
#define EIG_PQ_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SUP.h"

//basisklasse
class EIG_PQ{

   friend ostream &operator<<(ostream &,const EIG_PQ &);

   public:
      
      //constructor
      EIG_PQ(int M,int N);

      //copy constructor
      EIG_PQ(EIG_PQ &);

      //constructor met initialisatie op 
      EIG_PQ(SUP &);

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

      double **eig;

      int N;//nr of particles
      int M;//dim sp space
      int n_tp;//dim tp space

      int dim;//totale dimensie

};

#endif
