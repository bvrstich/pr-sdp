#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SUP.h"

class EIG{

   friend ostream &operator<<(ostream &,const EIG &);

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

      int gn_tp();

#ifndef PQ

      int gn_ph();

#endif

      //overload equality operator
      EIG &operator=(EIG &);

      double *operator[](int);

      //acces to the numbers
      double &operator()(int block,int index);

      double lsfunc(double );

      double min();

   private:

      double **eig;

      int N;//nr of particles
      int M;//dim sp space
      int n_tp;//dim tp space

      int n;//totale dimensie

#ifndef PQ

      int n_ph;//dim ph space

#endif

};

#endif
