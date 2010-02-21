#ifndef SUP_PQ_H
#define SUP_PQ_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "../TPM.h"

class EIG;

//basisklasse
class SUP_PQ{
  
   friend ostream &operator<<(ostream &,SUP_PQ &);

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

      TPM **SZ_tp;

      int M;
      int N;

      int n_tp;

      int dim;

};

#endif
