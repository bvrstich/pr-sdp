#ifndef SUP_PQGT1_H
#define SUP_PQGT1_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "../TPM.h"
#include "../DPM.h"
#include "SUP_PQG.h"

class EIG;

//neemt veel over van SUP_PQG
class SUP_PQGT1 : public SUP_PQG {

   friend ostream &operator<<(ostream &,SUP_PQGT1 &);
  
   public:

      //constructor
      SUP_PQGT1(int M,int N);

      //copy constructor
      SUP_PQGT1(SUP_PQGT1 &);

      //destructor
      ~SUP_PQGT1();

      //overload += operator
      SUP_PQGT1 &operator+=(SUP_PQGT1 &);

      //overload -= operator
      SUP_PQGT1 &operator-=(SUP_PQGT1 &);

      //overload equality operator
      SUP_PQGT1 &operator=(SUP_PQGT1 &);

      //overload equality operator
      SUP_PQGT1 &operator=(double &);

      SUP_PQGT1 operator*(SUP_PQGT1 &);

      int gn_dp();

      DPM &dpm();

      double ddot(SUP_PQGT1 &);

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP_PQGT1 &,SUP_PQGT1 &);

      void daxpy(double alpha,SUP_PQGT1 &);

      double trace();

      void fill(TPM &);

   protected:

      DPM *SZ_dp;

      int n_dp;

};

#endif
