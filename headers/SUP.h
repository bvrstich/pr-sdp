#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "TPM.h"

#ifndef PQ

#include "PHM.h"

#endif

#ifdef T_1

#include "DPM.h"

#endif

class EIG;

class SUP{
  
   friend ostream &operator<<(ostream &,const SUP &);

   public:

      //constructor
      SUP(int M,int N);

      //copy constructor
      SUP(SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(SUP &);

      //overload -= operator
      SUP &operator-=(SUP &);

      //overload equality operator
      SUP &operator=(SUP &);

      //overload equality operator
      SUP &operator=(double &);

      SUP operator*(SUP &);

      TPM &tpm(int i);

      int gN();

      int gM();

      int gn_tp();

#ifndef PQ

      PHM &phm();

      int gn_ph();

#endif

#ifdef T_1
      
      DPM &dpm();

      int gn_dp();

#endif

      double ddot(SUP &);

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP &,SUP &);

      void daxpy(double alpha,SUP &);

      double trace();

      void fill(TPM &);

   private:

      TPM **SZ_tp;

      int M;
      int N;

      int n_tp;

#ifndef PQ

      PHM *SZ_ph;

      int n_ph;

#endif

#ifdef T_1
   
      DPM *SZ_dp:

      int n_dp;

#endif

};

#endif
