#ifndef SUP_PQG_H
#define SUP_PQG_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "TPM.h"
#include "PHM.h"
#include "SUP_PQ.h"

class EIG;

//neemt veel over van SUP_PQ
class SUP_PQG : public SUP_PQ {
  
   public:

      //constructor
      SUP_PQG(int M,int N);

      //copy constructor
      SUP_PQG(SUP_PQG &);

      //destructor
      ~SUP_PQG();

      //overload += operator
      SUP_PQG &operator+=(SUP_PQG &);

      //overload -= operator
      SUP_PQG &operator-=(SUP_PQG &);

      //overload equality operator
      SUP_PQG &operator=(SUP_PQG &);

      //overload equality operator
      SUP_PQG &operator=(double &);

      SUP_PQG operator*(SUP_PQG &);

      PHM &phm();

      int gn_ph();

      double ddot(SUP_PQG &);

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP_PQG &,SUP_PQG &);

      void daxpy(double alpha,SUP_PQG &);

      double trace();

      void fill(TPM &);

   protected:

      PHM *SZ_ph;

      int n_ph;

};

#endif
