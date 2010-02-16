#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;
using std::string;

#include "Matrix.h"

class SPM;
class SUP;

#ifndef PQ

class PHM;

#endif

#ifdef T_1

class DPM;

#endif

class TPM : public Matrix {

   friend ostream &operator<<(ostream &,const TPM &);

   public:
      
      //constructor
      TPM(int M,int N);

      //copy constructor
      TPM(TPM &);

      //destructor
      virtual ~TPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d) const;

      //geef N terug
      int gN();

      //geef N terug
      int gM();

      //geef n terug
      int gn();

      void hubbard(double U);

      void Q(TPM &);

#ifndef PQ

      void G(PHM &);

#endif

#ifdef T_1

      void bar(DPM &);
      
      void T(DPM &);

#endif

      void init();

      void proj_Tr();

      void constr_grad(double t,TPM &,SUP &);

      int solve(double t,SUP &,TPM &);

      double line_search(double t,SUP &P,TPM &ham);

      double line_search(double t,TPM &,TPM &);

      void H(double t,TPM &b,SUP &P);

   private:

      static int **t2s;
      static int **s2t;

      static int counter;

      int N;//nr of particles
      int M;//dim sp space
      int n;//dim tp space

};

#endif
