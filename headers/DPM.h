#ifndef DPM_H
#define DPM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

class DPM : public Matrix {

   friend ostream &operator<<(ostream &,const DPM &);

   public:
      
      //constructor
      DPM(int M,int N);

      //copy constructor
      DPM(DPM &);

      //destructor
      virtual ~DPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int f) const;

      //geef N terug
      int gN();

      //geef M terug
      int gM();

      //geef dim terug
      int gn();

      //maak een DPM van een TPM
      void T(TPM &);

   private:

      //counter of objects
      static int counter;

      //relationship bewtween dp and sp basis
      static int **dp2s;
      static int ***s2dp;

      int N;//nr of particles
      int M;//dim sp space
      int n;//dim dpm space

};

#endif
