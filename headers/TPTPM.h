#ifndef TPTPM_H
#define TPTPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class PHM;

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class TPTPM is a class written for matrices of two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */
class TPTPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param hess_p the TPTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPTPM &hess_p);

   public:
      
      //constructor
      TPTPM();

      //copy constructor
      TPTPM(const TPTPM &);

      //destructor
      virtual ~TPTPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int z,int t,int h) const;

      //access to the numbers in tp mode
      double operator()(int I,int J,int K,int L) const;

      void dp(const PHM &);

      void dpt2(const DPM &);

      void dpw2(const PPHM &);

      void dptw(const PPHM &);

      void dpt2(const PPHM &);

      void I(const TPM &);

      void Q(const TPM &);

      void G(const PHM &);

      void T(const DPM &);

      void T(const PPHM &);

      static int gn();

      static int gtpmm2t(int,int);

      static int gt2tpmm(int,int);
      
      static void init(int,int);

      static void clear();

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!list relating the single-particle space to the TPTPM basis
      static vector< vector<int> > tpmm2t;

      //!list relating the single-particle space to the TPTPM basis
      static int **t2tpmm;

};

#endif
