#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"
#include "PPHM.h"

class PPHM;

/**
 * @author Brecht Verstichel
 * @date 22-02-2010\n\n
 * This class SPM is a class written for single particle matrices, it inherits alle the function from its mother 
 * Matrix and some member function which are special for SPM's.
 */
class SPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << spm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << spm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM you want to print
    */
   friend ostream &operator<<(ostream &output,SPM &spm_p);

   public:

      //constructor
      SPM(int M,int N);

      //copy constructor
      SPM(SPM &);

      //destructor
      virtual ~SPM();

      using Matrix::operator=;

      int gN();

      int gM();

      /**
       * function constructs SPM from TPM or PHM, this function is SPM::bar times 1/(N - 1)
       * @param MT input TPM or PHM
       */
      template<class MatrixType>
         void constr(MatrixType &MT){

            double ward;

            ward = 1.0/(N - 1.0);

            for(int a = 0;a < M;++a)
               for(int b = a;b < M;++b){

                  (*this)(a,b) = 0.0;

                  for(int l = 0;l < M;++l)
                     (*this)(a,b) += MT(a,l,b,l);

                  (*this)(a,b) *= ward;

                  (*this)(b,a) = (*this)(a,b);
               }

         }

      /**
       * Constructor of SPM with initialization on SPM::constr of input TPM or PHM
       * @param MT input TPM or PHM
       */
      template<class MatrixType>
         SPM(MatrixType &MT) : Matrix(MT.gM()) {

            this->M = MT.gM();
            this->N = MT.gN();

            this->constr(MT);

         }

      /**
       * construct bar matrix from TPM or PHM: e.g. SPM(a,c) = sum_{b} TPM(a,b,c,b)
       * @param MT input TPM or PHM
       */
      template<class MatrixType> void bar(MatrixType &MT){

            for(int a = 0;a < M;++a)
               for(int b = a;b < M;++b){

                  (*this)(a,b) = 0.0;

                  for(int l = 0;l < M;++l)
                     (*this)(a,b) += MT(a,l,b,l);

                  (*this)(b,a) = (*this)(a,b);
               }

         }

   private:

      //!dimension of the single particle space, and dimension of the Matrix
      int M;

      //!number of particles
      int N;

};

// template specialization of bar for PPHM matrices.
template<> void SPM::bar(PPHM &MT);
template<> void SPM::bar(T2PM &MT);

#endif

/* vim: set ts=3 sw=3 expandtab :*/
