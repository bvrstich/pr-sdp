#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;
using std::string;

#include "Matrix.h"

class SPM;
class SUP;
class PHM;
class DPM;

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class TPM is a class written for two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */
class TPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPM &tpm_p);

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

      //geef M terug
      int gM();

      //geef n terug
      int gn();

      void hubbard(double U);

      void Q(TPM &);

#ifdef __G_CON

      void G(PHM &);

#endif

#ifdef __T1_CON

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

      //!static list of dimension [n_tp][2] that takes in a tp index i and returns two sp indices: a = t2s[i][0] and b = t2s[i][1]
      static int **t2s;

      //!static list of dimension [M][M] that takes two sp indices a,b and returns a tp index i: i = s2t[a][b]
      static int **s2t;

      //!static counter that counts the number of TPM objects running in the program
      static int counter;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

      //!dimension of tp hilbert space and of the matrix
      int n;

};

#endif
