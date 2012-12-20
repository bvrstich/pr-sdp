#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

class SPM;
class SUP;
class PHM;
class DPM;
class PPHM;

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
      TPM(int,int);

      //copy constructor
      TPM(const TPM &);

      //file constructor
      TPM(const char *);

      //destructor
      virtual ~TPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      void hubbard(int option,double U);

      //Q afbeelding en zijn inverse
      void Q(int option,const TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,const TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,const TPM &);

      void unit();

      void proj_Tr();

      //de G down en inverse G up
      void G(int option,const PHM &);

      //trace one pair of indices of DPM
      void bar(const DPM &);

      //trace one pair of indices of PPHM
      void bar(const PPHM &);

      //T1 down
      void T(int option,const DPM &);

      //T2 down
      void T(const PPHM &);

      void min_unit(double scale);

      void min_qunit(double scale);

      void collaps(int option,const SUP &);

      void out(const char *filename);

      void sp_pairing(double );

      void in_sp(const char *);

      void constr_grad(double t,const TPM &,const SUP &);

      int solve(double t,const SUP &,TPM &);

      double line_search(double t,SUP &P,const TPM &ham);

      double line_search(double t,const TPM &,const TPM &);

      void H(double t,const TPM &b,const SUP &P);

      double S_2() const;

      void set_S_2();

      static int gn();

      static int gs2t(int,int);

      static int gt2s(int,int);

      static void init(int,int);

      static void clear();

      void constr_sp_diag(int);

   private:

      //!static list of dimension [n_tp][2] that takes in a tp index i and returns two sp indices: a = t2s[i][0] and b = t2s[i][1]
      static int **t2s;

      //!static list of dimension [M][M] that takes two sp indices a,b and returns a tp index i: i = s2t[a][b]
      static int **s2t;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!dimension of tp hilbert space and of the matrix
      static int n;

};

#endif
