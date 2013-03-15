#ifndef cPHM_H
#define cPHM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, cPHM, is a class written for particle-hole matrices, it inherits all the functions from its mother class
 * Matrix, some special member functions and two lists that give the relationship between the sp and the ph basis.
 */
class cPHM : public Matrix {

    /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << cphm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << cphm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param cphm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const cPHM &cphm_p);

   public:
      
      //constructor
      cPHM(int D);

      //copy constructor
      cPHM(const cPHM &);

      //destructor
      virtual ~cPHM();

      using Matrix::operator=;

      using Matrix::operator();

      //access the numbers in sp mode
      double operator()(int a,int b,int c,int d) const;

      int gD() const;

      void fill(const PHM &);

   private:

      //!static counter that counts the number of cPHM objects running in the program
      static int counter;

      //!static list of dimension [n_ph][2] that takes in a ph index i and returns two sp indices: a = ph2s[i][0] and b = ph2s[i][1]
      static int **ph2s;

      //!static list of dimension [M][M] that takes two sp indices a,b and returns a ph index i: i = s2t[a][b]
      static int **s2ph;

      int D;

};

#endif
