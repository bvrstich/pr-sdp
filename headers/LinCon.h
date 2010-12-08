#ifndef LINCON_H
#define LINCON_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class TPM;

/**
 * @author Brecht Verstichel
 * @date 08-12-2010\n\n
 * This is a class that contains the information about a single linear constraint.
 */

class LinCon{

   /**
    * output stream operator overloaded, will print the constraint matrix, the value of the minimal projection, and the current projection.
    * @param output The stream to which you are writing (e.g. cout)
    * @param LinCon_p de LinCon object you want to print
    */
   friend ostream &operator<<(ostream &output,LinCon &lc_p);

   public:

      //constructor
      LinCon(TPM &,double);

      //copy constructor
      LinCon(LinCon &);

      //destructor
      virtual ~LinCon();

      TPM &gI();

      double gi();

      void init(TPM &);

      double gtpm_I();

   private:

      //!Constraint matrix
      TPM *I_c;

      //!minimal projection on the constraint matrix, such that Tr(\Gamma I_C) \geq i_c
      double i_c;

      //!projection of the tpm object with wich the LinCon has been initialized on the constraint matrix.
      double tpm_I;

};

#endif
