#ifndef SUP_PQG_H
#define SUP_PQG_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "../TPM.h"
#include "../PHM.h"
#include "SUP_PQ.h"

class EIG;

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, SUP_PQG is a blockmatrix over the carrierspace's of the P, Q and G conditions. It inherits form its motherclass
 * SUP_PQ. This class will expand its mother by a pointer to a PHM object (which is independent of the TPM objects, by which I mean that
 * PHM::G( SUP_PQG::tpm (0)) is not neccesarily equal to SUP_PQG::phm()), and it will also redefine its mothers functions to include the PHM contribution.
 */
class SUP_PQG : public SUP_PQ {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param sup_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,SUP_PQG &sup_p);
  
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

      //!Pointer to the PHM object, will containt the G-space matrix of SUP
      PHM *SZ_ph;

      //!dimension of ph space
      int n_ph;

};

#endif
