#ifndef SUP_PQGT2_H
#define SUP_PQGT2_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "../TPM.h"
#include "../DPM.h"
#include "../PPHM.h"
#include "SUP_PQG.h"

class EIG;

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, SUP_PQGT2, is a blockmatrix over the carrierspace's of the P, Q, G and T1 conditions. It inherits form its direct motherclass
 * SUP_PQG and its mother SUP_PQ . This class will expand its mother by a pointer to a DPM object (which is independent of the other objects, by which I mean that
 * DPM::T( SUP_PQGT2::tpm (0) ) is not neccesarily equal to SUP_PQGT2::dpm() ), and it will also redefine its mothers functions to include the DPM contribution.
 */
class SUP_PQGT2 : public SUP_PQG {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param sup_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,SUP_PQGT2 &sup_p);

   public:

      //constructor
      SUP_PQGT2(int M,int N);

      //copy constructor
      SUP_PQGT2(SUP_PQGT2 &);

      //destructor
      ~SUP_PQGT2();

      //overload += operator
      SUP_PQGT2 &operator+=(SUP_PQGT2 &);

      //overload -= operator
      SUP_PQGT2 &operator-=(SUP_PQGT2 &);

      //overload equality operator
      SUP_PQGT2 &operator=(SUP_PQGT2 &);

      //overload equality operator
      SUP_PQGT2 &operator=(double &);

      SUP_PQGT2 operator*(SUP_PQGT2 &);

      int gn_pph();

      PPHM &pphm();

      double ddot(SUP_PQGT2 &);

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP_PQGT2 &,SUP_PQGT2 &);

      void daxpy(double alpha,SUP_PQGT2 &);

      double trace();

      void fill(TPM &);

   protected:

      //!pointer to a DPM object, will containt the T1 carrier space contribution of SUP
      PPHM *SZ_pph;

      //!dimension of three particle space
      int n_pph;

};

#endif
