#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Vector.h"
#include "SUP.h"

//definitions:
#ifdef PQG

#define __G_CON

#endif

#ifdef PQGT1

#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT

#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

/**
 * @author Brecht Verstichel
 * @date 09-03-2010\n\n
 * This class, EIG is a "block"-vector over the carrierspace's of the active condtions. It contains room
 * to store the eigenvalues and special member function that work with these eigenvalues.
 * This class should only be used when a SUP matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 */
class EIG{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << eig_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << eig_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,const EIG &eig_p);

   public:

      //default constructor
      EIG();
   
      //constructor met initialisatie op 
      EIG(SUP &);
      
      //copy constructor
      EIG(const EIG &);

      //destructor
      ~EIG();

      void diagonalize(SUP &);

      int gN() const;

      int gM() const;

      int gn_tp() const;

      int gdim() const;

      double centerpot(double,const EIG &,double,double) const;

      //overload equality operator
      EIG &operator=(const EIG &);

      Vector<TPM> &tpv(int);

      const Vector<TPM> &tpv(int) const;

#ifdef __G_CON

      int gn_ph() const;

      Vector<PHM> &phv();

      const Vector<PHM> &phv() const;

#endif

#ifdef __T1_CON

      int gn_dp() const;

      Vector<DPM> &dpv();

      const Vector<DPM> &dpv() const;

#endif

#ifdef __T2_CON

      int gn_pph() const;

      Vector<PPHM> &pphv();

      const Vector<PPHM> &pphv() const;

#endif

      double min() const;

      double max() const;

      double center_dev() const;

      double lsfunc(double) const;

   private:

      //!variable that tells if the memory has been allocated (flag = 1) or not (flag = 0)
      int flag;

      //!double pointer to a Vector<TPM> object, the eigenvalues of the P and Q part of a SUP matrix will be stored here.
      Vector<TPM> **v_tp;

      //!number of particles
      int N;

      //!dimension of sp space
      int M;

      //!dimension of tp space
      int n_tp;

#ifdef __G_CON

      //!pointer to a Vector<PHM> object that will contain the eigenvalues of the G part of a SUP matrix
      Vector<PHM> *v_ph;
      
      //!dimension of ph space
      int n_ph;

#endif

#ifdef __T1_CON

      //!pointer to a Vector<DPM> object that will contain the eigenvalues of the T1 part of a SUP matrix
      Vector<DPM> *v_dp;
      
      int n_dp;

#endif

#ifdef __T2_CON

      //!pointer to a Vector<PPHM> object that will contain the eigenvalues of the T1 part of a SUP matrix
      Vector<PPHM> *v_pph;

      int n_pph;

#endif

      //!total dimension of the EIG object
      int dim;

};

#endif
