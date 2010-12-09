#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * Constructor of a LinCon object
 * @param M The constraint matrix
 * @param N the minimal projection
 */
LinCon::LinCon(int M,int N){

   I_c = new TPM(M,N);

   this->M = M;
   this->N = N;

}

/**
 * copy constructor
 * @param lc_copy The LinCon object to be copied
 */
LinCon::LinCon(const LinCon &lc_copy){

   I_c = new TPM(lc_copy.gI());

   i_c = lc_copy.gi();

   tpm_I = lc_copy.gtpm_I();

   this->M = lc_copy.gM();
   this->N = lc_copy.gN();

}

/**
 * destructor
 */
LinCon::~LinCon(){

   delete I_c;

}

/**
 * @return the constraint TPM object
 */
TPM &LinCon::gI() const{

   return *I_c;

}

/**
 * @return The minimal projection
 */
double LinCon::gi() const{

   return i_c;

}

/**
 * set the constraint value
 * @param i the value that the minimal projection will be set to.
 */
void LinCon::si(double i){

   i_c = i;

}

/**
 * set the constraint Matrix
 * @param I the input constraint Matrix
 */
void LinCon::sI(const TPM &I){

   *I_c = I;

}

/**
 * Initialize the tpm_I value to a TPM object
 * @param tpm input TPM object
 */
void LinCon::stpm_I(const TPM &tpm){

   tpm_I = tpm.ddot(*I_c);

}

/**
 * @return the projection of the input (LinCon::init) TPM on the constraint matrix
 */
double LinCon::gtpm_I() const{

   return tpm_I;

}

ostream &operator<<(ostream &output,const LinCon &lc_p){

   cout << "The current projection and the minimal projection" << endl;
   cout << lc_p.gtpm_I() << "\t" << lc_p.gi() << endl;
   cout << endl;

   cout << "The Constraint matrix:" << endl;
   cout << endl;

   cout << lc_p.gI() << endl;

   return output;

}

/**
 * @return nr of sp orbitals
 */
int LinCon::gM() const{

   return M;

}

/**
 * @return nr of particles
 */
int LinCon::gN() const{

   return N;

}
