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
 * constructor of a LinCon object
 * @param I The constraint matrix
 * @param i the minimal projection
 */
LinCon::LinCon(TPM &I,double i){

   I_c = new TPM(I);

   i_c = i;

}

/**
 * copy constructor
 * @param lc_copy The LinCon object to be copied
 */
LinCon::LinCon(LinCon &lc_copy){

   I_c = new TPM(lc_copy.gI());

   i_c = lc_copy.gi();

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
TPM &LinCon::gI(){

   return *I_c;

}

/**
 * @return The minimal projection
 */
double LinCon::gi(){

   return i_c;

}

void LinCon::init(TPM &tpm){

   tpm_I = tpm.ddot(*I_c);

}

/**
 * @return the projection of the input (LinCon::init) TPM on the constraint matrix
 */
double LinCon::gtpm_I(){

   return tpm_I;

}

ostream &operator<<(ostream &output,LinCon &lc_p){

   cout << "The current projection and the minimal projection" << endl;
   cout << lc_p.gtpm_I() << "\t" << lc_p.gi() << endl;
   cout << endl;

   cout << "The Constraint matrix:" << endl;
   cout << endl;

   cout << lc_p.gI() << endl;

   return output;

}