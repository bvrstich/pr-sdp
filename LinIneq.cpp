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
 * constructor of a LinIneq object
 * @param M nr of sp orbitals
 * @param N nr of particles
 * @param nr nr of constraints
 */
LinIneq::LinIneq(int M,int N,int nr){

   this->M = M;
   this->N = N;

   this->nr = nr;

   li = new LinCon * [nr];

   for(int i = 0;i < nr;++i)
      li[i] = new LinCon(M,N);

}

/**
 * copy constructor
 * @param lc_copy The LinIneq object to be copied
 */
LinIneq::LinIneq(const LinIneq &lc_copy){

}

/**
 * destructor
 */
LinIneq::~LinIneq(){

   for(int i = 0;i < nr;++i)
      delete li[i];

   delete [] li;

}

/**
 * @return the number of currently applied linear constraints
 */
int LinIneq::gnr(){

   return nr;

}

/**
 * read and write access to your LinCon object
 * @param i row number
 * @return the entry on index i
 */
LinCon &LinIneq::operator[](int i){

   return *li[i];

}
