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
