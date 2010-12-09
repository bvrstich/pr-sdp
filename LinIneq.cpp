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
LinIneq::LinIneq(int M,int N){

   this->M = M;
   this->N = N;

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

}

/**
 * @return the number of currently applied linear constraints
 */
int LinIneq::gnr(){

   return li.size();

}

/**
 * Add a LinCon object to the vector, so add a constraint to the program.
 * @param lc The LinCon object to be added
 */
void LinIneq::add(const LinCon &lc){

   li.push_back(lc);

}

/** 
 * remove a LinCon object from the array
 * @param index The index at which the object to be removed is located.
 */
void LinIneq::rm(int index){

   li.erase(li.begin() + index - 1);

}
