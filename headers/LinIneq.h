#ifndef LININEQ_H
#define LININEQ_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;
using std::vector;

#include "LinCon.h"

/**
 * @author Brecht Verstichel
 * @date 08-12-2010\n\n
 * This is a class that contains the information about all the constraints active at a certain point in the program.
 */

class LinIneq{

   /**
    * output stream operator overloaded, will print all of the LinCon objects
    * @param output The stream to which you are writing (e.g. cout)
    * @param LinIneq_p de LinIneq object you want to print
    */
   friend ostream &operator<<(ostream &output,const LinIneq &li_p);

   public:

      //constructor
      LinIneq(int,int);

      //copy constructor
      LinIneq(const LinIneq &);

      //destructor
      virtual ~LinIneq();

      void add(const LinCon &);

      void rm(int);

      int gnr();

   private:

      //!std::vector object containing the different LinCon objects
      vector<LinCon> li;

      //!nr of sp orbitals
      int M;

      //!nr of particles
      int N;

};

#endif
