#ifndef LININEQ_H
#define LININEQ_H

#include <iostream>
#include <cstdlib>

using std::ostream;

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
    * @param li_p the LinIneq object you want to print
    */
   friend ostream &operator<<(ostream &output,const LinIneq &li_p);

   public:

      //constructor
      LinIneq(int,int,int);

      //copy constructor
      LinIneq(const LinIneq &);

      //destructor
      virtual ~LinIneq();

      int gnr();

      //easy to access and change the LinCon objects
      LinCon &operator[](int i);

   private:

      //!LinCon array containing the different LinCon objects
      LinCon **li;

      //!nr of linear constraints
      int nr;

      //!nr of sp orbitals
      int M;

      //!nr of particles
      int N;

};

#endif
