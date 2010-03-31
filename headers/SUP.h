/**
 * @file 
 * This is a wrapper class around the different SUP_PQ(GT1) classes. It is decided at compile time from which 
 * SUP_* file this class inherits. Compile with PQ to inherit from SUP_PQ, compile with PQG to inherit from SUP_PQG, etc. .
 * This way, you can use the SUP object everywhere in the program without having to worry about which conditions are used.
 */
#ifndef SUP_H
#define SUP_H

//if PQ is defined, inherit from SUP_PQ
#ifdef PQ

#include "SUP/SUP_PQ.h"

class SUP : public SUP_PQ{

   public :

      SUP(int M,int N) : SUP_PQ(M,N) { }

      SUP(SUP &SZ) : SUP_PQ(SZ) { }

      ~SUP(){ }

};

#endif

//if PQG is defined, inherit from SUP_PQG
#ifdef PQG

#include "SUP/SUP_PQG.h"

class SUP : public SUP_PQG{

   public :

      SUP(int M,int N) : SUP_PQG(M,N) { }

      SUP(SUP &SZ) : SUP_PQG(SZ) { }

      ~SUP(){ }

};

#endif

//if PQGT1 is defined, inherit from SUP_PQGT1
#ifdef PQGT1

#include "SUP/SUP_PQGT1.h"

class SUP : public SUP_PQGT1{

   public :

      SUP(int M,int N) : SUP_PQGT1(M,N) { }

      SUP(SUP &SZ) : SUP_PQGT1(SZ) { }

      ~SUP(){ }

};

#endif


//if PQGT1 is defined, inherit from SUP_PQGT1
#ifdef PQGT2

#include "SUP/SUP_PQGT2.h"

class SUP : public SUP_PQGT2{

   public :

      SUP(int M,int N) : SUP_PQGT2(M,N) { }

      SUP(SUP &SZ) : SUP_PQGT2(SZ) { }

      ~SUP(){ }

};

#endif

#endif
