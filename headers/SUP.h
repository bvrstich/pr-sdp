#ifndef SUP_H
#define SUP_H

#include "SUP/SUP_PQ.h"

#ifdef PQ

class SUP : public SUP_PQ{

   public :

      SUP(int M,int N) : SUP_PQ(M,N) { }

      SUP(SUP &SZ) : SUP_PQ(SZ) { }

      ~SUP(){ }

};

#endif

#ifdef PQG

#include "SUP/SUP_PQG.h"

class SUP : public SUP_PQG{

   public :

      SUP(int M,int N) : SUP_PQG(M,N) { }

      SUP(SUP &SZ) : SUP_PQG(SZ) { }

      ~SUP(){ }

};

#endif

#endif
