#ifndef EIG_H
#define EIG_H

#ifdef PQ

#include "EIG_PQ.h"

class EIG : public EIG_PQ { 

   public :

      EIG(int M,int N) : EIG_PQ(M,N) { }

      EIG(EIG &eig) : EIG_PQ(eig) { }

      EIG(SUP &SZ) : EIG_PQ(SZ) { }

      ~EIG(){ }
   
};

#endif

#ifdef PQG

#include "EIG_PQG.h"

class EIG : public EIG_PQG { 

   public :

      EIG(int M,int N) : EIG_PQG(M,N) { }

      EIG(EIG &eig) : EIG_PQG(eig) { }

      EIG(SUP &SZ) : EIG_PQG(SZ) { }

      ~EIG(){ }
   
};

#endif

#endif
