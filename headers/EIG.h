#ifndef EIG_H
#define EIG_H

#ifdef PQ

#include "EIG/EIG_PQ.h"

class EIG : public EIG_PQ { 

   public :

      EIG(int M,int N) : EIG_PQ(M,N) { }

      EIG(EIG &eig) : EIG_PQ(eig) { }

      EIG(SUP_PQ &SZ) : EIG_PQ(SZ) { }

      ~EIG(){ }
   
};

#endif

#ifdef PQG

#include "EIG/EIG_PQG.h"

class EIG : public EIG_PQG { 

   public :

      EIG(int M,int N) : EIG_PQG(M,N) { }

      EIG(EIG &eig) : EIG_PQG(eig) { }

      EIG(SUP_PQG &SZ) : EIG_PQG(SZ) { }

      ~EIG(){ }
   
};

#endif

#ifdef PQGT1

#include "EIG/EIG_PQGT1.h"

class EIG : public EIG_PQGT1 { 

   public :

      EIG(int M,int N) : EIG_PQGT1(M,N) { }

      EIG(EIG &eig) : EIG_PQGT1(eig) { }

      EIG(SUP_PQGT1 &SZ) : EIG_PQGT1(SZ) { }

      ~EIG(){ }
   
};

#endif

#endif
