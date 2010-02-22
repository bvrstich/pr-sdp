#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "../headers/EIG/EIG_PQG.h"
#include "../headers/SUP/SUP_PQG.h"
#include "../headers/lapack.h"

//constructor:
EIG_PQG::EIG_PQG(int M,int N) : EIG_PQ(M,N){
   
   this->n_ph = M*M;

   dim += n_ph;

   eig_ph = new double [n_ph];

}

//copy constructor:
EIG_PQG::EIG_PQG(EIG_PQG &eig_c) : EIG_PQ(eig_c){

   this->n_ph = eig_c.n_ph;

   dim += n_ph;

   eig_ph = new double [n_ph];

   int inc = 1;

   dcopy_(&n_ph,eig_c.eig_ph,&inc,eig_ph,&inc);

}

//constructor met initialisatie door SUP matrix:
EIG_PQG::EIG_PQG(SUP_PQG &SZ) : EIG_PQ(SZ){

   this->n_ph = SZ.gn_ph();
   dim += n_ph;

   eig_ph = new double [n_ph];

   SZ.phm().diagonalize(eig_ph);

}

EIG_PQG &EIG_PQG::operator=(EIG_PQG &eig_c){

   this->EIG_PQ::operator=(eig_c);

   int inc = 1;

   dcopy_(&n_ph,eig_c.eig_ph,&inc,eig_ph,&inc);

   return *this;

}

double *EIG_PQG::operator[](int i){

   if(i < 2)
      return this->EIG_PQ::operator[](i);
   else
      return eig_ph;

}

int EIG_PQG::gn_ph(){

   return n_ph;

}
 
//destructor
EIG_PQG::~EIG_PQG(){

   delete [] eig_ph;

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,EIG_PQG &eig_p){

   for(int i = 0;i < eig_p.gn_tp();++i)
      output << i << "\t" << eig_p[0][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_tp();++i)
      output << i << "\t" << eig_p[1][i] << std::endl;

   output << std::endl;

   for(int i = 0;i < eig_p.gn_ph();++i)
      output << i << "\t" << eig_p[2][i] << std::endl;

   return output;

}

double EIG_PQG::operator()(int block,int index){

   if(block < 2)
      return this->EIG_PQ::operator()(block,index);
   else
      return eig_ph[index];

}

double EIG_PQG::lsfunc(double alpha){

   double ward = this->EIG_PQ::lsfunc(alpha);

   for(int i = 0;i < n_ph;++i)
      ward += eig_ph[i]/(1.0 + alpha*eig_ph[i]);

   return ward;

}

double EIG_PQG::min(){

   double ward = this->EIG_PQ::min();

   if(ward > eig_ph[0])
      ward = eig_ph[0];

   return ward;

}
