#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

#include "headers/lapack.h"
#include "headers/Matrix.h"
#include "headers/TPM.h"
#include "headers/SPM.h"

#ifdef PQG

#include "headers/PHM.h"

#endif

#ifdef PQGT1

#include "headers/DPM.h"

#endif

#include "headers/SUP.h"
#include "headers/EIG.h"

int main(void){

   cout.precision(10);

   int M = 8;//dim sp hilbert space
   int N = 4;//nr of particles

   //hamiltoniaan
   TPM ham(M,N);

   ham.hubbard(1.0);

   double norm_ham = sqrt(ham.ddot(ham));

   ham /= norm_ham;

   TPM rdm(M,N);
   rdm.init();

   TPM backup_rdm(rdm);

   double t = 1.0;
   double tolerance = 1.0e-5;

   //outer iteration: scaling of the potential barrier
   while(t > 1.0e-12){

      cout << t << "\t" << rdm.trace() << "\t" << rdm.ddot(ham)*norm_ham << endl;

      double convergence = 1.0;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      while(convergence > tolerance){

         SUP P(M,N);

         P.fill(rdm);

         P.invert();

         //eerst -gradient aanmaken:
         TPM grad(M,N);

         grad.constr_grad(t,ham,P);

         //dit wordt de stap:
         TPM delta(M,N);

         //los het hessiaan stelsel op:
         cout << delta.solve(t,P,grad) << endl;

         //line search
         double a = delta.line_search(t,P,ham);

         //rdm += a*delta;
         rdm.daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);

      }

      t /= 2.0;

      //what is the tolerance for the newton method?
      tolerance = 1.0e-5*t;

      if(tolerance < 1.0e-12)
         tolerance = 1.0e-12;

      //extrapolatie:
      TPM extrapol(rdm);

      extrapol -= backup_rdm;

      //overzetten voor volgende stap
      backup_rdm = rdm;

      double a = extrapol.line_search(t,rdm,ham);

      rdm.daxpy(a,extrapol);

   }

   SUP P(M,N);
   P.fill(rdm);

   EIG eig(P);

   cout << eig;

   return 0;

}
