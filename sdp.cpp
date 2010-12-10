/**
 * @mainpage 
 * This is an implementation of the dual only, potential reduction interior point method
 * for optimizing the second order density matrix using the P, Q, G , T_1, T_2 and T_2' N-representability conditions.
 * This has an extra class for adding linear inequality constraints.
 * Compiling can be done with the options PQ, PQG, PQGT1, PQGT2, PQGT2P and PQGT (for all conditions active) with logical consequences for the program.
 * @author Brecht Verstichel, Ward Poelmans
 * @date 08-12-2010
 */

#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity density matrix normed on the particle number and minimize the 
 * ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(void){

   cout.precision(10);

   const int M = 8;//dim sp hilbert space
   const int N = 4;//nr of particles

   //hamiltoniaan
   TPM ham(M,N);

   //the zero is for pbc's
   ham.hubbard(0,1.0);

   TPM rdm(M,N);
   rdm.unit();

   TPM backup_rdm(rdm);

   double t = 1.0;
   double tolerance = 1.0e-5;

   LinIneq li(M,N,1);

   li[0].sI(ham);
   li[0].si(-3.0);

   //outer iteration: scaling of the potential barrier
   while(t > 1.0e-12){

      cout << t << "\t" << rdm.trace() << "\t" << rdm.ddot(ham) << "\t";

      int nr_cg_iter = 0;
      int nr_newton_iter = 0;

      double convergence = 1.0;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      while(convergence > tolerance){

         ++nr_newton_iter;

         SUP P(M,N);

         P.fill(rdm);

         li.fill(rdm);

         P.invert();

         //eerst -gradient aanmaken:
         TPM grad(M,N);

         grad.constr_grad(t,ham,P,li);

         //dit wordt de stap:
         TPM delta(M,N);

         //los het hessiaan stelsel op:
         nr_cg_iter += delta.solve(t,P,grad,li);

         //line search
         double a = delta.line_search(t,P,ham,li);

         //rdm += a*delta;
         rdm.daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);

      }

      cout << nr_newton_iter << "\t" << nr_cg_iter << endl;

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

      double a = extrapol.line_search(t,rdm,ham,li);

      rdm.daxpy(a,extrapol);

   }

   cout << endl;
   
   cout << "Final Energy:\t" << ham.ddot(rdm) << endl;
   cout << endl;
   cout << "Final Spin:\t" << rdm.S_2() << endl;

   return 0;
}
