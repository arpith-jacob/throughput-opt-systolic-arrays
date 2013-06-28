
//  throughput-ilp.hpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 7 2009
//
//  Create ILP to find k_max, the maximum number of points projected onto any
//  processing element by projection vector u.
//
//  Input is an input polyhedron Ax <= b in PIP matrix format
//  Output is the ILP in PIP matrix format

#ifndef __THROUGHPUT_ILP_H__
#   define __THROUGHPUT_ILP_H__

#include <cstdio>

// pip includes
#include <piplib/piplibMP.h>

class ThroughputILP
{

 public:

   // constructor
   ThroughputILP(PipMatrix *polyhedron, PipMatrix *context,
                 unsigned int dimensions, unsigned int parameters,
                 ublas::vector<int> *pv)
   {
//     pip_matrix_print(stdout, polyhedron);
//     pip_matrix_print(stdout, context);

     // generate throughput ilp
     GenThroughputILP(polyhedron, context, dimensions, parameters, pv);

//     cout << endl << "Projection vector: " << *pv << endl;
//     pip_matrix_print(stdout, throughputilp);
//     pip_matrix_print(stdout, contextilp);
     
     return;
   }

   // destructor
   ~ThroughputILP()
   {
     pip_matrix_free(contextilp);
     pip_matrix_free(throughputilp);
   }
   
   // get throughput ilp
   PipMatrix *getILP()
   {
     return throughputilp;
   }
   
   // get context
   PipMatrix *getContext()
   {
     return contextilp;
   }
   
   // get big parameter position
   int getBigParamPos()
   {
     return bigParamPos;
   }
   
private:

   PipMatrix *throughputilp, *contextilp;
   int bigParamPos;

   void GenThroughputILP(PipMatrix *polyhedron, PipMatrix *context,
                         unsigned int dimensions, unsigned int parameters,
                         ublas::vector<int> *pv)
   {
     int no_constraints = polyhedron->NbRows;

     // allocate memory for throughput ilp constraints
     //  Number of constraints = #orig_constraints * 2 + dimensions
     //        A x1 <= b; A x2 <= b; x1 - x2 = ku
     //
     //  Number of columns = 1 + dimensions * 2 + parameters + 3 (const, B, k')
     throughputilp = pip_matrix_alloc(
                       no_constraints * 2 + dimensions,
                       1 + dimensions * 2 + parameters + 3
                     );

     // format of columns is as follows
     //       equality?  k'  i1  ... k1 i2 ... k2 N1 ... N3 B const 

     // copy polyhedron constraints: A x1 <= b
     for (int i = 0; i < no_constraints; i++) {
       // equality/inequality?
       entier_assign (throughputilp->p[i][0], polyhedron->p[i][0]);

       // set k' column = 0
       entier_set_si (throughputilp->p[i][1], 0);

       // copy dimensions (unknowns) to A x1 <= b
       for (unsigned int j = 1; j <= dimensions; j++) {
         entier_assign (
                        throughputilp->p[i][j+1],
                        polyhedron->p[i][j]
                        );
       }

       // set dimensions (unknowns) for A x2 <= b to zero
       for (unsigned int j = 1; j <= dimensions; j++) {
         entier_set_si (throughputilp->p[i][j+dimensions+1], 0);
       }

       // set parameters
       for (unsigned int j = 1; j <= parameters; j++) {
         entier_assign (
                        throughputilp->p[i][j+2*dimensions+1],
                        polyhedron->p[i][j+dimensions]
                        );
       }

       // set B, const column
       entier_set_si (throughputilp->p[i][dimensions*2+parameters+2], 0);
       entier_assign (
                      throughputilp->p[i][dimensions*2+parameters+3],
                      polyhedron->p[i][dimensions+parameters+1]
                      );
     }

     // duplicate polyhedron constraints: A x2 <= b
     for (int i = no_constraints; i < no_constraints*2; i++) {
       // equality/inequality?
       entier_assign (
                      throughputilp->p[i][0],
                      polyhedron->p[i-no_constraints][0]
                      );

       // set k' column = 0
       entier_set_si (throughputilp->p[i][1], 0);

       // copy dimensions (unknowns) to A x1 <= b to zero
       for (unsigned int j = 1; j <= dimensions; j++) {
         entier_set_si (throughputilp->p[i][j+1], 0);
       }

       // set dimensions (unknowns) for A x2 <= b
       for (unsigned int j = 1; j <= dimensions; j++) {
         entier_assign (
                        throughputilp->p[i][j+dimensions+1],
                        polyhedron->p[i-no_constraints][j]
                        );
       }

       // set parameters
       for (unsigned int j = 1; j <= parameters; j++) {
         entier_assign (
                        throughputilp->p[i][j+2*dimensions+1],
                        polyhedron->p[i-no_constraints][j+dimensions]
                        );
       }

       // set B, const column
       entier_set_si (throughputilp->p[i][dimensions*2+parameters+2], 0);
       entier_assign (
                      throughputilp->p[i][dimensions*2+parameters+3],
                      polyhedron->p[i-no_constraints][dimensions+parameters+1]
                      );
     }

     // now specify the constraints: x1 - x2 = ku
     // to maximize k, we use k' = B - k where B is a big parameter and minimize
     // k'.  We specify k' as the first unknown to minimize it
     //
     for (unsigned int i = 0; i < dimensions; i++) {
        // equality
        entier_set_si (throughputilp->p[i+no_constraints*2][0], 0);

        // set k' column = projection_vector[]
        entier_set_si (
                       throughputilp->p[i+no_constraints*2][1],
                       (*pv)(i)
                       );

        // zero dimensions (unknowns) for A x1 <= b and A x2 <= b to zero
        for (unsigned int j = 1; j <= 2*dimensions; j++) {
          entier_set_si (throughputilp->p[i+no_constraints*2][j+1], 0);
        }

        // set x1 and -x2 for this dimension (i/j/k)
        entier_set_si (throughputilp->p[i+no_constraints*2][2+i], 1);
        entier_set_si (
                       throughputilp->p[i+no_constraints*2][2+dimensions+i],
                       -1
                       );

        // set parameters to zero
        for (unsigned int j = 1; j <= parameters; j++) {
          entier_set_si (
                         throughputilp->p[i+no_constraints*2][j+2*dimensions+1],
                         0
                         );
        }

        // set B, const column
        entier_set_si (
              throughputilp->p[i+no_constraints*2][dimensions*2+parameters+2],
              -1 * (*pv)(i)
              );

        entier_set_si (
              throughputilp->p[i+no_constraints*2][dimensions*2+parameters+3],
              0
              );
     }
    
     // column position of big parameter in constraint row
     // first column (equality/inequality?) starts at index 0
     bigParamPos = dimensions*2 + parameters + 2;

     //
     //  generate context for parameters
     //
     // allocate memory for parameter context constraints
     //  Number of constraints = #orig_constraints
     //
     //  Number of columns = original_columns + 1 (for new parameter B)
     contextilp = pip_matrix_alloc(
                       context->NbRows,
                       context->NbColumns + 1
                     );

     // format of columns is as follows
     //       equality?  N1 ... N3 B const 
     for (unsigned int i = 0; i < context->NbRows; i++) {
       for (unsigned int j = 0; j < context->NbColumns - 1; j++) {
         entier_assign (
                        contextilp->p[i][j],
                        context->p[i][j]
                        );
       }

       // B parameter
       entier_set_si (contextilp->p[i][context->NbColumns-1], 0);
       // const
       entier_assign (
                      contextilp->p[i][context->NbColumns],
                      context->p[i][context->NbColumns-1]
                      );
     }
   }

};

#endif // __THROUGHPUT_ILP_H__
