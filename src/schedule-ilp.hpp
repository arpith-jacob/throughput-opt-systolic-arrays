
//  schedule-ilp.hpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 10 2009
//
//  Create ILP to find schedule for a given projection vector.
//  The schedule is constrained to respect dependencies.
//  Objective function: minimize array utilization (\lambda u) and latency
//
//  Input is the projection vector (projection solution object), dependencies
//  and vertices in PIP matrix format
//  Output is the ILP in PIP matrix format

#ifndef __SCHEDULE_ILP_H__
#   define __SCHEDULE_ILP_H__

#include <cstdio>

// pip includes
#include <piplib/piplibMP.h>

class ScheduleILP
{

 public:

   // constructor
   ScheduleILP(unsigned int dimensions, unsigned int parameters,
               PipMatrix *dependencies, PipMatrix *vertices,
               unsigned int pepipelinestages,
               ProjectionSolution *ps)
   {
//     pip_matrix_print(stdout, dependencies);
//     pip_matrix_print(stdout, vertices);

     // generate schedule ilp
     GenScheduleILP(dimensions, parameters, dependencies, vertices,
                    pepipelinestages, ps);

//     pip_matrix_print(stdout, scheduleilp);

     return;
   }

   // destructor
   ~ScheduleILP()
   {
     pip_matrix_free(scheduleilp);
     pip_matrix_free(contextilp);
   }
   
   // regenerate ilp.  called after negating projection vector
   void regenILP(unsigned int dimensions, unsigned int parameters,
                 PipMatrix *dependencies, PipMatrix *vertices,
                 unsigned int pepipelinestages,
                 ProjectionSolution *ps)
   {
     GenScheduleILP(dimensions, parameters, dependencies, vertices,
                    pepipelinestages, ps);
   }

   // get schedule ilp
   PipMatrix *getILP()
   {
     return scheduleilp;
   }

   // get schedule context
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

   PipMatrix *scheduleilp, *contextilp;
   int bigParamPos;

   void GenScheduleILP(unsigned int dimensions, unsigned int parameters,
                       PipMatrix *dependencies, PipMatrix *vertices,
                       unsigned int pepipelinestages,
                       ProjectionSolution *ps)
   {
     int no_dependencies = dependencies->NbRows;
     int no_vertices     = vertices->NbRows;

     // temporary (big num) for manipulation
     Entier tmp, tmp2;
     entier_init (tmp);
     entier_init (tmp2);


     // allocate memory for schedule ilp constraints
     //  Number of constraints = #dependencies + 
     //                          (#vertices * #vertices - #vertices) + 3
     //
     //  ld <= -1 ; t >= lu ; lu >= 1 ; lv_d <= s ; q >= 2t + s
     //
     //  Number of columns = 1 + dimensions + 5 (q, t, s, const, B)
     //    We are using big parameter B so that l can be negative
     scheduleilp = pip_matrix_alloc(
                       no_dependencies + 
                       (no_vertices * no_vertices - no_vertices) + 3,
                       1 + dimensions + 5
                     );

     // format of columns is as follows
     //       equality?  q  t  s  l1  ... ln  B  const
     // each of l1 ... ln is of format l1 = l1' - G 

     //
     // minimize t = array utilization
     // constraint: t >= lu    t - lu >= 0    t - l'u + uG >= 0
     //
     entier_set_si (scheduleilp->p[0][0], 1);  // inequality
     entier_set_si (scheduleilp->p[0][1], 0);  // q
     entier_set_si (scheduleilp->p[0][2], 1);  // t
     entier_set_si (scheduleilp->p[0][3], 0);  // s

     // l1 ... ln
     int pv_sum = 0;
     for (unsigned int i = 0; i < dimensions; i++) {
       entier_oppose (tmp, ps->projection_vector->p[0][i]);
       entier_assign (scheduleilp->p[0][4+i], tmp);
       pv_sum += VALUE_TO_INT( ps->projection_vector->p[0][i] );
     }

     entier_set_si (scheduleilp->p[0][4+dimensions], pv_sum);  // u1 + ... + un
     entier_set_si (scheduleilp->p[0][4+dimensions+1], 0);       // const

     //
     // lu != 0.  Here we try lu >= 1; if this doesn't yield a solution,
     // negate u and rerun ILP.  Only one ILP can produce a solution.
     // constraint: lu >= 1   lu - 1 >= 0   lu - 1 - uG >= 0
     //
     entier_set_si (scheduleilp->p[1][0], 1);  // inequality
     entier_set_si (scheduleilp->p[1][1], 0);  // q
     entier_set_si (scheduleilp->p[1][2], 0);  // t
     entier_set_si (scheduleilp->p[1][3], 0);  // s

     // l1 ... ln
     for (unsigned int i = 0; i < dimensions; i++) {
       entier_assign (scheduleilp->p[1][4+i], ps->projection_vector->p[0][i]);
     }

     entier_set_si (scheduleilp->p[1][4+dimensions], -pv_sum);  // u1 + ... + un
     entier_set_si (scheduleilp->p[1][4+dimensions+1], -1);     // const

     //
     // Objective function
     // constraint: q >= 2048t + s   q - 2048t - s >= 0
     //
     // This is the final objective function to minimize.
     //  t - array utilization.  modify weight as needed
     //  s - latency
     //
     entier_set_si (scheduleilp->p[2][0], 1);     // inequality
     entier_set_si (scheduleilp->p[2][1], 1);     // q
     entier_set_si (scheduleilp->p[2][2], -2048); // t
     entier_set_si (scheduleilp->p[2][3], -1);    // s

     // l1 ... ln
     for (unsigned int i = 0; i < dimensions; i++) {
       entier_set_si (scheduleilp->p[2][4+i], 0);
     }

     entier_set_si (scheduleilp->p[2][4+dimensions], 0);    // u1 + ... + un
     entier_set_si (scheduleilp->p[2][4+dimensions+1], 0);  // const

     //
     // Constraints for all dependencies
     //
     // constraint: ld <= -1    -ld -1 >= 0    -ld -1 + dB >= 0
     //
     for (int i = 0; i < no_dependencies; i++) {
       entier_set_si (scheduleilp->p[3+i][0], 1);   // inequality
       entier_set_si (scheduleilp->p[3+i][1], 0);   // q
       entier_set_si (scheduleilp->p[3+i][2], 0);   // t
       entier_set_si (scheduleilp->p[3+i][3], 0);   // s

       // l1 ... ln
       int dep_sum = 0;
       for (unsigned int j = 0; j < dimensions; j++) {
         entier_oppose (tmp, dependencies->p[i][j]);
         entier_assign (scheduleilp->p[3+i][4+j], tmp);
         dep_sum += VALUE_TO_INT( dependencies->p[i][j] );
       }

        // d1 + ... + dn
       entier_set_si (scheduleilp->p[3+i][4+dimensions], dep_sum);
       entier_set_si (scheduleilp->p[3+i][4+dimensions+1], - (int) pepipelinestages);    // const
     }

     //
     // Constraints for all vertices (used to minimize latency)
     //
     // constraint: lV_d <= s   s - lV_d >= 0    s - lV_d +_ V_dB >= 0
     //             where V_d \in { V - V' | V, V' are vertices }
     //
     int cpos = 3 + no_dependencies;
     for (int i = 0; i < no_vertices; i++) {
       for (int j = 0; j < no_vertices; j++) {
         // for each pair of vertices
         if (i == j)
           continue;

         entier_set_si (scheduleilp->p[cpos][0], 1);   // inequality
         entier_set_si (scheduleilp->p[cpos][1], 0);   // q
         entier_set_si (scheduleilp->p[cpos][2], 0);   // t
         entier_set_si (scheduleilp->p[cpos][3], 1);   // s

         // l1 ... ln
         int vert_sum = 0;
         for (unsigned int k = 0; k < dimensions; k++) {
           entier_subtract (tmp, vertices->p[i][k], vertices->p[j][k]);
           vert_sum += VALUE_TO_INT( tmp );

           entier_oppose (tmp2, tmp);
           entier_assign (scheduleilp->p[cpos][4+k], tmp2);
         }

         // d1 + ... + dn
         entier_set_si (scheduleilp->p[cpos][4+dimensions], vert_sum);
         entier_set_si (scheduleilp->p[cpos][4+dimensions+1], 0);    // const

         // go to next constraint position
         cpos++;
       }
     }

     // column position of big parameter in constraint row
     // first column (equality/inequality?) starts at index 0
     bigParamPos = dimensions + 4;


     //
     // context: B >= 0
     //
     contextilp = pip_matrix_alloc(
                       1,
                       3
                     );
     entier_set_si (contextilp->p[0][0], 1);   // inequality
     entier_set_si (contextilp->p[0][1], 1);   // B
     entier_set_si (contextilp->p[0][2], 0);   // 0

     // free memory
     entier_clear (tmp);
     entier_clear (tmp2);
   }

};

#endif // __SCHEDULE_ILP_H__
