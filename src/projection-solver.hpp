
//  projection-solver.hpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 9 2009
//
//  Solve throughput and schedule ILP for given projection vector

#ifndef __PROJECTION_SOLVER_H__
#   define __PROJECTION_SOLVER_H__

#include <cstdio>
#include <string>
#include <cmath>

// pip includes
#include <piplib/piplibMP.h>

// barvinok enumeration library
#include <barvinok/barvinok.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>

#include <barvinok/basis_reduction.h>

// local includes
#include "projection-solution.hpp"
#include "throughput-ilp.hpp"
#include "schedule-ilp.hpp"

class ProjectionSolver
{

 public:

   // constructor
   ProjectionSolver(int _dimensions, int _parameters,
                    vector< int > *_parameterinstantiations,
                    vector< string > *_parameternames,
                    int _pepipelinestages,
                    string polyinputfile, string dependenciesfile,
                    string verticesfile) :
     dimensions (_dimensions),
     parameters (_parameters),
     parameterinstantiations (_parameterinstantiations),
     parameternames (_parameternames),
     pepipelinestages (_pepipelinestages)

   {
     //
     // open polyhedron input file
     //
     FILE *fp = fopen (polyinputfile.c_str(), "r");

     if (!fp) {
       cerr << "Failed to open " << polyinputfile << endl;
       exit (-1);
     }

     // read in domain (unknowns) and context (parameter inequalities)
     domain = pip_matrix_read(fp);
     context = pip_matrix_read(fp);

     fclose (fp);

     //
     // open input dependencies file
     //
     fp = fopen (dependenciesfile.c_str(), "r");

     if (!fp) {
       cerr << "Failed to open " << dependenciesfile << endl;
       exit (-1);
     }

     // read in dependencies matrix
     // # rows = # dependencies
     // # columns = # dimensions
     dependencies = pip_matrix_read(fp);

     fclose (fp);

     if (dependencies->NbColumns != dimensions) {
       cerr << "Number of columns in dependencies file should equal number of dimensions" << endl;
       exit(-1);
     }

     //
     // open input vertices file
     //
     fp = fopen (verticesfile.c_str(), "r");

     if (!fp) {
       cerr << "Failed to open " << verticesfile << endl;
       exit (-1);
     }

     // read in vertices matrix
     // # rows = # vertices
     // # columns = # dimensions
     vertices = pip_matrix_read(fp);

     fclose (fp);

     if (vertices->NbColumns != dimensions) {
       cerr << "Number of columns in vertices file should equal number of dimensions" << endl;
       exit(-1);
     }

     // temporary variables for counting PEs
     COB = Matrix_Alloc( dimensions + parameters + 1,
                         dimensions + parameters + 1 );
     COBI = Matrix_Alloc( dimensions + parameters + 1,
                          dimensions + parameters + 1 );

     // parameter instantiations, used to count PEs for an instance of
     // parameters
     parameter_inst_pecount = (Value * ) malloc (sizeof (Value) * _parameters);
     for (unsigned int i = 0; i < parameters; i++) {
       value_init (parameter_inst_pecount[i]);
       entier_set_si ( parameter_inst_pecount[i],
                       (*parameterinstantiations)[i] );
     }
   }

   // destructor
   ~ProjectionSolver()
   {
     Matrix_Free (COB);
     Matrix_Free (COBI);

     for (unsigned int i = 0; i < parameters; i++) {
       value_clear (parameter_inst_pecount[i]);
     }
     free (parameter_inst_pecount);

     pip_matrix_free(domain);
     pip_matrix_free(context);
   }
   
   // find throughput (block pipelining period) for given projection vector
   ProjectionSolution *findThroughput(ublas::vector<int> *pv)
   {
     PipOptions *options;
     PipQuast   *solution;
     ProjectionSolution *ps;


     //
     // generate parameterized ILP to compute throughput for a fixed
     // projection vector
     //
     ThroughputILP ilp(domain, context, dimensions, parameters,
                       pv);

     //
     // solve throughput ILP
     //
     // initialize options for PIP solver
     options = pip_options_init();

     // call solver
     solution = pip_solve(ilp.getILP(), ilp.getContext(),
                          ilp.getBigParamPos(), options);

     // extract throughput solution for this projection
     ps = extractThroughputSolution(solution, pv);

     // print QUAST of solution
//     pip_quast_print(stdout, solution, 0);

     // free memory
     pip_options_free(options);
     pip_quast_free(solution);
     pip_close();

     return ps;
   }

   // find schedule compatible with projection vector
   void findSchedule(ProjectionSolution *ps)
   {
     PipOptions *options;
     PipQuast   *solution;

     // temporary (big num) for manipulation
     Entier tmp;
     entier_init (tmp);


     //
     // generate ILP to compute schedule compatible with projection vector
     // minimizing array utilization and latency
     //
     ScheduleILP ilp(dimensions, parameters, dependencies, vertices,
                     pepipelinestages, ps);

     //
     // solve throughput ILP
     //
     // initialize options for PIP solver
     options = pip_options_init();

     // call solver
     solution = pip_solve(ilp.getILP(), ilp.getContext(),
                          ilp.getBigParamPos(), options);

     // extract schedule solution
     int res = extractScheduleSolution(solution, ps);

     // there was no solution to the ILP.  Negate projection vector and
     // retry
     if (res < 0) {
//       cout << "No solution, trying to negate projection vector" << endl;

       // invert projection vector
       for (unsigned int i = 0; i < dimensions; i++) {
         entier_oppose (tmp, ps->projection_vector->p[0][i]);
         entier_assign (ps->projection_vector->p[0][i], tmp);
       }

       // regenerate ilp with new projection vector
       ilp.regenILP(dimensions, parameters, dependencies, vertices,
                    pepipelinestages, ps);

       // call solver
       solution = pip_solve(ilp.getILP(), ilp.getContext(),
                            ilp.getBigParamPos(), options);

       // extract schedule solution
       res = extractScheduleSolution(solution, ps);

       if (res < 0) {
         cerr << "Unable to find schedule for projection vector" << endl;
         exit (-1);
       }
     }

     // print QUAST of solution
//     pip_quast_print(stdout, solution, 0);

     // free memory
     entier_clear (tmp);
     pip_options_free(options);
     pip_quast_free(solution);
     pip_close();
   }

   //
   // compute the schedule network for a projection
   //
   // this function computes the delay of the longest communication link
   // and the sum of delays of all communication links.
   //
   void computeScheduleNetwork(ProjectionSolution *ps)
   {
     int max_delay = 0;
     int sum_delays = 0;

     // matrix multiply schedule and dependencies
     // number of dependencies is stored in dependencies->NbRows
     int prod;
     for (unsigned int i = 0; i < dependencies->NbRows; i++) {
       prod = 0;
       for (unsigned int j = 0; j < dimensions; j++) {
         prod += VALUE_TO_INT( ps->schedule->p[0][j] ) *
                   VALUE_TO_INT( dependencies->p[i][j] );
       }

       // negate delay
       prod = -prod;

       // find longest delay
       if (prod > max_delay) max_delay = prod;

       // find sum of delays of all links
       sum_delays += prod;
     }

     ps->network_sum_delays = sum_delays;
     ps->network_max_delay  = max_delay;
     ps->network_avg_delay  = (float) sum_delays / dependencies->NbRows;
   }

   //
   // compute n-1 x n allocation matrix from projection vector
   //   - n is the dimension of the projection vector (and polyhedron)
   //
   // we simply find the nullspace basis of the projection vector
   // which become the n-1 rows of the allocation matrix.
   //
   // the allocation matrix with the schedule as the last row is the
   // change of basis matrix.
   //
   // the change of basis matrix must be unimodular (determinant = +1 or -1)
   //
   void computeAllocation(ProjectionSolution *ps)
   {
     Matrix *allocation;

     // compute the integer kernel (nullspace) of the projection vector
     allocation = int_ker ( ps->projection_vector );

     if (!allocation) {
       cerr << "Failed to find nullspace of projection vector" << endl;
       exit (-1);
     }

     // nullspace must be of dimension n x n-1
     if (allocation->NbRows != dimensions ||
           allocation->NbColumns != dimensions-1) {
       cerr << "Nullspace of projection vector is of invalid dimension" << endl;
       exit (-1);
     }

     // transpose the nullspace
     allocation = Transpose ( allocation );

#if 0
     // add a row to copy schedule
     changeofbasis = AddANullRow ( changeofbasis );
     
     // copy schedule to last row
     for (unsigned int i = 0; i < dimensions; i++) {
       changeofbasis->p[dimensions-1][i] = ps->schedule->p[0][i];
     }
#endif

     // assign to projection solution
     ps->allocation = allocation;
   }

   //
   // compute the interconnection network for a projection
   //
   // this function computes the length of the longest communication link
   // and the sum of lengths of all communication links generated by the
   // allocation matrix.
   //
   void computeInterconnectionNetwork(ProjectionSolution *ps)
   {
     int max_length = 0;
     int sum_lengths = 0;

     Matrix *allocation = ps->allocation;

     // matrix multiply allocation and dependencies
     // number of dependencies is stored in dependencies->NbRows
     int prod;
     for (unsigned int i = 0; i < allocation->NbRows; i++) {
       for (unsigned int j = 0; j < dependencies->NbRows; j++) {
         prod = 0;
         for (unsigned int k = 0; k < allocation->NbColumns; k++) {
           prod += VALUE_TO_INT( allocation->p[i][k] ) * 
                     VALUE_TO_INT( dependencies->p[j][k] );
         }

         // check absolute value
         prod = abs (prod);

         // find length of the furthest link
         if (prod > max_length) max_length = prod;

         // find sum of lengths of all links
         sum_lengths += prod;
       }
     }

     ps->network_max_length  = max_length;
     ps->network_avg_length  = (float) sum_lengths / dependencies->NbRows;
   }

   //
   // count number of PEs induced by an allocation
   //
   // this function should be called after an allocation and schedule has
   // been found.
   // we first compute the change of basis matrix (allocation & schedule)
   // and transform the original polyhedron by applying the change of basis
   // transformation.
   //
   // To transform the original polyhedron P, we find Preimage (P, COB_Inverse)
   //
   // we then find the number of points in the integer projection of the
   // transformed polyhedron.
   // we can alternatively find the number of integer points in the projection
   // of the transformed polyhedron (not enabled by default).
   //
   void countPEs(ProjectionSolution *ps)
   {
     // dimension of COB matrix = # dimensions in polyhedron + # parameters +
     //                           one for the constant
     unsigned int COB_dimensions = dimensions + parameters + 1;

     // Initialize matrix to zero (matrix COB and COB inverse have been
     // allocated in constructor)
     for (unsigned int i = 0; i < COB_dimensions; i++) {
       for (unsigned int j = 0; j < COB_dimensions; j++) {
         value_set_si (COB->p[i][j], 0);
       }
     }

     // copy allocation matrix
     for (unsigned int i = 0; i < dimensions - 1; i++) {
       for (unsigned int j = 0; j < dimensions; j++) {
         value_assign (COB->p[i][j], ps->allocation->p[i][j]);
       }
     }

     // copy schedule as final dimension which will be the existential
     // variable
     for (unsigned int i = 0; i < dimensions; i++) {
       value_assign (COB->p[dimensions - 1][i], ps->schedule->p[0][i]);
     }
     
     // parameters
     for (unsigned int i = 0; i < parameters; i++) {
       value_set_si (COB->p[dimensions + i][dimensions + i], 1);
     }
     
     // constant
     value_set_si (COB->p[dimensions + parameters][dimensions + parameters], 1);

     //
     // now we need to transform the original polyhedron by applying
     // PreImage (P, COB_inverse)
     //

     // first find inverse of COB matrix
     Matrix_Inverse (COB, COBI);
//     Matrix_Print ( stdout, P_VALUE_FMT, COBI );

     Polyhedron *dom = Constraints2Polyhedron ( (Matrix *) domain, 256);

#if 0
     struct barvinok_options *bo = barvinok_options_new_with_defaults();

     // read in domain (unknowns) and context (parameter inequalities)
     static Matrix *d1 = NULL;
     
     if (!d1) d1 = Matrix_Read();
     Polyhedron *A = Constraints2Polyhedron(d1, bo->MaxRays);

     Matrix_Print ( stdout, P_VALUE_FMT, d1 );
     Polyhedron_Print(stdout, P_VALUE_FMT, A);
     Matrix *m1 = Polyhedron_Reduced_Basis (A, bo);


     // Initialize matrix to zero (matrix COB and COB inverse have been
     // allocated in constructor)
     for (unsigned int i = 0; i < COB_dimensions; i++) {
       for (unsigned int j = 0; j < COB_dimensions; j++) {
         value_set_si (COB->p[i][j], 0);
       }
     }

     // copy allocation matrix
     for (unsigned int i = 0; i < dimensions; i++) {
       for (unsigned int j = 0; j < dimensions; j++) {
         value_assign (COB->p[i][j], m1->p[i][j]);
       }
     }

     // parameters
     for (unsigned int i = 0; i < parameters; i++) {
       value_set_si (COB->p[dimensions + i][dimensions + i], 1);
     }
     
     // constant
     value_set_si (COB->p[dimensions + parameters][dimensions + parameters], 1);

     Matrix_Inverse (COB, COBI);

     // find preimage by transforming polyhedron
     Polyhedron *cobd1 = Polyhedron_Preimage (dom,
                                              COBI, 256);


     Matrix_Print ( stdout, P_VALUE_FMT, m1 );

     Polyhedron_Print(stdout, P_VALUE_FMT, cobd1);

     domain = (PipMatrix *) Polyhedron2Constraints(cobd1);

#endif

     // find preimage by transforming polyhedron
     Polyhedron *cobdom = Polyhedron_Preimage (dom,
                                               COBI, 256);

     //
     // now count number of points in the integer projection of the
     // transformed domain
     //
     // we use the barvinok library treating the final dimension (time)
     // as an existential variable
     //
     // we have 1 existential variable and #parameters
     ps->pe_count = barvinok_enumerate_e (cobdom, 1, parameters, 256);
#if 0
     Polyhedron *con = Constraints2Polyhedron ( (Matrix *) context, 256);
     ps->pe_count = barvinok_enumerate_ev (dom, con, 256);
#endif

     // compute number of pes for an instance of parameters
     ps->instance_pe_count = (int) compute_evalue (ps->pe_count,
                                                   parameter_inst_pecount);

#if 0
     Value *parameter_inst_tmp = (Value * ) malloc (sizeof (Value) * parameters);
     for (unsigned int i = 0; i < parameters; i++) {
       value_init (parameter_inst_tmp[i]);
     }

     int pesreq = 9999;
     int maxN = 10000;
     while (pesreq > 1680) {
       maxN--;
       entier_set_si ( parameter_inst_tmp[0], maxN );

       pesreq = (int) compute_evalue (ps->pe_count,
                                      parameter_inst_tmp);
     }
     ps->maxN = maxN;
#endif

     // free polyhedron
     Polyhedron_Free (dom);
     Polyhedron_Free (cobdom);
   }

private:

   int extractScheduleUnknowns(
                       PipNewparm *newparm,
                       PipVector  *pv,   // unknown's PIP solution
                       boost::rational<int> *bigparmcoeff
                      )
   {
     int unknown;

     // constant coefficient
     if (VALUE_TO_INT( pv->the_deno[1] ) != 1) {
       cerr << "Schedule solution is not integral" << endl;
       exit(-1);
     } else {
       unknown = VALUE_TO_INT( pv->the_vector[1] );
     }

     // rational for the BIG PARAMETER used internally for the maximization ILP
     bigparmcoeff->assign(
               VALUE_TO_INT( pv->the_vector[0] ),
               VALUE_TO_INT( pv->the_deno[0] )
           );

     // we should not have newparm
     if (newparm) {
       cerr << "Cannot handle newparm in schedule solution" << endl;
       exit(-1);
     }

     return unknown;     
   }

   // extract schedule and array utilization from the QUAST returned by
   // the ILP solver
   //
   // return 0 if solution extracted successfully
   //       -1 if no solution
   int extractScheduleSolution(PipQuast *solution, ProjectionSolution *ps)
   {
     //
     // no solution to schedule ILP, try negating projection vector
     //
     if (solution->list == NULL) {
       return -1;
     }


     //
     // extract solution
     //
     if (solution->condition) {
       cerr << "Cannot handle conditions" << endl;
       exit(-1);
     } else {
       rational<int> rone(1, 1);
       rational<int> rzero(0, 1);

       // first element in list is solution to the minimization problem
       PipList *pl = solution->list;
       // second element in list is solution for array utilization
       pl = pl->next;

       boost::rational<int> bigparmcoeff;

       //
       // find solution for utilization
       //
       ps->utilization = extractScheduleUnknowns(
                                        solution->newparm,
                                        pl->vector,
                                        &bigparmcoeff
                                       );

       // ensure that the BIG PARAMETER has been cancelled out
       if (bigparmcoeff != rzero) {
         cerr << "Big parameter was not eliminated in schedule ILP\n";
         exit (-1);
       }


       //
       // find solution for latency
       //
       pl = pl->next;
       bigparmcoeff.assign(0, 1);
       ps->latency = extractScheduleUnknowns(
                                        solution->newparm,
                                        pl->vector,
                                        &bigparmcoeff
                                       );

       // ensure that the BIG PARAMETER has been cancelled out
       if (bigparmcoeff != rzero) {
         cerr << "Big parameter was not eliminated in schedule ILP\n";
         exit (-1);
       }


       //
       // find solution for schedule
       //
       int sched_element;
       for (unsigned int i = 0; i < dimensions; i++) {
         // go to next solution in the list
         pl = pl->next;

         bigparmcoeff.assign(0, 1);
         sched_element = extractScheduleUnknowns(
                                        solution->newparm,
                                        pl->vector,
                                        &bigparmcoeff
                                       );

         entier_set_si (ps->schedule->p[0][i], sched_element);

         // ensure that the BIG PARAMETER has been cancelled out
         if (bigparmcoeff - rone != rzero) {
           cerr << "Big parameter was not eliminated in schedule ILP\n";
           exit (-1);
         }
       }
     }

     return 0;
   }

   void extractThroughputUnknowns(
                       PipNewparm *newparm,
                       PipVector  *pv,   // unknown's PIP solution
                       boost::rational<int> *bigparmcoeff,
                       ublas::vector< boost::rational<int> > *unknown
                      )
   {
     // rational for the constant coefficient
     rational<int> r(
               VALUE_TO_INT( pv->the_vector[pv->nb_elements - 1] ),
               VALUE_TO_INT( pv->the_deno[pv->nb_elements - 1] )
           );
     (*unknown)(parameters) = r;

     // rational for the BIG PARAMETER used internally for the maximization ILP
     bigparmcoeff->assign(
               VALUE_TO_INT( pv->the_vector[parameters] ),
               VALUE_TO_INT( pv->the_deno[parameters] )
           );

     // set parameter multipliers
     for (unsigned int i = 0; i < parameters; i++) {
       (*unknown)(i).assign(
               VALUE_TO_INT( pv->the_vector[i] ),
               VALUE_TO_INT( pv->the_deno[i] )
           );
     }

     // set multipliers of all new parameters
     int multiplier, divider;
     while (newparm) {
       // multiplier co-efficient at rank new parameter
       multiplier = VALUE_TO_INT( pv->the_vector[newparm->rank] );

       // divider taken from newparm structure
       divider = VALUE_TO_INT( newparm->deno );

       for (unsigned int i = 0; i < parameters; i++) {
         rational<int> r(
                         VALUE_TO_INT( newparm->vector->the_vector[i] )
                            * multiplier,
                         VALUE_TO_INT( newparm->vector->the_deno[i] )
                            * divider
                        );
         (*unknown)(i) += r;
       }

       // compute constant coefficient
       int coeffrank = newparm->vector->nb_elements - 1;
       rational<int> r2(
                       VALUE_TO_INT( newparm->vector->the_vector[coeffrank] )
                          * multiplier,
                       VALUE_TO_INT( newparm->vector->the_deno[coeffrank] )
                          * divider
                      );
       (*unknown)(parameters) += r2;

       // compute BIG PARAMETER coefficient
       int bigparmrank = newparm->vector->nb_elements - 2;
       rational<int> r1(
                       VALUE_TO_INT( newparm->vector->the_vector[bigparmrank] )
                          * multiplier,
                       VALUE_TO_INT( newparm->vector->the_deno[bigparmrank] )
                          * divider
                      );
       *bigparmcoeff += r1;

       // go to next new parameter
       newparm = newparm->next;
     }
     
   }

   // extract the BPP, x1, x2 (projection solution) from the QUAST
   // returned by the ILP solver
   ProjectionSolution *extractThroughputSolution(PipQuast *solution,
                                                 ublas::vector<int> *pv)
   {
     //
     // no solution?  don't see how this is possible :(
     //
     if (solution->list == NULL) {
       cerr << "Throughput ILP yielded solution with condition, cannot handle" << endl;
       pip_quast_print(stdout, solution, 0);
       exit(-1);
     }

     // Projection solution
     // Solution is in terms of parameters and const
     ProjectionSolution *ps = new ProjectionSolution(dimensions,
                                                     parameters,
                                                     parameterinstantiations,
                                                     parameternames
                                                     );

     //
     // assign projection vector
     //
     for (unsigned int i = 0; i < dimensions; i++) {
       entier_set_si ( ps->projection_vector->p[0][i], (*pv)[i] );
     }

     //
     // extract solution
     //
     if (solution->condition) {
       cerr << "Cannot handle conditions" << endl;
       exit(-1);
     } else {
       rational<int> rone(1, 1);
       rational<int> rzero(0, 1);

       PipList *pl = solution->list;
       boost::rational<int> bigparmcoeff;

       //
       // find solution for k_max from k'
       //
       extractThroughputUnknowns(
                      solution->newparm,
                      pl->vector,
                      &bigparmcoeff,
                      ps->bpp
                     );

       // negate coefficients
       for (unsigned int i = 0; i <= parameters; i++) {
         (*ps->bpp)(i) = -(*ps->bpp)(i);
       }

       // ensure that the BIG PARAMETER has been cancelled out
       if (rone - bigparmcoeff != rzero) {
         cerr << "Big parameter was not eliminated in throughput ILP\n";
         exit (-1);
       }

       //
       // find solution for x1
       //
       for (unsigned int i = 0; i < dimensions; i++) {
         // advance to next unknown solution in list
         pl = pl->next;

         // extract X1[i]
         bigparmcoeff.assign(0, 1);
         extractThroughputUnknowns(
                        solution->newparm,
                        pl->vector,
                        &bigparmcoeff,
                        ps->x1[i]
                       );

         // ensure that the BIG PARAMETER coefficient is zero
         if (bigparmcoeff != rzero) {
           cerr << "Big parameter was not eliminated in throughput ILP\n";
           exit (-1);
         }
       }


       //
       // find solution for x2
       //
       for (unsigned int i = 0; i < dimensions; i++) {
         // advance to next unknown solution in list
         pl = pl->next;

         // extract X1[i]
         bigparmcoeff.assign(0, 1);
         extractThroughputUnknowns(
                        solution->newparm,
                        pl->vector,
                        &bigparmcoeff,
                        ps->x2[i]
                       );

         // ensure that the BIG PARAMETER coefficient is zero
         if (bigparmcoeff != rzero) {
           cerr << "Big parameter was not eliminated in throughput ILP\n";
           exit (-1);
         }
       }

     }
     
     return ps;
   }
   
  // number of dimensions and parameters in the input polyhedron
  unsigned int dimensions;
  unsigned int parameters;
  vector< int > *parameterinstantiations;
  vector< string > *parameternames;

  // # of pipeline stages.  Min. delay on each dependency link
  unsigned int pepipelinestages;

  // store polyhedron constraints
  PipMatrix *domain, *context;
  
  // store input dependencies
  PipMatrix *dependencies;

  // store polyhedron vertices
  PipMatrix *vertices;

  // temporary store for change of basis matrix and its inverse
  // used for counting number of points (processing elements) in a projected
  // domain
  Matrix *COB, *COBI;

  // temporary variables to hold parameter instantiations in data type
  // for computing number of PEs using the barvinok library
  Value *parameter_inst_pecount;

};

#endif // __PROJECTION_SOLVER_H__
