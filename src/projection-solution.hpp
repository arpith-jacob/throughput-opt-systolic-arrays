
//  projection-solution.hpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 8 2009
//
//  Projection solution container.  Contains projection vector, best block
//  pipelining period (BPP) possible (k_max), index points x1, x2 that produce
//  this BPP.

#ifndef __PROJECTION_SOLUTION_H__
#   define __PROJECTION_SOLUTION_H__

#include <polylib/polylibgmp.h>

#include <barvinok/evalue.h>

#include <boost/rational.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;

class ProjectionSolution
{

 public:

   // constructor
   ProjectionSolution(int _dimensions, int _parameters,
                      vector< int > *_parameterinstantiations,
                      vector< string > *_parameternames) :
     projection_vector ( Matrix_Alloc (1, _dimensions) ),
     bpp ( new ublas::vector< boost::rational<int> > (_parameters + 1) ),
     utilization (0),
     latency (0),
     schedule ( Matrix_Alloc (1, _dimensions) ),
     allocation ( NULL),
     network_max_length (0),
     network_avg_length (0.),
     pe_count (NULL),
     instance_pe_count (0),
     dimensions (_dimensions),
     parameters (_parameters),
     parameterinstantiations (_parameterinstantiations),
     parameternames (_parameternames)
   {
     // initialize x1, x2
     x1 = new ublas::vector< boost::rational<int> > * [_dimensions];
     x2 = new ublas::vector< boost::rational<int> > * [_dimensions];
     
     for (unsigned int i = 0; i < dimensions; i++) {
       x1[i] = new ublas::vector< boost::rational<int> > (_parameters + 1);
       x2[i] = new ublas::vector< boost::rational<int> > (_parameters + 1);
     }
   }

   // destructor
   ~ProjectionSolution()
   {
     if (parameternames)
       delete parameternames;

     if (parameterinstantiations)
       delete parameterinstantiations;

     if (pe_count)
       free_evalue_refs (pe_count);

     if (allocation)
       Matrix_Free (allocation);

     if (schedule)
       Matrix_Free (schedule);

     for (unsigned int i = 0; i < dimensions; i++) {
       delete x1[i];
       delete x2[i];
     }
     delete x1;
     delete x2;

     if (bpp)
       delete bpp;

     if (projection_vector)
       Matrix_Free (projection_vector);
   }

   // compute BPP for an instance of the parameters
   void computeInstanceBPP()
   {
     float bpp_i = 0;

     // multiply instance of each parameter with rational coefficient of BPP
     for (unsigned int i = 0; i < parameters; i++) {
       bpp_i += (*parameterinstantiations)[i] *
                (float)  (*bpp)(i).numerator() / (float) (*bpp)(i).denominator();
     }
     
     // add constant coefficient
     bpp_i += (float)  (*bpp)(parameters).numerator() /
                (float) (*bpp)(parameters).denominator();
     
     instance_bpp = (int) bpp_i;
     
     // round up
     if ((float) instance_bpp != bpp_i)
       instance_bpp++;
   }

   // print entire solution
   void print()
   {
     printProjectionVector();
     printBPP();
       //printInstanceBPP();
     printPECount();
     printInstancePECount();
//     printX1();
//     printX2();
     printSchedule();
     printUtil();
     printScheduleDelays();
     printLatency();
     printAllocation();
     printNetwork();
   }

   // print projection vector
   void printProjectionVector()
   {
     cout << "\"";
     for (unsigned int i = 0; i < projection_vector->NbColumns; i++) {
       cout << VALUE_TO_INT (projection_vector->p[0][i]) << " ";
     }
     cout << "\",";
   }

   // print BPP
   void printBPP()
   {
     cout << "\"";
     for (unsigned int i = 0; i < parameters; i++) {
       cout << (*bpp)(i) << (*parameternames)[i] << " + ";
     }
     cout << (*bpp)(parameters) << " + 1";
     cout << "\",";
   }

   // print instance BPP
   void printInstanceBPP()
   {
     cout << instance_bpp + 1 << ",";
   }

   // print number of PEs
   void printPECount()
   {
     const char **param_name;
     // TODO: this is inefficient, but we don't really care much about
     // printing this value
     param_name = (const char **) malloc ( sizeof (char *) * parameters );

     for (unsigned int i = 0; i < parameters; i++) {
       param_name[i] = (*parameternames)[i].c_str();
     }

     cout << "\"";
     print_evalue(stdout, pe_count, param_name);
     cout << "\",";
     
     free (param_name);
   }

   // print number of PEs for 
   void printInstancePECount()
   {
     cout << instance_pe_count << ", ";
//     cout << " (maxN=" << maxN << "), ";
   }

   // print x1
   void printX1()
   {
     cout << "X1: \n";
     for (unsigned int i = 0; i < dimensions; i++) {
       cout << "  " << i << ": ";
       for (unsigned int j = 0; j < parameters; j++) {
         cout << (*parameternames)[j] << (*x1[i])(j) << " + ";
       }
       cout << (*x1[i])(parameters) << endl;
     }
   }

   // print x2
   void printX2()
   {
     cout << "X2: \n";
     for (unsigned int i = 0; i < dimensions; i++) {
       cout << "  " << i << ": ";
       for (unsigned int j = 0; j < parameters; j++) {
         cout << (*parameternames)[j] << (*x2[i])(j) << " + ";
       }
       cout << (*x2[i])(parameters) << endl;
     }
   }

   // print utilization
   void printUtil()
   {
     cout << utilization << ",";
   }

   // print latency
   void printLatency()
   {
     cout << latency << ", ";
   }

   // print schedule
   void printSchedule()
   {
     cout << "\"";
     for (unsigned int i = 0; i < schedule->NbColumns; i++) {
       cout << VALUE_TO_INT (schedule->p[0][i]) << " ";
     }
     cout << "\",";
   }

   // print network delays
   void printScheduleDelays()
   {
     cout << network_sum_delays << ", " << network_avg_delay << ", " << network_max_delay << ", ";
   }

   // print allocation matrix
   void printAllocation()
   {
     cout << "\"";
     for (unsigned int i = 0; i < dimensions - 1; i++) {
       cout << "[ ";
       for (unsigned int j = 0; j < dimensions; j++) {
         cout << VALUE_TO_INT (allocation->p[i][j]) << " ";
       }
       cout << "]";
     }
     cout << "\",";
   }

   // print interconnection network info
   void printNetwork()
   {
     //cout << network_avg_length;
     cout << network_avg_length << ", " << network_max_length;
   }

   Matrix *projection_vector;
   ublas::vector< boost::rational<int> > *bpp;
   ublas::vector< boost::rational<int> > **x1;
   ublas::vector< boost::rational<int> > **x2;
   unsigned int instance_bpp;
   unsigned int utilization;
   unsigned int latency;
   Matrix *schedule;
   unsigned int network_sum_delays;
   unsigned int network_max_delay;
   float network_avg_delay;
   Matrix *allocation;
   unsigned int network_max_length;
   float network_avg_length;
   evalue *pe_count;
   unsigned int instance_pe_count;
//   unsigned int maxN;

private:

   unsigned int dimensions;
   unsigned int parameters;
   vector< int > *parameterinstantiations;
   vector< string > *parameternames;

};

#endif // __PROJECTION_SOLUTION_H__
