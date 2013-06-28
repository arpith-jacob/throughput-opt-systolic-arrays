
//  index-enumerator.hpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 9 2009
//
//  Enumerate indices

#ifndef __INDEX_ENUMERATOR_H__
#   define __INDEX_ENUMERATOR_H__

#include <boost/math/common_factor.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;

class IndexEnumerator
{

 public:

   // constructor
   IndexEnumerator(int _dimensions, int _maxval) :
     index ( new ublas::vector<int> (_dimensions) ),
     dimensions (_dimensions),
     maxval (_maxval),
     index_lowerbound ( new ublas::vector<int> (_dimensions) ),
     index_upperbound ( new ublas::vector<int> (_dimensions) )
   {
     // initialize indices
     init();
   }

   // destructor
   ~IndexEnumerator()
   {
     delete index;
     delete index_upperbound;
     delete index_lowerbound;
   }
   
   // initialize indices
   void init()
   {
     // initialize upper bound and the first index value
     for (unsigned int i = 0; i < dimensions; i++) {
       (*index_lowerbound)(i) = - maxval;
       (*index_upperbound)(i) = maxval;

       (*index)(i) = 0;
     }
     (*index)(dimensions - 1) = 1;
   }

   // increment index
   void incr()
   {
     //
     // increment all indices
     //

     // increment last index
     (*index)(dimensions - 1)++;

     // see if we must increment the other indices as well
     for (int i = dimensions-2; i >= 0; i--) {

       // if inner index has reached upper bound, reset and increment
       // the next highest index
       if ((*index)(i + 1) > (*index_upperbound)(i + 1)) {
         (*index)(i + 1) = (*index_lowerbound)(i + 1);
         (*index)(i)++;
       } else {
         break;
       }
     }
   }

   // signal when we've reached the end of index enumeration
   bool end()
   {
     // done when the first index has reached upper bound
     if ((*index)(0) > (*index_upperbound)(0))
       return true;
     else
       return false;
   }

   // return gcd of index elements (ignore 0 elements)
   int gcd()
   {
     int vecgcd = 0;
     for (unsigned int i = 0; i < dimensions; i++) {
       if ((*index)(i) != 0)
         vecgcd = boost::math::gcd(vecgcd, (*index)(i));
     }

     return vecgcd;
   }

   // is index magnitude greater than bound?
   bool isOverBound()
   {
     int bound = maxval * maxval;
     int mag = 0;
     for (unsigned int i = 0; i < dimensions; i++) {
       mag += (*index)(i) * (*index)(i);
     }

     if (mag > bound)
       return true;
     else
       return false;
   }

   ublas::vector<int> *index;

private:

  unsigned int dimensions;
  unsigned int maxval;
  ublas::vector<int> *index_lowerbound;
  ublas::vector<int> *index_upperbound;

};

#endif // __INDEX_ENUMERATOR_H__
