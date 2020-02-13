/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   Files: mpi_main.cpp clusters.cpp  clusters.h utils.h utils.cpp          */
/*        dbscan.cpp dbscan.h kdtree2.cpp kdtree2.hpp          */
/*      geometric_partitioning.h geometric_partitioning.cpp  */
/*                                 */
/*   Description: an mpi implementation of dbscan clustering algorithm       */
/*        using the disjoint set data structure        */
/*                                                                           */
/*   Author:  Md. Mostofa Ali Patwary                                        */
/*            EECS Department, Northwestern University                       */
/*            email: mpatwary@eecs.northwestern.edu                          */
/*                                                                           */
/*   Copyright, 2012, Northwestern University                                */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/*   Please cite the following publication if you use this package       */
/*                       */
/*   Md. Mostofa Ali Patwary, Diana Palsetia, Ankit Agrawal, Wei-keng Liao,  */
/*   Fredrik Manne, and Alok Choudhary, "A New Scalable Parallel DBSCAN      */
/*   Algorithm Using the Disjoint Set Data Structure", Proceedings of the    */
/*   International Conference on High Performance Computing, Networking,     */
/*   Storage and Analysis (Supercomputing, SC'12), pp.62:1-62:11, 2012.      */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "utils.h"

// Find Kth element without recusion
float findKMedian(vector<float>& A, int K) {
  int l,m;
  l=0;
  m=A.size()-1;
  while (l < m) {
    // NOTE: These variables need to stay here, because of the swap() lower down
    float x=A[K];
    int i=l;
    int j=m;
    do {
      while (A[i] < x) // start from the front, and move towards the back
        i++;
      while (x < A[j]) // start from the back, and move towards the front
        j--;
      if (i <= j) {
        /*
          Exchanges the values of a and b.
          TODO: moved from <algorithm> to <utility> in C++11
        */
        swap(A[i], A[j]);
        i++; 
        j--;
      }
    } while (i <= j);

    if (j<K) l=i;
    if (K<i) m=j;
  }

  return A[K];
}


