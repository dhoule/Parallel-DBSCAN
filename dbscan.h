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

#ifndef _DBSCAN_
#define _DBSCAN_

#include "utils.h"
#include "clusters.h"

namespace NWUClustering {

  class ClusteringAlgo : public Clusters {
  public:
    ClusteringAlgo(){ }
    virtual ~ClusteringAlgo();

    void set_dbscan_params(double eps, int minPts, double seed_percentage);
    void  get_clusters_distributed();
    void  writeCluster_distributed(string outfilename);

    void  trivial_compression(vector <int>* data, vector < vector <int> >* parser);
    void  trivial_decompression(vector <int>* data);

  public:
    
    double  m_epsSquare; // AKA radius. It is the square of the "radius" given by the user
    int   m_minPts; // The minimum number of points, given by the user, to start a cluster
    int   m_compression;
    double m_perc_of_dataset; // The percentage of the points each node is to use, given by the user. Default value is 1.0(all points)

    vector <int> m_parents; // Elements hold the pointers of the clustering tree
    vector <int> m_parents_pr; // Elements hold the pointers for which node the point is in
    vector <int> m_child_count; // number of "children" a possible centroid has

    vector <int> m_member; // Values are either 0 or 1. It's size = size_of(m_pts.m_i_num_points). Used to determine if a border point or not.
    vector <int> m_corepoint; // Values are either 0 or 1. It's size = size_of(m_pts.m_i_num_points). Used to determine centroids.
  };  

  void run_dbscan_algo_uf_mpi_interleaved(ClusteringAlgo& dbs); // union find dbscan algorithm using mpi with interleaved communication
};

#endif
