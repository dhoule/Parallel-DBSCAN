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

#include <unistd.h>
#include "dbscan.h"
#include "utils.h"
#include "kdtree2.hpp"
#include "geometric_partitioning.h"

static void usage(char *argv0) {
  const char *params =
    "Usage: %s [switches] -i filename -b -m minpts -e epsilon -o output [-k seed_percentage]\n"
    " -i filename : file containing input data to be clustered\n"
    " -b    : input file is in binary format (default, binary and currently the only supported format)\n"
    " -m minpts : input parameter of DBSCAN, min points to form a cluster, e.g. 2\n"
    " -e epsilon  : input parameter of DBSCAN, radius or threshold on neighbourhoods retrieved, e.g. 0.8\n"
    " -o output : clustering results, format, (points coordinates, cluster id)\n"
    " -k seed_percentage : the percentage of points for each node to use; range is [0.0,1.0), with default value of 1.0\n\n";
  
  fprintf(stderr, params, argv0);
  exit(-1);
}

// TODO splicing compression technique, using "Rem's union-find technique": https://algocoding.wordpress.com/2015/05/13/simple-union-find-techniques/#more-1395
int main(int argc, char** argv) {
  int   opt;
  int   minPts;
  double  eps, start, seed_percentage, preprocessing_start, actual_start;
  char*   outfilename;
  int     isBinaryFile;
  char*   infilename;
  int rank, nproc;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // some default values
  minPts = -1;
  eps = -1;
  seed_percentage = 1.0;
  isBinaryFile = 1; // default binary file
  outfilename = NULL;
  infilename = NULL;

  // determine command line options 
  while ((opt=getopt(argc,argv,"i:m:e:o:k:?b"))!= EOF) {
    switch (opt) {
      case 'i':
        infilename = optarg;
        break;
      case 'b':
        isBinaryFile = 1;
        break;
      case 'm':
        minPts = atoi(optarg);
        break;
      case 'e':
        eps  = atof(optarg);
        break;
      case 'o':
        outfilename = optarg;
        break;
      case 'k':
        seed_percentage = atof(optarg);
        break;
      case '?':
        usage(argv[0]);
        break;
      default:
        usage(argv[0]);
        break;
    }
  }
  // filter command line input
  if(infilename == NULL || minPts < 0 || eps < 0) {
    if(rank == proc_of_interest)
      usage(argv[0]);
    MPI_Finalize();
    exit(-1);
  }
  // 'proc_of_interest' defined in utils.h as process 0, the master node
  if(rank == proc_of_interest) cout << "Number of process cores " << nproc << endl;

  // check if `nproc` is NOT multiple of TWO
  unsigned int proc_count = nproc;
  // filter command line input, again
  while (((proc_count % 2) == 0) && proc_count > 1) // While x is even and > 1
    proc_count /= 2;
  
  if(proc_count != 1) {
    if(rank == proc_of_interest) cout << "\n\nPlease use the number of process cores as a multiple of TWO" << endl;
    MPI_Finalize();
    return 0;
  }

  // make sure `seed_percentage` is within the specified range
  if((0.0 >= seed_percentage) || (1.0 < seed_percentage)) {
    if(rank == proc_of_interest) cout << "\n\nFor -k, please use a percentage greater than 0.0 and less than or equal to 1.0." << endl;
    MPI_Finalize();
    return 0;
  }
  // declaring the ClusteringAlgo object 'dbs'
  NWUClustering::ClusteringAlgo dbs;
  // initialize some paramaters
  dbs.set_dbscan_params(eps, minPts, seed_percentage); // TODO need to modify this for SNG Alg

  if(rank == proc_of_interest) cout << "Epsilon: " << eps << " MinPts: " << minPts << endl;
  // Make ALL of the nodes/processes wait till they ALL get to this point
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  
  if(rank == proc_of_interest) cout << "Reading points from file: " << infilename << endl;
  // determine if there was an error reading from the binary file
  if(dbs.read_file(infilename, isBinaryFile) == -1) {
    MPI_Finalize();
    exit(-1);
  }

  if(rank == proc_of_interest) cout << "Reading the input data file took " << MPI_Wtime() - start << " seconds [pre_processing]"<< endl;
  // Make ALL of the nodes/processes wait till they ALL get to this point
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  preprocessing_start = start; // used to calculate total preprocessing time
  // parttition the data file geometrically: preprocessing_start, actual_start step
  start_partitioning(dbs);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == proc_of_interest) cout << "Partitioning the data geometrically took " << MPI_Wtime() - start << " seconds [pre_processing]" << endl;
  
  // gather extra(outer) points that falls within the eps radius from the boundary: preprocessing_start, actual_start step
  start = MPI_Wtime();
  get_extra_points(dbs);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == proc_of_interest) cout << "Gathering extra point took " << MPI_Wtime() - start << " seconds [pre_processing]" << endl;
    
  // build the kdtrees: preprocessing_start, actual_start step 
  start = MPI_Wtime();
  dbs.build_kdtree();
  dbs.build_kdtree_outer();
  // TODO call GetSeeds() from here
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == proc_of_interest) cout << "Build kdtree took " << MPI_Wtime() - start << " seconds [pre_processing]\n" << endl;
  if(rank == proc_of_interest) cout << "\nTotal preprocessing time: " << MPI_Wtime() - preprocessing_start << " seconds\n" << endl;
  //run the DBSCAN algorithm
  start = MPI_Wtime();
  actual_start = start;
  run_dbscan_algo_uf_mpi_interleaved(dbs);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == proc_of_interest) cout << "Parallel DBSCAN (init, local computation, and merging) took " << MPI_Wtime() - start << " seconds\n"<< endl;

  // assign cluster IDs to points
  start = MPI_Wtime();
  dbs.get_clusters_distributed(); 
  if(rank == proc_of_interest) cout << "Assigning cluster IDs to points " << MPI_Wtime() - start << " seconds [post_processing]" << endl;
  if(rank == proc_of_interest) cout << "\nTotal time for actual algorithm (including assigning cluster IDs): " << MPI_Wtime() - actual_start << " seconds\n" << endl;
  if(rank == proc_of_interest) cout << "\nTotal time for evrything: " << MPI_Wtime() - preprocessing_start << " seconds\n" << endl;


  if(outfilename != NULL) {
    start = MPI_Wtime();  

    // activate the following line to write the cluster to file
    dbs.writeCluster_distributed(outfilename);
    if(rank == proc_of_interest) cout << "Writing clusterIDs to disk took " << MPI_Wtime() - start << " seconds [pre_processing]"<< endl;
  }
    
  MPI_Finalize();
  return 0;
}
