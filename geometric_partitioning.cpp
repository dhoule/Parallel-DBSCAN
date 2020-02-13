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

#include "geometric_partitioning.h"

namespace NWUClustering {
  /*
    Called from mpi_main.cpp
    Gathers extra(outer) points that falls within the eps radius from the boundary
  */
  void get_extra_points(ClusteringAlgo& dbs) {

    int k, i, j, rank, nproc, thisNode, thatNode, dims = dbs.m_pts->m_i_dims; // `thisNode` & `thatNode` are used to speed up computation by only doing it once
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    // `local_box` size is the number of dimensions of the dataset
    interval* local_box = new interval[dims];
    compute_local_bounding_box(dbs, local_box);

    /* 
      Extend the bounding box of current node by a distance of `eps` in every direction,
      in each dimension and query other processors with the extended bounding box to 
      return their local points that falls within it.
    */
    float eps = sqrt(dbs.m_epsSquare);

    for(i = 0; i < dims; i++) {
      local_box[i].upper += eps;
      local_box[i].lower -= eps; 
    }
  
    // all together all the extending bounding box
    interval* gather_local_box = new interval[dims * nproc];

    /*
      gather the local bounding box first
      Gathers data from all processes and distributes it to all processes
      int MPI_Allgather(const void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
    */
    int sendcount = sizeof(interval) * dims; // instead of doing the calculation twice, just do it once...
    MPI_Allgather(local_box, sendcount, MPI_BYTE, gather_local_box, sendcount, MPI_BYTE, MPI_COMM_WORLD);

    bool if_inside, overlap;
    int count = 0, gcount;

    vector <float> empty;
    vector <vector <float> > send_buf;
    vector <vector <float> > recv_buf;
    send_buf.resize(nproc, empty);
    recv_buf.resize(nproc, empty);

    vector <int> empty_i;
    vector <vector <int> > send_buf_ind;
    vector <vector <int> > recv_buf_ind;
    send_buf_ind.resize(nproc, empty_i);
    recv_buf_ind.resize(nproc, empty_i);

    /*
      Looping over `gather_local_box`, which is the size of the number of nodes in the system.
      This is to get the "out box" points from other nodes.
      This is still considered a preprocessing step.
    */
    thisNode = rank * dims; // only needs to be computed once, but is used often
    for(k = 0; k < nproc; k++) {
      // makes no sense for a node to compare itself to itself. Move to next itteration.
      if(k == rank) 
        continue;

      // check the two extended bounding box of proc `rank` and `k`. If they don't overlap, there must be no points      

      overlap = true; // just a flag
      /* 
        Looping over each dimension, to determine if there is any overlap in the "bounding boxes"
        If there is no overlap, the `overlap` flag is switched to FALSE.
      */
      thatNode = k * dims; // Is used several times, so only need to calculate it once
      for(j = 0; j < dims; j++) {
        if(gather_local_box[thisNode + j].lower < gather_local_box[thatNode +j].lower) {
          if(gather_local_box[thisNode + j].upper - gather_local_box[thatNode + j].lower < eps) {
            overlap = false;
            break;
          }
        } else {
          if(gather_local_box[thatNode + j].upper - gather_local_box[thisNode + j].lower < eps) {
            overlap = false;
            break;
          }
        }
      }

      // the two bouding boxes are different, so continue to the next processors
      if(overlap == false)
        continue;

      // get the overlapping regions
      for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
        if_inside = true; // just a flag
        for(j = 0; j < dims; j++) {
          // determine of the current point, is within another bounding box or not.
          if(dbs.m_pts->m_points[i][j] < gather_local_box[thatNode + j].lower || dbs.m_pts->m_points[i][j] > gather_local_box[thatNode + j].upper) {
            if_inside = false;
            break;
          }
        }
        // if the flag was never tripped, put the points into the send buffer
        if(if_inside == true) {
          for(j = 0; j < dims; j++)
            send_buf[k].push_back(dbs.m_pts->m_points[i][j]);
          // TODO no idea what this does
          send_buf_ind[k].push_back(i);
          count++;
        }
      }
    }

    // now send buf have all the points. Send the size first to everyone
    vector <int> send_buf_size, recv_buf_size;
    send_buf_size.resize(nproc, 0);
    recv_buf_size.resize(nproc, 0);

    for(i = 0; i < nproc; i++)
      send_buf_size[i] = send_buf[i].size();

    /*
      All processes send data to all processes.
      Every node shares the number of "outter box" points they have found.
      int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
    */
    MPI_Alltoall(&send_buf_size[0], 1, MPI_INT, &recv_buf_size[0], 1, MPI_INT, MPI_COMM_WORLD);

    //return;
    int tag = 200, tagPlusOne = 201, send_count, recv_count;
    MPI_Request req_send[2 * nproc], req_recv[2 * nproc];
    MPI_Status stat_send[2 * nproc], stat_recv;

    recv_count = 0;
    for(i = 0; i < nproc; i++) {
      // only continue if the current element actually has elements in it
      if(recv_buf_size[i] > 0) {
        recv_buf[i].resize(recv_buf_size[i], 0);
        recv_buf_ind[i].resize(recv_buf_size[i] / dims, -1);
        /*
          Starts a standard-mode, nonblocking receive.
          int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
        */
        MPI_Irecv(&recv_buf[i][0], recv_buf_size[i], MPI_FLOAT, i, tag, MPI_COMM_WORLD, &req_recv[recv_count++]);
        MPI_Irecv(&recv_buf_ind[i][0], recv_buf_size[i] / dims, MPI_INT, i, tagPlusOne, MPI_COMM_WORLD, &req_recv[recv_count++]);
      }
    }

    send_count = 0;
    for(i = 0; i < nproc; i++) {
      if(send_buf_size[i] > 0) {
        /*
          Starts a standard-mode, nonblocking send.
          int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
        */
        MPI_Isend(&send_buf[i][0], send_buf_size[i], MPI_FLOAT, i, tag, MPI_COMM_WORLD, &req_send[send_count++]);
        MPI_Isend(&send_buf_ind[i][0], send_buf_size[i] / dims, MPI_INT, i, tagPlusOne, MPI_COMM_WORLD, &req_send[send_count++]);
      }
    }

    int rtag, rsource, rpos;

    // Sets up the Points_Outer object
    dbs.allocate_outer(dims);

    for(i = 0; i < recv_count; i++) {
      /*
        Waits for any specified send or receive to complete.
        int MPI_Waitany(int count, MPI_Request array_of_requests[], int *index, MPI_Status *status)
      */
      MPI_Waitany(recv_count, &req_recv[0], &rpos, &stat_recv);
    
      rtag = stat_recv.MPI_TAG;
      rsource = stat_recv.MPI_SOURCE;
    
      if(rtag == tag) {
        // process the request
        dbs.addPoints(rsource, recv_buf_size[rsource], dims, recv_buf[rsource]);

        recv_buf[rsource].clear();
      } else if(rtag == tagPlusOne) {
        // postpond this computation and call update points later
        // processing immediately might lead to invalid computation
      }
    } 
    
    // MAY NOT NEED THIS
    if(send_count > 0)
      MPI_Waitall(send_count, &req_send[0], &stat_send[0]);

    // got all the points
    // now update the indices of the outer points

    dbs.updatePoints(recv_buf_ind);
    /*
      TODO need to look more into this. `count` & `gcount` are local variables.
      Reduces values on all processes within a group.
      int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
    */
    MPI_Reduce(&count, &gcount, 1, MPI_INT, MPI_SUM, proc_of_interest, MPI_COMM_WORLD);
    // Clear the vectors manually to free up memory
    empty.clear();
    send_buf.clear();
    recv_buf.clear();
    send_buf_size.clear();
    recv_buf_size.clear();
    send_buf_ind.clear();
    recv_buf_ind.clear();

    delete [] gather_local_box;
    delete [] local_box;
  }

  /*
    Called from mpi_main.cpp
    Parttition the data file geometrically: preprocessing step
    Sets up new communication systems for each node.
    The "partners" for each node are specified when they are set up.
  */
  void start_partitioning(ClusteringAlgo& dbs) {
    int r_count, s_count, rank, nproc, i, j, k, dims = dbs.m_pts->m_i_dims;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // compute the local bouding box for each dimention
    interval* box = new interval[dims];
    
    for(i = 0; i < dims; i++) { 
      box[i].upper = dbs.m_pts->m_box[i].upper;
      box[i].lower = dbs.m_pts->m_box[i].lower;
    }

    // compute the global bouding box for each dimention
    interval* gbox = new interval[dims];
    compute_global_bounding_box(dbs, box, gbox, nproc);

    // find the loop count for nproc processors
    int internal_nodes, partner_rank, loops, b, color, sub_rank, d, max, sub_nprocs;

    MPI_Status status;

    loops = 0;
    i = nproc;
    internal_nodes = 1;
    
    while((i = i >> 1) > 0) { 
      loops++; 
      internal_nodes = internal_nodes << 1;
    }
    internal_nodes = internal_nodes << 1;
    /*
      From the prior while loop above:
      16 nodes:               8 nodes:
      loops = 4               loops = 3
      internal_nodes = 32     internal_nodes = 16

      4 nodes:                2 nodes:
      loops = 2               loops = 1
      internal_nodes = 8      internal_nodes = 4
    */

    //gbox for each node in the tree [ONLY upto to reaching each processor]
    interval** nodes_gbox = new interval*[internal_nodes];
    for(i = 0; i < internal_nodes; i++)
      nodes_gbox[i] = new interval[dims];
    
    copy_global_box_to_each_node(dbs, nodes_gbox, gbox, internal_nodes);
  
    vector <float> send_buf;
    vector <int>   invalid_pos_as;
    vector <float> recv_buf;
        
    int pow2_i, powColor;
    float median;

    for(i = 0; i < loops; i++) { 
      /*
        These are the values of the below variables, given the number of nodes.
        Nodes 2
        rank: 0           rank: 1
        i: 0              i: 0
        pow2_i: 1         pow2_i: 1
        b: 0              b: 0
        color: 0          color: 0
        partner_rank: 1   partner_rank: 0
        ==========================
        Nodes 4
        rank: 0          rank: 0          rank: 1          rank: 1          
        i: 0             i: 1             i: 0             i: 1             
        pow2_i: 1        pow2_i: 2        pow2_i: 1        pow2_i: 2        
        b: 0             b: 2             b: 0             b: 2             
        color: 0         color: 0         color: 0         color: 0         
        partner_rank: 2  partner_rank: 1  partner_rank: 3  partner_rank: 0  
        _______________________________________________________________________________________________________________________________________________
        rank: 2          rank: 2          rank: 3          rank: 3
        i: 0             i: 1             i: 0             i: 1
        pow2_i: 1        pow2_i: 2        pow2_i: 1        pow2_i: 2
        b: 0             b: 2             b: 0             b: 2
        color: 0         color: 1         color: 0         color: 1
        partner_rank: 0  partner_rank: 3  partner_rank: 1  partner_rank: 2
        ==========================
        Nodes 8
        rank: 0          rank: 0          rank: 0          rank: 1          rank: 1          rank: 1          
        i: 0             i: 1             i: 2             i: 0             i: 1             i: 2             
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 1        pow2_i: 2        pow2_i: 4        
        b: 0             b: 4             b: 6             b: 0             b: 4             b: 6             
        color: 0         color: 0         color: 0         color: 0         color: 0         color: 0         
        partner_rank: 4  partner_rank: 2  partner_rank: 1  partner_rank: 5  partner_rank: 3  partner_rank: 0  
        _______________________________________________________________________________________________________________________________________________
        rank: 2          rank: 2          rank: 2          rank: 3          rank: 3          rank: 3          
        i: 0             i: 1             i: 2             i: 0             i: 1             i: 2             
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 1        pow2_i: 2        pow2_i: 4        
        b: 0             b: 4             b: 6             b: 0             b: 4             b: 6             
        color: 0         color: 0         color: 1         color: 0         color: 0         color: 1         
        partner_rank: 6  partner_rank: 0  partner_rank: 3  partner_rank: 7  partner_rank: 1  partner_rank: 2  
        _______________________________________________________________________________________________________________________________________________
        rank: 4          rank: 4          rank: 4          rank: 5          rank: 5          rank: 5
        i: 0             i: 1             i: 2             i: 0             i: 1             i: 2
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 1        pow2_i: 2        pow2_i: 4
        b: 0             b: 4             b: 6             b: 0             b: 4             b: 6
        color: 0         color: 1         color: 2         color: 0         color: 1         color: 2
        partner_rank: 0  partner_rank: 6  partner_rank: 5  partner_rank: 1  partner_rank: 7  partner_rank: 4
        _______________________________________________________________________________________________________________________________________________
        rank: 6          rank: 6          rank: 6          rank: 7          rank: 7          rank: 7
        i: 0             i: 1             i: 2             i: 0             i: 1             i: 2
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 1        pow2_i: 2        pow2_i: 4
        b: 0             b: 4             b: 6             b: 0             b: 4             b: 6
        color: 0         color: 1         color: 3         color: 0         color: 1         color: 3
        partner_rank: 2  partner_rank: 4  partner_rank: 7  partner_rank: 3  partner_rank: 5  partner_rank: 6
        ==========================
        Nodes 16
        rank: 0          rank: 0          rank: 0          rank: 0          rank: 1          rank: 1          rank: 1          rank: 1
        i: 0             i: 1             i: 2             i: 3             i: 0             i: 1             i: 2             i: 3
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8
        b: 0             b: 8             b: 12            b: 14            b: 0             b: 8             b: 12            b: 14
        color: 0         color: 0         color: 0         color: 0         color: 0         color: 0         color: 0         color: 0
        partner_rank: 8  partner_rank: 4  partner_rank: 2  partner_rank: 1  partner_rank: 9  partner_rank: 5  partner_rank: 3  partner_rank: 0
        _______________________________________________________________________________________________________________________________________________
        rank: 2          rank: 2          rank: 2          rank: 2          rank: 3          rank: 3          rank: 3          rank: 3
        i: 0             i: 1             i: 2             i: 3             i: 0             i: 1             i: 2             i: 3
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8 
        b: 0             b: 8             b: 12            b: 14            b: 0             b: 8             b: 12            b: 14 
        color: 0         color: 0         color: 0         color: 1         color: 0         color: 0         color: 0         color: 1         
        partner_rank: 10 partner_rank: 6  partner_rank: 0  partner_rank: 3  partner_rank: 11 partner_rank: 7  partner_rank: 1  partner_rank: 2
        _______________________________________________________________________________________________________________________________________________
        rank: 4          rank: 4          rank: 4          rank: 4          rank: 5          rank: 5          rank: 5          rank: 5
        i: 0             i: 1             i: 2             i: 3             i: 0             i: 1             i: 2             i: 3
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8
        b: 0             b: 8             b: 12            b: 14            b: 0             b: 8             b: 12            b: 14
        color: 0         color: 0         color: 1         color: 2         color: 0         color: 0         color: 1         color: 2 
        partner_rank: 12 partner_rank: 0  partner_rank: 6  partner_rank: 5  partner_rank: 13 partner_rank: 1  partner_rank: 7  partner_rank: 4
        _______________________________________________________________________________________________________________________________________________
        rank: 6          rank: 6          rank: 6          rank: 6          rank: 7          rank: 7          rank: 7          rank: 7
        i: 0             i: 1             i: 2             i: 3             i: 0             i: 1             i: 2             i: 3
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8
        b: 0             b: 8             b: 12            b: 14            b: 0             b: 8             b: 12            b: 14
        color: 0         color: 0         color: 1         color: 3         color: 0         color: 0         color: 1         color: 3
        partner_rank: 14 partner_rank: 2  partner_rank: 4  partner_rank: 7  partner_rank: 15 partner_rank: 3  partner_rank: 5  partner_rank: 6
        _______________________________________________________________________________________________________________________________________________
        rank: 8          rank: 8          rank: 8          rank: 8          rank: 9          rank: 9          rank: 9          rank: 9
        i: 0             i: 1             i: 2             i: 3             i: 0             i: 1             i: 2             i: 3
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8
        b: 0             b: 8             b: 12            b: 14            b: 0             b: 8             b: 12            b: 14
        color: 0         color: 1         color: 2         color: 4         color: 0         color: 1         color: 2         color: 4
        partner_rank: 0  partner_rank: 12 partner_rank: 10 partner_rank: 9  partner_rank: 1  partner_rank: 13 partner_rank: 11 partner_rank: 8
        _______________________________________________________________________________________________________________________________________________
        rank: 10         rank: 10         rank: 10         rank: 10         rank: 11         rank: 11         rank: 11         rank: 11
        i: 0             i: 1             i: 2             i: 3             i: 0             i: 1             i: 2             i: 3
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8
        b: 0             b: 8             b: 12            b: 14            b: 0             b: 8             b: 12            b: 14
        color: 0         color: 1         color: 2         color: 5         color: 0         color: 1         color: 2         color: 5 
        partner_rank: 2  partner_rank: 14 partner_rank: 8  partner_rank: 11 partner_rank: 3  partner_rank: 15 partner_rank: 9  partner_rank: 10
        _______________________________________________________________________________________________________________________________________________
        rank: 12         rank: 12         rank: 12         rank: 12         rank: 13         rank: 13         rank: 13         rank: 13
        i: 0             i: 1             i: 2             i: 3             i: 0             i: 1             i: 2             i: 3
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8
        b: 0             b: 8             b: 12            b: 14            b: 0             b: 8             b: 12            b: 14
        color: 0         color: 1         color: 3         color: 6         color: 0         color: 1         color: 3         color: 6
        partner_rank: 4  partner_rank: 8  partner_rank: 14 partner_rank: 13 partner_rank: 5  partner_rank: 9  partner_rank: 15 partner_rank: 12
        _______________________________________________________________________________________________________________________________________________
        rank: 14         rank: 14         rank: 14         rank: 14         rank: 15         rank: 15         rank: 15         rank: 15
        i: 0             i: 1             i: 2             i: 3             i: 0             i: 1             i: 2             i: 3
        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8        pow2_i: 1        pow2_i: 2        pow2_i: 4        pow2_i: 8
        b: 0             b: 8             b: 12            b: 14            b: 0             b: 8             b: 12            b: 14
        color: 0         color: 1         color: 3         color: 7         color: 0         color: 1         color: 3         color: 7 
        partner_rank: 6  partner_rank: 10 partner_rank: 12 partner_rank: 15 partner_rank: 7  partner_rank: 11 partner_rank: 13 partner_rank: 14
        _______________________________________________________________________________________________________________________________________________
      */
      // POW2() defined in utils.h as: POW2(x) (1 << (x))
      pow2_i = POW2(i); // even (1 << 0) == 1, so no division by zero. Possible values are: 1,2,4,8
      b  = nproc - (int) (nproc / pow2_i); 
      color = (int)((rank & b) / POW2(loops - i )); // Bitwise AND
      partner_rank = rank ^ (int)(nproc/POW2(i + 1)); // Bitwise exclusive OR

      MPI_Comm new_comm;
      /*
        Creates new communicators based on colors and keys.
        int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
        color - Control of subset assignment (nonnegative integer).
        key   - Control of rank assignment (integer).
      */
      MPI_Comm_split(MPI_COMM_WORLD, color, rank, &new_comm);
      MPI_Comm_rank(new_comm, &sub_rank);

      powColor = pow2_i + color; // this computation is done multiple times. Just do it once, and reuse it.

      if(sub_rank == 0) {
        d = 0;
        for(j = 1; j < dims; j++) {
          // if the delta of J > delta of D, assign d to j
          if((nodes_gbox[powColor][j].upper - nodes_gbox[powColor][j].lower) > (nodes_gbox[powColor][d].upper - nodes_gbox[powColor][d].lower))
            d = j;
        }
      } 
      /*
        Broadcasts a message from the process with rank root to all other processes of the group.
        int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
      */
      MPI_Bcast(&d, 1, MPI_INT, 0, new_comm);

      // compute the median in this dimension `d`, for the new communication system
      float median  = get_median(dbs, d, new_comm);   
      // dtermine what points, and how many; `s_count`; for the current node to send to its `partner_rank` node, within the new communication system
      s_count = get_points_to_send(dbs, send_buf, invalid_pos_as, median, d, rank, partner_rank);

      int recvcount;
      int sendcount = s_count * dims;

      if (rank < partner_rank) {
        /*
          Sends and receives a message.
          int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, 
                           int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, 
                           int source, int recvtag, MPI_Comm comm, MPI_Status *status)
        */
        MPI_Sendrecv(&s_count, 1, MPI_INT, partner_rank, 4, &r_count, 1, MPI_INT, partner_rank, 5, MPI_COMM_WORLD, &status);
        recvcount = r_count * dims;
        recv_buf.resize(recvcount, 0.0);
        MPI_Sendrecv(&send_buf[0], sendcount, MPI_FLOAT, partner_rank, 2, &recv_buf[0], recvcount, MPI_FLOAT, partner_rank, 3, MPI_COMM_WORLD, &status);
        send_buf.clear();
      } else {
        MPI_Sendrecv(&s_count, 1, MPI_INT, partner_rank, 5, &r_count, 1, MPI_INT, partner_rank, 4, MPI_COMM_WORLD, &status);
        recvcount = r_count * dims;
        recv_buf.resize(recvcount, 0.0);
        MPI_Sendrecv(&send_buf[0], sendcount, MPI_FLOAT, partner_rank, 3, &recv_buf[0], recvcount, MPI_FLOAT, partner_rank, 2, MPI_COMM_WORLD, &status);
        send_buf.clear();
      }

      update_points(dbs, s_count, invalid_pos_as, recv_buf);
      recv_buf.clear(); // free the memory
      
      // LOWER() is defined in utils.h as: LOWER(i) (i<<1)
      copy_box(dbs, nodes_gbox[LOWER(powColor)], nodes_gbox[powColor]);
      nodes_gbox[LOWER(pow2_i+color)][d].upper =  median;
      // UPPER() is defined in utils.h as: UPPER(i) ((i<<1)+1)
      copy_box(dbs, nodes_gbox[UPPER(powColor)], nodes_gbox[powColor]);
      nodes_gbox[UPPER(powColor)][d].lower =  median; 
      /*
        Done with the new communication system, so need to get rid of it...
        Mark a communicator object for deallocation.
        int MPI_Comm_free(MPI_Comm *comm)
      */
      MPI_Comm_free(&new_comm);
    }

    // free the allocated memory
    for(i = 0; i < nproc; i++)
      delete [] nodes_gbox[i];

    delete [] nodes_gbox;
    delete [] gbox;
    delete [] box;

  }
  // TODO need to investigate this further, using `proc_of_interest` for cout() statements
  void update_points(ClusteringAlgo& dbs, int s_count, vector <int>& invalid_pos_as, vector <float>& recv_buf) {
    int i, j, k, l, r_count = recv_buf.size() / dbs.m_pts->m_i_dims;
    int newPtCt = dbs.m_pts->m_i_num_points + r_count - s_count; // this is used in multiple places, but none of the values change

    if(r_count >= s_count) {
      
      invalid_pos_as.resize(newPtCt, 1);

      //allocate memory for the points
      dbs.m_pts->m_points.resize(newPtCt);
      for(int ll = 0; ll < newPtCt; ll++)
        dbs.m_pts->m_points[ll].resize(dbs.m_pts->m_i_dims);

      j = 0;
      for(i = 0; i < invalid_pos_as.size(); i++) {
        if(invalid_pos_as[i] == 1) {
          for(k = 0; k < dbs.m_pts->m_i_dims; k++)
            dbs.m_pts->m_points[i][k] = recv_buf[j++];
        }
      }     

      dbs.m_pts->m_i_num_points = newPtCt;
    } else {
      j = 0;
      i = 0;  
      if(recv_buf.size() > 0) {
        for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
          if(invalid_pos_as[i] == 1) {
            for(k = 0; k < dbs.m_pts->m_i_dims; k++)
              dbs.m_pts->m_points[i][k] = recv_buf[j++];
          
            if(j == recv_buf.size()) {
              i++;
              break;
            }
          }
        }
      }
      
      l = dbs.m_pts->m_i_num_points;
      for( ; i < invalid_pos_as.size(); i++) {
        if(invalid_pos_as[i] == 1) {
          while(l > i) {
            l--;
            if(invalid_pos_as[l] == 0)
              break;
          }

          if(invalid_pos_as[l] == 0)  
            for(k = 0; k < dbs.m_pts->m_i_dims; k++)
              dbs.m_pts->m_points[i][k] = dbs.m_pts->m_points[l][k];
        }
      }

      //allocate memory for the points
      dbs.m_pts->m_points.resize(newPtCt);
      for(int ll = 0; ll < newPtCt; ll++)
        dbs.m_pts->m_points[ll].resize(dbs.m_pts->m_i_dims);

      dbs.m_pts->m_i_num_points = newPtCt;
    }   
  }

  /*
    Finds the points, and how many, for a node to send to its `partner_rank` node, 
    within the new communication system.
  */
  int get_points_to_send(ClusteringAlgo& dbs, vector <float>& send_buf, vector <int>& invalid_pos_as, float median, int d, int rank, int partner_rank) {
    int i, count = 0, j;
    bool lessRank = rank < partner_rank; // instead of calculating this for every itteration of the loop, just do it once.
    send_buf.reserve(dbs.m_pts->m_i_num_points * dbs.m_pts->m_i_dims);
    invalid_pos_as.clear();
    invalid_pos_as.resize(dbs.m_pts->m_i_num_points, 0);

    for(i = 0; i < dbs.m_pts->m_i_num_points; i++) {
      if (lessRank) {
        if(dbs.m_pts->m_points[i][d] > median) {
          invalid_pos_as[i] = 1;
          count++;
          for(j = 0; j < dbs.m_pts->m_i_dims; j++)
            send_buf.push_back(dbs.m_pts->m_points[i][j]);
        }
      } else {
        if(dbs.m_pts->m_points[i][d] <= median) {
          invalid_pos_as[i] = 1;
          count++;
          for(j = 0; j < dbs.m_pts->m_i_dims; j++)
            send_buf.push_back(dbs.m_pts->m_points[i][j]);
        }
      }
    }

    return count;
  } 

  /*
    Returns the median of a given attribute/dimension `d`, within the new
    communication system.
  */
  float get_median(ClusteringAlgo& dbs, int d, MPI_Comm& new_comm) { 

    float median;
    
    vector <float> data;
    data.reserve(dbs.m_pts->m_i_num_points);
    data.resize(dbs.m_pts->m_i_num_points, 0);
    // gather all values for a single attribute
    for (int k=0; k < dbs.m_pts->m_i_num_points; k++)
      data[k] = dbs.m_pts->m_points[k][d];
    // Find Kth element without recusion
    median = findKMedian(data, data.size()/2);
    data.clear(); // clear some memory

    int proc_count;
    MPI_Comm_size(new_comm, &proc_count);

    vector <float> all_medians;
    all_medians.resize(proc_count, 0);
    /*
      Gathers data from all processes and distributes it to all processes.
      int MPI_Allgather(const void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf, 
                        int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
    */
    MPI_Allgather(&median, sizeof(int), MPI_BYTE, &all_medians[0], sizeof(int), MPI_BYTE, new_comm);  
    // Find Kth element without recusion
    median = findKMedian(all_medians, all_medians.size()/2); 
    all_medians.clear(); // clear some memory
    
    return median;  
  }

  /*
    Determines the upper & lower values of each dimension, in the current
    node. 
  */
  void compute_local_bounding_box(ClusteringAlgo& dbs, interval* box) {
    int i, j;

    //we assume each processor has at least one point
    for(i = 0; i < dbs.m_pts->m_i_dims; i++) {
      box[i].upper = dbs.m_pts->m_points[0][i];
      box[i].lower = dbs.m_pts->m_points[0][i];
    }
    // find the upper and lower bound for each dimension
    for(i = 0; i < dbs.m_pts->m_i_dims; i++) {
      for(j = 1; j < dbs.m_pts->m_i_num_points; j++) {
        if(box[i].lower > dbs.m_pts->m_points[j][i])
          box[i].lower = dbs.m_pts->m_points[j][i];
        else if(box[i].upper < dbs.m_pts->m_points[j][i])
          box[i].upper = dbs.m_pts->m_points[j][i];
      }
    }
  }

  /*
    Determines the upper & lower values of each dimension, for the entire
    system.
  */
  void compute_global_bounding_box(ClusteringAlgo& dbs, interval* box, interval* gbox, int nproc) {
    int i, j, k;
  
    interval* gather_local_box = new interval[dbs.m_pts->m_i_dims * nproc];
  
    /* 
      gather the local bounding box first.
      Gathers data from all processes and distributes it to all processes.
      int MPI_Allgather(const void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf, 
                        int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
    */
    int sndRcvCt = sizeof(interval) * dbs.m_pts->m_i_dims; // calculation is done a couple of times. Do it once, and use the variable.
    MPI_Allgather(box, sndRcvCt, MPI_BYTE, gather_local_box, sndRcvCt, MPI_BYTE, MPI_COMM_WORLD);

    // compute the global bounding box, for each dimension
    for(i = 0; i < dbs.m_pts->m_i_dims; i++) {
      gbox[i].lower = gather_local_box[i].lower;
      gbox[i].upper = gather_local_box[i].upper;
      
      k = i;
      for(j = 0; j < nproc; j++, k += dbs.m_pts->m_i_dims) {
        if(gbox[i].lower > gather_local_box[k].lower)
          gbox[i].lower = gather_local_box[k].lower;
        
        if(gbox[i].upper < gather_local_box[k].upper)
          gbox[i].upper = gather_local_box[k].upper;
      }
    }
    
    delete [] gather_local_box;
  }

  /*
    Copies the values of `gbox` into every 1st vector element of `nodes_gbox`.
    The 1st vector element of `nodes_gbox` is of size `internal_nodes`.
    The 2nd vector element of `nodes_gbox`, and the size of `gbox`, is the number of dimensions in every point.
  */
  void copy_global_box_to_each_node(ClusteringAlgo& dbs, interval** nodes_gbox, interval* gbox, int internal_nodes) {
    int i, j;
    for(i = 0; i < internal_nodes; i++) {
      for(j = 0; j < dbs.m_pts->m_i_dims; j++) {
        nodes_gbox[i][j].upper = gbox[j].upper;
        nodes_gbox[i][j].lower = gbox[j].lower;
      }
    }
  }
  
  /*
    Copies the values of one vector of structs, into another
  */
  void copy_box(ClusteringAlgo& dbs, interval* target_box, interval* source_box) {
    for(int j = 0; j < dbs.m_pts->m_i_dims; j++) {
      target_box[j].upper = source_box[j].upper;
      target_box[j].lower = source_box[j].lower;
    }
  }
};

