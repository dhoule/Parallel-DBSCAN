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

#include "dbscan.h"

namespace NWUClustering {

  /*
    eps = epsilon/radious
    minPts = minimum number of points need to make a cluster
    seed_percentage = percentage of points each node is to use
  */
  void ClusteringAlgo::set_dbscan_params(double eps, int minPts, double seed_percentage) {
    m_epsSquare =  eps * eps;
    m_minPts =  minPts;
    m_compression = 0; // can set to 1 if want to compress specailly in the first round of communication
    m_perc_of_dataset = seed_percentage;
  }

  // Destructor
  ClusteringAlgo::~ClusteringAlgo() {
    m_parents.clear();
    m_parents_pr.clear();
    m_child_count.clear();
    m_corepoint.clear();
    m_member.clear();
  }

  /*
    Function keeps track of the time it takes to complete for some reason.
    called in run_dbscan_algo_uf_mpi_interleaved()
  */
  void ClusteringAlgo::trivial_decompression(vector <int>* data) {
    double start = MPI_Wtime();
    
    vector <int> parser;
    // allocates a MINIMUM amount of memory, the size of 'data'
    parser.reserve((*data).size());

    // declaring an INT automatically initializes it to 0... I believe the same goes for Vector objects as well
    int pid_count = parser[0], pos, i, j, pid, npid, npid_count;
    
    pos++;// the value becomes 1
    // TODO I have no idea what this does. Need to have outputs to determine what is going on;
      // besides rebuilding the 'data' variable
    while(pid_count > 0) {
      pid = parser[pos++];
      npid_count = parser[pos++];
      for(j = 0; j < npid_count; j++) {
        // push_back() adds elements to the "back" of the vector
        (*data).push_back(pid); 
        (*data).push_back(parser[pos++]);
      }

      pid_count--;
    }
    // removes all elements, destroying them. The size becomes 0.
    parser.clear();
  }

  // called in run_dbscan_algo_uf_mpi_interleaved()
  void ClusteringAlgo::trivial_compression(vector <int>* data, vector < vector <int> >* parser) {
    // get the starting time before doing anything in this function
    double start = MPI_Wtime();
    double org = 0;
    int pairs, pid, npid, i, j, pid_count, npid_count, data_size = (*data).size(), temp;

    // The number of "pairs" in the "data" vector
    pairs = data_size/2; 
    // TODO don't know what "org" is supposed to stand for, 'original' maybe
    org = data_size;   
    
    // loop over 'data' and add elements to the back of a vector 'parser'[pid]
    for(i = 0; i < pairs; i++) {
      temp = 2 * i;
      pid = (*data)[temp];
      npid = (*data)[temp + 1];
      (*parser)[pid].push_back(npid);
    }
    // empty the 'data' vector, and set the size to 0
    (*data).clear();
    // inititalize 'pid_count' to 0, and add it to the back of the 'data' vector
    pid_count = 0;
    (*data).push_back(pid_count); // uniques pids, should update later
    // 'm_pts' is the current cluster's struct object. Initialized in clusters.cpp read_file().
    // Loop rebuilds the 'data' vector and clears out dimensions of the 'parser' vector
    temp = m_pts->m_i_num_points;
    for(i = 0; i < temp; i++) {
      npid_count = (*parser)[i].size();
      if(npid_count > 0) {
        (*data).push_back(i);
        (*data).push_back(npid_count);
        for(j = 0; j < npid_count; j++) {
          (*data).push_back((*parser)[i][j]);       
        }
        pid_count++;
        (*parser)[i].clear();
      }
    }
  }

  // called in mpi_main.cpp.
    // The function merges Points from other nodes
  void ClusteringAlgo::get_clusters_distributed() {
    // Determine the current node's rank within the cluster, and the size of the cluster itself
    int rank, nproc, i, numPts, pid; // pid == process ID???
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    /* 
      get the root all the local points first
      store in a message buffer in case need to ask other processors
    */

    vector < vector <int > > merge_received;
    vector < vector <int > > merge_send1;
    vector < vector <int > > merge_send2;
    vector <int> init;
    // resize the vectors to the size of the number of nodes,
      //the new elements are initialized as copies of val; init; otherwise, they are value-initialized.
    merge_received.resize(nproc, init); 
    merge_send1.resize(nproc, init);
    merge_send2.resize(nproc, init); 
    init.clear(); // clear the unused variable ASAP

    // loop over the other dimensions of the vectors, 
      // resizing them to the minimum length of the number of points that 'm_pts' has.
    numPts = m_pts->m_i_num_points;
    for(pid = 0; pid < nproc; pid++) {
      merge_received[pid].reserve(numPts);
      merge_send1[pid].reserve(numPts);
      merge_send2[pid].reserve(numPts);
    }
    // These vectors, together, seem to be used to perform a bubble sort type of operation
      // 'p_cur_send' and 'p_cur_insert' have actual uses.
    vector < vector <int > >* pswap;
    vector < vector <int > >* p_cur_send;
    vector < vector <int > >* p_cur_insert;
    // assign the memory addresses of the variables???
    p_cur_send = &merge_send1; // TODO merge_send1 is only used for this purpose
    p_cur_insert = &merge_send2; // TODO merge_send2 is only used for this purpose
    // cannot clear `merge_send1` OR `merge_send2`, as the & gives the "address of". This seems stupid, as it just renames variables
    // `m_child_count` is declared in dbscan.h
    
    m_child_count.resize(numPts, 0);   
    
    int root, local_continue_to_run = 0, global_continue_to_run;
    
    /*
     loop over the points
    */
    for(i = 0; i < numPts; i++) {
      // find the point containing i
      root = i; // root is x in the paper
      
      // continue the loop as long as the node's ID matches the parent's ID???
      // TODO this while loop occures twice in this function.
      // TODO move it to another function, passing it the starting value of `root` and return the found value of `root`
      while(m_parents_pr[root] == rank) {
        if(m_parents[root] == root) { 
          // tree root is satisfied: p(x) == x
          break;
        }
        root = m_parents[root];
      }
      
      if(m_parents[root] == root && m_parents_pr[root] == rank) { // root is a local root
        // set the root of i directly to root
        m_parents[i] = root; // creating a new set for each element x is acheived by setting p(x) to x
        m_child_count[root] = m_child_count[root] + 1; // increase the child count by one
        //m_parents_pr[i] = rank; // NO NEED TO SET THIS AS IT // as it what??? Ans: it should already be in the current node
      } else {
        // set up info to request data from other nodes
        (*p_cur_insert)[m_parents_pr[root]].push_back(0); // flag: 0 means query and 1 means a reply
        (*p_cur_insert)[m_parents_pr[root]].push_back(m_parents[root]); // The point that belongs to another node
        (*p_cur_insert)[m_parents_pr[root]].push_back(i); // the element index of that point
        (*p_cur_insert)[m_parents_pr[root]].push_back(rank); // the MPI_SOURCE
        local_continue_to_run++; // increase msg counter
      }     
    }

    /*
      pos = used in MPI_Waitany(), for the `index` attribute
      quadraples = used to determine how many points have actually been sent to a node
      scount = used to determine the number of messages that have actually been sent
    */
    int pos, quadraples, scount, tid, tag = 0, tagPlusOne = 1, rtag, rsource, rcount, isend[nproc], irecv[nproc], flag;
    
    MPI_Request s_req_recv[nproc], s_req_send[nproc], d_req_send[nproc], d_req_recv[nproc]; 
    MPI_Status  d_stat_send[nproc], d_stat; // d_stat_send[nproc]: for MPI_Waitall; d_stat: for MPI_Waitany
    int target_point, source_point, source_pr;
    
    while(1) {
      global_continue_to_run = 0; // A flag to count the number of messages that still need to be passed
      // Combines values from all processes and distributes the result back to all processes
      // int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
      MPI_Allreduce(&local_continue_to_run, &global_continue_to_run, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      // All msg passing is done, break out of the loop
      if(global_continue_to_run == 0)
        break;
      // TODO bubble sort operation, that is somehow connected to Rem's union-find technique for Kruskal's Alg
      pswap = p_cur_insert;
      p_cur_insert = p_cur_send;
      p_cur_send = pswap;
      // wipe out the current vectors of `p_cur_insert`, because there are more msgs that need to be created
      for(tid = 0; tid < nproc; tid++)
        (*p_cur_insert)[tid].clear();
      // The number of msgs that were sent 
      scount = 0;
      // loop over the messages that need to be sent, and use a non-blocking operation
      for(tid = 0; tid < nproc; tid++) {
        isend[tid] = (*p_cur_send)[tid].size();
        // if the intended node, has at least 1 msg to be sent to it
        if(isend[tid] > 0) {
          /* 
            Starts a standard-mode, nonblocking send
            int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
          */  
          MPI_Isend(&(*p_cur_send)[tid][0], isend[tid], MPI_INT, tid, tagPlusOne, MPI_COMM_WORLD, &d_req_send[scount]);
          scount++;
        }
      }
      /*
        All nodes exchange data about how many messages they sent to each other.
        Each node calls MPI_Send() & MPI_Recv() for every node in the system, including itself.
        All processes send data to all processes, blocking call.
        int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
      */
      MPI_Alltoall(&isend[0], 1, MPI_INT, &irecv[0], 1, MPI_INT, MPI_COMM_WORLD);
      
      rcount = 0;
      for(tid = 0; tid < nproc; tid++) {
        if(irecv[tid] > 0) {
          merge_received[tid].clear();
          merge_received[tid].assign(irecv[tid], -1); // resize the vector and assign the values to -1
          /*
            Starts a standard-mode, nonblocking receive.
            int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
          */
          MPI_Irecv(&merge_received[tid][0], irecv[tid], MPI_INT, tid, tagPlusOne, MPI_COMM_WORLD, &d_req_recv[rcount]);
          rcount++;
        }
      }
      // reset the counter for the number of sent local msgs
      local_continue_to_run = 0;
      // loop over the received messages
      for(tid = 0; tid < rcount; tid++) { 
        // Waits for any specified send or receive to complete. 
        // int MPI_Waitany(int count, MPI_Request array_of_requests[], int *index, MPI_Status *status)
        MPI_Waitany(rcount, &d_req_recv[0], &pos, &d_stat);
        // Need to get the tag and the source of the msg
        rtag = d_stat.MPI_TAG;
        rsource = d_stat.MPI_SOURCE;
        // determine if the tag received is the same as the tag sent
        if(rtag == tagPlusOne) {
          quadraples = merge_received[rsource].size()/4; // divided by 4 because of how `(*p_cur_insert)` is built
          // `pid` is just a simple counter and not used for anything else, in this scope
          for(pid = 0; pid < quadraples; pid++) {
            /* 
              get the quadraple
              back() returns the last element of the vector
              pop_back() removes the last element of the vector
              Need to get the last value of the vector, then discard it to move onto the next
            */
            source_pr = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            source_point = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            target_point = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            flag = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            
            if(flag == 0) { // flag: 0 means query and 1 means a reply 
              root = target_point;
              // Need to move up the minimum spanning tree, to find the current root pt
              // TODO this while loop occures twice in this function.
              // TODO move it to another function, passing it the starting value of `root` and return the found value of `root`
              while(m_parents_pr[root] == rank) {
                // when the root pt is found, break out of the loop
                if(m_parents[root] == root)
                  break;
                root = m_parents[root];
              }
              if(m_parents[root] == root && m_parents_pr[root] == rank) { // root is a local root
                // have to return the child pt to the node it belongs to
                m_child_count[root] = m_child_count[root] + 1; // increase the child count by one
                
                (*p_cur_insert)[source_pr].push_back(1); // flag: 0 means query and 1 means a reply
                (*p_cur_insert)[source_pr].push_back(source_point); // the element index of that point
                (*p_cur_insert)[source_pr].push_back(m_parents[root]); // The root of the tree
                (*p_cur_insert)[source_pr].push_back(m_parents_pr[root]); // The node the point belongs to
                local_continue_to_run++;
              } else {
                // Still need to request more info about the pt from the node
                (*p_cur_insert)[m_parents_pr[root]].push_back(0); // flag: 0 means query and 1 means a reply
                (*p_cur_insert)[m_parents_pr[root]].push_back(m_parents[root]); // The root of the tree
                (*p_cur_insert)[m_parents_pr[root]].push_back(source_point); // the element index of that point
                (*p_cur_insert)[m_parents_pr[root]].push_back(source_pr); // The node the point belongs to
                local_continue_to_run++;
              }
            } else {
              // got a reply, so just set the parent
              m_parents[target_point] = source_point; // The root of the minimum spanning tree
              m_parents_pr[target_point] = source_pr; // The node the point belongs to
            }
          }
        }
      }

      tag++; // increment `tag` for the next msg set
      tagPlusOne++; // increment the next one to keep it in sync
      if(scount > 0) //  TODO is this block acting the same as MPI_Barrier(MPI_COMM_WORLD); ???
        MPI_Waitall(scount, &d_req_send[0], &d_stat_send[0]); // wait for all the sending operation
    }

    /* 
      `final_cluster_root` = number of clusters in the current node
      `total_final_cluster_root` = total number of clusters
      `points_in_cluster_final` = points that are in clusters, in the current node
      `total_points_in_cluster_final` = total number of points in clusters
    */
    int final_cluster_root = 0, total_final_cluster_root = 0;
    int points_in_cluster_final = 0, total_points_in_cluster_final = 0;
    // determine the number of points in the current cluster
    for(i = 0; i < numPts; i++) {
      if(m_parents[i] == i && m_parents_pr[i] == rank && m_child_count[i] > 1) {
        points_in_cluster_final += m_child_count[i];
        final_cluster_root++;
      }
    } 
    /*
      These MPI_Allreduce() are used to share metadata between all nodes.
      Combines values from all processes and distributes the result back to all processes.
      int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
    */
    // Determine the total number of pts that are in clusters, and not noise pts, in the entire system
    MPI_Allreduce(&points_in_cluster_final, &total_points_in_cluster_final, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // Determine the number of clusters in the entire system
    MPI_Allreduce(&final_cluster_root, &total_final_cluster_root, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // Determine the total number of points in entire system
    int total_points = 0;
    MPI_Allreduce(&m_pts->m_i_num_points, &total_points, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if(rank == proc_of_interest) cout << "Points in clusters " << total_points_in_cluster_final << " Noise " << (total_points - total_points_in_cluster_final) << " Total points " << total_points << endl;
    if(rank == proc_of_interest) cout << "Total number of clusters " << total_final_cluster_root << endl;
    
    vector<int> global_roots; // used to determine the roots in each node
    global_roots.resize(nproc, 0); // set the size and initialize to 0
    /*
      Need to determine the roots that are in other nodes.
      Gathers data from all processes and distributes it to all processes. 
      int MPI_Allgather(const void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf[starting address], int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
    */
    MPI_Allgather(&final_cluster_root, sizeof(int), MPI_BYTE, &global_roots[0], sizeof(int), MPI_BYTE, MPI_COMM_WORLD); 
    
    /*
      Determine the range of available IDs for each node by determining the largest number of roots on a single node
      `cluster_offset` = used to assign cluster IDs
    */
    int cluster_offset = 0;
    // every node only goes so far as its `rank` in the system
    for(i = 0; i <= rank; i++) 
      cluster_offset += global_roots[i];
    // `cluster_offset` becomes the MAX value for cluster IDs on a node
    // declared in cluster.h: vector <int>  m_pid_to_cid;
    m_pid_to_cid.clear(); // point ID to cluster ID. 1st use of this var in this function.
    m_pid_to_cid.resize(numPts, -1);
    // assign for the global roots only
    for(i = 0; i < numPts; i++) {
      // if current tree node is the current itterable AND the current pointer is for the current node
      if(m_parents[i] == i && m_parents_pr[i] == rank) {
        // if the number of points is greater than 1
        if(m_child_count[i] > 1) { 
          m_pid_to_cid[i] = cluster_offset;
          cluster_offset--;
        } else {
          m_pid_to_cid[i] = 0; // noise point
        }
      }
    }
    for(i = 0; i < numPts; i++) {
      if(m_parents_pr[i] == rank) {
        if(m_parents[i] != i) { //skip the noise points
          m_pid_to_cid[i] = m_pid_to_cid[m_parents[i]]; // assign local roots
        }
      } else {
        // ask the outer to to send back the clusterID
        (*p_cur_insert)[m_parents_pr[i]].push_back(0); // flag: 0 means query and 1 means a reply
        (*p_cur_insert)[m_parents_pr[i]].push_back(m_parents[i]);
        (*p_cur_insert)[m_parents_pr[i]].push_back(i);
        (*p_cur_insert)[m_parents_pr[i]].push_back(rank);
      }
    }
    // MAY BE REMOVED
    //MPI_Barrier(MPI_COMM_WORLD);

    /*
      TODO The following section begins to merge the info from other nodes
      TODO Why is the following for-loop only ran twice?
    */

    tag++;
    tagPlusOne++;
    int later_count;
    for(later_count = 0; later_count < 2; later_count++) {
      pswap = p_cur_insert;
      p_cur_insert = p_cur_send;
      p_cur_send = pswap;

      for(tid = 0; tid < nproc; tid++)
        (*p_cur_insert)[tid].clear();
  
      scount = 0;
      
      for(tid = 0; tid < nproc; tid++) {
        isend[tid] = (*p_cur_send)[tid].size();
        if(isend[tid] > 0) {
          // send the outer points to their corresponding nodes???
          // int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
          MPI_Isend(&(*p_cur_send)[tid][0], isend[tid], MPI_INT, tid, tagPlusOne, MPI_COMM_WORLD, &d_req_send[scount]);
          scount++;
        }
      }
      // All processes send data to all processes.
      // int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
      MPI_Alltoall(&isend[0], 1, MPI_INT, &irecv[0], 1, MPI_INT, MPI_COMM_WORLD);
      
      rcount = 0;
      for(tid = 0; tid < nproc; tid++) {
        if(irecv[tid] > 0) {
          merge_received[tid].clear();
          merge_received[tid].assign(irecv[tid], -1);
          // This is the corresponding Irec() to the Isend above???
          // int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
          MPI_Irecv(&merge_received[tid][0], irecv[tid], MPI_INT, tid, tagPlusOne, MPI_COMM_WORLD, &d_req_recv[rcount]);
          rcount++;
        }
      }
      for(tid = 0; tid < rcount; tid++) {
        MPI_Waitany(rcount, &d_req_recv[0], &pos, &d_stat);
        
        rtag = d_stat.MPI_TAG;
        rsource = d_stat.MPI_SOURCE;

        if(rtag == tagPlusOne) {
          quadraples = merge_received[rsource].size()/4; // divided by 4 because of how `(*p_cur_insert)` is built
          // break up `merge_received[rsource]`, to rebuild `(*p_cur_insert)[source_pr]`,
            // depending on the `flag` value, that is changed in the rebuild process
          for(pid = 0; pid < quadraples; pid++) {
            // get the quadraple
            source_pr = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            source_point = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            target_point = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            flag = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            if(flag == 0) {         
              (*p_cur_insert)[source_pr].push_back(1);
              (*p_cur_insert)[source_pr].push_back(source_point);
              (*p_cur_insert)[source_pr].push_back(m_pid_to_cid[m_parents[target_point]]);
              (*p_cur_insert)[source_pr].push_back(-1); // One extra INT, may be needed in future
            } else {
              // got a reply, so just set the parent
              m_pid_to_cid[target_point] = source_point; // this assigns the clusterID
            }
          }
        }
      }

      if(scount > 0)
        MPI_Waitall(scount, &d_req_send[0], &d_stat_send[0]); // wait for all the sending operation
    
      //MPI_Barrier(MPI_COMM_WORLD); // MAY NEED TO ACTIVATE THIS
      tag++;
      tagPlusOne++;
    }

    merge_received.clear();
    merge_send1.clear();
    merge_send2.clear();
    global_roots.clear();
  }

  // called in mpi_main.cpp.
  // writes data to file
  void ClusteringAlgo::writeCluster_distributed(string outfilename) {
    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    int i;

    // get the total number of points
    int total_points = 0;
    MPI_Allreduce(&m_pts->m_i_num_points, &total_points, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    vector<int> point_count;
    point_count.resize(nproc, 0);
    MPI_Allgather(&m_pts->m_i_num_points, sizeof(int), MPI_BYTE, &point_count[0], sizeof(int), MPI_BYTE, MPI_COMM_WORLD);

    int ret, ncfile;
    
    string outfilename_dis = outfilename;
    outfilename_dis = outfilename_dis; //"_clusters.nc";
  
    // create the file, if exists, open the file.
    ret = ncmpi_create(MPI_COMM_WORLD, outfilename_dis.c_str(), NC_CLOBBER, MPI_INFO_NULL, &ncfile);  
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    } 

    MPI_Offset num_particles = total_points;
    int num_particles_id;

    ret = ncmpi_def_dim(ncfile, "num_particles", num_particles, &num_particles_id);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }

    string column_name_initial = "position_col_X";
    stringstream column_name_id;
    string column_name;

    int j, ncolumn = 1, col_id = 0, varid[m_pts->m_i_dims + 1]; //number of column is 1, col_id is 0 as we use the first one

    // write the columns
    int temp = m_pts->m_i_dims;
    for(j = 0; j < temp; j++) {
      column_name_id.str("");
      column_name_id << j;
  
      column_name = column_name_initial + column_name_id.str(); 
          
      ret = ncmpi_def_var(ncfile, column_name.c_str(), NC_FLOAT, ncolumn, &col_id, &varid[j]);          
      if (ret != NC_NOERR) {
        handle_error(ret, __LINE__);
        return;
      }        
    }

    column_name = "cluster_id";

    ret = ncmpi_def_var(ncfile, column_name.c_str(), NC_INT, ncolumn, &col_id, &varid[j]);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }

    ret = ncmpi_enddef(ncfile);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }

    MPI_Offset start[2], count[2];
      
    start[0] = 0;
    for(i = 0; i < rank; i++) {
      start[0] += point_count[i]; 
    }

    count[0] = point_count[rank];
    start[1] = 0; // this to satisfy PnetCDF requirement
    count[1] = 1;//dim_sizes[dimids[1]];
    
    // allocate memory  
    float *data = new float[count[0] * count[1]];

    // write the data columns
    temp = m_pts->m_i_dims;
    int temp2 = m_pts->m_i_num_points;
    for(j = 0; j < temp; j++) {
      // get the partial column data
      for(i = 0; i < temp2; i++)
        data[i] = m_pts->m_points[i][j];

      // write the data
      ret = ncmpi_put_vara_float_all(ncfile, varid[j], start, count, data);
      if (ret != NC_NOERR) {
        handle_error(ret, __LINE__);
        return;
      }
    }
    delete [] data;

    int *data_id = new int[count[0] * count[1]];    

    //write the cluster_ids
    for(i = 0; i < temp2; i++)
      data_id[i] = m_pid_to_cid[i];

    ret = ncmpi_put_vara_int_all(ncfile, varid[m_pts->m_i_dims], start, count, data_id);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }

    delete [] data_id;  

    // close the file
    ret = ncmpi_close(ncfile);
    if (ret != NC_NOERR) {
      handle_error(ret, __LINE__);
      return;
    }
  }

  // "uf" == "Union Find"
  // called in mpi_main.cpp
    // Function gets the union of 2 tress, to create a larger cluster
  void run_dbscan_algo_uf_mpi_interleaved(ClusteringAlgo& dbs) {
    double start = MPI_Wtime();     
    int i, pid, j, k, npid;
    int rank, nproc;
    kdtree2_result_vector ne;
    kdtree2_result_vector ne_outer;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // initialize some parameters
    int numPts = dbs.m_pts->m_i_num_points; 
    // assign parent to itestf
    dbs.m_parents.resize(numPts, -1);
    dbs.m_parents_pr.resize(numPts, -1);

    int total_points = 0, points_per_pr[nproc], start_pos[nproc];

    // getting the total number of local points and assigning postions  
    MPI_Allgather(&dbs.m_pts->m_i_num_points, 1, MPI_INT, &points_per_pr[0], 1, MPI_INT, MPI_COMM_WORLD);
    
    for(i = 0; i < nproc; i++) {
      start_pos[i] = total_points;
      total_points += points_per_pr[i];
    }

    // assign proc IDs
    vector <int> vec_prID;
    vec_prID.resize(total_points, -1);

    k = 0;
    for(i = 0; i < nproc; i++) {
      for(j = 0; j < points_per_pr[i]; j++)
        vec_prID[k++] = i;
    }
    
    // resetting the membership and corepoints values
    dbs.m_member.resize(numPts, 0);
    dbs.m_corepoint.resize(numPts, 0);
    // returns the starting address of the `the_data`
    vector<int>* ind = dbs.m_kdtree->getIndex();
    vector<int>* ind_outer = dbs.m_kdtree_outer->getIndex();

    // setting parents to itself and corresponding proc IDs
    for(i = 0; i < numPts; i++) { // TODO this is the looping over all of the points
      pid = (*ind)[i]; // TODO need to change how the Point ID is retrieved
      dbs.m_parents[pid] = pid;
      dbs.m_parents_pr[pid] = rank;
    }

    vector < vector <int > > merge_received;
    vector < vector <int > > merge_send1;
    vector < vector <int > > merge_send2;
    vector <int> init;
    int rtag, rsource, tag = 0, tagPlusOne = 1, pos = 0, scount, rcount, isend[nproc], irecv[nproc];
    
    merge_received.resize(nproc, init);
    merge_send1.resize(nproc, init);
    merge_send2.resize(nproc, init);
    
    // reserving communication buffer memory
    int tempPlusNodes = numPts * nproc;
    for(pid = 0; pid < nproc; pid++) {
      merge_received[pid].reserve(tempPlusNodes);
      merge_send1[pid].reserve(tempPlusNodes);
      merge_send2[pid].reserve(tempPlusNodes);
    }

    int root, root1, root2, tid;

    vector < vector <int > >* pswap;
    vector < vector <int > >* p_cur_send;
    vector < vector <int > >* p_cur_insert;

    p_cur_send = &merge_send1;
    p_cur_insert = &merge_send2;

    if(rank == proc_of_interest) cout << "Init time " << MPI_Wtime() - start << endl; 

    MPI_Barrier(MPI_COMM_WORLD);
    
    // the main part of the DBSCAN algorithm (called local computation)
    start = MPI_Wtime();
    for(i = 0; i < numPts; i++) { // TODO this is the looping over all of the points
      pid = (*ind)[i]; // TODO need to change how the Point ID is retrieved
      // getting the local neighborhoods of local point
      ne.clear();
      dbs.m_kdtree->r_nearest_around_point(pid, 0, dbs.m_epsSquare, ne);
      
      ne_outer.clear();
      vector<float> qv(dbs.m_pts->m_i_dims);
      // `qv` stands for Query Vector. It is a vector of the current point's dimensions.
      for (int u = 0; u < dbs.m_pts->m_i_dims; u++) {
        qv[u] = dbs.m_kdtree->the_data[pid][u];
      }

      // getting the remote neighborhood of the local point
      if(dbs.m_pts_outer->m_i_num_points > 0)
        dbs.m_kdtree_outer->r_nearest(qv, dbs.m_epsSquare, ne_outer);
    
      qv.clear();
      int ne_outer_size = ne_outer.size();
      if(ne.size() + ne_outer_size >= dbs.m_minPts) {
        // pid is a core point
        root = pid;
        dbs.m_corepoint[pid] = 1;
        dbs.m_member[pid] = 1;
        
        // traverse the remote neighbors and add in the communication buffers  
        for(j = 0; j < ne_outer_size; j++) {
          npid = ne_outer[j].idx;
          int outer_parentIds = dbs.m_pts_outer->m_prIDs[npid];
          (*p_cur_insert)[outer_parentIds].push_back(pid);
          (*p_cur_insert)[outer_parentIds].push_back(dbs.m_pts_outer->m_ind[npid]);
        }
        
        //traverse the local neighbors and perform union operation
        for (j = 0; j < ne.size(); j++) {
          npid = ne[j].idx;
          
          // get the root containing npid
          root1 = npid;
          root2 = root;
          if(dbs.m_corepoint[npid] == 1 || dbs.m_member[npid] == 0) {
            dbs.m_member[npid] = 1;

            // REMS algorithm to (union) merge the trees
            while(dbs.m_parents[root1] != dbs.m_parents[root2]) {
              if(dbs.m_parents[root1] < dbs.m_parents[root2]) {
                if(dbs.m_parents[root1] == root1) {
                  dbs.m_parents[root1] = dbs.m_parents[root2];
                  root = dbs.m_parents[root2];
                  break;
                }

                // splicing comression technique
                int z = dbs.m_parents[root1];
                dbs.m_parents[root1] = dbs.m_parents[root2];
                root1 = z;
              } else {
                if(dbs.m_parents[root2] == root2) {
                  dbs.m_parents[root2] = dbs.m_parents[root1];
                  root = dbs.m_parents[root1];
                  break;
                }

                // splicing compressio technique
                int z = dbs.m_parents[root2];
                dbs.m_parents[root2] = dbs.m_parents[root1];                  
                root2 = z;
              }
            }
          }
        }
      }
    }
      
    MPI_Barrier(MPI_COMM_WORLD);

    int v1, v2, par_proc, triples, local_count, global_count;
    double temp_inter_med, inter_med;

    if(rank == proc_of_interest) cout << "Local computation took " << MPI_Wtime() - start << endl;

    inter_med = MPI_Wtime();

    i = 0;
      
    MPI_Request s_req_recv[nproc], s_req_send[nproc], d_req_send[nproc], d_req_recv[nproc]; 
    MPI_Status  d_stat_send[nproc], d_stat;

    start = MPI_Wtime();

    local_count = 0;

    // performing additional compression for the local points that are being sent 
    // this steps identifies the points that actually going to connect the trees in other processors
    // this step will eventually helps further compression before the actual communication happens
    for(tid = 0; tid < nproc; tid++) {
      triples = (*p_cur_insert)[tid].size()/2;
      local_count += triples;
      for(pid = 0; pid < triples; pid++) {
        int twoPid = 2 * pid;
        v1 = (*p_cur_insert)[tid][twoPid];
        root1 = v1;
        while(dbs.m_parents[root1] != root1) {
          root1 = dbs.m_parents[root1];
        }

        while(dbs.m_parents[v1] != root1) {
          int tmp = dbs.m_parents[v1];
          dbs.m_parents[v1] = root1;
          v1 = tmp;
        }

        (*p_cur_insert)[tid][twoPid] = root1;
      }
    }

    local_count = local_count/nproc;  
    global_count = 0;
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //message_per_round
    int uv, uf, um, ul, ucount;
    int local_continue_to_run, global_continue_to_run;

    vector <vector <int> > parser;
    vector <int> init_ex;
    parser.resize(numPts, init_ex);
    
    while(1) {
      pswap = p_cur_insert;
      p_cur_insert = p_cur_send;  
      p_cur_send = pswap;
      for(tid = 0; tid < nproc; tid++)
        (*p_cur_insert)[tid].clear();
    
      scount = 0;
      for(tid = 0; tid < nproc; tid++) {
        if(dbs.m_compression == 1 && i == 0 && (*p_cur_send)[tid].size() > 0) {
          dbs.trivial_compression(&(*p_cur_send)[tid], &parser);
        }

        isend[tid] = (*p_cur_send)[tid].size();
        if(isend[tid] > 0) {
          MPI_Isend(&(*p_cur_send)[tid][0], isend[tid], MPI_INT, tid, tagPlusOne, MPI_COMM_WORLD, &d_req_send[scount]);
          scount++;
        }
      }
      
      MPI_Alltoall(&isend[0], 1, MPI_INT, &irecv[0], 1, MPI_INT, MPI_COMM_WORLD);

      rcount = 0;
      for(tid = 0; tid < nproc; tid++) {
        if(irecv[tid] > 0) {
          merge_received[tid].clear();
          merge_received[tid].assign(irecv[tid], -1);
          MPI_Irecv(&merge_received[tid][0], irecv[tid], MPI_INT, tid, tagPlusOne, MPI_COMM_WORLD, &d_req_recv[rcount]);
          rcount++;
        }
      }

      local_count = 0; 
      //get the data and process them
      for(tid = 0; tid < rcount; tid++) {
        MPI_Waitany(rcount, &d_req_recv[0], &pos, &d_stat);
        rtag = d_stat.MPI_TAG;
        rsource = d_stat.MPI_SOURCE;  
        if(rtag == tagPlusOne) {
          // process received the data now
          if(i == 0) {
            if(dbs.m_compression == 1) {
              // call the decompression function
              dbs.trivial_decompression(&merge_received[rsource]);

              triples = merge_received[rsource].size()/2;
              par_proc = rsource;
            } else {
              triples = merge_received[rsource].size()/2;
              par_proc = rsource;
            }
          } else {
            triples = merge_received[rsource].size()/3;
          }

          for(pid = 0; pid < triples; pid++) {
            // get the pair
            v1 = merge_received[rsource].back();
            merge_received[rsource].pop_back();
            
            if(i > 0) {
              par_proc = merge_received[rsource].back();
              merge_received[rsource].pop_back();
            }
        
            v2 = merge_received[rsource].back();
            merge_received[rsource].pop_back();
        
            int con = 0;
            if(i > 0)
              con = 1;
            else if (i == 0 && (dbs.m_corepoint[v1] == 1 || dbs.m_member[v1] == 0)) { 
              dbs.m_member[v1] = 1;
              con = 1;
            }

            if(con == 1) {
        
              root1 = v1;
              // this will find the boundary vertex or the root if the root is in this processor
              while(dbs.m_parents_pr[root1] == rank) {
                if(dbs.m_parents[root1] == root1)
                  break;

                root1 = dbs.m_parents[root1];
              }
        
              // compress the local path
              while(v1 != root1 && vec_prID[v1] == rank) {
                int tmp = dbs.m_parents[v1];
                dbs.m_parents[v1] = root1;
                v1 = tmp;
              }

              if(dbs.m_parents[root1] == v2 && dbs.m_parents_pr[root1] == par_proc) {
                continue;
              }           
                
              if(par_proc == rank) {
                if(dbs.m_parents[root1] == dbs.m_parents[v2])
                  continue;
              }
                
              if(dbs.m_parents[root1] == root1 && dbs.m_parents_pr[root1] == rank) { // root1 is a local root
                if(start_pos[rank] + root1 < start_pos[par_proc] + v2) {
                  // do union
                  dbs.m_parents[root1] = v2;
                  dbs.m_parents_pr[root1] = par_proc;
                  continue;
                } else {
                  // ask the parent of v2
                  (*p_cur_insert)[par_proc].push_back(root1);
                  (*p_cur_insert)[par_proc].push_back(dbs.m_parents_pr[root1]);
                  (*p_cur_insert)[par_proc].push_back(v2);
                  local_count++;
                }
              } else {
                // root1 is not local
                if(start_pos[dbs.m_parents_pr[root1]] + root1 < start_pos[par_proc] + v2) {
                  // ask the parent of root1
                  (*p_cur_insert)[dbs.m_parents_pr[root1]].push_back(v2);
                  (*p_cur_insert)[dbs.m_parents_pr[root1]].push_back(par_proc);
                  (*p_cur_insert)[dbs.m_parents_pr[root1]].push_back(dbs.m_parents[root1]);
                  local_count++;
                } else {
                  // ask the parent of v2
                  (*p_cur_insert)[par_proc].push_back(dbs.m_parents[root1]);
                  (*p_cur_insert)[par_proc].push_back(dbs.m_parents_pr[root1]);
                  (*p_cur_insert)[par_proc].push_back(v2);
                  local_count++;
                }
              }
            }
          }
          
          merge_received[rsource].clear();
        } else {
          cout << "rank " << rank << " SOMETHING IS WRONG" << endl;
        }
      }
        
      if(scount > 0)
        MPI_Waitall(scount, &d_req_send[0], &d_stat_send[0]); 

      tag += 2; // change the tag value although not important
      tagPlusOne = tag + 1;
      
      local_continue_to_run = 0;
      local_count = 0;
      for(tid = 0; tid < nproc; tid++) {
        local_count += (*p_cur_insert)[tid].size()/3; 
        if((*p_cur_insert)[tid].size() > 0)
          local_continue_to_run = 1;
      }
      
      local_count = local_count / nproc;
      
      global_count = 0;
      global_continue_to_run = 0;
            
      MPI_Allreduce(&local_continue_to_run, &global_continue_to_run, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      if(global_continue_to_run == 0)
        break;
      i++;
    }

    if(rank == proc_of_interest) cout << "Merging took " << MPI_Wtime() - start << endl;

    pswap = NULL;
    p_cur_insert = NULL;
    p_cur_send = NULL;

    for(tid = 0; tid < nproc; tid++) {
      merge_received[tid].clear();
      merge_send1[tid].clear();
      merge_send2[tid].clear();
    }

    merge_received.clear();
    merge_send1.clear();
    merge_send2.clear();
    ind = NULL;
    ind_outer = NULL;

    vec_prID.clear();
    ne.clear();
    ne_outer.clear();
    parser.clear();
    init_ex.clear();
    init.clear();
  }
};

