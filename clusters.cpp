
#include "clusters.h"

namespace NWUClustering {

  // destructor 
    // can have 4 parameters: m_pts,m_kdtree,m_pts_outer,m_kdtree_outer
    // Clears, and deletes Points, Points_Outer, & kdtree2 objects
  Clusters::~Clusters() { // Pure virtual destructor 
    if(m_pts) {
      m_pts->m_points.clear();
      delete m_pts;
      m_pts = NULL;
    }

    if(m_kdtree) {
      delete m_kdtree;
      m_kdtree = NULL;
    }

    if(m_pts_outer) {
      m_pts_outer->m_prIDs.clear();
      m_pts_outer->m_ind.clear();
      m_pts_outer->m_points.clear();
      delete m_pts_outer;
      m_pts_outer = NULL;
    }

    if(m_kdtree_outer) {
      delete m_kdtree_outer;
      m_kdtree_outer = NULL;
    }
  }

  // Sets up the Points_Outer objects
    // Called in geometric_partitioning.cpp
  bool Clusters::allocate_outer(int dims) {
    if(m_pts_outer == NULL) {
      m_pts_outer = new Points_Outer;
      m_pts_outer->m_prIDs.clear();
      m_pts_outer->m_ind.clear();
      m_pts_outer->m_i_dims = dims;
      m_pts_outer->m_i_num_points = 0;
    }
  }

  // Adds points to a cluster object
    // Called in geometric_partitioning.cpp
  bool Clusters::addPoints(int source, int buf_size, int dims, vector<float>& raw_data) {

    // used as itterables
    int i, j, k; 
    // pos = original num of points, then used as a key in 'm_pts_outer->m_points'
    int pos;
    // num_points = new number of points for the cluster
    int num_points = buf_size / dims;
    
    // incorrect dimension
      // the number of dimensions/attributes can never change
    if(m_pts_outer->m_i_dims != dims)
      return false;

    // Save the original number of points, then increase by the new "num_points" value
    pos = m_pts_outer->m_i_num_points;
    m_pts_outer->m_i_num_points += num_points;
    
    //allocate memory for the points
      // resize() - Resizes the container so that it contains n elements
    m_pts_outer->m_points.resize(m_pts_outer->m_i_num_points);
    // resize each new point object
    int temp = m_pts_outer->m_i_num_points;
    for(int ll = 0; ll < temp; ll++) 
      m_pts_outer->m_points[ll].resize(dims);
    // resize the point IDs
    m_pts_outer->m_prIDs.resize(m_pts_outer->m_i_num_points, -1);

    k = 0;
    // loop over the new number of points
    for(i = 0; i < num_points; i++) {
      // loop over the number of dimensions/attributes
      for(j = 0; j < dims; j++) {
        // assign the new 'raw_data' elements to the newest points
        m_pts_outer->m_points[pos][j] = raw_data[k++]; // TODO possible buffer overflow????????
        
      }
      // assign the Node ID for the point
      m_pts_outer->m_prIDs[pos] = source; // `source` is the Node that the point is on
      // increment the counter for the key...
      pos++;
    }

    return true;
  }

  // Updates OUTER points' cluster IDs
    // Called in geometric_partitioning.cpp
  bool Clusters::updatePoints(vector< vector<int> >& raw_ind) {
    
    // used as itterables
    int i, j = -1;
    // I believe these are cluster IDs
    int source = -1,  prev_source = -1;
    // resize 'm_ind', by 'm_i_num_points' spaces, new elements are initialized as copies of -1
    m_pts_outer->m_ind.resize(m_pts_outer->m_i_num_points, -1);
    // loop over the Outer points, and update cluster IDs
    int temp = m_pts_outer->m_i_num_points;
    for(i = 0; i < temp; i++) {
      source = m_pts_outer->m_prIDs[i];
      
      if(source != prev_source)
        j = 0;

      m_pts_outer->m_ind[i] = raw_ind[source][j++];
      prev_source = source;
    }

    return true;
  }

  /*
    Called from mpi_main.cpp
    Reads the binary file. 
    Input data points are equally partioned, each core reading their corresponding part of the file. 
    Later the points will be partioned geometrically.
    Points are created here, based off of the data points.
    According to the README, 'num_points' & 'dims' HAVE to be the 1st 2 things in the file (each 4 bytes).
  */
  int Clusters::read_file(char* infilename, int isBinaryFile) {
    
    // used as itterables
    int i, j;
    // rank = current node's ID, nproc = total number of nodes in the system
    int rank, nproc;
    int num_points, dims;
    // Get the current node's 'rank'
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get the total number of nodes in the system
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    // only supports binary data file
    if(isBinaryFile == 1) {
      // NOTE: the input data points are equally partioned and each core read their corresponding part
      // later the points will be partioned geometrically

      ifstream file (infilename, ios::in|ios::binary);
      if(file.is_open()) {
        // size of the dataset
        file.read((char*)&num_points, sizeof(int));
        // number of dimensions for each data point/record
        file.read((char*)&dims, sizeof(int));

        // compute the respective segments of the file to read.
          // sch = section or search?
          // lower = the starting position in the file to start reading from
          // upper = the last position in the file to read from
        long long sch, lower, upper;
        // determine the size of the section of file to read
        if(num_points % nproc == 0)
          sch = num_points/nproc; // the 2 are evenly divisible
        else
          sch = num_points/nproc + 1; // taking care of the remainder, by adding a whole 1


        lower = sch * rank; // the pointer for the starting location
        upper = sch * (rank + 1); // the pointer for the end location
        // if the "end location" exceeds the number of points, reduce it to the number of points
        if(upper > num_points)
          upper = num_points;
        
        // allocate memory for points
        m_pts = new Points;
        m_pts->m_i_dims = dims; // number of dimensions/attributes
        m_pts->m_i_num_points = upper - lower; // the number of points
        // interval struct defined in kdtree2.hpp
          // m_box, is of type interval*
        m_pts->m_box = new interval[m_pts->m_i_dims]; 
        
        //allocate memory for the points
        m_pts->m_points.resize(m_pts->m_i_num_points);
        int temp = m_pts->m_i_num_points;
        for(int ll = 0; ll < temp; ll++)
          m_pts->m_points[ll].resize(dims);

        // point_coord_type = typedef float point_coord_type; in 'utils.h'     
        // initializes 'pt' variable
        point_coord_type* pt = new point_coord_type[dims * sizeof(point_coord_type)];

        // fseek to the respective position of the file
        file.seekg(lower * dims * sizeof(point_coord_type), ios::cur);
        // loop over the area of the file for the node
        int delta = upper - lower;
        for (i = 0; i < delta; i++) {
          // signature: istream& read (char* s, streamsize n);
          // Extracts n characters from the stream and stores them in the array pointed to by s.
          file.read((char*)pt, dims * sizeof(point_coord_type));
          
          for (j = 0; j < dims; j++) {
            m_pts->m_points[i][j] = pt[j];
            temp = pt[j];
            if(i == 0) { 
              m_pts->m_box[j].upper = temp;
              m_pts->m_box[j].lower = temp;
            } else {
              if(m_pts->m_box[j].lower > temp)
                m_pts->m_box[j].lower = temp;
              else if(m_pts->m_box[j].upper < temp)
                m_pts->m_box[j].upper = temp;
            }
          }
        }

        free(pt);
        pt = NULL;

        file.close();
      } else {
        cout << "rank " << rank << " Error: no such file: " << infilename << endl;
        return -1;
      }
    } else {
      cout << "Only supports binary data: Failed to read data" << endl;
      return -1;
    }
    return 0;
  }
  
  /*
   Called from mpi_main...
   Builds the kd-tree containing the "main" points for the node.
  */
  int Clusters::build_kdtree() {

    // if the Point objects weren't created, don't bother going further...
    if(m_pts == NULL) {
      cout << "Point set is empty" << endl;
      return -1;
    }

    //m_kdtree = new kdtree2(m_pts->m_points, true);
    m_kdtree = new kdtree2(m_pts->m_points, false);

    if(m_kdtree == NULL) {
      cout << "Falied to allocate new kd tree for orginial points" << endl;
      return -1;
    }

    return 0;   
  } 
  /*
   Called from mpi_main...
   Builds the kd-tree containing the outer points for the node.
  */
  int Clusters::build_kdtree_outer() {
    
    if(m_pts_outer == NULL) {
      cout << "Outer point set is empty" << endl;
      return -1;
    }

    if(m_pts_outer->m_i_num_points > 0) {

      //m_kdtree_outer = new kdtree2(m_pts->m_points_outer, true);
      m_kdtree_outer = new kdtree2(m_pts_outer->m_points, false);

      if(m_kdtree_outer == NULL) {
        cout << "Falied to allocate new kd tree for outer points" << endl;
        return -1;
      }
    }

    return 0;   
  } 
}
