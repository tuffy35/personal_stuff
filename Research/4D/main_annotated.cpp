#include <iostream> // write to the standard input/output streams
#include <iomanip> // manipulate the format of input and output of program
#include <fstream> // ability to read from and write to files
#include <cmath>  // basic math stuff
#include <array> // allows more easy-to-work-with arrays and array methods
#include <vector> // allows for easier work with vectors

struct neighbor_data { // structure for neighbor data, containing the number of neighbors, and a vector holding the identities of the neighbors
   int number_of_neighbors;
   std::vector<int> neighbors;
};

double three_determinant(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3) // function that computes a 3x3 determinant
{
    double fval = 0;

    fval = a1*(b2*c3-c2*b3) -a2*(b1*c3-c1*b3) + a3*(b1*c2-c1*b2);

    return fval;
}


double exact_solution(double x, double y, double z, double w) // function that exponentiates the sum of x,y,z,w
{
    //double fval = x + 2.*y + 3.*z + 4.*w;
    double fval = exp(x + y + z + w);

    return fval;
}

double forcing_term(double x, double y, double z, double w) // function that is just 4*exact_solution
{
    double fval = 4.0*exp(x + y + z + w);

    return fval;
}


bool four_by_four_inverse(const double m[16], double invOut[16]) // m[16] is a square array 16 elements (4x4) and invOut is the same size but now the inverse
{
    double inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12]; // computing the determinant of the 4x4 matrix

    if (det == 0) // if the determinant is 0, the matrix isn't invertible, so we just return false
    {
        return false;
    }

    det = 1.0 / det;

    for (i = 0; i < 16; i++) // actually computing the inverse now that we know it's invertible
    {
        invOut[i] = inv[i] * det;
    }

    return true;
}

// the & after the vector definitions tell us that we're manipulating the original vector, NOT a copy
void residual_computation(int n_points, int n_edges, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& w, 
                          std::vector<int>& number_of_neighbors_1D, std::vector<int>& offsets_1D, std::vector<int>& neighbors_1D,  
                          std::vector< std::vector< std::array<double,4> > >& coefficient_vector, std::vector< std::array<int,2> >& global_edges, 
                          std::vector< std::array<double,4> >& n01, std::vector<double>& dual_hypervolumes,
                          std::vector<double>& U, std::vector<double>& S, std::vector<double>& residual)
{
    for(int i = 0; i<n_points; ++i)
    {
        residual[i] = 0.; // initialize the residual to be 0 at each node
    }
    
    //*********************
    // Compute Gradients
    //*********************

    std::vector< std::array<double,4> > gradient; // vector of 4 valued arrays called gradient, which has as many arrays as there are points
    gradient.reserve(n_points);

    for(int i = 0; i<n_points; ++i) // loop through the points
    {
        for(int j=0; j<4; ++j)
        {
            gradient[i][j] = 0.; // make each value in each array 0 to start
        }

        int number_of_neighbors_local = number_of_neighbors_1D[i]; // number of neighbors of the i'th node
        int offset_local = offsets_1D[i]; // the offset for this node so we can locate the identities of these neighbors

        for(int j=0; j<number_of_neighbors_local; ++j) // loop through the neighbors for this node
        {
            int index_neighbor = neighbors_1D[offset_local+(j+1)]; // the current neighbor that we're dealing with
            for(int k=0; k<4; ++k) // looping through the x,y,z,w
            {
                gradient[i][k] += coefficient_vector[i][j][k]*(U[index_neighbor] - U[i]);  
            } // the k'th direction in gradient is updated with the coefficient_vector value at that point times the difference in solution value between the neighbor and this node
        }
    }
    
    /*
    for(int i = 0; i<n_points; ++i)
    {
        std::cout << "Gradient info " << std::endl;
        for(int j=0; j<4; ++j)
        {
            std::cout << " " << gradient[i][j];
        }
        std::cout << " " << std::endl;
    }
    */

    //************************
    // Edge Loop on Residuals
    //************************

    for(int i = 0; i<n_edges; ++i) // loop through the edges
    {
        int n0, n1;
        n0 = global_edges[i][0]; //n0 is the lower node, n01 is the upper node
        n1 = global_edges[i][1];

        double U0, U1;
        U0 = U[n0]; //U0 is the solution at n0, and U1 is the solution at n1
        U1 = U[n1];

        std::vector<double> n01_unit(4);
        
        double mag = 0.;

        for(int j=0; j<4; ++j)
        {
            mag += n01[i][j]*n01[i][j];
        }
        mag = sqrt(mag);

        for(int j=0; j<4; ++j)
        {
            n01_unit[j] = n01[i][j]/mag; // the unit vector of n01 for this edge
        }
        
        //XXX Needs to be modified if PDE changes
        std::vector<double> advection_direction = {1., 1., 1., 1.};

        double projected_direction = 0.;
        for(int j=0; j<4; ++j)
        {
           projected_direction += advection_direction[j]*n01_unit[j]; // dotting the advection direction with the unit normal vector
        }

        double flux_temp;
        //flux_temp = 0.5*projected_direction*(U0+U1)*mag;
        double UL, UR;
        
        UL = U0 + gradient[n0][0]*(x[n1]-x[n0])*0.5 //UL is the gradient at n0 times the difference between n1 and n0, divided by 2
                + gradient[n0][1]*(y[n1]-y[n0])*0.5
                + gradient[n0][2]*(z[n1]-z[n0])*0.5
                + gradient[n0][3]*(w[n1]-w[n0])*0.5;

        UR = U1 - gradient[n1][0]*(x[n1]-x[n0])*0.5
                - gradient[n1][1]*(y[n1]-y[n0])*0.5
                - gradient[n1][2]*(z[n1]-z[n0])*0.5
                - gradient[n1][3]*(w[n1]-w[n0])*0.5;

        flux_temp = (0.5*projected_direction*(UL + UR) - 0.5*abs(projected_direction)*(UR - UL))*mag; // calculation for the flux temp THE SECOND TERM IS THE UPWINDING TERM

        residual[n0] += flux_temp;
        residual[n1] -= flux_temp; // add that stuff to the residuals
    }

    //************************
    // Node Loop (forcing term)
    //************************

    for(int i=0; i<n_points; ++i)
    {
       residual[i] -= S[i]*dual_hypervolumes[i]; // each residual is updated with the forcing term at that point times the hypervolume of that point
    }

    return;
}




int main(void)
{
        //*******************
        //Read in coordinates
        //*******************

        std::ifstream coord_file;
        coord_file.open("coordinates_mesh_4.txt", std::ios_base::in);

        double a;
        int total_entries = 0;
        while(coord_file >> a) // reads one item at a time (NOT lines at a time) and returns true if there's any data there. if there isn't, then it reached the end and the loop stops
        {
            total_entries++;
        }
        coord_file.close();
        
        int n_points = total_entries/4; // it's divided by 4 because total_entries has 4 numbers per point (4 coordinates to define each)
        std::cout << "N points " << n_points << std::endl;
        // establishing vectors x,y,z,w and then making sure they have as many spots as there are points
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        std::vector<double> w;

        x.reserve(n_points);
        y.reserve(n_points);
        z.reserve(n_points);
        w.reserve(n_points);

        coord_file.open("coordinates_mesh_4.txt", std::ios_base::in);

        for(int i=0; i<n_points; i++) // writing into the vectors x,y,z,w. i corresponds to the line, and within each line, x 1st, y 2nd, z 3rd, w 4th
        {
            coord_file >> x[i];
            coord_file >> y[i];
            coord_file >> z[i];
            coord_file >> w[i];
        }
        coord_file.close();

        //********************
        //Read in connectivity
        //********************

        std::ifstream cell_file;
        cell_file.open("cells_mesh_4.txt", std::ios_base::in);

        int b;
        int total_indexes = 0;
        while(cell_file >> b)
        {
            total_indexes++;
        }
        cell_file.close();
        
        int n_pentatopes = total_indexes/5; // divide by 5 because each line has 5 numbers (5 nodes bound each pentatope)
        std::cout << "N cells " << n_pentatopes << std::endl;

        std::vector< std::array<int,5> > pent; // defining a vector called pent which contains elements which are arrays with 5 integer components (the nodes that bound each pent)
        pent.reserve(n_pentatopes); // size is: number of pentatopes x number of vertices (5)

        cell_file.open("cells_mesh_4.txt", std::ios_base::in);

        for(int i=0; i<n_pentatopes; i++) // filling in the pent vector with the correct values
        {
            cell_file >> pent[i][0];
            cell_file >> pent[i][1];
            cell_file >> pent[i][2];
            cell_file >> pent[i][3];
            cell_file >> pent[i][4];
        }
        cell_file.close();
 
        //Assumes base 1 indexing
        for(int i=0; i<n_pentatopes; i++) // this is assuming that the txt files we're reading from starts counting at 1, so we decrement everything by 1 to make it consistent with c++ (start at 0)
        {
            pent[i][0]--;
            pent[i][1]--;
            pent[i][2]--;
            pent[i][3]--;
            pent[i][4]--;
        }

        //pent, x, y, z, w

        //*********************
        //Create neighbor data
        //*********************

        int max_neighbors = 100; // is this arbitrary?

        std::vector< std::array<int,4> > local_neighbor = { {1, 2, 3, 4}, {2, 0, 3, 4}, {0, 1, 3, 4}, {0, 1, 4, 2}, {0, 1, 2, 3} }; // the LOCAL neighbors are the same for everything (local neighbors for node 0, then 1, then 2, then 3, then 4)
        std::vector< std::array<int,2> > local_edge = { {0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4} }; // the LOCAL edges are the same for everything

        std::vector<neighbor_data> total_neighbor_data; // creating a vector called total neighbor data that contains data type: neighbor_data (has # of neighbors and vector of those neighbors) 
        total_neighbor_data.reserve(n_points); // each point has a corresponding "total neighbor data" which contains this information -----------------------------------^

        for(int i = 0; i<n_points; ++i)
        {
            total_neighbor_data[i].number_of_neighbors = 0; // initializing each node to have 0 neighbors
            total_neighbor_data[i].neighbors.reserve(max_neighbors); // initializing each node to have an available list of max_neighbors = 100 neighbors
        }

        //Loop over pentatopes
        for(int i = 0; i<n_pentatopes; ++i)
        {
           //Loop over vertices
           for(int k = 0; k<5; ++k)
           {
               int n1 = pent[i][k]; // define n1 to be the node of interest that we're currently working with. it's the k'th vertex on the i'th pentatope

               //Loop over other nodes
               for(int j = 0; j<4; ++j)
               {
                  int n2 = pent[i][local_neighbor[k][j]]; // define n2 as the node on the i'th pentatope which is a neighbor to n1. this works because local_neighbor is defined such that the local neighbor will exclude the current node

                  //If no neighbors found yet
                  if( total_neighbor_data[n1].number_of_neighbors == 0) //if the node of interest n1 has no neighbors documented at this point, make the num of neighbors 1, and add node 2 to the list of neighbors
                  {
                      total_neighbor_data[n1].number_of_neighbors = 1;
                      total_neighbor_data[n1].neighbors[0] = n2;
                  }
                  else
                  {
                      //Has neighbor been added? // AKA: has this neighbor already been added to the list? if yes, then skip cuz we good
                      bool neighbor_found = false;
                      for(int m = 0; m<total_neighbor_data[n1].number_of_neighbors; ++m)
                      {
                          if(n2 == total_neighbor_data[n1].neighbors[m])
                          {
                              neighbor_found = true;
                              break;
                          }
                      }

                      //New neighbor found
                      if(!neighbor_found) // if the neighbor hadn't been identified yet,
                      {
                          total_neighbor_data[n1].number_of_neighbors++; // increment the list
                          total_neighbor_data[n1].neighbors[total_neighbor_data[n1].number_of_neighbors-1] = n2; // in the list of n1's neighbors, add n2 in the spot corresponding to the number of neighbors we have (we subtract 1 because of incrementing starting at 0)
                      }
                  }
               }
           }
        }


        //*********************************
        //Store neighbor data in 1D arrays
        //*********************************

        int total_number_of_neighbors = 0; // simply contains the total number of neighbors

        for(int i=0; i < n_points; ++i)
        {
           total_number_of_neighbors += total_neighbor_data[i].number_of_neighbors;
        }

        std::vector<int> neighbors_1D(total_number_of_neighbors); // a long 1D array that lists out, in order, the neighbors to node 0, then those for node 1, etc.
        std::vector<int> offsets_1D(n_points); // the i'th entry in offsets_1D tells us how far into neighbors_1D we go to see the neighbors for node i
        std::vector<int> number_of_neighbors_1D(n_points); // how many neighbors node i has

        int itemp = 0;
        offsets_1D[0] = -1;

        for(int i=0; i<n_points; ++i)
        {
            number_of_neighbors_1D[i] = total_neighbor_data[i].number_of_neighbors; // number_of_neighbors_1D has n_points elements. element i is how many neighbors the i'th node has
            
            for(int k = 0; k<number_of_neighbors_1D[i]; ++k) // loop from 0 to the number of neighbors that the i'th node has
            {
                itemp = offsets_1D[i] + (k+1); // itemp starts at 0, and we scribe neighbors_1D at 0 with node i's k'th neighbor. within this loop, itemp just increases by 1 each time
                neighbors_1D[itemp] = total_neighbor_data[i].neighbors[k];
            }

            if(i< (n_points-1))
            {
                offsets_1D[i+1] = itemp; // setting the next value for offsets_1D. this happens at the end of the loop so that on the next iteration, i+1 will be the new i. this allows us to track where each node's neighbors starts in this huge array
            }
        }

        //*********************************
        //Edge numbering
        //*********************************
        int  n_edges = 0;

        for(int i = 0; i<n_points; ++i)
        {
            for(int k = 0; k < number_of_neighbors_1D[i]; ++k)
            {
                 if(i < neighbors_1D[offsets_1D[i]+(k+1)])
                 {
                     n_edges++; // increment the number of edges if node i is lower number than its neighbor node k
                 }
             }
        }
 
        //global edge
        std::vector< std::array<int,2> > global_edges; // vector called global edges. each element is an 2 element array corresponding to each edge. one array for each edge
        global_edges.reserve(n_edges);

        std::vector< std::array<int,100> > global_edge_numbers; // vector of arrays with 100 spots reserved for each. one array for each node
        global_edge_numbers.reserve(n_points);

        int edge_count = -1;
        for(int i = 0; i<n_points; ++i) // loop through each node
        {
            for(int k = 0; k < number_of_neighbors_1D[i]; ++k) // loop through the neighbors of each of these nodes
            {
                 if(i < neighbors_1D[offsets_1D[i]+(k+1)]) // same if statement as before - identify an edge
                 {
                     edge_count++; // increment the edge count
                     global_edges[edge_count][0] = i;
                     global_edges[edge_count][1] = neighbors_1D[offsets_1D[i]+(k+1)]; // add the two nodes (i and the neighbor corresponding with hk)

                     global_edge_numbers[i][k] = edge_count; // add the edge number to global edge numbers. 
                 }
            }
        }   


        //*********************************
        //Cell to Edge
        //*********************************
        std::vector< std::array<int,10> > cell2edge; // vector that has an array for each pentatope
        cell2edge.reserve(n_pentatopes);


        for(int i = 0; i<n_pentatopes; ++i) // lets loop through all the pentatopes
        {
            int index_edge;
            int n1, n2;
                            // manually going through the 10 edges on the pent
            //edge 0-1
            index_edge = 0;
            
            if( pent[i][0] < pent[i][1]) // for edge connecting local node # 0 and 1, if global node # for 0 is less than that of 1,
            {
                n1 = pent[i][0]; // assign n1 as the 0th node and n2 as the 1st node
                n2 = pent[i][1];
            }
            else // otherwise, do it the other way around
            {
                n1 = pent[i][1];
                n2 = pent[i][0];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2) // search through the list of neighbors_1D for this node 2
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k]; // identify the edge that connects these 2 and add it to cell2edge
                } // we need to do this little loop to identify k because that's how we defined global_edge_numbers. also this allows us to define cell2edge in accordance with local edge numbers
            }

            //edge 0-2
            index_edge = 1;
            
            if( pent[i][0] < pent[i][2])
            {
                n1 = pent[i][0];
                n2 = pent[i][2];
            }
            else
            {
                n1 = pent[i][2];
                n2 = pent[i][0];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }

            //edge 0-3
            index_edge = 2;
            
            if( pent[i][0] < pent[i][3])
            {
                n1 = pent[i][0];
                n2 = pent[i][3];
            }
            else
            {
                n1 = pent[i][3];
                n2 = pent[i][0];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }

            //edge 0-4
            index_edge = 3;
            
            if( pent[i][0] < pent[i][4])
            {
                n1 = pent[i][0];
                n2 = pent[i][4];
            }
            else
            {
                n1 = pent[i][4];
                n2 = pent[i][0];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }

            //edge 1-2
            index_edge = 4;
            
            if( pent[i][1] < pent[i][2])
            {
                n1 = pent[i][1];
                n2 = pent[i][2];
            }
            else
            {
                n1 = pent[i][2];
                n2 = pent[i][1];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }

            //edge 1-3
            index_edge = 5;
            
            if( pent[i][1] < pent[i][3])
            {
                n1 = pent[i][1];
                n2 = pent[i][3];
            }
            else
            {
                n1 = pent[i][3];
                n2 = pent[i][1];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }

            //edge 1-4
            index_edge = 6;
            
            if( pent[i][1] < pent[i][4])
            {
                n1 = pent[i][1];
                n2 = pent[i][4];
            }
            else
            {
                n1 = pent[i][4];
                n2 = pent[i][1];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }

            //edge 2-3
            index_edge = 7;
            
            if( pent[i][2] < pent[i][3])
            {
                n1 = pent[i][2];
                n2 = pent[i][3];
            }
            else
            {
                n1 = pent[i][3];
                n2 = pent[i][2];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }

            //edge 2-4
            index_edge = 8;
            
            if( pent[i][2] < pent[i][4])
            {
                n1 = pent[i][2];
                n2 = pent[i][4];
            }
            else
            {
                n1 = pent[i][4];
                n2 = pent[i][2];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }

            //edge 3-4
            index_edge = 9;
            
            if( pent[i][3] < pent[i][4])
            {
                n1 = pent[i][3];
                n2 = pent[i][4];
            }
            else
            {
                n1 = pent[i][4];
                n2 = pent[i][3];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
            }
        }        


        //**************************
        //Directed hyperarea vectors
        //**************************

        std::vector< std::array<double,4> > n01; // n01 is a vector for the directed hyperarea vectors which correspond with each edge
        n01.reserve(n_edges);

        for(int i = 0; i<n_edges; ++i) // initialize them all to be zero
        {
            n01[i][0] = 0.;
            n01[i][1] = 0.;
            n01[i][2] = 0.;
            n01[i][3] = 0.;
        }

        for(int i = 0; i<n_pentatopes; ++i) // loop through all the pentatopes
        {
            for(int j=0; j<10; ++j) // loop through the 10 edges on each pentatope
            {
                int m;

                int index_edge = cell2edge[i][j]; // call index edge the current edge we're looking at based on pent i and local edge j
                if( global_edges[index_edge][0] == pent[i][local_edge[j][0]] ) // if the first node on index edge (based on global edges) matches the first node of i'th pent (based on local_edge)
                {
                    m = local_edge[j][0]; // call this first node node m
                }
                else
                {
                    m = local_edge[j][1]; // otherwise, call the other node on this local edge node m
                }

                //Pentatope facet opposite to the node m
                int v0, v1, v2, v3; // identify the other 4 points on this pentatope that AREN'T node m. this is why we use the local_neighbor construct
                v0 = pent[i][local_neighbor[m][0]];
                v1 = pent[i][local_neighbor[m][1]];
                v2 = pent[i][local_neighbor[m][2]];
                v3 = pent[i][local_neighbor[m][3]];
                // defining variables for the x, y, z, and w coordinates for these 4 other points on the pent that aren't m.
                double ax, ay, az, aw;
                double bx, by, bz, bw;
                double cx, cy, cz, cw;
                double dx, dy, dz, dw;
                // we're classifying the points as a, b, c, and d
                ax = x[v0];
                bx = x[v1];
                cx = x[v2];
                dx = x[v3];

                ay = y[v0];
                by = y[v1];
                cy = y[v2];
                dy = y[v3];

                az = z[v0];
                bz = z[v1];
                cz = z[v2];
                dz = z[v3];

                aw = w[v0];
                bw = w[v1];
                cw = w[v2];
                dw = w[v3];

                double bax, bay, baz, baw;
                double cax, cay, caz, caw;
                double dax, day, daz, daw;
                // defining vectors that connect a to b, a to c, and a to d
                bax = bx - ax;
                bay = by - ay;
                baz = bz - az;
                baw = bw - aw;

                cax = cx - ax;
                cay = cy - ay;
                caz = cz - az;
                caw = cw - aw;

                dax = dx - ax;
                day = dy - ay;
                daz = dz - az;
                daw = dw - aw;
                // determinant of these vectors, which gives us the normal vector to the facet opposite of node m
                //Normal vector components
                double i_comp =  three_determinant(bay, baz, baw, cay, caz, caw, day, daz, daw)*(1./6.);
                double j_comp = -three_determinant(bax, baz, baw, cax, caz, caw, dax, daz, daw)*(1./6.);
                double k_comp = three_determinant(bax, bay, baw, cax, cay, caw, dax, day, daw)*(1./6.);
                double l_comp = -three_determinant(bax, bay, baz, cax, cay, caz, dax, day, daz)*(1./6.);

                //std::cout << "Normal components " << i_comp << " " << j_comp << " " << k_comp << " " << l_comp << std::endl;
                // e0 and e1 bound the index edge of interest
                //Edge vector
                int e0, e1;
                e0 = global_edges[index_edge][0];
                e1 = global_edges[index_edge][1];
                // we dot this with the normal vectors we calculated
                double edge_normal_dot = i_comp*(x[e1] - x[e0]) + j_comp*(y[e1] - y[e0]) 
                                       + k_comp*(z[e1] - z[e0]) + l_comp*(w[e1] - w[e0]);

                // if it's positive, our vector's oriented properly
                if(edge_normal_dot > 0)
                {
                    n01[index_edge][0] += i_comp;
                    n01[index_edge][1] += j_comp;
                    n01[index_edge][2] += k_comp;
                    n01[index_edge][3] += l_comp;
                } // if not, we have to flip it
                else
                {
                    n01[index_edge][0] -= i_comp;
                    n01[index_edge][1] -= j_comp;
                    n01[index_edge][2] -= k_comp;
                    n01[index_edge][3] -= l_comp;
                }
            }
        }

        for(int i = 0; i<n_edges; ++i) // doing the proper correction at the end for each of our directed hyperarea vectors
        {
            n01[i][0] *= 1./10.; 
            n01[i][1] *= 1./10.;
            n01[i][2] *= 1./10.;
            n01[i][3] *= 1./10.;

            //std::cout << "Directed area vectors " << i << ": " << n01[i][0] << " " << n01[i][1] << " " << n01[i][2] << " " << n01[i][3] << std::endl;
        }

        //*********************************
        //Check directed hyperarea vectors
        //*********************************

        std::vector< std::array<double, 4> > vector_residual; // vector residual is like the aj vector i had for 2D
        vector_residual.reserve(n_points);

        for(int i = 0; i<n_points; ++i)
        {
            vector_residual[i][0] = 0.;
            vector_residual[i][1] = 0.;
            vector_residual[i][2] = 0.;
            vector_residual[i][3] = 0.;
        }

        std::cout << "N points " << n_points << std::endl;
        for(int i = 0; i<n_edges; ++i) // loop through the edges
        {
            int e0 = global_edges[i][0]; // identify e0 as the lower node and e1 as the higher node
            int e1 = global_edges[i][1];

            for(int j=0; j<4; ++j)
            {
                vector_residual[e0][j] += n01[i][j]; // add the vector to e0 and subtract it from e1 because n01 points from e0 to e1
                vector_residual[e1][j] -= n01[i][j];
            }
        }


        //************************************
        //Identify Boundary facets
        //************************************

        std::vector< std::array<int,4> > boundary_facets; // array called boundary_facets that has 2*n_pentatopes elements. this corresponds to boundary edges in 2D
        boundary_facets.reserve(2*n_pentatopes);

        std::vector<int> boundary_facet_2_cell; // the nth entry in this array is the pent number that corresponds to the nth boundary facet
        boundary_facet_2_cell.reserve(2*n_pentatopes);

        int boundary_facet_counter = 0;

        for(int i=0; i<n_pentatopes; ++i) // loop through all of the pentatopes
        {
            int ia, ib, ic, id, ie;
            // call ia-ie the 5 nodes that bound this current pent i
            ia = pent[i][0];
            ib = pent[i][1];
            ic = pent[i][2];
            id = pent[i][3];
            ie = pent[i][4];
            // call these the 5 facets that bound (equivelant to edges in 2D or faces in 3D)
            bool facet_abcd_found = false;
            bool facet_bcde_found = false;
            bool facet_acde_found = false;
            bool facet_abde_found = false;
            bool facet_abce_found = false;

            for(int j=0; j<n_pentatopes; ++j) // loop again through all the pentatopes, but now under a different index j
            {
                if(i!=j) // if statement is here so we can see each pent's relationship to each other pent, not itself
                {
                    int aflag, bflag, cflag, dflag, eflag;
                    // the flags tell us which of the nodes on i match the nodes on j
                    aflag = 0;
                    bflag = 0;
                    cflag = 0;
                    dflag = 0;
                    eflag = 0;
                    // the same is ia-ie but now on pent j
                    int ja, jb, jc, jd, je;
                    ja = pent[j][0];
                    jb = pent[j][1];
                    jc = pent[j][2];
                    jd = pent[j][3];
                    je = pent[j][4];
                    
                    if((ia==ja) || (ia==jb) || (ia==jc) || (ia==jd) || (ia==je)) // if the 1st node on i matches any of the nodes on j
                    {
                        aflag = 1;
                    }
                    if((ib==ja) || (ib==jb) || (ib==jc) || (ib==jd) || (ib==je)) // if the 2nd node on i matches any of the nodes on j
                    {
                        bflag = 1;
                    }
                    if((ic==ja) || (ic==jb) || (ic==jc) || (ic==jd) || (ic==je)) // etc. etc. etc.
                    {
                        cflag = 1;
                    }
                    if((id==ja) || (id==jb) || (id==jc) || (id==jd) || (id==je))
                    {
                        dflag = 1;
                    }
                    if((ie==ja) || (ie==jb) || (ie==jc) || (ie==jd) || (ie==je))
                    {
                        eflag = 1;
                    }
                    // check EVERY OTHER PENT to see if the facet is found anywhere in the entire mesh.
                    //Facet found
                    if((aflag+bflag+cflag+dflag)==4)
                    {
                        facet_abcd_found = true;   
                    }
                    if((bflag+cflag+dflag+eflag)==4)
                    {
                        facet_bcde_found = true;   
                    }
                    if((aflag+cflag+dflag+eflag)==4)
                    {
                        facet_acde_found = true;   
                    }
                    if((aflag+bflag+dflag+eflag)==4)
                    {
                        facet_abde_found = true;   
                    }
                    if((aflag+bflag+cflag+eflag)==4)
                    {
                        facet_abce_found = true;   
                    }
                }
            }
            // if not all 5 are found after scouring EVERY PENT, that means we have a boundary pent. this is the equivelant of saying if we only have 1 bordering cell to an edge in 2D, it's on the boundary
            if(facet_abcd_found==false)
            {
                boundary_facets[boundary_facet_counter][0] = ia;  
                boundary_facets[boundary_facet_counter][1] = ib;  
                boundary_facets[boundary_facet_counter][2] = ic;  
                boundary_facets[boundary_facet_counter][3] = id;  
                boundary_facet_2_cell[boundary_facet_counter] = i;
                boundary_facet_counter++;
            }
            if(facet_bcde_found==false)
            {
                boundary_facets[boundary_facet_counter][0] = ib;  
                boundary_facets[boundary_facet_counter][1] = ic;  
                boundary_facets[boundary_facet_counter][2] = id;  
                boundary_facets[boundary_facet_counter][3] = ie;  
                boundary_facet_2_cell[boundary_facet_counter] = i;
                boundary_facet_counter++;
            }
            if(facet_acde_found==false)
            {
                boundary_facets[boundary_facet_counter][0] = ia;  
                boundary_facets[boundary_facet_counter][1] = ic;  
                boundary_facets[boundary_facet_counter][2] = id;  
                boundary_facets[boundary_facet_counter][3] = ie;  
                boundary_facet_2_cell[boundary_facet_counter] = i;
                boundary_facet_counter++;
            }
            if(facet_abde_found==false)
            {
                boundary_facets[boundary_facet_counter][0] = ia;  
                boundary_facets[boundary_facet_counter][1] = ib;  
                boundary_facets[boundary_facet_counter][2] = id;  
                boundary_facets[boundary_facet_counter][3] = ie;  
                boundary_facet_2_cell[boundary_facet_counter] = i;
                boundary_facet_counter++;
            }
            if(facet_abce_found==false)
            {
                boundary_facets[boundary_facet_counter][0] = ia;  
                boundary_facets[boundary_facet_counter][1] = ib;  
                boundary_facets[boundary_facet_counter][2] = ic;  
                boundary_facets[boundary_facet_counter][3] = ie;  
                boundary_facet_2_cell[boundary_facet_counter] = i;
                boundary_facet_counter++;
            }
        }

       
        std::cout << "Boundary facets " << boundary_facet_counter << std::endl;
        /*
        for(int i = 0; i< boundary_facet_counter; ++i)
        {
             for(int k=0; k<4; ++k)
             {
                 std::cout << " " << boundary_facets[i][k];
             }
             std::cout << " " << std::endl;
        }
        */
        // now that we have all the boundary facets, we can find the boundary normals to do the correction
        std::vector< std::array<double,4> > boundary_normal;
        boundary_normal.reserve(boundary_facet_counter); // boundary_normal contains a vector for each of the boundary facets

        for(int i = 0; i < boundary_facet_counter; ++i) // loop through all of the boundary facets
        {   // the v's are the nodes on the boundary facet (equivelant of the nodes on the boundary edges)
            int v0, v1, v2, v3, v4;
            v0 = boundary_facets[i][0];
            v1 = boundary_facets[i][1];
            v2 = boundary_facets[i][2];
            v3 = boundary_facets[i][3];

            int index_cell = boundary_facet_2_cell[i]; // index cell is the cell that contains the i'th boundary facet
            // the i's are the nodes that bound this index cell
            int i0, i1, i2, i3, i4;
            i0 = pent[index_cell][0];
            i1 = pent[index_cell][1];
            i2 = pent[index_cell][2];
            i3 = pent[index_cell][3];
            i4 = pent[index_cell][4];
            // v4 is the one node in the cell that isn't on the boundary facet
            if((i0!=v0) && (i0!=v1) && (i0!=v2) && (i0!=v3))
            {
                v4 = i0;
            }
            if((i1!=v0) && (i1!=v1) && (i1!=v2) && (i1!=v3))
            {
                v4 = i1;
            }
            if((i2!=v0) && (i2!=v1) && (i2!=v2) && (i2!=v3))
            {
                v4 = i2;
            }
            if((i3!=v0) && (i3!=v1) && (i3!=v2) && (i3!=v3))
            {
                v4 = i3;
            }
            if((i4!=v0) && (i4!=v1) && (i4!=v2) && (i4!=v3))
            {
                v4 = i4;
            }
            // a corresonds with v0, b with v1, c with v2, and d with v3
            double ax, ay, az, aw;
            double bx, by, bz, bw;
            double cx, cy, cz, cw;
            double dx, dy, dz, dw;

            ax = x[v0];
            bx = x[v1];
            cx = x[v2];
            dx = x[v3];

            ay = y[v0];
            by = y[v1];
            cy = y[v2];
            dy = y[v3];

            az = z[v0];
            bz = z[v1];
            cz = z[v2];
            dz = z[v3];

            aw = w[v0];
            bw = w[v1];
            cw = w[v2];
            dw = w[v3];
            // bax is from a to b, and so on and so forth
            double bax, bay, baz, baw;
            double cax, cay, caz, caw;
            double dax, day, daz, daw;

            bax = bx - ax;
            bay = by - ay;
            baz = bz - az;
            baw = bw - aw;

            cax = cx - ax;
            cay = cy - ay;
            caz = cz - az;
            caw = cw - aw;

            dax = dx - ax;
            day = dy - ay;
            daz = dz - az;
            daw = dw - aw;

            //Normal vector components
            double i_comp =  three_determinant(bay, baz, baw, cay, caz, caw, day, daz, daw)*(1./6.);
            double j_comp = -three_determinant(bax, baz, baw, cax, caz, caw, dax, daz, daw)*(1./6.);
            double k_comp = three_determinant(bax, bay, baw, cax, cay, caw, dax, day, daw)*(1./6.);
            double l_comp = -three_determinant(bax, bay, baz, cax, cay, caz, dax, day, daz)*(1./6.);

            double x_cent, y_cent, z_cent, w_cent;
            // the centroid of the boundary facet
            x_cent = (ax + bx + cx + dx)*(1./4.);
            y_cent = (ay + by + cy + dy)*(1./4.);
            z_cent = (az + bz + cz + dz)*(1./4.);
            w_cent = (aw + bw + cw + dw)*(1./4.);
            // we dot the normal vector with the vector connecting the opposing node v4 and the centroid
            double normal_centroid_dot = (x_cent - x[v4])*i_comp + (y_cent - y[v4])*j_comp
                                         + (z_cent - z[v4])*k_comp + (w_cent - w[v4])*l_comp;
            // a positive dot product means it's oriented correctly, and a negative dot product means you flip it
            if(normal_centroid_dot > 0)
            {
                boundary_normal[i][0] = i_comp;
                boundary_normal[i][1] = j_comp;
                boundary_normal[i][2] = k_comp;
                boundary_normal[i][3] = l_comp;
            }
            else
            {
                boundary_normal[i][0] = -i_comp;
                boundary_normal[i][1] = -j_comp;
                boundary_normal[i][2] = -k_comp;
                boundary_normal[i][3] = -l_comp;
            }

            //correct the directed hyperarea vectors
            for(int m=0; m<10; ++m)
            {
               int index_edge = cell2edge[index_cell][m];
            // n0 and n1 are the 2 nodes that bound the index edge (these are all of the edges that comprise the boundary cell)
               int n0 = global_edges[index_edge][0];
               int n1 = global_edges[index_edge][1];

               //Check if on boundary
               if(((n0==v0) || (n0==v1) || (n0==v2) || (n0==v3)) &&
                  ((n1==v0) || (n1==v1) || (n1==v2) || (n1==v3))) // if this edge is on the boundary, perform the correction
               {
                    for(int q = 0; q<4; q++)
                    {
                        n01[index_edge][q] += boundary_normal[i][q]*(1./20.);
                    }
               }
            }
        }

        //**********************
        //Check boundary
        //**********************

        std::vector<double> boundary_residual = {0.0, 0.0, 0.0, 0.0}; // boundary residual[j] is the sum of all the boundary normals

        double surface_hyperarea = 0.0;
        for(int i = 0; i<boundary_facet_counter; ++i)
        {
             double mag = 0.0; // mag is the sum of the magnitudes of all the boundary nromals

             //std::cout << "Boundary normal " << i << std::endl;
             for(int j=0; j<4; ++j)
             {
                 //std::cout << " " << boundary_normal[i][j];
                 boundary_residual[j] += boundary_normal[i][j];
                 mag += (boundary_normal[i][j]*boundary_normal[i][j]);
             }
             //std::cout << " " << std::endl;

             mag = sqrt(mag);
             surface_hyperarea += mag;
        }

        std::cout << "Boundary residual " << boundary_residual[0] << " " << boundary_residual[1] << " " << boundary_residual[2] << " " << boundary_residual[3] << std::endl;

        std::cout << "Surface hyperarea " << std::setprecision(16) << surface_hyperarea << std::endl;

        //**************************
        //Recheck directed hyperarea vectors
        //**************************
        std::vector< std::array<double, 4> > vector_residual_corrected; // this is the equivelant of aj that we had in the 2D version
        vector_residual_corrected.reserve(n_points);

        for(int i = 0; i<n_points; ++i)
        {
            vector_residual_corrected[i][0] = 0.;
            vector_residual_corrected[i][1] = 0.;
            vector_residual_corrected[i][2] = 0.;
            vector_residual_corrected[i][3] = 0.;
        }

        std::cout << "N points " << n_points << std::endl;
        for(int i = 0; i<n_edges; ++i)
        {
            int e0 = global_edges[i][0];
            int e1 = global_edges[i][1];

            for(int j=0; j<4; ++j)
            {
                vector_residual_corrected[e0][j] += n01[i][j]; // we add n01 to the lower node on the edge
                vector_residual_corrected[e1][j] -= n01[i][j]; // and subtract it from the upper node
            }
        }

        for(int i = 0; i < boundary_facet_counter; ++i)
        {
            int v0, v1, v2, v3, v4;
            v0 = boundary_facets[i][0];
            v1 = boundary_facets[i][1];
            v2 = boundary_facets[i][2];
            v3 = boundary_facets[i][3];

            for(int k=0; k<4; ++k)
            {
                vector_residual_corrected[v0][k] += boundary_normal[i][k]*(1./4.); // we add the boundary normal to all 4 of the nodes on the boundary facet
                vector_residual_corrected[v1][k] += boundary_normal[i][k]*(1./4.);
                vector_residual_corrected[v2][k] += boundary_normal[i][k]*(1./4.);
                vector_residual_corrected[v3][k] += boundary_normal[i][k]*(1./4.);
            }
        }

        //std::cout << "Vector residuals (corrected)" << std::endl;
        for(int i = 0; i<n_points; ++i)
        {
            for(int k=0; k<4; k++)
            {
                //std::cout << " " << vector_residual_corrected[i][k];
                if(abs(vector_residual_corrected[i][k])>0.00000000001)
                {
                     std::cout << "WARNING: UNSAFE VALUE" << std::endl;
                }
            }
            //std::cout << " " << std::endl;
        }

        //**************************
        //Compute dual hypervolumes
        //**************************

        std::vector<double> dual_hypervolumes(n_points); // dual hypervolumes has one slot per point

        for(int i = 0; i<n_points; ++i)
        {
            dual_hypervolumes[i] = 0.; // initialize them all to be 0
        }


        for(int i = 0; i<n_edges; ++i)
        {
            int e0 = global_edges[i][0];
            int e1 = global_edges[i][1];

            dual_hypervolumes[e0] += (1./8.)*( (x[e1] - x[e0])*n01[i][0]
                                              +(y[e1] - y[e0])*n01[i][1]
                                              +(z[e1] - z[e0])*n01[i][2]
                                              +(w[e1] - w[e0])*n01[i][3] ); // add up the dual hypervolumes on each point based on the normal vectors and connectivity to neighbor
            
            dual_hypervolumes[e1] += (1./8.)*( (x[e1] - x[e0])*n01[i][0]
                                              +(y[e1] - y[e0])*n01[i][1]
                                              +(z[e1] - z[e0])*n01[i][2]
                                              +(w[e1] - w[e0])*n01[i][3] );
        }
 

        double total_hypervolume = 0.0;

        for(int i = 0; i<n_points; ++i)
        {
            total_hypervolume += dual_hypervolumes[i]; // add them all up
        }
        
        std::cout << "Hypervolume " << std::setprecision(16) << total_hypervolume << std::endl;

        //*********************
        // Least-Squares Coefficients (for Gradients)
        //*********************
            
        std::vector< std::vector< std::array<double,4> > > coefficient_vector; // this is a vector containing vectors which contain arrays
        coefficient_vector.reserve(n_points); // coefficient vector has n_points elements. each of these elements is a vector containing number_of_neighbors_1D elements. and each of those elements is a 1x4 array

        for(int i = 0; i<n_points; ++i) // loop through the outermost contents of coefficient_vector
        {
            int number_of_neighbors_local = number_of_neighbors_1D[i]; // just the number of neighbors that the i'th point has
            int offset_local = offsets_1D[i]; // the offset which helps us find the specific neighbor node numbers in the neighbors_1D data structure

            coefficient_vector[i].reserve(number_of_neighbors_local);
            std::vector< std::array<double,4> > diff_coordinate(number_of_neighbors_local); // new vector containing local # of neighbors elements. each element is a 4 element vector

            for(int j=0; j<number_of_neighbors_local; ++j) // looping through the number of local neighbors (loop through the second-tier of the 3-tier vector coefficient_vector)
            {
               int index_neighbor = neighbors_1D[offset_local+(j+1)]; // the node number of the neighbor

               diff_coordinate[j][0] = x[index_neighbor] - x[i]; // the j'th element in diff coordinate contains the vector pointing from the i'th node to its neighbor j
               diff_coordinate[j][1] = y[index_neighbor] - y[i];
               diff_coordinate[j][2] = z[index_neighbor] - z[i];
               diff_coordinate[j][3] = w[index_neighbor] - w[i];
            }           

            std::array< std::array<double,4>, 4 > A_matrix; // A_matrix is a 4x4 array

            for(int j=0; j<4; ++j)
            {
               for(int k=0; k<4; ++k)
               {
                  A_matrix[j][k] = 0.; // initialize its elements to be 0 everywhere
               }
            }

            for(int j=0; j<number_of_neighbors_local; ++j) // loop through the neighbors of the i'th node
            {
                for(int k=0; k<4; ++k)
                {                                                                                       // [delta_x^2        delta_x*delta_y   delta_x*delta_z   delta_x*delta_w]
                    for(int m=0; m<4; ++m)                                                              // [delta_y*dela_x   delta_y^2         delta_y*delta_z   delta_y*delta_w]
                    {                                                                                   // [delta_z*delta_x  delta_z*delta_y   delta_z^2         delta_z*delta_w]
                        A_matrix[m][k] += diff_coordinate[j][m]*diff_coordinate[j][k];                  // [delta_w*delta_x  delta_w*delta_y   delta_w*delta_z   delta_w^2      ]
                    }
                }
            }

            double A_matrix_1D[16]; // 1D-ifying the A_matrix
            int mat_counter = 0;

            for(int k=0; k<4; ++k)
            {
                for(int m=0; m<4; ++m)
                {
                    A_matrix_1D[mat_counter++] = A_matrix[k][m];
                }
            }
            
            double A_matrix_inverse_1D[16];

            //Invert 4x4
            bool inverse_flag = four_by_four_inverse(A_matrix_1D, A_matrix_inverse_1D);
            
            std::array< std::array<double,4>, 4 > A_matrix_inverse;

            mat_counter = 0;
            for(int k=0; k<4; ++k)
            {
                for(int m=0; m<4; ++m)
                {
                    A_matrix_inverse[k][m] = A_matrix_inverse_1D[mat_counter++]; // this becomes the inverse of A_matrix
                }
            }

            for(int j=0; j<number_of_neighbors_local; ++j)
            {
                for(int k=0; k<4; ++k)
                {
                   coefficient_vector[i][j][k] = 0.; // initialize coefficient vector to be 0 everywhere
                }
            }
 
            for(int j=0; j<number_of_neighbors_local; ++j) // update coefficient vector with A^(-1) * diff
            {
                for(int k=0; k<4; ++k)
                {
                    for(int m=0; m<4; ++m)
                    {
                        coefficient_vector[i][j][k] += A_matrix_inverse[k][m]*diff_coordinate[j][m];
                    }
                }
            }
        }     

        //*********************
        // Initialization
        //*********************

        std::vector<double> U(n_points); // vector called U that has n_points elements
        std::vector<double> Uinit(n_points); // vector called Uinit taht has n_points elements. this is the initial values
        std::vector<double> S(n_points); // vector called S which contains the forcing
        std::vector<double> residual(n_points); // vector called residual

        for(int i = 0; i<n_points; ++i)
        {
            U[i] = exact_solution(x[i],y[i],z[i],w[i]); // initializing values for U and S based on the definitions of exact_solution and forcing_term
            S[i] = forcing_term(x[i],y[i],z[i],w[i]);
        }
        
        
        int number_of_iterations = 1000;
        std::vector<double> local_time_steps(n_points); // local_time_steps is a vector that has an element for each point
        double CFL = 0.5;

        double residual_threshold = 0.000000000001;

        for(int m=0; m<number_of_iterations; ++m) // m corresponds to the current iteration that we're on
        {
             for(int i=0; i<n_points; ++i)
             {
                local_time_steps[i] = 0.; // set local_time_steps to 0 at each point
             }

             for(int i=0; i<n_edges; ++i) // loop through all the edges
             {
                 int n0 = global_edges[i][0]; // n0 is the lower node, n1 is the upper node
                 int n1 = global_edges[i][1];

                 std::vector<double> n01_unit(4); //n01_unit has 4 elements and is the unit vector of n01
                 
                 double mag = 0.;

                 for(int j=0; j<4; ++j)
                 {
                     mag += n01[i][j]*n01[i][j]; // mag is the magnitude of n01 for this edge
                 }
                 mag = sqrt(mag);

                 for(int j=0; j<4; ++j)
                 {
                     n01_unit[j] = n01[i][j]/mag; // n01 unit is the unit vector in the direction of n01
                 }
                 
                 //XXX Needs to be modified if PDE changes
                 std::vector<double> advection_direction = {1., 1., 1., 1.}; // advection direction is initialized to be 1 in all diretictions

                 double projected_direction = 0.;
                 for(int j=0; j<4; ++j)
                 {
                    projected_direction += advection_direction[j]*n01_unit[j]; // "projected direction" is the advection direction dotted with the unit normal direction
                 }
                 
                 local_time_steps[n0] += 0.5*abs(projected_direction)*mag; // the magnitude of projected direction is multiplied by the magnitude of the normal vector and added to node 0 and node 1
                 local_time_steps[n1] += 0.5*abs(projected_direction)*mag;
             }

             for(int i=0; i<n_points; ++i)
             {
                 local_time_steps[i] = CFL*dual_hypervolumes[i]/local_time_steps[i]; // local time step at the i'th node is updated to be CFL (const) times the dual hypervolume divided by the last local time step
                 //std::cout << local_time_steps[i] << std::endl;
             }

             //Stage 1
             for(int i=0; i<n_points; ++i)
             {
                 Uinit[i] = U[i]; // U initial is initialized as just U
             }

             residual_computation(n_points, n_edges, x, y, z, w, number_of_neighbors_1D, offsets_1D, neighbors_1D,  
                             coefficient_vector, global_edges, n01, dual_hypervolumes, U, S, residual); // perform the residual computation (see above)

             for(int i=0; i<n_points; ++i) // loop again through the points
             {
                 double xval = x[i]; // call x,y,z,wval the x,y,z,w coordinate at the current point
                 double yval = y[i];
                 double zval = z[i];
                 double wval = w[i];

                 double ethresh = 0.00000000001;
                // establish max/min values for x,y,z,w
                 double xmax = 1.0;
                 double xmin = 0.0;

                 double ymax = 1.0;
                 double ymin = 0.0;

                 double zmax = 1.0;
                 double zmin = 0.0;

                 double wmax = 1.0;
                 double wmin = 0.0;

                 if( (abs(xval-xmax) < ethresh) || (abs(xval-xmin) <ethresh) 
                    || (abs(yval-ymax) < ethresh) || (abs(yval-ymin) <ethresh) 
                    || (abs(zval-zmax) < ethresh) || (abs(zval-zmin) <ethresh) 
                    || (abs(wval-wmax) < ethresh) || (abs(wval-wmin) <ethresh) ) // if we're not so close to the max/min that we are less than this incredibly small threshold, we're chill
                 {
                    //Do nothing
                 }
                 else // if for some reason we are super close to the max/min
                 {
                     U[i] -= (local_time_steps[i]/dual_hypervolumes[i])*residual[i]; // perform a correction of sorts
                 }
             } 
 
             //Stage 2
             residual_computation(n_points, n_edges, x, y, z, w, number_of_neighbors_1D, offsets_1D, neighbors_1D,  
                             coefficient_vector, global_edges, n01, dual_hypervolumes, U, S, residual); // do the same computation again
            
             for(int i=0; i<n_points; ++i)
             {
                 double xval = x[i];
                 double yval = y[i];
                 double zval = z[i];
                 double wval = w[i];

                 double ethresh = 0.00000000001;

                 double xmax = 1.0;
                 double xmin = 0.0;

                 double ymax = 1.0;
                 double ymin = 0.0;

                 double zmax = 1.0;
                 double zmin = 0.0;

                 double wmax = 1.0;
                 double wmin = 0.0;

                 if( (abs(xval-xmax) < ethresh) || (abs(xval-xmin) <ethresh) 
                    || (abs(yval-ymax) < ethresh) || (abs(yval-ymin) <ethresh) 
                    || (abs(zval-zmax) < ethresh) || (abs(zval-zmin) <ethresh) 
                    || (abs(wval-wmax) < ethresh) || (abs(wval-wmin) <ethresh) ) // same as last time
                 {
                    //Do nothing
                 }
                 else
                 {
                     U[i] = 0.5*(U[i] + Uinit[i]) - 0.5*(local_time_steps[i]/dual_hypervolumes[i])*residual[i]; // perform different correction
                 }
             } 

             double max_residual = -1.;
             for(int i=0; i<n_points; ++i)
             {
                 if(abs(residual[i]) > max_residual)
                 {
                     double xval = x[i];
                     double yval = y[i];
                     double zval = z[i];
                     double wval = w[i];

                     double ethresh = 0.00000000001;

                     double xmax = 1.0;
                     double xmin = 0.0;

                     double ymax = 1.0;
                     double ymin = 0.0;

                     double zmax = 1.0;
                     double zmin = 0.0;

                     double wmax = 1.0;
                     double wmin = 0.0;

                     if( (abs(xval-xmax) < ethresh) || (abs(xval-xmin) <ethresh) 
                        || (abs(yval-ymax) < ethresh) || (abs(yval-ymin) <ethresh) 
                        || (abs(zval-zmax) < ethresh) || (abs(zval-zmin) <ethresh) 
                        || (abs(wval-wmax) < ethresh) || (abs(wval-wmin) <ethresh) ) // do it one more time
                     {
                        //Do nothing
                     }
                     else
                     {
                         max_residual = abs(residual[i]); // perform another correction
                     }
                 }
             }            

             if(max_residual < residual_threshold)
             {
                 std::cout << "Final iteration " << m << std::endl;
                 break;
             }
 
             std::cout << "Iteration number " << m << std::endl;
             std::cout << "Maximum residual " << max_residual << std::endl;
        }

        //*************************
        // Check interior residuals
        //*************************
        /*
        for(int i=0; i<n_points; ++i)
        {
            double xval = x[i];
            double yval = y[i];
            double zval = z[i];
            double wval = w[i];

            double ethresh = 0.00000000001;

            double xmax = 1.0;
            double xmin = 0.0;

            double ymax = 1.0;
            double ymin = 0.0;

            double zmax = 1.0;
            double zmin = 0.0;

            double wmax = 1.0;
            double wmin = 0.0;

            if( (abs(xval-xmax) < ethresh) || (abs(xval-xmin) <ethresh) 
               || (abs(yval-ymax) < ethresh) || (abs(yval-ymin) <ethresh) 
               || (abs(zval-zmax) < ethresh) || (abs(zval-zmin) <ethresh) 
               || (abs(wval-wmax) < ethresh) || (abs(wval-wmin) <ethresh) )
            {
               //Do nothing
            }
            else
            {
               std::cout << "Residual " << i << ": " << residual[i] << std::endl;
            }
        }
        */
}


