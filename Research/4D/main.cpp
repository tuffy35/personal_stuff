#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>

struct neighbor_data {
   int number_of_neighbors;
   std::vector<int> neighbors;
};

double three_determinant(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3)
{
    double fval = 0;

    fval = a1*(b2*c3-c2*b3) -a2*(b1*c3-c1*b3) + a3*(b1*c2-c1*b2);

    return fval;
}


double exact_solution(double x, double y, double z, double w)
{
    //double fval = x + 2.*y + 3.*z + 4.*w;
    double fval = exp(x + y + z + w);

    return fval;
}

double forcing_term(double x, double y, double z, double w)
{
    double fval = 4.0*exp(x + y + z + w);

    return fval;
}


bool four_by_four_inverse(const double m[16], double invOut[16])
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

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
    {
        return false;
    }

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
    {
        invOut[i] = inv[i] * det;
    }

    return true;
}


void residual_computation(int n_points, int n_edges, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& w, 
                          std::vector<int>& number_of_neighbors_1D, std::vector<int>& offsets_1D, std::vector<int>& neighbors_1D,  
                          std::vector< std::vector< std::array<double,4> > >& coefficient_vector, std::vector< std::array<int,2> >& global_edges, 
                          std::vector< std::array<double,4> >& n01, std::vector<double>& dual_hypervolumes,
                          std::vector<double>& U, std::vector<double>& S, std::vector<double>& residual)
{
    for(int i = 0; i<n_points; ++i)
    {
        residual[i] = 0.;
    }
    
    //*********************
    // Compute Gradients
    //*********************

    std::vector< std::array<double,4> > gradient;
    gradient.reserve(n_points);

    for(int i = 0; i<n_points; ++i)
    {
        for(int j=0; j<4; ++j)
        {
            gradient[i][j] = 0.;
        }

        int number_of_neighbors_local = number_of_neighbors_1D[i];
        int offset_local = offsets_1D[i];

        for(int j=0; j<number_of_neighbors_local; ++j)
        {
            int index_neighbor = neighbors_1D[offset_local+(j+1)];
            for(int k=0; k<4; ++k)
            {
                gradient[i][k] += coefficient_vector[i][j][k]*(U[index_neighbor] - U[i]);
            }
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

    for(int i = 0; i<n_edges; ++i)
    {
        int n0, n1;
        n0 = global_edges[i][0];
        n1 = global_edges[i][1];

        double U0, U1;
        U0 = U[n0];
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
            n01_unit[j] = n01[i][j]/mag;
        }
        
        //XXX Needs to be modified if PDE changes
        std::vector<double> advection_direction = {1., 1., 1., 1.};

        double projected_direction = 0.;
        for(int j=0; j<4; ++j)
        {
           projected_direction += advection_direction[j]*n01_unit[j];
        }

        double flux_temp;
        //flux_temp = 0.5*projected_direction*(U0+U1)*mag;
        double UL, UR;
        
        UL = U0 + gradient[n0][0]*(x[n1]-x[n0])*0.5
                + gradient[n0][1]*(y[n1]-y[n0])*0.5
                + gradient[n0][2]*(z[n1]-z[n0])*0.5
                + gradient[n0][3]*(w[n1]-w[n0])*0.5;

        UR = U1 - gradient[n1][0]*(x[n1]-x[n0])*0.5
                - gradient[n1][1]*(y[n1]-y[n0])*0.5
                - gradient[n1][2]*(z[n1]-z[n0])*0.5
                - gradient[n1][3]*(w[n1]-w[n0])*0.5;

        flux_temp = (0.5*projected_direction*(UL + UR) - 0.5*abs(projected_direction)*(UR - UL))*mag;

        residual[n0] += flux_temp;
        residual[n1] -= flux_temp;
    }

    //************************
    // Node Loop (forcing term)
    //************************

    for(int i=0; i<n_points; ++i)
    {
       residual[i] -= S[i]*dual_hypervolumes[i];
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
        while(coord_file >> a)
        {
            total_entries++;
        }
        coord_file.close();
        
        int n_points = total_entries/4;        
        std::cout << "N points " << n_points << std::endl;

        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        std::vector<double> w;

        x.reserve(n_points);
        y.reserve(n_points);
        z.reserve(n_points);
        w.reserve(n_points);

        coord_file.open("coordinates_mesh_4.txt", std::ios_base::in);

        for(int i=0; i<n_points; i++)
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
        
        int n_pentatopes = total_indexes/5;        
        std::cout << "N cells " << n_pentatopes << std::endl;

        std::vector< std::array<int,5> > pent;
        pent.reserve(n_pentatopes);

        cell_file.open("cells_mesh_4.txt", std::ios_base::in);

        for(int i=0; i<n_pentatopes; i++)
        {
            cell_file >> pent[i][0];
            cell_file >> pent[i][1];
            cell_file >> pent[i][2];
            cell_file >> pent[i][3];
            cell_file >> pent[i][4];
        }
        cell_file.close();
 
        //Assumes base 1 indexing
        for(int i=0; i<n_pentatopes; i++)
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

        int max_neighbors = 100;

        std::vector< std::array<int,4> > local_neighbor = { {1, 2, 3, 4}, {2, 0, 3, 4}, {0, 1, 3, 4}, {0, 1, 4, 2}, {0, 1, 2, 3} };
        std::vector< std::array<int,2> > local_edge = { {0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4} };

        std::vector<neighbor_data> total_neighbor_data;
        total_neighbor_data.reserve(n_points);

        for(int i = 0; i<n_points; ++i)
        {
            total_neighbor_data[i].number_of_neighbors = 0;
            total_neighbor_data[i].neighbors.reserve(max_neighbors);
        }

        //Loop over pentatopes
        for(int i = 0; i<n_pentatopes; ++i)
        {
           //Loop over vertices
           for(int k = 0; k<5; ++k)
           {
               int n1 = pent[i][k];

               //Loop over other nodes
               for(int j = 0; j<4; ++j)
               {
                  int n2 = pent[i][local_neighbor[k][j]];

                  //If no neighbors found yet
                  if( total_neighbor_data[n1].number_of_neighbors == 0)
                  {
                      total_neighbor_data[n1].number_of_neighbors = 1;
                      total_neighbor_data[n1].neighbors[0] = n2;
                  }
                  else
                  {
                      //Has neighbor been added?
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
                      if(!neighbor_found)
                      {
                          total_neighbor_data[n1].number_of_neighbors++;
                          total_neighbor_data[n1].neighbors[total_neighbor_data[n1].number_of_neighbors-1] = n2;
                      }
                  }
               }
           }
        }


        //*********************************
        //Store neighbor data in 1D arrays
        //*********************************

        int total_number_of_neighbors = 0;

        for(int i=0; i < n_points; ++i)
        {
           total_number_of_neighbors += total_neighbor_data[i].number_of_neighbors;
        }

        std::vector<int> neighbors_1D(total_number_of_neighbors);
        std::vector<int> offsets_1D(n_points);
        std::vector<int> number_of_neighbors_1D(n_points);

        int itemp = 0;
        offsets_1D[0] = -1;

        for(int i=0; i<n_points; ++i)
        {
            number_of_neighbors_1D[i] = total_neighbor_data[i].number_of_neighbors;
            
            for(int k = 0; k<number_of_neighbors_1D[i]; ++k)
            {
                itemp = offsets_1D[i] + (k+1);
                neighbors_1D[itemp] = total_neighbor_data[i].neighbors[k];
            }

            if(i< (n_points-1))
            {
                offsets_1D[i+1] = itemp;
            }
        }

        //*********************************
        //Edge numbering
        //*********************************
        int n_edges = 0;

        for(int i = 0; i<n_points; ++i)
        {
            for(int k = 0; k < number_of_neighbors_1D[i]; ++k)
            {
                 if(i < neighbors_1D[offsets_1D[i]+(k+1)])
                 {
                     n_edges++;
                 }
             }
        }
 
        //global edge
        std::vector< std::array<int,2> > global_edges;
        global_edges.reserve(n_edges);

        std::vector< std::array<int,100> > global_edge_numbers;
        global_edge_numbers.reserve(n_points);

        int edge_count = -1;
        for(int i = 0; i<n_points; ++i)
        {
            for(int k = 0; k < number_of_neighbors_1D[i]; ++k)
            {
                 if(i < neighbors_1D[offsets_1D[i]+(k+1)])
                 {
                     edge_count++;
                     global_edges[edge_count][0] = i;
                     global_edges[edge_count][1] = neighbors_1D[offsets_1D[i]+(k+1)];

                     global_edge_numbers[i][k] = edge_count;
                 }
            }
        }   


        //*********************************
        //Cell to Edge
        //*********************************
        std::vector< std::array<int,10> > cell2edge;
        cell2edge.reserve(n_pentatopes);


        for(int i = 0; i<n_pentatopes; ++i)
        {
            int index_edge;
            int n1, n2;

            //edge 0-1
            index_edge = 0;
            
            if( pent[i][0] < pent[i][1])
            {
                n1 = pent[i][0];
                n2 = pent[i][1];
            }
            else
            {
                n1 = pent[i][1];
                n2 = pent[i][0];
            }
 
            for(int k = 0; k < number_of_neighbors_1D[n1]; ++k)
            {
                if(neighbors_1D[offsets_1D[n1]+(k+1)] == n2)
                {
                    cell2edge[i][index_edge] = global_edge_numbers[n1][k];
                }
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

        std::vector< std::array<double,4> > n01;
        n01.reserve(n_edges);

        for(int i = 0; i<n_edges; ++i)
        {
            n01[i][0] = 0.; 
            n01[i][1] = 0.;
            n01[i][2] = 0.;
            n01[i][3] = 0.;
        }

        for(int i = 0; i<n_pentatopes; ++i)
        {
            for(int j=0; j<10; ++j)
            {
                int m;

                int index_edge = cell2edge[i][j];
                if( global_edges[index_edge][0] == pent[i][local_edge[j][0]] )
                {
                    m = local_edge[j][0];
                }
                else
                {
                    m = local_edge[j][1];
                }

                //Pentatope facet opposite to the node m
                int v0, v1, v2, v3;
                v0 = pent[i][local_neighbor[m][0]];
                v1 = pent[i][local_neighbor[m][1]];
                v2 = pent[i][local_neighbor[m][2]];
                v3 = pent[i][local_neighbor[m][3]];

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

                //std::cout << "Normal components " << i_comp << " " << j_comp << " " << k_comp << " " << l_comp << std::endl;

                //Edge vector
                int e0, e1;
                e0 = global_edges[index_edge][0];
                e1 = global_edges[index_edge][1];

                double edge_normal_dot = i_comp*(x[e1] - x[e0]) + j_comp*(y[e1] - y[e0]) 
                                       + k_comp*(z[e1] - z[e0]) + l_comp*(w[e1] - w[e0]);


                if(edge_normal_dot > 0)
                {
                    n01[index_edge][0] += i_comp;
                    n01[index_edge][1] += j_comp;
                    n01[index_edge][2] += k_comp;
                    n01[index_edge][3] += l_comp;
                }
                else
                {
                    n01[index_edge][0] -= i_comp;
                    n01[index_edge][1] -= j_comp;
                    n01[index_edge][2] -= k_comp;
                    n01[index_edge][3] -= l_comp;
                }
            }
        }

        for(int i = 0; i<n_edges; ++i)
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

        std::vector< std::array<double, 4> > vector_residual;
        vector_residual.reserve(n_points);

        for(int i = 0; i<n_points; ++i)
        {
            vector_residual[i][0] = 0.;
            vector_residual[i][1] = 0.;
            vector_residual[i][2] = 0.;
            vector_residual[i][3] = 0.;
        }

        std::cout << "N points " << n_points << std::endl;
        for(int i = 0; i<n_edges; ++i)
        {
            int e0 = global_edges[i][0];
            int e1 = global_edges[i][1];

            for(int j=0; j<4; ++j)
            {
                vector_residual[e0][j] += n01[i][j];
                vector_residual[e1][j] -= n01[i][j];
            }
        }


        //************************************
        //Identify Boundary facets
        //************************************

        std::vector< std::array<int,4> > boundary_facets;
        boundary_facets.reserve(2*n_pentatopes);

        std::vector<int> boundary_facet_2_cell;
        boundary_facet_2_cell.reserve(2*n_pentatopes);

        int boundary_facet_counter = 0;

        for(int i=0; i<n_pentatopes; ++i)
        {
            int ia, ib, ic, id, ie;

            ia = pent[i][0];
            ib = pent[i][1];
            ic = pent[i][2];
            id = pent[i][3];
            ie = pent[i][4];

            bool facet_abcd_found = false;
            bool facet_bcde_found = false;
            bool facet_acde_found = false;
            bool facet_abde_found = false;
            bool facet_abce_found = false;

            for(int j=0; j<n_pentatopes; ++j)
            {
                if(i!=j)
                {
                    int aflag, bflag, cflag, dflag, eflag;
                    
                    aflag = 0;
                    bflag = 0;
                    cflag = 0;
                    dflag = 0;
                    eflag = 0;

                    int ja, jb, jc, jd, je;
                    ja = pent[j][0];
                    jb = pent[j][1];
                    jc = pent[j][2];
                    jd = pent[j][3];
                    je = pent[j][4];
                    
                    if((ia==ja) || (ia==jb) || (ia==jc) || (ia==jd) || (ia==je))
                    {
                        aflag = 1;
                    }
                    if((ib==ja) || (ib==jb) || (ib==jc) || (ib==jd) || (ib==je))
                    {
                        bflag = 1;
                    }
                    if((ic==ja) || (ic==jb) || (ic==jc) || (ic==jd) || (ic==je))
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

        std::vector< std::array<double,4> > boundary_normal;
        boundary_normal.reserve(boundary_facet_counter);

        for(int i = 0; i < boundary_facet_counter; ++i)
        {
            int v0, v1, v2, v3, v4;
            v0 = boundary_facets[i][0];
            v1 = boundary_facets[i][1];
            v2 = boundary_facets[i][2];
            v3 = boundary_facets[i][3];

            int index_cell = boundary_facet_2_cell[i];

            int i0, i1, i2, i3, i4;
            i0 = pent[index_cell][0];
            i1 = pent[index_cell][1];
            i2 = pent[index_cell][2];
            i3 = pent[index_cell][3];
            i4 = pent[index_cell][4];

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
  
            x_cent = (ax + bx + cx + dx)*(1./4.);
            y_cent = (ay + by + cy + dy)*(1./4.);
            z_cent = (az + bz + cz + dz)*(1./4.);
            w_cent = (aw + bw + cw + dw)*(1./4.);

            double normal_centroid_dot = (x_cent - x[v4])*i_comp + (y_cent - y[v4])*j_comp
                                         + (z_cent - z[v4])*k_comp + (w_cent - w[v4])*l_comp;

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

               int n0 = global_edges[index_edge][0];
               int n1 = global_edges[index_edge][1];

               //Check if on boundary
               if(((n0==v0) || (n0==v1) || (n0==v2) || (n0==v3)) &&
                  ((n1==v0) || (n1==v1) || (n1==v2) || (n1==v3)))
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

        std::vector<double> boundary_residual = {0.0, 0.0, 0.0, 0.0};

        double surface_hyperarea = 0.0;
        for(int i = 0; i<boundary_facet_counter; ++i)
        {
             double mag = 0.0;

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
        std::vector< std::array<double, 4> > vector_residual_corrected;
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
                vector_residual_corrected[e0][j] += n01[i][j];
                vector_residual_corrected[e1][j] -= n01[i][j];
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
                vector_residual_corrected[v0][k] += boundary_normal[i][k]*(1./4.);
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

        std::vector<double> dual_hypervolumes(n_points);

        for(int i = 0; i<n_points; ++i)
        {
            dual_hypervolumes[i] = 0.;
        }


        for(int i = 0; i<n_edges; ++i)
        {
            int e0 = global_edges[i][0];
            int e1 = global_edges[i][1];

            dual_hypervolumes[e0] += (1./8.)*( (x[e1] - x[e0])*n01[i][0]
                                              +(y[e1] - y[e0])*n01[i][1]
                                              +(z[e1] - z[e0])*n01[i][2]
                                              +(w[e1] - w[e0])*n01[i][3] );
            
            dual_hypervolumes[e1] += (1./8.)*( (x[e1] - x[e0])*n01[i][0]
                                              +(y[e1] - y[e0])*n01[i][1]
                                              +(z[e1] - z[e0])*n01[i][2]
                                              +(w[e1] - w[e0])*n01[i][3] );
        }
 

        double total_hypervolume = 0.0;

        for(int i = 0; i<n_points; ++i)
        {
            total_hypervolume += dual_hypervolumes[i];
        }
        
        std::cout << "Hypervolume " << std::setprecision(16) << total_hypervolume << std::endl;

        //*********************
        // Least-Squares Coefficients (for Gradients)
        //*********************
            
        std::vector< std::vector< std::array<double,4> > > coefficient_vector;
        coefficient_vector.reserve(n_points);

        for(int i = 0; i<n_points; ++i)
        {
            int number_of_neighbors_local = number_of_neighbors_1D[i];
            int offset_local = offsets_1D[i];

            coefficient_vector[i].reserve(number_of_neighbors_local);
            std::vector< std::array<double,4> > diff_coordinate(number_of_neighbors_local);

            for(int j=0; j<number_of_neighbors_local; ++j)
            {
               int index_neighbor = neighbors_1D[offset_local+(j+1)];

               diff_coordinate[j][0] = x[index_neighbor] - x[i];
               diff_coordinate[j][1] = y[index_neighbor] - y[i];
               diff_coordinate[j][2] = z[index_neighbor] - z[i];
               diff_coordinate[j][3] = w[index_neighbor] - w[i];
            }           

            std::array< std::array<double,4>, 4 > A_matrix;

            for(int j=0; j<4; ++j)
            {
               for(int k=0; k<4; ++k)
               {
                  A_matrix[j][k] = 0.;
               }
            }

            for(int j=0; j<number_of_neighbors_local; ++j)
            {
                for(int k=0; k<4; ++k)
                {
                    for(int m=0; m<4; ++m)
                    {
                        A_matrix[m][k] += diff_coordinate[j][m]*diff_coordinate[j][k]; 
                    }
                }
            }

            double A_matrix_1D[16];
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
                    A_matrix_inverse[k][m] = A_matrix_inverse_1D[mat_counter++];
                }
            }

            for(int j=0; j<number_of_neighbors_local; ++j)
            {
                for(int k=0; k<4; ++k)
                {
                   coefficient_vector[i][j][k] = 0.;
                }
            }
 
            for(int j=0; j<number_of_neighbors_local; ++j)
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

        std::vector<double> U(n_points);
        std::vector<double> Uinit(n_points);
        std::vector<double> S(n_points);
        std::vector<double> residual(n_points);

        for(int i = 0; i<n_points; ++i)
        {
            U[i] = exact_solution(x[i],y[i],z[i],w[i]);
            S[i] = forcing_term(x[i],y[i],z[i],w[i]);
        }
        
        
        int number_of_iterations = 1000;
        std::vector<double> local_time_steps(n_points);
        double CFL = 0.5;

        double residual_threshold = 0.000000000001;

        for(int m=0; m<number_of_iterations; ++m)
        {
             for(int i=0; i<n_points; ++i)
             {
                local_time_steps[i] = 0.;
             }

             for(int i=0; i<n_edges; ++i)
             {
                 int n0 = global_edges[i][0];
                 int n1 = global_edges[i][1];

                 std::vector<double> n01_unit(4);
                 
                 double mag = 0.;

                 for(int j=0; j<4; ++j)
                 {
                     mag += n01[i][j]*n01[i][j];
                 }
                 mag = sqrt(mag);

                 for(int j=0; j<4; ++j)
                 {
                     n01_unit[j] = n01[i][j]/mag;
                 }
                 
                 //XXX Needs to be modified if PDE changes
                 std::vector<double> advection_direction = {1., 1., 1., 1.};

                 double projected_direction = 0.;
                 for(int j=0; j<4; ++j)
                 {
                    projected_direction += advection_direction[j]*n01_unit[j];
                 }
                 
                 local_time_steps[n0] += 0.5*abs(projected_direction)*mag;
                 local_time_steps[n1] += 0.5*abs(projected_direction)*mag;
             }

             for(int i=0; i<n_points; ++i)
             {
                 local_time_steps[i] = CFL*dual_hypervolumes[i]/local_time_steps[i];
                 //std::cout << local_time_steps[i] << std::endl;
             }

             //Stage 1
             for(int i=0; i<n_points; ++i)
             {
                 Uinit[i] = U[i];
             }

             residual_computation(n_points, n_edges, x, y, z, w, number_of_neighbors_1D, offsets_1D, neighbors_1D,  
                             coefficient_vector, global_edges, n01, dual_hypervolumes, U, S, residual);

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
                     U[i] -= (local_time_steps[i]/dual_hypervolumes[i])*residual[i];
                 }
             } 
 
             //Stage 2
             residual_computation(n_points, n_edges, x, y, z, w, number_of_neighbors_1D, offsets_1D, neighbors_1D,  
                             coefficient_vector, global_edges, n01, dual_hypervolumes, U, S, residual);

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
                     U[i] = 0.5*(U[i] + Uinit[i]) - 0.5*(local_time_steps[i]/dual_hypervolumes[i])*residual[i];
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
                        || (abs(wval-wmax) < ethresh) || (abs(wval-wmin) <ethresh) )
                     {
                        //Do nothing
                     }
                     else
                     {
                         max_residual = abs(residual[i]);
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


