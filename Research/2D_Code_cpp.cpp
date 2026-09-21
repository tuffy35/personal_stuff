#include <iostream> // write to the standard input/output streams
#include <iomanip> // manipulate the format of input and output of program
#include <fstream> // ability to read from and write to files
#include <cmath>  // basic math stuff
#include <array> // allows more easy-to-work-with arrays and array methods
#include <vector> // allows for easier work with vectors
#include <algorithm> // allows for the "find" command

int main(void) 
{
    // read in coordinates
    std::ifstream coord_file;
    coord_file.open("mesh0_coords.txt",std::ios_base::in);

    double a;
    int total_entries = 0;
    while(coord_file >> a)
    {
        total_entries ++ ;
    }
    coord_file.close();
    
    int n_nodes = total_entries/2;

    std::vector<double> x;
    std::vector<double> y;

    x.reserve(n_nodes);
    y.reserve(n_nodes);

    coord_file.open("mesh0_coords.txt",std::ios_base::in);
    for(int i; i<n_nodes; i++)
    {
        coord_file >> x[i];
        coord_file >> y[i];
    }
    coord_file.close();
    
    
    // read in connectivity
    std::ifstream cell_file;
    cell_file.open("mesh0_connectivity.txt",std::ios_base::in); // reading in the connectivity file

    int b;
    int total_indices = 0;
    while(cell_file >> b)
    {
        total_indices++;
    }
    cell_file.close();

    int n_tet = total_indices/3; // divide by 3 since there's 3 coordinates per triangles

    std::vector< std::array<int,3> > meshing; // creating a template for the meshing array
    meshing.reserve(n_tet); // there are as many little arrays in meshing as there are triangles in the mesh

    cell_file.open("mesh0_connectivity.txt",std::ios_base::in);
    for(int i=0; i<n_tet; i++)
    {
        cell_file >> meshing[i][0];
        cell_file >> meshing[i][1];
        cell_file >> meshing[i][2];
    }
    cell_file.close();
    
    // NOTE: this assumes the indexing starts at 0 in the mesh file. If it starts at 1, you need to subtract 1 from each of the indices


    // Neighbor Data Structure

    std::vector<std::vector<int>> neighbor_nodes;

    for(int i=0; i<n_nodes; i++) // loop through each of the nodes in the mesh
    {
        std::vector<int> local_neighbors;
        for(int j=0; j<n_tet; j++)
        {
            if( std::find(meshing[j].begin(), meshing[j].end(), i) != meshing[j].end() ) // if the current node i is in the j'th triangle
            {
                for(int k=0; k<3; k++)
                {
                    if( meshing[j][k] != i && std::find(local_neighbors.begin(), local_neighbors.end(), meshing[j][k]) == local_neighbors.end())
                    {
                        local_neighbors.push_back(meshing[j][k]);
                    }
                }
            }
        }

        neighbor_nodes.push_back(local_neighbors);
    }

    std::cout << "There are " << n_nodes << "nodes";
    std:: cout << "\n";
    for (int i = 0; i < n_nodes; i++)
    {
        std::cout << "Node " << i << ": ";
        for (int val : neighbor_nodes[i])
        {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }
    
}