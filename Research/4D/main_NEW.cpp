#include <iostream> // write to the standard input/output streams
#include <iomanip> // manipulate the format of input and output of program
#include <fstream> // ability to read from and write to files
#include <cmath>  // basic math stuff
#include <array> // allows more easy-to-work-with arrays and array methods
#include <vector> // allows for easier work with vectors

struct neighbor_data {
    int number_of_neighbors;
    std::vector <int> neighbor_nodes;
};

double three_determinant(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3){

    double func_value = 0;
    
    func_value = a1*(b2*c3 - b3*c2) - a2*(b1*c3 - b3*c1) + a3*(b1*c2 - b2*c1);

    return func_value;
}

int main(void) {

    // read in coordinates file
    std::ifstream coord_file;
    coord_file.open("coordinates_mesh_4.txt",std::ios_base::in);

    double a;
    int total_node_entries = 0;

    while(coord_file >> a)
    {
        total_node_entries ++;
    }
    coord_file.close();

    int n_nodes = total_node_entries/4; // coordinates file is n_nodes x 4, so the total amount of entries is 4x the n_nodes

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> w;

    x.reserve(n_nodes);
    y.reserve(n_nodes);
    z.reserve(n_nodes);
    w.reserve(n_nodes);

    coord_file.open("coordinates_mesh_4.txt",std::ios_base::in);

    for(int i=0; i<n_nodes; i++) // x,y,z,w are populated
    {
        coord_file >> x[i];
        coord_file >> y[i];
        coord_file >> z[i];
        coord_file >> w[i];
    }

    // read in connectivity

    std::ifstream cell_file;
    cell_file.open("cells_mesh_4.txt",std::ios_base::in); // reading in the connectivity file

    int b;
    int total_indices = 0;
    while(cell_file >> b)
    {
        total_indices++;
    }
    cell_file.close();

    int n_pents = total_indices/5; // divide by 5 since there are 5 coordinates per pentatope

    // establish pentatope data structure
    std::vector< std::array<int,5> > pent;
    pent.reserve(n_pents);

    cell_file.open("cells_mesh_4.txt",std::ios_base::in); // reading in the connectivity file
    for(int i =0; i<n_pents; i++)
    {
        cell_file >> pent[i][0];
        cell_file >> pent[i][1];
        cell_file >> pent[i][2];
        cell_file >> pent[i][3];
        cell_file >> pent[i][4];
    }
    cell_file.close();

    // if the mesh file is using base 1 indexing, deduct 1 from each index

    for(int i=0; i < n_pents; i++)
    {
        pent[i][0] --;
        pent[i][1] --;
        pent[i][2] --;
        pent[i][3] --;
        pent[i][4] --;
    }

    // define neighbor nodes structure

    int max_neighbors = 100;

    // local_neighbor[i] is the list of which local node numbers are the neighbors to the i'th node
    std::vector < std::array < int, 4 > > local_neighbor = {{1, 2, 3, 4}, {2, 0, 3, 4}, {0, 1, 3, 4}, {0, 1, 4, 2}, {0, 1, 2, 3} };
    
    // local_edge[i] are the node numbers of the i'th local edge of a certain pentatope
    std::vector < std::array < int, 2 > > local_edge = {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};


    std::vector < neighbor_data > total_neighbor_data;
    total_neighbor_data.reserve(n_nodes);
    
    // initialize total_neighbor_data to be zeros everywhere
    for(int i = 0; i < n_nodes; i++)
    {
        total_neighbor_data[i].number_of_neighbors = 0;
        total_neighbor_data[i].neighbor_nodes.reserve(max_neighbors);
    }

    for(int i = 0; i < n_pents; i++) // loop over all the pents
    {
        for(int k = 0; k < 5; k++) // loop over the nodes on a given pent
        {
            int n1 = pent[i][k]; // call current node "n1"
            
            for(int j = 0; j < 4; j++) // loop over the other nodes besides k on pent i
            {
                int n2 = pent[i][local_neighbor[k][j]]; // call current neighbor "n2"
                
                if(total_neighbor_data[n1].number_of_neighbors == 0) // what to do if it's the first node being added to this node's data
                {
                    total_neighbor_data[n1].number_of_neighbors = 1;
                    total_neighbor_data[n1].neighbor_nodes[0] = n2;
                }

                bool neighbor_found = false;

                for(int m = 0; m < total_neighbor_data[n1].number_of_neighbors; m++)
                {
                    if(total_neighbor_data[n1].neighbor_nodes[m] == n2)
                    {
                        neighbor_found = true;
                        break;
                    }
                }

                if(!neighbor_found)
                {
                    total_neighbor_data[n1].number_of_neighbors ++;
                    total_neighbor_data[n1].neighbor_nodes[total_neighbor_data[n1].number_of_neighbors-1] = n2;
                }
            }
        }
    }

    // write a 1-D data structure for neighbor_data
    int total_number_of_neighbors = 0;
    for(int i=0; i < n_nodes; i ++)
    {
        total_number_of_neighbors += total_neighbor_data[i].number_of_neighbors;
    }

    std::vector < int > number_of_neighbors_1D(n_nodes);
    std::vector < int > offsets_1D(n_nodes);
    std::vector < int > neighbors_1D(total_number_of_neighbors);

    // establish in a separate data structure how many neighbors a given node has

    int itemp = 0;
    offsets_1D[0] = -1;

    for(int i=0; i<n_nodes; i ++)
    {
        number_of_neighbors_1D[i] = total_neighbor_data[i].number_of_neighbors;

        for(int j = 0; j<number_of_neighbors_1D[i]; j++)
        {
            itemp = offsets_1D[i] + (j+1);
            neighbors_1D[itemp] = total_neighbor_data[i].neighbor_nodes[j];
        }
        
        if(i< (n_nodes-1)) // this ensures that we don't try to index offsets_1D at an index that doesn't exist
        {
            offsets_1D[i+1] = itemp;
        }
    }

    // global edge numbering system

    int n_edges = 0; // first count the number of edges

    for(int i =  0; i<n_nodes; i++)
    {
        for(int k = 0; k <number_of_neighbors_1D[i]; k++) // loop through each of the neighbors of each of the nodes
        {
            if(i < neighbors_1D[offsets_1D[i] + (k+1)])
            {
                n_edges ++; // only increment the edge counter if node i is lower number than its kth neighbor
            }
        }
    }

    // now we want to establish two data structures.

    // the first takes edge number as an input and returns the nodes as outputs
    std::vector < std::array < int, 2 > > global_edges;
    global_edges.reserve(n_edges);

    // the second takes node-neighbor pairs as inputs and returns the edge number as the output
    std::vector < std::array < int, 100 > > global_edge_numbers;
    global_edge_numbers.reserve(n_nodes);

    int edge_counter = -1;

    for(int i = 0; i < n_nodes; i++)
    {
        for(int k = 0; k <number_of_neighbors_1D[i]; k++) // again loop through all node-edge pairs
        {
            if(i < neighbors_1D[offsets_1D[i] + (k+1)])
            {
                edge_counter ++; // increment edge-counter
                global_edges[edge_counter][0] = i;
                global_edges[edge_counter][1] = neighbors_1D[offsets_1D[i] + (k+1)]; // when i < neighbor, add both nodes to current edge_counter

                global_edge_numbers[i][k] = edge_counter; // assign this combo in global_edge_numbers
            }
        }
    }

    // cell2edge
    // idea: loop over every pentatope and then manually through each of the 10 edges on the pent (use local_edges)
    // take that combo of node-neighbor pairs and extract the global edge number from global_edge_numbers
    // update that in the cell2edge entry for that pent

    std::vector < std::array <int,10> > cell2edge; // 10 edges per pentatope
    cell2edge.reserve(n_pents);

    for(int i = 0; i < n_pents; i++)
    {
        int index_edge;
        int n1;
        int n2;

        // edge 0-1
        index_edge = 0;
        if(pent[i][0] < pent[i][1])
        {
            n1 = pent[i][0];
            n2 = pent[i][1];
        }
        else
        {
            n2 = pent[i][0];
            n1 = pent[i][1];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }

        // edge 0-2
        index_edge = 1;
        if(pent[i][0] < pent[i][2])
        {
            n1 = pent[i][0];
            n2 = pent[i][2];
        }
        else
        {
            n2 = pent[i][0];
            n1 = pent[i][2];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }

        // edge 0-3
        index_edge = 2;
        if(pent[i][0] < pent[i][3])
        {
            n1 = pent[i][0];
            n2 = pent[i][3];
        }
        else
        {
            n2 = pent[i][0];
            n1 = pent[i][3];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }

        // edge 0-4
        index_edge = 3;
        if(pent[i][0] < pent[i][4])
        {
            n1 = pent[i][0];
            n2 = pent[i][4];
        }
        else
        {
            n2 = pent[i][0];
            n1 = pent[i][4];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }

        // edge 1-2
        index_edge = 4;
        if(pent[i][1] < pent[i][2])
        {
            n1 = pent[i][1];
            n2 = pent[i][2];
        }
        else
        {
            n2 = pent[i][1];
            n1 = pent[i][2];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }

        // edge 1-3
        index_edge = 5;
        if(pent[i][1] < pent[i][3])
        {
            n1 = pent[i][1];
            n2 = pent[i][3];
        }
        else
        {
            n2 = pent[i][1];
            n1 = pent[i][3];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }

        // edge 1-4
        index_edge = 6;
        if(pent[i][1] < pent[i][4])
        {
            n1 = pent[i][1];
            n2 = pent[i][4];
        }
        else
        {
            n2 = pent[i][1];
            n1 = pent[i][4];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }

        // edge 2-3
        index_edge = 7;
        if(pent[i][2] < pent[i][3])
        {
            n1 = pent[i][2];
            n2 = pent[i][3];
        }
        else
        {
            n2 = pent[i][2];
            n1 = pent[i][3];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }
        // edge 2-4
        index_edge = 8;
        if(pent[i][2] < pent[i][4])
        {
            n1 = pent[i][2];
            n2 = pent[i][4];
        }
        else
        {
            n2 = pent[i][2];
            n1 = pent[i][4];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }

        // edge 3-4
        index_edge = 9;
        if(pent[i][3] < pent[i][4])
        {
            n1 = pent[i][3];
            n2 = pent[i][4];
        }
        else
        {
            n2 = pent[i][3];
            n1 = pent[i][4];
        }

        for(int k = 0; k < number_of_neighbors_1D[n1]; k++)
        {
            if(neighbors_1D[offsets_1D[n1] + (k+1)] == n2)
            {
                cell2edge[i][index_edge] = global_edge_numbers[n1][k];
            }
        }
    }

    // directed hyperarea vectors

    // idea: 
    // 1) loop through pentatopes
    // 2) loop through local edges
    // 3) define how this pentatope contributes to the njk vector of that local edge
    // 3.1) list out the 4 nodes that AREN'T the node from which the edge is eminating
    // 3.2) define the vectors that connect those nodes
    // 3.3) compute the 4D cross product of all those vectors
    // 3.4) check it's orientation by dotting it with the vector j-to-k
    // 3.5) scale it and add it to njk contribution

    std::vector < std::array < double, 4 > > n01;
    n01.reserve(n_edges);

    for(int i = 0; i < n_edges; i++) // first, initialize all the vectors to be 0
    {
        n01[i][0] = 0.;
        n01[i][1] = 0.;
        n01[i][2] = 0.;
        n01[i][3] = 0.;
    }

    for(int i = 0; i < n_pents; i++) // loop through the pentatopes
    {
        for(int j = 0; j < 10; j++) // loop through the local edges on that pentatope
        {
            int m = 0; // m is the local node number (0-4) from with njk eminates from for this current iteration
            int index_edge = cell2edge[i][j];
            
            // which of the 5 nodes (0,1,2,3,4) on this pent is the node that we find from the global_edge data structure?
            if( global_edges[index_edge][0] == pent[i][local_edge[j][0]]) // check: is the node number the 1st item in local_edges (m is that first local node number in local_edges) or is it the other one?
            {
                m = local_edge[j][0];
            }
            else
            {
                m = local_edge[j][1];
            }

            // list out the node numbers that define the facet opposite of node m
            int a,b,c,d;
            a = pent[i][local_neighbor[m][0]];
            b = pent[i][local_neighbor[m][1]];
            c = pent[i][local_neighbor[m][2]];
            d = pent[i][local_neighbor[m][3]];

            // x,y,z,w coordinates of these four points
            double ax,ay,az,aw;
            double bx,by,bz,bw;
            double cx,cy,cz,cw;
            double dx,dy,dz,dw;

            ax = x[a];
            ay = y[a];
            az = z[a];
            aw = w[a];

            bx = x[b];
            by = y[b];
            bz = z[b];
            bw = w[b];

            cx = x[c];
            cy = y[c];
            cz = z[c];
            cw = w[c];

            dx = x[d];
            dy = y[d];
            dz = z[d];
            dw = w[d];

            // now we need vectors to define the connection of these points. for simplicity, let's say all vectors eminate from a.

            double bax, bay, baz, baw;
            double cax, cay, caz, caw;
            double dax, day, daz, daw;

            bax = bx-ax;
            bay = by-ay;
            baz = bz-az;
            baw = bw-aw;

            cax = cx-ax;
            cay = cy-ay;
            caz = cz-az;
            caw = cw-aw;

            dax = dx-ax;
            day = dy-ay;
            daz = dz-az;
            daw = dw-aw;

            // now we need to take the cross product of these 3 vectors

            double i_comp = three_determinant(bay,baz,baw,cay,caz,caw,day,daz,daw)*(1./6.);
            double j_comp = -three_determinant(bax,baz,baw,cax,caz,caw,dax,daz,daw)*(1./6.);
            double k_comp = three_determinant(bax,bay,baw,cax,cay,caw,dax,day,daw)*(1./6.);
            double l_comp = -three_determinant(bax,bay,baz,cax,cay,caz,dax,day,daz)*(1./6.);

            // check if this one dotted with the j-to-k vector is positive or negative

            int n1 = global_edges[index_edge][0];
            int n2 = global_edges[index_edge][1];

            double edge_vector_i = x[n2]-x[n1];
            double edge_vector_j = y[n2]-y[n1];
            double edge_vector_k = z[n2]-z[n1];
            double edge_vector_l = z[n2]-z[n1];

            double dot_prod = 0.;

            dot_prod = i_comp*edge_vector_i + j_comp*edge_vector_j + k_comp*edge_vector_k + l_comp*edge_vector_l;

            if(dot_prod < 0) // flip the vector orientation if it's not pointing from node j to node k
            {
                i_comp = -i_comp;
                j_comp = -j_comp;
                k_comp = -k_comp;
                l_comp = -l_comp;
            }

            // add all the contributions to the respective n01 entry
            n01[index_edge][0] += i_comp;
            n01[index_edge][1] += j_comp;
            n01[index_edge][2] += k_comp;
            n01[index_edge][3] += l_comp;

        }
    }
    
    for(int i = 0; i < n_edges; i++) // scaling factor
    {
        n01[i][0] = n01[i][0]*(1./10.);
        n01[i][1] = n01[i][1]*(1./10.);
        n01[i][2] = n01[i][2]*(1./10.);
        n01[i][3] = n01[i][3]*(1./10.);
    }

    // boundary considerations

    std::vector < std::array < int,4 > > boundary_facets; // data structure that takes the boundary facet number as an input and returns the 4 nodes
    boundary_facets.reserve(2*n_pents);

    std::vector <int> boundary_facet_2_cell; // data structure that takes the boundary facet number as an input and returns its parent cell number
    boundary_facet_2_cell.reserve(2*n_pents);

    int boundary_facet_counter = 0;
    
    for(int i = 0; i < n_pents; i++)
    {
        int ia, ib, ic, id, ie; // the nodes on pent i 

        ia = pent[i][0];
        ib = pent[i][1];
        ic = pent[i][2];
        id = pent[i][3];
        ie = pent[i][4];

        for(int j = 0; j < n_pents; j++)
        {
            if(i!=j)
            {
                int ja, jb, jc, jd, je; // the nodes on the other pent j
                
                ja = pent[j][0];
                jb = pent[j][1];
                jc = pent[j][2];
                jd = pent[j][3];
                je = pent[j][4];

                int aflag = 0;
                int bflag = 0;
                int cflag = 0;
                int dflag = 0;
                int eflag = 0;

                bool facet_abcd_found = false;
                bool facet_acde_found = false;
                bool facet_abde_found = false;
                bool facet_abce_found = false;
                bool facet_bcde_found = false;

                // check which nodes on i are found on this second pent j
                if((ia == ja) || (ia == jb) || (ia == jc) || (ia == jd) || (ia == je))
                {
                    aflag = 1;
                }

                if((ib == ja) || (ib == jb) || (ib == jc) || (ib == jd) || (ib == je))
                {
                    bflag = 1;
                }

                if((ic == ja) || (ic == jb) || (ic == jc) || (ic == jd) || (ic == je))
                {
                    cflag = 1;
                }

                if((id == ja) || (id == jb) || (id == jc) || (id == jd) || (id == je))
                {
                    dflag = 1;
                }

                if((ie == ja) || (ie == jb) || (ie == jc) || (ie == jd) || (ie == je))
                {
                    eflag = 1;
                }

                
                if((aflag + bflag + cflag + dflag) == 4)
                {
                    facet_abcd_found = true;
                }

                if((aflag + cflag + dflag + eflag) == 4)
                {
                    facet_acde_found = true;
                }

                if((aflag + bflag + dflag + eflag) == 4)
                {
                    facet_abde_found = true;
                }

                if((aflag + bflag + cflag + eflag) == 4)
                {
                    facet_abce_found = true;
                }

                if((bflag + cflag + dflag + eflag) == 4)
                {
                    facet_bcde_found = true;
                }

                // at this point, we now would have identified all of the boundary facets that could have been found

                if(facet_abcd_found == false)
                {
                    boundary_facets[boundary_facet_counter][0] = ia;
                    boundary_facets[boundary_facet_counter][1] = ib;
                    boundary_facets[boundary_facet_counter][2] = ic;
                    boundary_facets[boundary_facet_counter][3] = id;
                    
                    boundary_facet_2_cell[boundary_facet_counter] = i;
                    
                    boundary_facet_counter ++;
                }

                if(facet_acde_found == false)
                {
                    boundary_facets[boundary_facet_counter][0] = ia;
                    boundary_facets[boundary_facet_counter][1] = ic;
                    boundary_facets[boundary_facet_counter][2] = id;
                    boundary_facets[boundary_facet_counter][3] = ie;
                    
                    boundary_facet_2_cell[boundary_facet_counter] = i;
                    
                    boundary_facet_counter ++;
                }

                if(facet_abde_found == false)
                {
                    boundary_facets[boundary_facet_counter][0] = ia;
                    boundary_facets[boundary_facet_counter][1] = ib;
                    boundary_facets[boundary_facet_counter][2] = id;
                    boundary_facets[boundary_facet_counter][3] = ie;
                    
                    boundary_facet_2_cell[boundary_facet_counter] = i;
                    
                    boundary_facet_counter ++;
                }

                if(facet_acde_found == false)
                {
                    boundary_facets[boundary_facet_counter][0] = ia;
                    boundary_facets[boundary_facet_counter][1] = ic;
                    boundary_facets[boundary_facet_counter][2] = id;
                    boundary_facets[boundary_facet_counter][3] = ie;
                    
                    boundary_facet_2_cell[boundary_facet_counter] = i;
                    
                    boundary_facet_counter ++;
                }

                if(facet_bcde_found == false)
                {
                    boundary_facets[boundary_facet_counter][0] = ib;
                    boundary_facets[boundary_facet_counter][1] = ic;
                    boundary_facets[boundary_facet_counter][2] = id;
                    boundary_facets[boundary_facet_counter][3] = ie;
                    
                    boundary_facet_2_cell[boundary_facet_counter] = i;
                    
                    boundary_facet_counter ++;
                }
            }
        }
    }

    // now we have found all the boundary facets, so we need to loop over all of those and compute the normal vectors orthogonal to the facets (nB)

    std::vector < std::array <int,4> > boundary_normal;
    boundary_normal.reserve(boundary_facet_counter);

    for(int i = 0; i < boundary_facet_counter; i++)
    {
        // identify which of the 5 nodes on index_cell is NOT on the boundary facet. this one is called v4

        int v0,v1,v2,v3,v4;
        v0 = boundary_facets[i][0];
        v1 = boundary_facets[i][1];
        v2 = boundary_facets[i][2];
        v3 = boundary_facets[i][3];

        int index_cell = boundary_facet_2_cell[i]; // the cell that this boundary facet is attached to

        int i0, i1, i2, i3, i4;
        
        i0 = pent[index_cell][0];
        i1 = pent[index_cell][1];
        i2 = pent[index_cell][2];
        i3 = pent[index_cell][3];
        i4 = pent[index_cell][4];

        if((i0 != v0) && (i0 != v1) && (i0 != v2) && (i0 != v3))
        {
            v4 = i0;
        }

        if((i1 != v0) && (i1 != v1) && (i1 != v2) && (i1 != v3))
        {
            v4 = i1;
        }

        if((i2 != v0) && (i2 != v1) && (i2 != v2) && (i2 != v3))
        {
            v4 = i2;
        }

        if((i3 != v0) && (i3 != v1) && (i3 != v2) && (i3 != v3))
        {
            v4 = i3;
        }

        if((i4 != v0) && (i4 != v1) && (i4 != v2) && (i4 != v3))
        {
            v4 = i4;
        }

        // now we know that v4 is the node number that isn't on the boundary facet. The next step is to compute the 3-cross product between
        // the 4 nodes that are on the boundary.

        double ax, ay, az, aw;
        double bx, by, bz, bw;
        double cx, cy, cz, cw;
        double dx, dy, dz, dw;

        ax = x[v0];
        ay = y[v0];
        az = z[v0];
        aw = w[v0];

        bx = x[v1];
        by = y[v1];
        bz = z[v1];
        bw = w[v1];

        cx = x[v2];
        cy = y[v2];
        cz = z[v2];
        cw = w[v2];

        dx = x[v3];
        dy = y[v3];
        dz = z[v3];
        dw = w[v3];

        double bax, bay, baz, baw;
        double cax, cay, caz, caw;
        double dax, day, daz, daw;

        bax = bx-ax;
        bay = by-ay;
        baz = bz-az;
        baw = bw-aw;

        cax = cx-ax;
        cay = cy-ay;
        caz = cz-az;
        caw = cw-aw;

        dax = dx-ax;
        day = dy-ay;
        daz = dz-az;
        daw = dw-aw;

        // now we need to take the cross product of these 3 vectors

        double i_comp = three_determinant(bay,baz,baw,cay,caz,caw,day,daz,daw)*(1./6.);
        double j_comp = -three_determinant(bax,baz,baw,cax,caz,caw,dax,daz,daw)*(1./6.);
        double k_comp = three_determinant(bax,bay,baw,cax,cay,caw,dax,day,daw)*(1./6.);
        double l_comp = -three_determinant(bax,bay,baz,cax,cay,caz,dax,day,daz)*(1./6.);

        // now need to check orientation, but to do that we need to find the vector from v4 to the centroid of the boundary facet

        double x_cent, y_cent, z_cent, w_cent;

        x_cent = (ax + bx + cx + dx)*(1./4.);
        y_cent = (ay + by + cy + dy)*(1./4.);
        z_cent = (az + bz + cz + dz)*(1./4.);
        w_cent = (aw + bw + cw + dw)*(1./4.);

        double dot_prod = ( (x_cent-x[v4]) *i_comp) + ( (y_cent-y[v4]) *j_comp) + ( (z_cent-z[v4]) *k_comp) + ( (w_cent-w[v4]) *l_comp);

        if(dot_prod < 0) // flip it if it's not oriented properly
        {
            i_comp = -i_comp;
            j_comp = -j_comp;
            k_comp = -k_comp;
            l_comp = -l_comp;
        }

        boundary_normal[i][0] = i_comp;
        boundary_normal[i][1] = j_comp;
        boundary_normal[i][2] = k_comp;
        boundary_normal[i][3] = l_comp;

        // now we need to identify which edges on index_cell are on the boundary facet because those are the ones onto which we need to
        // add this correction

        for(int j = 0; j < 10; j++)
        {
            int index_edge = cell2edge[index_cell][j];
            int n0 = global_edges[index_edge][0];
            int n1 = global_edges[index_edge][1];

            // if both the nodes on this edge are on the boundary facet, the edge is on the boundary facet and we must update

            if(( (n0 == v0) || (n0 == v1) || (n0 == v2) || (n0 = v3) ) && ( (n1 == v0) || (n1 == v1) || (n1 == v2) || (n1 == v3) ) )
            {
                for(int q = 0; q < 4; q++)
                {
                    n01[index_edge][q] += boundary_normal[i][q] * (1./20.);
                }
            }
        }
    }

    // check boundary information

    std::vector<double> boundary_residual = {0.0, 0.0, 0.0, 0.0};
    double boundary_hyperarea = 0;

    for(int i = 0; i < boundary_facet_counter; i++)
    {
        double mag = 0;
        for(int j = 0; j<4; j++)
        {
            boundary_residual[j] += boundary_normal[i][j];
            mag += boundary_normal[i][j] * boundary_normal[i][j];
        }

        mag = sqrt(mag);

        boundary_hyperarea += mag;
    }

    std::cout << "Boundary Hyperarea: " << std::setprecision(16) << boundary_hyperarea << std::endl;

    std::cout << "Boundary Residual: " << boundary_residual[0] << " " << boundary_residual[1] << " " << boundary_residual[2] << " " << boundary_residual[3] << std::endl;
}
