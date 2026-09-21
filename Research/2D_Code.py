import math

# Residual Computation Function

def residual_computation(n_nodes,n_edges,x,y,neighbor_nodes,coefficient_vector,edge_list,n_jk_vectors,dual_hypervolumes,U,S,residual):

    for i in range(n_nodes):
        residual[i] = 0

    ## computing gradients
    # gradient is computed as a collection of the gradient values from all of node i's neighbors
    gradient = [[0,0] for i in range(n_nodes)]

    for i in range(n_nodes):

        number_of_neighbors_local = len(neighbor_nodes[i])
        
        for j in range(number_of_neighbors_local):
            
            index_neighbor = neighbor_nodes[i][j]
            
            for k in range(2):
                gradient[i][k] += coefficient_vector[i][j][k]*(U[index_neighbor]-U[i]) 
    
    ## edge loop on residuals
    for i in range(n_edges):

        n0 = edge_list[i][0]
        n1 = edge_list[i][1]

        U0 = U[n0]
        U1 = U[n1]

        n01_unit = [0,0]
        mag = 0

        for j in range(2):
            mag += n_jk_vectors[i][j]*n_jk_vectors[i][j]
        mag = math.sqrt(mag)

        for j in range(2):
            n01_unit[j] = n_jk_vectors[i][j]/mag  # unit vector of n_jk for this edge
        
        advection_direction = [1.,1.]
        projected_direction = 0

        for j in range(2):
            projected_direction += advection_direction[j]*n01_unit[j]

        UL = U0 + gradient[n0][0]*(x[n1]-x[n0])*0.5 + gradient[n0][1]*(y[n1]-y[n0])*0.5
        UR = U1 - gradient[n1][0]*(x[n1]-x[n0])*0.5 - gradient[n1][1]*(y[n1]-y[n0])*0.5

        flux_temp = (0.5*projected_direction*(UL + UR) - 0.5*abs(projected_direction)*(UR - UL))*mag

        residual[n0] += flux_temp
        residual[n1] -= flux_temp

    ## node loop (forcing term)

    for i in range(n_nodes):
        residual[i] -= S[i]*dual_hypervolumes[i]

    #for i in range(n_nodes):
        #for k in range(2):
            #xval = x[i]
            #yval = y[i]
            #ethresh = 0.0000000001
            #xmax = 1.0
            #xmin = 0.0
            #ymax = 1.0
            #ymin = 0.0
            #if abs(xval-xmax) < ethresh or abs(xval-xmin) < ethresh or abs(yval-ymax) < ethresh or abs(yval-ymin) < ethresh:
                #pass
            #else:
                #print(gradient[i][k])
    return


# MAIN CODE

# reading a sample mesh
text = open("mesh5_connectivity.txt","r")
lines = text.readlines()
meshing = []
for line in lines:
    nodes = line.split()
    for i in range(len(nodes)):
        nodes[i] = int(nodes[i])
    meshing.append(nodes)
#meshing = [[1,3,4,5],[2,6,5,4],[2,9,10,3],[7,8,11,4]] # sample mesh connectivity array
n_tet = len(meshing) # number of tets in the mesh

# computing the number of nodes
n_nodes = 0
for rows in meshing:
    for node in rows:
        if node > n_nodes:
            n_nodes = node
n_nodes += 1
# added 1 because the node count starts at 0

# the nth list in neighbor_nodes contains the nth node's neighbors
neighbor_nodes = [[] for i in range(n_nodes)]
for row in meshing:
    # For each node in the current row
    for node in row:
        # Add all other nodes in the row to the current node's neighbor list
        for neighbor in row:
            if neighbor != node and neighbor not in neighbor_nodes[node]:
                neighbor_nodes[node].append(neighbor)

# figure out how to connect all the edges together on our shape
current_edge = 1;
edge_list = [];
for i in range(n_nodes):
    for k in neighbor_nodes[i]:
        if i < k:
            current_edge = [i,k]
            if current_edge not in edge_list:
                edge_list.append([i,k])
# the ith entry in edge_list is the ith edge number, and the 2 numbers in that entry are the nodes that bound the edge
# if the 4th entry is [2,6], that means that the 4th global edge number is bounded by nodes 2 and 6


# cell2edge has the following format:
# cell 1 [bounding_edge_#1, bounding_edge_#2, bounding_edge_#3, ...]
# cell 2 [bounding_edge_#1, bounding_edge_#2, bounding_edge_#3, ...]
# cell 3 [bounding_edge_#1, bounding_edge_#2, bounding_edge_#3, ...]
edge_number = 0;
cell2edge = [[] for i in range(len(meshing))]

for k in range(len(meshing)): # looping through each cell
    for edges in edge_list: # looping through the list of edges
        if (edges[0] in meshing[k]) and (edges[1] in meshing[k]): # if both bounding nodes also bound a cell...
            cell2edge[k].append(edge_list.index(edges)) # that edge is a bounding edge of the cell.

# Now, bring in the coordinates for the connectivity matrix

text = open("mesh5_coords.txt","r")
lines = text.readlines()
coordinates = []
for line in lines:
    xy_coords = line.split()
    coordinates.append(xy_coords)
    for i in range(len(xy_coords)):
        xy_coords[i] = float(xy_coords[i])

coordinates = [node_coords[:-1] for node_coords in coordinates]

# now we have:
# meshing: the jth item in meshing contains the nodes that bound the jth cell
# neighbor_nodes: the jth item in neighbor_nodes contains the jth node's neighbors
# edge_list: the jth item in edge_list contains the nodes that connect the jth edge
# cell2edge: the jth item in cell2edge contains the edges that bound the jth cell
# coordinates: the jth item in coordinates contains the coordinates of the jth cell

n_jk_vectors = []
for edge_number in range(len(edge_list)):
    edge = edge_list[edge_number]
    current_node = edge[0]
    current_neighbor = edge[1]
    # step 1: find the vector for the edge that connects j and k
    j_to_k = []
    j_to_k.append(coordinates[current_neighbor][0] - coordinates[current_node][0])
    j_to_k.append(coordinates[current_neighbor][1] - coordinates[current_node][1])

    # step 2: identify which cells border current_edge
    current_cells = []
    for cells in cell2edge:
        if edge_number in cells:
            current_cells.append(cell2edge.index(cells))

    # NEED TO ACCOUNT FOR BORDERS WHERE THERE'S ONLY 1 CELL TOUCHING AN EDGE
    if len(current_cells) == 1: # if there's only 1 cell touching the edge (it's an edge-case)
        # THIS IS AN EDGE CASE AND WE NEED TO USE THE CORRECTION 

        border_cell = current_cells[0]
        border_nodes = meshing[border_cell][:] # making a shallow copy
        border_nodes.remove(current_node)
        
        # this is the vector that points from current neighbor to the other node attached to the cell
        border_vector = []
        border_node_0 = border_nodes[0]
        border_node_1 = border_nodes[1]
        border_vector.append(coordinates[border_node_0][0] - coordinates[border_node_1][0])
        border_vector.append(coordinates[border_node_0][1] - coordinates[border_node_1][1])

        # rotating it so that it's the normal vector that points away from current_node
        n_jk1 = []
        n_jk1.append(-border_vector[1])
        n_jk1.append(border_vector[0])

        # enforcing that it points away from current_node
        dot_prod = 0
        for i in range(len(n_jk1)):
            dot_prod += n_jk1[i]*j_to_k[i]
        if dot_prod < 0:
            n_jk1 = [-a for a in n_jk1]

        # computing the vector that is orthogonal to j_to_k
        n_jk2 = []
        n_jk2.append(-j_to_k[1])
        n_jk2.append(j_to_k[0])

        # this vector should point away from n_jk1, so we want to make sure the dot product between these vectors is negative
        dot_prod = 0
        for i in range(len(n_jk2)):
            dot_prod += n_jk2[i]*n_jk1[i]
        if dot_prod > 0: # note: this time, we flip the sign if it's > because we want the dot product to be negative
            n_jk2 = [-a for a in n_jk2]

        # now we have n_jk1, which points across the cell, and n_jk2, which points perpendicular to the edge on the end of the mesh
        
        # calculating the n_jk vector for current_node and current_neighbor by doing njk = 1/3(n1+n2/2)
        n_jk = []
        for i in range(len(n_jk2)):
            n_jk.append((n_jk1[i]+n_jk2[i]/2)/3)
        n_jk_vectors.append(n_jk)
    

    else: # it's not an edge case, so compute it the regular way
        # now we have current_cells which are the 2 cells that share the current_edge connecting node j and node k
        left_cell = current_cells[0]
        right_cell = current_cells[1]
        left_nodes = meshing[left_cell][:]  # Make a shallow copy
        right_nodes = meshing[right_cell][:]  # Make a shallow copy

        # for the left and right cells, we just want the nodes that AREN'T the current node, that way we can find the vectors connecting them

        left_nodes.remove(current_node)
        right_nodes.remove(current_node)

        # now left_nodes and right_nodes are just the nodes that bound the edge opposite of the current node
        # left_vector and right_vector are the vectors along the edge opposite of current_node
        left_vector = []
        left_node_0 = left_nodes[0]
        left_node_1 = left_nodes[1]
        left_vector.append(coordinates[left_node_0][0] - coordinates[left_node_1][0])
        left_vector.append(coordinates[left_node_0][1] - coordinates[left_node_1][1])

        right_vector = []
        right_node_0 = right_nodes[0]
        right_node_1 = right_nodes[1]
        right_vector.append(coordinates[right_node_0][0] - coordinates[right_node_1][0])
        right_vector.append(coordinates[right_node_0][1] - coordinates[right_node_1][1])

        # rotating left and right vector by 90 degrees to get njks

        n_jk1 = []
        n_jk1.append(-left_vector[1])
        n_jk1.append(left_vector[0])

        n_jk2 = []
        n_jk2.append(-right_vector[1])
        n_jk2.append(right_vector[0])

        # checking to see if n_jk is oriented correctly and adjusting accordingly
        dot_prod = 0
        for i in range(len(n_jk1)):
            dot_prod += n_jk1[i]*j_to_k[i]
        if dot_prod < 0:
            n_jk1 = [-a for a in n_jk1]

        dot_prod = 0
        for i in range(len(n_jk2)):
            dot_prod += n_jk2[i]*j_to_k[i]
        if dot_prod < 0:
            n_jk2 = [-a for a in n_jk2]

        # calculating the n_jk vector for current_node and current_neighbor by doing njk = 1/3(n1+n2)
        n_jk = []
        for i in range(len(n_jk2)):
            n_jk.append((n_jk1[i]+n_jk2[i])/3)
        n_jk_vectors.append(n_jk)


# n_jk vectors has a vector corresponding to each edge in edge_list. each edge in edge_list is oriented [lower_node_# , higher_node_#]
# so if the 10th entry in edge_list is [15,20]
# and the 10th entry in n_jk_vectors is [0.1,0.1]
# then aj[15] would have the vector [0.1,0.1] added to it, and aj[20] would have the vector [-0.1,-0.1] added to it
# after looping through every edge, all of aj's should be 0
# 
# slight adjustment for boundary nodes. you add 1/(d+1) * the outward facing normal to the edge

# N_JK VALIDATION FOR INTERIOR EDGES

aj = [[0,0] for i in range(n_nodes)] # aj is an array that will eventually contain the sum of the vectors pointing out of each node j
for edge_number in range(len(edge_list)):
    
    edge = edge_list[edge_number]
    current_node = edge[0]
    current_neighbor = edge[1]
    # check how many cells border this edge
    current_cells = []
    for cells in cell2edge:
        if edge_number in cells:
            current_cells.append(cell2edge.index(cells))
    # find the vector for the edge that connects j and k
    j_to_k = []
    j_to_k.append(coordinates[current_neighbor][0] - coordinates[current_node][0])
    j_to_k.append(coordinates[current_neighbor][1] - coordinates[current_node][1])
    
    # CASE 1: INTERIOR EDGES
    if len(current_cells) > 1:
        node1 = edge_list[edge_number][0]
        node2 = edge_list[edge_number][1]
        aj[node1][0] += n_jk_vectors[edge_number][0]
        aj[node1][1] += n_jk_vectors[edge_number][1]

        aj[node2][0] -= n_jk_vectors[edge_number][0]
        aj[node2][1] -= n_jk_vectors[edge_number][1]

    
    # CASE 2: BOUNDARY CELLS
    # in the end, there should be 2 different vectors being added to node 1 AND to node 2: n_jk, and the correction value corresponding to n_B
    #
    else:
        # grab border cell normal vector code from above and use the formula from the paper
        
        # exact same code as interior edges to add on the n_jk that points from j to k
        node1 = edge_list[edge_number][0]
        node2 = edge_list[edge_number][1]
        aj[node1][0] += n_jk_vectors[edge_number][0]
        aj[node1][1] += n_jk_vectors[edge_number][1]

        aj[node2][0] -= n_jk_vectors[edge_number][0]
        aj[node2][1] -= n_jk_vectors[edge_number][1]

        # but now, we also have to incorporate the correction factor (copying code from line 96)
        # literally copied everything just so we can re-calculate the value of n_jk2 (orthogonal to the boundary edge)
        border_cell = current_cells[0]
        border_nodes = meshing[border_cell][:] # making a shallow copy
        border_nodes.remove(current_node)

        #if edge_number == 0:
            #print("the nodes that border node 0 are", border_nodes)

        # this is the vector that points from current neighbor to the other node attached to the cell
        border_vector = []
        border_node_0 = border_nodes[0]
        border_node_1 = border_nodes[1]
        border_vector.append(coordinates[border_node_0][0] - coordinates[border_node_1][0])
        border_vector.append(coordinates[border_node_0][1] - coordinates[border_node_1][1])

        # rotating it so that it's the normal vector that points away from current_node
        n_jk1 = []
        n_jk1.append(-border_vector[1])
        n_jk1.append(border_vector[0])

        # enforcing that it points away from current_node
        dot_prod = 0
        for i in range(len(n_jk1)):
            dot_prod += n_jk1[i]*j_to_k[i]
        if dot_prod < 0:
            n_jk1 = [-a for a in n_jk1]

        # computing the vector that is orthogonal to j_to_k
        n_jk2 = []
        n_jk2.append(-j_to_k[1])
        n_jk2.append(j_to_k[0])


        # this vector should point away from n_jk1, so we want to make sure the dot product between these vectors is negative
        dot_prod = 0
        for i in range(len(n_jk2)):
            dot_prod += n_jk2[i]*n_jk1[i]
        if dot_prod > 0: # note: this time, we flip the sign if it's > because we want the dot product to be negative
            n_jk2 = [-a for a in n_jk2] 
        
        # now that we have n_jk2 again, we can add that to aj to correct the boundary nodes
        
        aj[current_node][0] += n_jk2[0]/2
        aj[current_node][1] += n_jk2[1]/2

        aj[current_neighbor][0] += n_jk2[0]/2
        aj[current_neighbor][1] += n_jk2[1]/2

# cleaning up aj so values that are negligible (xE-10 for example) are rewritten as 0's
for i in range(len(aj)):
    for j in range(len(aj[i])):
        if abs(aj[i][j]) < 0.00001:
            aj[i][j] = 0.0

njk_flag = 0
for i in range(n_nodes):
    for j in range(2):
        if aj[i][j] !=0:
            njk_flag = 1;
            break

if njk_flag == 0:
    print("njk check: PASS")
else:
    print("njk check: FAIL")


# LITLE BIT OF CODE TO CHECK IF A NODE IS ON THE BOUNDARY OR NOT

is_boundary = [0 for i in range(n_nodes)]

for edge_number in range(len(edge_list)):    
    # check how many cells border this edge
    current_cells = []
    for cells in cell2edge:
        if edge_number in cells:
            current_cells.append(cell2edge.index(cells))
    
    if len(current_cells) == 1:
        node1 = edge_list[edge_number][0]
        node2 = edge_list[edge_number][1]
        is_boundary[node1] = 1
        is_boundary[node2] = 1

# COMPUTING DUAL HYPERVOLUMES
# theorem 5 in the paper
## code is based off of the 4D code

# creating new variables x and y to make it more similar to the 4D code provided by Dr. Williams
x = []
y = []

for i in range(n_nodes):
    x.append(coordinates[i][0])
    y.append(coordinates[i][1])

dual_hypervolumes = [0 for i in range(n_nodes)]

# loop through all of the edges
n_edges = len(edge_list)
for i in range(n_edges):
    e0 = edge_list[i][0]
    e1 = edge_list[i][1]

    dual_hypervolumes[e0] += (1/4)*( (x[e1]-x[e0])*n_jk_vectors[i][0]
                                  + (y[e1]-y[e0])*n_jk_vectors[i][1])
    
    dual_hypervolumes[e1] += (1/4)*( (x[e1]-x[e0])*n_jk_vectors[i][0]
                                  + (y[e1]-y[e0])*n_jk_vectors[i][1])

# you can just add to both nodes because both n_jk and the coordinates thing will have flipped signs so it'd be the same

# add them up - it should be 1 because that is the total area of the region
total_hypervolume = 0
for i in range(n_nodes):
    total_hypervolume += dual_hypervolumes[i]

if ( (total_hypervolume > 1.000000001) | (total_hypervolume < 0.99999999) ):
    print("dual hypervolume check: FAIL")
else:
    print("dual hypervolume check: PASS")


# least-squares coefficients (for gradients)

coefficient_vector = [[[0 for k in range(2)] for j in range(1)] for i in range(n_nodes)] # 3-dimensional array that has n_nodes vectors. for each node, we start with 1 list of [0,0]. we'll be adding more based on how many
# neighbors there are

for i in range(n_nodes):
    
    number_of_neighbors_local = len(neighbor_nodes[i]) # how many neighbors this current node has
    
    for j in range(number_of_neighbors_local-1):
        coefficient_vector[i].append([0,0])
    
    # now coefficient_vector has: n_node elements which each contain number_of_neighbors_local copies of [0,0]
    
    diff_coordinate = [[0,0] for j in range(number_of_neighbors_local)]

    for j in range(number_of_neighbors_local):
        index_neighbor = neighbor_nodes[i][j]

        diff_coordinate[j][0] = x[index_neighbor]-x[i]
        diff_coordinate[j][1] = y[index_neighbor]-y[i]

    A_matrix = [[0,0],[0,0]] # initializing the A matrix as a 2x2 (as opposed to the 4x4 that it is in the 4D code)

    for j in range(number_of_neighbors_local):
        for k in range(2):
            for m in range(2):
                A_matrix[m][k] += diff_coordinate[j][m]*diff_coordinate[j][k]

    # compute the inverse of the A_matrix
    det_A = A_matrix[0][0]*A_matrix[1][1]-A_matrix[0][1]*A_matrix[1][0]
    A_matrix_inverse = [[0,0],[0,0]]
    A_matrix_inverse[0][0] += 1/det_A * A_matrix[1][1]
    A_matrix_inverse[0][1] += 1/det_A * -1 * A_matrix[0][1]
    A_matrix_inverse[1][0] += 1/det_A * -1 * A_matrix[1][0]
    A_matrix_inverse[1][1] += 1/det_A * A_matrix[0][0]

    #print(coefficient_vector[i])

    for j in range(number_of_neighbors_local):
        for k in range(2):
            for m in range(2):
                coefficient_vector[i][j][k] += A_matrix_inverse[k][m]*diff_coordinate[j][m]
    #print(coefficient_vector)

## Initialization
U = []
S = []
residual = [0 for i in range(n_nodes)]
Uinit = [0 for i in range(n_nodes)]
kfactor = 0.1
pi = math.pi
for i in range(n_nodes):
    U.append(3.0*x[i]*x[i] + 5.0*y[i]*y[i])
    S.append(6.0*x[i] + 10.0*y[i])

# initial residual calculation
residual_init = [0 for i in range(n_nodes)]
residual_computation(n_nodes,n_edges,x,y,neighbor_nodes,coefficient_vector,edge_list,n_jk_vectors,dual_hypervolumes,U,S,residual_init)

max_residual_init = 0
for i in range(n_nodes):
    if abs(residual_init[i]) > max_residual_init:
        max_residual_init = abs(residual_init[i])

#print("max residual init is ")
#print(max_residual_init)
#print("\n")

number_of_iterations = 1000
CFL = 0.5
residual_threshold = 0.000000000001
local_time_steps = [0 for i in range(n_nodes)]

for m in range(number_of_iterations):
    
    for i in range(n_nodes):
        local_time_steps[i] = 0

    # this part essentially just finds the unit njk vector for each edge and dots that with the advection direction
    for i in range(n_edges):
        n0 = edge_list[i][0]
        n1 = edge_list[i][1]

        n01_unit = [0,0]
        mag = 0
        
        for j in range(2):
            mag += n_jk_vectors[i][j]*n_jk_vectors[i][j]

        mag = math.sqrt(mag)

        for j in range(2):
            n01_unit[j] = n_jk_vectors[i][j]/mag

        advection_direction = [1,1]
        projected_direction = 0

        for j in range(2):
            projected_direction += advection_direction[j] * n01_unit[j]

        local_time_steps[n0] += 0.5 * abs(projected_direction) * mag
        local_time_steps[n1] += 0.5 * abs(projected_direction) * mag

    for i in range(n_nodes):
        local_time_steps[i] = CFL * dual_hypervolumes[i] / local_time_steps[i]

        
    # Stage 1
    for i in range(n_nodes):
        Uinit[i] = U[i]

    residual_computation(n_nodes,n_edges,x,y,neighbor_nodes,coefficient_vector,edge_list,n_jk_vectors,dual_hypervolumes,U,S,residual)

    for i in range(n_nodes):
        xval = x[i]
        yval = y[i]

        ethresh = 0.00000000001
        
        xmax = 1.0
        xmin = 0.0
        ymax = 1.0
        ymin = 0.0

        if abs(xval-xmax) < ethresh or abs(xval-xmin) < ethresh or abs(yval-ymax) < ethresh or abs(yval-ymin) < ethresh: # IF the node isn't on the boundary...
            pass
        else:
            U[i] -= (local_time_steps[i]/dual_hypervolumes[i])*residual[i] # we update the solution at the node
    
    
    # Stage 2

    residual_computation(n_nodes,n_edges,x,y,neighbor_nodes,coefficient_vector,edge_list,n_jk_vectors,dual_hypervolumes,U,S,residual) # now do another pass but with this updated version of U

    for i in range(n_nodes):
        xval = x[i]
        yval = y[i]

        ethresh = 0.00000000001
        
        xmax = 1.0
        xmin = 0.0
        ymax = 1.0
        ymin = 0.0

        if abs(xval-xmax) < ethresh or abs(xval-xmin) < ethresh or abs(yval-ymax) < ethresh or abs(yval-ymin) < ethresh: # IF we're not on the boundary
            pass
        else:
            U[i] = 0.5*(U[i] + Uinit[i]) - 0.5*(local_time_steps[i]/dual_hypervolumes[i])*residual[i] # update the solution at the node
    
    max_residual = -1.

    for i in range(n_nodes):

        if abs(residual[i]) > max_residual:
            xval = x[i]
            yval = y[i]

            ethresh = 0.0000000001
            xmax = 1.0
            xmin = 0.0
            ymax = 1.0
            ymin = 0.0
            if abs(xval-xmax) < ethresh or abs(xval-xmin) < ethresh or abs(yval-ymax) < ethresh or abs(yval-ymin) < ethresh:
                pass
            else:
                max_residual = abs(residual[i])
    
    normalized_max_residual = max_residual/max_residual_init
    residual_threshold = 0.0000000001

    #print("iteration number: ", m)
    #print("norm max residual is: ", normalized_max_residual, ", ")
    #print("max residual is ", max_residual)

    if normalized_max_residual < residual_threshold:
        break


max_error = 0
for i in range(n_nodes):
    exact_sol = 3.0*x[i]*x[i] + 5.0*y[i]*y[i]
    point_error = abs(U[i]-exact_sol)
    if point_error > max_error:
        max_error = point_error
print("\n")
print("maximum error is ", max_error)

with open("mesh5_nonsymmetric_quadratic_solution.txt","w") as file:
    for i in U:
        file.write(str(i) + "\n")
