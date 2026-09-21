import math


# reading a sample mesh
text = open("elms_orientation_test.txt","r")
lines = text.readlines()
meshing = []
for line in lines:
    line_numbers = line.split()
    nodes = line_numbers[-4:]
    for i in range(len(nodes)):
        nodes[i] = int(nodes[i])
    meshing.append(nodes)

#meshing = [[1,3,4,5],[2,6,5,4],[2,9,10,3],[7,8,11,4]] # sample mesh connectivity array
n_tet = len(meshing) # number of tets in the mesh
n_nodes = max(max(meshing)) # computing the total number of nodes by finding largest node number

neighbor_nodes = [[] for i in range(n_nodes)]


# for every node that exists in our meshing
for node in range(n_nodes):
    # go through row by row in the meshing array
    for row in range(len(meshing)):
        # go through each node in that row
        for item in range(len(meshing[1])):
            # if our current node is in the row, add the row to the node's neighbor list
            if meshing[row][item] == node+1:
                # once we identify that the node is in the row, loop through each item in the row
                for i in range(len(meshing[row])):
                     # as long as the item list isn't equal to the node number itself, add it to the list of neighbors
                    if meshing[row][i] != node+1:
                        neighbor_nodes[node].append(meshing[row][i]);

# figure out how to connect all the edges together on our shape
current_edge = 1;
edge_list = [];
for i in range(n_nodes):
    for k in neighbor_nodes[i]:
        if i+1 < k:
            current_edge = [i+1,k]
            if current_edge not in edge_list:
                edge_list.append([i+1,k])
# the ith entry in edge_list is the ith edge number, and the 2 numbers in that entry are the nodes that bound the edge
# if the 4th entry is [2,6], that means that the 4th global edge number is bounded by nodes 2 and 6


# cell2edge has the following format:
# cell 1 [bounding_edge_#1, bounding_edge_#2, bounding_edge_#3, ...]
# cell 2 [bounding_edge_#1, bounding_edge_#2, bounding_edge_#3, ...]
# cell 3 [bounding_edge_#1, bounding_edge_#2, bounding_edge_#3, ...]
edge_number = 0;
cell2edge = [[] for i in range(len(meshing))]

for k in range(len(meshing)): # looping through each cell
    for i in range(len(edge_list)): # looping through the list of edges
        if (edge_list[i][0] in meshing[k]) & (edge_list[i][1] in meshing[k]): # if both bounding nodes also bound a cell...
            cell2edge[k].append(i+1) # that edge is a bounding edge of the cell.
print(cell2edge)

# note: this is using 1 as the first index rather than 0 ^^

# notes from meeting 12/20/2024

# for each edge, we need a vector (directed area vector) and have an association between each edge and each vertex
# array that contains all the edges and the directed area vectors along the vectors *********
# create data structures from Hiro video

# either use the tetrahedra as the CVs, or build "dual-volumes" that surround each node itself. Just different ways to
# preserve the same conservation laws (energy, mass, momentum, etc)
# primary = tetrahedra, dual = around each node
# need to measure the flux through the surface of the CV (dual volume)
#   every time the edge goes through the surface, just take the sum of the vectors associated with each of the DV face
#   and that gives you the directed area vector
#   only add the ones that intersect with the edge that sticks out
#   but how do we get the vectors to add together? theres a formula
#   normal vectors from j to the face opposite of j and then use formula 15