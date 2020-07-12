# Subtractive Vertex Magic Labelings

## Centre College CSC 350 Parallel Computing, MAT 490/491 Research in Magic Graphs

## Jordan Turley, Jason Pinto, Matthew Ko, Aaron Davis

This program generates and finds subtractive vertex magic labelings for directed graphs. We are focusing on two identical directed cycles joined by a certain number of vertices in the middle. We have been able to find all labelings for eight graphs and hope to be able to generalize this to any labeling of this style.

Read more about the algorithm and results in [CSC_350_Subtractive_Vertex_Magic_Paper.pdf](https://github.com/ztizzlegaming/parallel-graph-labelings/blob/master/CSC_350_Subtractive_Vertex_Magic_Paper.pdf).

### Usage

First, clone the repository:

    git clone https://github.com/ztizzlegaming/parallel-graph-labelings
    
**Sequential**

    g++ subtractive_vertex_magic.c -O3 -o vertex_magic
    ./vertex_magic 3 2
    ./vertex_magic 3 1
    ./vertex_magic 4 3
    ...
    
**OpenMP**

    g++ subtractive_vertex_magic_openmp.c -o vertex_magic_openmp -O3 -fopenmp
    ./vertex_magic_openmp 3 2
    ./vertex_magic_openmp 3 1
    ./vertex_magic_openmp 4 3
    ...

**MPI**

    mpic++ subtractive_vertex_magic_mpi.c -o vertex_magic_mpi -O3
    mpirun -n [num processes] -machinefile [machinefile] ./vertex_magic_openmp 3 2
    mpirun -n [num processes] -machinefile [machinefile] ./vertex_magic_openmp 3 1
    mpirun -n [num processes] -machinefile [machinefile] ./vertex_magic_openmp 4 3
    ...
