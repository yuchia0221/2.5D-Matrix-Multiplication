## 2.5D Matrix Multiplication

#### Table of Contents

-   [Performance](#performance)
-   [Folder Structure](#folder-structure)
-   [Commands](#commands)


---

### Performance
![performance](/images/result.png)

### Folder Structure

    .
    ├── multiply.cpp            # Parallelized Version
    ├── Makefile                # Recipes for building and running your program
    └── README.md

## Commands

Makefile:

> a Makefile that includes recipes for building and running your program.

```bash
make                # builds your code
make view           # views an execution trace with Jumpshot
make run-hpc        # creates a HPCToolkit database for performance measurements
make clean          # removes all executable files
make clean-hpc      # removes all HPCToolkit-related files
```

To run an MPI program on the cluster with the default settings
(number of ranks = number of nodes x ntasks per node):

```bash
srun ./multiply <matrix_size> <c>
```

To run an MPI program on the cluster with a specified number of ranks
(which must be <= number of nodes x ntasks per node):

```bash
srun -n <numranks> ./multiply <matrix_size> <c>
```
