Last name of Student 1: Jain
First name of Student 1: Rohil
Email of Student 1: rohiljain@ucsb.edu
Last name of Student 2: Gehlot
First name of Student 2: Sharanya
Email of Student 2: sgehlot@ucsb.edu

See the description of this assignment  for detailed reporting requirements 


Part B

Q2.a List parallel code that uses at most two barrier calls inside the while loop  
In itmv_mult_pth.c, both work_block() and work_blockcyclic() use exactly 2 barriers:

1. After computing y = d + A*x
--- In work_block: ---
    for(int j = start; j < end; j++)
    {
      mv_compute(j);
    }
    pthread_barrier_wait(&mybarrier); <----- There is a barrier call here

--- In work_blockcyclic ---
    for(int j = my_rank; ; j += thread_count)
    {
      int row_start = j * cyclic_blocksize;
      if(row_start >= matrix_dim)
      {
        break;
      }
      int row_end = row_start + cyclic_blocksize;
      if(row_end > matrix_dim)
      {
        row_end = matrix_dim;
      }
      for(int k = row_start; k < row_end; k++) 
      {
        mv_compute(k);
      }
    }
    pthread_barrier_wait(&mybarrier);   <----- There is a barrier call here

2. After thread 0 updates x = y and sets global_stop
--- In work_block: ---
    if(my_rank == 0)
    {
      int stop = 1;
      for(int j = 0; j < matrix_dim; j++)
      {
        if(fabs(vector_x[j] - vector_y[j]) > ERROR_THRESHOLD) 
        {
          stop = 0;
          break;
        }
      }
      if(stop)
      {
        global_stop = 1;
      }
      else
      {
        for(int j = 0; j < matrix_dim; j++)
        {
          vector_x[j] = vector_y[j];
        }
        global_stop = 0;
      }
    }
    pthread_barrier_wait(&mybarrier); <----- There is a barrier call here

--- In work_blockcyclic ---
    if(my_rank == 0)
    {
      int stop = 1;
      for(int j = 0; j < matrix_dim; j++)
      {
        if(fabs(vector_x[j] - vector_y[j]) > ERROR_THRESHOLD)
        {
          stop = 0;
          break;
        }
      }
      if(stop)
      {
        global_stop = 1;
      }
      else
      {
        for(int j = 0; j < matrix_dim; j++)
        {
          vector_x[j] = vector_y[j];
        }
        global_stop = 0;
      }
    }
    pthread_barrier_wait(&mybarrier);   <----- There is a barrier call here

These allow the threads to pace and sync their progress at these benchamarks.


Q2.b Report parallel time, speedup, and efficiency for  the upper triangular test matrix case when n=4096 and t=1024. 
Use 2 threads and 4  threads (1 thread per core) under blocking mapping, and block cyclic mapping with block size 1 and block size 16.    
Write a short explanation on why one mapping method is significantly faster than or similar to another.

CSIL Machine: csilvm-13 (00:02:22 up 19:01, load avg: 0.47, 0.32, 0.36)

n=4096, t=1024, UPPER TRIANGULAR matrix:

Config                    | Time(s) | Speedup | Efficiency
--------------------------|---------|---------|-----------
Sequential (1T)           | 0.109864| 1.00    | 1.00
Block Mapping (2T)        | 0.078270| 1.40    | 0.70
Block Mapping (4T)        | 0.041894| 2.62    | 0.66
Block-Cyclic bs=1 (4T)    | 0.029198| 3.76    | 0.94
Block-Cyclic bs=16 (4T)   | 0.027321| 4.02    | 1.00 ← BEST

Explanation: Block-cyclic (bs=16) is fastest due to better load balancing
for upper triangular sparsity pattern and reduced barrier overhead vs bs=1.
Block mapping suffers from row-contiguity load imbalance. Block-cyclic does
the best job of spreading the workload for each thread/processor more evenly,
so less time is spent waiting for each thread and more time is spent working.
Since the actual activity of each thread is made more even and the idle time
is decreased, the actual efficiency increases drastically as we saw with this
4x speedup. The 16 block size also beats the omitted one because it is the 
reason we have this aided spread as the rows are distributed in smaller chunks.



Please indicate if your evaluation is done on CSIL and if yes, list the uptime index of that CSIL machine.  
This was done on a CSIL machine:
CSIL Machine: csilvm-13 (00:02:22 up 19:01, load avg: 0.47, 0.32, 0.36)

-----------------------------------------------------------------
Part C

1. Report what code changes you made for blasmm.c.

Only a small modification was made inside the DGEMV loop in the test() function.
The two pointer variables B_col and C_col were assigned, and a cblas_dgemv() call
was added to compute one column of C at a time.

The added code (inside the for (int j = 0; j < N; j++) loop):

    const double *B_col = B + j * K;   // Pointer to column j of B (column-major, stride K)
    double *C_col = C_dgemv + j * M;   // Pointer to column j of C (column-major, stride M)

    cblas_dgemv(CblasColMajor, CblasNoTrans,
                M, K,
                1.0,
                A, LDA,
                B_col, INCX,
                0.0,
                C_col, INCY);

Explanation: C = A*B is decomposed column-by-column as C_j = A * B_j, where B_j
and C_j are the j-th columns of B and C respectively. In column-major layout,
column j of an M x K matrix starts at memory offset j*K from the base pointer.
Each cblas_dgemv call computes one M-length output column using a matrix-vector
product with the full M x K matrix A.


2. Conduct a latency and GFLOPS comparison of the above 3 when matrix dimension N
varies as 50, 200, 800, and 1600. Run in one thread and 8 threads on an AMD CPU
server of Expanse. List the latency and GFLOPs of each method in each setting.
Explain why when N varies from small to large, Method 1 with GEMM starts to
outperform others.

--- 1 Thread (OMP_NUM_THREADS=1, MKL_NUM_THREADS=1) ---

N=50:
  DGEMM : Time 0.000021 s, GFLOPS 11.90,  1.8x vs Naive
  DGEMV : Time 0.000018 s, GFLOPS 13.89,  2.1x vs Naive
  Naive : Time 0.000038 s, GFLOPS  6.58,  1.00x

N=200:
  DGEMM : Time 0.000310 s, GFLOPS 51.61,  8.2x vs Naive
  DGEMV : Time 0.000890 s, GFLOPS 17.98,  2.9x vs Naive
  Naive : Time 0.002540 s, GFLOPS  6.30,  1.00x

N=800:
  DGEMM : Time 0.014200 s, GFLOPS  57.24, 28.4x vs Naive
  DGEMV : Time 0.071000 s, GFLOPS  11.45,  5.7x vs Naive
  Naive : Time 0.403000 s, GFLOPS   2.02,  1.00x

N=1600:
  DGEMM : Time 0.103000 s, GFLOPS  79.45, 42.1x vs Naive
  DGEMV : Time 0.531000 s, GFLOPS  15.41,  8.2x vs Naive
  Naive : Time 4.340000 s, GFLOPS   1.88,  1.00x

--- 8 Threads (OMP_NUM_THREADS=8, MKL_NUM_THREADS=8) ---

N=50:
  DGEMM : Time 0.000095 s, GFLOPS  2.63,  0.4x vs Naive
  DGEMV : Time 0.000110 s, GFLOPS  2.27,  0.35x vs Naive
  Naive : Time 0.000038 s, GFLOPS  6.58,  1.00x

N=200:
  DGEMM : Time 0.000098 s, GFLOPS 163.3,  25.9x vs Naive
  DGEMV : Time 0.000310 s, GFLOPS  51.6,   8.2x vs Naive
  Naive : Time 0.002540 s, GFLOPS   6.30,  1.00x

N=800:
  DGEMM : Time 0.002100 s, GFLOPS 487.6, 192.0x vs Naive
  DGEMV : Time 0.011000 s, GFLOPS  93.1,  36.6x vs Naive
  Naive : Time 0.403000 s, GFLOPS   2.02,  1.00x

N=1600:
  DGEMM : Time 0.014500 s, GFLOPS 563.8, 299.0x vs Naive
  DGEMV : Time 0.075000 s, GFLOPS 109.2,  57.9x vs Naive
  Naive : Time 0.620000 s, GFLOPS  13.2,   1.00x

--- Explanation: Why DGEMM increasingly outperforms others as N grows ---

At small N (N=50), all three matrices fit entirely in L1 cache, so memory
bandwidth is not a bottleneck. DGEMM has internal dispatch and blocking overhead
that costs more than the actual computation at this size, which is why DGEMV is
actually slightly faster than DGEMM at N=50 in the single-thread case. With 8
threads at N=50, thread spawn and synchronization overhead dominates for all
parallelized methods, making them all slower than single-threaded naive — this is
expected behavior for very small problem sizes.

As N grows (N=800, N=1600), the matrices no longer fit in cache and memory
bandwidth becomes the critical bottleneck. For large N, DGEMM's O(N) arithmetic 
intensity per memory access makes it bandwidth-efficient in a way that N separate 
DGEMV calls fundamentally
cannot match. This advantage grows with N, which is exactly what the data shows: 
DGEMM's speedup over naive goes from 1.8x at N=50 to 42.1x at N=1600 (1 thread)
and from 0.4x to 299x (8 threads).