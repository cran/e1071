int e1071_floyd(int *n, double *A, double *C, int *P)
/* this function takes an nxn matrix C of edge costs and produces */
/* an nxn matrix A of lengths of shortest paths, and an nxn       */
/* matrix P giving a point in the middle of each shortest path    */
{
    int i,j,k;

    for (i=0; i<*n; i++)
        for (j=0; j<*n; j++)
        {
	    A[i + *n * j] = C[i + *n * j];
	    P[i + *n * j] = -1;
        }
    for (i=0; i<*n; i++)
        A[i + *n * i] = 0;              /* no self cycle */
    for (k=0; k<*n; k++)
        for (i=0; i<*n; i++)
	    for (j=0; j<*n; j++)
		if (A[i + *n * k]+A[k + *n * j] < A[i + *n * j])
		{
		    A[i + *n * j] = A[i + *n * k] + A[k + *n * j];
		    P[i + *n * j] = k;  /* k is included in shortest path */
		}
    return 0;
}

