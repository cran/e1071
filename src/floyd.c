int e1071_floyd(int *n, double A[*n][*n], double C[*n][*n], int P[*n][*n])
/* this function takes an nxn matrix C of edge costs and produces */
/* an nxn matrix A of lengths of shortest paths, and an nxn       */
/* matrix P giving a point in the middle of each shortest path    */
{
    int i,j,k;

    for (i=0; i<*n; i++)
        for (j=0; j<*n; j++)
        {
	    A[i][j] = C[i][j];
	    P[i][j] = -1;
        }
    for (i=0; i<*n; i++)
        A[i][i] = 0;              /* no self cycle */
    for (k=0; k<*n; k++)
        for (i=0; i<*n; i++)
	    for (j=0; j<*n; j++)
		if (A[i][k]+A[k][j] < A[i][j])
		{
		    A[i][j] = A[i][k] + A[k][j];
		    P[i][j] = k;  /* k is included in the shortest path */
		}
    return 0;
}

