#include "mex.h"
#include <math.h>
#include <string.h>

// Function to perform a binary search of the distance among distance bins.
int binary_search(double distance, double distance_bin_edges[], int N_edges)
{
    // Initialize indices
    int mid_ind;
    int start_ind = 0;
    int end_ind = N_edges;

    // Return out-of-bound index if the distance is larger than the largest bin edge.
    if (distance < distance_bin_edges[0]){
        return 0;
    }
    if (distance > distance_bin_edges[N_edges-1]){
        return N_edges;
    }
    // Perform the binary search.
    while (start_ind + 1 < end_ind){
        mid_ind = (start_ind + end_ind)/2;
        if (distance <= distance_bin_edges[mid_ind])
            end_ind = mid_ind;
        else {
            start_ind = mid_ind;
            }
//        printf("%d %d\n",start_ind,end_ind);
    }
    return end_ind;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *datap,*BinEdges_p,Binwidth_log,temp;
    double *XI, *XJ, *XI0;
    double squareddist,Y;
    int N_data,N_bins,i,j,k,l,MidBinInd;
    int N_dim,N_points;

    /*  get the dimensions of the matrix input y */
    N_dim = mxGetM(prhs[0]);
    N_points = mxGetN(prhs[0]);

    // Get the data.
    XI = mxGetPr(prhs[0]);
    N_data = mxGetNumberOfElements(prhs[0]);

    // Initialize the bin counts and square the bin edges for faster binning.
    N_bins = mxGetNumberOfElements(prhs[1]);
    double BinEdges_squared[N_bins];
    double BinCounts[N_bins];
    BinEdges_p = mxGetPr(prhs[1]);
    for (i = 0; i < N_bins; i++) {
        BinEdges_squared[i] = pow(*(BinEdges_p + i),2);
        BinCounts[i] = 0;
    }

    // Determine if periodic boundary conditions should be used.
    int periodic = nrhs > 2;
    double Dim_sizes[N_dim];

    if (periodic) {
        // Check the shape of the input boundaries.
        int N_boundaries = mxGetM(prhs[2]);
        int N_boundaries_dim = mxGetN(prhs[2]);
        if (N_boundaries != 2 || N_boundaries_dim !=N_dim) {
            mexErrMsgIdAndTxt("corrdim2:Boundaries",
                    "The shape of the input boundaries is incorrect. A shape of 2 x n is expected.");
        }

        // Copy the input boundaries into an array.
        //double boundaries[2][N_dim] = (double) (*(mxGetPr(prhs[2])));
        double boundaries[2][N_dim];
        double *boundaries_p = &boundaries[0][0];
        double *input_boundaries_p = mxGetPr(prhs[2]);

        for (i=0; i<N_boundaries; i++) {
            for (j=0; j<N_dim; j++){
                boundaries[i][j] = *(input_boundaries_p + i + j*N_boundaries);
            }
        }

        // Calculate the size of each dimension.
        for (i=0; i<N_dim; i++) {
            Dim_sizes[i] = boundaries[1][i] - boundaries[0][i];
            if (Dim_sizes[i] <=0 ){
                char err_msg[75];
                snprintf(err_msg,sizeof(err_msg),"The size of the boundary in dimension %i is negative.",i+1);
                mexErrMsgIdAndTxt("corrdim2:Boundary_sizes",err_msg);
            }
            //printf("%f\n", Dim_sizes[i]);
        }
    }  
    
    //Binwidth_log=log(BinEdges_p[N_bins-1]/BinEdges_p[0])/(N_bins-1);
    MidBinInd = (int) (N_bins/2);

    for (i=0; i<N_points; i++) {
            XI0 = XI;
            XJ = XI+N_dim;
            for (j=i+1; j<N_points; j++) {
                /* XI = x + i*n; XJ = x + j*n; */
                XI = XI0;
                squareddist = 0;
                for (k=0; k<N_dim; k++,XI++,XJ++){
                    // Calculate the squared distance between x_i and x_j in the kth dimension.
                    Y = fabs((*XI)-(*XJ));

                    // If periodic boundary conditions are used,
                    // find the smallest distance between the two points in the kth dimension.
                    if (periodic && Y > Dim_sizes[k]/2){
                        Y = Dim_sizes[k] - Y;
                    }
                    squareddist += Y*Y;
                }

                // Find a bin search start index that is closest to the squared distance value.
//                 l = 0;
//                 l = (squareddist > BinEdges_squared[MidBinInd]) ? MidBinInd : l;
//                 l = (squareddist > BinEdges_squared[N_bins-1]) ? N_bins : l;
                
                // Find where the largest bin edge that is smaller or equal than the distance between XI and XJ.
                l = binary_search(squareddist,BinEdges_squared,N_bins);
                
//                 // Scan across the bins to find the bin index.
//                 while (l < N_bins && squareddist >= BinEdges_squared[l]) {
//                     l++;
//                 }
                if (l < N_bins){
                    BinCounts[l]++;
                }

//                 while (l < N_bins){
//                     //printf("%f\n",*(BinEdges_p+l));
//                     if(squareddist < BinEdges_squared[l])){
//                         BinCounts[l]++;
//                         l=N_bins;
//                     }
//                     l++;
//                 }
            }
        }

    // Format the output.
    double *BinCounts_output_p;
    plhs[0] = mxCreateDoubleMatrix(N_bins,1,mxREAL);
    BinCounts_output_p = mxGetPr(plhs[0]);
    for (i=0; i<N_bins; i++) {
        *(BinCounts_output_p + i) = (double) (BinCounts[i]);
    }
}