#include <mex.h>
#include "kdtree.h"
#include <string.h>

#include "kdtree.cpp"

void mexFunction(int nlhs, mxArray * plhs[],
		 int nrhs, const mxArray * prhs[])
{

  if(mxIsClass(prhs[0],"kdtree")==0) {
    mexPrintf("First argument must be a kdtree class\n");
    return;
  }
  KDTree *tree = KDTree::unserialize(mxGetPr(mxGetFieldByNumber(prhs[0],0,0)));

  // Verify the point array
  if (mxGetNumberOfDimensions(prhs[1]) != 2) {
    mexPrintf("Invalid point array passed in.\n");
    return;
  }
  int npoints = (int) mxGetM(prhs[1]);
  int ndim = (int) mxGetN(prhs[1]);
  if (ndim != tree->ndim) {
    mexPrintf("Points have wrong number of dimensions.\n");
    mexPrintf("Tree dimension = %d\n",tree->ndim);
    mexPrintf("Input array dimension = %d\n",ndim);
    delete tree;
    return;
  }
  // Check the format of the input
  bool isDouble = false;
  mxClassID id = mxGetClassID(prhs[1]);
  double *dPtr = (double *)0;
  float  *sPtr = (float *)0;
  if (id == mxDOUBLE_CLASS) {
    dPtr = (double *) mxGetPr(prhs[1]);
    isDouble = true;
  } else if (id == mxSINGLE_CLASS) {
    sPtr = (float *) mxGetPr(prhs[1]);
    isDouble = false;
  } else {
    mexPrintf("Input points must be either sinlgle or double\n");
    delete tree;
    return;
  }


  // Create an output array of indices
  //mexPrintf("A.\n");
  plhs[0] = mxCreateNumericMatrix(npoints, 1, mxDOUBLE_CLASS, mxREAL);
  //mexPrintf("A.\n");
  double *idxptr = (double *) mxGetPr(plhs[0]);
  //mexPrintf("A.\n");

  // Create an output array of points, if needed
  float *fpntptr = (float *) 0;
  double *dpntptr = (double *) 0;
  if (nlhs > 1) {
    if (isDouble) {
  //mexPrintf("B.\n");
      plhs[1] =
	  mxCreateNumericMatrix(npoints, ndim, mxDOUBLE_CLASS, mxREAL);
  //mexPrintf("B.\n");
      dpntptr = (double *) mxGetPr(plhs[1]);
  //mexPrintf("B.\n");
    } else {
  //mexPrintf("B.\n");
      plhs[1] =
	  mxCreateNumericMatrix(npoints, ndim, mxSINGLE_CLASS, mxREAL);
  //mexPrintf("B.\n");
      fpntptr = (float *) mxGetPr(plhs[1]);
  //mexPrintf("B.\n");
    }
  }
  // Allocate the point to check
  //mexPrintf("C.\n");
  float *curPoint = new float[ndim];
  //mexPrintf("C.\n");
  for (int i = 0; i < npoints; i++) {
    // Extract the point in the correct format
    if (isDouble) {
      for (int j = 0; j < ndim; j++) {
        // i*ndims+j
        // j*npoints + i;
        // MATLAB stores the transpose of the normal order
         curPoint[j] = (float) dPtr[j * npoints + i];
      }
    } else {
      for (int j = 0; j < ndim; j++) 
        curPoint[j] = sPtr[j * npoints + i];
    }

    // Check the point
//mexPrintf("D.\n");
    int idx;
    tree->closest_point(curPoint, idx);
////mexPrintf("D.\n");
  
    // Set the output arrays
    // First the index
    // adding 1 since MATLAB arrays are referenced to 
    // 1, not zero
    idxptr[i] = (double) (idx + 1);

    // Then, the actual point -- if requested
    if (nlhs > 1) {
      if (isDouble) {
          for (int j = 0; j < ndim; j++)
            dpntptr[j * npoints + i] = tree->points[idx*tree->ndim+j];
      } 
      else {
        for (int j = 0; j < ndim; j++)
          fpntptr[j * npoints + i] = tree->points[idx*tree->ndim+j];
      }
    }
  }
  // Deallocate the point to check
  delete curPoint;
  
  // Get rid of the tree
  // (This doesn't delete the serialized data)
  delete tree;
  return;
  
} // end of kdtree_closestpoint
