/*********************************************************************
* SCCfull.c 
* Created by: M. Heinz 27Sept2004
* Modified from InnerSACmex.c
* 
* This is a MEX-file for MATLAB.
*
* Computes Shuffled Cross Correlogram from two Spike Matrices
* 
* Call:  [SCC,delay] = SCCfull(SpikeMAT1,NUMspikes1,TOTALspikes1,SpikeMAT2,NUMspikes2,TOTALspikes2);
* 		SpikesMAT1(2): (Kmax x NUMreps) matrix with spike trains (NaNs fill end of each column for shorter spike trains)
*		NUMspikes1(2): (1x NUMreps) vector with number of spikes in each train
*		TOTALspikes1(2): needed to allocate SCC
*
*		SCC: Shuffled Cross Correlogram.
*		delay: delay corresponding to SCC. 
*		
*********************************************************************/

#include "mex.h"

void SCCfull(double *SpikeMAT1, double *NUMspikes1, long NUMspikeREPS1, long Kmax1, double *SpikeMAT2, double *NUMspikes2, long NUMspikeREPS2, long Kmax2, double *DelayBinWidth, long nzmax, double *ints)
{
   long REPindREF,SpikeIND,REPindCOMP,COMPSpikeIND,intIND=0;
   
   for (REPindREF=0; REPindREF<NUMspikeREPS1; REPindREF++) {
	   for (SpikeIND=0; SpikeIND<NUMspikes1[REPindREF]; SpikeIND++) {
		   for (REPindCOMP=0; REPindCOMP<NUMspikeREPS2; REPindCOMP++) {
 				for (COMPSpikeIND=0; COMPSpikeIND<NUMspikes2[REPindCOMP]; COMPSpikeIND++) {
					intIND=nzmax+floor((SpikeMAT2[REPindCOMP*Kmax2+COMPSpikeIND]-SpikeMAT1[REPindREF*Kmax1+SpikeIND])/DelayBinWidth[0]);
					
					if (intIND>=0 && intIND<=2*nzmax)
					{
					   ints[intIND]+=1;
					}
      }  }	}  }
}


// the gateway function
// [SCC,delay] = SCCfull(SpikeMAT1,NUMspikes1,TOTALspikes1,SpikeMAT2,NUMspikes2,TOTALspikes2);
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   double   *SpikeMAT1,*SCC,*NUMspikes1, *DelayMax, *DelayBinWidth; // input-output arrays
   double   *SpikeMAT2,*NUMspikes2; // input-output arrays
	long     TOTALspikes1,TOTALspikes2; // input scalar
   double   *delay;  // output scalar
   long     NUMspikeREPS1,Kmax1;
   long     NUMspikeREPS2,Kmax2, nzmax;
   

	int      i,j;

//	__asm int 3h ; // Setting a breakpoint for debugging (compile with the -g flag!)
							// For "Nel_matlab", run 'doallmex' with '(0,1)' flags.
							// Run Nel software as usual, and Matlab will automatically
							//  break at this line and open this file in Visual Studio
							//  debugger environment.

   /*  check for proper number of arguments */
   if(nrhs!=8)
    mexErrMsgTxt("Eight inputs required.");
   if(nlhs!=2)     mexErrMsgTxt("Two outputs required.");
  
  	/*  create a pointer to the input matrix SpikeMAT1(2) */
   // SpikeMAT1 is stored in MATLAB as: each row is one spike train
	//   in MEX/C, linear indexing goes down columns first, so we pass the transpose to 
	//   allow easier indexing
   SpikeMAT1 = mxGetPr(prhs[0]);  // each column is one spike train
	/*  get the dimensions of the matrix input SpikeMAT1 */
   NUMspikeREPS1 = mxGetN(prhs[0]);
   Kmax1 = mxGetM(prhs[0]);
   SpikeMAT2 = mxGetPr(prhs[3]);
   NUMspikeREPS2 = mxGetN(prhs[3]);
   Kmax2 = mxGetM(prhs[3]);
  
	/*  create a pointer to the input vector NUMspikes1(2) */
   NUMspikes1 = mxGetPr(prhs[1]);
   NUMspikes2 = mxGetPr(prhs[4]);

   /*  get the scalar input TOTALspikes1(2) */
   TOTALspikes1 = mxGetScalar(prhs[2]);
   TOTALspikes2 = mxGetScalar(prhs[5]);

   /* check to make sure the TOTALspikes1(2) input arguments are scalars */
   if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
      mexErrMsgTxt("Input TOTALspikes1 must be a scalar.");
   }
   if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
      mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
      mexErrMsgTxt("Input TOTALspikes2 must be a scalar.");
   }

      /*  get the scalar input DelayMax */
   DelayMax = mxGetPr(prhs[6]);
   /* check to make sure the last input argument is a scalar */
   if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
      mxGetN(prhs[6])*mxGetM(prhs[6])!=1 ) {
      mexErrMsgTxt("Input DelayMax must be a scalar.");
   }

    /*  get the scalar input DelayBinWidth  */
   DelayBinWidth = mxGetPr(prhs[7]);
   /* check to make sure the last input argument is a scalar */
   if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
      mxGetN(prhs[7])*mxGetM(prhs[7])!=1 ) {
      mexErrMsgTxt("Input DelayBinWidth must be a scalar.");
   }
   
   /*  set the output pointer to the output matrix SCC */	
	nzmax = floor(DelayMax[0]/DelayBinWidth[0]);
	plhs[0] = mxCreateDoubleMatrix(1,2*nzmax+1, mxREAL);
   /*  create a C pointer to a copy of the output matrix SCC */
   SCC = mxGetPr(plhs[0]);
  
   /*  set the output pointer to the output scalar delay */
   plhs[1] = mxCreateDoubleMatrix(1,2*nzmax+1, mxREAL);
   /*  create a C pointer to a copy of the output variable delay */
   delay = mxGetPr(plhs[1]);
   
   for (i=0; i<=2*nzmax; i++)
   {
	   delay[i]=-DelayMax[0]+DelayBinWidth[0]*i;
   }
   
   /*  call the C subroutine */
   SCCfull(SpikeMAT1,NUMspikes1,NUMspikeREPS1,Kmax1,SpikeMAT2,NUMspikes2,NUMspikeREPS2,Kmax2,DelayBinWidth, nzmax,SCC);
	
	}
