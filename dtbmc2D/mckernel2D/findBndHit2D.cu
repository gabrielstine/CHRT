//
// CUDA PTX kernel to find the first index of hitting upper bound in DTB 
// 2D Monte Carlo simulation.
//

// Copyright 2015 Jian Wang

__global__ void withoutPos(const double* randpool, 
		const double* bup, const double* blo, double* iup, double* ilo,                                                        
		const int ntime, const int ntrials, 
		const double y0, const double mu, const double cov)
{
	int thisThread = blockIdx.x * blockDim.x + threadIdx.x;

	if (thisThread < ntrials) {                
		double rup = y0, rlo = y0; // Upper & lower competing races. 
		double thisBUp, thisBLow; // Upper & lower bound height at this time.

		int m = thisThread * ntime * 3.0; // Index of rand number.
		int k = thisThread * ntime; // Index of boundary.

		for (int n = 0; n < ntime-1; ++n) {
			rup += ((1-abs(cov))*randpool[m] + abs(cov)*randpool[m+1]) + mu;
			rlo += ((1-abs(cov))*randpool[m+2] + cov*randpool[m+1])    - mu;

			thisBUp = bup[n+1];
			thisBLow = blo[n+1];

			// Check whether hitting lower boundary.
			if (rup < thisBLow) { 
				rup = thisBLow;
			}

			if (rlo < thisBLow) {
				rlo = thisBLow;
			}

			// Check whether hitting upper boundary.
			if (rup > thisBUp) {
				iup[k+1] = 1.0;
				break;
			}
			else if (rlo > thisBUp) {
				ilo[k+1] = 1.0;
				break;
			}

			m += 3;
			k++;
		}
	}
}


//
// GPU PTX kernel to find the index of 1st up or lower bound hitting
// and the possibility of in positive half axis.
// 

// Copyright 2014 Jian Wang

//__global__ void withPos(const bool* upBnd, double* iup,		
//                        const bool* loBnd, double* ilo,
//                        const int ntime, const int ntrials,
//                        const double* dft, double* ipos)
//{
//    int thisThread = blockIdx.x * blockDim.x + threadIdx.x;
//    
//    if (thisThread < ntrials) {
//        int m = thisThread * ntime;
//        for (int n = 0; n < ntime; ++n) {
//        	if (dft[m] > 0.0) {        		
//        	    ipos[m] = 1.0;        	    
//        	}
//        	
//            if (upBnd[m]) {
//                iup[m] = 1.0;
//                break;
//            }
//            else if (loBnd[m]) {
//                ilo[m] = 1.0;
//                break;
//            }
//                        
//            ++m;
//        }
//    }
//}

