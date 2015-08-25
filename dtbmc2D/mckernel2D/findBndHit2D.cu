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
// CUDA PTX kernel to find the first index of hitting upper bound with pos in DTB 
// 2D Monte Carlo simulation.
//

// Copyright 2015 Jian Wang

__global__ void withPos(const double* randpool, 
		const double* bup, const double* blo, double* iup, double* ilo,                                                        
		const int ntime, const int ntrials, 
		const double y0, const double mu, const double cov, double* pos)
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
			
			// Check positive non-absorptive possibility.
			if (rup > rlo) {
				pos[k+1] = 1.0;
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



