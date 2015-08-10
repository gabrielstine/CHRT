//
// GPU PTX kernel to find the index of 1st up or lower bound hitting
// 
//

// Copyright 2014 Jian Wang


__global__ void find1HitBndPos(const bool* upBnd, double* iup,
                            const bool* loBnd, double* ilo,
                            const int ntime, const int ntrials,
                            const double* dft, double* ipos)
{
    int thisThread = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (thisThread < ntrials) {
        int m = thisThread * ntime;
        for (int n = 0; n < ntime; ++n) {
        	if (dft[m] > 0.0) {        		
        	    ipos[m] = 1.0;        	    
        	}
        	
            if (upBnd[m]) {
                iup[m] = 1.0;
                break;
            }
            else if (loBnd[m]) {
                ilo[m] = 1.0;
                break;
            }
                        
            ++m;
        }
    }
}