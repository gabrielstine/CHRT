//
// GPU PTX kernel to find the index of 1st up or lower bound hitting
//

// Copyright 2014 Jian Wang


__global__ void find1HitBnd(const bool* upBnd, double* iup,
                            const bool* loBnd, double* ilo,
                            const int ntime, const int ntrials)
{
    int thisThread = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (thisThread < ntrials) {
        int m = thisThread * ntime;
        for (int n = 0; n < ntime; ++n) {
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