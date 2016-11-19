//Getting started with OpenMP
#include <iostream>
#include <omp.h>
int main()
{
#pragma omp parallel
	{
	int NumThreads,threadID;
    NumThreads = omp_get_num_threads();
	threadID = omp_get_thread_num();
#pragma omp critical
{
    std::cout << "This is thread " << threadID << std::endl;
	}
    }

return 0;
}

