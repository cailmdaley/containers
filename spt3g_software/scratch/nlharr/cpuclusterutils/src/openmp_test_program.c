#include <omp.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>

int main(){
#pragma omp parallel
	{
		time_t timer;
		struct tm* tm_info;
		time(&timer);
		printf("Start %s", asctime(localtime(&timer)));
		sleep(3);
		time(&timer);
		printf("End %s", asctime(localtime(&timer)));

		
	}	
}
