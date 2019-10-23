#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main( int argc, char **argv )
{
  int rank;
  int size;
  int retval=0;
  
  if(argc < 2){
    printf("Usage mpibatch <command_0> <command_1> <command_2> ...\n");
    exit(1);
  }
  
  MPI_Init( 0, 0 );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("rank %d size %d\n",rank, size);  

  int currentLine=0;
  while(1){
    if (currentLine + 1 >= argc) break;
    if((currentLine % size) == rank){
      int ret = system(argv[currentLine+1]);
    }
    currentLine++;
  }
  MPI_Finalize();
  printf("Exit code %d\n", retval);
  return retval;
}
