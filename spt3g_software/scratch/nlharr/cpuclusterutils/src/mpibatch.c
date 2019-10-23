#define _GNU_SOURCE
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>

#include <ctype.h>

int main( int argc, char **argv )
{
  int rank;
  int size;
  int retval=0;
  
  if(argc < 2){
    printf("Usage mpibatch <batchfile>\n");
    exit(1);
  }
  
  MPI_Init( 0, 0 );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("rank %d size %d\n",rank,size);  
  
  FILE *input_file=fopen(argv[1], "r");
  if(!input_file){
    printf("Failed to open file \"%s\"\n", argv[1]);
    exit(1);
  }
  
  int currentLine=0;
  
  while(1){
    //printf(".\n");
    //char *command=readline(input_file, &retval);
    char *command=0;
    size_t buffsize=0;
    size_t rv=getline(&command, &buffsize, input_file);  
    command[strlen(command)-1]='\0';
    if(rv == -1){ 
       break;
    }
    
    if((currentLine % size) == rank){
      //printf("%d %d: Executing command: %s\n", rv, rank, command);
      int ret=system(command);
      
      //printf("PE: %d, Command: %s\n", rank, 
      //     command);
    }
    free(command);
    command=0;
    buffsize=0;
    currentLine++;
  }
  
  
  
  MPI_Finalize();
  printf("Exit code %d\n", retval);

  fclose(input_file);

  return retval;
  
}
