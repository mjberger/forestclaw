/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#ifdef FCLAW_ENABLE_MPI

#include <mpi.h>
#include <unistd.h>    //  To get process ids


void fclaw2d_mpi_debug()
{
  int i;

  /* Find out process rank */
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Find out number of processes */
  int num;
  MPI_Comm_size(MPI_COMM_WORLD, &num);

  /* Don't do anything until we are done with this part */
  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0)
  {
      printf("\n");
      printf("Getting setup for parallel debugging\n");
      printf("------------------------------------\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(i = 0; i < num; i++)
  {
      if (my_rank == i)
      {
          printf("Proc %d with process %d is waiting to be attached\n",my_rank,getpid());
          fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
  }

  int ii = 0;
  while (ii == 0)  /* (gdb) set ii=1 */
  {
    /* Code will stop here;  set ii=1 to continue in gdb */
  }
}

#else
void fclaw2d_mpi_debug()
{
}
#endif
