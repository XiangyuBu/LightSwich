program MC
USE global_parameters
USE utility_routines    

implicit none
INCLUDE 'mpif.h'
call MPI_INIT(ierr)
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

call initialize()

call var_s()

call MPI_FINALIZE(ierr)
print*, "end of program SCMFT"

End program MC







 



    









     















