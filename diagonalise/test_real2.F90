!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
!>
!> Fortran test programm to demonstrates the use of
!> ELPA 2 real case library.
!> If "HAVE_REDIRECT" was defined at build time
!> the stdout and stderr output of each MPI task
!> can be redirected to files if the environment
!> variable "REDIRECT_ELPA_TEST_OUTPUT" is set
!> to "true".
!>
!> By calling executable [arg1] [arg2] [arg3] [arg4]
!> one can define the size (arg1), the number of
!> Eigenvectors to compute (arg2), and the blocking (arg3).
!> If these values are not set default values (4000, 1500, 16)
!> are choosen.
!> If these values are set the 4th argument can be
!> "output", which specifies that the EV's are written to
!> an ascii file.
!>
program test_real_example

!-------------------------------------------------------------------------------
! Standard eigenvalue problem - REAL version
!
! This program demonstrates the use of the ELPA module
! together with standard scalapack routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
!-------------------------------------------------------------------------------
   
   use iso_c_binding

   use elpa

#ifdef HAVE_MPI_MODULE
   use mpi
   implicit none
#else
   implicit none
   include 'mpif.h'
#endif

   
   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! na:   System size
   ! nev:  Number of eigenvectors to be calculated
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer           :: nblk,reclen
   integer                          :: na, nev

   integer                          :: np_rows, np_cols, na_rows, na_cols

   integer                          :: myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer                          :: i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol

   integer, external                :: numroc

   real(kind=c_double), allocatable :: a(:,:), z(:,:), ev(:)
   
   real(kind=c_double) :: aaa
   

   !integer                          :: iseed(4096) ! Random seed, size should be sufficient for every generator

   integer                          :: STATUS
   integer                          :: success
   character*10                     :: ch
   integer                          :: j,il,iu,jl,ju

   integer, parameter               :: error_units = 0

   class(elpa_t), pointer           :: e
   !-------------------------------------------------------------------------------


   ! hard-coded parameters 
   na = 5001*3     ! number of modes = 3 * number of atoms
   nev = 5001*3    ! number of eigenvectors
   nblk = 3       ! number of blocks, has to be carefully choosen, 
   ! it seems that na must be divisble by nblk**2
   ! it seems that it runs only on number of processors = nblk**2
   ! we need to pad the dynamical matrix in order for all block sizes to have the same dimensions
   ! the block sizes are written in the output, see np_rows,np_cols (line 144)
   ! probably, for a 648 atom model, nblk=4, nblk=8 will work, but nblk=10 will not work
   

   call mpi_init(mpierr)
   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
   call mpi_comm_size(mpi_comm_world,nprocs,mpierr)

   do np_cols = NINT(SQRT(REAL(nprocs))),2,-1
     if(mod(nprocs,np_cols) == 0 ) exit
   enddo
   ! at the end of the above loop, nprocs is always divisible by np_cols

   np_rows = nprocs/np_cols
   
   if(myid==0) then
     write(*,*)'np_rows,np_cols:', np_rows,np_cols
     write(*,*)'nblock:', nblk
   endif

   ! initialise BLACS

   my_blacs_ctxt = mpi_comm_world
   call BLACS_Gridinit(my_blacs_ctxt, 'C', np_rows, np_cols)
   if (myid==0) then
     print '(a)','| Past BLACS_Gridinit.'
   end if
   
   call BLACS_Gridinfo(my_blacs_ctxt, nprow, npcol, my_prow, my_pcol)

   if (myid==0) then
     print '(a)','| Past BLACS_Gridinfo.'
   end if
   ! determine the neccessary size of the distributed matrices,
   ! we use the scalapack tools routine NUMROC

   na_rows = numroc(na, nblk, my_prow, 0, np_rows)
   na_cols = numroc(na, nblk, my_pcol, 0, np_cols)
   
   
     write(*,'(A4, 5i6)')'my: ', myid,na_rows,na_cols,my_prow,my_pcol
   

   ! set up the scalapack descriptor for the checks below
   ! For ELPA the following restrictions hold:
   ! - block sizes in both directions must be identical (args 4 a. 5)
   ! - first row and column of the distributed matrix must be on
   !   row/col 0/0 (arg 6 and 7)

   if (np_rows .ne. na_cols) then
      print *, "na_rows must be equal to na_cols for elpa to work"
      stop
    endif

   call descinit(sc_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info)

   if (info .ne. 0) then
     write(error_units,*) 'Error in BLACS descinit! info=',info
     write(error_units,*) 'Most likely this happend since you want to use'
     write(error_units,*) 'more MPI tasks than are possible for your'
     write(error_units,*) 'problem size (matrix size and blocksize)!'
     write(error_units,*) 'The blacsgrid can not be set up properly'
     write(error_units,*) 'Try reducing the number of MPI tasks...'
     call MPI_ABORT(mpi_comm_world, 1, mpierr)
   endif

   if (myid==0) then
     print '(a)','| Past scalapack descriptor setup.'
   end if

   allocate(a (na_rows,na_cols))
   allocate(z (na_rows,na_cols))
   
   if (myid==0) then
     print '(a)','| Mem allocated.'
   end if
   
   	 inquire(iolength=reclen)z (:,1)
   
    il=my_prow*na_rows+1
	iu=(my_prow+1)*na_rows
	jl=my_pcol*na_cols+1
	ju=(my_pcol+1)*na_cols
	
	 !write(*,'(A4,5i5)')'2 ',myid,il,iu,jl,ju
   
   a=0.0d0
   z=0.0d0

   allocate(ev(na))
   
   ev=0.0d0

   print '(a)','| Matrix reading...name must be dmat.dat'
   
   open(111,file='dmat.dat',status='old')
 123  read(111,*,end=666)i,j,aaa
      if(i>=il .and. i<=iu .and. j>=jl .and. j<=ju) then
        a(i-il+1,j-jl+1) = aaa
	  end if
     goto 123
666   close(111)
 

   
    print '(a)','| Matrix has been set up.'
 

   !-------------------------------------------------------------------------------

   if (elpa_init(20171201) /= elpa_ok) then
     print *, "ELPA API version not supported"
     stop
   endif
   e => elpa_allocate()

   ! set parameters decribing the matrix and it's MPI distribution
   call e%set("na", na, success)
   call e%set("nev", nev, success)
   call e%set("local_nrows", na_rows, success)
   call e%set("local_ncols", na_cols, success)
   call e%set("nblk", nblk, success)
   call e%set("mpi_comm_parent", mpi_comm_world, success)
   call e%set("process_row", my_prow, success)
   call e%set("process_col", my_pcol, success)

   success = e%setup()
   
   call e%set("omp_threads", 2, success)

   call e%set("solver", elpa_solver_2stage, success)
   
   ! set the AVX BLOCK2 kernel, otherwise ELPA_2STAGE_REAL_DEFAULT will
 ! be used
    !call e%set("real_kernel", ELPA_2STAGE_REAL_AVX_BLOCK2, success)


   ! Calculate eigenvalues/eigenvectors

   if (myid==0) then
     print '(a)','| Entering two-step ELPA solver ... '
     print *
   end if

   call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
   call e%eigenvectors(a, ev, z, success)
  
    if(myid<10000) write(ch,'(i4)')myid
	if(myid<1000) write(ch,'(i3)')myid
	if(myid<100) write(ch,'(i2)')myid
	if(myid<10) write(ch,'(i1)')myid
 
    open(111,file='eigen_'//trim(ch)//'.bin',form='unformatted',access='direct',status='unknown',recl=reclen)
 
     do j=1,na_cols
       write(111,rec=j)z(:,j)
     enddo
 
    close(111)
   
   if (myid==0) then
     open(222,file='frequencies.dat',status='unknown')
     do i=1,na
	   write(222,*)i,ev(i)
	 enddo
	 close(222)
   end if

   if (myid==0) then
     print '(a)','| Two-step ELPA solver complete.'
     print *
   end if

   call elpa_deallocate(e)
   call elpa_uninit()

   call blacs_gridexit(my_blacs_ctxt)
   call mpi_finalize(mpierr)

end

