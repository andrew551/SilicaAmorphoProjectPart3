!   it seems that na must be divisble by nblk**2
! it seems that it runs only on number of processors = nblk**2
! number of files generated and read has to match the number of processors

program gather
implicit none
integer, parameter :: blocks=4, bsize=3664*3/blocks

integer :: i,j,k,k1,reclen,allreclen,myid,ul,ur,kpos,cnt
real*8 :: bvec(bsize),allvec(bsize*blocks)

character*8 :: ch

inquire(iolength=reclen)bvec(:)
inquire(iolength=allreclen)allvec(:)

open(9222,file='eigen_all.bin',form='unformatted', &
                    access='direct', status='unknown',recl=allreclen)
 
    do j=1,blocks
      do i=1,blocks
      
        myid=(i-1)+(j-1)*blocks
		
		write(*,*)myid
    
        if (myid<10000)write(ch,'(i4)')myid
        if (myid<1000)write(ch,'(i3)')myid
        if (myid<100)write(ch,'(i2)')myid
        if (myid<10)write(ch,'(i1)')myid
    
        open(111+myid,file='eigen_'//trim(ch)//'.bin',form='unformatted', &
                        access='direct', status='old',recl=reclen)
                                  
      enddo
	enddo
	
	cnt=1
	do k1=1,bsize/blocks
	  do j=1,blocks
	    do k=1,blocks
	      kpos=k+(k1-1)*blocks
		  !write(12345,*)cnt,k1,kpos
          do i=1,blocks
    	    myid=(i-1)+(j-1)*blocks
		  
    	     read(111+myid,rec=kpos)bvec(:)
    	    ul=(i-1)*bsize+1
    	    ur=i*bsize
    	   allvec(ul:ur)=bvec(:)
    	  enddo
    	
    	 write(9222,rec=cnt)allvec(:)
		cnt=cnt+1
		enddo
      enddo
    enddo
    close(9222)
 
end program gather
