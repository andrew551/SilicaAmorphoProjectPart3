      Program convert_to_1d
      
      implicit none
       
      real*8 :: var
      integer :: ind(6),reclen,reclen1,cnt
	  character*100 fname
	  
	  fname=''
	  call getarg(1,fname)
	  if(fname=='') then
	    write(*,*)'Enter d3mat filename to convert'
		stop
	  endif
	  
      cnt=1
      open(unit=26,file='d3mataaa.dat',status='unknown')
      
      inquire(iolength=reclen)ind(:)
      inquire(iolength=reclen1)var
      open(unit=251,file='d3mataac_ind.bin',status='unknown',access='direct',form='unformatted',recl=reclen)
      open(unit=252,file='d3mataac_val.bin',status='unknown',access='direct',form='unformatted',recl=reclen1)
   
      open(113,file=trim(fname),status='old')
      do while(.true.)
        read(113,*,end=123)ind(:),var
		
		if(ind(1)==ind(3) .and. ind(3)==ind(5)) write(26,'(6I6,E20.10)')ind(1:6),var
		
		if(dabs(var)< 1d-10) cycle
		
		write(251,rec=cnt)ind(:)
        write(252,rec=cnt)var
		cnt=cnt+1
    
      enddo
123	  close(113)

      close(251)
      close(252)
      close(26)
      open(111,file='nlines',status='unknown')	  
      write(111,*)(cnt-1)
	  close(111)
 
      End program convert_to_1d