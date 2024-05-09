      Program Grun
       
      implicit none
      include 'mpif.h' 
      
      real*8, parameter :: pi=4.0d0*datan(1.0d0), eJ=1.602176487d-19, a_to_m=1.0d-10, &
                     amu=1.660599d-27, freq_to_wthz=dsqrt(eJ/amu/a_to_m**2)/1.0d12, &
                     togpa=eJ/a_to_m**3/1.0d9, &
                     hbar = 1.054571d-34,& ! J*s
                     kb=1.380649d-23, &
                     thtocm=1.0d0/0.0299792458d0, &
                     mm(1:2)=[15.9994d0, 28.0855d0]
      
      integer ::mj,i,j,k,ii,jj,kk,dd,nu,natom,nlines, typ, ind(6),itp,reclen,reclen1,ga,gisto(-50:50,4), iii,cnt,length
      real*8 :: dum,box(3),x1,x2,vol,qq,ssum,dx(3),ssum1,v(6),var,T, arg,eng,smin,smax,smin1,smax1,rr(3)
      
      real*8, allocatable :: mass(:),eigen(:),phir(:),rrr(:),& 
                             gamma1(:),allgamma(:),nj(:),CV(:),ww(:),tempa(:,:)
      integer, allocatable :: ia(:),ib(:)
 
      integer :: m_rank, m_size, m_ierror
      character*2 ch
      character*200 path_model,path_w,path_eigen,path_fc3,isrmean
 
      m_rank=0
      m_size=1
      
      call MPI_INIT(m_ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, m_size, m_ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, m_rank, m_ierror)
  
      call getarg(1,isrmean)
    
      if(trim(isrmean)/='MAIN' .and. trim(isrmean)/='IS') then 
        write(*,*)'unspecified type:', trim(isrmean)
        stop
      endif
	  
	  open(11,file='inputs',status='old')
      read(11,*)path_model
      if(m_rank==0)write(*,*)path_model
      read(11,*)path_w
      if(m_rank==0)write(*,*)path_w
      read(11,*)path_eigen
      if(m_rank==0)write(*,*)path_eigen
      read(11,*)path_fc3
      if(m_rank==0)write(*,*)path_fc3
  
      close(11)
      
      open(unit=45,file=trim(path_model),status='old')
 
      read(45,*)natom
    
      allocate(mass(3*natom),eigen(3*natom))
      allocate(gamma1(3*natom),allgamma(3*natom), rrr(3*natom),ww(3*natom),nj(natom*3),CV(3*natom))
 
      read(45,*)   !box(:)
      !vol=box(1)*box(2)*box(3)
 
      do i=1,natom
        typ=0
        read(45,*)ch,rrr(3*i-2:3*i)
        if(trim(ch)=='O')typ=1
        if(trim(ch)=='Si')typ=2
        mass(3*i-2:3*i) = dsqrt(mm(typ))
      enddo
      close(45)
  
      if(trim(isrmean)=='IS') then 
        open(333,file='rmean.dat',status='old')
        do i=1,natom
          read(333,*)j,rrr(i*3-2:i*3)
        enddo
        close(333)
        write(*,*)'IS case'
      endif
 
      open(55,file=trim(path_w),status='old')
      do mj=1,3*natom
        read(55,*)i,ww(mj)
      enddo
      close(55)

      inquire(iolength=reclen) eigen(:)
      open(1000+m_rank,file=trim(path_eigen),form='unformatted',access='direct',status='old',recl=reclen) 
 
      if(m_rank==0) then
	    allocate(tempa(3*natom,3*natom))    
        tempa=0.0d0      

        inquire(iolength=reclen)ind(1:6)
        inquire(iolength=reclen1)var
        open(100,file=trim(path_fc3)//'fc3_ind.bin',form='unformatted',action='read',access='direct',recl=reclen,status='old')
        open(200,file=trim(path_fc3)//'fc3_val.bin',form='unformatted',action='read',access='direct',recl=reclen1,status='old')
         
        cnt=1
        do while(.true.)
          read(100,rec=cnt,err=12)ind(:)
          read(200,rec=cnt,err=12)var
	      ii=ind(1)*3-3+ind(2)
          jj=ind(3)*3-3+ind(4)
          kk=ind(5)*3-3+ind(6)
          tempa(ii,jj)=tempa(ii,jj)+var*rrr(kk)
          cnt=cnt+1
        enddo
12      close(100)
        close(200)
        cnt=0
        do ii=1,3*natom
          do jj=1,3*natom
            if(dabs(tempa(ii,jj))>1d-6) cnt=cnt+1
          enddo
        enddo
		
	    nlines=cnt
         
        allocate(phir(nlines),ia(nlines),ib(nlines))
        cnt=0
        do ii=1,3*natom
          do jj=1,3*natom
            if(dabs(tempa(ii,jj))>1d-6) then
              cnt=cnt+1
              phir(cnt)=tempa(ii,jj)
              ia(cnt)=ii
              ib(cnt)=jj
            endif
          enddo
        enddo 
	    deallocate(tempa)
      endif
	  
      call MPI_BCAST(nlines,1,MPI_INTEGER,0,MPI_COMM_WORLD,m_ierror)
	  
      if(m_rank /= 0) allocate(phir(nlines),ia(nlines),ib(nlines))
 
      call MPI_BCAST(phir,nlines,MPI_DOUBLE,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(ia,nlines,MPI_INTEGER,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(ib,nlines,MPI_INTEGER,0,MPI_COMM_WORLD,m_ierror)
   
      do mj=1+3*natom/m_size*m_rank,3*natom/m_size*(m_rank+1)
        if(mj<4)cycle
        write(*,*)'mj=',mj
        read(1000+m_rank,rec=mj)eigen(:)
 
        eigen(:)=eigen(:)/mass(:)
      
        var=0.0
        do cnt=1,nlines
          var=var+phir(cnt)*eigen(ia(cnt))*eigen(ib(cnt))
        enddo
        gamma1(mj)=-var/6.0d0/ww(mj)
      enddo
  
      call MPI_ALLREDUCE(gamma1, allgamma, 3*natom, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, m_ierror)

      if(m_rank==0) then
        if(trim(isrmean)=='IS') then 
          open(10,file='gruneisen_IS.dat',status='unknown')
        else  
          open(10,file='gruneisen_MAIN.dat',status='unknown')
        endif
        do mj=1+3,3*natom
          write(10,*)freq_to_wthz*sqrt(dabs(ww(mj)))/2/pi*thtocm,allgamma(mj)
        enddo
        close(10)
         
        
      endif
	  
	  deallocate(phir,ia,ib)
 
      call MPI_FINALIZE(m_ierror)
 
 
      End program Grun
   
