      Program Grun
       
      implicit none
      include 'mpif.h' 
  
       real*8, parameter :: pi=4.0d0*datan(1.0d0), eJ=1.602176487d-19, a_to_m=1.0d-10, &
                     amu=1.660599d-27, freq_to_wthz=dsqrt(eJ/amu/a_to_m**2)/1.0d12, &
                     togpa=eJ/a_to_m**3/1.0d9, &
                     hbar = 1.054571d-34,& ! J*s
                     kb=1.380649d-23, &
                     thtocm=1.0d0/0.0299792458d0
      
      integer :: i,j,k,natom,mj,cnt,ia,al,be,ic,ga,length0,length1,reclen,typ,ii,jj,kk,length_dm,ind(6),reclen1
      integer :: m_rank, m_size, m_ierror,itp
      real*8 :: var,qq,r(3),mm(10),eng,T,arg,box(3),x1,x2,vol,dx(3),tmp
      real*8, allocatable :: dphir(:),phir(:,:)
      integer, allocatable :: v1_1(:),v2_1(:),v1_2(:),v2_2(:),&
               v1_3(:),v2_3(:),vd1(:),vd2(:)
      
      real*8, allocatable :: eigen(:),mass(:),rr(:),ww(:),gamma1(:),allgamma(:),nj(:),CV(:)
 
      character*100 path_model,path_w,path_eigen,path_fc3,isrmean
 
      m_rank=0
      m_size=1
      
      call MPI_INIT(m_ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, m_size, m_ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, m_rank, m_ierror)
       
      open(11,file='inputs',status='old')
      read(11,*)isrmean
      write(*,*)isrmean
      if(trim(isrmean)/='MAIN' .and. trim(isrmean)/='IS') then 
        write(*,*)'unspecified type:', trim(isrmean)
        stop
      endif
      read(11,*)path_model
      write(*,*)path_model
      read(11,*)path_w
      write(*,*)path_w
      read(11,*)path_eigen
      write(*,*)path_eigen
      read(11,*)path_fc3
      write(*,*)path_fc3
  
      close(11)
      
      open(unit=45,file=trim(path_model),status='old')
      read(45,*)
      read(45,*)
      read(45,*)natom
      read(45,*)itp
      read(45,*)
      
      allocate(mass(3*natom),rr(natom*3),&
           gamma1(3*natom),allgamma(3*natom),ww(3*natom),eigen(natom*3),nj(natom*3),CV(3*natom))
 
      read(45,*)x1,x2
      box(1)=x2-x1
      read(45,*)x1,x2
      box(2)=x2-x1
      read(45,*)x1,x2
      box(3)=x2-x1
      vol=box(1)*box(2)*box(3)   
      
      read(45,*)
      read(45,*)
      read(45,*)
      
      do i=1,itp
        read(45,*)j,mm(i)
      enddo
       
      do i=1,8
        read(45,*)
      enddo
      do i=1,natom
        read(45,*)j,typ,qq,r(:)
        rr(j*3-2:j*3)=r(:)
        mass(3*j-2:3*j) = dsqrt(mm(typ))
      enddo
      close(45)
      if(trim(isrmean)=='IS') then 
      open(333,file='rmean.dat',status='old')
        do i=1,natom
          read(333,*)j,rr(i*3-2:i*3)
        enddo
        close(333)
        write(*,*)'IS case'
      endif
  
      open(111,file='nlines',status='old')
  
      read(111,*)length0,length1      ! wc in d3mataac.bin
      close(111)
  
      open(55,file=trim(path_w),status='old')
      do mj=1,3*natom
        read(55,*)i,ww(mj)
      enddo
      close(55)
      
      inquire(iolength=reclen) eigen(:)
      open(551,file=trim(path_eigen),form='unformatted',access='direct',status='old',recl=reclen)
   
      allocate(dphir(length0),vd1(length0),vd2(length0))
      allocate(v1_1(length1),v2_1(length1),phir(length1,3))
      allocate(v1_2(length1),v2_2(length1))
      allocate(v1_3(length1),v2_3(length1))
 
      if(m_rank==0) then

        inquire(iolength=reclen)ind(:)
        inquire(iolength=reclen1)var 
        
        open(22111,file=trim(path_fc3)//'d3mataaa_ind.bin',action='read',form='unformatted',access='direct',status='old',recl=reclen)
        open(55666,file=trim(path_fc3)//'d3mataaa_val.bin',action='read',form='unformatted',access='direct',status='old',recl=reclen1)
      
        cnt=1
        do while(.true.)
          read(22111,rec=cnt)ind(:)
          read(55666,rec=cnt)var
          dphir(cnt)=var*rr((ind(1)-1)*3+ind(6))
          vd1(cnt)=(ind(1)-1)*3+ind(2)
          vd2(cnt)=(ind(1)-1)*3+ind(4)
         
          cnt=cnt+1
          if(cnt>length0) goto 123
        enddo
123     close(22111)    
        close(55666)
      
        open(22111,file=trim(path_fc3)//'d3mataac_ind.bin',action='read',form='unformatted',access='direct',status='old',recl=reclen)
        open(55666,file=trim(path_fc3)//'d3mataac_val.bin',action='read',form='unformatted',access='direct',status='old',recl=reclen1)
        cnt=1
        do while(.true.)
          read(22111,rec=cnt)ind(:)
          read(55666,rec=cnt)var 
   
          phir(cnt,1)=var*rr((ind(5)-1)*3+ind(6))  ! V(a,a,b)*R(b)
          phir(cnt,2)=var*rr((ind(3)-1)*3+ind(4))  ! V(a,b,a)*R(a)
          phir(cnt,3)=var*rr((ind(3)-1)*3+ind(4))  ! V(b,a,a)*R(a)
 
          v1_1(cnt)=(ind(1)-1)*3+ind(2)  ! a\alpha
          v2_1(cnt)=(ind(3)-1)*3+ind(4)  ! a\alpha
          
          v1_2(cnt)=(ind(1)-1)*3+ind(2)  ! a\alpha
          v2_2(cnt)=(ind(5)-1)*3+ind(6)  ! b\beta
           
          v1_3(cnt)=(ind(5)-1)*3+ind(6)  ! b\beta
          v2_3(cnt)=(ind(1)-1)*3+ind(2)  ! a\alpha

          cnt=cnt+1
          if(cnt>length1) goto 1234
        enddo
1234    close(22111)    
        close(55666)
        write(*,*)'Ready to send'
      endif
      
      call MPI_BCAST(dphir,length0,MPI_DOUBLE,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(vd1,length0,MPI_INT,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(vd2,length0,MPI_INT,0,MPI_COMM_WORLD,m_ierror)
      
      call MPI_BCAST(phir(:,1),length1,MPI_DOUBLE,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(phir(:,2),length1,MPI_DOUBLE,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(phir(:,3),length1,MPI_DOUBLE,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(v1_1,length1,MPI_INT,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(v2_1,length1,MPI_INT,0,MPI_COMM_WORLD,m_ierror)
 
      call MPI_BCAST(v1_2,length1,MPI_INT,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(v2_2,length1,MPI_INT,0,MPI_COMM_WORLD,m_ierror)
 
      call MPI_BCAST(v1_3,length1,MPI_INT,0,MPI_COMM_WORLD,m_ierror)
      call MPI_BCAST(v2_3,length1,MPI_INT,0,MPI_COMM_WORLD,m_ierror)
      
      
      gamma1=0.0d0
      
      do mj=1+3*natom/m_size*m_rank,3*natom/m_size*(m_rank+1) ! modes
        write(*,*)mj,m_rank
      
        read(551,rec=mj)eigen(:)
 
        eigen(:)=eigen(:)/mass(:)

        tmp=0.0d0
        do i=1,length0
          tmp=tmp+dphir(i)*eigen(vd1(i))*eigen(vd2(i))
        enddo
         
        do i=1,length1
          tmp = tmp + phir(i,1) * eigen(v1_1(i)) * eigen(v2_1(i))
          tmp = tmp + phir(i,2) * eigen(v1_2(i)) * eigen(v2_2(i))
          tmp = tmp + phir(i,3) * eigen(v1_3(i)) * eigen(v2_3(i))
        enddo
       
        gamma1(mj)=-tmp/6.0d0/ww(mj)
      enddo
  
 
      call MPI_ALLREDUCE(gamma1(:), allgamma(:), 3*natom, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, m_ierror)
     
      if(m_rank==0) then
        write(*,*)'m_ierror',m_ierror
        if(trim(isrmean)=='IS') then 
          open(333,file='gruneisen_IS.dat',status='unknown')
        else  
          open(333,file='gruneisen_MAIN.dat',status='unknown')

        endif
        do mj=1+3,3*natom
          write(333,*)freq_to_wthz*sqrt(dabs(ww(mj)))/2/pi*thtocm,allgamma(mj)
        enddo
        close(333)
      endif
                                                          
 
      call MPI_FINALIZE(m_ierror)
 
      End program Grun
      
