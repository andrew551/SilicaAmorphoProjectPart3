      Program Grun
       
      implicit none
 
  
       real*8, parameter :: pi=4.0d0*datan(1.0d0), eJ=1.602176487d-19, a_to_m=1.0d-10, &
                     amu=1.660599d-27, freq_to_wthz=dsqrt(eJ/amu/a_to_m**2)/1.0d12, &
                     togpa=eJ/a_to_m**3/1.0d9, &
                     hbar = 1.054571d-34,& ! J*s
                     kb=1.380649d-23, &
                     thtocm=1.0d0/0.0299792458d0
      
      integer :: i,j,k,natom,mj,cnt,ia,al,be,ic,ga,length,length0,reclen,typ,ii,jj,kk,length_dm
      integer :: m_rank, m_size, m_ierror,itp
      real*8 :: var,qq,r(3),mm(10),eng,T,arg,box(3),x1,x2,vol,bulk,shift(3),rtemp(3),dx(3),rad(648,648,3),dum,var2(2),fc
      real*8, allocatable :: vaac(:),vaaa(:),rrr(:,:)
      integer, allocatable :: v1_1(:),v2_1(:),v3_1(:),v1_2(:),v2_2(:),v3_2(:),&
               v1_3(:),v2_3(:),v3_3(:),vd1(:),vd2(:),vd3(:),bb(:)
      
      real*8, allocatable :: eigen(:),mass(:),rr(:),rmean(:),ww(:),dm(:),gamma1(:),allgamma(:,:),nj(:),CV(:)
 
      character*100 path,str,factor, path_model,path_w,path_eigen,isrmean
	  
	  call getarg(1,factor)
	  read(factor,*)fc
	  write(*,*)fc
  
    open(11,file='inputs',status='old')
    read(11,*)isrmean
    write(*,*)isrmean
  
    read(11,*)path_model
    write(*,*)path_model
    read(11,*)path_w
    write(*,*)path_w
    read(11,*)path_eigen
    write(*,*)path_eigen
     
    close(11)

      open(unit=45,file=trim(path_model),status='old')
      read(45,*)
      read(45,*)
      read(45,*) natom

      write(*,*)natom
  
      allocate(allgamma(natom*3,2),ww(natom*3),cv(natom*3),nj(natom*3))
      
      
      open(55,file='frequencies.dat',status='old')
      do mj=1,3*natom
        read(55,*)i,ww(mj)
      enddo
      close(55)
       
!    
 
      open(333,file='gruneisen_MAIN.dat',status='old')
      do mj=1+3,3*natom
        read(333,*)dum,allgamma(mj,1)
      enddo
      close(333)
	  
	  open(333,file='gruneisen_IS.dat',status='old')
      do mj=1+3,3*natom
        read(333,*)dum,allgamma(mj,2)
      enddo
      close(333)
    
  
        open(666,file='gruneisen_T_'//trim(factor)//'_.dat',status='unknown')
        do i=1,1000
          T=i
          nj=0.0d0 ! boson factor at mode freq []
          do mj=1+3,3*natom
            eng=dsqrt(ww(mj))*freq_to_wthz*1d12*hbar ! in [J]
            arg=eng/(T*kb) ! dimensionless
            if(dabs(arg)>100)arg=100*arg/dabs(arg)
            nj(mj)=1.0d0/(dexp(arg)-1.0d0)
         !   write(*,*)nj(mj)
          enddo

          CV=0.0d0 ! capacity at mode freq positions [J/K]
          do mj=4,3*natom
            CV(mj)=(hbar*dsqrt(ww(mj))*freq_to_wthz*1d12)**2/ &
                   (kb*T*T)*nj(mj)*(1.0d0+nj(mj))  ! [J/K]
        !     write(*,*)cv(mj)
          enddo
          
       
          write(666,*)i,sum(CV(4:3*natom)*(allgamma(4:3*natom,1)+fc*allgamma(4:3*natom,2)))/sum(CV(4:3*natom))
         enddo
    
                                                          
 
 
      End program Grun
      
