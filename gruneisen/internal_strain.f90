      Program istress
      implicit none
 
    
      integer :: i,j,k,natom,nmode,reclen,step
      real*8 dum      
      real*8, allocatable :: xx(:),xx0(:),emode(:),smode(:),mass(:)
      character*6 ch
      character*100 path_model,path_w,path_eigen,isrmean
  
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
   
      nmode=3*natom
 
   
      allocate(xx(nmode),xx0(nmode),emode(nmode),smode(nmode),mass(nmode))
  
      ch=''
      do while(trim(ch) /= 'Atoms')
        read(45,*)ch
      enddo
      read(45,*)
  
      do i=1,natom
        read(45,*)j,k,dum,xx0(3*j-2:3*j)
        if(k==1)mass(3*j-2:3*j)=dsqrt(15.9994d0)
        if(k==2)mass(3*j-2:3*j)=dsqrt(28.086d0)
      enddo
      close(45)

      inquire(iolength=reclen)emode(:)  
      open(13,file=trim(path_eigen),form='unformatted',access='direct',status='old', &
                    recl=reclen,action='read')
      
      smode=0.0d0
      do j=1,nmode
        !write(*,*)j
        read(13,rec=j)emode(:)
        smode(j)=dot_product(emode(:)*mass(:),xx0(:))
      enddo
 
      xx=0.0d0
      do j=1,nmode
        if(j<4)cycle
        read(13,rec=j)emode(:)
        xx(:)=xx(:)-emode(:)/mass(:)*smode(j)
      enddo
  
      close(13)
 
      open(11,file='rmean.dat', status='unknown')
      do i=1,natom
        write(11,'(I8,6E15.5)')i,xx(i*3-2:i*3),xx0(i*3-2:i*3)
      enddo
      close(11)
    
 
      End program istress
     
  
  
      
      
