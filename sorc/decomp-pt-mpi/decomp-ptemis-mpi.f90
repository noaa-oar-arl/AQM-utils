  program decomp_ptemis_mpi
!-------------------------------------------------------------------------------------
! split point source emission file into sub-domain based on 
! runtime configure setting.
!
! by Youhua Tang, Sep 2022

  use netcdf
  use mpi
  
  integer, allocatable, dimension (:) :: ixt, jyt, mtime, mindx
  real, allocatable, dimension (:) :: tlat,tlon,stkdm,stkout
  real, allocatable, dimension (:,:) :: xlat, xlon, ain, bout
  integer idlists(200),idlists2(200),dimids(2),ibuff(4)
  character*200 afile,topofile
  character*20  argval,vname
  integer, parameter :: kcompress=4
  
  call get_environment_variable('NX',argval,lenval,istatus,.true.)
  if(istatus.ne.0) then
   print*,'failed to get NX'
   stop
  endif
  read(argval,*)nx

  call get_environment_variable('NY',argval,lenval,istatus,.true.)
  if(istatus.ne.0) then
   print*,'failed to get NY'
   stop
  endif
  read(argval,*)ny

  call get_environment_variable('LAYOUT_X',argval,lenval,istatus,.true.)
  if(istatus.ne.0) then
   print*,'failed to get LAYOUT_X'
   stop
  endif
  read(argval,*)layout_x

  call get_environment_variable('LAYOUT_Y',argval,lenval,istatus,.true.)
  if(istatus.ne.0) then
   print*,'failed to get LAYOUT_Y'
   stop
  endif
  read(argval,*)layout_y

  call get_environment_variable('TOPO',topofile,lenval,istatus,.true.)
  if(istatus.ne.0) then
   print*,'failed to get TOPO'
   stop
  endif

  call get_environment_variable('PT_IN',afile,lenval,istatus,.true.)
  if(istatus.ne.0) then
   print*,'failed to get PT_IN'
   stop
  endif
  
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, npe, ierr)
  
  if(my_rank.eq.0) then
    ibuff(1)=nx
    ibuff(2)=ny
    ibuff(3)=layout_x
    ibuff(4)=layout_y
  endif
  call MPI_BCAST(ibuff,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if(ierr.ne.MPI_SUCCESS) then
   print*,'mpi_bcast error'
   stop
  endif
  if(my_rank.gt.0) then
   nx=ibuff(1)
   ny=ibuff(2)
   layout_x=ibuff(3)
   layout_y=ibuff(4)
  endif    
  
  if(npe.ne.layout_x*layout_y) then
    print*,'inconsistent npe ',npe,layout_x,layout_y
    stop
  endif
  
  call MPI_BCAST(topofile,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  if(ierr.ne.MPI_SUCCESS) stop
  call MPI_BCAST(afile,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  if(ierr.ne.MPI_SUCCESS) stop
  
! read topo file lat/lon info
!  print*,'read topofile',my_rank,trim(topofile) 
  call check(nf90_open(trim(topofile),nf90_nowrite,ncid2))
  call check(nf90_inq_dimid(ncid2,'grid_xt',iddim_lon))
  call check(nf90_inquire_dimension(ncid2,iddim_lon,len=nx2))
  call check(nf90_inq_dimid(ncid2,'grid_yt',iddim_lat))
  call check(nf90_inquire_dimension(ncid2,iddim_lat,len=ny2))
  if(nx.ne.nx2.or.ny.ne.ny2) then
   print*,'inconsistent topo file dimension'
   stop
  endif 
  call check(nf90_inq_varid(ncid2,'grid_latt',idvar_geolat))
  call check(nf90_inq_varid(ncid2,'grid_lont',idvar_geolon))
  allocate(xlon(nx,ny),xlat(nx,ny))
  call check(nf90_get_var(ncid2,idvar_geolon,xlon))
  call check(nf90_get_var(ncid2,idvar_geolat,xlat))
  call check(nf90_close(ncid2))

  distnear=0.7071*haversine(xlat(1,1),xlon(1,1),xlat(2,1),xlon(2,1))  
  if(distnear.le.1) then
    print*,'distnear too small ',distnear,xlat(1,1),xlat(2,1),xlon(1,1),xlon(2,1)
    stop
  endif
    
  iperde=nx/layout_x
  iextra=mod(nx,layout_x)
  iextra1=iextra/2+mod(iextra,2)  ! index of extra from boundary
  iextra2=iextra/2
  
  jperde=ny/layout_y
  jextra=mod(ny,layout_y)
  jextra1=jextra/2+mod(jextra,2)
  jextra2=jextra/2  
  
!  print*,'iperde,iextra,jperde,jextra=',iperde,iextra,jperde,jextra
  
  nowpe=0
  find_pe: do j=1, layout_y
   if(j.eq.1) then
    jstart=1
   else 
    jstart=jend+1 
   endif
    
   if (j.le.jextra1.or.j.ge.(layout_y-jextra2+1)) then
    jend=jstart+jperde
   else
    jend=jstart+jperde-1
   endif 
   
   do i=1, layout_x
    if(i.eq.1) then
     istart=1
    else
     istart=iend+1 
    endif
     
    if (i.le.iextra1.or.i.ge.(layout_x-iextra2+1)) then
     iend=istart+iperde
    else
     iend=istart+iperde-1
    endif
    
    if(my_rank.eq.nowpe) exit find_pe
    nowpe=nowpe+1
   enddo
  enddo find_pe
  if(nowpe.ge.npe) then
   print*,'can not match PE ',my_rank, nowpe
   stop
  endif 
! print*,'pe,is,ie,js,je=',nowpe,istart,iend,jstart,jend 
  
!--read input PT file
  call check(nf90_open(trim(afile),nf90_nowrite, ncid))
  call check(nf90_inq_dimid(ncid,'nlocs',iddim_stack))
  call check(nf90_inquire_dimension(ncid,iddim_stack,len=nstack))
  call check(nf90_inq_dimid(ncid,'time',iddim_time))
  call check(nf90_inquire_dimension(ncid,iddim_time,len=ntime))
  
  allocate(tlat(nstack),tlon(nstack),ixt(nstack),mindx(nstack),jyt(nstack), &
   stkdm(nstack),ain(nstack,ntime),mtime(ntime))
  
  call check(nf90_inq_varid(ncid,'LATITUDE',idvar))
  call check(nf90_get_var(ncid,idvar,tlat))
  call check(nf90_inq_varid(ncid,'LONGITUDE',idvar))
  call check(nf90_get_var(ncid,idvar,tlon))

! -widen two grid cells  
  jstart2=max(1,jstart-2)
  jend2=min(ny,jend+2)
  istart2=max(1,istart-2)
  iend2=min(nx,iend+2)
  
  m=0
  do n=1,nstack
    search_loop: do j = jstart2, jend2
      do i = istart2, iend2
	if(haversine(xlat(i,j),xlon(i,j),tlat(n),tlon(n)).le.distnear)  exit search_loop
      enddo
    enddo search_loop
    if(i.le.iend2.and.j.le.jend2) then
      m=m+1
      mindx(m)=n
    endif  
  enddo  
  mtotal=m
  if(mtotal.eq.0) then
   mtotal=1
   mindx(1)=1
  endif
  print*,'pe,mtotal,distnear,istart2,iend2,jstart2,jend2=',my_rank,mtotal,distnear,istart2,iend2,jstart2,jend2
  
!--- create output file
  allocate(stkout(mtotal),bout(mtotal,ntime))
    
  write(afile,"('PT/pt-',i4.4,'.nc')")my_rank
  call check(nf90_create(trim(afile),cmode=IOR(IOR(NF90_CLOBBER,NF90_NETCDF4),NF90_CLASSIC_MODEL),ncid=ncid2))
  
  call check(nf90_def_dim(ncid2,'nlocs',len=mtotal,dimid=iddim2_stack))
  call check(nf90_def_dim(ncid2,'time',len=ntime,dimid=iddim2_time))
     
  call check(nf90_inq_varids(ncid,nvars,idlists))
  
  do n=1,nvars
    call check(nf90_inquire_variable(ncid,idlists(n),vname,xtype=itype,ndims=ndims,&
      dimids=dimids,natts=natts))
    if(ndims.eq.1) then  
     if(dimids(1).eq.iddim_stack) then
      call check(nf90_def_var(ncid2,vname,itype,[iddim2_stack],idlists2(n),deflate_level=kcompress))
     else if (dimids(1).eq.iddim_time) then
      call check(nf90_def_var(ncid2,vname,itype,[iddim2_time],idlists2(n),deflate_level=kcompress))
     else
      print*,'unknown dimids ',dimids,vname,idlists(n)
      stop
     endif
    else if(ndims.eq.2) then
      if(dimids(2).eq.iddim_time.and.dimids(1).eq.iddim_stack) then
       call check(nf90_def_var(ncid2,vname,itype,[iddim2_stack,iddim2_time],idlists2(n),shuffle=.true.,deflate_level=kcompress))
      else
       print*,'1, please check dimids ',dimids,iddim_time,iddim_stack,vname,idlists(n)
       stop
      endif
    else
       print*,'2, please check dimids ',dimids,vname,idlists(n)
       stop
    endif
    do nat=1,natts  ! copy attributes
     call check(nf90_inq_attname(ncid,idlists(n),nat,afile))
     call check(nf90_copy_att(ncid,idlists(n),trim(afile),ncid2,idlists2(n)))
    enddo 
  enddo
  call check(nf90_enddef(ncid2)) 

! start output  
  do n=1,nvars
    call check(nf90_inquire_variable(ncid,idlists(n),vname,xtype=itype,ndims=ndims,&
      dimids=dimids,natts=natts))    
    if(ndims.eq.1) then      
     if(dimids(1).eq.iddim_stack) then  ! stack variables
       if(itype.ne.NF90_REAL) then
         print*,'wrong type  ',trim(vname),itype
	 stop
       endif
       call check(nf90_get_var(ncid,idlists(n),stkdm))
       do m=1,mtotal
        stkout(m)=stkdm(mindx(m))
       enddo
       call check(nf90_put_var(ncid2,idlists2(n),stkout))
     else if (dimids(1).eq.iddim_time) then
       if(itype.ne.NF90_INT) then
         print*,'wrong type  ',trim(vname),itype
	 stop
       endif
       call check(nf90_get_var(ncid,idlists(n),mtime))
       call check(nf90_put_var(ncid2,idlists2(n),mtime))
     endif
    else if (ndims.eq.2) then
      if(dimids(2).eq.iddim_time.and.dimids(1).eq.iddim_stack) then 
        call check(nf90_get_var(ncid,idlists(n),ain))
	do m=1,mtotal
	  bout(m,:ntime)=ain(mindx(m),:ntime)
	enddo
	call check(nf90_put_var(ncid2,idlists2(n),bout))
      else
        print*,'wrong dimids ',vname,dimids(:2)
	stop
      endif
    else
      print*,'ndims wrong',ndims,vname
      stop
    endif
   enddo
   call check(nf90_close(ncid2))
   call MPI_Finalize(ierr)  
   
   contains
       subroutine check(status)
       integer, intent ( in) :: status

       if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
       end if
      end subroutine check
     end

      function to_radian(degree) result(rad)
          ! degrees to radians
          real,intent(in) :: degree
          real, parameter :: deg_to_rad = atan(1.0)/45 ! exploit intrinsic atan to generate pi/180 runtime constant
          real :: rad
 
          rad = degree*deg_to_rad
      end function to_radian
 
      function haversine(deglat1,deglon1,deglat2,deglon2) result (dist)
          real,intent(in) :: deglat1,deglon1,deglat2,deglon2
          real :: a,c,dist,dlat,dlon,lat1,lat2
          real,parameter :: radius = 6372.8 ! in km
 
          dlat = to_radian(deglat2-deglat1)
          dlon = to_radian(deglon2-deglon1)
          lat1 = to_radian(deglat1)
          lat2 = to_radian(deglat2)
          a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
          c = 2*asin(sqrt(a))
          dist = radius*c
      end function haversine
