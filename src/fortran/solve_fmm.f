
	subroutine solve_fmm(TOL,pot,ltree,npbox,norder,
	1    nboxes,nlevels,itree,iptr,fvals,boxsize,centers)
	implicit real *8 (a-h,o-z)
	integer iptr(9)
	integer npols, nd

	integer itree(ltree)
	real *8 boxsize(0:nlevels)
	real *8 fvals(1,npbox,nboxes),centers(2,nboxes)
	real *8 pot(1,npbox,nboxes)
	real *8, allocatable :: xref(:,:)
	real *8 xyztmp(3),rintl(0:200),umat,vmat,wts
	real *8 timeinfo(6),tprecomp(3)

	complex *16 zpars

c       real *8 potmat(np,nboxes)

	
	real *8, allocatable :: potex(:,:,:)
	complex *16 ima,zz,ztmp,zk

	real *8 alpha,beta,targ(2)

	character *1 type
	data ima/(0.0d0,1.0d0)/

	external fgaussn,fgauss1
	logical flag
	
	eps = TOL

	npols = norder*norder
	nd = 1
	call prini(6,13)
	zk = ima
	done = 1
	pi = atan(done)*4
c       
c       initialize function parameters
c       
	delta = 4d-4
	boxlen = 1.0d0
	

	nd = 1
c       Gaussian source centers and variance
	
	norder = 8
	iptype = 1
	eta = 1.0d0
	
	eps = TOL
	call prin2('eps=*',eps,1)
	call cpu_time(t1)
C       $      t1 = omp_get_wtime()
	
	call prinf('nboxes=*',nboxes,1)
	call prinf('nlevels=*',nlevels,1)
	call prinf('ltree=*',ltree,1)
	
	call prin2('boxsize=*',boxsize(0),1)
	call prin2('centers=*',centers(1,1),2)
	
c       
c       convert values to coefs
c       
	
cccc    npols = norder*(norder+1)*(norder+2)/6

	do i=1,nboxes
	   do j=1,npbox
	      pot(1,j,i) = 0
	   enddo
	enddo

	type = 'f'
	iperiod = 0
	call cpu_time(t1) 
C       $     t1 = omp_get_wtime()
	
c       print *, "boxsize = ", boxsize
c       print *, "npbox = ", npbox
	call lbfmm2d(nd,eps,iperiod,nboxes,nlevels,ltree,
	1    itree,iptr,norder,npols,type,fvals,centers,boxsize,npbox,
	2    pot,timeinfo)
	call cpu_time(t2) 
C       $     t2 = omp_get_wtime()      
	call prin2('time taken in bfmm=*',t2-t1,1)
	
	nlfbox = 0
	do ilevel=1,nlevels
	   do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
	      if(itree(iptr(4)+ibox-1).eq.0) nlfbox = nlfbox+1
	   enddo
	enddo
	call prinf('nlfbox=*',nlfbox,1)
	call prinf('ntotal=*',nlfbox*npbox,1)
ccc     d = 0
ccc     do i = 1,6
ccc     d = d + timeinfo(i)
ccc     enddo
	
	call prin2('speed in pps=*',
	1    (npbox*nlfbox+0.0d0)/(t2-t1),1)
	
	end
c       
c       
c       
c       



	subroutine exactnd1(targ,dpars,pot)

	implicit real*8 (a-h,o-z)
	real*8 targ(2)
	real*8 pot
	real*8 gf(2),dpars(*)
c       

	pot=0.0d0


	ng=4


	do i=1,ng
	   idp = (i-1)*3

	   dx = targ(1)-dpars(idp+1)
	   dy = targ(2)-dpars(idp+2)
	   r2 = dx*dx+dy*dy

	   sigma = dpars(idp+3)
	   
	   pot=pot+exp(-r2/sigma)
	enddo


	return
	end
