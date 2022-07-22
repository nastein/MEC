program em_2body
   use mympi
   use mc_module
   use dirac_matrices
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmpi=139.5d0,xmd=1236.0d0,xmn=939.0d0
   real*8, parameter :: xmmu=105.658357,q2max_v=4.0d0*1.e6
   integer*8, allocatable :: irn(:),irn0(:)
   integer*4 :: nw,nev,i,xA,i_fg,np,ne,nq,nenu,iq,ie,i_fsi
   integer*4 :: nwlk=128
   real*8 :: wmax,qval,xpf,hw,w
   real*8 :: ee,ef,mlept,thetalept,cos
   real*8 :: sig,sigA,sig_err,sigA_err
   real*8 :: enu_max,q2max,he,q2max_c
   real*8 :: ti,tf
   real*8 :: ti1, ti2 !isospins of initial nucleon pair
   real*8, allocatable :: enu_v(:)
   real*8 :: sig_pi_pi, sig_pi_ni
   real*8 :: sig_pi_pi_err, sig_pi_ni_err
   character*40 :: nk_fname
   logical :: pwia
! initialize mpi
   call init0()

! read input
    open(unit=7, file='C12_2700_15p0_2b_pp+nn.out')
    open(unit=8, file='C12_2700_15p0_2b_pn+np.out')

   if (myrank().eq.0) then
      read(5,*) nev
      read(5,*) ee,thetalept
      !read(5,*) enu,cos
      read(5,*) wmax,nw      
      read(5,*) xpf
      read(5,*) xA
      read(5,*) i_fg
      read(5,*) nk_fname
      read(5,*) np,ne
      read(5,*) i_fsi
      read(5,*) pwia
   endif
   call bcast(nev)
   call bcast(ee)
   call bcast(thetalept)
!   call bcast(cos)   
   call bcast(wmax)
   call bcast(nw)
   call bcast(xpf)
   call bcast(xA)
   call bcast(i_fg)
   call bcast(nk_fname)
   call bcast(np)
   call bcast(ne)
   call bcast(i_fsi)
   call bcast(pwia)
   ti=MPI_Wtime()


    allocate(irn0(nwlk))
    do i=1,nwlk
       irn0(i)=19+i
    enddo
    if (myrank().eq.0) then
       write (6,'(''number of cpus ='',t50,i10)') nproc()
       if (mod(nwlk,nproc()).ne.0) then
          write(6,*)'Error: nwalk must me a multiple of nproc'
          stop
       endif
    endif
    nwlk=nwlk/nproc()
    
    allocate(irn(nwlk))
    irn(:)=irn0(myrank()*nwlk+1:myrank()*nwlk+nwlk)
    thetalept=thetalept/180.0d0*pi
    !thetalept=acos(cos)

   !...this is for the total xsec
   hw=(wmax)/dble(nw)

   call dirac_matrices_in(xmd,xmn,xmpi)

   !Cross section for pp and nn initial pairs
   ti1=0.5 !+ 1/2 isospin
   ti2=0.5 !+ 1/2 isospin
   call set_isospins(ti1,ti2)
   call mc_init(i_fg,pwia,i_fsi,irn,nev,nwlk,xpf,thetalept,xmpi,xmd,xmn,xA,np,ne,nk_fname)

   do i=1,nw
      w=+dble(i-0.5d0)*hw
      ef=ee-w
      call mc_eval(ee,w,sig_pi_pi,sig_pi_pi_err)

      if (myrank().eq.0) then
         write(6,*)'Computing for pp and nn initial states'
         write(6,*)w,ee-w, sig_pi_pi*2.0d0, sig_pi_pi_err !multiply cross section by two because pp is same as nn
         write(7,*)w,ee-w, sig_pi_pi*2.0d0, sig_pi_pi_err
         flush(7)
      endif
   enddo
   close(7)

      !Cross section for pp and nn initial pairs
   ti1=0.5 !+ 1/2 isospin
   ti2=-0.5 !- 1/2 isospin
   call set_isospins(ti1,ti2)
   
   do i=1,nw
      w=+dble(i-0.5d0)*hw
      ef=ee-w
      call mc_eval(ee,w,sig_pi_ni,sig_pi_ni_err)

      if (myrank().eq.0) then
         write(6,*)'Computing for pn and np initial states'
         write(6,*)w,ee-w, sig_pi_ni, sig_pi_ni_err
         write(8,*)w,ee-w, sig_pi_ni, sig_pi_ni_err
         flush(8)
      endif
   enddo
   close(8)

   close(9)
   tf=MPI_Wtime()
   if (myrank().eq.0) then
      write(6,*)'Elapsed time is',tf-ti
   endif
   call done()   

end program
