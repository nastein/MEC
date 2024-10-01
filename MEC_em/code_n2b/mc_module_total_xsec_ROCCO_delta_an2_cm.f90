module mc_module
   implicit none
   integer*4, private, save :: nev,xA,i_fg,i_fsi,np,ne,nwlk,npot,np_del,total_isospin
   integer*4, private, parameter :: neq=10000,nvoid=10
   real*8, private, save ::  xpf
   real*8, private, save:: xmpi,xmd,xmn,norm,norm0,norm1,thetalept
   real*8, private, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0,xmp=938.0d0,ppmax=1.0d0*1.e3
   real*8, private,parameter :: alpha=1.0d0/137.0d0
   real*8, private, allocatable :: pv(:),dp(:,:),ep(:),dp1(:,:),dp0(:,:)
   real*8, private, allocatable :: kin(:),pot(:),pdel(:),pot_del(:)
   integer*8, private, allocatable, save :: irn(:)
contains

subroutine mc_init(i_fg_in,pwia,i_fsi_in,irn_in,nev_in,nwlk_in,xpf_in,thetalept_in,xmpi_in,xmd_in,xmn_in,xA_in, &
     &   total_isospin_in,nk_fname_in)
  use mathtool
   implicit none
   integer*8 :: irn_in(nwlk_in)
   integer*4 :: nev_in,nwlk_in,xA_in,i_fg_in,np_in,i,j,ne_in,np0,ne0,ien,total_isospin_in
   integer*4 :: ipot,i_fsi_in
   real*8 :: xpf_in,xmpi_in,xmd_in,xmn_in,mlept_in,hp,he,thetalept_in
   character*40 :: nk_fname_in
   logical :: pwia

   nev=nev_in
   nwlk=nwlk_in
   xpf=xpf_in
   xmpi=xmpi_in
   xmd=xmd_in
   xmn=xmn_in
   thetalept=thetalept_in
   xA=xA_in
   total_isospin=total_isospin_in
   i_fg=i_fg_in
   i_fsi=i_fsi_in

   allocate(irn(nwlk))
   irn(:)=irn_in(:)
   if(i_fg.ne.1) then
      open(unit=8,file=nk_fname_in,status='unknown',form='formatted')
      read(8,*) np
      allocate(pv(np),dp(np,np),dp1(np,np),dp0(np,np))
      do i=1,np
         do j=1,np
           read(8,*) pv(i),pv(j),dp(i,j),dp1(i,j),dp0(i,j) 
           !print*, pv(i),pv(j),dp(i,j),dp1(i,j),dp0(i,j) 
         enddo  
      enddo
      close(8)
   endif
   pv=pv*hbarc       
   dp=dp/hbarc**6
   dp1=dp1/hbarc**6
   dp0=dp0/hbarc**6
  
   norm=0.0d0
   norm1=0.0d0  
   norm0=0.0d0  
   dp=dp/(2.0d0*pi)**6
   dp1=dp1/(2.0d0*pi)**6
   dp0=dp0/(2.0d0*pi)**6
   do i=1,np
      do j=1,np
         norm=norm+dp(i,j)*pv(i)**2*pv(j)**2*(4.0d0*pi*(pv(2)-pv(1)))**2
         norm1=norm1+dp1(i,j)*pv(i)**2*pv(j)**2*(4.0d0*pi*(pv(2)-pv(1)))**2
         norm0=norm0+dp0(i,j)*pv(i)**2*pv(j)**2*(4.0d0*pi*(pv(2)-pv(1)))**2
       enddo
   enddo    

   dp=dp/norm*(4.0d0*pi*xpf**3/3.0d0)**2
   dp1=dp1/norm*(4.0d0*pi*xpf**3/3.0d0)**2
   dp0=dp0/norm*(4.0d0*pi*xpf**3/3.0d0)**2
   

   norm=0.0d0
   norm1=0.0d0  
   norm0=0.0d0 
   do i=1,np
      do j=1,np
         norm=norm+dp(i,j)*pv(i)**2*pv(j)**2*(4.0d0*pi*(pv(2)-pv(1)))**2
         norm1=norm1+dp1(i,j)*pv(i)**2*pv(j)**2*(4.0d0*pi*(pv(2)-pv(1)))**2
         norm0=norm0+dp0(i,j)*pv(i)**2*pv(j)**2*(4.0d0*pi*(pv(2)-pv(1)))**2
       enddo
   enddo  
   ! this needs to be updated

   if(i_fsi.eq.1) then
      open(8, file='pke/realOP_12C_EDAI.dat')
      read(8,*) npot
      allocate(kin(npot),pot(npot))
      do ipot=1,npot
         read(8,*)kin(ipot),pot(ipot)
         pot(ipot)=pot(ipot)
         !....kin and pot are in MeV
      enddo
   endif

   open(10, file='rho_1.dat')
   read(10,*) np_del
   allocate(pdel(np_del),pot_del(np_del))
   do i=1,np_del
      read(10,*) pdel(i),pot_del(i)
   enddo
      
   return
end subroutine   


subroutine mc_eval(ee,w,sig_avg_tot,sig_err_tot)
  use mathtool
  use mympi
   implicit none
   integer*4 :: i,ie1,ie2,ne1,ne2,j
   real*8 :: w,emax,ee
   real*8 :: sig_o(nwlk)
   real*8 :: sig_avg,sig_err
   real*8 :: sig_avg_tot,sig_err_tot

   real*8 :: xpmax,qval_in,sig
   real*8 :: wmax,q2,q2max,nk(np,np),nk_norm
   integer*4 :: i_acc,i_avg,i_acc_tot,i_avg_tot
   integer*4 :: ip1_o(nwlk),ie1_o(nwlk),ip2_o(nwlk),ie2_o(nwlk)
   integer*4 :: ip1_n(nwlk),ip2_n(nwlk),ie1_n(nwlk),ie2_n(nwlk)
   real*8 ::  g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk) 

   sig_avg=0.0d0
   sig_err=0.0d0
   i_acc=0
   i_avg=0
   g_o=0.0d0
   !if (i_fg.eq.1) then
   !   xpmax=pv(np)
   !   emax=1.0d0
   !else
      xpmax=pv(np)!-pv(1)
   !   emax=ep(ne)-ep(1)
   !endif
   
   !

   if(total_isospin.eq.1) then 
      nk = dp1 
      nk_norm = norm1
      !print*,'PP and NN initial states'
   else if(total_isospin.eq.0) then
      nk = dp0
      nk_norm = norm0
      !print*,'NP and PN initial states'
   else
      nk = dp  
      nk_norm = norm
      !print*,'ALL initial states'
   endif
   
   do i=1,nwlk
      call setrn(irn(i))
      do while(g_o(i).le.0.0d0)
         ip1_o(i)=1+int(np*ran())
         ip2_o(i)=1+int(np*ran())
         call g_eval(pv(ip1_o(i)),pv(ip2_o(i)),nk(ip1_o(i),ip2_o(i)),xpmax,nk_norm,g_o(i))
      enddo
      call getrn(irn(i))
   enddo
   
      
   do i=1,nev
      do j=1,nwlk
         call setrn(irn(j))
         ip1_n(j)=nint(ip1_o(j)+0.05d0*np*(-1.0d0+2.0d0*ran()))
         ip2_n(j)=nint(ip2_o(j)+0.05d0*np*(-1.0d0+2.0d0*ran()))
         if (ip1_n(j).le.np.and.ip1_n(j).ge.1.and.ip2_n(j).le.np.and.ip2_n(j).ge.1) then
            call g_eval(pv(ip1_n(j)),pv(ip2_n(j)),nk(ip1_n(j),ip2_n(j)),xpmax,nk_norm,g_n(j))
         else
            g_n(j)=0.0d0
         endif

         if (g_n(j)/g_o(j).ge.ran()) then
            ip1_o(j)=ip1_n(j)
            ip2_o(j)=ip2_n(j)
            g_o(j)=g_n(j)
            i_acc=i_acc+1
         endif
         if (i.ge.neq.and.mod(i,nvoid).eq.0) then
            if(ip1_o(j).gt.np.or.ip2_o(j).gt.np) cycle
            call f_eval(ee,pv(ip1_o(j)),pv(ip2_o(j)),ip1_o(j),ip2_o(j),w,sig_o(j),nk(ip1_o(j),ip2_o(j)))    
            sig_o(j)=sig_o(j)/g_o(j)
            sig_avg=sig_avg+sig_o(j)
            sig_err=sig_err+sig_o(j)**2
            i_avg=i_avg+1
         endif
         call getrn(irn(j))
      enddo
   enddo
   call addall(sig_avg,sig_avg_tot) 
   call addall(sig_err,sig_err_tot)
   call addall(i_avg,i_avg_tot) 
   call addall(i_acc,i_acc_tot) 
   if (myrank().eq.0) then
      sig_avg_tot=sig_avg_tot/dble(i_avg_tot)
      sig_err_tot=sig_err_tot/dble(i_avg_tot)
      sig_err_tot=sqrt((sig_err_tot-sig_avg_tot**2)/dble(i_avg_tot-1))
      write(6,*)'acceptance=',dble(i_acc_tot)/dble(nev*nwlk*nproc())
   endif
 
   return
end subroutine   

subroutine f_eval(ee,p1,p2,ip1,ip2,w,sig,nk)
  use mathtool
  implicit none
  integer*4 :: ip1,ip2
  real*8 :: ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,ee,eef,w,q2
  real*8 :: qval,cos_theta,jac_c,tan2
  real*8 :: v_ll,v_t,nk
  real*8 :: tnl2,sig0,delta,rho,rhop
  real*8 :: r_now(5),sig,qp(4)

  ctpp1=-1.0d0+2.0d0*ran()
  ctp2=-1.0d0+2.0d0*ran()
  phip2=2.0d0*pi*ran()
  ctp1=-1.0d0+2.0d0*ran()
  phip1=2.0d0*pi*ran()
  sig=0.0d0
  eef=ee/hbarc
  cos_theta=cos(thetalept)
  tan2=(1.0d0-cos_theta)/(1.0d0+cos_theta)
  !.....compute sigma_mott [ fm^2 --> mb ]
  sig0=alpha**2/2.0d0/(1.0d0-cos_theta)/eef**2/tan2
  sig0=sig0*10.0d0
  q2=2.0d0*ee*(ee-w)*(1.0d0-cos_theta)
  if(q2.lt.0.0d0) stop
  qval=sqrt(q2+w**2)
  v_ll=(-q2/qval**2)**2
  v_t=(q2/qval**2/2.0d0+tan2)
  call int_eval(ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,ip1,ip2,w,qval,r_now,nk)
  r_now(:)=r_now(:)*2.0d0**3*(2.0d0*pi)**2
  sig=sig0*(v_ll*r_now(1)+v_t*r_now(4))*1.e9

  return

end subroutine f_eval

subroutine int_eval(ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,ip1,ip2,w,qval,r_now,nk)
   use dirac_matrices         
   use mathtool
   implicit none
   real*8, parameter :: lsq=0.71*1.e6,l3=3.5d0*1.e6,xma2=1.1025d0*1.e6
   real*8, parameter :: fstar=2.13d0,eps=10.0d0,e_gs=-92.16,e_bg=-64.75
   integer*4 :: ip1,ip2
   real*8 :: w,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,stpp1,stp1,stp2
   real*8 :: at,bt,vt,par1,par2,pp1,den,jac,arg,qval
   real*8 :: q2,rho,norm,ca5,cv3,gep,nk
   real*8 :: p1_4(4),p2_4(4),pp1_4(4),pp2_4(4),k2_4(4),k1_4(4),q_4(4),pp_4(4)
   real*8 :: k2e_4(4),k1e_4(4)
   real*8 :: pp1_4cm(4),pp2_4cm(4),phipp1_cm,ctpp1_cm
   real*8 :: vcm(3),vcm_mag,gammacm,uhatcm(3) 
   real*8 :: stpp1_cm,E_tot,p_tot(3),p_totmag,pp1_cm_mag,lorentz_jac
   real*8 :: r_cc_pi,r_cl_pi,r_ll_pi,r_t_pi,r_tp_pi
   real*8 :: r_cc_del,r_cl_del,r_ll_del,r_t_del,r_tp_del   
   real*8 :: r_cc_int,r_cl_int,r_ll_int,r_t_int,r_tp_int  
   real*8 :: dp1,dp2 
   real*8 :: tkin_pp1,tkin_pp2, u_pp1,u_pp2
   real*8 :: dir(5),exc(5),r_now(5)

   stpp1=sqrt(1.0d0-ctpp1**2)
   stp1=sqrt(1.0d0-ctp1**2)
   stp2=sqrt(1.0d0-ctp2**2)

   p1_4(1)=sqrt(p1**2+xmn**2)
   p1_4(2)=p1*stp1*cos(phip1)
   p1_4(3)=p1*stp1*sin(phip1)
   p1_4(4)=p1*ctp1
   p2_4(1)=sqrt(p2**2+xmn**2)
   p2_4(2)=p2*stp2*cos(phip2)
   p2_4(3)=p2*stp2*sin(phip2)
   p2_4(4)=p2*ctp2


   !Compute qtilde
   !if(i_fg.eq.1) then
   !   q_4(1)=w-40.0d0
   !else
     ! q_4(1)=w-p1_4(1)-p2_4(1)-ep(ie1)+xmn-ep(ie2)+xmn+60.0d0!-u_pp1-u_pp2  
   !   q_4(1)=w+e_gs-e_bg &!-sum(p1_4(2:4)+p2_4(2:4))**2/2.0d0/(10.0d0*xmn) &
   !    & -p1_4(1)-p2_4(1)+2.0d0*xmn!-u_pp1-u_pp2 
   !endif
   q_4(1) = w
   if (q_4(1).lt.0d0) then                                                                      
      r_now=0.0d0                                                                               
      return                                                                                    
   endif
   q_4(2:3)=0.0d0
   q_4(4)=qval

   !Compute the total energy and momentum in lab frame
   E_tot = p1_4(1) + p2_4(1) + q_4(1)! + 40.0d0
   p_tot = p1_4(2:4) + p2_4(2:4) + q_4(2:4)
   p_totmag = sqrt(sum(p_tot(1:3)**2))

   !Check that we have enough energy to create the two final state particles
   if((E_tot**2 - p_totmag**2 - 4.0d0*xmn**2).lt.0.0d0) then 
      r_now(:) = 0.0d0 
      return
   endif

   !Now we go to the CM frame of pp1 and pp2 to pick momenta
   !Choose angles
   phipp1_cm=2.0d0*pi*ran() 
   ctpp1_cm=-1.0d0+2.0d0*ran()
   stpp1_cm = sqrt(1.0d0 - ctpp1_cm**2)
   !Use lorentz invariance to get pp1_cm(1) and momentum
   pp1_4cm(1) = 0.5d0*sqrt(E_tot**2 - p_totmag**2)
   pp1_cm_mag = sqrt(pp1_4cm(1)**2 - xmn**2)
   pp1_4cm(2) = pp1_cm_mag*stpp1_cm*cos(phipp1_cm)
   pp1_4cm(3) = pp1_cm_mag*stpp1_cm*sin(phipp1_cm)
   pp1_4cm(4) = pp1_cm_mag*ctpp1_cm
   !pp1_3cm = - pp2_3cm
   pp2_4cm(2:4) = -pp1_4cm(2:4)
   pp2_4cm(1) = pp1_4cm(1)

   !Ok now I have 4vecs in cm frame
   !I want to boost back to lab frame
   !velocity of cm frame
   vcm(:) = p_tot(:)/E_tot 
   vcm_mag = sqrt(sum(vcm(1:3)**2))
   uhatcm(:) = vcm(:)/vcm_mag
   gammacm = 1.0d0/sqrt(1 - vcm_mag**2)

   !Lorentz transform
   pp1_4(2:4) = (gammacm*vcm_mag*pp1_4cm(1) + (gammacm - 1.0d0)*dot_product(pp1_4cm(2:4),uhatcm))*uhatcm(:) + pp1_4cm(2:4)
   pp2_4(2:4) = (gammacm*vcm_mag*pp2_4cm(1) + (gammacm - 1.0d0)*dot_product(pp2_4cm(2:4),uhatcm))*uhatcm(:) + pp2_4cm(2:4)
   pp1_4(1) = sqrt(xmn**2 + sum(pp1_4(2:4)**2))
   pp2_4(1) = sqrt(xmn**2 + sum(pp2_4(2:4)**2))

   !....Pauli blocking
   if(sqrt(sum(pp1_4(2:4)**2)).lt.xpf) then   
      r_now=0.0d0
      return
   endif        

   if(sqrt(sum(pp2_4(2:4)**2)).lt.xpf) then
      r_now=0.0d0
      return
   endif

   !Jacobian
   lorentz_jac = pp1_cm_mag/2.0d0/pp1_4cm(1)

   !...define pion momenta
   k1_4(:)=pp1_4(:)-p1_4(:)
   k2_4(:)=q_4(:)-k1_4(:)
   k1e_4(:)=pp2_4(:)-p1_4(:)
   k2e_4(:)=q_4(:)-k1e_4(:)

   !Define energy transfer for currents
   !q_4(1)= w +0.5d0*(e_gs-e_bg)+xmn-(p1_4(1)+p2_4(1))*0.5d0+20.0d0
   !if(q_4(1).lt.0.0d0) then
   !  r_now=0.0d0
   !  return
   !endif

   !......define constants and ff
   q2=q_4(1)**2-qval**2
   gep=1.0d0/(1.0d0-q2/lsq)**2 
   cv3=fstar/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
   ca5=0.0d0
   rho=xpf**3/(1.5d0*pi**2)

!.......currents
   call current_init(p1_4,p2_4,pp1_4,pp2_4,q_4,w,k1_4,k2_4,1,total_isospin)      
   call define_spinors()
   call det_Jpi(gep)
   call det_JpiJpi(r_cc_pi,r_cl_pi,r_ll_pi,r_t_pi,r_tp_pi)
   call det_JaJb_JcJd(e_gs,e_bg,cv3,ca5,np_del,pdel,pot_del)
   call det_JaJc_dir(r_cc_del,r_cl_del,r_ll_del,r_t_del,r_tp_del)
   call det_JpiJaJb(r_cc_int,r_cl_int,r_ll_int,r_t_int,r_tp_int)

   dir(1)=r_cc_pi+2.0d0*(r_cc_del+r_cc_int)
   dir(2)=r_cl_pi+2.0d0*(r_cl_del+r_cl_int)
   dir(3)=r_ll_pi+2.0d0*(r_ll_del+r_ll_int)
   dir(4)=r_t_pi+2.0d0*(r_t_del+r_t_int)
   dir(5)=r_tp_pi+2.0d0*(r_tp_del+r_tp_int)
   
   call current_init(p1_4,p2_4,pp2_4,pp1_4,q_4,w,k1e_4,k2e_4,2,total_isospin)
   call det_JaJb_JcJd(e_gs,e_bg,cv3,ca5,np_del,pdel,pot_del)
   call det_JaJc_exc(r_cc_del,r_cl_del,r_ll_del,r_t_del,r_tp_del)
   call det_JpiJaJb_exc(r_cc_int,r_cl_int,r_ll_int,r_t_int,r_tp_int)

   exc(1)=2.0d0*(r_cc_del+r_cc_int)
   exc(2)=2.0d0*(r_cl_del+r_cl_int)
   exc(3)=2.0d0*(r_ll_del+r_ll_int)
   exc(4)=2.0d0*(r_t_del+r_t_int)
   exc(5)=2.0d0*(r_tp_del+r_tp_int)
   
   

      r_now(:) =nk*p1**2*p2**2/(2.0d0*pi)**8*(dir(:)-exc(:))* &
   &      lorentz_jac/rho*dble(xA)/2.0d0/2.0d0 ! /2.0d0 for the electromagnetic piece

   return
end subroutine   




subroutine g_eval(pj1,pj2,dp12,pmax,norm,g)
    implicit none
    real*8, parameter :: pi=acos(-1.0d0)
    real*8 :: pj1,pj2,dp12,g,pmax,norm
    g=(4.0d0*pi)**2*pj1**2*pj2**2*dp12
    g=g/norm
    
    return
    end subroutine


end module        
