program test_iso

   use dirac_matrices
   implicit none

   integer*4 :: i,j,k,l
   real*8 :: test_p4(4),w 
   complex*16 :: p(2), n(2), result1, iso(2,2), testing
   complex*16:: result(4)
   complex*16 :: czero = (0.0d0,0.0d0)
   complex*16 :: cone  = (1.0d0,0.0d0)
   complex*16 :: ci    = (0.0d0,1.0d0)

   w = 0.0d0
   test_p4 = (/0.0d0, 0.0d0, 0.0d0, 0.0d0/)

   call dirac_matrices_in(938.0d0,938.0d0,100.0d0)
   call current_init(w,test_p4,test_p4,test_p4,test_p4,test_p4,test_p4,test_p4,2)
   call define_spinors()

   p(1) = cone
   p(2) = czero

   n(1) = czero
   n(2) = cone

   iso(1,:) = p
   iso(2,:) = n

   !JAJA
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaADag(iso(i,:),iso(j,:),p,p)*IDeltaA(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaADag(iso(i,:),iso(j,:),p,n)*IDeltaA(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaADag(iso(i,:),iso(j,:),n,p)*IDeltaA(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaADag(iso(i,:),iso(j,:),n,n)*IDeltaA(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JA*JA--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JAJB
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaADag(iso(i,:),iso(j,:),p,p)*IDeltaB(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaADag(iso(i,:),iso(j,:),p,n)*IDeltaB(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaADag(iso(i,:),iso(j,:),n,p)*IDeltaB(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaADag(iso(i,:),iso(j,:),n,n)*IDeltaB(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JA*JB--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JAJC
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaADag(iso(i,:),iso(j,:),p,p)*IDeltaC(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaADag(iso(i,:),iso(j,:),p,n)*IDeltaC(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaADag(iso(i,:),iso(j,:),n,p)*IDeltaC(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaADag(iso(i,:),iso(j,:),n,n)*IDeltaC(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JA*JC--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JAJD
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaADag(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaADag(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaADag(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaADag(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JA*JD--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))


   !JBJB
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaBDag(iso(i,:),iso(j,:),p,p)*IDeltaB(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaBDag(iso(i,:),iso(j,:),p,n)*IDeltaB(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaBDag(iso(i,:),iso(j,:),n,p)*IDeltaB(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaBDag(iso(i,:),iso(j,:),n,n)*IDeltaB(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JB*JB--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))


   !JBJC
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaBDag(iso(i,:),iso(j,:),p,p)*IDeltaC(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaBDag(iso(i,:),iso(j,:),p,n)*IDeltaC(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaBDag(iso(i,:),iso(j,:),n,p)*IDeltaC(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaBDag(iso(i,:),iso(j,:),n,n)*IDeltaC(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JB*JC--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JBJD
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaBDag(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaBDag(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaBDag(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaBDag(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JB*JD--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JCJC
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaCDag(iso(i,:),iso(j,:),p,p)*IDeltaC(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaCDag(iso(i,:),iso(j,:),p,n)*IDeltaC(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaCDag(iso(i,:),iso(j,:),n,p)*IDeltaC(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaCDag(iso(i,:),iso(j,:),n,n)*IDeltaC(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JC*JC--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JCJD
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaCDag(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaCDag(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaCDag(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaCDag(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JC*JD--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JDJD
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaDDag(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) + IDeltaDDag(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) + IDeltaDDag(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) + IDeltaDDag(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'JD*JD--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))
   
   !! All of the exchange matrix elements are below
   !JAJAexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaADag(iso(i,:),iso(j,:),p,p)*IDeltaA(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaADag(iso(i,:),iso(j,:),p,n)*IDeltaA(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaADag(iso(i,:),iso(j,:),n,p)*IDeltaA(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaADag(iso(i,:),iso(j,:),n,n)*IDeltaA(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JA*JAexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JAJBexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaADag(iso(i,:),iso(j,:),p,p)*IDeltaB(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaADag(iso(i,:),iso(j,:),p,n)*IDeltaB(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaADag(iso(i,:),iso(j,:),n,p)*IDeltaB(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaADag(iso(i,:),iso(j,:),n,n)*IDeltaB(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JA*JBexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JAJCexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaADag(iso(i,:),iso(j,:),p,p)*IDeltaC(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaADag(iso(i,:),iso(j,:),p,n)*IDeltaC(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaADag(iso(i,:),iso(j,:),n,p)*IDeltaC(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaADag(iso(i,:),iso(j,:),n,n)*IDeltaC(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JA*JCexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JAJDexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaADag(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaADag(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaADag(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaADag(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JA*JDexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))


   !JBJBexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaBDag(iso(i,:),iso(j,:),p,p)*IDeltaB(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaBDag(iso(i,:),iso(j,:),p,n)*IDeltaB(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaBDag(iso(i,:),iso(j,:),n,p)*IDeltaB(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaBDag(iso(i,:),iso(j,:),n,n)*IDeltaB(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JB*JBexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))


   !JBJCexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaBDag(iso(i,:),iso(j,:),p,p)*IDeltaC(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaBDag(iso(i,:),iso(j,:),p,n)*IDeltaC(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaBDag(iso(i,:),iso(j,:),n,p)*IDeltaC(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaBDag(iso(i,:),iso(j,:),n,n)*IDeltaC(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JB*JCexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JBJDexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaBDag(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaBDag(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaBDag(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaBDag(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JB*JDexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JCJCexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaCDag(iso(i,:),iso(j,:),p,p)*IDeltaC(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaCDag(iso(i,:),iso(j,:),p,n)*IDeltaC(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaCDag(iso(i,:),iso(j,:),n,p)*IDeltaC(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaCDag(iso(i,:),iso(j,:),n,n)*IDeltaC(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JC*JCexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JCJDexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaCDag(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaCDag(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaCDag(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaCDag(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JC*JDexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !JDJDexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) + IDeltaDDag(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) + IDeltaDDag(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) + IDeltaDDag(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) + IDeltaDDag(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'JD*JDexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !!!!!!!!!!!!! Pion Delta interference below
   !PiJA
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*IDeltaA(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*IDeltaA(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*IDeltaA(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*IDeltaA(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'Pi*JA--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !PiJB
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*IDeltaB(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*IDeltaB(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*IDeltaB(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*IDeltaB(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'Pi*JB--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !PiJC
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*IDeltaC(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*IDeltaC(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*IDeltaC(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*IDeltaC(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'Pi*JC--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !PiJD
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'Pi*JD--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))


   !!!!!!!!!!!!! Pion Delta exchange interference below
   !PiJAexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*IDeltaA(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*IDeltaA(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*IDeltaA(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*IDeltaA(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'Pi*JAexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !PiJBexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*IDeltaB(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*IDeltaB(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*IDeltaB(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*IDeltaB(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'Pi*JBexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !PiJCexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*IDeltaC(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*IDeltaC(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*IDeltaC(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*IDeltaC(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'Pi*JCexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !PiJDexc
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*IDeltaD(p,p,iso(j,:),iso(i,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*IDeltaD(p,n,iso(j,:),iso(i,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*IDeltaD(n,p,iso(j,:),iso(i,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*IDeltaD(n,n,iso(j,:),iso(i,:))
      enddo
   enddo
   write(6,*)'Pi*JDexc--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))

   !PiPi
   result=czero
   do i=1,2
      do j=1,2
         result(1) = result(1) - Iv(iso(i,:),iso(j,:),p,p)*Iv(p,p,iso(i,:),iso(j,:))
         result(2) = result(2) - Iv(iso(i,:),iso(j,:),p,n)*Iv(p,n,iso(i,:),iso(j,:))
         result(3) = result(3) - Iv(iso(i,:),iso(j,:),n,p)*Iv(n,p,iso(i,:),iso(j,:))
         result(4) = result(4) - Iv(iso(i,:),iso(j,:),n,n)*Iv(n,n,iso(i,:),iso(j,:))
      enddo
   enddo
   write(6,*)'Pi*Pi--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   write(6,*)'sum = ', REAL(sum(result(:)))



   
   
end program test_iso




