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

   call dirac_matrices_in(938.0d0,938.0d0,100.0d0, 0.0d0, 0.0d0)
   call current_init(w,test_p4,test_p4,test_p4,test_p4,test_p4,test_p4,test_p4,2)
   call define_spinors()

   p(1) = cone
   p(2) = czero

   n(1) = czero
   n(2) = cone

   iso(1,:) = p
   iso(2,:) = n

   !JAdag
   result=czero
   !do i=1,2
    !  do j=1,2
         result(1) = result(1) + IDeltaADag(n,p,p,p)
         result(2) = result(2) + IDeltaADag(n,n,p,n)
         result(3) = result(3) + IDeltaADag(p,p,n,p)
         result(4) = result(4) + IDeltaADag(p,n,n,n)
    !  enddo
   !enddo
   write(6,*)'JAdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))


   !JBdag
   result=czero
   !do i=1,2
   !   do j=1,2
         result(1) = result(1) + IDeltaBDag(n,p,p,p)
         result(2) = result(2) + IDeltaBDag(n,n,p,n)
         result(3) = result(3) + IDeltaBDag(p,p,n,p)
         result(4) = result(4) + IDeltaBDag(p,n,n,n)
   !   enddo
   !enddo
   write(6,*)'JBdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))


   !JCdag
   result=czero
   !do i=1,2
   !   do j=1,2
         result(1) = result(1) + IDeltaCDag(n,p,p,p)
         result(2) = result(2) + IDeltaCDag(n,n,p,n)
         result(3) = result(3) + IDeltaCDag(p,p,n,p)
         result(4) = result(4) + IDeltaCDag(p,n,n,n)
   !   enddo
   !enddo
   write(6,*)'JCdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))


   !JDdag
   result=czero
   !do i=1,2
   !   do j=1,2
         result(1) = result(1) + IDeltaDDag(n,p,p,p)
         result(2) = result(2) + IDeltaDDag(n,n,p,n)
         result(3) = result(3) + IDeltaDDag(p,p,n,p)
         result(4) = result(4) + IDeltaDDag(p,n,n,n)
   !   enddo
   !enddo
   write(6,*)'JDdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))
   
   !! All of the exchange matrix elements are below
   !JAexcdag
   result=czero
   !do i=1,2
    !  do j=1,2
         result(1) = result(1) + IDeltaADag(p,n,p,p)
         result(2) = result(2) + IDeltaADag(n,n,p,n)
         result(3) = result(3) + IDeltaADag(p,p,n,p)
         result(4) = result(4) + IDeltaADag(n,p,n,n)
   !   enddo
   !enddo
   write(6,*)'JAexcdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))


   !JBexcdag
   result=czero
   !do i=1,2
   !   do j=1,2
         result(1) = result(1) + IDeltaBDag(p,n,p,p)
         result(2) = result(2) + IDeltaBDag(n,n,p,n)
         result(3) = result(3) + IDeltaBDag(p,p,n,p)
         result(4) = result(4) + IDeltaBDag(n,p,n,n)
   !   enddo
   !enddo
   write(6,*)'JBexcdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))

   !JCexcdag
   result=czero
   !do i=1,2
   !   do j=1,2
         result(1) = result(1) + IDeltaCDag(p,n,p,p)
         result(2) = result(2) + IDeltaCDag(n,n,p,n)
         result(3) = result(3) + IDeltaCDag(p,p,n,p)
         result(4) = result(4) + IDeltaCDag(n,p,n,n)
   !   enddo
   !enddo
   write(6,*)'JCexcdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))

   !JDexcdag
   result=czero
   !do i=1,2
   !   do j=1,2
         result(1) = result(1) + IDeltaDDag(p,n,p,p)
         result(2) = result(2) + IDeltaDDag(n,n,p,n)
         result(3) = result(3) + IDeltaDDag(p,p,n,p)
         result(4) = result(4) + IDeltaDDag(n,p,n,n)
   !   enddo
  ! enddo
   write(6,*)'JDexcdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))

   !!!!!!!!!!!!! Pion 
   !Pidag
   result=czero
   !do i=1,2
   !   do j=1,2
         result(1) = result(1) - Ivminus(n,p,p,p)
         result(2) = result(2) - Ivminus(n,n,p,n)
         result(3) = result(3) - Ivminus(p,p,n,p)
         result(4) = result(4) - Ivminus(p,n,n,n)
   !   enddo
   !enddo
   write(6,*)'Pidag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))


   !Piexcdag
   result=czero
   !do i=1,2
   !   do j=1,2
         result(1) = result(1) - Ivminus(p,n,p,p)
         result(2) = result(2) - Ivminus(n,n,p,n)
         result(3) = result(3) - Ivminus(p,p,n,p)
         result(4) = result(4) - Ivminus(n,p,n,n)
   !   enddo
   !enddo
   write(6,*)'Piexcdag--------------------'
   write(6,*)'pp = ',REAL(result(1))
   write(6,*)'pn = ',REAL(result(2))
   write(6,*)'np = ',REAL(result(3))
   write(6,*)'nn = ',REAL(result(4))



   
   
end program test_iso




