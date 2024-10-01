module isospin
    implicit none

	contains

	function Iv(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: Iv

	    Iv = ci*(me(1,it1,itp1)*me(2,it2,itp2) - me(2,it1,itp1)*me(1,it2,itp2))

	    return
	end function Iv


	function IDeltaA(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: IDeltaA, c

	    c = me(3,it2,itp2)*iden(it1,itp1)

	    IDeltaA = (2.*c/3.) - (Iv(it1,it2,itp1,itp2)/3.)

	    return
	end function IDeltaA

	function IDeltaADag(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: IDeltaADag, c

	    c = me(3,it2,itp2)*iden(it1,itp1)

	    IDeltaADag = (2.*c/3.) + (Iv(it1,it2,itp1,itp2)/3.)


	    return
	end function IDeltaADag

	function IDeltaB(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: IDeltaB, c

	    c = me(3,it2,itp2)*iden(it1,itp1) 

	    IDeltaB = (2.*c/3.) + (Iv(it1,it2,itp1,itp2)/3.) 

	    return
	end function IDeltaB

	function IDeltaBDag(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: IDeltaBDag, c

	    c = me(3,it2,itp2)*iden(it1,itp1)

	    IDeltaBDag = (2.*c/3.) - (Iv(it1,it2,itp1,itp2)/3.) 

	    return
	end function IDeltaBDag

	function IDeltaC(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: IDeltaC, c

	    c = me(3,it1,itp1)*iden(it2,itp2) 

	    IDeltaC = (2.*c/3.) + (Iv(it1,it2,itp1,itp2)/3.)

	    return
	end function IDeltaC

	function IDeltaCDag(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: IDeltaCDag, c

	    c = me(3,it1,itp1)*iden(it2,itp2)  

	    IDeltaCDag = (2.*c/3.) - (Iv(it1,it2,itp1,itp2)/3.)

	    return
	end function IDeltaCDag

	function IDeltaD(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: IDeltaD, c

	    c = me(3,it1,itp1)*iden(it2,itp2)  

	    IDeltaD = (2.*c/3.) - (Iv(it1,it2,itp1,itp2)/3.) 

	    return
	end function IDeltaD

	function IDeltaDDag(it1,it2,itp1,itp2)
	    implicit none
	    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
	    complex*16 :: IDeltaDDag, c

	    c = me(3,it1,itp1)*iden(it2,itp2)  

	    IDeltaDDag = (2.*c/3.) + (Iv(it1,it2,itp1,itp2)/3.) 

	    return
	end function IDeltaDDag

	function me(i,it,itp)
	    implicit none
	    integer*4 :: i
	    complex*16 :: me, it(2),itp(2), matrix(2)

	    me = sum(itp(:)*matmul(sig(i,:,:),it))

	    return
	end function me


	function iden(it,itp)
	    implicit none 
	    complex*16:: iden, it(2),itp(2)

	    iden = sum(itp(:)*matmul(id,it))
	    return
	end function iden
end module
