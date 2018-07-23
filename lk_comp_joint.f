      subroutine lk_comp_joint(bev,xiv,k,piv,Pi,si2,
     c typ,n,M,npX,XM,YM,PM,btv,Pp1,tv,lk)

      integer k,n,M,npX,npW,i,mi,m1,PM(n,M)
      double precision bev(npX),xiv(k),piv(k),Pi(k)
      double precision si2,XM(n,M,npX),YM(n,M),btv(M)
      double precision Pp1(k,M,n),tv(n)
      double precision Wp(n),muv(k),lk,lsi2pi
      character (len=4) typ
      parameter (pig=3.1415926535897932)

c to compute the second and third components of the expect complete log-likelihood

c preliminaries
	  if (typ=="norm") then
        lsi2pi = log(si2)+log(2*pig)
	  end if
c compute functions
	  do i=1,n
		mi = count(btv.le.tv(i))
		do m1 = 1,mi
          if (PM(i,m1)==1) then
            muv = xiv+sum(XM(i,m1,:)*bev)
            if (typ=="norm") then
              lk = lk-0.5*sum((lsi2pi+(YM(i,m1)-muv)**2/si2)*
     c             Pp1(:,m1,i))
            end if
            if (typ=="bino") then
              muv = exp(muv)/(1+exp(muv))
              lk = lk+sum((YM(i,m1)*log(muv)+(1-YM(i,m1))*log(1-muv))*
     c             Pp1(:,m1,i))
            end if
          end if
		end do
      end do

      end
