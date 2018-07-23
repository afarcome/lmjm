      subroutine lk_comp_joint_surv(bev,lnu,xiv,phi,psiv,k,piv,Pi,si2,
     c typ,n,M,npX,XM,npW,W,YM,PM,btv,Pp1,dev,tv,lk)

      integer k,n,M,npX,npW,i,mi,m1,dev(n),PM(n,M)
      double precision bev(npX),lnu,xiv(k),phi,psiv(npW),piv(k),Pi(k)
      double precision si2,XM(n,M,npX),YM(n,M),W(n,npW),btv(M)
      double precision Pp1(k,M,n),tv(n),btvnu(M),dbtvnu(M-1)
      double precision l2,l3,tmp00(k),Wp(n),tmp0(k),muv(k),etmp0(k),lk
      double precision nu,lsi2pi
      character (len=4) typ
      parameter (pig=3.1415926535897932)      

c to compute the second and third components of the expect complete log-likelihood

c preliminaries
      nu = exp(lnu)
      btvnu = btv**nu
	  dbtvnu = btvnu(2:M)-btvnu(1:(M-1))
	  if (typ=="norm") then
        lsi2pi = log(si2)+log(2*pig)
	  end if
c compute functions
	  l2 = 0
	  l3 = 0
	  tmp00 = xiv*phi
	  Wp = matmul(W,psiv)
	  do i=1,n
      	tmp0 = tmp00+Wp(i)
		etmp0 = exp(tmp0)
		mi = count(btv.le.tv(i))
		do m1 = 1,mi
          if (PM(i,m1)==1) then
            muv = xiv+sum(XM(i,m1,:)*bev)
            if (typ=="norm") then
              l2 = l2-0.5*sum((lsi2pi+(YM(i,m1)-muv)**2/si2)*
     c             Pp1(:,m1,i))
            end if
            if (typ=="bino") then
              muv = exp(muv)/(1+exp(muv))
              l2 = l2+sum((YM(i,m1)*log(muv)+(1-YM(i,m1))*log(1-muv))*
     c             Pp1(:,m1,i))
            end if
          end if
		end do
		if(mi>1) then
		  do m1=2,mi
		    l3 = l3-dbtvnu(m1-1)*sum(etmp0*Pp1(:,m1-1,i))
		  end do
		endif
        if (dev(i)==1) then
          l3 = l3+sum((lnu+(nu-1)*log(tv(i))+tmp0-
     c                 etmp0*(tv(i)**nu-btvnu(mi)))*Pp1(:,mi,i))
        else
          l3 = l3-sum((etmp0*(tv(i)**nu-btvnu(mi)))*Pp1(:,mi,i))
        end if
      end do
      lk = l2+l3

      end
