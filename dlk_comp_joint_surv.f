      subroutine dlk_comp_joint_surv(bev,lnu,xiv,phi,psiv,k,piv,Pi,si2,
     c typ,n,M,npX,XM,npW,W,YM,PM,btv,Pp1,dev,tv,dlk)

      integer k,n,M,npX,npW,i,mi,m1,dev(n),PM(n,M)
      double precision bev(npX),lnu,xiv(k),phi,psiv(npW),piv(k),Pi(k)
      double precision si2,XM(n,M,npX),YM(n,M),W(n,npW),btv(M)
      double precision Pp1(k,M,n),tv(n),btvnu(M),dbtvnu(M-1)
      double precision l2,l3,tmp00(k),Wp(n),tmp0(k),muv(k),etmp0(k),lk
      double precision nu,dl2bev(npX),dl3nu,dl2xiv(k),dl3xiv(k),dl3phi
      double precision dl3psiv(npW),tmp1(k),tmp(k),dlk(npX+1+k+1+npW)
      character (len=4) typ
      
c preliminaries
      nu = exp(lnu)
      btvnu = btv**nu
	  dbtvnu = btvnu(2:M)-btvnu(1:(M-1))
c compute functions
      dl2bev = 0
      dl3nu = 0
      dl2xiv = 0
      dl3xiv = 0
      dl3phi = 0
      dl3psiv = 0
      do i=1,n
		tmp0 = exp(xiv*phi+sum(W(i,:)*psiv))
		mi = count(btv.le.tv(i))
		do m1=1,mi
		  if(PM(i,m1)==1) then
			muv = xiv+sum(XM(i,m1,:)*bev)
			if (typ=="norm") then
              dl2bev = dl2bev+XM(i,m1,:)*sum(((YM(i,m1)-muv)/si2)*
     c                 Pp1(:,m1,i))
              dl2xiv = dl2xiv+(YM(i,m1)-muv)/si2*Pp1(:,m1,i)
			end if
			if (typ =="bino") then
			  muv = exp(muv)/(1+exp(muv))
			  tmp1 = (YM(i,m1)-muv)*Pp1(:,m1,i)
			  dl2bev = dl2bev+XM(i,m1,:)*sum(tmp1)
			  dl2xiv = dl2xiv+tmp1
			end if
		  end if
		end do
		if (mi>1) then
		  do m1=2,mi
            if (btv(m1-1)>0) then
		      tmp = log(btv(m1-1))
		    else
		      tmp = 0
		    end if
		    dl3nu = dl3nu-sum((tmp0*(log(btv(m1))*btvnu(m1)-
     c   		    tmp*btvnu(m1-1)))*Pp1(:,m1-1,i))
		    dl3xiv = dl3xiv-(tmp0*dbtvnu(m1-1))*phi*Pp1(:,m1-1,i)
		    dl3phi = dl3phi-sum((tmp0*dbtvnu(m1-1)*xiv)*Pp1(:,m1-1,i))
		    dl3psiv = dl3psiv-W(i,:)*sum((tmp0*dbtvnu(m1-1))*
     c                Pp1(:,m1-1,i))
		  end do
		end if
		if (btv(mi)>0) then
		  tmp = log(btv(mi))
		else
		  tmp = 0
		end if
		if (dev(i)==1) then
          dl3nu = dl3nu+sum(((1/nu+log(tv(i)))-tmp0*(log(tv(i))*
     c            tv(i)**nu-tmp*btvnu(mi)))*Pp1(:,mi,i))
		  dl3xiv = dl3xiv+(phi-tmp0*(tv(i)**nu-btvnu(mi))*phi)*
     c             Pp1(:,mi,i)
		  dl3phi = dl3phi+sum((xiv-tmp0*(tv(i)**nu-btvnu(mi))*
     c             xiv)*Pp1(:,mi,i))
		  dl3psiv = dl3psiv+W(i,:)*sum((1-tmp0*(tv(i)**nu-
     c              btvnu(mi)))*Pp1(:,mi,i))
		else
         dl3nu = dl3nu-sum((tmp0*(log(tv(i))*tv(i)**nu-tmp*btvnu(mi)))*
     c            Pp1(:,mi,i))
		  dl3xiv = dl3xiv-(tmp0*(tv(i)**nu-btvnu(mi))*phi)*Pp1(:,mi,i)
		  dl3phi = dl3phi-sum((tmp0*(tv(i)**nu-btvnu(mi))*xiv)*
     c             Pp1(:,mi,i))
		  dl3psiv = dl3psiv-W(i,:)*sum((tmp0*(tv(i)**nu-btvnu(mi)))*
     c              Pp1(:,mi,i))
		end if
      end do
c output
      dlk(1:npX) = dl2bev
      dlk(npX+1) = dl3nu*nu
      dlk((npX+2):(k+npX+1)) = dl2xiv+dl3xiv
      dlk(k+npX+2) = dl3phi
      dlk((k+npX+3):(k+npX+2+npW)) = dl3psiv
     
      end