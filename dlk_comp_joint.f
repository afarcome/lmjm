      subroutine dlk_comp_joint(bev,xiv,k,piv,Pi,si2,
     c typ,n,M,npX,XM,YM,PM,btv,Pp1,tv,dlk)

      integer k,n,M,npX,i,mi,m1,PM(n,M)
      double precision bev(npX),xiv(k),piv(k),Pi(k)
      double precision si2,XM(n,M,npX),YM(n,M),btv(M)
      double precision Pp1(k,M,n),tv(n)
      double precision l2,muv(k),lk
      double precision dl2bev(npX),dl3nu,dl2xiv(k),dl3xiv(k),dl3phi
      double precision tmp1(k),dlk(npX+k)
      character (len=4) typ
      
c compute functions
      dl2bev = 0
      dl3nu = 0
      dl2xiv = 0
      dl3xiv = 0
      do i=1,n
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
      end do
c output
      dlk(1:npX) = dl2bev
      dlk((npX+1):(npX+k)) = dl2xiv
     
      end