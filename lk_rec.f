      subroutine lk_rec(k,piv,Pi,npX,bev,xiv,si2,
     c typ,n,XM,YM,M,PM,btv,tv,lk,Pp1,Pp2)

      integer k,npX,npW,n,M,mi,m1,i,PM(n,M)
      double precision piv(k),Pi(k,k),bev(npX),xiv(k)
      double precision si2,XM(n,M,npX),YM(n,M),btv(M)
      double precision tv(n),lk,Pp1(k,M,n),Pp2(k,k)
      double precision tmp0(k),Frec(k,M),muv(k)
      double precision frecv(k),tmp(k),Grec(k,M),grecv(k),cost
      double precision Tmpm(k,k),fi,Pc(k,M,n)
      character (len=4) typ
      parameter (pig=3.1415926535897932)

c preliminaries
      if (typ=="norm") then
        cost = sqrt(2*pig*si2)
      end if	
c recursion for each individual
      Pc = 1
      lk = 0
      do i = 1,n
        mi = count(btv.le.tv(i))
c FORWARD RECURSION	
c first time
		muv = xiv+sum(XM(i,1,:)*bev)
        if (typ=="norm") then
          Pc(:,1,i) = exp(-(YM(i,1)-muv)**2/(2*si2))/cost
        end if          
		if (typ=="bino") then
		  muv = exp(muv)/(1+exp(muv))
		  Pc(:,1,i) = YM(i,1)*muv+(1-YM(i,1))*(1-muv)
		end if
        frecv = piv*Pc(:,1,i)
		Frec(:,1) = frecv
c following time occasions
		if (mi>1) then
		  do m1 = 2,mi
            frecv = matmul(transpose(Pi),frecv)
			if (PM(i,m1)==1) then
			  muv = xiv+sum(XM(i,m1,:)*bev)
			  if(typ=="norm") then
                Pc(:,m1,i) = exp(-(YM(i,m1)-muv)**2/(2*si2))/cost
              end if
			  if (typ=="bino") then
			    muv = exp(muv)/(1+exp(muv))
			    Pc(:,m1,i) = YM(i,m1)*muv+(1-YM(i,m1))*(1-muv)
			  end if
              frecv = frecv*Pc(:,m1,i)
			end if
			Frec(:,m1) = frecv
          end do
		end if
		fi = sum(frecv)
c accumulate log-likelihood
		lk = lk+log(fi)

c BACKWARD RECURSION
c last time occasion
		grecv = 1
		Grec(:,mi) = grecv
c previous time occasions
        if(mi>1) then
          do m1 = mi-1,1,-1
			if (PM(i,m1+1)==1) then
			  grecv = grecv*Pc(:,m1+1,i)
			end if
			grecv = matmul(Pi,grecv)
			Grec(:,m1) = grecv
          end do	
		end if
c Posterior probabilities for the first time occasion
		Pp1(:,1:mi,i) = Frec(:,1:mi)*Grec(:,1:mi)/fi
c Joint Posterior probabilities for two occasions
		if (mi>1) then
		  do m1 = 1,mi-1
			Tmpm = spread(Frec(:,m1),2,k)
			Tmpm = Tmpm*spread(Grec(:,m1+1)*Pc(:,m1+1,i),1,k)*Pi/fi
			Pp2 = Pp2+Tmpm
          end do
        end if
      end do
      
      end