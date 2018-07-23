      subroutine lk_rec_surv(k,piv,Pi,npX,bev,nu,xiv,phi,npW,psiv,si2,
     c typ,n,XM,W,YM,M,PM,btv,dev,tv,lk,Pp1,Pp2)

      integer k,npX,npW,n,M,mi,m1,i,dev(n),PM(n,M)
      double precision piv(k),Pi(k,k),bev(npX),nu,xiv(k),phi,psiv(npW)
      double precision si2,XM(n,M,npX),W(n,npW),YM(n,M),btv(M)
      double precision tv(n),lk,Pp1(k,M,n),Pp2(k,k)
      double precision btvnu(M),dbtvnu(M-1),tmp0(k),Frec(k,M),muv(k)
      double precision frecv(k),tmp(k),Grec(k,M),grecv(k),cost
      double precision Tmpm(k,k),fi,Pc(k,M,n)
      character (len=4) typ
      parameter (pig=3.1415926535897932)      
	
c preliminaries
      btvnu = btv**nu
      dbtvnu = btvnu(2:M)-btvnu(1:(M-1))
      if (typ=="norm") then
        cost = sqrt(2*pig*si2)
      end if
c recursion for each individual
      Pc = 1
      lk = 0
      do i = 1,n
        mi = count(btv.le.tv(i))
		tmp0 = exp(xiv*phi+sum(W(i,:)*psiv))
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
            tmp = exp(-tmp0*dbtvnu(m1-1))
            frecv = matmul(transpose(Pi),frecv*tmp)
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
		tmp = exp(-tmp0*(tv(i)**nu-btvnu(mi)))
		if (dev(i)==1) then
		  tmp = tmp*nu*tv(i)**(nu-1)*tmp0
		end if
		frecv = frecv*tmp
		fi = sum(frecv)
c accumulate log-likelihood
		lk = lk+log(fi)

c BACKWARD RECURSION
c last time occasion
		grecv = exp(-tmp0*(tv(i)**nu-btvnu(mi)))
		if (dev(i)==1) then
		  grecv = grecv*nu*tv(i)**(nu-1)*tmp0
		end if
		Grec(:,mi) = grecv
c previous time occasions
		if(mi>1) then
		  do m1 = mi-1,1,-1
			tmp = exp(-tmp0*dbtvnu(m1))
			if (PM(i,m1+1)==1) then
				grecv = grecv*Pc(:,m1+1,i)
			end if
			grecv = tmp*matmul(Pi,grecv)
			Grec(:,m1) = grecv
          end do	
		end if
c Posterior probabilities for the first time occasion
		Pp1(:,1:mi,i) = Frec(:,1:mi)*Grec(:,1:mi)/fi
c Joint Posterior probabilities for two occasions
		if (mi>1) then
		  do m1 = 1,mi-1
            tmp = exp(-tmp0*dbtvnu(m1))
			Tmpm = spread(Frec(:,m1)*tmp,2,k)
			Tmpm = Tmpm*spread(Grec(:,m1+1)*Pc(:,m1+1,i),1,k)*Pi/fi
			Pp2 = Pp2+Tmpm
          end do
        end if
      end do

      end