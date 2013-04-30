      program localmass
************************************************************************
*       Program for calculating local mass from correlators            *
*       Calculates Jackknived local masses                             *
************************************************************************
      implicit double precision(a-h,m,o-z)
      parameter(nt=40,nt2=nt/2,nb=100,nb1=nb-1,nconf=100)
      parameter(ntm1=nt-1)
      dimension corj(0:ntm1,nb)
      dimension corb(0:ntm1,0:nb)
      dimension rmas(0:nb),rcorr(0:nb)
      dimension a(0:ntm1)
      CHARACTER::filenam*50
      do j=0,nt-1,1
         corb(j,0)=0.d0
         do k=1,nb
            corj(j,k)=0.d0
            corb(j,k)=0.d0
         enddo
         a(j)=0.d0
      enddo
************************************************************************
c     Reading file and evaluating the blocks
      
      open(unit=16,file='CORRLT40',status='old')
      ntherm=0
      ntr=(nconf-ntherm)/nb
      do i=1,nb
         do k=1,ntr
            do j=0,nt-1
               read(16,*)zz,a(j)
               corj(j,i)=corj(j,i)+a(j)
            enddo
         enddo
         do j=0,nt-1,1
            corj(j,i)=corj(j,i)/(ntr)
         enddo
      enddo
      close(unit=16,status='keep')
      do j=0,nt-1
         corb(j,0)=0.d0
         do i=1,nb
            corb(j,0)=corb(j,0)+corj(j,i)
         enddo
      enddo
      do i=1,nb
         do j=0,nt-1
            corb(j,i)=(corb(j,0)-corj(j,i))/nb1
         enddo
      enddo
      do j=0,nt-1
         corb(j,0)=corb(j,0)/nb
      enddo
      do j=0,nt-1
         do i=0,nb
            rcorr(i)=corb(j,i)
         enddo
         call jack(rcorr,val,err)
         write(28,*)j,val,err
      enddo
      nt21=nt2-1
      do j=0,nt21
         do i=0,nb
            rat=corb(j+1,i)/corb(j,i)
            ij=j
            if (j.eq.nt21) then
               rms=(1.0+sqrt(1.0-rat**2))/rat
            else
               rms=root(rat,ij)
            endif
            rmas(i)=log(rms)
         enddo
         call jack(rmas,val,err)
         print *,float(j)+0.5,val,err
      enddo
      end
************************************************************************


      function root(xf,iloc)
************************************************************************
*       Finds roots of eqn. func=0                                     *
*       By Newton-Raphson's method                                     *
************************************************************************
      double precision root,rtsafe,xf,accu,df,dx,f,pr
      character con
      external func,rtsafe
      max=200
      accu=1d-4
      pr=1d0/xf
      do iso=1,max
         call func(pr,f,df,xf,iloc)
         dx=f/df
         pr=pr-dx
         if (abs(dx).lt.accu) then
            root=pr
            return
         endif
      enddo
      root=pr
      write(*,801)iloc,dx
      write(*,*)'(c)ontinue with present value?(t)ry rtsafe?(q)uit?'
c 800  read (*,*)con
      con='t'
c        con='t'
      if (con.eq.'c') then
         return
      elseif (con.eq.'t') then
         root=rtsafe(xf,iloc)
      elseif (con.eq.'q') then
         stop
      else
         write(*,*)'please answer c/t/q'
      endif
      return
 801  format("No convergence at distance",i2,"accuracy obtained",f7.4)
      end
************************************************************************


      subroutine func(x,fn,der,xg,ilcn)
************************************************************************
*       Evaluates the function and its derivative                      *
*       For Newton-Raphson iteration                                   *
************************************************************************
      parameter(nt=60,nt2=nt/2)
      double precision x,fn,der,xg,trm1
      icf=nt2-ilcn
      trm1=x**icf
      fn=trm1*(xg-1d0/x)+(xg-x)/trm1
      der=xg*icf*(trm1/x-1/(x*trm1))-(icf-1)*(trm1/(x**2)-1d0/trm1)
      return
      end
************************************************************************


      subroutine jack(estm,val,err)
************************************************************************
*       Given the block values, forms the Jackknife estimators         *
*           And calculates values and errors                           *
************************************************************************
      parameter(nconf=50)
      double precision estm(0:nconf)
      double precision val,err,x0,x
      val=0.d0
      err=0.d0
      x0=nconf*estm(0)
      nb1=nconf-1
      do ik=1,nconf
         x=x0-nb1*estm(ik)        !The Jackknife estimators
         val=val+x
         err=err+x*x
      enddo
      val=val/nconf
      err=(err/nconf-val*val)/nb1
      if (err.le.0.) then
         err=0.
      else
         err=sqrt(err)
      endif
      return
      end
************************************************************************


      function rtsafe(st,ilc)
************************************************************************
*       Uses a combination of Newton-Raphson and Bisection             *
************************************************************************
      double precision rtsafe,accu,st,xup,xlr,f,df,fnl,x1,x2,dx,dol,tmp
      double precision xdn
      external func
      max=200
      accu=1d-4
      xlr=1d0/st
      call func(xlr,f,df,st,ilc)
      if (f.eq.0d0) then
         rtsafe=xlr
         return
      endif
      fnl=f
      xup=xlr
      xdn=xlr
      do iull=1,200
      xup=xup+1d-1
      xdn=xdn-1d-1
      call func(xup,f,df,st,ilc)
      if((f*fnl).le.0d0) go to 900
      call func(xdn,f,df,st,ilc)
      if((f*fnl).le.0d0) go to 910
      enddo
      write(*,*)'rtsafe could not find the limits'
      return
 900  if (f.eq.0d0)then
         rtsafe=xup
         return
      endif
      x1=xlr
      x2=xup
      go to 920
 910  if (f.eq.0d0) then
         rtsafe=xdn
         return
      endif
      x1=xdn
      x2=xlr
      go to 920
 920  rtsafe=.5*(x1+x2)
      dol=abs(x2-x1)
      dx=dol
      call func(rtsafe,f,df,st,ilc)
      do irs=1,max
         if(((rtsafe-x1)*df-f)*((rtsafe-x2)*df-f).gt.0d0
     #         .or.abs(2.*f).gt.abs(dol*df)) then
            dol=dx
            dx=.5*(x2-x1)
            rtsafe=x1+dx
            if(x1.eq.rtsafe)return
         else
            dol=dx
            dx=f/df
            tmp=rtsafe
            rtsafe=rtsafe-dx
            if(tmp.eq.rtsafe)return
         endif
         if(abs(dx).lt.accu)return
         call func(rtsafe,f,df,st,ilc)
         if (f.lt.0d0) then
            x1=rtsafe
         else
            x2=rtsafe
         endif
      enddo
      write(*,*)'rtsafe cannot converge either'
      stop
      end
************************************************************************
************************************************************************
