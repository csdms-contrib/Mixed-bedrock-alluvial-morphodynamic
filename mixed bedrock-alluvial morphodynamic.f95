Program Mixed_bedrock_alluvial



!******************** Mixed bedrock_alluvial morphodynamic model ************************
!************************ Copyright (C) 2019  Sadegh Jafarinik **************************



!****************************************************************************************
!****************************************************************************************
!*********** This program is free software: you can redistribute it and/or modify********
!*********** it under the terms of the GNU General Public License as published by********
!*********** the Free Software Foundation, either version 3 of the License  *************
!*********** , or (at your option) any later version.                       *************
!*********** This program is distributed in the hope that it will be useful, ************
!*********** but WITHOUT ANY WARRANTY; without even the implied warranty of *************
!*********** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the **************
!*********** GNU General Public License for more details.                  **************
!*********** You should have received a copy of the GNU General Public License **********
!*********** along with this program.  If not, see <https://www.gnu.org/licenses/>. *****
!****************************************************************************************
!****************************************************************************************

!----------------define parameter------------------------
implicit none
integer :: N=21,  M=1000000, i, j,jj,v,k

Real*8 :: dt,Sbase, etabo, etaao, lmr, b, qw,qww, Gbf, lp, nk
real*8 :: R, g, qbf, Ho,xt,dx,fr2,func,so,d50,dg
real*8 ::  cumulative,qbft=0
real*8 :: psibar=0,psibarf=0,sigmaf

!---------------D90,ks,kc,la-------------------
real*8 :: fad=0,fau=0,d90i,alphr
real*8, dimension(21):: d90,laprime=0.01,cf=0.15,ks,rh,dlaprime=0,fr,kbrs,kbr,kbf,ksc,kt

!---------------Side wall correction parameters---------------------
real*8 :: alph=0.2,la,reb=0,ab0=0,ab1=0,ratio=0,cont=0,fw0=0,fw1=0,rew=0,rw=0,aw=0
Real*8, dimension(21):: cfb,cfre,cfwrew,fwc,rewc,sigma1=.002,cfbs=0
Real*8, dimension(21):: re,ab,a,pw,cfw,rou=1000,taub,tauw,taut,rhw,fre
Real*8, dimension(21):: rbs=1.0,cfs=0,rb=0,rbs1=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, dimension(9):: db,psib,ff
real*8, dimension(8):: d,psi,f
real*8, dimension(8):: qbfi
real*8,dimension(21,8):: taustarc=0,qstarb=0,taustar=0,qb=0,fa,fl=.166,pb, fs=.166
real*8, dimension(21):: x=0,etab=0,etaa=0,pc=0,fc=0,h=0,slope=0,qbt=0,detaa=0 ,lprime =0,dgf=0 ,u

!--------------Output parameters-------------------
real*8, dimension(21,102) :: etaaout=0, etabout=0, qbout=0, hout=0, dgout=0,pcout=0,slopeout=0




!---------open results files-----------------------
open(1,file='etaa')
open(2,file='etab')
open(4,file='dgf')
open(5,file='water depth')
open(7,file='slope')
open(8,file='pc')
open(12,file='qbt')

!-----------------Define Parameters----------------


so= .002        ! Initial bed surface slope
alphr=8.1       ! Constant in friction coefficient relation
dt= 0.3155      ! sec (depending on your dx and other parameter this dt can be increased)
kbrs=0.0001     ! roughness of the bedrock surface (m)
sbase= 0.000    ! bedrock slope 
etabo= -0.12    ! downstream most node of the bedrock relative to the fixed water downstream water level (m)
etaao= -.07     ! initial value for downstream node of the bed surface relative to the downstream water level(m)
sigmaf=0.004    ! Standard deviation of the bedforms in fully alluvial reach (representative of the average bedform height) (m)
lmr=1.5*sigmaf  ! Minimum thickness of alluvial cover (the coefficient can be calibrated for different bedform heights)
dx=.2905        ! m  
b=.19           ! Channel width (m)
qw=.01          ! Flow discharge (m^3/S)
qww=qw/b        ! Flow discharge per unit channel width (m^2/S)
gbf=400         ! sediment supply rate (gr/min of sand)
qbf= gbf/(60*1000*1650*b)       !Sediment supply per unit channel width (m^2/S)
lp=0.4          ! sediment porosity (lambda)
r= 1.65         ! submerged specific gravity
g=9.81          ! Gravity acceleration (m/s^2)
nk=2.0          

!-------------Grain size distribution input------------------

!------Diameters-------       
db(1)=.297     !mm        
db(2)=.422 !mm
db(3)=.599    !mm
db(4)=.853   !mm
db(5)=1.0		   !mm
db(6)=1.4      !mm
db(7)=1.680      !mm
db(8)=2.38      !mm
db(9)=3.36      !mm
!-------Fractions-------
ff(1)=0  
ff(2)=0.072
ff(3)=0.246
ff(4)=0.452
ff(5)=0.738
ff(6)=.870
ff(7)=.928
ff(8)=.966
ff(9)=1.0				

                        
!----------------------sediment distribution-------------------
qbfi(1)=qbf*ff(2)
qbfi(2)=qbf*(ff(3)-ff(2))
qbfi(3)=qbf*(ff(4)-ff(3))
qbfi(4)=qbf*(ff(5)-ff(4))
qbfi(5)=qbf*(ff(6)-ff(5))
qbfi(6)=qbf*(ff(7)-ff(6))
qbfi(7)=qbf*(ff(8)-ff(7))
qbfi(8)=qbf*(ff(9)-ff(8))
do v=1,8
qbft=qbft+qbfi(v)
end do


!-----------------Geometric mean diameter of sediment input---------------
psibar=0
    
do i=1,9
  psib(i)=log(db(i))/log(2.0)
end do

do i=1,8
  psi(i)=.5*(psib(i)+psib(i+1))
  d(i)=(db(i)*db(i+1))**.5 
  f(i)=ff(i+1)-ff(i)
end do

do i=1,8
  psibar=psi(i)*f(i)+psibar
end do
dg=2**psibar
dg=dg/1000
d=d/1000

!------------D90 of the sediment input----------------------------

do v=1,8
  fau=f(v)+fau
  if (fau>=.90)then
    fad=fau-f(v)
    d90i=exp((log(db(v+1))-log(db(v)))*(.9-fad)/(fau-fad)+log(db(v)))/1000
    exit 
  end if
end do

!-----------Initial active layer thickness-----------------------

La=2*d90i

!------------Initial sediment fraction in avtive layer and substrate----------
do i=1,n
  do v=1,8
    fs(i,v)=f(v)
    fa(i,v)=f(v)
  end do
end do

              
                                             
!------------compute initial condition--------------------

do i=1,n
  x(i)= (i-1)                              ! node numbers starts from upstream
  etab(i)=etabo+sbase*(x(n)-x(i))*dx       ! bedrock surface elevation relative to the water level at downstream
end do

!---------- Initial Pc calculation--------------------
etaa(n) = etaao
 
do i=1,n
  if (i<n) then
    etaa(i)=etaa(n)+so*(x(n)-x(i))*dx
  end if
  if (etaa(i)<etab(i)) then
    etaa(i)=etab(i)
  end if
  fc(i)=(etaa(i)-etab(i))/lmr
  if (fc(i)>1.0) then
    pc(i)=1.0
  else if (fc(i)<0.0) then
    pc(i)=0.05
  else
    pc(i)=0.05+0.95*fc(i)
  end if
end do

!--------- Initial water depth, Cf and total shields number everywhere---------------

!---------the downstream most node---------     
h(n)= -etaa(n)
u(n) =qww/h(n)
rh(n)=b*h(n)/(b+2*h(n))
re(n)=rh(n)*u(n)/0.000001
ks(n)=2*d90i
fr(n)=u(n)/(9.81*h(n))**.5
cf(n)=(alphr**(-2))*(h(n)/ks(n))**(-1/3)


!--------the rest of the nodes-----------
do i=1,n
  if (i<n) then
    slope(i)=(etaa(i)-etaa(i+1))/dx
  else
    slope(i)=(etaa(i-1)-etaa(i))/dx
  end if  
end do
do i=1,n
  if (i<n) then
    u(n-i+1)=qww/h(n-i+1)
    fr2=(u(n-i+1)**2)/(g*h(n-i+1))
    func=(slope(n-i)-cf(n-i)*fr2)/(1-fr2)
    h(n-i)=h(n-i+1)-func*dx
    u(n-i)=qww/h(n-i)
    rh(n-i)=b*h(n-i)/(b+2*h(n-i))  
    ks(n-i)=2*d90i
    cf(n-i)=(alphr**(-2))*(h(n-i)/ks(n-i))**(-1/3)
    re(n-i)=rh(n-i)*u(n-i)/0.000001
    fr(n-i)=u(n-i)/(9.81*h(n-i))**.5
  
  end if
end do


       
!----------------------Initial Bedload calculation (Ashida Michiue)-----------------

do i=1,n
  do v=1,8

    if ((d(v)/dg)>0.4) then
      taustarc(i,v)=.05*(log(19.0)/log(19.0*(d(v)/dg)))**2
    else
      taustarc(i,v)=.05*.843*(d(v)/dg)**(-1)
    end if

    taustar(i,v)=cf(i)*(u(i)**2)/(r*g*d(v))

    if (taustar(i,v)>taustarc(i,v)) then
      qstarb(i,v)= 17*(taustar(i,v)-taustarc(i,v))*(taustar(i,v)**.5-taustarc(i,v)**.5)
    else 
      qstarb(i,v)=0
    end if

    qb(i,v)=fa(i,v)*pc(i)*qstarb(i,v)*(r*g*d(v))**.5*d(v)
    qbt(i)=qbt(i)+qb(i,v)
  end do
  do v=1,8
    if(qbt(i)==0)then
      pb(i,v)=0
    else  
    pb(i,v)=qb(i,v)/qbt(i)
    end if
  end do
end do


    
    

!------------------Printout Initial Conditoin--------------------
etaaout(:,2)=etaa
etabout(:,1)=x
etabout(:,2)=etab
slopeout(:,1)=x
slopeout(:,2)=slope
dgout(:,1)=x
dgout(:,2)=dg
qbout(:,1)=x
qbout(:,2)=qbt
hout(:,1)=x
hout(:,2)=h
pcout(:,1)=x
pcout(:,2)=pc

               
!*********************Start the calculations for time steps*******************
!*****************************************************************************
!*****************************************************************************
                

do jj=1,m
               
  do i=1,n      
!---------------------Compute new bed surface------------
     
     if (i==1) then
         etaa(i)=etaa(i)-((qbt(i)-qbf)/((1.0-lp)*dx*pc(i)))*dt

         detaa(i)=-(qbt(i)-qbf)/((1.0-lp)*dx*pc(i))
     else

         etaa(i)=etaa(i)-((qbt(i)-qbt(i-1))/((1.0-lp)*dx*pc(i)))*dt

         detaa(i)=-(qbt(i)-qbt(i-1))/((1.0-lp)*dx*pc(i))
     end if
!---------------------Compute new Pc-------------------

     fc(i)=(etaa(i)-etab(i))/lmr

     if (fc(i)>=1.0) then

          pc(i)=1.0
     else if (fc(i)<=0.0) then
          pc(i)=0.05
     else
          pc(i)=(0.05+0.95*fc(i))
     end if
      
  end do

!-------------------Compute fraction of grains at the interface of substracte and active layer (fl)-------------

  do i=1,n
     do v=1,8
        if ((detaa(i)-dlaprime(i))>0) then
            fl(i,v)=alph*fa(i,v)+(1-alph)*pb(i,v)

      else if ( (etaa(i)-laprime(i))<etab(i)) then
        fl(i,v)=0

      else
        fl(i,v)=f(v)

      end if
    end do
  end do

!-----------------Compute the fraction of grains in the active layer (fa)------------------
  do i=1,n
    do v=1,8
      if (i==1) then   

        fa(i,v)=((-qb(i,v)+qbfi(v))/(dx*(1.0-lp))-fl(i,v)*pc(i)*detaa(i)-fa(i,v)*dlaprime(i))*&
        &(dt/laprime(i))+fa(i,v)     

      else 

        fa(i,v)=((-qb(i,v)+qb(i-1,v))/(dx*(1.0-lp))-fl(i,v)*pc(i)*detaa(i)-fa(i,v)*dlaprime(i))*&
        &(dt/laprime(i))+fa(i,v)  

      end if
    
    
    end do

!--------------Normalize the fraction--------------------

      cumulative = 0
      do v=1,8

		 cumulative = cumulative + fa(i,v)
      end do

      do v= 1,8

	 	fa(i,v) = fa(i,v)/cumulative
      end do
   
  end do

!------------Geometric mean diamter of grains in the active layer------------
                          
  do i=1,n

    psibarf=0               
    do v=1,8

      psibarf=psi(v)*fa(i,v)+psibarf
    end do

    dgf(i)=2**psibarf
    dgf(i)=dgf(i)/1000
  end do

!-------------D90 of the grains in the active layer------------------

  do i=1,n
    fau=0
    fad=0
    do v=1,9
      fau=fa(i,v)+fau
      if (fau>=.90)then
        fad=fau-fa(i,v)
        d90(i)=exp((log(db(v+1))-log(db(v)))*(.9-fad)/(fau-fad)+log(db(v)))/1000
        exit 
      end if
    end do
  end do

!---------------Compute La, laprime and Ks-------------------

  do i=1,n
     dlaprime(i)=laprime(i)
     if ((etaa(i)-la-etab(i))>lmr) then
       laprime(i)=la
     else
       laprime(i)=2*sigma1(i)*dgf(i)
     end if    
     dlaprime(i)=laprime(i)-dlaprime(i)

     Ks(i)=nk*d90(i)

  end do

   
!----------------Computation of flow---------------------------
    
  do i= 1,n
     if (i<n) then
        slope(i)=(etaa(i)-etaa(i+1))/dx
     else
        slope(i)=(etaa(i-1)-etaa(i))/dx
     end if
  end do

  h(n)=-etaa(n)
  
  do i=1,n

    u(n-i+1)=qww/h(n-i+1)
    rh(n-i+1)=b*h(n-i+1)/(b+2*h(n-i+1))
    fr(n-i+1)=u(n-i+1)/(9.81*h(n-i+1))**.5
    sigma1(n-i+1)=(-25.968*fr(n-i+1)+23.505)          ! It should be calibrated for any specific case
    kt(n-i+1)=d90(n-i+1)*(.1739*exp(.3515*sigma1(n-i+1)))     ! It should be calibrated for any specific case
    cont=1

!------------------Sidewall correction and Einstein Decompistion----------------------

     do
       if (cont==1) then
         if (b>2*h(n-i+1)) then
            ab1=b*h(n-i+1)-(h(n-i+1)**2)
         else
            ab1=b*h(n-i+1)/2.0 
         end if 
       end if       
       rb(n-i+1)=ab1/b
       cfb(n-i+1)=(alphr**(-2))*(rb(n-i+1)/kt(n-i+1))**(-1.0/3)
       reb=rb(n-i+1)*(U(n-i+1))/.000001
       ratio=8*cf(n-i+1)/re(n-i+1)
       fw0=.301*(ratio**(1.0/5.0))
       rew=fw0/ratio
       do
          fw1=fw0
          fw0=(1/(0.86*log((4*rew)*(fw1**0.5))-.8))**(2.0)
          if ((abs(fw1-fw0)/fw0)<0.000000001)then 
             exit
          end if  
       end do
       fw0=(fw0+fw1)/2
       rew=fw0/ratio
       rw=rew*0.000001/(u(n-i+1))
       aw=rw*2*h(n-i+1)
       ab0=b*h(n-i+1)-aw
       cont=cont+1
       if((abs(ab1-ab0)/ab0)<.00000001) then
          exit
       end if
       ab1=ab0
       ab0=(ab1+ab0)/2.0
       rb(n-i+1)=ab0/b
       reb=rb(n-i+1)*(u(n-i+1)/0.000001)
       cfb(n-i+1)=(alphr**(-2))*(rb(n-i+1)/kt(n-i+1))**(-1.0/3)
       cf(n-i+1)=(cfb(n-i+1)*b+(fw0/8.0)*2.0*h(n-i+1))/(b+2.0*h(n-i+1))
       re(n-i+1)=rh(n-i+1)*u(n-i+1)/0.000001
    end do
    fr2=(u(n-i+1)**2.0)/(g*h(n-i+1))
    if (i<n) then  
       func=(slope(n-i)-cfb(n-i+1)*u(n-i+1)**2.0/(g*h(n-i+1)))/(1-fr2)    
       h(n-i)=h(n-i+1)-func*dx 
    end if
    ab(n-i+1)=ab1
  end do 

!--------- Calculation of friction coefficient due to skin friction--------------
  do i=1,n

    rb(i)=ab(i)/b
    ksc(i)=pc(i)*ks(i)+(1-pc(i))*kbrs(i)
    do 
       rbs1(i)=rbs(i)
       rbs(i)=(rb(i)/cfb(i))*(alphr**(-2))*(rbs1(i)/ksc(i))**(-1.0/3) 
       if (abs(rbs(i)-rbs1(i))<.000001) then
           exit
       end if
    end do

    cfbs(i)=(alphr**(-2.0))*(rbs(i)/ksc(i))**(-1.0/3)

  end do
          
!----------------------Bedload calculation (Ashida Michiue)-----------------
  do i=1,n
    qbt(i)=0
    do v=1,8
   
      if ((d(v)/dgf(i))>0.4) then
         taustarc(i,v)=.05*(log10(19.0)/log10(19.0*(d(v)/dgf(i))))**2
      else
         taustarc(i,v)=.05*.843*(d(v)/dgf(i))**(-1)
      end if
      taustar(i,v)=cfbs(i)*(u(i)**2)/(r*g*d(v))
      if (taustar(i,v)>taustarc(i,v)) then
         qstarb(i,v)= 17*(taustar(i,v)-taustarc(i,v))*(taustar(i,v)**.5-taustarc(i,v)**.5)
      else 
         qstarb(i,v)=0
      end if
      qb(i,v)=fa(i,v)*pc(i)*qstarb(i,v)*(r*g*d(v))**.5*d(v)
      qbt(i)=qbt(i)+qb(i,v)
    end do
    do v=1,8
      if (qbt(i)==0)then
        pb(i,v)=0
      else  
        pb(i,v)=qb(i,v)/qbt(i)
      end if
    end do
  end do


!------------End of time step loop----------------------
 end do

!------------Print out the results-----------------------

 etaaout(:,3)=etaa
 etabout(:,3)=etab
 slopeout(:,3)=slope
 dgout(:,3)=dgf
 qbout(:,3)=qbt
 hout(:,3)=h
 pcout(:,3)=pc
      

  do i=1,n

    write(1,*)  etaaout(i,1),etaaout(i,2),etaaout(i,3)
    write(2,*)  etabout(i,1),etabout(i,2),etabout(i,3)
    write(5,*)  hout(i,1), hout(i,2), hout(i,3)
    write(7,*)  slopeout(i,1),slopeout(i,2),slopeout(i,3)
    write(4,*)  dgout(i,1),dgout(i,2),dgout(i,3)
    write(8,*)  pcout(i,1),pcout(i,2),pcout(i,3)
    write(12,*) qbout(i,1),qbout(i,2),qbout(i,3)
       
  end do

 
end program Mixed_bedrock_alluvial
