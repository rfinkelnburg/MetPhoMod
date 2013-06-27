      SUBROUTINE ps2str(nlevel,zen,rsfc,tauu,omu,gu,
     $     dsdh, nid, delta,
     $     fdr, fup, fdn, edr, eup, edn,
     $     fdn0, cosradsoil)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Solve two-stream equations for multiple layers.  The subroutine is based =*
*=  on equations from:  Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.=*
*=  It contains 9 two-stream methods to choose from.  A pseudo-spherical     =*
*=  correction has also been added.                                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NLEVEL  - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  ZEN     - DOUBLE PRECISION, solar zenith angle (degrees)                          (I)=*
*=  RSFC    - DOUBLE PRECISION, surface albedo at current wavelength                  (I)=*
*=  TAUU    - DOUBLE PRECISION, unscaled optical depth of each layer                  (I)=*
*=  OMU     - DOUBLE PRECISION, unscaled single scattering albedo of each layer       (I)=*
*=  GU      - DOUBLE PRECISION, unscaled asymmetry parameter of each layer            (I)=*
*=  DSDH    - DOUBLE PRECISION, slant path of direct beam through each layer crossed  (I)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (I)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*=  DELTA   - LOGICAL, switch to use delta-scaling                        (I)=*
*=            .TRUE. -> apply delta-scaling                                  =*
*=            .FALSE.-> do not apply delta-scaling                           =*
*=  FDR     - DOUBLE PRECISION, contribution of the direct component to the total     (O)=*
*=            actinic flux at each altitude level                            =*
*=  FUP     - DOUBLE PRECISION, contribution of the diffuse upwelling component to    (O)=*
*=            the total actinic flux at each altitude level                  =*
*=  FDN     - DOUBLE PRECISION, contribution of the diffuse downwelling component to  (O)=*
*=            the total actinic flux at each altitude level                  =*
*=  EDR     - DOUBLE PRECISION, contribution of the direct component to the total     (O)=*
*=            spectral irradiance at each altitude level                     =*
*=  EUP     - DOUBLE PRECISION, contribution of the diffuse upwelling component to    (O)=*
*=            the total spectral irradiance at each altitude level           =*
*=  EDN     - DOUBLE PRECISION, contribution of the diffuse downwelling component to  (O)=*
*=            the total spectral irradiance at each altitude level           =*
*=  FDN     - The incoming diffuse radiation at the model top (Added by S.P.)=*
*=  cosradsoil - The cosine of the angle between the direct sun beam and the =*
*=            normal of the soil surface                                     =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/96  Added pseudo-spherical correction                                 =*
*=  05/94  Added various two-stream methods                                  =*
*=  Nov./Dez. 1998   Adaption for MetPhoMod by S.Perego                      =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER nrows
      PARAMETER(nrows=2*kz)

*******
* input:
*******
      INTEGER nlevel
      DOUBLE PRECISION zen, rsfc
      DOUBLE PRECISION tauu(kz), omu(kz), gu(kz)
      DOUBLE PRECISION dsdh(0:kz,kz)
      INTEGER nid(0:kz)
      LOGICAL delta

*******
* output:
*******
      DOUBLE PRECISION fup(kz),fdn(kz),fdr(kz)
      DOUBLE PRECISION eup(kz),edn(kz),edr(kz)

*******
* local:
*******
      DOUBLE PRECISION tausla(0:kz), tauc(0:kz)
      DOUBLE PRECISION mu2(0:kz), mu, sum

* internal coefficients and matrix
      DOUBLE PRECISION lam(kz),taun(kz),bgam(kz)
      DOUBLE PRECISION e1(kz),e2(kz),e3(kz),e4(kz)
      DOUBLE PRECISION cup(kz),cdn(kz),cuptn(kz),cdntn(kz)
      DOUBLE PRECISION mu1(kz)
      INTEGER row
      DOUBLE PRECISION a(nrows),b(nrows),d(nrows),e(nrows),y(nrows)

*******
* other:
*******
      DOUBLE PRECISION pifs, fdn0, cosradsoil
      DOUBLE PRECISION gi(kz), omi(kz), tempg
      DOUBLE PRECISION f, g, om
      DOUBLE PRECISION gam1, gam2, gam3, gam4

* For calculations of Associated Legendre Polynomials for GAMA1,2,3,4
* in delta-function, modified quadrature, hemispheric constant,
* Hybrid modified Eddington-delta function metods, p633,Table1.
* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440, 
* uncomment the following two lines and the appropriate statements further
* down.
C     DOUBLE PRECISION YLM0, YLM2, YLM4, YLM6, YLM8, YLM10, YLM12, YLMS, BETA0,
C    >     BETA1, BETAn, amu1, subd

      DOUBLE PRECISION expon, expon0, expon1, divisr, temp, up, dn
      DOUBLE PRECISION ssfc
      INTEGER nlayer, mrows, lev

      INTEGER i, j

* Some additional program constants:
      DOUBLE PRECISION eps, precis
      PARAMETER (eps = 1.E-3, precis = 1.E-7) 
*_______________________________________________________________________

* MU = cosine of solar zenith angle
* RSFC = surface albedo
* TAUU =  unscaled optical depth of each layer
* OMU  =  unscaled single scattering albedo
* GU   =  unscaled asymmetry factor
* KLEV = max dimension of number of layers in atmosphere
* NLAYER = number of layers in the atmosphere
* NLEVEL = nlayer + 1 = number of levels

* initial conditions:  pi*solar flux = 1;  diffuse incidence = 0

      pifs = 1.      

      nlayer = nlevel - 1

       mu = COS(zen*pi/180.)

************** compute coefficients for each layer:
* GAM1 - GAM4 = 2-stream coefficients, different for different approximations
* EXPON0 = calculation of e when TAU is zero
* EXPON1 = calculation of e when TAU is TAUN
* CUP and CDN = calculation when TAU is zero
* CUPTN and CDNTN = calc. when TAU is TAUN
* DIVISR = prevents division by zero

        do j = 0, kz
           tauc(j) = 0.
           tausla(j) = 0.
           mu2(j) = 1./SQRT(largest)

        end do

       IF( .NOT. delta ) THEN
         DO i = 1, nlayer
           gi(i) = gu(i)
           omi(i) = omu(i)
           taun(i) = tauu(i)
         ENDDO
       ELSE 
* delta-scaling. Have to be done for delta-Eddington approximation, 
* delta discrete ordinate, Practical Improved Flux Method, delta function,
* and Hybrid modified Eddington-delta function methods approximations
         DO i = 1, nlayer
           f = gu(i)*gu(i)
           gi(i) = (gu(i) - f)/(1 - f)
           omi(i) = (1 - f)*omu(i)/(1 - omu(i)*f)       
           taun(i) = (1 - omu(i)*f)*tauu(i)
         ENDDO
        END IF

*
* calculate slant optical depth at the top of the atmosphere when zen>90.
* in this case, higher altitude of the top layer is recommended which can 
* be easily changed in gridz.f.
*
         IF(zen .GT. 90.0) THEN
           IF(nid(0) .LT. 0) THEN
             tausla(0) = largest
           ELSE
             sum = 0.0
             DO j = 1, nid(0)
              sum = sum + 2.*taun(j)*dsdh(0,j)
             END DO
             tausla(0) = sum 
           END IF
         END IF
  
*
        DO 11, i = 1, nlayer

         g = gi(i)
         om = omi(i)
         tauc(i) = tauc(i-1) + taun(i)

* stay away from 1 by precision.  For g, also stay away from -1

         tempg = MIN(abs(g),1. - precis)
         g = SIGN(tempg,g)
         om = MIN(om,1.-precis)


* calculate slant optical depth
*              
          IF(nid(i) .LT. 0) THEN
            tausla(i) = largest
          ELSE
            sum = 0.0
            DO j = 1, MIN(nid(i),i)
               sum = sum + taun(j)*dsdh(i,j)
            ENDDO
            DO j = MIN(nid(i),i)+1,nid(i)
               sum = sum + 2.*taun(j)*dsdh(i,j)
            ENDDO
            tausla(i) = sum 
            IF(tausla(i) .EQ. tausla(i-1)) THEN
              mu2(i) = SQRT(largest)
            ELSE
              mu2(i) = (tauc(i)-tauc(i-1))/(tausla(i)-tausla(i-1))
              mu2(i) = SIGN( MAX(ABS(mu2(i)),1./SQRT(largest)),
     $                     mu2(i) )
            END IF
          END IF
*
*** the following gamma equations are from pg 16,289, Table 1
*** save mu1 for each approx. for use in converting irradiance to actinic flux

* Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):

        gam1 =  (7. - om*(4. + 3.*g))/4.
        gam2 = -(1. - om*(4. - 3.*g))/4.
        gam3 = (2. - 3.*g*mu)/4.
        gam4 = 1. - gam3
        mu1(i) = 0.5

* quadrature (Liou, 1973, JAS, 30, 1303-1326; 1974, JAS, 31, 1473-1475):

c          gam1 = 1.7320508*(2. - om*(1. + g))/2.
c          gam2 = 1.7320508*om*(1. - g)/2.
c          gam3 = (1. - 1.7320508*g*mu)/2.
c          gam4 = 1. - gam3
c          mu1(i) = 1./sqrt(3.)
         
* hemispheric mean (Toon et al., 1089, JGR, 94, 16287):

c          gam1 = 2. - om*(1. + g)
c          gam2 = om*(1. - g)
c          gam3 = (2. - g*mu)/4.
c          gam4 = 1. - gam3
c          mu1(i) = 0.5

* PIFM  (Zdunkovski et al.,1980, Conrib.Atmos.Phys., 53, 147-166):
c         GAM1 = 0.25*(8. - OM*(5. + 3.*G))
c         GAM2 = 0.75*OM*(1.-G)
c         GAM3 = 0.25*(2.-3.*G*MU)
c         GAM4 = 1. - GAM3
c         mu1(i) = 0.5

* delta discrete ordinates  (Schaller, 1979, Contrib.Atmos.Phys, 52, 17-26):
c         GAM1 = 0.5*1.7320508*(2. - OM*(1. + G))
c         GAM2 = 0.5*1.7320508*OM*(1.-G)
c         GAM3 = 0.5*(1.-1.7320508*G*MU)
c         GAM4 = 1. - GAM3
c         mu1(i) = 1./sqrt(3.)

* Calculations of Associated Legendre Polynomials for GAMA1,2,3,4
* in delta-function, modified quadrature, hemispheric constant,
* Hybrid modified Eddington-delta function metods, p633,Table1.
* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440
c      YLM0 = 2.
c      YLM2 = -3.*G*MU
c      YLM4 = 0.875*G**3*MU*(5.*MU**2-3.)
c      YLM6=-0.171875*G**5*MU*(15.-70.*MU**2+63.*MU**4)
c     YLM8=+0.073242*G**7*MU*(-35.+315.*MU**2-693.*MU**4
c    *+429.*MU**6)
c     YLM10=-0.008118*G**9*MU*(315.-4620.*MU**2+18018.*MU**4
c    *-25740.*MU**6+12155.*MU**8)
c     YLM12=0.003685*G**11*MU*(-693.+15015.*MU**2-90090.*MU**4
c    *+218790.*MU**6-230945.*MU**8+88179.*MU**10)
c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
c      YLMS=0.25*YLMS
c      BETA0 = YLMS
c
c         amu1=1./1.7320508
c      YLM0 = 2.
c      YLM2 = -3.*G*amu1
c      YLM4 = 0.875*G**3*amu1*(5.*amu1**2-3.)
c      YLM6=-0.171875*G**5*amu1*(15.-70.*amu1**2+63.*amu1**4)
c     YLM8=+0.073242*G**7*amu1*(-35.+315.*amu1**2-693.*amu1**4
c    *+429.*amu1**6)
c     YLM10=-0.008118*G**9*amu1*(315.-4620.*amu1**2+18018.*amu1**4
c    *-25740.*amu1**6+12155.*amu1**8)
c     YLM12=0.003685*G**11*amu1*(-693.+15015.*amu1**2-90090.*amu1**4
c    *+218790.*amu1**6-230945.*amu1**8+88179.*amu1**10)
c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
c      YLMS=0.25*YLMS
c      BETA1 = YLMS
c
c         BETAn = 0.25*(2. - 1.5*G-0.21875*G**3-0.085938*G**5
c    *-0.045776*G**7)


* Hybrid modified Eddington-delta function(Meador and Weaver,1980,JAS,37,630):
c         subd=4.*(1.-G*G*(1.-MU))
c         GAM1 = (7.-3.*G*G-OM*(4.+3.*G)+OM*G*G*(4.*BETA0+3.*G))/subd
c         GAM2 =-(1.-G*G-OM*(4.-3.*G)-OM*G*G*(4.*BETA0+3.*G-4.))/subd
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = (1. - g*g*(1.- mu) )/(2. - g*g)

*****
* delta function  (Meador, and Weaver, 1980, JAS, 37, 630):
c         GAM1 = (1. - OM*(1. - beta0))/MU
c         GAM2 = OM*BETA0/MU
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = mu
*****
* modified quadrature (Meador, and Weaver, 1980, JAS, 37, 630):
c         GAM1 = 1.7320508*(1. - OM*(1. - beta1))
c         GAM2 = 1.7320508*OM*beta1
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = 1./sqrt(3.)

* hemispheric constant (Toon et al., 1989, JGR, 94, 16287):
c         GAM1 = 2.*(1. - OM*(1. - betan))
c         GAM2 = 2.*OM*BETAn
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = 0.5

*****

* lambda = pg 16,290 equation 21
* big gamma = pg 16,290 equation 22
 
         lam(i) = sqrt(gam1*gam1 - gam2*gam2)
         bgam(i) = (gam1 - lam(i))/gam2

         expon = EXP(-lam(i)*taun(i))

* e1 - e4 = pg 16,292 equation 44
         
         e1(i) = 1. + bgam(i)*expon
         e2(i) = 1. - bgam(i)*expon
         e3(i) = bgam(i) + expon
         e4(i) = bgam(i) - expon

* the following sets up for the C equations 23, and 24
* found on page 16,290
* prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
* which is approx equiv to shifting MU by 0.5*EPS* (MU)**3

         expon0 = EXP(-tausla(i-1))
         expon1 = EXP(-tausla(i))
          
         divisr = lam(i)*lam(i) - 1./(mu2(i)*mu2(i))
         temp = MAX(eps,abs(divisr))
         divisr = SIGN(temp,divisr)

         up = om*pifs*((gam1 - 1./mu2(i))*gam3 + gam4*gam2)/divisr
         dn = om*pifs*((gam1 + 1./mu2(i))*gam4 + gam2*gam3)/divisr
         
* cup and cdn are when tau is equal to zero
* cuptn and cdntn are when tau is equal to taun

         cup(i) = up*expon0
         cdn(i) = dn*expon0
         cuptn(i) = up*expon1
         cdntn(i) = dn*expon1
 
   11 CONTINUE

***************** set up matrix ******
* ssfc = pg 16,292 equation 37  where pi Fs is one (unity).changed!
      ssfc = rsfc*cosradsoil*EXP(-tausla(nlayer))*pifs

* MROWS = the number of rows in the matrix

      mrows = 2*nlayer     
      
* the following are from pg 16,292  equations 39 - 43.
* set up first row of matrix:

      i = 1
      a(1) = 0.
      b(1) = e1(i)
      d(1) = -e2(i)
      e(1) = fdn0 - cdn(i)

      row=1

* set up odd rows 3 thru (MROWS - 1):

      i = 0
      DO 20, row = 3, mrows - 1, 2
         i = i + 1
         a(row) = e2(i)*e3(i) - e4(i)*e1(i)
         b(row) = e1(i)*e1(i + 1) - e3(i)*e3(i + 1)
         d(row) = e3(i)*e4(i + 1) - e1(i)*e2(i + 1)
         e(row) = e3(i)*(cup(i + 1) - cuptn(i)) + 
     $        e1(i)*(cdntn(i) - cdn(i + 1))
   20 CONTINUE

* set up even rows 2 thru (MROWS - 2): 

      i = 0
      DO 30, row = 2, mrows - 2, 2
         i = i + 1
         a(row) = e2(i + 1)*e1(i) - e3(i)*e4(i + 1)
         b(row) = e2(i)*e2(i + 1) - e4(i)*e4(i + 1)
         d(row) = e1(i + 1)*e4(i + 1) - e2(i + 1)*e3(i + 1)
         e(row) = (cup(i + 1) - cuptn(i))*e2(i + 1) - 
     $        (cdn(i + 1) - cdntn(i))*e4(i + 1)
   30 CONTINUE

* set up last row of matrix at MROWS:

      row = mrows
      i = nlayer
      
      a(row) = e1(i) - rsfc*e3(i)
      b(row) = e2(i) - rsfc*e4(i)
      d(row) = 0.
      e(row) = ssfc - cuptn(i) + rsfc*cdntn(i)

* solve tri-diagonal matrix:

      CALL tridag(a, b, d, e, y, mrows)

**** unfold solution of matrix, compute output fluxes:

      row = 1 
      lev = 1
      j = 1
      
* the following equations are from pg 16,291  equations 31 & 32

      fdr(lev) = EXP( -tausla(0) )
      edr(lev) = mu * fdr(lev)
      edn(lev) = fdn0
      eup(lev) =  y(row)*e3(j) - y(row + 1)*e4(j) + cup(j)
      fdn(lev) = edn(lev)/mu1(lev)
      fup(lev) = eup(lev)/mu1(lev)

      DO 60, lev = 2, nlayer + 1
         fdr(lev) = EXP(-tausla(lev-1))
         edr(lev) =  mu *fdr(lev)
         edn(lev) =  y(row)*e3(j) + y(row + 1)*e4(j) + cdntn(j)
         eup(lev) =  y(row)*e1(j) + y(row + 1)*e2(j) + cuptn(j)
         fdn(lev) = edn(lev)/mu1(j)
         fup(lev) = eup(lev)/mu1(j)

         row = row + 2
         j = j + 1
   60 CONTINUE
*_______________________________________________________________________

      RETURN
      END

************************************************************************

      SUBROUTINE tridag(a,b,c,r,u,n)
*_______________________________________________________________________
* solves tridiagonal system.  From Numerical Recipies, p. 40
*_______________________________________________________________________

      IMPLICIT NONE

* input:
      INTEGER n
      DOUBLE PRECISION a, b, c, r
      DIMENSION a(n),b(n),c(n),r(n)

* output:
      DOUBLE PRECISION u
      DIMENSION u(n)

* local:
      INTEGER j

      INCLUDE 'params'
      DOUBLE PRECISION bet, gam
      DIMENSION gam(2*kz)
*_______________________________________________________________________

      IF (b(1) .EQ. 0.) STOP 1001
      bet   = b(1)
      u(1) = r(1)/bet
      DO 11, j = 2, n   
         gam(j) = c(j - 1)/bet
         bet = b(j) - a(j)*gam(j)
         IF (bet .EQ. 0.) STOP 2002 
         u(j) = (r(j) - a(j)*u(j - 1))/bet
   11 CONTINUE
      DO 12, j = n - 1, 1, -1  
         u(j) = u(j) - gam(j + 1)*u(j + 1)
   12 CONTINUE
*_______________________________________________________________________

      RETURN
      END


