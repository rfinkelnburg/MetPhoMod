      SUBROUTINE rtlink(nz,z,
     $     ag, zen, 
     $     dsdh, nid,
     $     dtrl, 
     $     dto3, 
     $     dto2, 
     $     dtso2, 
     $     dtno2, 
     $     dtcld, omcld, gcld,
     $     dtaer,omaer,gaer,
     $     edir, edn, eup, fdir, fdn, fup, 
     $     f0z, cosradsoil)
*_______________________________________________________________________

      IMPLICIT NONE

      INCLUDE 'params'

* input:

      INTEGER nz
      DOUBLE PRECISION z(kz)
      DOUBLE PRECISION ag
      DOUBLE PRECISION zen
      DOUBLE PRECISION dtrl(kz)
      DOUBLE PRECISION dto3(kz), dto2(kz), dtso2(kz), dtno2(kz)
      DOUBLE PRECISION dtcld(kz), omcld(kz), gcld(kz)
      DOUBLE PRECISION dtaer(kz), omaer(kz), gaer(kz)
      DOUBLE PRECISION dsdh(0:kz,kz)
      DOUBLE PRECISION f0z, cosradsoil
      INTEGER nid(0:kz)


* output
      DOUBLE PRECISION edir(kz), edn(kz), eup(kz)
      DOUBLE PRECISION fdir(kz), fdn(kz), fup(kz)

* more program constants:
      DOUBLE PRECISION dr
      PARAMETER (dr = pi/180.)

* local:
      DOUBLE PRECISION dt(kz), om(kz), g(kz)
      DOUBLE PRECISION ediri(kz), edni(kz), eupi(kz)
      DOUBLE PRECISION fdiri(kz), fdni(kz), fupi(kz)
      DOUBLE PRECISION daaer, dtsct, dtabs, dsaer, dscld, dacld
      INTEGER i, ii
      LOGICAL delta

      DATA delta /.true./
*_______________________________________________________________________

* initialize:

      DO 5 i = 1, nz
         fdir(i) = 0.
         fup(i) = 0.
         fdn(i) = 0.
         edir(i) = 0.
         eup(i) = 0.
         edn(i) = 0.
 5    CONTINUE

*  set here any coefficients specific to rt scheme, 
* ----

      DO 10, i = 1, nz - 1

         dscld = dtcld(i)*omcld(i)
         dacld = dtcld(i)*(1.-omcld(i))

         dsaer = dtaer(i)*omaer(i)
         daaer = dtaer(i)*(1.-omaer(i))

         dtsct = dtrl(i) + dscld + dsaer
         dtabs = dto3(i) + dto2(i) + dtso2(i) + 
     >           dtno2(i) + dacld + daaer

 	 dtabs = MAX(dtabs,1./largest)
 	 dtsct = MAX(dtsct,1./largest)

* invert z-coordinate:

         ii = nz - i
         dt(ii) = dtsct + dtabs
         om(ii) = dtsct/(dtsct + dtabs)
           IF (dtsct .EQ. 1./largest) om(ii) = 1./largest
         g(ii) = (gcld(i)*dscld + gaer(i)*dsaer)/dtsct

   10 CONTINUE

*  call rt routine:

      CALL ps2str(nz,zen,ag,dt,om,g,
     $     dsdh, nid, delta,
     $     fdiri, fupi, fdni, ediri, eupi, edni,
     $     f0z, cosradsoil)

* put on upright z-coordinate

      DO 20, i = 1, nz
         ii = nz - i + 1
         fdir(i) = fdiri(ii)
         fup(i) = fupi(ii)
         fdn(i) = fdni(ii)
         edir(i) = ediri(ii)
         eup(i) = eupi(ii)
         edn(i) = edni(ii)
 20   CONTINUE
*_______________________________________________________________________

      RETURN
      END   
