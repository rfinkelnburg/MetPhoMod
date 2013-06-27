      SUBROUTINE seto2(nz,nw,wl,cz,zen,dto2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute equivalent optical depths for O2 absorption, including absorption=*
*=  in SR bands.                                                             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - DOUBLE PRECISION, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL      - DOUBLE PRECISION, vector of lower limits of wavelength intervals in     (I)=*
*=            working wavelength grid                                        =*
*=  CZ      - DOUBLE PRECISION, number of air molecules per cm^2 at each specified    (I)=*
*=            altitude layer                                                 =*
*=  ZEN     - DOUBLE PRECISION, solar zenith angle                                    (I)=*
*=  DTO2    - DOUBLE PRECISION, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer at each specified wavelength                    =*
*=  XSO2    - DOUBLE PRECISION, molecular absorption cross section in SR bands at     (O)=*
*=            each specified altitude and wavelength.  Includes Herzberg     =*
*=            continuum.                                                     =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  03/97  Fix dto2 problem at top level (nz)                                =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax))                                               =*
*=  08/96  Modified for early exit, no redundant read of data and smaller    =*
*=         internal grid if possible;  internal grid uses user grid points   =*
*=         whenever possible                                                 =*
*=  07/96  Modified to work on internal grid and interpolate final values    =*
*=         onto the user-defined grid                                        =*
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

      DOUBLE PRECISION wl(kw)
      DOUBLE PRECISION cz(kz)
      INTEGER nz, nw
      DOUBLE PRECISION zen
      DOUBLE PRECISION dto2(kz,kw)

* grid on which Kockarts' parameterization is defined
      INTEGER ngast
      PARAMETER (ngast = 17)
      DOUBLE PRECISION wlgast(ngast)
      SAVE wlgast
 

* O2 optical depth and equivalent cross section on Kockarts' grid
      DOUBLE PRECISION dto2k(kz,ngast-1)


* internal grid and O2 cross section on internal grid
      INTEGER kdata,nwint
      PARAMETER (kdata = 200)
      DOUBLE PRECISION wlint(kdata),xso2int(kdata)
      PARAMETER (nwint = 105)

* temporary one-dimensional storage for optical depth and cross section values
* XXtmp  - on internal grid
* XXuser - on user defined grid
      DOUBLE PRECISION dttmp(2*kw), xstmp(2*kw)
      DOUBLE PRECISION dtuser(kw)

      DOUBLE PRECISION o2col(kz)

      DOUBLE PRECISION secchi
      DOUBLE PRECISION fchap
      EXTERNAL fchap

* cross section data for use outside the SR-Bands (combined from
* Brasseur and Solomon and the JPL 1994 recommendation)
      INTEGER nosr
      PARAMETER (nosr = 105)

* auxiliaries
      DOUBLE PRECISION dr
      PARAMETER (dr = pi/180.)
      INTEGER i, iw, igast
      INTEGER iz

      LOGICAL call1
      SAVE call1
      DATA call1/.TRUE./
C  Add prepared tables calculated by the original Madronich-Code
      DATA (wlint(i),i=1,nwint)/
     * 0, 116.6499, 116.65, 117.3, 117.95, 118.65, 119.4, 120.15,
     * 120.85, 121.59, 121.6, 122.35, 123.1, 123.85, 124.6, 125.4,
     * 126.2, 127, 127.8, 128.6, 129.45, 130.3, 131.15, 132, 132.85,
     * 133.75, 134.65, 135.55, 136.5, 137.45, 138.4, 139.85, 141.8,
     * 143.85, 145.95, 148.1, 150.35, 152.65, 155, 157.45, 160, 162.6,
     * 165.3, 168.1, 170.95, 173.15, 174.65, 175.4386, 176.9911,
     * 178.5714, 180.1802, 181.8182, 183.4862, 185.1852, 186.9159,
     * 188.6792, 190.4762, 192.3077, 194.1748, 196.0784, 198.0198, 200,
     * 202.0202, 204.0816, 205, 206, 207, 208, 209, 210, 211, 212, 213,
     * 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225,
     * 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238,
     * 239, 240, 241, 242, 243, 243.0001, 1e+38/
      DATA (xso2int(i),i=1,nwint)/
     * 0, 1e-20, 6.35e-19, 7.525e-19, 1.425e-19, 2.025e-19,
     * 2.412e-18, 6.4e-18,
     * 4.18e-18, 3.05e-19, 4.175e-19, 3.425e-19, 3.925e-19,
     * 8.918e-18, 9.197e-18,
     * 6.625e-19, 2.7e-19, 1.575e-19, 3.24e-19, 4.99e-19,
     * 4.875e-19, 5.525e-19, 1.068e-18,
     * 1.85e-18, 2.275e-18, 3.425e-18, 5.89e-18, 8.365e-18,
     * 1.09e-17, 1.275e-17,
     * 1.34e-17, 1.38e-17, 1.44e-17, 1.445e-17, 1.35e-17,
     * 1.22e-17, 1.071e-17, 9.075e-18,
     * 7.41e-18, 5.775e-18, 4.21e-18, 2.765e-18, 1.655e-18,
     * 9.76e-19, 5.9e-19, 3.66e-19,
     * 2.061e-19, 3.892e-20, 4.092e-21, 2.273e-21, 1.256e-21,
     * 6.901e-22, 3.749e-22,
     * 2.097e-22, 1.226e-22, 6.508e-23, 4.579e-23, 3.22e-23,
     * 2.507e-23, 1.941e-23,
     * 1.609e-23, 1.319e-23, 9.908e-24, 7.419e-24, 7.24e-24,
     * 7.09e-24, 6.955e-24,
     * 6.77e-24, 6.595e-24, 6.375e-24, 6.145e-24, 5.97e-24,
     * 5.805e-24, 5.655e-24,
     * 5.47e-24, 5.24e-24, 5.005e-24, 4.76e-24, 4.55e-24,
     * 4.36e-24, 4.175e-24, 3.99e-24,
     * 3.78e-24, 3.56e-24, 3.33e-24, 3.095e-24, 2.875e-24,
     * 2.875e-24, 2.875e-24,
     * 2.7e-24, 2.53e-24, 2.34e-24, 2.175e-24, 2.02e-24,
     * 1.86e-24, 1.705e-24, 1.555e-24,
     * 1.41e-24, 1.28e-24, 1.16e-24, 1.055e-24, 9.55e-25, 4.5e-25, 0, 0/

*-------------------------------------------------------------------------------


* check, whether user grid is in the O2 absorption band at all...
* if not, set cross section and optical depth values to zero and return

      IF (wl(1) .GT. 243.) THEN
         DO iw = 1, nw-1
           DO i = 1, nz
             dto2(i,iw) = 0.
           ENDDO
         ENDDO
         RETURN
      ENDIF

* sec Xhi or Chapman calculation
      IF (zen .LE. 75.) THEN
         secchi = 1./COS(zen*dr)
      ELSEIF (zen .LE. 95. ) THEN
         secchi = fchap(zen)
      ELSE
         RETURN
      ENDIF

* O2 overhead columns calculation
      o2col(nz-1) = 0.2095 * cz(nz-1) * secchi
      DO i = nz-2, 1, -1
        o2col(i) = o2col(i+1) + 0.2095*cz(i)*secchi
      END DO

* read O2 cross section data outside SR-bands only....etc
* this section was replaced by the above DATA statements in order
* to prevent from having to read from as file
      IF (call1) THEN
        DO iw = 1, ngast
           wlgast(iw) = 1E7/(57000.-(iw-1)*500.)
        ENDDO

        IF (call1) call1 = .FALSE.

      ENDIF


* if necessary:
* do Kockarts' parameterization of the SR bands, output values of O2
* optical depth and O2 equivalent cross section are on his grid
      IF ((wl(1) .LT. wlgast(ngast)) .AND. 
     >    (wl(nw) .GT. wlgast(1))) THEN
        DO iw = 1, ngast-1
           CALL schu(nz,o2col,iw,secchi,dto2k)
        ENDDO
      ENDIF


* loop through the altitude levels 
      DO iz = 1, nz

         igast = 0

* loop through the internal wavelength grid
         DO iw = 1, nwint-1

* if outside Kockarts' grid, use the JPL/Brasseur+Solomon data, if inside
* Kockarts' grid, use the parameterized values from the call to SCHU
           IF ((wlint(iw) .LT. wlgast(1)) .OR.
     >         (wlint(iw) .GT. wlgast(ngast-1))) THEN
              IF (iz .EQ. nz) THEN
                dttmp(iw) = 0.
              ELSE
                dttmp(iw) = xso2int(iw) * 0.2095*cz(iz)
              ENDIF
              xstmp(iw) = xso2int(iw)
           ELSE
              igast = igast+1
              dttmp(iw) = dto2k(iz,igast)
           ENDIF

* compute the area in each bin (for correct interpolation purposes only!)
           dttmp(iw) = dttmp(iw) * (wlint(iw+1)-wlint(iw))

         ENDDO

* interpolate O2 optical depth from the internal grid onto the user grid
         CALL inter3(nw,wl,dtuser, nwint,wlint,dttmp, 0)
         DO iw = 1, nw-1
            dto2(iz,iw) = dtuser(iw)/(wl(iw+1)-wl(iw))
         ENDDO
      
      ENDDO

      RETURN
      END
