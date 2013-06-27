      FUNCTION fchap(zeta)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the Chapman function.  Chapman function is used when the solar =*
*=  zenith angle exceeds 75 degrees, but is not greater than 95 degrees.     =*
*=  The function value is calculated by interpolate between values given in, =*
*=  e.g., McCartney (1976).                                                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  ZETA - DOUBLE PRECISION, solar zenith angle (degrees)                             (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  01/95  Original                                                          =*
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

* input:
      DOUBLE PRECISION zeta

* function value:
      DOUBLE PRECISION fchap

* internal:
      DOUBLE PRECISION rm
      INTEGER j
      DOUBLE PRECISION y(21)
*_______________________________________________________________________

      DATA y /
     1     3.800,4.055,4.348,4.687,5.083,5.551,6.113,
     2     6.799,7.650,8.732,10.144,12.051,14.730,18.686,
     3     24.905,35.466,55.211,96.753,197.,485.,1476./
*_______________________________________________________________________

      j = MAX(INT(zeta)+1,75)
      rm = FLOAT(j)

      fchap = y(j-75) + (y(j-74) - y(j-75))*(zeta - (rm-1.))
*_______________________________________________________________________

      RETURN
      END
