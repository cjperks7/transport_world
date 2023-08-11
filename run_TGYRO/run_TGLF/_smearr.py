'''

Adapted from the D, V smooting function, smearr,

Options are:
    1) Simple smoothing to remove single point outliers
        that jump in magnitude or flip sign

    2) what's used in ASTRA (sglazh)--

!----------------------------------------------------------------------|
!  Subroutine minimizes the value of functional
!  INTEGRAL(alfa*(dU/dx)**2+(U-F)**2)*dx with respect to U(x).
!  f_in(1:NA1) is a given array on the grid X(1:NA1)
!  The result is a smoothed array f_out(1:NA1) given on the same grid
!	ALFA ~ 0.01*X(NA1)**2  is a regularizator
!	The target function f_out obeys the additional conditions:
!	     df_out/dx(x=0)=0 - cylindrical case
!	     f_out(XN(NA1))=f_in(XO(NA1))
!----------------------------------------------------------------------|


'''

# Modules
import numpy as np

__all__ = [
    'smearr',
    ]

##########################################################
#
#                   Main
#
##########################################################

def smearr(
    alfa = None, 
    f_in = None,
    rho = None,
    alg = None,
    ):

    # Runs smoothing
    if alg == 'sglazh':
        f_out = _sglazh(
            alfa=alfa,
            n0 = len(rho)-1,
            f_in = f_in,
            coor = rho,
            )

    elif alg == 'simple':
        f_out = _simple(
            alfa=alfa,
            f_in = f_in,
            coor = rho,
            )


    return f_out



##########################################################
#
#                   Utilities
#
##########################################################

def _simple(
    alfa = None,
    f_in = None,
    coor = None,
    ):

    # Initializes
    f_out = f_in.copy()

    # Loop over radial coordinate
    for jj in np.flip(np.arange(1,len(f_out)-1)):
        if np.sign(f_in[jj-1]) == np.sign(f_in[jj+1]):
            if np.sign(f_out[jj]) != np.sign(f_in[jj-1]):
                f_out[jj] = (f_in[jj-1]+ f_in[jj+1])/2
        
        if (abs(f_out[jj]) >= abs(alfa * (f_out[jj-1] + f_out[jj+1])/2)):
            print(f_out[jj])
            print((f_in[jj-1] + f_in[jj+1])/2)
            f_out[jj] = (f_in[jj-1]+ f_in[jj+1])/2

    return f_out 

def _sglazh(
    alfa=None,      # Usually ~ 0.01 *rho[-1]**2
    n0 = None,
    f_in = None,    # Quantity to smooth
    coor = None,    # Coordinate to smooth over
    ):

    # Intializes arrays
    f_out = np.zeros(len(f_in))
    PP = np.zeros(len(f_in))

    # Loop over coordinate
    for jj in np.arange(1,n0):
        PP[jj] = alfa/(coor[jj] - coor[jj-1])/coor[n0]**2
    PP[0] = 0

    # Fixed first value
    f_out[0] = f_in[0]

    # Initializes
    II = 0
    YF = (f_in[1]- f_in[0])/(coor[1]-coor[0])
    YX = 2/(coor[1]+coor[0])
    YP = 0.0
    YQ = 0.0

    for jj in np.arange(n0-1):
        if (coor[II] <= coor[jj]):
            II = II + 1
            II = min(II, n0)

            if (II == n0 or coor[II] >= coor[jj]):
                break

            YF = (f_in[II] - f_in[II-1])/(coor[II] - coor[II-1])

        FJ = f_in[II] + YF*(coor[jj]- coor[II])

        YD = 1 + YX*(YP + PP[jj+1])

        PP[jj] = YX*PP[jj+1]/YD

        f_out[jj] = (FJ + YX*YQ)/YD

        YX = 2/(coor[jj+2] - coor[jj])

        YP = (1 - PP[jj]) *PP[jj+1]

        YQ = f_out[jj] * PP[jj+1]


    f_out[n0] = f_in[n0]

    for jj in np.flip(np.arange(n0)):
        f_out[jj] = PP[jj] * f_out[jj+1] + f_out[jj]

    return f_out
