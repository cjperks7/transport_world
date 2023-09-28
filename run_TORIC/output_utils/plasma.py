#!/usr/bin/env python
#
#   Try reading g-eqdsk file with re (regex) module
#   instead of the non-existant fortran format code
#   python feature.
#
#   WARNING: this code has only been testing on the 
#   two files listed below and my regex skills are
#   quite poor (as is my python ) so you have been 
#   warned.
#
#   DLG - 14-Aug-08
#
#   JCW - 27-Jan-10
#
# some improvements to regular expression and implemented direct casting
# of string arrays to floating arrays.
# Regexp are a bit fragile here. Should take advantage of the known data
# ordering with n.fromfile
#
#   JCW - 02-May-12 
#
#   Make into a function call returning a structure with optional outputs
#
# some eqdsk files have data in 16.8, some in 16.9
def readGEQDSK(filename='eqdsk.dat', dointerior=False, doplot=None, width=9):
    import re
    import numpy as n
    import pylab as p

    file = open (filename)
    data    = file.read ()

    dimensionsRE    = re.compile ( ' {1,3}\d?\d?\d?\d\d' ) # Equivilant to i5 fortran code, JCW these should be i4
    dimensionsRE5    = re.compile ( ' {1,3}\d?\d?\d?\d' ) # Equivilant to i5 fortran code
    headerRE    = re.compile ( '^.*\\n') # First line
    if width==9:
        valuesRE   = re.compile ( '([ \-]\d\.\d{9}[eEdD][\+\-]\d\d)' )   # Equivilant to e16.9 fortran code
    else:
        valuesRE   = re.compile ( '([ \-]\d\.\d{8}[eEdD][\+\-]\d\d)' )   # Equivilant to e16.8 fortran code

#bbbsRE  = re.compile ( '( {1,3}\d?\d?\d?\d\d {1,3}\d?\d?\d?\d\d)' )   # Candidate dimension lines (2i5 fortran code)
    bbbsRE  = re.compile ( r'(?m)^.{10}\n' ) #there should be only one 10 character line

    dataStr     = valuesRE.findall ( data )
    headerStr   = headerRE.findall ( data )
    bbbStr  = bbbsRE.findall ( data )

    file.close ()
    if len(bbbStr) > 0:
        nbbbsStr    = dimensionsRE5.findall ( bbbStr[0] )
    else:
        print('no bounding box find. should be Line with 2 integers length of 10 characters')
        return -1
        
    nWnHStr = dimensionsRE.findall ( headerStr[0] )

    nW  = n.int ( nWnHStr[-2] )
    nH  = n.int ( nWnHStr[-1] )

    nbbbs   = n.int ( nbbbsStr[-2] )
    limitr   = n.int( nbbbsStr[-1] )
   
    rdim    = n.float ( dataStr[0] )
    zdim    = n.float ( dataStr[1] )
    rcentr  = n.float ( dataStr[2] )
    rleft   = n.float ( dataStr[3] )
    zmid    = n.float ( dataStr[4] )

    rmaxis  = n.float ( dataStr[5] )
    zmaxis  = n.float ( dataStr[6] )
    simag   = n.float ( dataStr[7] )
    sibry   = n.float ( dataStr[8] )
    bcentr  = n.float ( dataStr[9] )

    current = n.float ( dataStr[10] )

    fpol    = n.zeros ( nW )
    pres    = n.zeros ( nW )
    ffprim  = n.zeros ( nW )
    pprime  = n.zeros ( nW )
    psizr   = n.zeros ( ( nW, nH ) )
    qpsi    = n.zeros ( nW )
    rbbbs   = n.zeros ( nbbbs )
    zbbbs   = n.zeros ( nbbbs )
    rlim    = n.zeros ( limitr )
    zlim    = n.zeros ( limitr )


#   If you know how to cast a list of strings to
#   a numpy array without a loop please let me 
#   know, as these loops should not be required.

#   1D arrays

    for i in n.arange ( nW ) : 
    
        fpol[i] = dataStr[n.cast['int'](i+20)]
        pres[i] = dataStr[n.cast['int'](i+20+nW)]
        ffprim[i] = dataStr[n.cast['int'](i+20+2*nW)]
        pprime[i] = dataStr[n.cast['int'](i+20+3*nW)]
        qpsi[i] = dataStr[n.cast['int'](i+20+4*nW+nW*nH)]

    for i in n.arange ( nbbbs ) :
    
        rbbbs[i]    = dataStr[n.cast['int'](i*2+20+5*nW+nW*nH)]
        zbbbs[i]    = dataStr[n.cast['int'](i*2+1+20+5*nW+nW*nH)]
   
    for i in n.arange ( limitr ) :
       
        rlim[i] = dataStr[n.cast['int'](i*2+20+5*nW+nW*nH+2*nbbbs)] 
        zlim[i] = dataStr[n.cast['int'](i*2+1+20+5*nW+nW*nH+2*nbbbs)] 

#   2D array

    for i in n.arange ( nW ) :
        for j in n.arange ( nH ) :
            psizr[i,j] = dataStr[n.cast['int'](i+20+4*nW+j*nW)]

    rStep   = rdim / ( nW - 1 )
    zStep   = zdim / ( nH - 1 )
    fStep   = -( simag - sibry ) / ( nW - 1 )

    r   = n.arange ( nW ) * rStep + rleft
    z   = n.arange ( nH ) * zStep + zmid - zdim / 2.0

    fluxGrid    = n.arange ( nW ) * fStep + simag

#   Find indices of points inside and outside
#   the rbbbs/zbbbs boundary.
    import matplotlib.path as mplPath
    import numpy as np
    lcf=mplPath.Path( np.column_stack( (rbbbs,zbbbs) ) )
    iiInsideA   = n.zeros ( psizr.shape )
    iiInside = -1
    iiOutside = -1
    if (dointerior):
        for i in n.arange ( nW ) :
            for j in n.arange ( nH ) :
                if lcf.contains_point( (r[i],z[i]) ):
                    iiInsideA[i,j] = 1
                #q1  = n.size ( n.where ( ( r[i] - rbbbs > 0 ) & ( z[j] - zbbbs > 0 ) ) )
                #q2  = n.size ( n.where ( ( r[i] - rbbbs > 0 ) & ( z[j] - zbbbs <= 0 ) ) )
                #q3  = n.size ( n.where ( ( r[i] - rbbbs <= 0 ) & ( z[j] - zbbbs > 0 ) ) )
                #q4  = n.size ( n.where ( ( r[i] - rbbbs <= 0 ) & ( z[j] - zbbbs <= 0 ) ) )

                #if ( q1 > 0 ) & ( q2 > 0 ) & ( q3 > 0 ) & ( q4 > 0 ) :
                #    iiInsideA[i,j]  = 1
                
        iiInside    = n.where ( iiInsideA > 0 )
        iiOutside   = n.where ( iiInsideA == 0 )

#    print nW, nH, nbbbs, limitr
#    print rdim, zdim, rcentr, rleft, zmid
#    print rmaxis, zmaxis, simag, sibry, bcentr

#   Plot output
    fig='No figure'
    if (doplot):
        N=10
        if not isinstance(doplot,bool):
            if isinstance(doplot,int):
                 N=doplot
        fig = p.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        p.contour ( r, z, psizr.T, N ,aspect='equal')
        p.plot ( rbbbs, zbbbs, 'k', linewidth = 3 )
        p.plot ( rlim, zlim, 'g', linewidth = 4 )
        p.show ()

    #checks
    # rmaxis =/ rcentr
    eqdsk = {'nW':nW, 'nH':nH, 'nbbbs':nbbbs, 'limitr':limitr, 'rdim':rdim,
             'zdim':zdim, 'rcentr':rcentr, 'rleft':rleft, 'zmid':zmid, 
             'rmaxis':rmaxis, 'zmaxis':zmaxis, 'simag':simag, 'sibry':sibry,
             'bcentr':bcentr, 'current':current, 'fpol':fpol, 'pres':pres,
             'ffprim':ffprim, 'psizr':psizr, 'qpsi':qpsi, 'rbbbs':rbbbs,
             'zbbbs':zbbbs, 'rlim':rlim, 'zlim':zlim, 'r':r, 'z':z,
             'fluxGrid':fluxGrid, 'iiInside':iiInside}

    return eqdsk,fig


def getModB(eq):
    """
    Calculate the magnitude of the magnetic field on the RZ mesh.

        |B| = \sqrt(Fpol^2+d\Psi/dZ^2+d\Psi/dR^2)/R
    """
    import numpy as np
    from scipy import interpolate

    #poloidal component
    R=eq.get('r')
    Z=eq.get('z')
    Rv,Zv=np.meshgrid(R,Z) #these are R and Z on RZ mesh
    psiRZ=np.transpose(eq.get('psizr'))
    spline_psi = interpolate.RectBivariateSpline(R,Z,psiRZ.T,bbox=[np.min(R),np.max(R),np.min(Z),np.max(Z)],kx=5,ky=5)
    psi_int_r=spline_psi.ev(Rv,Zv,dx=1)
    psi_int_z=spline_psi.ev(Rv,Zv,dy=1)
    grad_psi=np.sqrt(psi_int_z**2+psi_int_r**2)
    
    #toroidal component
    #get Fpol and interpolate to RZ mesh to get fpolRZ
    fpol=eq.get('fpol')
    psi=eq.get('fluxGrid') #mesh for fpol
    #fpolRZ=fpol(psiRZ)
    # scipy 0.18    spline_fpol=interpolate.CubicSpline(psi,fpol,bc_type='natural')
    #fill value is B0*R0
    spline_fpol=interpolate.interp1d(psi,fpol,bounds_error=False,fill_value=fpol[-1],kind='cubic')
    fpolRZ=[] #np.zeros(psiRZ.shape)
    for psirow in psiRZ:
        fpolRZ.append( spline_fpol(psirow) )
    fpolRZ=np.array(fpolRZ) #Fpol numpy array on RZ mesh

    modB=np.sqrt(grad_psi**2+fpolRZ**2)/Rv
    #Add components
    return modB,grad_psi,fpolRZ,Rv,Zv


def getLCF(eq):
    #find which contour in LCF, same as rbbbs? 

    import matplotlib.path as mplPath
    import numpy as np
    import pylab as p

    R=eq.get('r')
    Z=eq.get('z')
    psiRZ=np.transpose(eq.get('psizr'))    
    CSlcf=p.contour(R,Z,psiRZ,levels=[eq['sibry']-.01])
    cntr=(eq['rmaxis'],eq['zmid'])
    lcf=(0,0)
    for p in CSlcf.collections[0].get_paths():
        v = p.vertices
        x = v[:,0]
        y = v[:,1]
        bbPath = mplPath.Path(np.column_stack( (x,y)))
        if bbPath.contains_point(cntr):
            lcf=(x,y)
            return lcf
