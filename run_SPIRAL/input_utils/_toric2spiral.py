'''

Script to create the netCDF file needed to use TORIC data
in SPIRAL

cjperks
Oct 31, 2024

'''

# Module
import sys, os
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from transport_world.run_TORIC.output_utils import toric_tools
from transport_world.run_profiletools import eqtools3 as eq

__all__ = [
    't2s'
    ]

#################################################
#
#           Main
#
#################################################

# Converts TORIC output to SPIRAL E-field input
def t2s(
    toric_file = None,
    outpath = None,
    device = 'CMOD',
    Prf_abs = 0.5*0.3,
    # Data trimming
    gfile = None,
    afile = None,
    ):

    # Data trimming outside rhop>1
    if gfile is not None:
        dedr = eq._get_eq(
            gfile=gfile,
            afile=afile,
            machine=device,
            )
    else:
        dedr = None

    # Organizes the needed TORIC data
    ddata = _get_toric(
        toric_file = toric_file,
        device = device,
        Prf_abs = Prf_abs,
        dedr = dedr,
        )

    # Writes output file
    _write_spiral(
        outpath=outpath,
        ddata=ddata,
        )

#################################################
#
#           Utilities
#
#################################################

# Organizes the needed TORIC data
def _get_toric(
    toric_file=None,
    device=None,
    Efield_type = 'slab',
    Prf_abs = None, # [MW]
    # Data trimming
    dedr = None,
    # Plotting
    plt_all = True,
    ):

    # Reads TORIC data
    toric = toric_tools.toric_analysis(
        toric_name=toric_file,
        mode='ICRF',
        path=os.path.dirname(toric_file)+'/'
        )

    # Machine geometry
    if device == 'CMOD':
        Rmaj = 68.491/100 # [m]
        ntor = 10
        freq = 80e6 # [Hz]
    elif device == 'SPARC':
        Rmaj = 189.553/100 # [m]
        ntor = 44
        freq = 120e6 # [Hz]

    # Init
    ddata = {}

    # Scalars
    ddata['ntor'] = ntor
    ddata['freq'] = freq

    # Prepares E-field data
    if Efield_type == 'slab':
        ddata = _Efield_slab(
            toric = toric,
            ddata = ddata,
            Rmaj = Rmaj,
            Prf_abs = Prf_abs,
        )

    # Data trimming if requested
    if dedr is not None:
        ddata['Er_re'] = eq._trim_hull(
            dedr=dedr,
            RR_2d = ddata['rgrid_2d'],
            ZZ_2d = ddata['zgrid_2d'],
            data_2d = ddata['Er_re']
            )
        ddata['Er_re'][np.isnan(ddata['Er_re'])] = 0.0
        ddata['Er_im'] = eq._trim_hull(
            dedr=dedr,
            RR_2d = ddata['rgrid_2d'],
            ZZ_2d = ddata['zgrid_2d'],
            data_2d = ddata['Er_im']
            )
        ddata['Er_im'][np.isnan(ddata['Er_im'])] = 0.0
        ddata['Ez_re'] = eq._trim_hull(
            dedr=dedr,
            RR_2d = ddata['rgrid_2d'],
            ZZ_2d = ddata['zgrid_2d'],
            data_2d = ddata['Ez_re']
            )
        ddata['Ez_re'][np.isnan(ddata['Ez_re'])] = 0.0
        ddata['Ez_im'] = eq._trim_hull(
            dedr=dedr,
            RR_2d = ddata['rgrid_2d'],
            ZZ_2d = ddata['zgrid_2d'],
            data_2d = ddata['Ez_im']
            )
        ddata['Ez_im'][np.isnan(ddata['Ez_im'])] = 0.0

    if plt_all:
        fig, ax = plt.subplots(2,2, figsize = (10,8))
        fig.subplots_adjust(hspace=0.4, wspace = 0.5)

        con = ax[0,0].contourf(
            ddata['rgrid_2d'].T, ddata['zgrid_2d'].T, ddata['Er_re'].T,
            levels = 20
            )
        ax[0,0].set_title('real(Er)')
        ax[0,0].set_aspect('equal')
        cbar = fig.colorbar(con, ax=ax[0,0])
        cbar.set_label('[V/m]')
        cbar.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax[0,0].set_xlabel('R [m]')
        ax[0,0].set_ylabel('Z [m]')

        con = ax[1,0].contourf(
            ddata['rgrid_2d'].T, ddata['zgrid_2d'].T, ddata['Er_im'].T,
            levels = 20
            )
        ax[1,0].set_title('imag(Er)')
        ax[1,0].set_aspect('equal')
        cbar = fig.colorbar(con, ax=ax[1,0])
        cbar.set_label('[V/m]')
        cbar.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax[1,0].set_xlabel('R [m]')
        ax[1,0].set_ylabel('Z [m]')

        con = ax[0,1].contourf(
            ddata['rgrid_2d'].T, ddata['zgrid_2d'].T, ddata['Ez_re'].T,
            levels = 20
            )
        ax[0,1].set_title('real(Ez)')
        ax[0,1].set_aspect('equal')
        cbar = fig.colorbar(con, ax=ax[0,1])
        cbar.set_label('[V/m]')
        cbar.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax[0,1].set_xlabel('R [m]')
        ax[0,1].set_ylabel('Z [m]')

        con = ax[1,1].contourf(
            ddata['rgrid_2d'].T, ddata['zgrid_2d'].T, ddata['Ez_im'].T,
            levels = 20
            )
        ax[1,1].set_title('imag(Ez)')
        ax[1,1].set_aspect('equal')
        cbar = fig.colorbar(con, ax=ax[1,1])
        cbar.set_label('[V/m]')
        cbar.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax[1,1].set_xlabel('R [m]')
        ax[1,1].set_ylabel('Z [m]')


    # Ouput
    return ddata

# Writes SPIRAL output file
def _write_spiral(
    outpath = None,
    ddata = None,
    ):


    # Create a new netCDF file
    with Dataset(
        outpath, 
        'w', format='NETCDF4'
        ) as ncfile:

        # Define dimensions
        ncfile.createDimension('nr', ddata['nr'])
        ncfile.createDimension('nz', ddata['nz'])
        
        # Define scalar variables
        ntor = ncfile.createVariable('n_tor', 'i4')  # int
        ntor.description = "toroidal mode number [#]"
        ntor[:] = ddata['ntor']
        
        freq = ncfile.createVariable('frequency', 'f8')  # double
        freq.description = "mode frequency [Hz]"
        freq[:] = ddata['freq']

        # Define 2D grid variables
        rr = ncfile.createVariable('rr', 'f8', ('nz', 'nr'))
        rr.description = "radial grid [m]"
        rr[:] = ddata['rgrid_2d']

        zz = ncfile.createVariable('zz', 'f8', ('nz', 'nr'))
        zz.description = "vertical grid [m]"
        zz[:] = ddata['zgrid_2d']

        # Define 2D field variables (real and imaginary components)
        field_vars = {
            'Er_re': "real radial electrical field [V/m]",
            'Er_im': "imaginary radial electrical field [V/m]",
            'Ez_re': "real vertical electrical field [V/m]",
            'Ez_im': "imaginary vertical electrical field [V/m]",
            'Ephi_re': "real toroidal electrical field [V/m]",
            'Ephi_im': "imaginary toroidal electrical field [V/m]",
            'Br_re': "real radial magnetic field [T]",
            'Br_im': "imaginary radial magnetic field [T]",
            'Bz_re': "real vertical magnetic field [T]",
            'Bz_im': "imaginary vertical magnetic field [T]",
            'Bphi_re': "real toroidal magnetic field [T]",
            'Bphi_im': "imaginary toroidal magnetic field [T]"
        }

        # Populates data
        for var_name, desc in field_vars.items():
            var = ncfile.createVariable(var_name, 'f8', ('nz', 'nr'))
            var.description = desc
            var[:] = ddata[var_name]



#################################################
#
#           Conversion
#
#################################################

def _Efield_slab(
    toric = None,
    ddata = None,
    Rmaj = None,
    Prf_abs = None, # [MW]
    ):

    # Gets data
    RR = toric.cdf_hdl.variables['Xplasma'].data/100 + Rmaj
    ZZ = toric.cdf_hdl.variables['Zplasma'].data/100
    EplusRe = toric.cdf_hdl.variables['Re2Eplus'].data *Prf_abs # [V/m]
    EplusIm = toric.cdf_hdl.variables['Im2Eplus'].data *Prf_abs # [V/m]
    EminusRe = toric.cdf_hdl.variables['Re2Eminus'].data *Prf_abs # [V/m]
    EminusIm = toric.cdf_hdl.variables['Im2Eminus'].data *Prf_abs # [V/m]

    # Prepares regular mesh to interpolate data on
    nR = RR.shape[0]
    nZ = RR.shape[1]

    rlim = [
        np.min(RR.flatten()),
        np.max(RR.flatten())
        ]
    zlim = [
        np.min(ZZ.flatten()),
        np.max(ZZ.flatten())
        ]
    ZZ2, RR2 = np.meshgrid(
        np.linspace(zlim[0], zlim[1], nZ),
        np.linspace(rlim[0], rlim[1], nR),
        )

    # Interpolates E-fields
    nEpRe = LinearNDInterpolator(
        (RR.flatten(), ZZ.flatten()),
        EplusRe.flatten(),
        fill_value = 0.0
        )(
            (RR2.flatten(), ZZ2.flatten())
            ).reshape(RR2.shape)
    nEpIm = LinearNDInterpolator(
        (RR.flatten(), ZZ.flatten()),
        EplusIm.flatten(),
        fill_value = 0.0
        )(
            (RR2.flatten(), ZZ2.flatten())
            ).reshape(RR2.shape)

    nEmRe = LinearNDInterpolator(
        (RR.flatten(), ZZ.flatten()),
        EminusRe.flatten(),
        fill_value = 0.0
        )(
            (RR2.flatten(), ZZ2.flatten())
            ).reshape(RR2.shape)
    nEmIm = LinearNDInterpolator(
        (RR.flatten(), ZZ.flatten()),
        EminusIm.flatten(),
        fill_value = 0.0
        )(
            (RR2.flatten(), ZZ2.flatten())
            ).reshape(RR2.shape)

    # Stores grid ouput
    ddata['nr'] = nR
    ddata['nz'] = nZ
    ddata['rgrid_2d'] = RR2.T # [m], dim(nR, nZ)
    ddata['zgrid_2d'] = ZZ2.T # [m], dim(nR, nZ)

    # Calculates (Er,Ez) in slab limit
    Er = np.sqrt(2)/2*(
        nEpRe + 1j * nEpIm
        + (nEmRe +1j * nEmIm)
        ) # [V/m], dim(nR, nZ)
    Ez = np.sqrt(2)/2 * 1j*(
        nEpRe + 1j * nEpIm
        - (nEmRe +1j * nEmIm)
        ) # [V/m], dim(nR, nZ)

    # Stores E-field output
    ddata['Er_re'] = np.real(Er.T)
    ddata['Er_im'] = np.imag(Er.T)
    ddata['Ez_re'] = np.real(Ez.T)
    ddata['Ez_im'] = np.imag(Ez.T)

    # Remaining data
    left = [
        'Ephi_re', 'Ephi_im',
        'Br_re', 'Br_im',
        'Bz_re', 'Bz_im',
        'Bphi_re', 'Bphi_im'
        ]

    for ll in left:
        ddata[ll] = np.zeros_like(ddata['Er_re'])

    # Output
    return ddata