'''

Recreates a dout structure from an input file

cjperks
Jan 3., 2024


'''

import numpy as np

__all__ = [
    'reload_dout'
    ]

######################################
#
#           Main
#
######################################

# Method to create dout structure from input file
def reload_dout(
    source = 'gacode',
    ddata = None,
    ):

    # Determines which type of input file
    if source == 'gacode':
        dout = reload_from_gacode(
            ddata=ddata,
            )
    else:
        print('NOT IMPLEMENTED YET!!!')

    # Output
    return dout



######################################
#
#           Utilities
#
######################################

def reload_from_gacode(
    ddata=None,
    ):

    # Modules
    from omfit_classes import omfit_eqdsk, omfit_gapy

    # Init
    dout = {}

    # Loads GACODE data
    ga = omfit_gapy.OMFITgacode(ddata['path'])

    # Input path
    dout['paths'] = {}
    dout['paths']['gacode'] = ddata['path']

    # Radial coordinates
    dout['nexp'] = str(ga['N_EXP'])
    dout['rhop'] = np.sqrt(ga['polflux']/ga['polflux'][-1])
    dout['rhot']= ga['rho']
    dout['rmin_m'] = ga['rmin']

    # Kinetic profiles
    dout['ne_19m3'] = ga['ne']
    dout['Te_keV'] = ga['Te']
    dout['Ti_keV'] = ga['Ti_1']
    dout['ptot_Pa'] = ga['ptot']
    dout['Zeff'] = ga['z_eff']
    dout['omega0_rad/s'] = ga['omega0']
    dout['vtor_m/s'] = ga['vtor_1']
    dout['vpol_m/s'] = ga['vpol_1']

    # Electron data
    dout['e'] = {}
    dout['e']['mass'] = '%0.7E'%(5.4488741e-4)
    dout['e']['charge'] = '%0.7E'%(-1)

    # Ion profiles
    dout['ions'] = {}
    dout['ions']['nion'] = str(ga['N_ION'])
    name =''
    itype=''
    mass=''
    charge=''
    for key in ga['IONS'].keys():
        ion = ga['IONS'][key][0]
        name += '%s '%(ion)
        charge += '%0.7E '%(ga['IONS'][key][1])
        mass += '%0.7E '%(ga['IONS'][key][2])
        itype += '[%s] '%(ga['IONS'][key][3])

        dout['ions'][ion] = {}
        dout['ions'][ion]['Z'] = ga['IONS'][key][1]
        dout['ions'][ion]['M'] = ga['IONS'][key][2]
        dout['ions'][ion]['ni_tot_19m3'] = ga['ni_'+str(key)]

    dout['ions']['name'] = name
    dout['ions']['itype'] = itype
    dout['ions']['mass'] = mass
    dout['ions']['charge'] = charge

    # Power profiles
    dout['powers'] = {}
    dout['powers']['sync'] = {}
    dout['powers']['sync']['prof_MW/m3'] = ga['qsync']
    dout['powers']['cont'] = {}
    dout['powers']['cont']['prof_MW/m3'] = ga['qbrem']
    dout['powers']['line'] = {}
    dout['powers']['line']['prof_MW/m3'] = ga['qline']
    dout['powers']['rad'] = {}
    dout['powers']['rad']['prof_MW/m3'] = (
        ga['pow_e_sync']
        + ga['pow_e_line']
        + ga['pow_e_brem']
        )
    dout['powers']['ohm'] = {}
    dout['powers']['ohm']['prof_MW/m3'] = ga['qohme']
    dout['powers']['rfe'] = {}
    dout['powers']['rfe']['prof_MW/m3'] = ga['qrfe']
    dout['powers']['rfi'] = {}
    dout['powers']['rfi']['prof_MW/m3'] = ga['qrfi']
    dout['powers']['rf'] = {}
    dout['powers']['rf']['prof_MW/m3'] = (
        ga['qrfe']
        + ga['qrfi']
        )
    dout['powers']['alpi'] = {}
    dout['powers']['alpi']['prof_MW/m3'] = ga['qfusi']
    dout['powers']['alpe'] = {}
    dout['powers']['alpe']['prof_MW/m3'] = ga['qfuse']
    dout['powers']['fus'] = {}
    dout['powers']['fus']['prof_MW/m3'] = (
        ga['qfuse']
        + ga['qfusi']
        )
    dout['powers']['ei'] = {}
    dout['powers']['ei']['prof_MW/m3'] = ga['qei']

    dout['powers']['beame'] = {}
    dout['powers']['beame']['prof_MW/m3'] = ga['qbeame']
    dout['powers']['beami'] = {}
    dout['powers']['beami']['prof_MW/m3'] = ga['qbeami']
    dout['powers']['ione'] = {}
    dout['powers']['ione']['prof_MW/m3'] = ga['qione']
    dout['powers']['ioni'] = {}
    dout['powers']['ioni']['prof_MW/m3'] = ga['qioni']
    dout['powers']['cxi'] = {}
    dout['powers']['cxi']['prof_MW/m3'] = ga['qcxi']

    dout['powers']['par_beam'] = {}
    dout['powers']['par_beam']['prof_MW/m3'] = ga['qpar_beam']
    dout['powers']['par_wall'] = {}
    dout['powers']['par_wall']['prof_MW/m3'] = ga['qpar_wall']
    dout['powers']['mom'] = {}
    dout['powers']['mom']['prof_MW/m3'] = ga['qmom']

    # Magnetic geometry
    dout['polflux_Wb/rad'] = ga['polflux']
    dout['q'] = ga['q']
    dout['rmaj_m'] = ga['rmaj']
    dout['zmag_m'] = ga['zmag']
    dout['kappa'] = ga['kappa']
    dout['delta'] = ga['delta']
    dout['zeta'] = ga['zeta']
    dout['shape_cos0'] = ga['shape_cos0']
    dout['shape_cos1'] = ga['shape_cos1']
    dout['shape_cos2'] = ga['shape_cos2']
    dout['shape_cos3'] = ga['shape_cos3']
    dout['shape_sin3'] = ga['shape_sin3']
    dout['johm_MA/m2'] = ga['johm']
    dout['jbs_MA/m2'] = ga['jbs']
    dout['jbstor_MA/m2'] = ga['jbstor']
    dout['jrf_MA/m2'] = ga['jrf']
    dout['jnb_MA/m2'] = ga['jnb']
    dout['vol_m3'] = ga['vol']
    dout['Bt_T'] = ga['bt0'] # ??????????

    # Scalars
    dout['torfluxa_Wb/rad'] = '%0.7E'%(0) # ?????????
    dout['bcentr_T'] = '%0.7E'%(ga['BT_EXP'])
    dout['rcentr_m'] = '%0.7E'%(0) # ?????????
    dout['current_MA'] = '%0.7E'%(ga['IP_EXP'])

    # Documentation
    dout['t0_s'] = ga['TIME']
    dout['device'] = 'UNKNOWN'
    dout['shot'] = ga['SHOT']

    # Ouptut
    return dout


