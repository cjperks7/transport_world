'''

_profs_cmod is a fucntion to handle obtain C-Mod 
experimental kinetic profiles

'''

# Modules
import MDSplus
from scipy.interpolate import interp1d
import numpy as np

__all__ = [
    'profs_cmod',
    ]

###########################################
#
#               Main
#
###########################################

def profs_cmod(
    shot=None,
    ):

    # Initializes output dictionary
    dkin = {}

    # Obtains electron temperature
    dkin['Te'] = _get_Te(
        shot=shot,
        )

    # Obtains ion temperature
    dkin['Ti'] = _get_Ti(
        shot=shot,
        )

    # Obtains electron density
    dkin['ne'] = _get_ne(
        shot=shot,
        )

    # Output
    return dkin

###########################################
#
#               Extra
#
###########################################

def _get_ne(
    shot=None,
    ):
    '''
    NOTE: does not include edge probes or edge Thomson,
    just core diagnostics
    see w_dens
    '''

    # Obtains data tree
    tree = MDSplus.Tree('electrons', shot)

    # Initializes output
    dne = {}

    # ------------------
    # Obtains TCI data
    # ------------------
    dne['TCI'] = {}

    try:
        # line-of-sight characteristic length
        LOS = 0.6 # [m]

        # Obtain line-average electron density data
        dne['TCI']['val_1e20m3'] = tree.getNode('\electrons::top.tci.results:nl_04').data()/LOS # [1e20 1/m2], dim(t,)

        # Time domain
        dne['TCI']['t_s'] = tree.getNode('\electrons::top.tci.results:nl_04').dim_of(0).data() # [s], dim(t,)

    # In case this diag is not available
    except:
        print('TCI for nebar not available')

    # ------------------
    # Obtains Thomson data
    # ------------------

    '''
    NOTE: Be aware that there was an old YAG system
    before shot 1070807001
    '''

    dne['Thomson'] = {}

    try:
        # Obtains electron density data
        dne['Thomson']['val_1e20m3'] = tree.getNode('\yag_new.results.profiles:ne_rz').data() # [1e20 1/m3], dim(nchan,t)

        # Obtains error
        dne['Thomson']['err_1e20m3'] = tree.getNode('\yag_new.results.profiles:ne_err').data() # [1e20 1/m3], dim(nchan,t)

        # Radial domain
        dne['Thomson']['r_m'] = tree.getNode('\yag_new.results.profiles:r_mid_t').data() # [m], dim(nchan,t)

        # Time domain
        dne['Thomson']['t_s'] = tree.getNode('\yag_new.results.profiles:ne_rz').dim_of(0).data() # [s], dim(t,)


    # In case this diag is not available
    except:
        print('Thomson for ne not available')

    # ------------------
    # Obtains edge Thomson data
    # ------------------

    dne['Edge'] = {}

    try:
        # Obtains electron density data
        dne['Edge']['val_1e20m3'] = tree.getNode('\electrons::top.yag_edgets.results:ne').data()/1e20 # [1e20 1/m3], dim(nchan,t)

        # Radial domain
        dne['Edge']['psin'] = tree.getNode('\electrons::top.yag_edgets.results:psinorm').data() # dim(nchan,t)

        # Time domain
        dne['Edge']['t_s'] = tree.getNode('\electrons::top.yag_edgets.results:ne').dim_of(0).data() # [s], dim(t,)

    # In case this diag is not available
    except:
        print('Edge Thomson for ne not available')


    #Output
    return dne

    
def _get_Ti(
    shot=None,
    ):

    # Obtains data tree
    tree = MDSplus.Tree('spectroscopy', shot)

    # Initializes output
    dTi = {}

    # ------------------
    # Obtains HIREXSR z line data
    # ------------------
    dTi['HIREXSR_z'] = {}

    try:
        # Obtain ion temp data
        dTi['HIREXSR_z']['val_keV'] = tree.getNode(
            r'\spectroscopy::top.hirexsr.analysis.helike.profiles.z:pro'
            ).data()[3,:,:].T # [keV], dim(nchan,t)

        # Obtain ion temp errorbar data
        dTi['HIREXSR_z']['err_keV'] = tree.getNode(
            r'\spectroscopy::top.hirexsr.analysis.helike.profiles.z:proerr'
            ).data()[3,:,:].T # [keV], dim(nchan,t)

        # Radial domain
        dTi['HIREXSR_z']['psin'] = tree.getNode(
            r'\spectroscopy::top.hirexsr.analysis.helike.profiles.z:rho'
            ).data().T # [], dim(mchan,t)

        # Time domain
        dTi['HIREXSR_z']['t_s'] = tree.getNode(
            r'\spectroscopy::top.hirexsr.analysis.helike.profiles.z:rho'
            ).dim_of(0).data() # [s], dim(t,)

        # Error check
        if dTi['HIREXSR_z']['val_keV'].shape[0] > dTi['HIREXSR_z']['psin'].shape[0]:
            ss = dTi['HIREXSR_z']['psin'].shape[0]
            dTi['HIREXSR_z']['val_keV'] = dTi['HIREXSR_z']['val_keV'][:ss,:]
            dTi['HIREXSR_z']['err_keV'] = dTi['HIREXSR_z']['err_keV'][:ss,:]

    # In case this diag is not available
    except:
        print('HIREXSR z line for Ti not available')

    # ------------------
    # Obtains HIREXSR Lya1 line data
    # ------------------
    dTi['HIREXSR_Lya1'] = {}

    try:
        # Obtain ion temp data
        dTi['HIREXSR_Lya1']['val_keV'] = tree.getNode(
            r'\spectroscopy::top.hirexsr.analysis.hlike.profiles.lya1:pro'
            ).data()[3,:,:].T # [keV], dim(nchan,t)

        # Obtain ion temp errorbar data
        dTi['HIREXSR_Lya1']['err_keV'] = tree.getNode(
            r'\spectroscopy::top.hirexsr.analysis.hlike.profiles.lya1:proerr'
            ).data()[3,:,:].T # [keV], dim(nchan,t)

        # Radial domain
        dTi['HIREXSR_Lya1']['psin'] = tree.getNode(
            r'\spectroscopy::top.hirexsr.analysis.hlike.profiles.lya1:rho'
            ).data().T # [], dim(mchan,t)

        # Time domain
        dTi['HIREXSR_Lya1']['t_s'] = tree.getNode(
            r'\spectroscopy::top.hirexsr.analysis.hlike.profiles.lya1:rho'
            ).dim_of(0).data() # [s], dim(t,)

    # In case this diag is not available
    except:
        print('HIREXSR Lya1 line for Ti not available')


    # Output
    return dTi
    

def _get_Te(
    shot=None,
    ):
    '''
    NOTE: does not include edge probes or edge Thomson,
    just core diagnostics
    see w_temp
    '''

    # Obtains data tree
    tree = MDSplus.Tree('electrons', shot)

    # Initializes output
    dTe = {}

    # ------------------
    # Obtains GPC data
    # ------------------
    dTe['GPC'] = {}

    try:
        # Number of channels
        nchan = 9

        # Time domain
        dTe['GPC']['t_s'] = tree.getNode(
            r'\electrons::GPC_TE1'
            ).dim_of(0).data() # dim(t,)

        # Initiaizes data array
        dTe['GPC']['val_keV'] = np.zeros((nchan, len(dTe['GPC']['t_s']))) # [keV], dim(nchan,t)
        dTe['GPC']['r_m'] = np.zeros((nchan, len(dTe['GPC']['t_s']))) # [m], dim(nchan,t)

        # Loop over channels
        for ii in np.arange(nchan):
            # Obtains measured electron temp
            dTe['GPC']['val_keV'][ii,:] = tree.getNode(
                r'\electrons::GPC_TE'+str(ii+1)
                ).data() # [keV]

            # Obtains measured radial location
            r_tmp = tree.getNode(
                r'\electrons::top.ece.gpc_results.rad:r'+str(ii+1)
                ).data() # [m]
            t_tmp = tree.getNode(
                r'\electrons::top.ece.gpc_results.rad:r'+str(ii+1)
                ).dim_of(0).data() # [s]

            # Interpolates radial domain
            dTe['GPC']['r_m'][ii,:] = interp1d(
                t_tmp,
                r_tmp,
                bounds_error = False,
                fill_value = (r_tmp[0], r_tmp[-1])
                )(dTe['GPC']['t_s'])

    # In case this diag is not available
    except:
        print('GPC for Te not available')

    # ------------------
    # Obtains GPC2 data
    # ------------------

    dTe['GPC2'] = {}

    try:
        # Obtains measured electron temp
        dTe['GPC2']['val_keV'] = tree.getNode('\gpc2_te').data() # [keV], dim(nchan, t)

        # Time domain
        dTe['GPC2']['t_s'] = tree.getNode('\gpc2_te').dim_of(0).data() # [s], dim(t,)

        # Radial domain
        r_tmp = tree.getNode('\gpc2_r').data() # [m], dim(nchan, t_tmp)
        t_tmp = tree.getNode('\gpc2_r').dim_of(0).data() # [s], dim(t_tmp,)

        # Interpolates onto time domain
        dTe['GPC2']['r_m'] = interp1d(
            t_tmp,
            r_tmp,
            axis = 1,
            bounds_error = False,
            fill_value = (r_tmp[0,0], r_tmp[0,-1])
            )(dTe['GPC2']['t_s'])

    # In case this diag is not available
    except:
        print('GPC2 for Te not available')

    # ------------------
    # Obtains FRC data
    # ------------------

    dTe['FRC'] = {}

    try:
        # Number of channels
        nchan = 32

        # Time domain
        dTe['FRC']['t_s'] = tree.getNode(
            r'\electrons::te_hrece01'
            ).dim_of(0).data() # dim(t,)

        # Initiaizes data array
        dTe['FRC']['val_keV'] = np.zeros((nchan, len(dTe['FRC']['t_s']))) # [keV], dim(nchan,t)
        dTe['FRC']['r_m'] = np.zeros((nchan, len(dTe['FRC']['t_s']))) # [m], dim(nchan,t)

        # Loop over channels
        for ii in np.arange(nchan):
            if ii < 9:
                num = '0'+str(ii+1)
            else:
                num = str(ii+1)
            # Obtains measured electron temp
            dTe['FRC']['val_keV'][ii,:] = tree.getNode(
                r'\electrons::te_hrece'+num
                ).data() # [keV]

            # Obtains measured radial location
            r_tmp = tree.getNode(
                r'\electrons::rmid_hrece'+num
                ).data()[:,0] # [m]
            t_tmp = tree.getNode(
                r'\electrons::rmid_hrece'+num
                ).dim_of(0).data() # [s]

            # Interpolates radial domain
            dTe['FRC']['r_m'][ii,:] = interp1d(
                t_tmp,
                r_tmp,
                bounds_error = False,
                fill_value = (r_tmp[0], r_tmp[-1])
                )(dTe['FRC']['t_s'])

    # In case this diag is not available
    except:
        print('FRC for Te not available')

    # ------------------
    # Obtains Thomson data
    # ------------------

    '''
    NOTE: Be aware that there was an old YAG system
    before shot 1070807001
    '''

    dTe['Thomson'] = {}

    try:
        # Obtains electron temp data
        dTe['Thomson']['val_keV'] = tree.getNode('\yag_new.results.profiles:te_rz').data() # [keV], dim(nchan,t)

        # Obtains error
        dTe['Thomson']['err_keV'] = tree.getNode('\yag_new.results.profiles:te_err').data() # [keV], dim(nchan,t)

        # Radial domain
        dTe['Thomson']['r_m'] = tree.getNode('\yag_new.results.profiles:r_mid_t').data() # [m], dim(nchan,t)

        # Time domain
        dTe['Thomson']['t_s'] = tree.getNode('\yag_new.results.profiles:te_rz').dim_of(0).data() # [s], dim(t,)

    # In case this diag is not available
    except:
        print('Thomson for Te not available')


    # ------------------
    # Obtains edge Thomson data
    # ------------------

    dTe['Edge'] = {}

    try:
        # Obtains electron density data
        dTe['Edge']['val_keV'] = tree.getNode('\electrons::top.yag_edgets.results:te').data()/1e3 # [keV], dim(nchan,t)

        # Radial domain
        dTe['Edge']['psin'] = tree.getNode('\electrons::top.yag_edgets.results:psinorm').data() # dim(nchan,t)

        # Time domain
        dTe['Edge']['t_s'] = tree.getNode('\electrons::top.yag_edgets.results:ne').dim_of(0).data() # [s], dim(t,)

    # In case this diag is not available
    except:
        print('Edge Thomson for ne not available')


    # ------------------
    # Obtains Michelson data
    # ------------------

    dTe['Michelson'] = {}

    try:
        # Obtain electron temp data
        dTe['Michelson']['val_keV'] = tree.getNode('\electrons::ece_te').data().T # [keV] #dim(nchan,t)

        # Radial domain
        dTe['Michelson']['r_m'] = tree.getNode('\electrons::ece_te').dim_of(0).data() # [keV] #dim(nchan,)

        # Time domain
        dTe['Michelson']['t_s'] = tree.getNode('\electrons::ece_te').dim_of(1).data() # [keV] #dim(t,)

    # In case this diag is not available
    except:
        print('Michelson for Te not available')


    # Output
    return dTe
