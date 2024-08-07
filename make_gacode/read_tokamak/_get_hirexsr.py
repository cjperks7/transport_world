'''

_get_hirexsr.py is a module to help faciliate
loading spectroscopy data from HIREX SR on C-Mod

'''

# Modules
import MDSplus
import matplotlib.pyplot as plt
import numpy as np
import copy
from scipy.interpolate import interp1d

plt.rcParams.update({'font.size': 16})

__all__ = [
    'get_hirexsr',
    'plt_hirexsr',
    ]

###########################################
#
#               Main
#
###########################################

def get_hirexsr(
    shot=None,
    tht = 0,
    quants=None, # None -> ['int', 'spectra', 'moments']
    plt_all=None,
    plt_ch=None,
    ):

    # Initializes output
    dout = {}
    dout['shot'] = shot
    dout['tht'] = tht

    # Defaults
    if quants is None:
        quants = ['int', 'spectra', 'moments', 'profiles']

    # Gets intensity time traces
    if 'int' in quants:
        dout = _get_int(
            dout=dout,
            )

    # Gets spectra 
    if 'spectra' in quants:
        dout = _get_spectra(
            dout=dout,
            )

    # Gets moments of the signal
    if 'moments' in quants:
        dout = _get_moments(
            dout=dout,
            )

    # Get inversion profiles
    if 'profiles' in quants:
        dout = _get_profs(
            dout=dout,
            )

    # Plotting
    if plt_all:
        plt_hirexsr(dout=dout, plt_ch=plt_ch)

    return dout

###########################################
#
#               Plotting
#
###########################################

def plt_hirexsr(
    dout=None,
    plt_ch=None,
    dplt_mom={
        't1': [0.75, 0.95],
        't2': [1.15,1.35],
        },
    ):

    # Plots intensity versus time
    if 'int' in list(dout.keys()):
        fig1, ax1 = plt.subplots(1,2)
        fig1.suptitle('Midplane Intensity')
        fig1.tight_layout(pad=1.0)

        # Loop over lines
        for ll in dout['int'].keys():
            if ll in ['A', 'M']:
                num = 0
                data = dout['int'][ll]['data']
            elif ll in ['W','Z']:
                num = 1
                data = dout['int'][ll]['data'][0,:]

            ax1[num].plot(
                dout['int'][ll]['t_s'],
                data,
                label = ll
                )  
        
        ax1[0].set_xlim(0.5,1.5)
        ax1[1].set_xlim(0.5,1.5)

        leg1 = ax1[0].legend()
        leg2 = ax1[1].legend()
        leg1.set_draggable('on')
        leg2.set_draggable('on')

        ax1[0].grid('on')
        ax1[1].grid('on')

        ax1[0].set_xlabel('t [s]')
        ax1[1].set_xlabel('t [s]')

        ax1[0].set_ylabel('Intensity [arb]')
        ax1[1].set_ylabel('Intensity [arb]')

    # Plots averaged spectra
    if 'spectra' in list(dout.keys()):
        from transport_world.plot_utils import plt_slider
        #fig2, ax2 = plt.subplots(1,2)

        # Loop over branches
        for bb in dout['spectra'].keys():
            # Determines which channel to plot
            if plt_ch is None:
                # If not selected, plot midplane
                plt_ch_tmp = np.argmin(
                    abs(dout['spectra'][bb]['pos'][98,:,3])
                    )
            else:
                if plt_ch > dout['spectra'][bb]['pos'].shape[1] -1:
                    plt_ch_tmp = dout['spectra'][bb]['pos'].shape[1] -1 
                else:
                    plt_ch_tmp = copy.deepcopy(plt_ch)
            plt_slider(
                xxx=dout['spectra'][bb]['lam_AA'][:,0,plt_ch_tmp],
                yyy=dout['spectra'][bb]['t_s'],
                dzzz=dout['spectra'][bb]['avespec'][:,:,plt_ch_tmp].T,
                xxlabel=r"$\lambda$ [$\AA$]",
                yylabel="t [s] (ch="+str(plt_ch_tmp)+")",
                zzlabel=r"Counts [arb]",
                plt_sum=False,
                yscale='log',
                )

    # Plots brightness ratio
    if 'mom' in list(dout.keys()) and dplt_mom is not None:
        fig3, ax3 = plt.subplots(1, len(dout['mom'].keys()))
        fig3.tight_layout(pad=1.0)

        # Loop over lines
        for ll, line in enumerate(list(dout['mom'].keys())):
            # Time window to average data
            ind1 = np.where(
                (dout['mom'][line]['t_s'] >= dplt_mom['t1'][0])
                & (dout['mom'][line]['t_s'] <= dplt_mom['t1'][-1])
                )[0]
            ind2 = np.where(
                (dout['mom'][line]['t_s'] >= dplt_mom['t2'][0])
                & (dout['mom'][line]['t_s'] <= dplt_mom['t2'][-1])
                )[0]

            # Counts data
            data1 = np.mean(dout['mom'][line]['data'][0,ind1,:], axis=-2) # [], dim(nchan,)
            data2 = np.mean(dout['mom'][line]['data'][0,ind2,:], axis=-2) # [], dim(nchan,)

            psin1 = np.mean(dout['mom'][line]['psin'][ind1,:], axis=-2) # [], dim(nchan,)
            psin2 = np.mean(dout['mom'][line]['psin'][ind2,:], axis=-2) # [], dim(nchan,)

            ax3[ll].plot(
                psin1,
                1- (data2/data1),
                'r*'
                )

            ax3[ll].set_xlabel(r'$Psi_n$')
            ax3[ll].set_ylabel('Brightness reduction')
            ax3[ll].set_title(line)
            ax3[ll].set_ylim(0,1)
            ax3[ll].grid('on')

    # Plots emissivity ratio
    if 'profs' in list(dout.keys()) and dplt_mom is not None:
        fig4, ax4 = plt.subplots(1, len(dout['profs'].keys()))
        fig4.tight_layout(pad=1.0)

        # Loop over lines
        for ll, line in enumerate(list(dout['profs'].keys())):
            # Time window to average data
            ind1 = np.where(
                (dout['profs'][line]['t_s'] >= dplt_mom['t1'][0])
                & (dout['profs'][line]['t_s'] <= dplt_mom['t1'][-1])
                )[0]
            ind2 = np.where(
                (dout['profs'][line]['t_s'] >= dplt_mom['t2'][0])
                & (dout['profs'][line]['t_s'] <= dplt_mom['t2'][-1])
                )[0]

            # Counts data
            data1 = np.mean(
                dout['profs'][line]['emis']['data'][:,ind1], 
                axis=-1
                ) # [], dim(nrho,)
            data2 = np.mean(
                dout['profs'][line]['emis']['data'][:,ind2], 
                axis=-1
                ) # [], dim(nrho,)

            psin1 = np.mean(dout['profs'][line]['psin'][:,ind1], axis=-1) # [], dim(nrho,)
            psin2 = np.mean(dout['profs'][line]['psin'][:,ind2], axis=-1) # [], dim(nrho,)

            ax4[ll].plot(
                psin1,
                1- (data2/data1),
                'r*'
                )

            ax4[ll].set_xlabel(r'$\sqrt{\Psi_n}$')
            ax4[ll].set_ylabel('Emissivity reduction')
            ax4[ll].set_title(line)
            ax4[ll].set_ylim(0,1)
            ax4[ll].grid('on')


###########################################
#
#               Utilities
#
###########################################

# Gets inversion profiles (emissivity, Ti, vtor) v. psin
def _get_profs(
    dout=None,
    ):

    # Initializes output
    dout['profs'] = {}

    # Available data
    lines = ['W', 'X', 'Z', 'LYA1', 'J', 'MO4D']

    # Loop over lines
    for line in lines:
        if line in ['W', 'X', 'Z']:
            bnch = 'HELIKE'
        elif line in ['LYA1', 'J', 'MO4D']:
            bnch = 'HLIKE'

        # MDSplus tree
        _, branchNode = _cmod_tree(dout=dout, branch=bnch)

        try:
            # get time basis
            all_times = branchNode.getNode(
                'profiles.'+line+':rho'
                ).dim_of(0).data() # [s], dim(t,)

            dout['profs'][line] = {}
            dout['profs'][line]['emis'] = {}
            dout['profs'][line]['Ti'] = {}
            dout['profs'][line]['vtor'] = {}

            mask = all_times > -1
            dout['profs'][line]['t_s'] = all_times[mask]

            # Obtain inversion data
            dout['profs'][line]['emis']['data'] = branchNode.getNode(
                'profiles.'+line+':pro'
                ).data()[0,mask,:].T # dim(nrho,t)
            dout['profs'][line]['Ti']['data'] = branchNode.getNode(
                'profiles.'+line+':pro'
                ).data()[3,mask,:].T # [keV], dim(nrho,t)
            dout['profs'][line]['vtor']['data'] = branchNode.getNode(
                'profiles.'+line+':pro'
                ).data()[1,mask,:].T # [kHz], dim(nrho,t)

            # Obtain inversion errorbar data
            dout['profs'][line]['emis']['err'] = branchNode.getNode(
                'profiles.'+line+':proerr'
                ).data()[0,mask,:].T # dim(nrho,t)
            dout['profs'][line]['Ti']['err'] = branchNode.getNode(
                'profiles.'+line+':proerr'
                ).data()[3,mask,:].T # [keV], dim(nrho,t)
            dout['profs'][line]['vtor']['err'] = branchNode.getNode(
                'profiles.'+line+':proerr'
                ).data()[1,mask,:].T # [kHz], dim(nrho,t)

            # Radial domain
            dout['profs'][line]['psin'] = branchNode.getNode(
                'profiles.'+line+':rho'
                ).data()[mask,:].T # [], dim(nrho,t); rho = norm. pol. flux

            # Error check
            if dout['profs'][line]['Ti']['err'].shape[0] > dout['profs'][line]['psin'].shape[0]:
                ss = dout['profs'][line]['psin'].shape[0]
                dout['profs'][line]['emis']['data'] = dout['profs'][line]['emis']['data'][:ss,:]
                dout['profs'][line]['emis']['err'] = dout['profs'][line]['emis']['err'][:ss,:]
                dout['profs'][line]['Ti']['data'] = dout['profs'][line]['Ti']['data'][:ss,:]
                dout['profs'][line]['Ti']['err'] = dout['profs'][line]['Ti']['err'][:ss,:]
                dout['profs'][line]['vtor']['data'] = dout['profs'][line]['vtor']['data'][:ss,:]
                dout['profs'][line]['vtor']['err'] = dout['profs'][line]['vtor']['err'][:ss,:]

        except:
            print('No inversion data for line: '+line)

    # Output
    return dout

# Gets moments time traces v. psin
def _get_moments(
    dout=None,
    ):

    # Initializes output
    dout['mom'] = {}

    # Available data
    lines = ['X', 'Z', 'LYA1']

    # Loop over lines
    for line in lines:
        if line == 'X' or line == 'Z':
            bnch = 'HELIKE'
        elif line == 'LYA1':
            bnch = 'HLIKE'

        # MDSplus tree
        _, branchNode = _cmod_tree(dout=dout, branch=bnch)

        # MDSplus node
        mom_nd = branchNode.getNode('moments.'+line+':mom')

        # get time basis
        all_times = np.asarray(mom_nd.dim_of(1))
        mask = all_times > -1

        # Loads channel mapping
        chmap = branchNode.getNode('binning:chmap').data()[0]
        maxChan = np.max(chmap)
        ch_num = np.arange(maxChan) + 1

        data = mom_nd.data()[:,mask,:]
        psin = mom_nd.dim_of(0).data()[mask,:] # rho = norm. pol. flux

        # Stores data
        dout['mom'][line] = {}
        dout['mom'][line]['data'] = data[:,:,ch_num-1]          # [], dim(mom, nt, nchan)
        dout['mom'][line]['psin'] = psin[:,ch_num-1]            # [], dim(nt, nchan)
        dout['mom'][line]['t_s'] = all_times[mask]              # [s], dim(nt,)

    # Output
    return dout

# Gets intensity time traces
def _get_int(
    dout=None,
    ):

    # Initializes output
    dout['int'] = {}

    # Available data
    lines = ['A', 'M', 'W', 'Z']

    # MDSplus tree
    spec, _ = _cmod_tree(dout=dout)

    # Loop over available lines
    for ll in lines:
        # MDSplus node
        int_nd = spec.getNode(r'\spectroscopy::top.hirex_sr.analysis.'+ll+':int') 

        # Stores data
        dout['int'][ll] = {}
        dout['int'][ll]['data'] = int_nd.data()
        dout['int'][ll]['t_s'] = int_nd.dim_of(0).data()

        # Stores branch
        if ll == 'A' or ll == 'M':
            dout['int'][ll]['branch'] = 'HLIKE'
        elif ll == 'W' or ll == 'Z':
            dout['int'][ll]['branch'] = 'HELIKE'

    # Output
    return dout

# Loads averaged-spectra
def _get_spectra(
    dout=None,
    ):

    # Initializes output
    dout['spectra'] = {}

    # Available data
    branch = ['HELIKE', 'HLIKE']

    # Loop over modules
    for bb in branch:
        # MDSplus tree
        spec, branchNode = _cmod_tree(dout=dout, branch=bb)

        # MDSplus node
        br_nd = branchNode.getNode('spec:specbr')
        lam_nd = branchNode.getNode('spec:lam')

        # Initializes output
        dout['spectra'][bb] = {}

        # Stores time axis
        dout['spectra'][bb]['t_s'] = br_nd.dim_of(1).data() # dim(t,)

        # Stores averaged-spectra
        dout['spectra'][bb]['avespec'] = br_nd.data() # dim(hor_pxl, t, vert_ch)

        # Stores wavleength mesh
        dout['spectra'][bb]['lam_AA'] = lam_nd.data() # dim(hor_pxl, t, vert_ch)

        # Stores branch
        dout['spectra'][bb]['branch'] = bb

        # Gets POS vector for the spectra
        dout['spectra'][bb]['pos'] = _get_pos(
            dout=dout,
            branch=bb,
            avespec=True,
            ) # dim(hor_pxl, vert_ch, 4)

        # Removes padding
        numch = dout['spectra'][bb]['pos'].shape[1]
        dout['spectra'][bb]['avespec'] = dout['spectra'][bb]['avespec'][:,:,:numch]
        dout['spectra'][bb]['lam_AA'] = dout['spectra'][bb]['lam_AA'][:,:,:numch]

        ind_t = np.where(dout['spectra'][bb]['t_s']>0)

        dout['spectra'][bb]['avespec'] = dout['spectra'][bb]['avespec'][:,ind_t[0],:]
        dout['spectra'][bb]['lam_AA'] = dout['spectra'][bb]['lam_AA'][:,ind_t[0],:]
        dout['spectra'][bb]['t_s'] = dout['spectra'][bb]['t_s'][ind_t[0]]

    # Output
    return dout

# Obtains POS vectors
def _get_pos(
    dout=None,
    branch=None,
    avespec=None,
    ):
    '''
    Definition --
        POS = (R0, Z0, RT, Psi)
            R0 -- radial=coordinate of crystal center
            Z0 -- height-coordinate of crystal center
            RT -- tangency radial-coordinate of LOS 
            Psi -- angle in constant Z-plane of LOS
                tan(Psi) = ( Z0 -ZT ) / sqrt( R0**2 - RT**2 )

    '''

    # MDSplus tree
    spec, _ = _cmod_tree(dout=dout, branch=branch)

    # Obtains POS vector data
    if branch == 'HLIKE':
        pos = spec.getNode(
            r'\spectroscopy::top.hirexsr.calib.mod4:pos'
            ).data() # dim(hor_pxl, vert_pxl, 4)
    
    elif branch == 'HELIKE':
        pos = np.hstack([
            spec.getNode(
                r'\spectroscopy::top.hirexsr.calib.mod1:pos'
                ).data(), 
            spec.getNode(
                r'\spectroscopy::top.hirexsr.calib.mod2:pos'
                ).data(), 
            spec.getNode(
                r'\spectroscopy::top.hirexsr.calib.mod3:pos'
                ).data()
            ]) # dim(hor_pxl, vert_pxl, 4)

    # If binning vertical channels
    if avespec:
        pos = _cmod_pxl2bin(
            dout=dout,
            branch=branch,
            quant=pos,
            ) # dim(hor_pxl, vert_ch, 4)

    # Output
    return pos



###########################################
#
#               Extra
#
###########################################

# Gets MDSplus paths
def _cmod_tree(
    dout=None,
    branch=None,
    ):

    # Obtains data
    specTree = MDSplus.Tree('spectroscopy', dout['shot'])
    rootPath = r'\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'

    # Loads THACO branch
    if dout['tht'] > 0:
        rootPath += str(tht)
    rootPath += '.'

    # Analysis branch
    if branch is not None:
        branchNode = specTree.getNode(rootPath+branch)
    else:
        branchNode = None

    # Output
    return specTree, branchNode


# Gets map for spatial binning
def _cmod_pxl2bin(
    dout=None,
    branch=None,
    quant=None,
    ):

    # Obtains shot data
    _, branchNode = _cmod_tree(dout=dout, branch=branch)

    # Obtains mapping from pixels to bins
    # NOTE: axis=0 is just padding
    chmap = branchNode.getNode('BINNING:CHMAP').data() # dim(hor_pxl, vert_pxl * all modules)
    pxl_2_ch = chmap[0,:] # dim(vert_pxl * all modules)

    # Initializes binned quantity
    numBin = len(np.unique(pxl_2_ch[pxl_2_ch>0]))
    quant_bin = np.zeros((quant.shape[0], numBin) +quant.shape[2:]) # dim(hor_pxl, vert_ch, etc...)

    # Loop over bins
    for ii, chd in enumerate(np.unique(pxl_2_ch[pxl_2_ch>0])):
        # Averages quantity v. spatial pixel over the bin window
        quant_bin[:,ii] = np.mean(quant[:, pxl_2_ch == chd], axis = 1)

    # dim(hor_pxl, vert_ch, etc...)
    return quant_bin
