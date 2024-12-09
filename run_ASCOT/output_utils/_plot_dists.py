'''

Functions to plot distributions from ASCOT

cjperks
Sep 28, 2024

'''

# Modules
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cnt
from scipy.interpolate import griddata, interp1d
from matplotlib.ticker import ScalarFormatter
from matplotlib import gridspec

plt.rcParams.update({'font.size': 16})

formatter = ScalarFormatter()
formatter.set_powerlimits((-2, 2))  # Change the limits as needed

from omfit_classes import omfit_eqdsk


__all__ = [
    '_plot_rhoDist',
    '_plot_rzDist',
    '_plot_rzVDist'
    ]


#############################################
#
#               Main
#
#############################################

# Plots (vpara, vperp) distribution
def _plot_rzVDist(
    # Intial condition data
    idata = None,
    idf = None,
    # Output data
    ddata = None,
    df = None,
    # rho controls
    rhop_min = 0,
    rhop_max = 1,
    # Plot controls
    plt_t = -1,
    plt_all = True,
    # Extra
    gfile = None,
    ):

    # Gets data for trapped/passing boundary
    if gfile is not None:
        geq = omfit_eqdsk.OMFITgeqdsk(gfile)
        rmin_g = geq['fluxSurfaces']['geo']['a']
        rhop_g = np.sqrt(
            geq['fluxSurfaces']['geo']['psin']
            /geq['fluxSurfaces']['geo']['psin'][-1]
            )

        rmin = interp1d(
            rhop_g, rmin_g
            )((rhop_min+rhop_max)/2)

        # Critical angle
        cos_thetaC = np.sqrt(2*rmin/geq['RMAXIS']) # [v_||/v]
        thetaC = np.arccos(cos_thetaC)


    # Distribution data
    if idf is None:
        if idata is not None:
            idf = _process_dist(
                ddata=idata,
                key = 'rVDist',
                rhop_min = rhop_min,
                rhop_max = rhop_max
                )

    if df is None:
        df = _process_dist(
            ddata = ddata,
            key = 'rzVDist',
            rhop_min = rhop_min,
            rhop_max = rhop_max,
            )

    if not plt_all:
        return df

    #### --- PLot df/dVpe vs. Vperp for whoel plasma volume --- ####
    # Units [s/m^4]


    fig1, ax1  = plt.subplots()

    if idf is not None:
        '''
        ax1.plot(
            idf['whole']['Vpe_grid']['vals'],
            idf['whole'][0]['Vpe_dist']['vals'],
            '*-',
            label= r'init, $\rho_p$=[0,1]'
            )
        '''
        ax1.plot(
            idf['rhop']['Vpe_grid']['vals'],
            idf['rhop'][0]['Vpe_dist']['vals'],
            '*-',
            label= (
                r'init, $\rho_p$=[%0.2f, %0.2f]'%(rhop_min, rhop_max)
                )
            )

    for itime in np.arange(df['rhop']['t_grid']['size']):
        '''
        ax1.plot(
            df['whole']['Vpe_grid']['vals'],
            df['whole'][itime]['Vpe_dist']['vals'],
            '*-',
            label = r'$\rho_p$=[0,1]; $t_{RF,on}$ =%0.2d - %0.2d ms'%(
                df['whole']['t_grid']['vals'][itime], 
                df['whole']['t_grid']['vals'][itime+1]
                )
            )
        '''
        ax1.plot(
            df['rhop']['Vpe_grid']['vals'],
            df['rhop'][itime]['Vpe_dist']['vals'],
            '*-',
            label = r'$\rho_p$=[%0.2f,%0.2f]; $t_{RF,on}$ =%0.2d - %0.2d ms'%(
                rhop_min, rhop_max,
                df['rhop']['t_grid']['vals'][itime], 
                df['rhop']['t_grid']['vals'][itime+1]
                )
            )
        

    ax1.grid('on')
    ax1.set_yscale('log')
    #ax1.set_xscale('log')
    leg = ax1.legend(labelcolor='linecolor')
    leg.set_draggable('on')

    ax1.set_xlabel(r'$v_\perp$ [$m/s$]')
    ax1.set_ylabel(r'$f_{Ar16+}$ [$s/m^4$]')
 
    #ax1.set_ylim(1e7, 1e14)
    ax1.set_xlim(0, 1.5e4)
    ax1.xaxis.set_major_formatter(formatter)


    #### --- Plot log10(df/dVpa/dVpe) on (Vpara, Vperp) grid for given rho region --- ####
    # Units [s^2/m^5]

    fig2 = plt.figure(figsize=(16,10))
    nrows = int(np.ceil((df['rhop']['t_grid']['size']+1)/3))
    ncols = 3
    gs = gridspec.GridSpec(nrows, ncols, figure=fig2)

    # Plot initial distribution
    ax = fig2.add_subplot(gs[0,0])
    ax.set_ylabel(r'$v_{\perp}$ [$m/s$]')
    if idf is not None:
        vpaGrid_1d = idf['rhop']['Vpa_grid']['vals']
        vpeGrid_1d = idf['rhop']['Vpe_grid']['vals']
        vpeGrid_2d, vpaGrid_2d = np.meshgrid(vpeGrid_1d, vpaGrid_1d) # dim(nvpara, nvperp)

        tmp = np.log10(idf['rhop'][0]['VpaVpe_dist']['vals']) # dim(nvpara, nvperp)
        tmp[np.isinf(tmp)] = np.nan

        pcm = ax.pcolormesh(
            vpaGrid_2d,
            vpeGrid_2d,
            tmp
            )
        cbar = plt.colorbar(pcm, ax=ax)
        pcm.set_clim(-2,6.5)

        ax.set_title('init')
        #ax.set_yscale('log')
        #ax.set_xscale('log')
        ax.set_ylim(0,1.5e4)
        ax.set_xlim(-1.5e5, 1.5e5)
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_formatter(formatter)

        if gfile is not None:
            ax.plot(
                np.r_[0, 1.5e4/np.tan(thetaC)],
                np.r_[0,1.5e4],
                'k--'
                )
            ax.plot(
                np.r_[0, -1.5e4/np.tan(thetaC)],
                np.r_[0,1.5e4],
                'k--'
                )


    # Init
    key = 0
    vpaGrid_1d = df['rhop']['Vpa_grid']['vals']
    vpeGrid_1d = df['rhop']['Vpe_grid']['vals']
    vpeGrid_2d, vpaGrid_2d = np.meshgrid(vpeGrid_1d, vpaGrid_1d) # dim(nvpara, nvperp)

    # Loop over rows
    for rr in np.arange(nrows):
        # Loop over columns
        for cc in np.arange(ncols):
            if rr == 0 and cc == 0:
                continue
            if key == df['rhop']['t_grid']['size']:
                continue
            ax = fig2.add_subplot(gs[rr,cc])

            tmp = np.log10(df['rhop'][0]['VpaVpe_dist']['vals']) # dim(nvpara, nvperp)
            tmp[np.isinf(tmp)] = np.nan

            pcm = ax.pcolormesh(
                vpaGrid_2d,
                vpeGrid_2d,
                tmp
                )
            cbar = plt.colorbar(pcm, ax=ax)
            pcm.set_clim(-2,6.5)
            ax.grid('on')
            ax.yaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_formatter(formatter)

            ax.set_title(r'$t_{RF,on}$ =%0.2d - %0.2d ms'%(
                df['rhop']['t_grid']['vals'][key], 
                df['rhop']['t_grid']['vals'][key+1]
                ))
            #ax.set_yscale('log')
            #ax.set_xscale('log')
            ax.set_ylim(0,1.5e4)
            ax.set_xlim(-1.5e5, 1.5e5)

            if rr == nrows -1:
                ax.set_xlabel(r'$v_{||}$ [$m/s$]')
            if cc == ncols -1:
                cbar.set_label(r'$log_{10}(f_{Ar16+})$ [$log_{10}(s^2/m^5)$]')
            if cc == 0:
                ax.set_ylabel(r'$v_{\perp}$ [$m/s$]')

            if gfile is not None:
                ax.plot(
                    np.r_[0, 1.5e4/np.tan(thetaC)],
                    np.r_[0,1.5e4],
                    'k--'
                    )
                ax.plot(
                    np.r_[0, -1.5e4/np.tan(thetaC)],
                    np.r_[0,1.5e4],
                    'k--'
                    )

            # Increase counter
            key +=1

    fig2.suptitle(
        r'$\rho_p$=[%0.2f,%0.2f]'%(
                rhop_min, rhop_max
            )
        )   


    # Output
    return df


# Plots (pitch, E) distribution
def _plot_rzDist(
    # Intial condition data
    idata = None,
    idf = None,
    # Output data
    ddata = None,
    df = None,
    # rho controls
    rhop_min = 0,
    rhop_max = 1,
    # Plot controls
    plt_all = True,
    plt_whole = True,
    plt_t = -1,
    # Extra
    gfile = None,
    ):

    # Gets data for trapped/passing boundary
    if gfile is not None:
        geq = omfit_eqdsk.OMFITgeqdsk(gfile)
        rmin_g = geq['fluxSurfaces']['geo']['a']
        rhop_g = np.sqrt(
            geq['fluxSurfaces']['geo']['psin']
            /geq['fluxSurfaces']['geo']['psin'][-1]
            )

        rmin = interp1d(
            rhop_g, rmin_g
            )((rhop_min+rhop_max)/2)

        # Critical angle
        cos_thetaC = np.sqrt(2*rmin/geq['RMAXIS']) # [v_||/v]

        

    # Distribution data
    if idf is None:
        if idata is not None:
            idf = _process_dist(
                ddata=idata,
                key = 'rzDist',
                rhop_min = rhop_min,
                rhop_max = rhop_max
                )

    if df is None:
        df = _process_dist(
            ddata = ddata,
            key = 'rzDist',
            rhop_min = rhop_min,
            rhop_max = rhop_max,
            )

    if not plt_all:
        return df

    #### --- PLot df/dE vs. E for whoel plasma volume --- ####
    # Units [1/m3/eV]


    fig1, ax1  = plt.subplots()

    if idf is not None:
        if plt_whole:
            ax1.plot(
                idf['whole']['E_grid']['vals'],
                idf['whole'][0]['E_dist']['vals'],
                '*-',
                label= r'init, $\rho_p$=[0,1]'
                )
        ax1.plot(
            idf['rhop']['E_grid']['vals'],
            idf['rhop'][0]['E_dist']['vals'],
            '*-',
            label= (
                r'init, $\rho_p$=[%0.2f, %0.2f]'%(rhop_min, rhop_max)
                )
            )

    for itime in np.arange(df['rhop']['t_grid']['size']):
        label = r'$t_{RF,on}$ =%0.2d - %0.2d ms'%(
                df['whole']['t_grid']['vals'][itime], 
                df['whole']['t_grid']['vals'][itime+1]
                )

        if plt_whole:
            ax1.plot(
                df['whole']['E_grid']['vals'],
                df['whole'][itime]['E_dist']['vals'],
                '*-',
                label = r'$\rho_p$=[0,1]; '+label
                )

            label = r'$\rho_p$=[%0.2f,%0.2f]; '%(rhop_min, rhop_max) + label
        else:
            ax1.set_title(r'$\rho_p$=[%0.2f,%0.2f]; '%(rhop_min, rhop_max))

        ax1.plot(
            df['rhop']['E_grid']['vals'],
            df['rhop'][itime]['E_dist']['vals'],
            '*-',
            label = label
            )
        

    ax1.grid('on')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    leg = ax1.legend(labelcolor='linecolor')
    leg.set_draggable('on')

    ax1.set_xlabel('Energy [keV]')
    ax1.set_ylabel(r'$f_{Ar16+}$ [$1/m^3/eV$]')
 
    #ax1.set_ylim(1e7, 1e14)
    #ax1.set_xlim(1e-1, 1e2)
    #ax1.set_xlim(0,4e5)
    #ax1.xaxis.set_major_formatter(ScalarFormatter(useMathText=True)) 
    #ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    #### --- Plot log10(df/dpitch/dE) on (pitch, E) grid for given rho region --- ####
    # Units [1/m3/eV]

    fig2 = plt.figure(figsize=(14,8))
    nrows = int(np.ceil((df['rhop']['t_grid']['size']+1)/3))
    ncols = 3
    gs = gridspec.GridSpec(nrows, ncols, figure=fig2)

    # Plot initial distribution
    ax = fig2.add_subplot(gs[0,0])
    ax.set_ylabel('Energy [keV]')
    if idf is not None:
        pGrid_1d = idf['rhop']['pitch_grid']['vals']
        eGrid_1d = idf['rhop']['E_grid']['vals']
        eGrid_2d, pGrid_2d = np.meshgrid(eGrid_1d, pGrid_1d) # dim(npitch, nE)

        tmp = np.log10(idf['rhop'][0]['pitchE_dist']['vals']) # dim(npitch, nE)
        tmp[np.isinf(tmp)] = np.nan

        pcm = ax.pcolormesh(
            pGrid_2d,
            eGrid_2d,
            tmp
            )
        cbar = plt.colorbar(pcm, ax=ax)
        pcm.set_clim(1,13)

        ax.set_title('init')
        ax.set_yscale('log')
        ax.set_ylim(1e-1,1e3)

        if gfile is not None:
            ax.plot(
                np.r_[cos_thetaC, cos_thetaC],
                np.r_[1e-1,1e3],
                'k--'
                )
            ax.plot(
                np.r_[cos_thetaC, cos_thetaC]*-1,
                np.r_[1e-1,1e3],
                'k--'
                )



    # Init
    key = 0
    pGrid_1d = df['rhop']['pitch_grid']['vals']
    eGrid_1d = df['rhop']['E_grid']['vals']
    eGrid_2d, pGrid_2d = np.meshgrid(eGrid_1d, pGrid_1d) # dim(npitch, nE)

    # Loop over rows
    for rr in np.arange(nrows):
        # Loop over columns
        for cc in np.arange(ncols):
            if rr == 0 and cc == 0:
                continue
            if key == df['rhop']['t_grid']['size']:
                continue
            ax = fig2.add_subplot(gs[rr,cc])

            tmp = np.log10(df['rhop'][key]['pitchE_dist']['vals']) # dim(npitch, nE)
            tmp[np.isinf(tmp)] = np.nan

            pcm = ax.pcolormesh(
                pGrid_2d,
                eGrid_2d,
                tmp
                )
            cbar = plt.colorbar(pcm, ax=ax)
            pcm.set_clim(1,13)
            ax.grid('on')
            #ax.yaxis.set_major_formatter(formatter)

            ax.set_title(r'$t_{RF,on}$ =%0.2d - %0.2d ms'%(
                df['rhop']['t_grid']['vals'][key], 
                df['rhop']['t_grid']['vals'][key+1]
                ))
            ax.set_yscale('log')
            ax.set_ylim(1e-1,1e3)

            if rr == nrows -1:
                ax.set_xlabel(r'pitch [$v_{||}/v$]')
            if cc == ncols -1:
                cbar.set_label(r'$log_{10}(f_{Ar16+})$ [$log_{10}(1/m^3/eV)$]')
            if cc == 0:
                ax.set_ylabel('Energy [keV]')

            if gfile is not None:
                ax.plot(
                    np.r_[cos_thetaC, cos_thetaC],
                    np.r_[1e-1,1e3],
                    'k--'
                    )
                ax.plot(
                    np.r_[cos_thetaC, cos_thetaC]*-1,
                    np.r_[1e-1,1e3],
                    'k--'
                    )

            # Increase counter
            key +=1

    fig2.suptitle(
        r'$\rho_p$=[%0.2f,%0.2f]'%(
                rhop_min, rhop_max
            )
        )   


    # Output
    return df

# Plot 1D profile data
def _plot_rhoDist(
    ind_dist = 1,
    ddata = None,
    ylabel = None,
    # Extra
    dext = None, # External data
    shot = None,
    ):

    # Init
    tmp = ddata['rhoDist']['dists'][ind_dist]
    rhoD = ddata['rhoDist']

    fig, ax = plt.subplots()

    if dext is not None:
        ax.plot(
            dext['rhop'],
            dext['vals'],
            'k-',
            label = dext['label'],
            linewidth = 2
            )

    for itime in np.arange(rhoD['cents']['dim2']['size']):
        ax.plot(
            rhoD['cents']['dim1']['vals'],
            tmp['vals'][:,itime,0]/rhoD['cents']['dim2']['dbin'],
            label = r'$t_{RF,on}$ =%0.2d - %0.2d ms'%(
                rhoD['edges']['dim2']['vals'][itime]*1e3, 
                rhoD['edges']['dim2']['vals'][itime+1]*1e3
                )
        )

    ax.set_xlabel(r'$\sqrt{\Psi_n}$')
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel(r'[$'+tmp['units'].decode("utf-8")+'$]')
    
    leg = ax.legend(labelcolor='linecolor')
    leg.set_draggable('on')
    ax.grid('on')
    ax.set_title(shot)
    ax.yaxis.set_major_formatter(formatter)


#############################################
#
#               Utilities
#
#############################################

# Processes (R,Z) distribution data
def _process_dist(
    ddata = None,
    key = None,
    rhop_min = None,
    rhop_max = None,
    ):
    print(rhop_min)
    print(rhop_max)

    # Init
    rzD = ddata[key]
    df = {}
    df['whole'] = {} # Sum over whole plasma volume
    df['rhop'] = {} # Sum over given rho interval

    if key == 'rzDist':
        key_2d = 'pitchE_dist'
        units_2d = r'$1/m^3/eV$'
        key_1d = 'E_dist'
        units_1d = units_2d
        scalef = cnt.e

    elif key == 'rzVDist':
        key_2d = 'VpaVpe_dist'
        units_2d = r'$s^2/m^5$'
        key_1d = 'Vpe_dist'
        units_1d = r'$s/m^4$'
        scalef = 1

    # Total distribution data
    df['tot'] = {}
    df['tot']['vals'] = rzD['dist']['vals'] *scalef # dim(nR, nZ, npitch, nE, nt), [1/m3/eV]
    df['tot']['units'] = units_2d

    # If we need time-averaging
    if rzD['cents']['dim5']['size'] > 1:
        df['tot']['vals'] /= rzD['cents']['dim5']['dbin']

    # Finds rhop values of each data point
    RR_bf_1d = ddata['bfield']['RR']
    ZZ_bf_1d = ddata['bfield']['ZZ']
    ZZ_bf_2d, RR_bf_2d = np.meshgrid(ZZ_bf_1d, RR_bf_1d)

    RR_rz_1d = rzD['cents']['dim1']['vals'] # dim(nR,)
    ZZ_rz_1d = rzD['cents']['dim2']['vals'] # dim(nZ, )
    ZZ_rz_2d, RR_rz_2d = np.meshgrid(ZZ_rz_1d,RR_rz_1d) # dim(nR, nZ)
    rhop_vals = griddata(
        (RR_bf_2d.flatten(), ZZ_bf_2d.flatten()),
        ddata['bfield']['rhop'].T.flatten(),
        (RR_rz_2d.flatten(), ZZ_rz_2d.flatten()),
        method = 'linear'
        ).reshape(RR_rz_2d.shape) # dim(nR, nZ)

    # Loop over time
    for itime in np.arange(df['tot']['vals'].shape[-1]):
        # Dictionary for this time slice
        df['whole'][itime] = {}
        df['rhop'][itime] = {}

        # Time slice data
        df['whole'][itime]['time'] = {}
        df['rhop'][itime]['time'] = {}
        df['whole'][itime]['time']['vals'] = df['rhop'][itime]['time']['vals'] = [
            rzD['edges']['dim5']['vals'][itime],
            rzD['edges']['dim5']['vals'][itime+1]
            ]
        df['whole'][itime]['time']['dbin'] = df['rhop'][itime]['time']['dbin'] =  rzD['edges']['dim5']['dbin']

        # 4D distribution at this time
        df['whole'][itime]['4D_dist'] = {}
        df['rhop'][itime]['4D_dist'] = {}
        df['whole'][itime]['4D_dist']['vals'] = df['rhop'][itime]['4D_dist']['vals'] = (
            df['tot']['vals'][...,itime]
            ) # dim(nR, nZ, npitch, nE), [1/m3/eV]
        df['whole'][itime]['4D_dist']['units'] = df['rhop'][itime]['4D_dist']['units'] = df['tot']['units']


        # Integrates 4D dist over 
        # ... whole plasma volume
        df['whole'][itime][key_2d] = {}
        df['whole'][itime][key_2d]['units'] = units_2d
        df['whole'][itime][key_2d]['vals'] = np.zeros((
            rzD['cents']['dim3']['size'], rzD['cents']['dim4']['size']
            )) # dim(npitch, nE); [1/m3/eV]

        # ... rho interval
        df['rhop'][itime][key_2d] = {}
        df['rhop'][itime][key_2d]['units'] = units_2d
        df['rhop'][itime][key_2d]['vals'] = np.zeros((
            rzD['cents']['dim3']['size'], rzD['cents']['dim4']['size']
            )) # dim(npitch, nE); [1/m3/eV]

        # Running counters
        vol_tot = 0
        vol_rhop = 0

        for rr in np.arange(rzD['cents']['dim1']['size']):
            for zz in np.arange(rzD['cents']['dim2']['size']):

                # Annulus volume
                voxel = (
                    np.pi *(
                        rzD['edges']['dim1']['vals'][rr+1]**2
                        - rzD['edges']['dim1']['vals'][rr]**2
                        )
                    * (
                        rzD['edges']['dim2']['vals'][zz+1]
                        -rzD['edges']['dim2']['vals'][zz]
                        )
                    ) # [m3]
                vol_tot += voxel

                # Intergrates over plasma volume
                df['whole'][itime][key_2d]['vals'] += (
                    df['tot']['vals'][rr,zz,:,:,itime]*voxel
                    ) # [1/eV]

                # If volume is in rhop interval
                if rhop_min <= rhop_vals[rr,zz] <= rhop_max:
                    vol_rhop += voxel
                    df['rhop'][itime][key_2d]['vals'] += (
                        df['tot']['vals'][rr,zz,:,:,itime]*voxel
                        ) # [1/eV]

                

        # Averages over volume
        df['whole'][itime][key_2d]['vals'] /= vol_tot
        df['rhop'][itime][key_2d]['vals'] /= vol_rhop
        df['whole'][itime][key_2d]['units'] = df['rhop'][itime][key_2d]['units'] = df['tot']['units']
        
        # 4D dist averaged over whole plasma volume and integrated over pitch
        df['whole'][itime][key_1d] = {}
        df['whole'][itime][key_1d]['vals'] = np.sum(
            df['whole'][itime][key_2d]['vals'],
            axis = 0
            )*rzD['edges']['dim3']['dbin'] # dim(nE); [1/m3/eV]

        df['rhop'][itime][key_1d] = {}
        df['rhop'][itime][key_1d]['vals'] = np.sum(
            df['rhop'][itime][key_2d]['vals'],
            axis = 0
            )*rzD['edges']['dim3']['dbin'] # dim(nE); [1/m3/eV]
        df['whole'][itime][key_1d]['units'] = df['rhop'][itime][key_1d]['units'] = units_1d

        # Store rhop grid
        df['rhop']['grid'] = {}
        df['rhop']['grid']['RR'] = RR_rz_2d
        df['rhop']['grid']['ZZ'] = ZZ_rz_2d
        df['rhop']['grid']['rhop'] = rhop_vals


        if key == 'rzDist':
            # Store energy grid
            df['whole']['E_grid'] = {}
            df['rhop']['E_grid'] = {}
            df['whole']['E_grid']['vals'] = df['rhop']['E_grid']['vals'] = rzD['cents']['dim4']['vals']/cnt.e/1e3 # dim(nE,), [keV]
            df['whole']['E_grid']['units'] = df['rhop']['E_grid']['units'] = 'keV'

            # Store pitch grid
            df['whole']['pitch_grid'] = {}
            df['rhop']['pitch_grid'] = {}
            df['whole']['pitch_grid']['vals'] = df['rhop']['pitch_grid']['vals'] = rzD['cents']['dim3']['vals'] # dim(npitch,), []
            df['whole']['pitch_grid']['units'] = df['rhop']['pitch_grid']['units'] = 'v_{||}/v'
        elif key == 'rzVDist':
            # Store parallel energy grid
            df['whole']['Vpa_grid'] = {}
            df['rhop']['Vpa_grid'] = {}
            df['whole']['Vpa_grid']['vals'] = df['rhop']['Vpa_grid']['vals'] = rzD['cents']['dim3']['vals'] # dim(nVpa,), [m/s]
            df['whole']['Vpa_grid']['units'] = df['rhop']['Vpa_grid']['units'] = 'm/s'

            # Store perpendicular energy grid
            df['whole']['Vpe_grid'] = {}
            df['rhop']['Vpe_grid'] = {}
            df['whole']['Vpe_grid']['vals'] = df['rhop']['Vpe_grid']['vals'] = np.flipud(rzD['cents']['dim4']['vals']) # dim(nVpe,), [m/s]
            df['whole']['Vpe_grid']['units'] = df['rhop']['Vpe_grid']['units'] = 'm/s'

        # Store time grid
        df['whole']['t_grid'] = {}
        df['rhop']['t_grid'] = {}
        df['whole']['t_grid']['vals'] = df['rhop']['t_grid']['vals'] = rzD['edges']['dim5']['vals']*1e3 # dim(nt,), [ms]
        df['whole']['t_grid']['units'] = df['rhop']['t_grid']['units'] = 'ms'
        df['whole']['t_grid']['size'] = df['rhop']['t_grid']['size'] = rzD['cents']['dim5']['size']


    # Output
    return df
