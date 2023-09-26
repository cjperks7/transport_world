'''

_plt_slider.py is a multi-purpose tool to plot
2D lines with a slider for a third dimension 

'''

import aurora
import numpy as np

__all__ = [
    'plt_slider'
    ]

###############################################
#
#               Main
#
###############################################

def plt_slider(
    xxx=None,
    yyy=None,
    dzzz=None,
    xxlabel=None,
    yylabel=None,
    zzlabel=None,
    dzlabels=None,
    plt_sum=None,
    x_line=None,
    y_line=None,
    xscale='linear',
    yscale='linear',
    xlim=None,
    ylim=None,
    ):
    '''
    NOTE: I assume my data is in the form of (yyy,xxx)
        i.e., spectra in the shape (rhop, lambda)
                and I want to slide over rhop

        so whatever you want to slide over is first dimension

    '''

    # Initializes
    zzz = None
    labels = []

    # Reorganizes input data
    if isinstance(dzzz, dict):
        for key1 in dzzz.keys():
            if isinstance(dzzz[key1], dict):
                for key2 in dzzz[key1].keys():
                    if zzz is None:
                        zzz = (dzzz[key1][key2].T)[None,:,:]
                    else:
                        zzz = np.concatenate(
                            (zzz, (dzzz[key1][key2].T)[None,:,:]),
                            axis = 0
                            )

                    labels.append(
                        str(key1)
                        + ' +'
                        + str(key2)
                        )

            else:
                if zzz is None:
                    zzz = (dzzz[key1].T)[None,:,:]

                else:
                    zzz = np.concatenate(
                        (zzz, (dzzz[key1].T)[None,:,:]),
                        axis = 0
                        )

                labels.append(str(key1))

    elif dzzz.ndim == 2:
        zzz = (dzzz.T)[None,:,:] # di(1, xxx, yyy)

        labels = ['total']

    else:
        zzz = dzzz.transpose(0,2,1)  # dim(nlabels, xxx, yyy)
        labels = dzlabels

    # Uses Aurora slider tools
    aurora.plot_tools.slider_plot(
        xxx,
        yyy,
        zzz,
        xlabel=xxlabel,
        ylabel=yylabel,
        zlabel=zzlabel,
        labels=labels,
        plot_sum=plt_sum,
        x_line=x_line,
        y_line=y_line,
        xscale=xscale,
        yscale=yscale,
        xlim=xlim,
        ylim=ylim,
        )


