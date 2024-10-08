'''

Defaults settings for various CMOD shots


'''

# Modules
import numpy as np

__all__ = [
    'def_dshot'
]

########################################################
#
#           Main
#
########################################################

def def_dshot(
    shot = None,
    t0 = None,
    ):


    dshot = {
        1140221013:{
            't1':{ # [850, 950] ms
                'H2D': 0.51,
                'D_con': 0.6,
                'Ar_con': 1e-3,
                'Mo_con': 7e-4,
                'Prf_in': 0,
                'dt': [0.85, 0.95],
                'fimprad': None,
                },
            't2':{ # [1150, 1300] ms
                'H2D' : 0.51,
                'D_con': 0.6,
                'Ar_con': 0.5e-3,
                'Mo_con': 10e-4,
                'Prf_in': 0.5,
                'dt': [1.15, 1.3],
                'fimprad': None,
            },
        },

        1140221012:{
            't1':{ # [860, 980] ms
                'H2D': 0.51,
                'D_con': 0.6,
                'Ar_con': 1e-3,
                'Mo_con': 7e-4,
                'Prf_in': 0,
                'dt': [0.86, 0.98],
                'fimprad': 'imprad_1140221012.npy',
                },
            't2':{ # [1170, 1350] ms
                'H2D' : 0.51,
                'D_con': 0.6,
                'Ar_con': 0.5e-3,
                'Mo_con': 10e-4,
                'Prf_in': 0.5,
                'dt': [1.17, 1.35],
                'fimprad': None,
                },
            },

        1140408020:{
            't1':{
                'H2D': 0.05,
                'D_con': 0.6,
                'Ar_con': 1e-3,
                'Mo_con': 7e-4,
                'Prf_in': 0,
                'dt': 0.01,
                'fimprad': 'imprad_1140408020.npy',
                },
            't2':{
                'H2D' : 0.05,
                'D_con': 0.6,
                'Ar_con': 0.5e-3,
                'Mo_con': 10e-4,
                'Prf_in': 0.5,
                'dt': 0.05,
                'fimprad': None,
                },
            },
        }

    # Ion concentration modeling
    dmodel = {
        'option': 'concentration',
        'ions': {
            'D': {
                'con': dshot[shot][t0]['D_con'],
                },
            'H': {
                'con': dshot[shot][t0]['H2D']
                },
            'Ar': {
                'con': dshot[shot][t0]['Ar_con'],
                },
            'Mo': {
                'con': dshot[shot][t0]['Mo_con'],
                },
            },
        'rotation':{
            'option': 'zero',
            },
        }
    dmodel['ions']['H']['con'] *= dmodel['ions']['D']['con']


    # Output
    return dmodel, dshot[shot][t0]['dt'], dshot[shot][t0]['Prf_in'], dshot[shot][t0]['fimprad']


