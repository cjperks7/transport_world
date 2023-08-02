'''

plot_kin.py plots time traces of kinetic profs

'''

# Module
from transport_world.make_gacode import read_tokamak as rTok
from transport_world.plot_utils import plt_slider

shot = 1140221013

# Obtains experimental kinetic profiles
dkin = rTok.profs_cmod(shot=shot)

kin = 'ne'
diag = 'Thomson'
quant = 'val_1e20m3'

t0 = 1.2
t_ind = np.argmin(abs(dkin[kin][diag]['t_s'] - t0))

plt_slider(
    xxx=dkin[kin][diag]['r_m'][:,t_ind],
    yyy=dkin[kin][diag]['t_s'],
    dzzz=dkin[kin][diag][quant].T,
    xxlabel=r"$r_V$ [$m$]",
    yylabel=r"$t$ [$s$]",
    zzlabel=r"",
    plt_sum=False
    )