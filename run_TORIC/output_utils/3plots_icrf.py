#!/usr/bin/env python
"""
Make the three most common plots for toric and exit
"""

from toric_tools import *
import matplotlib.pyplot as plt
plt.ioff()
Run=toric_analysis("toric.ncdf",mode='ICRF')
Run.threeplots()
exit

