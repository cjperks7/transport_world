


from itertools import islice


filename = 'test_ascot_31289655.orbits.test'
skip_rows = 24
x = int(1e7)

orbits = {}

quants = [
    'energy',
    'pitch',
    'vmagn',
    'vparal',
    'vperp',
    'vphi',
    'phi',
    'R',
    'z',
    'id',
    'rho',
    'time',
    'endcond'
    ]

with open(filename, 'r') as file:
    for line in islice(file, skip_rows, skip_rows +x):
        data = line.split()

        if len(data) <= 1:
            continue
        mid = data[9]

        if mid not in orbits.keys():
            orbits[mid] = {}
            
            for qq in quants:
                orbits[mid][qq] = []

        for ii, qq in enumerate(quants):
            orbits[mid][qq].append(float(data[ii]))


key = '59'

import matplotlib.pyplot as plt
fig, ax = plt.subplots()

ax.plot(
    orbits[key]['R'],
    orbits[key]['z']
)

ax.set_xlim(0.4, 0.9)
ax.set_ylim(-0.4, 0.4)