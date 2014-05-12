#!/usr/bin/python

from gpe import LatticeGPE, GPE

import os

REPULSIVE = 0
ATTRACTIVE = 1
text = ["repulsive", "attractive"]

# Settings
geometry = REPULSIVE
sites = 3

ite_steps_init = 2000
energy_steps = 200
ite_steps = 1000

depth = 80

a = 1.0
add = 3.0

# End Settings

if geometry == REPULSIVE:
    direction = 0
    width = [0.5, 4, 0.5]
else: # ATTRACTIVE
    direction = 2
    width = [0.5, 4, 0.5]

text_dir = ["x", "y", "z"]

title = "single_{geometry}_{add}".format(add=add, geometry = text[geometry])
print title
sim = LatticeGPE(title, sites, direction, depth, width)

# Grid parameters
sizeLat = 64
sizePerp = 64
if geometry == REPULSIVE:
    sim.setGridSize(sizeLat, sizePerp, sizePerp)
    sim.set("maxX", 2.5)
    sim.set("maxY", 2.5)
    sim.set("maxZ", 1.5)
else:
    sim.setGridSize(sizePerp, sizePerp, sizeLat)
    sim.set("maxX", 1.5)
    sim.set("maxY", 1.5)
    sim.set("maxZ", 2.5)

# Contact and dipolar interaction
sim.set("atomNumber",       2)
sim.set("scatteringLength", a)
sim.set("dipolarLength",    add)

sim.setBool("evolutionContact", False)
sim.setBool("evolutionDipolar", False)

sim.initialize()
sim.write("initial")

sim.ite(ite_steps_init, energy_steps)
sim.write("non_interacting")

sim.setBool("evolutionContact", True)
sim.ite(ite_steps_init, energy_steps)
sim.write("contact_interacting")

sim.setBool("evolutionDipolar", True)
sim.ite(2000, energy_steps)

i = 0
while a > 0.7:
    sim.set("scatteringLength", a)
    sim.ite(100)
    sim.write("state_{i}".format(i=i))
    a = a - 0.01
    i = i + 1
