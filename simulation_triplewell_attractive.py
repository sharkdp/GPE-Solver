#!/usr/bin/python

from gpe import LatticeGPE, GPE

sites = 3
direction = 2

ite_steps = 15000
energy_steps = 100

depth = 80

width = [30, 30, 0.5]

#for a in [0.01, 0.1, 1, 2, 4, 8]:
for a in [1, 2, 4, 8]:
    d_add = a / 10
    add = -a
    while add <= a:
        title = "TripleWell_attractive_{a}_{add}".format(a=a, add=add)
        print title
        sim = LatticeGPE(title, sites, direction, depth, width)

        #TODO
        sim.set("maxX", 7)
        sim.set("maxY", 7)

        # Grid parameters
        #sim.setGridSize(128, 128, 256)
        sim.setGridSize(64, 64, 128)

        # Contact and dipolar interaction
        sim.set("atomNumber",       2)

        sim.set("scatteringLength", a)
        sim.set("dipolarLength",    add)

        sim.setBool("evolutionContact", False)
        sim.setBool("evolutionDipolar", False)

        sim.initialize()
        sim.write("initial")

        sim.ite(ite_steps, energy_steps)
        sim.write("non_interacting")

        sim.setBool("evolutionContact", True)
        #sim.ite(ite_steps, energy_steps)
        #sim.write("contact_interacting")

        sim.setBool("evolutionDipolar", True)
        sim.ite(ite_steps, energy_steps)
        sim.write("dipole_interacting")

        sim.setBool("evolutionContact", False)
        sim.setBool("evolutionDipolar", False)
        sim.ite(ite_steps, energy_steps)
        sim.write("non_interacting2")

        del sim

        add = add + d_add
