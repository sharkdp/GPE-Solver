#!/usr/bin/python

from gpe import LatticeGPE, GPE

sites = 3
direction = 0

depth = 80

width = [0.5, 0.5, 0.35]

for a in [0.1, 1, 2]:
    for add in [0, 0.05, 0.1, 0.5, 1.0]:
        if add > a:
            break
        title = "TripleWell_Dynamics_{a}_{add}".format(a = a, add = add)
        sim = LatticeGPE(title, sites, direction, depth, width)

        # Grid parameters
        #    sim.setGridSize(256, 128, 128)
        sim.setGridSize(128, 64, 64)
        #    sim.setGridSize(128, 128, 128)

        # Contact and dipolar interaction
        sim.set("atomNumber",       2)

        sim.set("scatteringLength", a)
        sim.set("dipolarLength",    add)

#        sim.setBool("evolutionContact", False)
#        sim.setBool("evolutionDipolar", False)

        # place in "left" well
        sim.set("offsetX", -1.0)
        sim.set("sigmaX", 0.25)

        # only the "left" well in the beginning
        sim.set("latticeSites", 1)
        sim.set("latticeOffsetX", -1.0)

        sim.initialize()
        sim.write("initial")
        sim.ite(5000, 100)
        sim.write("groundstate")

        # enable remaining wells
        sim.set("latticeSites", 3)
        sim.set("latticeOffsetX", 0.0)

        sim.initialize(False)

        sim.rte(50000, 500, 500)
        del sim
