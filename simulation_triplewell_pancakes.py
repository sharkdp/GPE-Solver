#!/usr/bin/python

from gpe import LatticeGPE, GPE

sites = 3
direction = 0

ite_steps = 10000
energy_steps = 100

depth = 80


width = [0.5, 0.5, 0.35]

#for a in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]:
for a in [0.1]:
  for add in [0.52]:
    # Load GPE
    #sim = LatticeGPE("TripleWell_June3_{a}_{add}_{width}_{depth}".format(a=a, add=add, width=width, depth=depth), sites, direction, depth, width)
    title = "TripleWell_Pancake_{a}_{add}".format(w = width[0], w2 = width[2], a=a, add=add)
    print title
    sim = LatticeGPE(title, sites, direction, depth, width)

    # Grid parameters
#    sim.setGridSize(256, 128, 128)
    sim.setGridSize(128, 64, 64)

    # Contact and dipolar interaction
    sim.set("atomNumber",       2)

    sim.set("scatteringLength", a)
    sim.set("dipolarLength",    add)

    sim.setBool("evolutionContact", False)
    sim.setBool("evolutionDipolar", False)

    sim.initialize()

    sim.writePotential("potential")
    sim.write("initial")

    sim.ite(ite_steps, energy_steps)
    sim.write("non_interacting")

    sim.setBool("evolutionContact", True)
    sim.ite(ite_steps, energy_steps)
    sim.write("contact_interacting")

    sim.setBool("evolutionDipolar", True)
    sim.ite(ite_steps, energy_steps)
    sim.write("dipole_interacting")

    sim.setBool("evolutionContact", False)
    sim.setBool("evolutionDipolar", False)
    sim.ite(ite_steps, energy_steps)
    sim.write("non_interacting2")

    del sim
