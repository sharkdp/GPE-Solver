from gpe import LatticeGPE

dipolar = False
contact = False

name = "DoubleWell"
if dipolar:
    name = name + "_dipolar"

sim = LatticeGPE(name, 2, 2, 50, 0.5)

# Grid parameters
#sim.setGridSize(128, 128, 128)
#sim.setGridSize(64, 64, 64)
sim.setGridSize(128, 256, 256)

# Contact and dipolar interaction
sim.set("atomNumber",       2)
sim.set("scatteringLength", 2)
sim.set("dipolarLength",    1.5)

sim.setBool("evolutionDipolar", dipolar)
sim.setBool("evolutionContact", contact)

sim.initialize()

sim.ite(3000, 500, 1000)
sim.setBool("evolutionPotential", False)
sim.set("timeStep", 0.00025)
sim.rte(1600, 100, 100)
