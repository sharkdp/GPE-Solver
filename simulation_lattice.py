from gpe import LatticeGPE, GPE

sites = 3
direction = 0
depth = 50
width = 0.5

a = 10
add = 10

# Load GPE
sim = LatticeGPE("Lattice", sites, direction, depth, width)

# Grid parameters
sim.setGridSize(128, 32, 32)

# Contact and dipolar interaction
sim.set("atomNumber",       2)
sim.set("scatteringLength", a)
sim.set("dipolarLength",    add)

# load only one well
#sim.set("sigmaX", 0.2)
#sim.set("sigmaY", 0.2)
#sim.set("sigmaZ", 0.2)
#sim.set("maxX", 3)
#sim.set("maxY", 3)
#sim.set("maxZ", 3)

# pancake shaped
#sim.set("latticeWidthZ", 0.2)
#sim.set("sigmaZ", 0.2)
#sim.set("maxZ", 0.6)

sim.initialize()

sim.ite(5000, 50, 500)

#for i in range(0, 50):
#    sim.ite(500, 50, 500)
#    sim.set("latticeDepth", 50 + i * 50)
#    sim.initialize(False)

#sim.set("latticeSites", 3)
#sim.initialize(False)
#sim.rte(400000, 50, 500)

