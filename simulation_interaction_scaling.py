from gpe import HarmonicTrapGPE, GPE

# Load GPE
sim = HarmonicTrapGPE("U_scaling", 1.0, 1.0, 1.0)

# Grid parameters
sim.setGridSize(128)
sim.set("maxX", 10)
sim.set("maxY", 10)
sim.set("maxZ", 10)


# Contact
sim.set("scatteringLength", 1)
sim.set("atomNumber",       1) # start non-interacting

sim.set("dipolarLength", 0)
sim.setBool("evolutionDipolar", False)

sim.initialize()

#sim.ite(1000)

sim.set("scatteringLength", 1)
print "0 0"

#dataf = file("U_scaling/int_energy", "w")

#for n in [2, 3, 4, 5, 10, 20, 50, 100, 150, 200, 300, 400, 500]:
for n in [100]:
    sim.set("atomNumber", n)
    sim.ite(10000, 10)
    sim.write("n_{n}".format(n = n))
    #dataf.write("{n} {e}\n".format(n = n, e = sim.energyInt()))
    print "{n} {e}\n".format(n = n, e = sim.energyInt())


