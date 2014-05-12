from gpe import HarmonicTrapGPE, GPE

# Load GPE
sim = HarmonicTrapGPE("Vortices", 1.0, 1.0, 9.0)

# Grid parameters
sim.setGridSize(64, 64, 16)

# Contact and dipolar interaction
sim.set("atomNumber",       2)
sim.set("scatteringLength", 1)
sim.set("dipolarLength",    0)

#sim.setEvaluation({"dispX" : GPE.DISPERSION_X,
#                   "dispZ" : GPE.DISPERSION_Z,
#                   "ratio" : lambda sim: sim.evaluate(GPE.DISPERSION_X) / sim.evaluate(GPE.DISPERSION_Z),
#                   "add" : "dipolarLength",
#                   "time" : lambda sim: sim.time
#                   })

#sim.setThreads(10)

sim.initialize()

#sim.setEvaluation([GPE.DISPERSION_X, GPE.DISPERSION_Y, GPE.DISPERSION_Z, GPE.PSI_POWER_2, GPE.PSI_POWER_6])

#sim.ite(1000, 100)
sim.applyRotation(1.0)

sim.rte(10000, 500, 500)
