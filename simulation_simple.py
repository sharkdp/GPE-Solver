from gpe import HarmonicTrapGPE, GPE

# Load GPE
sim = HarmonicTrapGPE("Test", 1.0, 1.0, 1.0)

# Grid parameters
sim.setGridSize(32)

# Contact and dipolar interaction
sim.set("atomNumber",       2)
sim.set("scatteringLength", 1)
sim.set("dipolarLength",    1)

#sim.setEvaluation({"dispX" : GPE.DISPERSION_X,
#                   "dispZ" : GPE.DISPERSION_Z,
#                   "ratio" : lambda sim: sim.evaluate(GPE.DISPERSION_X) / sim.evaluate(GPE.DISPERSION_Z),
#                   "add" : "dipolarLength",
#                   "time" : lambda sim: sim.time
#                   })

sim.setThreads(6)

sim.initialize()

sim.setEvaluation([GPE.DISPERSION_X, GPE.DISPERSION_Y, GPE.DISPERSION_Z, GPE.PSI_POWER_2, GPE.PSI_POWER_6])

sim.ite(1000, 100)

sim.set("lossRate", 5)
sim.setBool("evolutionThreeBodyLosses", True)

sim.set("modulationAmplitude", 0.1)
sim.set("modulationFreq", 0.1)
sim.setString("externalPotential", "modulatedHarmonic")

sim.rte(10000, 100, 100)
