from gpe import HarmonicTrapGPE, GPE

# Load GPE
sim = HarmonicTrapGPE("Simple", 1.0, 1.0, 1.0)

# Grid parameters
sim.setGridSize(32)

# Contact interaction
# Only the product (atomNumber - 1) * scatteringLength is important
sim.set("atomNumber",       2)
sim.set("scatteringLength", 1)

# Dipolar interaction (switched off)
sim.set("dipolarLength", 0)
sim.set("evolutionDipolar", False)

# Number of threads. Equals number of cores if not set
# sim.setThreads(4)

sim.initialize()

# Imaginary time evolution
sim.ite(steps=10000, monitorSteps=500, plotSteps=500)
