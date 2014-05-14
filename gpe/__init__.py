import math
import os
from gpe.gpepython import GPEPython #@UnresolvedImport


class GPE(GPEPython):
    """Wrapper class to control the GPE program more easily"""

    # evaluation quantities
    DISPERSION_X = 1
    DISPERSION_Y = 2
    DISPERSION_Z = 3
    PSI_POWER_2  = 4
    PSI_POWER_6  = 5
    PSI_0RE      = 6
    PSI_0IM      = 7

    def __init__(self, name):
        """Initialize the GPE Simulation object"""

        # initialize the C++ GPEPython object
        GPEPython.__init__(self)

        self.name = name

        # simulation details
        self.evolutionType = None
        self.totalSteps = 0
        self.time = 0.0

        # data storage
        self.simulationFolder = None
        self.dataFiles = {}
        self.quantities = {}

        # set default parameters
        self.set("timeStep", 0.001)

        # evolution types
        self.setBool("evolutionKinetic",         True)
        self.setBool("evolutionPotential",       True)
        self.setBool("evolutionContact",         True)
        self.setBool("evolutionDipolar",         True)
        self.setBool("evolutionThreeBodyLosses", False)

    def __del__(self):
        """Close all data files"""

        for file in self.dataFiles.values():
            file.close()

    def getSimulationFolder(self):
        """returns the path of the current simulation folder (and creates it, if not existent)"""

        if not self.simulationFolder:
            self.simulationFolder = os.path.join(os.path.abspath(os.curdir), self.name)
            try:
                os.mkdir(self.simulationFolder)
            except:
                pass

        return self.simulationFolder

    def getDataFile(self, name):
        """returns the file handle of the [name].dat file within the current simulation folder"""

        if name not in self.dataFiles.keys():
            self.dataFiles[name] = open(os.path.join(self.getSimulationFolder(), name + ".dat"), "w")

        return self.dataFiles[name]

    def setGridSize(self, sizeX, sizeY=None, sizeZ=None):
        """initializes the grid"""

        if not sizeY:
            sizeY = sizeX
        if not sizeZ:
            sizeZ = sizeY

        self.set("sizeX", sizeX)
        self.set("sizeY", sizeY)
        self.set("sizeZ", sizeZ)

    def setMax(self, maxX, maxY=None, maxZ=None):
        """set maximum x, y and z position on the grid"""

        if not maxY:
            maxY = maxX
        if not maxZ:
            maxZ = maxY

        self.set("maxX", maxX)
        self.set("maxY", maxY)
        self.set("maxZ", maxZ)

    def evolutionCmd(self, steps):
        if self.evolutionType == "ite":
            GPEPython.ite(self, steps)
        else:
            GPEPython.rte(self, steps)

        self.totalSteps = self.totalSteps + steps
        self.time = self.time + self.get("timeStep") * steps

    def evolution(self, steps=1, monitorSteps=None, plotSteps=None):
        """run the simulation for a certain number of steps
        monitor the energy every 'monitorSteps' steps and
        write projections every 'plotSteps' steps"""

        if self.evolutionType not in ["ite", "rte"]:
            raise Exception("Wrong evolution type: " + type)

        if plotSteps and (plotSteps % monitorSteps != 0):
            raise Exception("plotSteps should be a multiple of monitorSteps")

        s = 0
        if monitorSteps is not None and monitorSteps > 0:
            self.monitor()
            if plotSteps is not None and plotSteps > 0:
                self.plot()

            while s < steps:
                self.evolutionCmd(monitorSteps)
                s = s + monitorSteps

                self.monitor()

                if plotSteps is not None and (s % plotSteps == 0):
                    self.plot()

        # simulate the remaining steps
        if s < steps:
            self.evolutionCmd(steps - s)

    def ite(self, steps=1, monitorSteps=None, plotSteps=None):
        """run imaginary time evolution"""

        self.evolutionType = "ite"
        self.evolution(steps, monitorSteps, plotSteps)

    def rte(self, steps=1, monitorSteps=None, plotSteps=None):
        """run real time evolution"""

        self.evolutionType = "rte"
        self.evolution(steps, monitorSteps, plotSteps)

    def energy(self):
        """returns the total energy"""
        Ekin = self.energyKinetic()
        Epot = self.energyPotential()
        Econ = self.energyContact()
        Edip = self.energyDipolar()

        return Ekin + Epot + Econ + Edip

    def energyInt(self):
        """returns the contact + dipolar interaction energy"""
        Econ = self.energyContact()
        Edip = self.energyDipolar()

        return abs(Econ) + abs(Edip)

    def virial(self):
        """returns viral"""
        Ekin = self.energyKinetic()
        Epot = self.energyPotential()
        Econ = self.energyContact()
        Edip = self.energyDipolar()

        return 2 * (Ekin - Epot) + 3 * (Econ + Edip)

    def energyString(self):
        Ekin = self.energyKinetic()
        Epot = self.energyPotential()
        Econ = self.energyContact()
        Edip = self.energyDipolar()
        Etot = Ekin + Epot + Econ + Edip
        Vir = 2 * (Ekin - Epot) + 3 * (Econ + Edip)

        energies = [Etot, Ekin, Epot, Econ, Edip, Vir]

        string = "{totalSteps:08}\t".format(totalSteps=self.totalSteps)
        string = string + "\t".join([format(E, " 10.9f") for E in energies]) + "\n"
        return string

    def monitor(self):
        """print energy values and evaluation quantities to monitoring files"""

        # monitor energy
        string = self.energyString()

        print(string.strip())

        file = self.getDataFile("energies_" + self.evolutionType)
        file.write(string)
        file.flush()

        # monitor quantities
        if len(self.quantities) > 0:
            string = self.evaluationString()
            file = self.getDataFile("quantities_" + self.evolutionType)
            file.write(string)
            file.flush()

    def printEnergy(self):
        """print a short summary of all energies"""

        print(self.energyString())

    def write(self, filename):
        """write 1D and 2D projections"""

        path = os.path.join(self.getSimulationFolder(), filename)

        GPEPython.write(self, path)

    def writePotential(self, filename):
        """write 1D and 2D projections of the external potential"""

        path = os.path.join(self.getSimulationFolder(), filename)

        GPEPython.writePotential(self, path)

    def plot(self):
        """write 1D and 2D projections with pre-defined name (step_XXXX_...)"""

        self.write("step_{totalSteps:08}_{type}".format(totalSteps=self.totalSteps, type=self.evolutionType))

    def evaluationString(self):
        """Generates a string with all evaluated quantities"""

        string = "{totalSteps:08}\t".format(totalSteps=self.totalSteps)
        values = []
        for quantity in self.quantities:
            if type(quantity) == str:
                value = self.get(quantity)
            elif type(quantity) == int:
                value = self.evaluate(quantity)
            elif hasattr(quantity, '__call__'):
                value = quantity(self)
            else:
                raise Exception("Wrong type of quantity")

            values.append(value)

        string = string + "\t".join([format(v, " 8.5f") for v in values]) + "\n"
        return string

    def setEvaluation(self, quantities):
        self.quantities = quantities


class HarmonicTrapGPE(GPE):
    """Defines a single harmonic trap simulation with oscillatory units defined by the frequency in z direction."""

    def __init__(self, name, omegaX=1, omegaY=None, omegaZ=None):
        GPE.__init__(self, name)

        # external potential
        self.setString("externalPotential", "harmonic")

        if not omegaY:
            omegaY = omegaX
        if not omegaZ:
            omegaZ = omegaX

        self.omegaX = omegaX
        self.omegaY = omegaY
        self.omegaZ = omegaZ

        length = math.sqrt(1.0 / omegaZ)
        self.length = length

        self.hoX = math.sqrt(1.0 / omegaX) / length
        self.hoY = math.sqrt(1.0 / omegaY) / length
        self.hoZ = 1.0

        self.setTrap()
        factor = 5
        self.setMax(factor * self.hoX,
                    factor * self.hoY,
                    factor * self.hoZ)

    def setTrap(self):
        """set all omega parameters and appropriate sigma's (harmonic oscillator)"""

        self.set("omegaX", self.omegaX)
        self.set("omegaY", self.omegaY)
        self.set("omegaZ", self.omegaZ)

        self.set("sigmaX", self.hoX)
        self.set("sigmaY", self.hoY)
        self.set("sigmaZ", self.hoZ)


class LatticeGPE(GPE):
    """Defines a lattice simulation"""

    def __init__(self, name, lSites, lDirection, lDepth, lWidth=0.5):
        """
            lSites:     number of lattice sites
            lDirection: 0, 1, 2 for x, y or z direction
            lDepth:     prefactor of the potential V0: - V0 * exp(...)
            lWidth:     either a single number or a vector of 3 different widths of the single well
        """

        GPE.__init__(self, name)

        # external potential
        self.setString("externalPotential", "lattice")

        self.set("latticeSites", lSites)
        self.set("latticeDirection", lDirection)
        self.set("latticeDepth", lDepth)

        # Geometry of a single well
        if not hasattr(lWidth, '__getitem__'):
            lWidth = 3 * [lWidth]

        self.set("latticeWidthX", lWidth[0])
        self.set("latticeWidthY", lWidth[1])
        self.set("latticeWidthZ", lWidth[2])

        # Initial state
        sigma = [0, 0, 0]
        for i in range(0, 3):
            sigma[i] = math.sqrt(lWidth[i]/(2.0 * math.sqrt(lDepth)))

        if lSites > 1:
            sigmaLat = 1.2 * (lSites - 1) / 2  # enlarge sigma in lattice direction (to load into the lattice)
            sigma[lDirection] = sigmaLat

        self.set("sigmaX", sigma[0])
        self.set("sigmaY", sigma[1])
        self.set("sigmaZ", sigma[2])

        # Grid
        maxP = [0, 0, 0]
        for i in range(0, 3):
            maxP[i] = 6 * lWidth[i]

        maxPLat = 2.0 + (float(lSites) - 1) / 2  # enlarge the grid in lattice direction
        maxP[lDirection] = maxPLat

        self.setMax(maxP[0], maxP[1], maxP[2])
