#!/usr/bin/python

from gpe import LatticeGPE, GPE

import os

REPULSIVE = 0
ATTRACTIVE = 1
text = ["repulsive", "attractive"]

# Settings
geometry = REPULSIVE 
sites = 3

ite_steps_init = 2000
energy_steps = 200
ite_steps = 1000

depth = 80
d_add = 0.1
d_a = d_add

start_a = 1.0
start_add = 0.0
stop_add = 5.0

# End Settings

def getAtomnumbers(filen):
    f = file(filen)
    lines = f.readlines()
    n = [0.0, 0.0, 0.0]
    lastx = 0
    dx = 0
    for l in lines:
        [x, dens] = l.split(" ")
        dens = float(dens.strip())
        x = float(x)
        
        if x < -0.5:
            n[0] = n[0] + dens
        else:
            if x < 0.5:
                n[1] = n[1] + dens
            else: 
                n[2] = n[2] + dens
        
        if dx == 0:
            if lastx == 0:
                lastx = x
            else:
                dx = x - lastx

    for i in range(0, len(n)):
        n[i] = n[i] * dx

    f.close()
    return n

if geometry == REPULSIVE:
    direction = 0
    width = [0.5, 4, 0.5]
else: # ATTRACTIVE
    direction = 2
    width = [0.5, 4, 0.5]

text_dir = ["x", "y", "z"]

add = start_add
while add <= stop_add:
    title = "{geometry}_{add}".format(add=add, geometry = text[geometry])
    print title
    sim = LatticeGPE(title, sites, direction, depth, width)

    # Grid parameters
    sizeLat = 128
    sizePerp = 32
    if geometry == REPULSIVE:
        sim.setGridSize(sizeLat, sizePerp, sizePerp)
        sim.set("maxX", 2.5)
        sim.set("maxY", 2.5)
        sim.set("maxZ", 1.5)
    else:
        sim.setGridSize(sizePerp, sizePerp, sizeLat)
        sim.set("maxX", 1.5)
        sim.set("maxY", 1.5)
        sim.set("maxZ", 2.5)

    # Contact and dipolar interaction
    sim.set("atomNumber",       2)
    sim.set("scatteringLength", start_a)
    sim.set("dipolarLength",    add)

    sim.setBool("evolutionContact", False)
    sim.setBool("evolutionDipolar", False)

    sim.initialize()
    sim.write("initial")

    sim.ite(ite_steps_init, energy_steps)
    sim.write("non_interacting")

    sim.setBool("evolutionContact", True)
    sim.ite(ite_steps_init, energy_steps)
    sim.write("contact_interacting")

    if add > 0:
        sim.setBool("evolutionDipolar", True)

    file_ratio = file(title + "/ratio", "a")
    collapsed = False
    a = start_a
    while not collapsed:
        print "a = {a}\tadd = {add}\n".format(a = a, add = add)
        filen = "state_{a}".format(a = a)
        filep = title + "/" + filen + "_project_" + text_dir[direction] + ".data"
        
        if os.path.exists(filep):
            print "State exists already... skipping"
            break

        if a == start_a:
            sim.ite(ite_steps_init, energy_steps)
        else:
            sim.ite(ite_steps, energy_steps)

        sim.write(filen)
        
        [n1, n2, n3] = getAtomnumbers(filep)
        print "n1 = {n1}, n2 = {n2}, n3 = {n3}".format(n1 = n1, n2 = n2, n3 = n3)
        ratio = (n1 + n3)/(2 * n2)
        ratio_n1n3 = n1/n3

        if ratio_n1n3 < 0.9 or ratio_n1n3 > 1.1 or ratio < 0.001:
            # Collapse
            ratio = 0.0 # per definition
            collapsed = True
        
        file_ratio.write("{a} {add} {ratio} {n1} {n2} {n3}\n".format(a = a, add = add, ratio = ratio, n1 = n1, n2 = n2, n3 = n3))
        
        a = a - d_a
        if abs(a) < 1e-10:
            a = 0.0
        sim.set("scatteringLength", a)

    del sim
    add = add + d_add

    file_ratio.close()
