#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BRANCHING AND ANNIHILATING RANDOM WALK TEST SIMULATION

Uses the modules defined in the scripts "branching_rules" and "branching_simulation"

@author: mucar
"""
import numpy as np
import branching_simulation as sim
import matplotlib.pyplot as plt

""" 
Decide on branching probability, external guidance and self-interaction strengths:
"""
prob = 0.05  # branching probability
fc = 0.1    # external guidance strength
fs = -0.1   # self-avoidance strength
tmax = 200  # maximal simulation time

# Run simulation:
test_run = sim.simulation_loop(prob,fc,fs,tmax)

# Assign coordinates, active tip lengths over time, and angle values:
coordinates, evolve, angles = test_run['coordinates'], test_run['evolve'], test_run['angles']

# Save data to analyze:
np.save(f'data/test_simulation_pb_{prob}_fc_{fc}_fs_{fs}_coordinates.npy',coordinates)
np.save(f'data/test_simulation_pb_{prob}_fc_{fc}_fs_{fs}_angles.npy',angles)
np.save(f'data/test_simulation_pb_{prob}_fc_{fc}_fs_{fs}_evolve.npy',evolve)

# Plot the entire network:
plt.plot(coordinates[:,0],coordinates[:,1],'o',ms=2)
