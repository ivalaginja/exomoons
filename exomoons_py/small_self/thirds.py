"""
This program calculated the ring radii in a system in which the gaps get filled with a lane of dust in their middle third. Then transferred by hand.
"""

import numpy as np

inradii = [23.80214119,24.65023232,37.55720901,39.52087784,60.06858826,62.88878632,78.04033661]

gap1 = inradii[1]-inradii[0]
gap2 = inradii[3]-inradii[2]
gap3 = inradii[5]-inradii[4]

third1 = gap1/3.
third2 = gap2/3.
third3 = gap3/3.

rad1 = inradii[0]
rad2 = inradii[0] + third1
rad3 = inradii[0] + third1*2
rad4 = inradii[1]
rad5 = inradii[2]
rad6 = inradii[2] + third2
rad7 = inradii[2] + third2*2
rad8 = inradii[3]
rad9 = inradii[4]
rad10 = inradii[4] + third3
rad11 = inradii[4] + third3*2
rad12 = inradii[5]
rad13 = inradii[6]

radii = [rad1,rad2,rad3,rad4,rad5,rad6,rad7,rad8,rad9,rad10,rad11,rad12,rad13]

print radii
