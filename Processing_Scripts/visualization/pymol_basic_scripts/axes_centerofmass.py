# axes.py
from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain

# for pymol.cgo reference, refer to http://pymol.sourceforge.net/newman/user/S0500cgo.html
# create the axes object, draw axes with cylinders coloured red, green, blue for X, Y and Z

cx, cy, cz = cmd.centerofmass()

# Note for load_cgo:
# cylSpec, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2
radius = 0.2


# let x be red, y be green, z be blue
obj = [
   CYLINDER, cx, cy, cz, 10.0+cx, cy, cz, radius, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
   CYLINDER, cx, cy, cz, cx, 10.0+cy, cz, radius, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
   CYLINDER, cx, cy, cz, cx, cy, 10.0+cz, radius, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
   ]

cmd.load_cgo(obj,'axes')
