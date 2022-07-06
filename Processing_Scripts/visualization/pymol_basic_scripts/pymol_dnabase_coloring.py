from pymol import cmd
# pymol script: dna base pair coloring

# -1- set background color:
cmd.bg_color("white")

# -2- set dna bases to desired color:
cmd.color("tv_red", "resn dA in chain A")
cmd.color("tv_red", "resn dT in chain B")

cmd.color("tv_blue", "resn dT in chain A")
cmd.color("tv_blue", "resn dA in chain B")

cmd.color("tv_green", "resn dG in chain A")
cmd.color("tv_green", "resn dC in chain B")

cmd.color("tv_orange", "resn dC in chain A")
cmd.color("tv_orange", "resn dG in chain B")

