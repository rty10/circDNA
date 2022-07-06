from pymol import cmd

cmd.bg_color('white')
cmd.cartoon('dumbbell')
cmd.set('cartoon_ring_mode', 0)
cmd.set('cartoon_ring_finder', 1)
cmd.set('cartoon_ring_color', 'gray50')
cmd.set('cartoon_ladder_mode', 1)
cmd.set('cartoon_ladder_color', 'black')
cmd.set('cartoon_nucleic_acid_mode', 2)
cmd.set('cartoon_dumbbell_radius', 0.8)

cmd.color('tv_red', 'resn dA in chain A')
cmd.color('tv_blue', 'resn dT in chain A')
cmd.color('tv_green', 'resn dG in chain A')
cmd.color('tv_orange', 'resn dC in chain A')
cmd.color('gray80', 'chain B')
cmd.show('cartoon')

