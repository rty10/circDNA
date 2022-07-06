from pymol import cmd

cmd.bg_color('white')
cmd.cartoon('rectangle')
cmd.set('cartoon_ring_mode', 1)
cmd.set('cartoon_ring_finder', 1)
cmd.set('cartoon_ring_color', 'gray50')
cmd.set('cartoon_ladder_mode', 1)
cmd.set('cartoon_ladder_color', 'black')
cmd.set('cartoon_nucleic_acid_mode', 1)
cmd.show('cartoon')

