from pymol import *

cmd.load("hueco.pdb")
cmd.load("1mtn.pdb")
cmd.load("bbox.pdb")
cmd.load("in_bbox.pdb")

cmd.color("salmon", "in_bbox")
cmd.show("spheres", "in_bbox")
