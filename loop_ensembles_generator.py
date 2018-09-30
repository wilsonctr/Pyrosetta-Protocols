import rosetta
import pyrosetta
import math
from pyrosetta.toolbox import mutate_residue
from itertools import groupby
from operator import itemgetter
pyrosetta.init()

def resi_within(pose, angstrom, the_sun):
    
    def resi_atom_coords(res_in_pose): # input = the pose.residue(number) object
        xyz_coords_list = [float(i) for i in str(res_in_pose.atoms())[31:-1].split(',')] # All atom coords in a residue
        atoms_coords_list = [xyz_coords_list[j:j+3] for j in range(0, len(xyz_coords_list), 3)] # List of lists. xyz of individual residue
        return atoms_coords_list
    
    def distance_between_two_xyz_coords(a,b):
        distance = math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
        return distance
    
    resi_within = []
    center = resi_atom_coords(pose.residue(the_sun))
    for resi in range(1, pose.total_residue()+1): # loop over protein
        planets = resi_atom_coords(pose.residue(resi))
        for coords in center:
            for other_coords in planets:
                distance = distance_between_two_xyz_coords(coords,other_coords)
                if distance <= angstrom:
                    resi_within.append(resi)
                
    return sorted(list(set(resi_within)))[0:-1]

pose = pyrosetta.pose_from_file(<PDB>)

sfxn = pyrosetta.get_fa_scorefxn()

dssp = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)
ss = dssp.get_dssp_secstruct()

#Active site loops

loops_at_interface = []
all_loops = []
for j,k in zip(ss, range(1,len(pose.sequence()))):
    if j == 'L':
        all_loops.append(k)

active_site_resi = resi_within(pose, 10, len(pose.sequence()))

for a,b in groupby(enumerate(all_loops), lambda (i, x): i-x):
    data = map(itemgetter(1), b)
    if len(data) > 3 and len(data) < 8:
        for resi in data:
            if resi in active_site_resi:
                loops_at_interface.append(data)

kic = pyrosetta.rosetta.protocols.kinematic_closure.KicMover()
loops = pyrosetta.rosetta.protocols.loops.Loops()
for loop in loops_at_interface:
    loops.add_loop(loop[0], loop[int(len(loop)/2)], loop[1], 0, False)
kic.set_loops(loops)

mc = rosetta.protocols.simple_moves.GenericMonteCarloMover()
mc.set_drift(True) # this sets drift for the maxtrials (not technically mc anymore)
mc.set_maxtrials(20)  # CHANGE THIS TO 10 for real runs
mc.set_sampletype('low')
mc.set_temperature(0.6)
mc.set_mover(kic)
mc.set_scorefxn(sfxn)

mc.apply(pose)
