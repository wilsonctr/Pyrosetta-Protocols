import pyrosetta
import math
pyrosetta.init('')

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


def CA_distance(res_num_1, res_num_2):

    res_1 = {}
    res_2 = {}
    for lig_atom in str(pose.residue(res_num_1)).split('Atom Coordinates:')[1].replace(' ','').split('\n')[1:-2]:
        res_1[lig_atom.split(':')[0]] = lig_atom.split(':')[1]

    for lig_atom in str(pose.residue(res_num_2)).split('Atom Coordinates:')[1].replace(' ','').split('\n')[1:-2]:
        res_2[lig_atom.split(':')[0]] = lig_atom.split(':')[1]

    def distance_between_two_xyz_coords(a,b):
        distance = math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
        return distance

    CA_coords = [float(i) for i in res_1['CA'].split(',')]

    distances = []
    for phosphate_oxygen in phosphate_atoms:
        PO_coords = [float(i) for i in res_2[phosphate_oxygen].split(',')]
        distances.append(distance_between_two_xyz_coords(CA_coords, PO_coords))

    return distances
