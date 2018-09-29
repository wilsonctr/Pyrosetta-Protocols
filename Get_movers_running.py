#Loop modeling

loop_modeling = pyrosetta.rosetta.protocols.loop_modeling.LoopModeler()
loops = pyrosetta.rosetta.protocols.loops.Loop (38, 42, 40, 0, False)
loop_modeling.set_loop(loops)
loop_modeling.apply(pose)

#Docking using transform mover

gm = pyrosetta.rosetta.protocols.qsar.scoring_grid.GridSet()
classic = pyrosetta.rosetta.protocols.qsar.scoring_grid.ClassicGrid()
classic.set_chain('X')
gm.width(40.0)
gm.add_grid('Classic', classic)
#(gm object above, chain: unicode, box_size: float, move_distance: float, angle: float, cycles: int, temp: float)                                   
transform = pyrosetta.rosetta.protocols.ligand_docking.Transform(gm,'X',20.0,0.1,1.0,1000000,5.0)
transform.apply(pose)

#H-bond count filter between residues

hbond_filter = pyrosetta.rosetta.protocols.protein_interface_design.filters.HbondsToResidueFilter()
hbond_filter.set_scorefxn(sfxn)
hbond_filter.set_partners(1)
hbond_filter.set_resnum(resi)
hbond_filter.set_energy_cutoff(-1.0)
hbond_filter.set_selector(resi_selector)

# Rotamer explosion

extra_rotamers = pyrosetta.rosetta.core.pack.task.operation.ExtraRotamersGeneric()
extra_rotamers.ex1(True)
extra_rotamers.ex2(True)
extra_rotamers.ex3(True)
extra_rotamers.ex4(True)
extra_rotamers.ex1_sample_level(pyrosetta.rosetta.core.pack.task.ExtraRotSample(4))                                                                                    
extra_rotamers.ex2_sample_level(pyrosetta.rosetta.core.pack.task.ExtraRotSample(4))                                                                                    
extra_rotamers.ex3_sample_level(pyrosetta.rosetta.core.pack.task.ExtraRotSample(4))                                                                                    
extra_rotamers.ex4_sample_level(pyrosetta.rosetta.core.pack.task.ExtraRotSample(4))                                                                                    
extra_rotamers.extrachi_cutoff(1)
tf.push_back(extra_rotamers)
