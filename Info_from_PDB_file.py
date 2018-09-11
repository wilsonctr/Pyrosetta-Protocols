def data_from_pdb(pdb_file):
    
    sfxn = pyrosetta.get_fa_scorefxn()
    emeth_opt = sfxn.energy_method_options().clone()
    emeth_opt.hbond_options().decompose_bb_hb_into_pair_energies( True )
    sfxn.set_energy_method_options( emeth_opt )
    pose = pyrosetta.pose_from_file(pdb_file)
    seq = pose.sequence()
    total_energy = sfxn(pose)
    name = str(pdb_file)
    
    info_dict = {}
    for i in range(1,300):
        info_dict[str(pyrosetta.rosetta.core.scoring.ScoreType(i)).split('.')[1]] = pyrosetta.rosetta.core.scoring.ScoreFunction.score_by_scoretype(sfxn, pose, pyrosetta.rosetta.core.scoring.ScoreType(i))
    
    info_dict['Sequence'] = str(seq)
    info_dict['total_energy'] = sfxn(pose)
    info_dict['File_name'] = name
    
    return info_dict
