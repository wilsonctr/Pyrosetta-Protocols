def score_pdb(pdb_file):
    
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


def energies_and_seq_extract(scored_pdb):
    
    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}

    with open(scored_pdb, 'r') as f:
        
        sequence = []
    
        prev = '-1'
        for line in f:
            toks = line.split()
            if len(toks)<1: continue
            if toks[0] != 'ATOM': continue
            if toks[5] != prev:
                sequence.append('%c' % letters[toks[3]])
            prev = toks[5]

        protein_seq = ''.join(sequence)
    
    with open(scored_pdb, 'r') as f:
        
        data_dict = {}

        for i in f.readlines():
            if 'pose' in i:
                scores = i.split()
            if 'label' in i:
                labels = i.split()

        for j,k in zip(labels, scores):
            data_dict[j] = k
            
        data_dict['Name'] = scored_pdb
        data_dict['Sequence'] = protein_seq
        
    return data_dict
