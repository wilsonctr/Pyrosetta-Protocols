import operator
from random import randint
import pyrosetta
import math
import os
from random import randint
from pyrosetta.toolbox import mutate_residue
import argparse

parser = argparse.ArgumentParser(description='Generic docking')
parser.add_argument('--input_pdb', help='Your input pdb.',required=True)
parser.add_argument('--constraints_file', help='User-defined, enzdes style constraints file.',required=True)
parser.add_argument('--output_folder_name', help='Output folder name for docked pdbs outputs.',required=True)
pyrosetta.init('-extra_res_fa X00.params Y00.params -preserve_header True -output_virtual T -enzdes_out T -score:weights ref2015_cst')

args = parser.parse_args()

class My_New_Mover(pyrosetta.rosetta.protocols.moves.Mover):
    def __init__(self,sfxn,pt):
        print( 'My_New_Mover.__init__...' )
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        self.sfxn = sfxn
        self.pt = pt
        
    def get_name(self): return 'My_New_Mover'
    
    def apply(self, p):
        #print( 'My_New_Mover.apply:', type(p) )
        #print( 'This My_New_Mover apply...' )
            
        cstopt = pyrosetta.rosetta.protocols.enzdes.EnzdesBaseProtocol()
        cstopt.set_scorefxn( self.sfxn )
        cstopt.set_minimize_options(True, False, True, True) # check this fn signature for details
            
        # this actuall runs the minimizer !
        cstopt.cst_minimize(p, self.pt, True)

pose = pyrosetta.pose_from_file(args.input_pdb)

sfxn = pyrosetta.get_fa_scorefxn()
pyrosetta.rosetta.protocols.enzdes.enzutil.enable_constraint_scoreterms(sfxn)

soft_rep = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("soft_rep")
pyrosetta.rosetta.basic.options.set_real_option('enzdes:cut1',6.0)
pyrosetta.rosetta.basic.options.set_real_option('enzdes:cut2',8.0)
pyrosetta.rosetta.basic.options.set_real_option('enzdes:cut3',10.0)
pyrosetta.rosetta.basic.options.set_real_option('enzdes:cut4',12.0)
pyrosetta.rosetta.basic.options.set_boolean_option('enzdes:detect_design_interface',True)
tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
dp = pyrosetta.rosetta.protocols.enzdes.DetectProteinLigandInterface()
dp.init_from_options()
dp.set_design(False)

extra_rotamers = pyrosetta.rosetta.core.pack.task.operation.ExtraRotamersGeneric()
extra_rotamers.ex1(True)
extra_rotamers.ex2(True)
extra_rotamers.extrachi_cutoff(1)
tf.push_back(dp)
tf.push_back(extra_rotamers)

predock = pyrosetta.rosetta.protocols.enzdes.PredesignPerturbMover()
predock.trans_magnitude(0.1)
predock.rot_magnitude(2)
## This set ligand must be set to the ligand, here it is hard coded for the example
predock.set_ligand(len(pose.sequence()))

enzdes_wbb = pyrosetta.rosetta.protocols.enzdes.EnzRepackMinimize()
enzdes_wbb.set_scorefxn_minimize(sfxn)
enzdes_wbb.set_min_lig(True)
enzdes_wbb.set_min_rb(True)
enzdes_wbb.set_min_sc(True)
enzdes_wbb.set_design(False)
enzdes_wbb.set_scorefxn_repack(soft_rep)
enzdes_wbb.set_min_bb(True)
enzdes_wbb.task_factory( tf )
        
addcsts = pyrosetta.rosetta.protocols.enzdes.AddOrRemoveMatchCsts()
addcsts.cstfile(args.constraints_file)
addcsts.set_cst_action( pyrosetta.rosetta.protocols.enzdes.CstAction.ADD_NEW )
addcsts.apply(pose)

# Create a packer task, specifically for the cstopt mover to work

pt = tf.create_task_and_apply_taskoperations(pose)
cstopt = My_New_Mover(sfxn,pt)
        
parsed = pyrosetta.rosetta.protocols.moves.SequenceMover()
parsed.add_mover(predock)
parsed.add_mover(enzdes_wbb)
parsed.add_mover(cstopt)

gm = pyrosetta.rosetta.protocols.qsar.scoring_grid.GridSet()
classic = pyrosetta.rosetta.protocols.qsar.scoring_grid.ClassicGrid()
classic.set_chain('X')
gm.width(30.0)
gm.add_grid('Classic', classic)
#(gm object above, chain: unicode, box_size: float, move_distance: float, angle: float, cycles: int, temp: float)
transform = pyrosetta.rosetta.protocols.ligand_docking.Transform(gm,'X',10.0,0.3,3.0,1000000,5.0)
transform.apply(pose)

mc = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover()
mc.set_drift(True) # This sets drift for the maxtrials (not technically mc anymore)
mc.set_maxtrials(20)  # CHANGE THIS TO 10 for real runs
mc.set_sampletype('low')
mc.set_temperature(0.6)
mc.set_mover(parsed)
mc.set_max_accepted_trials(20)
mc.set_scorefxn(sfxn)
mc.apply(pose)
sfxn(pose)

constraint_energy = pyrosetta.rosetta.protocols.enzdes.enzutil.sum_constraint_scoreterms(pose, -1)
random_number = randint(1, 100000)
end = len(pose.sequence())
os.system('mkdir %s/docking_results'%args.output_folder_name)
pose.dump_pdb('%s/docking_results/%s_%s'%(args.output_folder_name, randint(1, 100000), args.input_pdb.split('/')[1]))
