import pyrosetta
pyrosetta.init()

mm = pyrosetta.MoveMap()
for i in resi_within(pose, 6, 262):
    mm.set_bb(i, True)
    mm.set_chi(i,True)
backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
backrub.set_movemap(mm)

for i in range(1,20):
    pose = pyrosetta.pose_from_file('<PDB>)
    for j in range(1,20):
        backrub.apply(pose)
