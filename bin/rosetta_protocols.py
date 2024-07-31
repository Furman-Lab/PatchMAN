

from pyrosetta import *
from pyrosetta.rosetta import *

#Core
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation

#Protocols
from pyrosetta.rosetta.protocols import minimization_packing as pack_min

ONE_LETTER_AA = ('G', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'R',
                 'C', 'N', 'Q', 'T', 'S', 'P', 'H', 'K', 'D', 'E')

ONE_TO_THREE_AA = {'G': 'GLY',
                   'A': 'ALA',
                   'V': 'VAL',
                   'L': 'LEU',
                   'I': 'ILE',
                   'P': 'PRO',
                   'C': 'CYS',
                   'M': 'MET',
                   'H': 'HIS',
                   'F': 'PHE',
                   'Y': 'TYR',
                   'W': 'TRP',
                   'N': 'ASN',
                   'Q': 'GLN',
                   'S': 'SER',
                   'T': 'THR',
                   'K': 'LYS',
                   'R': 'ARG',
                   'D': 'ASP',
                   'E': 'GLU'}

THREE_TO_ONE_AA = {'GLY': 'G',
                   'ALA': 'A',
                   'VAL': 'V',
                   'LEU': 'L',
                   'ILE': 'I',
                   'PRO': 'P',
                   'CYS': 'C',
                   'MET': 'M',
                   'HIS': 'H',
                   'PHE': 'F',
                   'TYR': 'Y',
                   'TRP': 'W',
                   'ASN': 'N',
                   'GLN': 'Q',
                   'SER': 'S',
                   'THR': 'T',
                   'LYS': 'K',
                   'ARG': 'R',
                   'ASP': 'D',
                   'GLU': 'E'}


def fixbb_design(lig_selection, filename, pepseq, scrfxn):
    """Receives indices of the residues to mutate, filename, peptide sequence and score function,
    run fixbb and dumps designed pdb if successful"""
    complex_pose = pose_from_pdb(filename)

    task_factory = core.pack.task.TaskFactory()
    packer_task = task_factory.create_packer_task(complex_pose)
    packer_task.restrict_to_residues(lig_selection)

    lig_nums = list(map(str, core.select.residue_selector.selection_positions(lig_selection))) # list of str

    with open(filename+'_resfile', 'w') as resfile:
        resfile.write('NATRO\nSTART\n')
        for i, res in enumerate(pepseq):
            if complex_pose.chain_sequence(2)[i] != res:  # to keep the original rotamers
                line = '{orig_aa} {chain} PIKAA {new_aa} EX 1 EX 2\n'.format(orig_aa=lig_nums[i],
                                                                             chain=complex_pose.pdb_info().chain(int(lig_nums[i])),
                                                                             new_aa=res)
            else:
                line = '{orig_aa} {chain} NATRO \n'.format(orig_aa=lig_nums[i],
                                                           chain=complex_pose.pdb_info().chain(int(lig_nums[i])))
            resfile.write(line)

    read_resfile = core.pack.task.operation.ReadResfile()
    read_resfile.filename(filename+'_resfile')
    read_resfile.apply(complex_pose, packer_task)

    fixbb = protocols.minimization_packing.PackRotamersMover(scrfxn, packer_task)

    try:
        fixbb.apply(complex_pose)  # If failing with unrecognized residue
    except RuntimeError:
        return 1
    complex_pose.dump_pdb(os.path.splitext(filename)[0] + '_0001.pdb')

    os.remove(filename + '_resfile')


def prepack(pose):

    new_pose = pose.clone()

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())
    tf.push_back(operation.IncludeCurrent())
    tf.push_back(operation.NoRepackDisulfides())

    packer = pack_min.PackRotamersMover()
    packer.task_factory(tf)

    packer.apply(new_pose)

    return new_pose

def prepack_min(pose):
    new_pose = pose.clone

    prepack(new_pose)

    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(True)

    minmover = pack_min.MinMover()
    minmover.movemap(mm)

    minmover.apply(new_pose)

    return new_pose
