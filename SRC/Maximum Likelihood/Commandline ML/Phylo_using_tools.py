from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from io import StringIO
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline

P10 = open("P10_strain.fasta", "r")
R50 = open("R50_strain.fasta", "r")
AR465 = open("AR465_strain.fasta", "r")
M48 = open("M48_strain.fasta", "r")
NRS1 = open("NRS1_strain.fasta", "r")
NCTC_8325 = open("NCTC_8325_strain.fasta", "r")

multifile = open("Multi-fasta.fasta", "w+")


def merging_sequences(file):
    """This function creates a Multi-Fasta file
        in case user enters multiple files of genes of interest"""
    for line in file:
        line1 = line.rstrip()
        if line1.startswith(">") or line1.endswith("\t"):
            var = multifile.write(str(line1)+"\n")
        elif line1[1::].startswith("A") or line1[1::].startswith("T") or \
                line1[1::].startswith("C") or line1[1::].startswith("G"):
            sequence = line[1:71]
            var2 = multifile.write(str(sequence))

    return multifile


merging_sequences(P10)
merging_sequences(R50)
merging_sequences(AR465)
merging_sequences(NCTC_8325)
merging_sequences(NRS1)
merging_sequences(M48)

multifile.close()

in_file = "Multi-fasta.fasta"   # test files!
out_file = "alignment.phylip"

# specify the location of your muscle exe file
muscle_exe = r"F:\\New folder (2)\\lab 7\\muscle3.8.425_win32.exe"


def align_v1(Fasta):
    muscle_cline = MuscleCommandline(muscle_exe, input=Fasta, out=out_file)
    stdout = muscle_cline(stdin=in_file)
    stderr = muscle_cline(stdin=in_file)
    MultipleSeqAlignment = AlignIO.read(out_file, "fasta")
    with AlignIO.as_handle("MultipleSeqAlignment.phylip","w") as msa:
        AlignIO.write(MultipleSeqAlignment, msa, "phylip")
    return msa


align_v1(in_file)


infile = "MultipleSeqAlignment.phylip"

# specify the location of your executable.
cmd_exe = r"F:\PhyML-3.1\PhyML-3.1\PhyML-3.1_win32.exe"


def phyml(PHYLIP):
    phyml_cline = PhymlCommandline(cmd_exe, input=PHYLIP)
    stdout = phyml_cline(stdin=infile)
    stderr = phyml_cline(stdin=infile)
    return stdout

phyml(infile)

tree = Phylo.read("MultipleSeqAlignment.phylip_phyml_tree.txt", "newick")
tree.rooted = false
Phylo.draw_ascii(tree)
