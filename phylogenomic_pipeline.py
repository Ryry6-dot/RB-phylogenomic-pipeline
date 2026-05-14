# Import needed modules
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo

####################################################
#Set our file names
indir = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
outdir = "/scratch/boldajir/BB485/06/outdir/"

# Start a blank list to store our topology results
topology_results = []

###Start loop here####
file_list = glob.glob(indir + "*")

for file_path in file_list:   
    #print(file_path)
    ########################################################
    # Create a new file path pointing to the output directory (this is how we tell mafft what to name the output)
    new_file_path = file_path.replace(indir, outdir)

    print(new_file_path)

    # Create a command string (this is what get called using the 'system call'.
    aln_cmd = 'mafft --auto --quiet '+file_path+' > '+new_file_path

    #Check the command
    print(aln_cmd)

    #Run the command
    os.system(aln_cmd) #Uncomment this once you've double-checked that it's looking good.

    ############################################################################################

    #Create the command. -nt 2 means two threads. If running this from within a job submission, you could use more threads to make it go faster.
    tree_command = f"iqtree -s {new_file_path} -m TEST -nt 24"

    #Check the command 
    print(tree_command)

    #Run the command using a 'system call'
    os.system(tree_command) #uncomment once you've check the command

    ##################################################################################################
    tree = new_file_path + ".treefile"


    #Read in the tree and store as phylo object
    temp_tree = Phylo.read(tree, "newick")

    #Loop through the tips in the tree to find which one contains Es (the outgroup)
    for tip in temp_tree.get_terminals():
        if "Es_" in tip.name:
            es_tip = tip
            #Stope the loop once we found the correct tip
            break
        
    #Root the tree by the outgroup taxon
    temp_tree.root_with_outgroup(es_tip)
        
    #Get a list of all terminal (aka tips) branches
    all_terminal_branches = temp_tree.get_terminals()
        
    #Loop through the branches and store the names of the tips of each
    for t in all_terminal_branches:
        if "Bs_" in t.name:
            Bs_temp = t 
        elif "Cr_" in t.name:
            Cr_temp = t
        elif "At_" in t.name:
            At_temp = t
        else:
            out_temp = t
            
    #Make lists of pairs of branches, so that we can ask which is monophyletic
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]
        

    #Use series of if/else statements to ask which pair in monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str = "Unknown"

    #append topo string to a list
    topology_results.append(topo_str)
    # print(topo_str)

###end loop here

from collections import Counter

counts = Counter(topology_results)
print(counts)
######################################################################################