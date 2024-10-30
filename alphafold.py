from pymol import cmd
import sys

# get the protein name in the argument

prot = sys.argv[1] 

file_crs = "data/proteins2fold/" + prot.upper() + "/" + "fold_" + prot +"_crs_model_1.cif"
file_deni = "data/proteins2fold/" + prot.upper() + "/" + "fold_" + prot +"_deni_model_1.cif"
f = open("data/proteins2fold/" + prot.upper() + "/" + prot.upper() + "_mut_resi.txt", "r")
mutations = f.readlines()
resi=mutations[0].replace(" ", "+")



# load the files
cmd.load(file_crs, "crs_"+prot)
cmd.load(file_deni, "deni_"+prot)

cmd.align("deni_"+prot, "crs_"+prot)
cmd.select("mutations", "resi " + resi + " and deni_" + prot)
cmd.select("wt", "resi " + resi + " and crs_" + prot)
cmd.show("licorice", "wt")
cmd.show("licorice", "mutations")
cmd.color("red", "mutations")
cmd.color("yellow", "wt")



