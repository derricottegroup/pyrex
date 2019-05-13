import os
import sys
import itertools
import subprocess
import json
import numpy as np

# What type of Linking Do you Want? #
link_type = "By Charge"
#link_type = "50-50"


# Get Fragments #
json_input = sys.argv[1]
json_data=open(json_input).read()
input_params = json.loads(json_data)
monomer_A_labels = input_params['fsapt']['monomer_A_labels']
monomer_B_labels = input_params['fsapt']['monomer_B_labels']
monomer_A_labels.append("All")
monomer_B_labels.append("All")
print(monomer_A_labels)
print(monomer_B_labels)
# Grab coordinates from the irc file
irc_filename = input_params['pyrex']['irc_filename']
step_size = input_params['pyrex']['irc_stepsize']
force_min = input_params['pyrex']['force_min']
force_max = input_params['pyrex']['force_max']
irc_file = open(irc_filename, "r")
coordinates = []
for line in irc_file:
    if("Full IRC Point" in line):
        coord_line = line.split()
        coord_num = int(coord_line[3])
        coordinates.append(float(coord_num*step_size))
unique_pairs_list = list(itertools.product(monomer_A_labels,monomer_B_labels))
unique_pairs = []
for i in range(len(unique_pairs_list)):
    unique_pairs.append("%s-%s" %(unique_pairs_list[i][0],unique_pairs_list[i][1]))
print(unique_pairs)

# Initialize empty dictionaries to hold pair energy, force(f), and work(w) data
pair_elst, pair_exch, pair_indab, pair_indba, pair_disp, pair_total = ({} for i in range(6))
f_elst, f_exch, f_indab, f_indba, f_disp, f_total = ({} for i in range(6))
w_elst_1, w_exch_1, w_indab_1, w_indba_1, w_disp_1, w_total_1 = ({} for i in range(6))
w_elst_2, w_exch_2, w_indab_2, w_indba_2, w_disp_2, w_total_2 = ({} for i in range(6))

# Initialize empty array for every pair property
for pair in unique_pairs:
    pair_elst[pair], pair_exch[pair], pair_indab[pair], pair_indba[pair], pair_disp[pair], pair_total[pair] = ([] for i in range(6))
    #f_elst[pair], f_exch[pair], f_indab[pair], f_indba[pair], f_disp[pair], f_total[pair] = ([] for i in range(6))

#print(pair_elst)

os.chdir("fsapt_output")
num_geoms = len(os.listdir())
num_ts = 135 #Hard coded number for transition state of this system, delete line to make general
num_geoms = 136 # NOTE: Use a few geometries past the transition state to get better agreement with pyrex data
for i in range(num_geoms):
    fsapt_dir = "fsapt%d" %i
    os.chdir(fsapt_dir)
    print("running in %s" %fsapt_dir)
    subprocess.call(["fsapt.py"])
    if(link_type=="By Charge"):
        do_store = True
    else:
        do_store = False
    fsapt_file = open('fsapt.dat', 'r')
    for line in fsapt_file:
        if("Reduced Analysis" in line):
            if(do_store):
                for j in range(3):
                    line = next(fsapt_file)
                for j in range(len(unique_pairs)):
                    pair_data = line.split()
                    pair = "%s-%s" %(pair_data[0],pair_data[1])
                    elst = pair_data[2]
                    exch = pair_data[3]
                    indab = pair_data[4]
                    indba = pair_data[5]
                    disp = pair_data[6]
                    total = pair_data[7]
                    pair_elst[pair].append(float(elst))
                    pair_exch[pair].append(float(exch))
                    pair_indab[pair].append(float(indab))
                    pair_indba[pair].append(float(indba))
                    pair_disp[pair].append(float(disp))
                    pair_total[pair].append(float(total))
                    #print(pair_elst)
                    line = next(fsapt_file) 
            if(link_type=="By Charge"):
                do_store = False
            else:
                do_store = True
    os.chdir("..")
pair_dicts = [pair_elst, pair_exch, pair_indab, pair_indba, pair_disp, pair_total]
f_dicts = [f_elst, f_exch, f_indab, f_indba, f_disp, f_total]
w_dicts_1 = [w_elst_1, w_exch_1, w_indab_1, w_indba_1, w_disp_1, w_total_1]
w_dicts_2 = [w_elst_2, w_exch_2, w_indab_2, w_indba_2, w_disp_2, w_total_2]
for dict_ in pair_dicts:
    index = pair_dicts.index(dict_)
    print(index)
    for pair in unique_pairs:
          f_dicts[index][pair] = -1.0*np.gradient(dict_[pair],step_size)


# Store energy data in .csv files
elst_csv = open("fsapt_elst.csv", "w+")
exch_csv = open("fsapt_exch.csv", "w+")
indab_csv = open("fsapt_indab.csv", "w+")
indba_csv = open("fsapt_indba.csv", "w+")
disp_csv = open("fsapt_disp.csv", "w+")
total_csv = open("fsapt_total.csv", "w+")

csv_files = [elst_csv,exch_csv,indab_csv,indba_csv,disp_csv,total_csv]

for csv_f in csv_files:
    index = csv_files.index(csv_f)
    prop_dict = pair_dicts[index]
    csv_f.write("Coordinate,")
    for pair in unique_pairs:
        if(pair==unique_pairs[-1]):
            csv_f.write("%s,\n" %pair)
        else:
            csv_f.write("%s," %pair)
    for i in range(num_geoms):
        csv_f.write("%f," %coordinates[i])
        for pair in unique_pairs:
            if(pair==unique_pairs[-1]):
                csv_f.write("%.3f,\n" %prop_dict[pair][i])
            else:
                csv_f.write("%.3f," %prop_dict[pair][i])
elst_csv.close(), exch_csv.close(), indab_csv.close(), indba_csv.close(), disp_csv.close(), total_csv.close()

# Store force data in .csv files
f_elst_csv = open("force_elst.csv", "w+")
f_exch_csv = open("force_exch.csv", "w+")
f_indab_csv = open("force_indab.csv", "w+")
f_indba_csv = open("force_indba.csv", "w+")
f_disp_csv = open("force_disp.csv", "w+")
f_total_csv = open("force_total.csv", "w+") 

f_csv_files = [f_elst_csv,f_exch_csv,f_indab_csv,f_indba_csv,f_disp_csv,f_total_csv]

for csv_f in f_csv_files:
    index = f_csv_files.index(csv_f)
    prop_dict = f_dicts[index]
    csv_f.write("Coordinate,")
    for pair in unique_pairs:
        if(pair==unique_pairs[-1]):
            csv_f.write("%s,\n" %pair)
        else:
            csv_f.write("%s," %pair)
    for i in range(num_geoms):
        csv_f.write("%f," %coordinates[i])
        for pair in unique_pairs:
            if(pair==unique_pairs[-1]):
                csv_f.write("%.3f,\n" %prop_dict[pair][i])
            else:
                csv_f.write("%.3f," %prop_dict[pair][i])
f_elst_csv.close(), f_exch_csv.close(), f_indab_csv.close(), f_indba_csv.close(), f_disp_csv.close(), f_total_csv.close()

# Calculate Work Contributions for each Pair
index_min = coordinates.index(force_min) 
f_dict_count = 0
for dict_ in f_dicts:
    for pair in unique_pairs:
          w_dicts_1[f_dict_count][pair] = -1.0*np.trapz(dict_[pair][:index_min],dx=step_size)
    f_dict_count = f_dict_count + 1
f_dict_count = 0
for dict_ in f_dicts:
    for pair in unique_pairs:
          w_dicts_2[f_dict_count][pair] = -1.0*np.trapz(dict_[pair][index_min-1:num_ts],dx=step_size)
    f_dict_count = f_dict_count + 1

work_values = open("work_values.dat", "w+")
work_values.write('\n\n--Reaction Work Decomposition (Region 1)--\n')
work_values.write('\n-----------------------------------------------------------------------------------------------------------------')
work_values.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('Pair(A-B)','W_elst','W_exch', 'W_indAB', 'W_indBA', 'W_disp', 'W_total'))
work_values.write('-----------------------------------------------------------------------------------------------------------------\n')
#print(w_dicts[0])
for pair in unique_pairs:
    work_values.write('\n{:>15s} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(pair,w_dicts_1[0][pair],w_dicts_1[1][pair], w_dicts_1[2][pair], w_dicts_1[3][pair], w_dicts_1[4][pair], w_dicts_1[5][pair]))
work_values.write('------------------------------------------------------------------------------------------------------------------\n')
work_values.write('\n\n--Reaction Work Decomposition (Region 2)--\n')
work_values.write('\n-----------------------------------------------------------------------------------------------------------------')
work_values.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('Pair(A-B)','W_elst','W_exch', 'W_indAB', 'W_indBA', 'W_disp', 'W_total'))
work_values.write('-----------------------------------------------------------------------------------------------------------------\n')
#print(w_dicts[0])
for pair in unique_pairs:
    work_values.write('\n{:>15s} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(pair,w_dicts_2[0][pair],w_dicts_2[1][pair], w_dicts_2[2][pair], w_dicts_2[3][pair], w_dicts_2[4][pair], w_dicts_2[5][pair]))
work_values.write('------------------------------------------------------------------------------------------------------------------\n')
work_values.close()
