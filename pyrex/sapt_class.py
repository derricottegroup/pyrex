import os
import psi4
import itertools
import subprocess
import numpy as np


class sapt(object):

    def __init__(self, geometries, data, method, basis, outfile):
        self.geometries = geometries
        self.method = method
        self.basis = basis
        self.outfile = outfile
        self.frag_A = data.frag_A
        self.frag_B = data.frag_B
        self.do_fsapt = data.do_fsapt
        self.keywords = data.keywords
        self.step_size = data.irc_stepsize
        self.coordinates = data.coordinates
        self.force_min = data.force_min
        self.nthreads = data.nthreads
        self.set_memory = data.set_memory
        self.memory_allocation = data.memory_allocation
        # Specifics for F-SAPT Calculations
        if(data.do_fsapt):
            self.monomer_A_frags = data.monomer_A_frags
            self.monomer_B_frags = data.monomer_B_frags
            self.monomer_A_labels = data.monomer_A_labels
            self.monomer_B_labels = data.monomer_B_labels

    def psi4_sapt(self):
        self.int_  = []
        self.elst_ = []
        self.exch_ = []
        self.ind_ =  []
        self.disp_ = []
        all_atoms = self.frag_A + self.frag_B
        au_to_kcal = 627.51
        count = 0
        output = open(self.outfile, "a")
        output.write('\n\n--SAPT Decomposition--\n')
        output.write('\n-------------------------------------------------------------------------------------------------')
        output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('IRC Point', 'E_int(kcal)', 'E_elst(kcal)','E_exch(kcal)', 'E_ind(kcal)', 'E_disp(kcal)'))
        output.write('-------------------------------------------------------------------------------------------------\n')
        output.close()
        scale_sapt = False
        # See if the user is requesting a scaled SAPT0 calculation
        if(self.keywords['ssapt0_scale'] or self.keywords['SSAPT0_SCALE']):
            scale_sapt = True
        for geometry in self.geometries:
            output = open(self.outfile, "a")
            psi4.core.set_output_file("psi4_output/irc_%d_sapt.out" %count, False)
            geometry += "symmetry c1"
            psi4.geometry(geometry)
            # Adds fsapt capabilities
            if(self.do_fsapt):
                try:
                    os.mkdir('fsapt_output')
                except FileExistsError:
                    pass
                try:
                    os.mkdir('s-fsapt_output')
                except FileExistsError:
                    pass
                self.keywords['FISAPT_FSAPT_FILEPATH'] = 'fsapt_output/fsapt%d' %count
                if(scale_sapt):
                    self.keywords['FISAPT_FSSAPT_FILEPATH'] = 's-fsapt_output/fsapt%d' %count
            psi4.set_options(self.keywords)
            #print("pyREX:SAPT Calculation on IRC Point %d" %(count))
            psi4.set_num_threads(self.nthreads)
            if(self.set_memory):
                psi4.set_memory(self.memory_allocation)
            e_int = psi4.energy(self.method)
            #TODO: Actually have this create the SAPT files necessary for FSAPT partitioning. Make user options as well
            #      for fragment definitions. Something in the sapt block.
            if(self.do_fsapt):
                fsaptA_outfile = open('fsapt_output/fsapt%d/fA.dat' %count, 'w+')
                if(scale_sapt):
                    s_fsaptA_outfile = open('s-fsapt_output/fsapt%d/fA.dat' %count, 'w+')
                for i in range(len(self.monomer_A_labels)):
                    fsaptA_outfile.write("%s" %self.monomer_A_labels[i])
                    if(scale_sapt):
                        s_fsaptA_outfile.write("%s" %self.monomer_A_labels[i])
                    for j in range(len(self.monomer_A_frags[i])):
                        if(j==len(self.monomer_A_frags[i])-1):
                            fsaptA_outfile.write(" %d \n" %(all_atoms.index(self.monomer_A_frags[i][j])+1))
                            if(scale_sapt):
                                s_fsaptA_outfile.write(" %d \n" %(all_atoms.index(self.monomer_A_frags[i][j])+1))
                        else:
                            fsaptA_outfile.write(" %d " %(all_atoms.index(self.monomer_A_frags[i][j])+1))
                            if(scale_sapt):
                                s_fsaptA_outfile.write(" %d " %(all_atoms.index(self.monomer_A_frags[i][j])+1))
                fsaptA_outfile.close()
                s_fsaptA_outfile.close()
                fsaptB_outfile = open('fsapt_output/fsapt%d/fB.dat' %count, 'w+')
                if(scale_sapt):
                    s_fsaptB_outfile = open('s-fsapt_output/fsapt%d/fB.dat' %count, 'w+')
                for i in range(len(self.monomer_B_labels)):
                    fsaptB_outfile.write("%s" %self.monomer_B_labels[i])
                    if(scale_sapt):
                        s_fsaptB_outfile.write("%s" %self.monomer_B_labels[i])
                    for j in range(len(self.monomer_B_frags[i])):
                        if(j==len(self.monomer_B_frags[i])-1):
                            fsaptB_outfile.write(" %d \n" %(all_atoms.index(self.monomer_B_frags[i][j])+1))
                            if(scale_sapt):
                                s_fsaptB_outfile.write(" %d \n" %(all_atoms.index(self.monomer_B_frags[i][j])+1))
                        else:
                            fsaptB_outfile.write(" %d " %(all_atoms.index(self.monomer_B_frags[i][j])+1))
                            if(scale_sapt):
                                s_fsaptB_outfile.write(" %d " %(all_atoms.index(self.monomer_B_frags[i][j])+1))
                fsaptB_outfile.close()
                s_fsaptB_outfile.close()
            e_int =  psi4.core.variable("SAPT TOTAL ENERGY")
            e_elst = psi4.core.variable("SAPT ELST ENERGY")
            e_exch = psi4.core.variable("SAPT EXCH ENERGY")
            e_ind  = psi4.core.variable("SAPT IND ENERGY")
            e_disp = psi4.core.variable("SAPT DISP ENERGY")

            self.int_.append(e_int)
            self.elst_.append(e_elst)
            self.exch_.append(e_exch)
            self.ind_.append(e_ind)
            self.disp_.append(e_disp)
            output.write('\n{:>15} {:>15.4f} {:>15.4f} {:>15.4f} {:>15.4f} {:>15.4f}\n'.format(count, e_int*au_to_kcal, e_elst*au_to_kcal, e_exch*au_to_kcal, e_ind*au_to_kcal, e_disp*au_to_kcal))
            count = count+1
        output.write('-------------------------------------------------------------------------------------\n')
        output.close()
        return self.int_, self.elst_, self.exch_, self.ind_, self.disp_

    def psi4_super(self):
        self.int_  = []
        all_atoms = self.frag_A + self.frag_B
        au_to_kcal = 627.51
        count = 0
        output = open(self.outfile, "a")
        output.write('\n\n--Supermolecular Interaction Energy--\n')
        output.write('\n-------------------------------------------------------------------------------------------------')
        output.write('\n{:>15} {:>15}\n'.format('IRC Point', 'E_int(kcal)'))
        output.write('-------------------------------------------------------------------------------------------------\n')
        output.close()
        for geometry in self.geometries:
            output = open(self.outfile, "a")
            psi4.core.set_output_file("psi4_output/irc_%d_sapt.out" %count, False)
            geometry += "symmetry c1"
            psi4.geometry(geometry)
            psi4.set_options({'reference': 'rhf', 'basis' : self.basis})
            #print("pyREX:SAPT Calculation on IRC Point %d" %(count))
            e_int = psi4.energy("%s/%s" %(self.method,self.basis), bsse_type='cp')
            output.write('\n{:>15} {:>15.4f}\n'.format(count, e_int*au_to_kcal))
            output.close()
            self.int_.append(e_int)
            count = count + 1
        return self.int_

    def fsapt_post_analysis(self):
        fsapt_script_dir = os.environ['PSIPATH'][:-1] + '/share/psi4/fsapt/fsapt.py'
        link_type = "50-50"
        step_size = self.step_size
        coordinates = self.coordinates
        force_min = self.force_min
        self.monomer_A_labels.append("All")
        self.monomer_B_labels.append("All")
        unique_pairs_list = list(itertools.product(self.monomer_A_labels,self.monomer_B_labels))
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
        os.chdir("fsapt_output")
        dirs_list = [name for name in os.listdir() if os.path.isdir(name)] 
        num_geoms = 0
        for directory in dirs_list:
            if("fsapt" in directory):
                num_geoms = num_geoms + 1
        print(num_geoms)
        for i in range(num_geoms):
            fsapt_dir = "fsapt%d" %i
            os.chdir(fsapt_dir)
            print("running in %s" %fsapt_dir)
            subprocess.call(['python', fsapt_script_dir])
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
                  w_dicts_2[f_dict_count][pair] = -1.0*np.trapz(dict_[pair][index_min-1:],dx=step_size)
            f_dict_count = f_dict_count + 1
        
        work_values = open("work_values.dat", "w+")
        work_values.write('\n\n--Reaction Work Decomposition (Region 1)--\n')
        work_values.write('\n-----------------------------------------------------------------------------------------------------------------')
        work_values.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('Pair(A-B)','W_elst','W_exch',
        'W_indAB', 'W_indBA', 'W_disp', 'W_total'))
        work_values.write('-----------------------------------------------------------------------------------------------------------------\n')
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
