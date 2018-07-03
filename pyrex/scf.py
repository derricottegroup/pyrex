import numpy as np
import psi4

def frag_opt_new(geometries,level_of_theory, outfile, natoms_A, natoms_B):
    frag_A_geom = ""
    geom_string = geometries[0][1]
    geom_string = geom_string.split('\n')[:(natoms_A+2)]
    for i in range(natoms_A+2):
        frag_A_geom += "%s\n" %geom_string[i]
    frag_B_geom = ""
    geom_string = geometries[0][2]
    geom_string = geom_string.split('\n')[:(natoms_B+2)]
    for i in range(natoms_B+2):
        frag_B_geom += "%s\n" %geom_string[i]
    output = open(outfile, "a")
    # Optimize Fragment A
    output.write("\n\n--Fragment A Geometry Optimization--\n\n")
    print("%s Geometry Optimization on Fragment A" %(level_of_theory))
    output.write("\nInitial Geometry of Fragment A:\n")
    output.write(frag_A_geom[frag_A_geom.find('\n')+4:])
    frag_A_geom += "symmetry c1"
    psi4.set_options({'reference': 'rhf'})
    psi4.geometry(frag_A_geom)
    psidump = "psi4_output/fragment_A_opt.out"
    psi4.core.set_output_file(psidump, False)
    e_A = psi4.optimize(level_of_theory)
    output.write("\nOutput Written To: %s\n" %psidump)
    output.write("Final Energy: %f\n" %e_A)
    # Optimize Fragment B
    output.write("\n\n--Fragment B Geometry Optimization--\n\n")
    print("%s Geometry Optimization on Fragment B" %(level_of_theory))
    output.write("\nInitial Geometry of Fragment B:\n")
    output.write(frag_B_geom[frag_B_geom.find('\n')+4:])
    frag_B_geom += "symmetry c1"
    psi4.set_options({'reference': 'rhf'})
    psi4.geometry(frag_B_geom)
    psidump = "psi4_output/fragment_B_opt.out"
    psi4.core.set_output_file(psidump, False)
    e_B = psi4.optimize(level_of_theory)
    output.write("\nOutput Written To: %s\n" %psidump)
    output.write("Final Energy: %f\n" %e_B)
    output.close()
    return e_A , e_B

def frag_opt(geometry_frag, level_of_theory, outfile, label):
    output = open(outfile, "a")
    output.write("\n\n--Fragment %s Geometry Optimization--\n\n" %label)
    print("%s Geometry Optimization on Fragment %s" %(level_of_theory, label))
    output.write("\nInitial Geometry of Fragment:\n")
    output.write(geometry_frag[geometry_frag.find('\n')+4:])
    geometry_frag += "symmetry c1"
    psi4.set_options({'reference': 'rhf'})
    psi4.geometry(geometry_frag)
    psidump = "psi4_output/fragment_%s_opt.out" %label
    psi4.core.set_output_file(psidump, False)
    e_frag = psi4.optimize(level_of_theory)
    output.write("\nOutput Written To: %s\n" %psidump)
    output.write("Final Energy: %f\n" %e_frag)
    output.close()
    return e_frag

def psi4_scf(geometries, level_of_theory, pol=False, do_eda=False):
    energies = [] # Energy array of tuples -> (Dimer_Energy, Fragment_A_Energy, Fragment_B_Energy)
    wavefunctions = []
    interaction_energies = []
    for i in range(len(geometries)):
        psi4.core.set_output_file("psi4_output/irc_%d.out" %i, False)
        geometry = geometries[i][0]
        geometry += "symmetry c1"
        psi4.geometry(geometry)
        #TODO Give the user control over these SCF options
        psi4.set_options({'reference': 'rhf', 'save_jk' : 'true'})
        print("pyREX:Single Point Calculation on IRC Point %d" %(i))
        dimer_energy, dimer_wfn = psi4.energy(level_of_theory, return_wfn=True)
        #e_1e = psi4.core.get_variable("ONE-ELECTRON ENERGY")
        #e_2e = psi4.core.get_variable("TWO-ELECTRON ENERGY")
        e_nuc = psi4.core.get_variable("NUCLEAR REPULSION ENERGY")
        if(do_eda==True):
            j = np.array(dimer_wfn.jk().J()[0].to_array())
            k = np.array(dimer_wfn.jk().K()[0].to_array())
            h = np.array(dimer_wfn.H().to_array()) 
            c = np.array(dimer_wfn.Ca_subset("AO", "OCC").to_array())
            D = np.dot(c, c.transpose())
        if(pol==True):
            # Fragment A SCF
            psi4.core.set_output_file("psi4_output/irc_%d_A_scf.out" %i, False)
            geometry_A = geometries[i][1]
            geometry_A += "symmetry c1"
            psi4.geometry(geometry_A)
            psi4.set_options({'reference': 'rhf', 'save_jk' : 'true'})
            print("pyREX:Single Point Calculation on IRC Point %d (Fragment A)" %(i))
            frag_A_energy, frag_A_wfn = psi4.energy(level_of_theory, return_wfn=True)
            #frag_A_e_1e = psi4.core.get_variable("ONE-ELECTRON ENERGY")
            #frag_A_e_2e = psi4.core.get_variable("TWO-ELECTRON ENERGY")
            frag_A_e_nuc = psi4.core.get_variable("NUCLEAR REPULSION ENERGY")
            if(do_eda==True):
                k_A = np.array(frag_A_wfn.jk().K()[0].to_array())
                c_A = np.array(frag_A_wfn.Ca_subset("AO", "OCC").to_array())
                D_A = np.dot(c_A, c_A.transpose())
            # Fragment B SCF
            psi4.core.set_output_file("psi4_output/irc_%d_B_scf.out" %i, False)
            geometry_B = geometries[i][2]
            geometry_B += "symmetry c1"
            psi4.geometry(geometry_B)
            psi4.set_options({'reference': 'rhf', 'save_jk' : 'true'})
            print("pyREX:Single Point Calculation on IRC Point %d (Fragment B)" %(i))
            frag_B_energy, frag_B_wfn = psi4.energy(level_of_theory, return_wfn=True)
            #frag_B_e_1e = psi4.core.get_variable("ONE-ELECTRON ENERGY")
            #frag_B_e_2e = psi4.core.get_variable("TWO-ELECTRON ENERGY")
            frag_B_e_nuc = psi4.core.get_variable("NUCLEAR REPULSION ENERGY")
            interaction_energy = dimer_energy - (frag_A_energy + frag_B_energy)
            if(do_eda==True):
                k_B = np.array(frag_B_wfn.jk().K()[0].to_array())
                c_B = np.array(frag_B_wfn.Ca_subset("AO", "OCC").to_array())
                D_B = np.dot(c_B, c_B.transpose())
                approx_energy_1 = 2.0*np.vdot(D , h) + 2.0*np.vdot(D , j) - (np.vdot(D_A, k_A) + np.vdot(D_B, k_B)) + e_nuc
                approx_energy_2 = 2.0*np.vdot(D , h) + 2.0*np.vdot(D , j) - np.vdot(D, k) + e_nuc
                #interaction_energy = dimer_energy - (frag_A_energy + frag_B_energy)
                electrostatic = (frag_A_energy  + frag_B_energy) - approx_energy_1
                exchange = approx_energy_1 - approx_energy_2
                relaxation = interaction_energy - electrostatic - exchange
        else:
            frag_A_energy = 0.0
            frag_A_wfn = 0.0
            frag_B_energy = 0.0
            frag_B_wfn = 0.0
            interaction_energy = 0.0
            electrostatic = 0.0
            exchange = 0.0
            relaxation = 0.0
        energies.append((dimer_energy,frag_A_energy,frag_B_energy))
        wavefunctions.append((dimer_wfn, frag_A_wfn, frag_B_wfn))
        if(do_eda==True):
            interaction_energies.append((interaction_energy, electrostatic,exchange,relaxation))
    return energies, wavefunctions, interaction_energies        
