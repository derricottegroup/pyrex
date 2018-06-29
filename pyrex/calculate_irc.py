"""
Script for Calculating Energies Along IRC in Psi4
"""
__authors__ = "Wallace D. Derricotte"
__credits__ = ["Wallace D. Derricotte"]

__copyright__ = "(c) 2018, Derricotte Research Group"
__license__ = "BSD-3-Clause"
__date__ = "2018-5-29"

import psi4
import numpy as np
import os
import datetime

def num_first_derivs(irc_energies,step_size):
    reaction_force = []
    for i in range(2,len(irc_energies)-2):
        current_force = -1.0*(-1.0*irc_energies[i+2] + 8.0*irc_energies[i+1] - 8.0*irc_energies[i-1] + 1.0*irc_energies[i-2])/(12.0*step_size)
        reaction_force.append(current_force)
    return reaction_force

def num_integrate(force_coordinates, reaction_force_values, lower_limit, upper_limit):
    integral_value = 0.0
    for i in range(lower_limit,upper_limit):
        integral_value += (force_coordinates[i+1] - force_coordinates[i])*((reaction_force_values[i] + reaction_force_values[i+1])/2.0)
    return integral_value

header ='''
-----------------------------------------------------------------------
     pyREX: Python Reaction Energy eXtension for Quantum Chemistry
                            
                         pyREX beta 0.1

                      Wallace D. Derricotte
                    Derricotte Research Group
                        Morehouse College
                         Atlanta,Georgia
-----------------------------------------------------------------------
'''


full_irc = open("full_irc.xyz", "r")
output = open("edacity_output.dat", "w+")
csv_file = open("raw_data.csv","w+")
irc = []
pid = os.getpid()
output.write(header)
datetime_now = datetime.datetime.now().strftime('pyREX Run on %A %B %d, %Y at %H:%M%p')
output.write("\n")
output.write(datetime_now)
output.write("\n")
output.write("Process ID: %d" %pid)
output.close()

charge_A = 0 #Specify total charge on Monomer B
mult_A = 1 #Specify multiplicity on Monomer B
frag_A_atom_list = [0,1,4]
charge_B = 0 #Specify total charge on Monomer B
mult_B = 1 #Specify multiplicity on Monomer B
frag_B_atom_list = [2,3,5]

charge_dimer = 0 #Specify total charge on the supermolecular complex
mult_dimer = 1 #Specify multiplicity of the supermolecular complex

irc_step_size = 0.2 #in units au*amu^(1/2), Psi4 default is 0.2
level_of_theory = "scf/3-21G" # Level of Theory for Total Energies

# Grab number of atoms (natoms) from the top of the XYZ file.
natoms = int(full_irc.readline())
coordinates = []

# Grab and store geometries from the IRC
for line in full_irc:
    if "Full IRC Point" in line:
        geom = []
        irc_num_line = line.split()
        irc_num = int(irc_num_line[3])
        for i in range(natoms):
            line = next(full_irc)
            geom.append(line)
        irc.append((irc_num, geom))
        coordinates.append(irc_num*irc_step_size)

irc_energies = []
reaction_force = []
strain_energies = []
chemical_potentials = []
chemical_potentials_A = []
chemical_potentials_B = []
reaction_electronic_flux = []
reaction_electronic_flux_A = []
reaction_electronic_flux_B = []

e_A = 0.0
e_B = 0.0

for i in range(len(irc)):
    if (i==0):
        output = open("edacity_output.dat", "a")
        geom_A_header ='''
----------------------------------
Geometry Optimization On Monomer A
----------------------------------
'''
        output.write(geom_A_header)
        print("%s Geometry Optimization on Monomer A" %level_of_theory)
        #geometry_A = ""
        geometry_A = "\n%d %d\n" %(charge_A, mult_A)
        for j in range(len(frag_A_atom_list)):
            geometry_A += irc[i][1][frag_A_atom_list[j]]
        output.write("\nInitial Geometry of Monomer A:\n")
        output.write(geometry_A[geometry_A.find('\n')+4:])
        geometry_A += "symmetry c1" 
        psi4.set_options({'reference': 'rhf'})
        #print(geometry_A)
        psi4.geometry(geometry_A)
        outfile = "monomer_A_opt.out"
        psi4.core.set_output_file(outfile, False)
        e_A = psi4.optimize(level_of_theory)
        output.write("\nOutput Written To: %s\n" %outfile)
        output.write("Final Energy: %f" %e_A)
        output.close()
        #print("%s Geometry Optimization on Monomer B" %level_of_theory)
        geom_B_header ='''
----------------------------------
Geometry Optimization On Monomer B
----------------------------------
'''
        output = open("edacity_output.dat", "a")
        output.write(geom_B_header)
        #geometry_B = ""
        geometry_B = "\n%d %d\n" %(charge_B, mult_B)
        for j in range(len(frag_B_atom_list)):
            geometry_B += irc[i][1][frag_B_atom_list[j]]
        psi4.set_options({'reference': 'rhf'})
        print(geometry_B)
        output.write("\nInitial Geometry of Monomer B:\n")
        output.write(geometry_B[geometry_B.find('\n')+4:])
        geometry_B += "symmetry c1"
        psi4.geometry(geometry_B)
        outfile = "monomer_B_opt.out"
        psi4.core.set_output_file(outfile, False)
        e_B = psi4.optimize(level_of_theory)
        output.write("\nOutput Written To: %s\n" %outfile)
        output.write("Final Energy: %f" %e_B)
        output.close()
    geometry = ""
    geometry += "\n%d %d\n" %(charge_dimer, mult_dimer)
    for j in range(len(irc[i][1])):
        geometry += irc[i][1][j]
    output.write("=================================================\n")
    output.write("Single Point Energy Calculation on IRC Point %d\n" %irc[i][0])
    output.write("=================================================\n")
    psi4.core.set_output_file("irc_%d.out" %irc[i][0], False)
    psi4.geometry(geometry)
    output.write("Current Geometry:\n")
    for j in range(natoms):
        output.write("%s  \n" %irc[i][1][j])
    psi4.set_options({'reference': 'rhf'})
    print("%s Single Point Calculation on IRC Point %d" %(level_of_theory,irc[i][0]))
    current_energy, current_wfn = psi4.energy(level_of_theory, return_wfn=True)
    # Run SCF on Both Monomers for Flux Polarization Calculation (TODO Make this optional)
    geometry_A = ""
    geometry_A += "\n%d %d\n" %(charge_A, mult_A)
    for j in range(len(frag_A_atom_list)):
        geometry_A += irc[i][1][frag_A_atom_list[j]]
    geometry_A += "symmetry c1"
    psi4.geometry(geometry_A)
    psi4.core.set_output_file("irc_%d_A_scf.out" %irc[i][0], False)
    psi4.set_options({'reference': 'rhf','mom_start': '5'})
    #monomer_A_energy, monomer_A_wfn = psi4.energy(level_of_theory, return_wfn=True)
    # Grab number of occupied orbitals
    #ndocc_A = monomer_A_wfn.doccpi()[0]
    # Coefficient Matrix
    #C_A = np.array(monomer_A_wfn.Ca())
    # Orbital energies
    #eps_A = monomer_A_wfn.epsilon_a()
    #eps_A = np.array([eps_A.get(x) for x in range(C_A.shape[0])])
    #homo_energy_A = eps_A[ndocc_A-1]
    #lumo_energy_A = eps_A[ndocc_A]
    #chemical_potential_A = 0.5*(homo_energy_A + lumo_energy_A)
    #chemical_potentials_A.append(chemical_potential_A)
    #geometry_B = ""
    #geometry_B += "\n%d %d\n" %(charge_B, mult_B)
    #for j in range(len(frag_B_atom_list)):
    #    geometry_B += irc[i][1][frag_B_atom_list[j]]
    #geometry_B += "symmetry c1"
    #psi4.geometry(geometry_B)
    #psi4.core.set_output_file("irc_%d_B_scf.out" %irc[i][0], False)
    #psi4.set_options({'reference': 'rhf','mom_start': '5'})
    #monomer_B_energy, monomer_B_wfn = psi4.energy(level_of_theory, return_wfn=True)
    # Grab number of occupied orbitals
    #ndocc_B = monomer_B_wfn.doccpi()[0]
    #print(ndocc)
    # Coefficient Matrix
    #C_B = np.array(monomer_B_wfn.Ca())
    # Orbital energies
    #eps_B = monomer_B_wfn.epsilon_a()
    #eps_B = np.array([eps_B.get(x) for x in range(C_B.shape[0])])
    #homo_energy_B = eps_B[ndocc_B-1]
    #lumo_energy_B = eps_B[ndocc_B]
    #chemical_potential_B = 0.5*(homo_energy_B + lumo_energy_B)
    #chemical_potentials_B.append(chemical_potential_B)
    # Grab number of occupied orbitals
    ndocc = current_wfn.doccpi()[0]
    print(ndocc)
    # Coefficient Matrix
    C = np.array(current_wfn.Ca())
    # Orbital energies
    eps = current_wfn.epsilon_a()
    eps = np.array([eps.get(x) for x in range(C.shape[0])])
    homo_energy = eps[ndocc-1]
    lumo_energy = eps[ndocc]
    chemical_potential = 0.5*(homo_energy + lumo_energy)
    output.write("Energy = %.10f Hartree\n" %current_energy)
    strain_energy = current_energy - (e_A + e_B)
    #print("Strain Energy for Point %d = %.4f" %(irc[i][0],strain_energy))
    irc_energies.append(current_energy)
    strain_energies.append(strain_energy)
    chemical_potentials.append(chemical_potential)
force_coordinates = coordinates[2:len(coordinates)-2]
reaction_force_values = num_first_derivs(irc_energies,irc_step_size)
reaction_electronic_flux = num_first_derivs(chemical_potentials,irc_step_size)
reaction_electronic_flux_A = num_first_derivs(chemical_potentials_A,irc_step_size)
reaction_electronic_flux_B = num_first_derivs(chemical_potentials_B,irc_step_size)

output.write("=================================================\n")
output.write("          Edacity Reaction Force Analysis        \n")
output.write("               Wallace D. Derricotte             \n")
output.write("=================================================\n")
output.write("Minimum Force =  %.10f\n" %(min(reaction_force_values)))
output.write("Maximum Force =  %.10f\n" %(max(reaction_force_values)))


for i in range(len(reaction_force_values)):
    reaction_force.append((force_coordinates[i],reaction_force_values[i]))
print(reaction_force)

index_min = np.argmin(np.asarray(reaction_force_values))
index_ts  = force_coordinates.index(0.0000)
index_max = np.argmax(np.asarray(reaction_force_values))

output.write("\nReactant Region:          %.3f ------> %.3f\n" %(coordinates[0], force_coordinates[index_min]))
output.write("\nTransition State Region:  %.3f ------> %.3f\n" %(force_coordinates[index_min], force_coordinates[index_max]))
output.write("\nProduct Region:            %.3f ------> %.3f\n" %(force_coordinates[index_max], coordinates[-1]))

# Calculate Work in Reactant Region
W_1 = -1.0*num_integrate(force_coordinates, reaction_force_values, 0, index_min)

#Calculate Work in Transition State Region 1
W_2 = -1.0*num_integrate(force_coordinates, reaction_force_values, index_min, index_ts)

#Calculate Work in Transition State Region 2
W_3 = -1.0*num_integrate(force_coordinates, reaction_force_values, index_ts, index_max)

#Calculate Work in Product Region
W_4 = -1.0*num_integrate(force_coordinates, reaction_force_values, index_max, len(reaction_force)-1)

output.write("===============\n")
output.write("Work Integrals\n")
output.write("===============\n")
output.write("%f kcal/mol\n" %(W_1*627.51))
output.write("%f kcal/mol\n" %(W_2*627.51))
output.write("%f kcal/mol\n" %(W_3*627.51))
output.write("%f kcal/mol\n" %(W_4*627.51))

output.writelines(["%f," % energy  for energy in irc_energies])
output.write("\n\n")
output.writelines(["%f," % coord  for coord in coordinates])
output.write("\n\n")
output.writelines(["%f," % force  for force in reaction_force_values])


# Calculate Energy Difference
min_energy_index = np.argmin(np.asarray(irc_energies))
Del_E_raw = []
for i in range(len(irc_energies)):
    Del_E_raw.append(irc_energies[i] - irc_energies[min_energy_index])

Del_E = Del_E_raw[2:len(coordinates)-2]

#Create CSV File

output.write("----------------------------------------------------------------------------------------\n")
output.write("   Coordinate(au amu^(1/2))            Energy_diff(au)                  Force(au)    \n")
output.write("----------------------------------------------------------------------------------------\n")
for i in range(len(reaction_force_values)):
    output.write("       %.3f                           %.8f                             %.8f\n" %(force_coordinates[i], Del_E[i], reaction_force_values[i]))

potentials_truncated = chemical_potentials[2:len(coordinates)-2]
csv_file.write("Coordinate,DeltaE,Force,Chemical Potential,Reaction Electronic Flux, Flux A, Flux B\n")
for i in range(len(reaction_force_values)):
    csv_file.write("%f,%f,%f,%f,%f,%f,%f\n" %(force_coordinates[i], Del_E[i], reaction_force_values[i], potentials_truncated[i],reaction_electronic_flux[i], reaction_electronic_flux_A[i], reaction_electronic_flux_B[i]))
