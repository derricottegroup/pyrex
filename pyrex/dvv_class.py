"""
Class for calculating IRC using Morokuma Algorithm
"""

import sys
import os
import psi4
import json
import md_helper
import numpy as np

# Global Constants (Atomic Units conversion)
fs_timeau = 41.34137314
amu2au = 1822.8884850
bohr2ang = 0.529177

def JSON2XYZ(input_params):
    charge = input_params['molecule']['molecular_charge']
    mult = input_params['molecule']['molecular_multiplicity']
    geom = ''
    geom += '\n%d %d\n' %(charge, mult)
    symbols = input_params['molecule']['symbols']
    geometry = input_params['molecule']['geometry']
    for i in range(len(symbols)):
        geom += "%s  " %symbols[i]
        for j in range(3):
            if(j==2):
                geom += "   %f   \n" %geometry[(i + 2*i) + j]
            else:
                geom += "   %f   " %geometry[(i + 2*i) + j]
    #print(geom)
    return geom

def normal_mode_reader(input_params):
    natoms = len(input_params['molecule']['symbols'])
    normal_mode_file = open(input_params['dvv']['normal_mode_file'], 'r')
    ts_vec = []
    for line in normal_mode_file:
        if("vibration %d" %input_params['dvv']['mode'] in line):
            print(line)
            line = next(normal_mode_file)
            print(line)
            for i in range(natoms):
                trj = line.split()
                print(trj)
                trj = list(map(float, trj))
                np_trj = np.array(trj)
                print(np_trj)
                #print(bohr2ang)
                #np_trj_ang = bohr2ang*np_trj
                #print(np_trj_ang)
                ts_vec.append(np_trj)
                line = next(normal_mode_file)
    print(ts_vec)
    return ts_vec
    
class ToolKit():
    # Placeholder for Common Data
    def __init__(self,name):
        self.natoms=0
        self.basis = '3-21G'
        self.method = 'SCF'
        self.level_of_theory = "%s/%s" %(self.method,self.basis)
        self.alpha=0.1
        self.basename='.'.join(name.split('.')[:-1])
        self.out=open("%s.log"%(self.basename),'a',1)
        self.delta=0.05
        self.direction=1
        self.energies=[]
        self.energy=0.00
        self.geos=[]
        self.restart=False
        self.geometry=[]
        self.guessfn=''
        self.damp=0.05
        self.algorithm=1
        self.autodamp=False
        self.prevgrad=0.0
        self.hessfn=''
        self.maxdispl=0.01
        self.mode=0
        self.dgrad=0.0
        self.npoints=25
        self.template=[]
        self.keywords={}
        self.tolerance=1.0e-04
        self.orcacmd='UNDEFINED'
        self.ReadInput(name)
        self.ComputeHessian()
        self.cons_vel = 0.04
        self.err_tol = 0.003
        #self.ts_vec = []
        self.normal_mode_file=''
        #self.symbols=[]
    def printPars(self):
        stmp="  %13s: %s\n"
        ftmp="  %13s: %7.4f\n"
        itmp="  %13s: %7d\n"
        self.out.write('')
        self.out.write(itmp%('Algorithm',self.algorithm))
        self.out.write(itmp%('N. Points',self.npoints))
        self.out.write(ftmp%('Grad. Tol.',self.tolerance))
        self.out.write('')
        self.out.write(stmp%('Hessian',self.hessfn))
        self.out.write(itmp%('Mode',self.mode))
        self.out.write(itmp%('Direction',self.direction))
        self.out.write('')
        self.out.write(ftmp%('Alpha',self.alpha))
        self.out.write(ftmp%('Delta',self.delta))
        self.out.write('')
        self.out.write(ftmp%('Damp Factor',self.damp))
        self.out.write(itmp%('Damp Update',self.autodamp))
        self.out.write('')
        self.out.write(ftmp%('Max. Displ.',self.maxdispl))
        self.out.write('')
        self.out.write(stmp%('Guess',self.guessfn))
        self.out.write("\n------------------------------------------------\n")
    def ReadInput(self,json_input):
        # Read in input from JSON input file
        json_data=open(json_input).read()
        input_params = json.loads(json_data)
        if input_params['molecule']['geometry']:
            self.geometry = JSON2XYZ(input_params)
        if input_params['model']['basis']:
            self.basis = input_params['model']['basis']
            self.level_of_theory = "%s/%s" %(self.method,self.basis)
            #print("IN HERE")
        if input_params['model']['method']:
            self.method = input_params['model']['method']
            self.level_of_theory = "%s/%s" %(self.method,self.basis)
        if input_params['molecule']['symbols']:
            self.symbols = input_params['molecule']['symbols'] 
            self.natoms = len(input_params['molecule']['symbols'])
            print(self.symbols)
            print(self.natoms)
        if input_params['keywords']:
            self.keywords = input_params['keywords']
        if input_params['dvv']['cons_vel']:
            self.cons_vel = input_params['dvv']['cons_vel']
        if input_params['dvv']['err_tol']:
            self.err_tol = input_params['dvv']['err_tol']
        if input_params['dvv']['mode']:
            self.mode = input_params['dvv']['mode']
        if input_params['dvv']['normal_mode_file']:
            self.normal_mode_file = input_params['dvv']['normal_mode_file']
            print(input_params['dvv']['normal_mode_file'])
            self.ts_vec = normal_mode_reader(input_params)
            print(self.ts_vec) 
        if "irc" in input_params:
            if input_params['irc']['hessfile']:
                self.hessfn = input_params['irc']['hessfile']
            if input_params['irc']['guess']:
                self.guessfn = input_params['irc']['guess']
            if input_params['irc']['alpha']:
                self.alpha = input_params['irc']['alpha']
            if input_params['irc']['delta']:
                self.delta = input_params['irc']['delta']
            if input_params['irc']['damp']:
                self.damp = input_params['irc']['damp']
            if input_params['irc']['restart']:
                if(input_params['irc']['restart']==1):
                    self.restart = True
                else:
                    self.restart = False
            if input_params['irc']['autodamp']:
                if(input_params['irc']['autodamp']==1):
                    self.autodamp=True
                else:
                    self.autodamp=False
            if input_params['irc']['tol']:
                self.tolerance= input_params['irc']['tol']
            if input_params['irc']['alg']:
                self.algorithm = input_params['irc']['alg']
            if input_params['irc']['dir']:
                self.direction = input_params['irc']['dir']
            if input_params['irc']['mode']:
                self.mode = input_params['irc']['mode']
            if input_params['irc']['maxd']:
                self.maxdispl = input_params['irc']['maxd']
            if input_params['irc']['pts']:
                self.npoints = input_params['irc']['pts']
        else: 
            pass
   # def JSON2XYZ(self, input_params):
   #     geom = ''
   #     geom += '%d ' %input_params['molecule']['molecular_charge']
   #     geom += '%d \n' %input_params['molecule']['molecular_multiplicity'] 
   #     symbols = input_params['molecule']['symbols']
   #     geometry = input_params['molecule']['geometry']
   #     for i in range(len(symbols)):
   #         geom += "%s  " %symbols[i]
   #         for j in range(3):
   #             if(j==2):
   #                 geom += "   %f   \n" %geometry[(i + 2*i) + j]
   #             else:
   #                 geom += "   %f   " %geometry[(i + 2*i) + j]
   #     self.geom = geom	

    def ComputeHessian(self):
        # Use PSI4 to calculate Hessian matrix
        psi4.core.set_output_file("hessian.out", False)
        psi4.geometry(self.geometry)
        #print(self.keywords)
        psi4.set_options(self.keywords)
        H = psi4.hessian(self.method)
        Hess = np.array(H)
        print(Hess.shape)
        self.displacement = Hess[:,self.mode]
        print(self.displacement.shape)
        self.grad = np.zeros(3*self.natoms) 


##########################
# Energies and Gradients #
##########################

def doEnergy(geom,pars):
    psi4.geometry(geom)
    E = psi4.energy(pars.level_of_theory)
    return E

def doGrad(geom,pars):
    psi4.geometry(geom)
    Grad, wfn = psi4.gradient(pars.level_of_theory, return_wfn=True)
    e = wfn.energy()
    G = np.array(Grad)
    print(G.shape)
    pars.grad = G
    return G, e


######################################
# Geometry manipulation and printing #
######################################

def geodisplace(geom, displace):
    new_geom = []
    displaced_geom = ''
    symbols = []
    #print(geom)
    lines = geom.split('\n')
    lines.pop(0)
    del lines[-1]
    charge_mult_line = lines[0].split() #Extract the Charge and Multiplicity from top of string
    charge = int(charge_mult_line[0])
    mult = int(charge_mult_line[1])
    for line in lines[1:]: # Loop Over lines that actually contain atoms
        atom_line = line.split()
        #print(atom_line)
        symbols.append(atom_line[0])
        atom_line.pop(0) # Get rid of atom symbol temporarily
        new_geom.append(list(map(float,atom_line)))
    i=0
    for atom in new_geom:
        atom[0] += displace[i]
        atom[1] += displace[i+1]
        atom[2] += displace[i+2]
        i += 3
#TODO (12/15/18) Finish this, all we have to do is transform the geometry from a list to a
#string like the geometry given as an input to this function.
    displaced_geom += "\n%d %d\n" %(charge, mult)
    for i in range(len(symbols)):
        displaced_geom += "%s  " %symbols[i]
        for j in range(3):
            if(j==2):
                #displaced_geom += "   %f   \n" %new_geom[(i + 2*i) + j]
                displaced_geom += "   %f   \n" %new_geom[i][j]
            else:
                #displaced_geom += "   %f   " %new_geom[(i + 2*i) + j]
                displaced_geom += "   %f   " %new_geom[i][j]
    print(displaced_geom)
    return displaced_geom

def printTrj(params, n):
    trajectory_file = open(params.basename+'.trj', 'a')
    trajectory_file.write("%d\npyREX IRC point %d E=%14.7f\n"%(params.natoms,n,params.energy))
    trajectory_file.write("%s\n" %params.geometry)
    trajectory_file.close()

########################
#Damped Velocity Verlet#
########################

def md_main(params):
#MD Options
    timestep =  1.0                       # Time step for each iteration in time atomic units
    max_md_step = 400                 # Number of MD iterations
    #veloc0 = np.zeros((2,3))            # Numpy array (natoms,3) with inital velocities
    trajec = True                       # Boolean: Save all trajectories in a single xyz file 
    int_alg = 'veloc_verlet'            # Algorithm to use as integrator
    opt    = False                      # Optimize geometry using F=ma
    E = doEnergy(params.geometry, params)
    geom = params.geometry
    # Initial positions, velocities, accelerations, and forces 
    mol = psi4.geometry(geom)
    energy, forces = md_helper.get_forces(params.level_of_theory)
    natoms = mol.natom()
    atom_mass = np.asarray([mol.mass(atom) for atom in range(natoms)])*amu2au
    veloc0 = np.zeros((natoms,3))            # Numpy array (natoms,3) with inital velocities
    veloc = params.ts_vec
    print(veloc0)
    print(veloc)
    accel = forces/(atom_mass.reshape((natoms,1)))

    # MD Loop
    pos = np.asarray(mol.geometry())

    # Save energy of each iteration on md_energy file
    md_energy = open('md_energy.dat','w')
    md_energy.write('File with the energy of each MD trajectory point\n\n')
    md_energy.write('Trajectory Number\tEnergy (Hartree)\n')
    md_energy.write('{0:>3d}\t\t\t{1:10.8f}\n'.format(1,energy))
    pos_vec = []
    vel_vec = []
    accel_vec = []
    energy_vec = []
    time_vec = []
    for i in range(1,max_md_step+1):
        # Saving energies and trajectory points
        md_energy.write('{0:>3d}\t\t\t{1:10.8f}\n'.format(i,energy))
        if trajec:
            mol.save_xyz_file('md_step_'+str(i)+'.xyz',False)
    
        # Updating positions velocities and accelerations using Velocity Verlet Integrator
        pos_new,vel_new,accel_new,energy_new = md_helper.integrator(int_alg,timestep,pos,veloc,accel,mol,params.level_of_theory,params.symbols)
        if(i>3):
            # Compute Displacement point x(prime) eq 6. in paper
            print(len(pos_vec))
            print(len(vel_vec))
            print(len(accel_vec))
            print(len(time_vec))
            i_minus_two = i-3 #Had to offset these since the loop starts at 1. 
            i_minus_one = i-2 #Had to offset these since the loop starts at 1.
            print(i-2)
            print(i-1)
            pos_disp = pos_vec[i_minus_two] + vel_vec[i_minus_two]*(time_vec[i_minus_one]+ timestep) + 0.5*accel_vec[i_minus_two]*((time_vec[i_minus_one] + timestep)**2.0)
            disp = pos_disp - pos # Error Estimate
            disp_norm = np.linalg.norm(disp)
            disp_matrix = np.asarray(disp)
            max_disp = disp_matrix.max()
            err_est = 0.0
            if(disp_norm > max_disp):
                err_est = disp_norm
            else:
                err_est = max_disp
            # Update timestep based on error estimation
            time_new = timestep*((params.err_tol/err_est)**(1.0/3.0))
            print("New Timestep = %f" %time_new)
            timestep = time_new
            pos = pos_new
            veloc = md_helper.damp_velocity(vel_new, params)
            accel = accel_new
            energy = energy_new
            pos_vec.append(pos)
            vel_vec.append(veloc)
            accel_vec.append(accel)
            energy_vec.append(energy)
            time_vec.append(timestep)
        else:
            pos = pos_new
            veloc = md_helper.damp_velocity(vel_new, params)
            accel = accel_new
            energy = energy_new
            pos_vec.append(pos)
            vel_vec.append(veloc)
            accel_vec.append(accel)
            energy_vec.append(energy)
            time_vec.append(timestep)
    md_energy.close()
    if trajec:
        md_helper.md_trajectories(max_md_step)
    print("Done with Molecular Dynamics Program!")
    print(energy)
    print(forces)
    print(natoms)

inpname = sys.argv[1]
params=ToolKit(inpname)
md_main(params)
