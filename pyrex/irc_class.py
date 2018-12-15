"""
Class for calculating IRC using Morokuma Algorithm
"""


import psi4
import json
import numpy


    
class Toolkit():
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
        self.tolerance=1.0e-04
        self.orcacmd='UNDEFINED'
        self.ReadInput(name)
        self.ReadHessian()
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
            JSON2XYZ(input_params)
        elif input_params['model']['basis']:
            self.basis = input_params['model']['basis']
            self.level_of_theory = "%s/%s" %(self.method,self.basis)
        elif input_params['model']['method']:
            self.method = input_params['model']['method']
            self.level_of_theory = "%s/%s" %(self.method,self.basis)
        elif input_params['molecule']['symbols']: 
            self.natoms = len(input_params['molecule']['symbols'])
        elif input_params['irc']['hessfile']:
            self.hessfn = input_params['irc']['hessfile']
        elif input_params['irc']['guess']:
            self.guessfn = input_params['irc']['guess']
        elif input_params['irc']['alpha']:
            self.alpha = input_params['irc']['alpha']
        elif input_params['irc']['delta']:
            self.delta = input_params['irc']['delta']
        elif input_params['irc']['damp']:
            self.damp = input_params['irc']['damp']
        elif input_params['irc']['restart']:
            if(input_params['irc']['restart']==1):
                self.restart = True
            else:
                self.restart = False
        elif input_params['irc']['autodamp']:
            if(input_params['irc']['autodamp']==1):
                self.autodamp=True
            else:
                self.autodamp=False
        elif input_params['irc']['tol']:
            self.tolerance= input_params['irc']['tol']
        elif input_params['irc']['alg']:
            self.algorithm = input_params['irc']['alg']
        elif input_params['irc']['dir']:
            self.direction = input_params['irc']['dir']
        elif input_params['irc']['mode']:
            self.mode = input_params['irc']['mode']
        elif input_params['irc']['maxd']:
            self.maxdispl = input_params['irc']['maxd']
        elif input_params['irc']['pts']:
            self.npoints = input_params['irc']['pts']
    def JSON2XYZ(self, input_params):
        geom = ''
        geom += '%d ' %input_params['molecule']['molecular_charge']
        geom += '%d \n' %input_params['molecule']['molecular_multiplicity'] 
        symbols = input_params['molecule']['symbols']
        geometry = input_params['molecule']['geometry']
        for i in range(len(symbols)):
            geom += "%s  " %symbols[i]
            for j in range(3):
                if(j==2):
                    geom += "   %f   \n" %geometry[(i*3) + j]
                else:
                    geom += "   %f   " %geometry[(i*3) + j]
        self.geom = geom
			

	def ComputeHessian(self):
		# Use PSI4 to calculate Hessian matrix
        psi4.geometry(self.geom)
		H = psi4.hessian(self.level_of_theory)
        Hess = np.array(H)
        self.displacement = Hess[:,self.mode]
        self.grad = np.zeros(3*self.natoms) 


##########################
# Energies and Gradients #
##########################

def doEnergy(geom):
    psi4.geometry(geom)
    E = psi4.energy(self.level_of_theory)
    return E

def doGrad(geom):
    psi4.geometry(geom)
    G = gradient(self.level_of_theory)
    return G
