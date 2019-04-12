import os
import psi4
import numpy as np


class sapt(object):

    def __init__(self, geometries, method, basis, outfile):
        self.geometries = geometries
        self.method = method
        self.basis = basis
        self.outfile = outfile

    def psi4_sapt(self):
        self.int_  = []
        self.elst_ = []
        self.exch_ = []
        self.ind_ =  []
        self.disp_ = []
        au_to_kcal = 627.51
        count = 0
        output = open(self.outfile, "a")
        output.write('\n\n--SAPT Decomposition--\n')
        output.write('\n-------------------------------------------------------------------------------------------------')
        output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('IRC Point', 'E_int(kcal)', 'E_elst(kcal)','E_exch(kcal)', 'E_ind(kcal)', 'E_disp(kcal)'))
        output.write('-------------------------------------------------------------------------------------------------\n')
        output.close()
        for geometry in self.geometries:
            output = open(self.outfile, "a")
            psi4.core.set_output_file("psi4_output/irc_%d_sapt.out" %count, False)
            geometry += "symmetry c1"
            psi4.geometry(geometry)
            # Adds fsapt capabilities 
            try:
                os.mkdir('fsapt_output')
            except FileExistsError:
                pass
            psi4.set_options({'reference': 'rhf', 'basis' : self.basis, 'exch_scale_alpha' : 'true', 'FISAPT_FSAPT_FILEPATH' :'fsapt_output/fsapt%d' %count})
            #print("pyREX:SAPT Calculation on IRC Point %d" %(count))
            e_int = psi4.energy(self.method)
            #TODO: Actually have this create the SAPT files necessary for FSAPT partitioning. Make user options as well
            #      for fragment definitions. Something in the sapt block. 
            fsaptA_outfile = open('fsapt_output/fsapt%d/fA.dat' %count, 'w+')
            fsaptB_outfile = open('fsapt_output/fsapt%d/fB.dat' %count, 'w+')
            e_elst = psi4.core.get_variable("SAPT ELST ENERGY")
            e_exch = psi4.core.get_variable("SAPT EXCH ENERGY")
            e_ind  = psi4.core.get_variable("SAPT IND ENERGY")
            e_disp = psi4.core.get_variable("SAPT DISP ENERGY")

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
