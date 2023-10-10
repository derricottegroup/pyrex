import sys
from . import rex
import argparse
from . import input_gen

pyrex_help = """PYREX is an open-source python toolkit for intrinsic reactivity analysis"""
parser = argparse.ArgumentParser(description=pyrex_help)
def read_options():
    parser.add_argument('--input', type=str, help="Specify input filename",default=None)
    #parser.add_argument('--output', type=str, help="Specify output filename",default='pyrex_output.dat')
    parser.add_argument('--rexplot', action="store_true", help="Use PYREX plotting utility")
    parser.add_argument('--join_irc', action="store_true", help='Combine forward and backward IRCs. (string, default = false)')
    parser.add_argument('--input_gen', help="Use Built in PYREX input generator", default=None)
    parser.add_argument('--connectivity_matrix', action="store_true", help="Produce Connectivity Matrix for a single structure")
    parser.add_argument('--fsapt_analyze', action="store_true", help="Utility to perform reaction F-SAPT decomposition of reaction force")

def main():
    read_options()
    args = vars(parser.parse_args())
    # These are all methods that need to skip the normal pyrex stuff
    if(args["input_gen"]=="irc"):
        input_gen.build_irc()
    elif(args["input_gen"]=="energy"):
        input_gen.build_standard()
    elif(args["input_gen"]=="frag"):
        input_gen.build_frag()
    elif(args["join_irc"]==True):
        input_gen.join_irc()
    elif(args["rexplot"]==True):
        json_input = args["input"]
        rex.rexplot_utility(json_input)
    elif(args["fsapt_analyze"]==True):
        params = rex.read_params(args)
        params.fsapt_analyze = True
        rex.geometry_builder(params)
        rex.sapt_geometry_build(params)
        rex.fsapt_analysis(params)
    elif(args["connectivity_matrix"]==True):
        params = rex.read_params(args)
        rex.fragility_spectrum.single_connectivity_matrix(params)

    # Normal pyrex stuff
    else:
        params = rex.read_params(args)
        rex.calculate_irc(params)
        rex.calculate_surface_scan(params)
        rex.geometry_builder(params)
        rex.energy_calculations(params)
        rex.reaction_force_analysis(params)
        rex.reaction_force_constant(params)
        rex.activation_strain_analysis(params)
        rex.conceptual_dft_analysis(params)
        rex.total_reaction_work_integrals(params)
        rex.sapt_geometry_build(params)
        rex.atomic_force_decomposition(params)
        rex.sapt_energy_calc(params)
        rex.reaction_fragility_spectrum(params)
        #rexplot_utility(params)
        #flux_polarization(params)
        rex.success(params)

if __name__ == '__main__':
    main()
