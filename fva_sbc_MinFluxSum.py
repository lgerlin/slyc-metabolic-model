# -*- coding:Utf8 -*-

#############################################################################
# Program Python type
# authors: Gerlin et al., 2021
#############################################################################

#############################################################################
# External functions
import cplex as cplex
from cplex.exceptions import CplexError
import pandas as pd
import numpy as np
import network_parser as mnet  # open the parsed network

#############################################################################

#                                MAIN PROGRAM

#############################################################################


#############################################################################

#                                GENERAL VARIABLES

#############################################################################

total_r = len(mnet.reactions_id) + mnet.nb_exchangereac

#############################################################################

#                DEFINE THE CONSTRAINTS ON THE 3 COMP MODEL

#############################################################################
# Upper bounds
ub = [10000] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id)) + [cplex.infinity] * len(mnet.reactions_id)


# Lower bounds. Set to 0 for exchange reactions to avoid infinite uptake of substrates other than those wished
lb = [-10000 if rev == 1 else 0 for rev in mnet.reactions_reversible] \
     + [0] * mnet.nb_exchangereac * 1 + [-cplex.infinity] * len(mnet.reactions_id)

# Right-hand side of the LP problem
rhs = [0] * mnet.nb_metab + [0] * len(mnet.reactions_id) * 2

# Senses of the LP problem
sense = "E" * mnet.nb_metab + "L" * len(mnet.reactions_id) * 2

# Objective function of the LP problem (minimizations sum fluxes)
obj = [0.0] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id)) + [1] * len(mnet.reactions_id)

# Variables' initialisation for additional constraints of the LP problem, not on lb or ub
minnames = ["MinFluxConstraint1_" + id for id in mnet.reactions_id] + \
           ["MinFluxConstraint2_" + id for id in mnet.reactions_id]

minstoichMat_rows = [index for index in
                     range(mnet.nb_metab, mnet.nb_metab + len(mnet.reactions_id))] + \
                    [index for index in
                     range(mnet.nb_metab, mnet.nb_metab + len(mnet.reactions_id))] + \
                    [index for index in range(mnet.nb_metab + len(mnet.reactions_id),
                                              mnet.nb_metab + len(mnet.reactions_id) * 2)] + \
                    [index for index in range(mnet.nb_metab + len(mnet.reactions_id),
                                              mnet.nb_metab + len(mnet.reactions_id) * 2)]
minstoichMat_columns = [index for index in range(len(mnet.reactions_id))] + \
                       [index for index in range(len(mnet.exchangereactions_id) + len(mnet.reactions_id),
                                                 len(mnet.exchangereactions_id) + len(mnet.reactions_id) * 2)] + \
                       [index for index in range(len(mnet.reactions_id))] + \
                       [index for index in range(len(mnet.exchangereactions_id) + len(mnet.reactions_id),
                                                 len(mnet.exchangereactions_id) + len(mnet.reactions_id) * 2)]

minstoichMat_values = [1] * len(mnet.reactions_id) + [-1] * len(mnet.reactions_id) * 3


#############################################################################

#                      FLUX BALANCE ANALYSIS
#                     1) Photon uptake minimization
#                     2) ABS(Flux) sum minimization

#############################################################################

try:
    #############################
    # FBA 1) Flux sum minimized
    #############################

    probFBA = cplex.Cplex()

    probFBA.set_problem_name("FBATomatoMinPhotons")

    probFBA.objective.set_sense(probFBA.objective.sense.minimize)

    probFBA.linear_constraints.add(rhs=rhs, senses=sense,
                                   names=["QSSA_" + id for id in mnet.metabolites_id] + minnames)

    probFBA.variables.add(obj=obj, lb=lb, ub=ub, names=mnet.reactions_id + mnet.exchangereactions_id + ["MinVar_" + id for id in mnet.reactions_id])

    probFBA.linear_constraints.set_coefficients(
        zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows + minstoichMat_rows,
            mnet.stoichMat_columns + mnet.exchangestoichMat_columns + minstoichMat_columns,
            mnet.stoichMat_values + mnet.exchangestoichMat_values + minstoichMat_values))

    probFBA.write("output/LP/" + "FBATomatoMinFluxSum" + "_probFBA_SBC.lp")

    probFBA.solve()

except CplexError as exc:
    print(exc)

# Print results
print("Solution status = ", probFBA.solution.get_status(), ":", )
print(probFBA.solution.status[probFBA.solution.get_status()])
print("Objective value  = ", probFBA.solution.get_objective_value())
print("Iteration count = ", probFBA.solution.progress.get_num_iterations())


#############################################################################

#                           FLUX VARIABILITY ANALYSIS
#                            On exchange fluxes only

#############################################################################

print()
print("Starting FVA")

#define new optimization problem by adding constraint on sum of min fluxes
rhsFVA = rhs + [probFBA.solution.get_objective_value()]
senseFVA = sense + 'E'

lbFVA = lb
ubFVA = ub

constraint_minSumFlux_rows = [len(senseFVA)-1]*len(mnet.reactions_id)
constraint_minSumFlux_columns = list(range(len(mnet.reactions_id) + len(mnet.exchangereactions_id), len(lb)))
constraint_minSumFlux_values = [1] * len(mnet.reactions_id)


#create variables to store FVA results
statusMin = ['NA'] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))
statusMax = ['NA'] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))
FluxesMin = ['NA'] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))
FluxesMax = ['NA'] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))


#perform FVA on list of reactions of FVAreactions.csv file
FVAreactions = pd.read_csv("input/FVAreactions.csv", sep=";").id

for FVAreaction in FVAreactions:
    index_reaction = mnet.reactions_id.index(FVAreaction)

    #objective function: min or max of each reaction flux
    objFVA = [0]*(len(mnet.reactions_id) + len(mnet.exchangereactions_id)) + [0] * len(mnet.reactions_id)
    objFVA[index_reaction] = 1

    try:

        ##################################################
        # FBA 2) et FBA 3) min et max flux of each fluxes
        #################################################

        probFVA = cplex.Cplex()

        probFVA.set_problem_name("FVATomatoMinFluxSBC")

        probFVA.linear_constraints.add(rhs=rhsFVA, senses=senseFVA,
                                        names=["QSSA_" + id for id in mnet.metabolites_id] + minnames + ["Sum flux constraint"])

        probFVA.variables.add(obj=objFVA, lb=lbFVA, ub=ubFVA,
                               names=mnet.reactions_id + mnet.exchangereactions_id + ["MinVar_" + id for id in mnet.reactions_id])
        probFVA.linear_constraints.set_coefficients(
            zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows + minstoichMat_rows + constraint_minSumFlux_rows,
                mnet.stoichMat_columns + mnet.exchangestoichMat_columns + minstoichMat_columns + constraint_minSumFlux_columns,
                mnet.stoichMat_values + mnet.exchangestoichMat_values + minstoichMat_values + constraint_minSumFlux_values))

        #no terminal output to save time
        probFVA.set_log_stream(None)
        probFVA.set_error_stream(None)
        probFVA.set_warning_stream(None)
        probFVA.set_results_stream(None)

        # solve FVA min
        probFVA.objective.set_sense(probFVA.objective.sense.minimize)
        probFVA.solve()
        probFVA.write("output/LP/" + "FVATomatoMinFlux" + "_probFVA.lp")

        # store solution
        statusMin[index_reaction]= probFVA.solution.get_status()
        FluxesMin[index_reaction]= probFVA.solution.get_objective_value()

        # Print results
        #print("Solution status = ", probFVA.solution.get_status(), ":", )
        #print(probFVA.solution.status[probFVA.solution.get_status()])
        #print("Objective value  = ", probFVA.solution.get_objective_value())


        # solve FVA max
        probFVA.objective.set_sense(probFVA.objective.sense.maximize)
        probFVA.solve()
        probFVA.write("output/LP/" + "FVATomatoMaxFlux" + "_probFVA.lp")

        # store solution
        statusMax[index_reaction] = probFVA.solution.get_status()
        FluxesMax[index_reaction] = probFVA.solution.get_objective_value()

        # Print results
        #print("Solution status = ", probFVA.solution.get_status(), ":", )
        #print(probFVA.solution.status[probFVA.solution.get_status()])
        #print("Objective value  = ", probFVA.solution.get_objective_value())

    except CplexError as exc:
        print(exc)


print("FVA finished")
print()

# save results
dfFVA = pd.DataFrame({
    'Reaction Name': mnet.reactions_name + mnet.exchangereactions_name,
    'Reaction Formula Name': mnet.reactions_formulasNames + mnet.exchangereactions_formulasNames,
    'Reaction Number': range(0, len(mnet.reactions_id + mnet.exchangereactions_id)),
    'Reaction Id': mnet.reactions_id + mnet.exchangereactions_id,
    'Reaction Formula Id': mnet.reactions_formulas + mnet.exchangereactions_formulas,
    'FVA min': FluxesMin,
    'FVA max': FluxesMax,
    'FVA status min': statusMin,
    'FVA status max': statusMax})

val_atp = str(mnet.cost)
dfFVA.to_excel('output/FVAMinFlux_SBC.xlsx',
               sheet_name='FVAMinFlux_SBC', index=False)