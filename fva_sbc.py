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
ub = [10000] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))


# Lower bounds. Set to 0 for exchange reactions to avoid infinite uptake of substrates other than those wished
lb = [-10000 if rev == 1 else 0 for rev in mnet.reactions_reversible] \
     + [0] * mnet.nb_exchangereac * 1

# Right-hand side of the LP problem
rhs = [0] * mnet.nb_metab

# Senses of the LP problem
sense = "E" * mnet.nb_metab

# Objective function of the LP problem
obj = [0.0] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))

#############################################################################

#                           FLUX VARIABILITY ANALYSIS
#                            On exchange fluxes only

#############################################################################

print()
print("Starting FVA")

#define optimization problem
rhsFVA = rhs
senseFVA = sense

lbFVA = lb
ubFVA = ub

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
    objFVA = [0]*(len(mnet.reactions_id) + len(mnet.exchangereactions_id))
    objFVA[index_reaction] = 1

    try:

        ##################################################
        # FBA 1) et FBA 2) min et max flux of each fluxes
        #################################################

        probFVA = cplex.Cplex()

        probFVA.set_problem_name("FVATomatoSBC")

        probFVA.linear_constraints.add(rhs=rhsFVA, senses=senseFVA,
                                        names=["QSSA_" + id for id in mnet.metabolites_id])

        probFVA.variables.add(obj=objFVA, lb=lbFVA, ub=ubFVA,
                               names=mnet.reactions_id + mnet.exchangereactions_id)
        probFVA.linear_constraints.set_coefficients(
            zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows,
                mnet.stoichMat_columns + mnet.exchangestoichMat_columns,
                mnet.stoichMat_values + mnet.exchangestoichMat_values))

        #no terminal output to save time
        probFVA.set_log_stream(None)
        probFVA.set_error_stream(None)
        probFVA.set_warning_stream(None)
        probFVA.set_results_stream(None)

        # solve FVA min
        probFVA.objective.set_sense(probFVA.objective.sense.minimize)
        probFVA.solve()
        #probFVA.write("output/LP/" + "FVATomato" + "_probFVA_SBC.lp")

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
        #probFVA.write("output/LP/" + "FVATomato" + "_probFVA_SBC.lp")

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
dfFVA.to_excel('output/FVA_SBC.xlsx',
               sheet_name='FVA_SBC', index=False)