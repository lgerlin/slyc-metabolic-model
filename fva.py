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
import network_parser as mnet  # import parsed network

#############################################################################

#                                MAIN PROGRAM

#############################################################################


#############################################################################

#                DEFINE THE CONSTRAINTS ON THE 3 COMP MODEL

#############################################################################

NH4_NO3 = mnet.NH4_NO3
LEAF_STEM_PHOTOSYNTHESIS = mnet.LEAF_STEM_PHOTOSYNTHESIS

cost = mnet.cost
xyl = mnet.xyl

# Upper bounds
ub = [cplex.infinity] * (
        len(mnet.reactions_id) + len(mnet.exchangereactions_id))

# Lower bounds
lb = [-cplex.infinity if rev == 1 else 0 for rev in mnet.reactions_reversible
      ] + [-cplex.infinity] * mnet.nb_exchangereac * 1

# Right-hand side
rhs = [0] * mnet.nb_metab + [0] + [0] + [0] + [0]

# Senses of the problem
sense = "E" * mnet.nb_metab + "L" + "E" + "E" + "E"

# Constraints
constraints = pd.read_csv("input/constraints.csv", sep=";")
for num, id in enumerate(constraints.id):
    lb[mnet.reactions_id.index(id)] = float(constraints.lb[num])
    ub[mnet.reactions_id.index(id)] = float(constraints.ub[num])

# Exchanges
exchanges = pd.read_csv("input/exchanges.csv", sep=";")
for organ in mnet.organs_id:
    for num, id in enumerate(exchanges.id):
        lb[
            len(
                mnet.reactions_id
            ) + mnet.exchangereactions_id.index(id + '_' + organ)
            ] = float(exchanges.lb[num])
        ub[
            len(
                mnet.reactions_id
                ) + mnet.exchangereactions_id.index(id + '_' + organ)
            ] = float(exchanges.ub[num])

# xylem fluxes

xylem = pd.read_csv("input/xylem.csv", sep=";")

if xyl == 1:
    for num, id in enumerate(xylem.id):
        lb[mnet.reactions_id.index(id)] = xylem.lb[num]
        ub[mnet.reactions_id.index(id)] = xylem.ub[num]

# Objective function to minimize photons import

obj = [0.0] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))

objective = pd.read_csv("input/objective.csv", sep=";")
num_obj = len(
    mnet.reactions_id
) + mnet.exchangereactions_id.index(objective.exchange[0])
ub[num_obj] = cplex.infinity
lb[num_obj] = -cplex.infinity

num_obj = len(
    mnet.reactions_id
) + mnet.exchangereactions_id.index(objective.exchange[1])
ub[num_obj] = cplex.infinity
lb[num_obj] = -cplex.infinity

obj[
    mnet.reactions_id.index(objective.reaction[0])
] = float(1.0 * mnet.LEAF_WEIGHT)
obj[
    mnet.reactions_id.index(objective.reaction[1])
] = float(1.0 * mnet.STEM_WEIGHT)


# MATRIX FOR RATIO CONSTRAINTS

addstoichMat_rows = []
addstoichMat_values = []
addstoichMat_columns = []

total_r = len(mnet.reactions_id) + mnet.nb_exchangereac

# photosynth red constraint

addstoichMat_values.append(LEAF_STEM_PHOTOSYNTHESIS)
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab)
addstoichMat_rows.append(mnet.nb_metab)

addstoichMat_columns.append(mnet.reactions_id.index('R_EX_photon_h_s'))
addstoichMat_columns.append(mnet.reactions_id.index('R_EX_photon_h_l'))

# nitrogen ratio constraint

addstoichMat_values.append(NH4_NO3)
addstoichMat_values.append(-1)
# ROWS
addstoichMat_rows.append(mnet.nb_metab + 1)
addstoichMat_rows.append(mnet.nb_metab + 1)
# COLUMNS
addstoichMat_columns.append(mnet.reactions_id.index('R_EX_no3_c_r'))
addstoichMat_columns.append(mnet.reactions_id.index('R_EX_nh4_c_r'))

# add constraint for growth

addstoichMat_values.append(0.17/0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab + 2)
addstoichMat_rows.append(mnet.nb_metab + 2)

addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_LEAF_l'))
addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_ROOT_r'))

addstoichMat_values.append(0.26/0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab + 3)
addstoichMat_rows.append(mnet.nb_metab + 3)

addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_LEAF_l'))
addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_STEM_s'))

objsolmax = []
objsolmin = []

#############################################################################

#                      FLUX BALANCE ANALYSIS
#                     1) Photon uptake minimization
#                     2) ABS(Flux) sum minimization

#############################################################################

try:
    #####################################
    # FBA 1) Photon uptake minimization
    #####################################

    probFBA = cplex.Cplex()

    probFBA.set_problem_name("FBATomatoMinPhotons")

    probFBA.objective.set_sense(probFBA.objective.sense.minimize)

    probFBA.linear_constraints.add(
        rhs=rhs, senses=sense,
        names=mnet.metabolites_id + ['other_constraints']*4)

    probFBA.variables.add(
        obj=obj, lb=lb, ub=ub,
        names=mnet.reactions_id + mnet.exchangereactions_id)

    probFBA.linear_constraints.set_coefficients(
        zip(mnet.stoichMat_rows +
            mnet.exchangestoichMat_rows + addstoichMat_rows,
            mnet.stoichMat_columns +
            mnet.exchangestoichMat_columns +
            addstoichMat_columns,
            mnet.stoichMat_values +
            mnet.exchangestoichMat_values +
            addstoichMat_values))

    probFBA.solve()

    # PRINT RESULTS

    print ("Solution status = ",
           probFBA.solution.get_status(), ":",)
    print (probFBA.solution.status[probFBA.solution.get_status()])
    print ("Objective value  = ",
           probFBA.solution.get_objective_value())
    print ("Iteration count = ",
           probFBA.solution.progress.get_num_iterations())

    x = probFBA.solution.get_values()

    sum_photons_min = probFBA.solution.get_objective_value()

    photons_constraint_leaf = x[
        mnet.reactions_id.index('R_EX_photon_h_l')]
    photons_constraint_stem = x[
        mnet.reactions_id.index('R_EX_photon_h_s')]

    no3_constraint_root = x[
        mnet.reactions_id.index('R_EX_no3_c_r')]
    nh4_constraint_root = x[
        mnet.reactions_id.index('R_EX_nh4_c_r')]

    print("NH4", x[
        mnet.reactions_id.index('R_EX_nh4_c_r')])
    print("NO3", x[
        mnet.reactions_id.index('R_EX_no3_c_r')])

    print("co2_l", x[
        mnet.reactions_id.index('R_EX_co2_e_l')])
    print("co2_s", x[
        mnet.reactions_id.index('R_EX_co2_e_s')])
    print("co2_r", x[
        mnet.reactions_id.index('R_EX_co2_e_r')])

    obj2 = [0] * (len(mnet.reactions_id) + mnet.nb_exchangereac) + \
           [1] * (len(mnet.reactions_id))

    lb2 = lb + [-cplex.infinity] * (len(mnet.reactions_id))

    ub2 = ub + [cplex.infinity] * (len(mnet.reactions_id))

    ub2[
        mnet.reactions_id.index('R_EX_photon_h_l')
    ] = sum_photons_min

    ub2[
        mnet.reactions_id.index('R_EX_photon_h_s')
    ] = sum_photons_min

    ub2[
        mnet.reactions_id.index('R_EX_no3_c_r')
    ] = no3_constraint_root

    ub2[
        mnet.reactions_id.index('R_EX_nh4_c_r')
    ] = nh4_constraint_root

    lb2[
        mnet.reactions_id.index('R_BIOMASS_LEAF_l')
    ] = x[
        mnet.reactions_id.index('R_BIOMASS_LEAF_l')
    ]

    lb2[
        mnet.reactions_id.index('R_BIOMASS_STEM_s')
    ] = x[
        mnet.reactions_id.index('R_BIOMASS_STEM_s')
    ]

    lb2[
        mnet.reactions_id.index('R_BIOMASS_ROOT_r')
    ] = x[
        mnet.reactions_id.index('R_BIOMASS_ROOT_r')
    ]

    rhs2 = [0] * (mnet.nb_metab + (len(mnet.reactions_id) * 2)) +\
           [0] + [0] + [sum_photons_min * 1.01]

    sense2 = "E" * mnet.nb_metab + \
             "L" * (len(mnet.reactions_id) * 2) + "L" + \
             "E" + "L"

    minstoichMat_rows = []
    minstoichMat_values = []
    minstoichMat_columns = []

    minstoichMat_values = [1] * len(mnet.reactions_id) + \
                          [-1] * len(mnet.reactions_id) * 3

    # ROWS

    for i in range(mnet.nb_metab - 1, mnet.nb_metab +
                   (len(mnet.reactions_id)) * 2 - 1):
        minstoichMat_rows.append(i + 3)

    for i in range(mnet.nb_metab - 1, mnet.nb_metab +
                   (len(mnet.reactions_id)) * 2 - 1):
        minstoichMat_rows.append(i + 3)

    # COLUMNS

    for i in range(-1, len(mnet.reactions_id) - 1):
        minstoichMat_columns.append(i + 1)

    for i in range(-1, len(mnet.reactions_id) - 1):
        minstoichMat_columns.append(i + 1)

    for i in range(
            total_r - 1, total_r + len(mnet.reactions_id) - 1):
        minstoichMat_columns.append(i + 1)

    for i in range(
            total_r - 1, total_r + len(mnet.reactions_id) - 1):
        minstoichMat_columns.append(i + 1)

    addstoichMat_values.append(mnet.STEM_WEIGHT)
    addstoichMat_values.append(mnet.LEAF_WEIGHT)

    addstoichMat_rows.append(minstoichMat_rows[-1] + 1)
    addstoichMat_rows.append(minstoichMat_rows[-1] + 1)

    addstoichMat_columns.append(
        mnet.reactions_id.index('R_EX_photon_h_s'))
    addstoichMat_columns.append(
        mnet.reactions_id.index('R_EX_photon_h_l'))

    print('PS stem', x[
        mnet.reactions_id.index('R_EX_photon_h_s')])
    print('PS leaf', x[
        mnet.reactions_id.index('R_EX_photon_h_l')])

    try:
        #####################################
        # FBA 2) ABS(Flux) sum minimization
        #####################################

        probFBA = cplex.Cplex()

        probFBA.set_problem_name("FBATomatoMinFluxSum")

        probFBA.objective.set_sense(probFBA.objective.sense.minimize)

        probFBA.linear_constraints.add(rhs=rhs2, senses=sense2,
                                       names=mnet.metabolites_id +
                                       mnet.reactions_id +
                                       mnet.reactions_id +
                                       ['additional constraint']*3)

        probFBA.variables.add(obj=obj2, lb=lb2, ub=ub2,
                              names=mnet.reactions_id +
                              mnet.exchangereactions_id +
                              mnet.reactions_id)
        print(len(rhs2))
        print(len(mnet.stoichMat_rows + mnet.exchangestoichMat_rows +
                  minstoichMat_rows + addstoichMat_rows))

        probFBA.linear_constraints.set_coefficients(
            zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows +
                minstoichMat_rows + addstoichMat_rows,
                mnet.stoichMat_columns +
                mnet.exchangestoichMat_columns +
                minstoichMat_columns +
                addstoichMat_columns,
                mnet.stoichMat_values +
                mnet.exchangestoichMat_values +
                minstoichMat_values + addstoichMat_values))

        probFBA.solve()

        # PRINT RESULTS

        print (
            probFBA.solution.status[
                probFBA.solution.get_status()])
        flux_sum = probFBA.solution.get_objective_value()
        print("flux sum", flux_sum)
        x2 = probFBA.solution.get_values()

    except CplexError as exc:
        print(exc)

except CplexError as exc:
    print(exc)

#############################################################################

#                           FLUX VARIABILITY ANALYSIS
#                            On exchange fluxes only

#############################################################################

variable = pd.read_csv("input/variable.csv", sep=";")

list_variable = variable.id

rhs3 = rhs2 + [flux_sum * 1.001]

sense3 = sense2 + "L"

fba_values_variable = []

for i in range(total_r - 1, total_r + len(mnet.reactions_id) - 1):
    addstoichMat_columns.append(i + 1)
    addstoichMat_values.append(1)
    addstoichMat_rows.append(minstoichMat_rows[-1] + 2)

for reac in list_variable:

    fba_values_variable.append(x2[mnet.reactions_id.index(reac)])

    obj3 = [0] * (len(mnet.reactions_id) + mnet.nb_exchangereac) + \
           [0] * (len(mnet.reactions_id))

    obj3[mnet.reactions_id.index(reac)] = 1

    try:

        # SOLVE MINIMIZATION OF FLUX PROBLEM

        probFBA = cplex.Cplex()

        probFBA.set_problem_name("FBATomatoMinFlux")

        probFBA.objective.set_sense(probFBA.objective.sense.minimize)

        info = ["1:1 ratio for NO3 and NH4 uptake"]

        probFBA.linear_constraints.add(rhs=rhs3, senses=sense3,
                                       names=mnet.metabolites_id +
                                       mnet.reactions_id +
                                       mnet. reactions_id +
                                       ['other constraints']*4)

        probFBA.variables.add(obj=obj3, lb=lb2, ub=ub2,
                              names=mnet.reactions_id +
                              mnet.exchangereactions_id + mnet.reactions_id)

        probFBA.linear_constraints.set_coefficients(
            zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows +
                minstoichMat_rows +
                addstoichMat_rows,
                mnet.stoichMat_columns +
                mnet.exchangestoichMat_columns + minstoichMat_columns +
                addstoichMat_columns,
                mnet.stoichMat_values +
                mnet.exchangestoichMat_values + minstoichMat_values +
                addstoichMat_values))

        probFBA.solve()
        numcols = probFBA.variables.get_num()
        objsolmin.append(probFBA.solution.get_objective_value())
        x = probFBA.solution.get_values()

        # PRINT RESULTS
        print(reac)
        print("Solution status = ", probFBA.solution.get_status(), ":",)
        print(probFBA.solution.status[probFBA.solution.get_status()])
        print("Objective value  = ", probFBA.solution.get_objective_value())
        print("Iteration count = ",
              probFBA.solution.progress.get_num_iterations())
        print("Solution =", x)

    except CplexError as exc:
        print(exc)
        objsolmin.append("error")

    try:

        # SOLVE MAXIMIZATION OF FLUX PROBLEM

        probFBA = cplex.Cplex()

        probFBA.set_problem_name("FBATomatoMaxFlux")

        probFBA.objective.set_sense(probFBA.objective.sense.maximize)

        probFBA.linear_constraints.add(rhs=rhs3, senses=sense3,
                                       names=mnet.metabolites_id +
                                       mnet.reactions_id +
                                       mnet.reactions_id +
                                       ['other constraints']*4)

        probFBA.variables.add(obj=obj3, lb=lb2, ub=ub2,
                              names=mnet.reactions_id +
                              mnet.exchangereactions_id +
                              mnet.reactions_id)

        probFBA.linear_constraints.set_coefficients(
            zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows +
                minstoichMat_rows + addstoichMat_rows,
                mnet.stoichMat_columns + mnet.exchangestoichMat_columns +
                minstoichMat_columns + addstoichMat_columns,
                mnet.stoichMat_values + mnet.exchangestoichMat_values +
                minstoichMat_values + addstoichMat_values))

        probFBA.solve()

        # PRINT RESULTS

        print("Solution status = ", probFBA.solution.get_status(), ":",)
        print(probFBA.solution.status[probFBA.solution.get_status()])
        print("Objective value  = ", probFBA.solution.get_objective_value())
        print("Iteration count = ",
              probFBA.solution.progress.get_num_iterations())

        x = probFBA.solution.get_values()

        numcols = probFBA.variables.get_num()
        print("Solution =", x)
        objsolmax.append(probFBA.solution.get_objective_value())

    except CplexError as exc:
        print(exc)
        objsolmax.append("error")

val_atp = str(cost)

# save results

df = pd.DataFrame({'Reaction Id': list_variable, 'min': objsolmin,
                   'max': objsolmax, 'fba': fba_values_variable})

df = df.reindex(
    columns=['Reaction Id', 'min', 'max', 'fba'])

df.to_excel('output/FVA'
            '_cost_{0}_ratio_{1}_Xyl_{2}.xlsx'.format(
                val_atp, LEAF_STEM_PHOTOSYNTHESIS, xyl),
            sheet_name='sheet1', index=False)
