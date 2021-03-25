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

addstoichMat_values.append(0.17 / 0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab + 2)
addstoichMat_rows.append(mnet.nb_metab + 2)

addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_LEAF_l'))
addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_ROOT_r'))

addstoichMat_values.append(0.26 / 0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab + 3)
addstoichMat_rows.append(mnet.nb_metab + 3)

addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_LEAF_l'))
addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_STEM_s'))

#############################################################################

#                               STEM PROPORTION

#############################################################################

photons_perc_list = []
perc_stem_list = []
carbon_prod_stem_list = []
sucrose_uptake_stem_list = []

PSIHL = []
PSIHS = []
PSIIHS = []
PSIIHL = []

SUCR_ROOT = []
SUCR_LEAF = []
SUCR_PRODS = []

LEAF_WEIGHT = mnet.calib.value[
    mnet.calib_names.index('leaf_relative_weight')]
ROOT_WEIGHT = mnet.calib.value[
    mnet.calib_names.index('root_relative_weight')]

STEM_WEIGHT_initial = mnet.calib.value[
    mnet.calib_names.index('stem_relative_weight')]

PLANT_WEIGHT_initial = LEAF_WEIGHT + ROOT_WEIGHT + STEM_WEIGHT_initial

perc_stem = [i for i in np.arange(10, 50, 0.05 * 10)]
solution_status = []
for perc in perc_stem:

    obj = []

    obj = [0.0] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))

    STEM_WEIGHT = ((perc / 100) * (LEAF_WEIGHT + ROOT_WEIGHT)) / \
                  (1 - (perc / 100))

    PLANT_WEIGHT = LEAF_WEIGHT + ROOT_WEIGHT + STEM_WEIGHT

    for plant_stem in mnet.index_plant_stem:
        mnet.stoichMat_values[
            plant_stem
        ] = mnet.stoichMat_values[plant_stem] * (
                STEM_WEIGHT_initial / PLANT_WEIGHT_initial
        ) * PLANT_WEIGHT / STEM_WEIGHT

    for organ_stem in mnet.index_organ_stem:
        mnet.stoichMat_values[
            organ_stem
        ] = mnet.stoichMat_values[
                        organ_stem
            ] * PLANT_WEIGHT_initial / PLANT_WEIGHT

    STEM_WEIGHT_initial = STEM_WEIGHT
    PLANT_WEIGHT_initial = PLANT_WEIGHT

    obj[
        mnet.reactions_id.index(objective.reaction[0])
    ] = float(1.0 * LEAF_WEIGHT)

    obj[
        mnet.reactions_id.index(objective.reaction[1])
    ] = float(1.0 * STEM_WEIGHT)

    #########################################################################

    #                      FLUX BALANCE ANALYSIS
    #                     1) Photon uptake minimization
    #                     2) ABS(Flux) sum minimization

    #########################################################################

    try:
        #####################################
        # FBA 1) Photon uptake minimization
        #####################################
        x = []

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
        x = probFBA.solution.get_values()
        sum_photons_min = probFBA.solution.get_objective_value()

        # PRINT RESULTS

        print ("Solution status = ",
               probFBA.solution.get_status(), ":",)
        print (probFBA.solution.status[probFBA.solution.get_status()])
        print ("Objective value  = ",
               probFBA.solution.get_objective_value())
        print ("Iteration count = ",
               probFBA.solution.progress.get_num_iterations())
        print("Sump photons min", sum_photons_min)
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

        print('PS stem', x[
            mnet.reactions_id.index('R_EX_photon_h_s')])
        print('PS leaf', x[
            mnet.reactions_id.index('R_EX_photon_h_l')])

        # add constraints for flux min

        photons_perc_list.append(sum_photons_min)
        solution_status.append(
            probFBA.solution.status[probFBA.solution.get_status()])

        photons_constraint_leaf = x[
            mnet.reactions_id.index('R_EX_photon_h_l')]
        photons_constraint_stem = x[
            mnet.reactions_id.index('R_EX_photon_h_s')]

        no3_constraint_root = x[
            mnet.reactions_id.index('R_EX_no3_c_r')]
        nh4_constraint_root = x[
            mnet.reactions_id.index('R_EX_nh4_c_r')]

        obj2 = []

        obj2 = [0] * (len(mnet.reactions_id) + mnet.nb_exchangereac) + \
               [1] * (len(mnet.reactions_id))
        lb2 = []

        lb2 = lb + [-cplex.infinity] * (len(mnet.reactions_id))

        ub2 = []

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

        rhs2 = []

        rhs2 = [0] * (mnet.nb_metab + (len(mnet.reactions_id) * 2)) +\
               [0] + [0] + [sum_photons_min * 1.01]

        sense2 = []

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

        addstoichMat_values_2 = []
        addstoichMat_rows_2 = []
        addstoichMat_columns_2 = []

        addstoichMat_values_2.append(STEM_WEIGHT)
        addstoichMat_values_2.append(mnet.LEAF_WEIGHT)

        addstoichMat_rows_2.append(minstoichMat_rows[-1] + 1)
        addstoichMat_rows_2.append(minstoichMat_rows[-1] + 1)

        addstoichMat_columns_2.append(
            mnet.reactions_id.index('R_EX_photon_h_s'))
        addstoichMat_columns_2.append(
            mnet.reactions_id.index('R_EX_photon_h_l'))

        try:
            #####################################
            # FBA 2) ABS(Flux) sum minimization
            #####################################
            x2 = []

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

            probFBA.linear_constraints.set_coefficients(
                zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows +
                    minstoichMat_rows + addstoichMat_rows +
                    addstoichMat_rows_2,
                    mnet.stoichMat_columns +
                    mnet.exchangestoichMat_columns +
                    minstoichMat_columns +
                    addstoichMat_columns +
                    addstoichMat_columns_2,
                    mnet.stoichMat_values +
                    mnet.exchangestoichMat_values +
                    minstoichMat_values + addstoichMat_values +
                    addstoichMat_values_2))

            probFBA.solve()
            x2 = probFBA.solution.get_values()
            flux_sum = probFBA.solution.get_objective_value()

            # PRINT RESULTS

            print("flux sum", flux_sum)

            perc_stem_list.append(perc)

            carbon_prod_stem = x2[mnet.reactions_id.index('R_RBPCh_s')]
            sucrose_uptake_stem = x2[
                mnet.reactions_id.index('Exch_phl_s_M_sucr_c')]
            carbon_prod_stem_list.append(carbon_prod_stem)
            sucrose_uptake_stem_list.append(sucrose_uptake_stem)
            PSIHL.append(x2[mnet.reactions_id.index('R_PSIh_l')])
            PSIHS.append(x2[mnet.reactions_id.index('R_PSIh_s')])
            PSIIHS.append(x2[mnet.reactions_id.index('R_PSII_h_s')])
            PSIIHL.append(x2[mnet.reactions_id.index('R_PSII_h_l')])
            SUCR_ROOT.append(
                x2[mnet.reactions_id.index('Exch_phl_r_M_sucr_c')])
            SUCR_PRODS.append(
                x2[mnet.reactions_id.index('Exch_s_phl_M_sucr_c')])
            SUCR_LEAF.append(
                x2[mnet.reactions_id.index('Exch_l_phl_M_sucr_c')])

        except CplexError as exc:
            print(exc)

    except CplexError as exc:
        print(exc)

# save results

df = pd.DataFrame({'Perc. STEM': perc_stem_list,
                   'Photon uptake': photons_perc_list,
                   'Carbon production stem': carbon_prod_stem_list,
                   'Sucrose uptake stem': sucrose_uptake_stem_list,
                   'R_PSIh_leaf': PSIHL,
                   'R_PSIh_stem': PSIHS,
                   'R_PSII_h_stem': PSIIHS,
                   'R_PSII_h_leaf': PSIIHL,
                   'Exchange_phloem_root_M_sucr_c': SUCR_ROOT,
                   'Exchange_leaf_phloem_M_sucr_c': SUCR_LEAF,
                   'Exchange_stem_phloem_M_sucr_c': SUCR_PRODS,
                   'Solution status': solution_status
                   })

df = df.reindex(columns=['Perc. STEM',
                         'Photon uptake',
                         'Carbon production stem',
                         'Sucrose uptake stem',
                         'R_PSIh_leaf',
                         'R_PSIh_stem',
                         'R_PSII_h_stem',
                         'R_PSII_h_leaf',
                         'Exchange_phloem_root_M_sucr_c',
                         'Exchange_leaf_phloem_M_sucr_c',
                         'Exchange_stem_phloem_M_sucr_c',
                         'Solution status'
                         ])

df.to_excel('output/perc_stem.xlsx',
            sheet_name='sheet1', index=False)
