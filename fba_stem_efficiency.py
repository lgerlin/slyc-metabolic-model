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

# photosynth red

addstoichMat_values.append(LEAF_STEM_PHOTOSYNTHESIS)
num_ratio_eff = len(addstoichMat_values)-1
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab)
addstoichMat_rows.append(mnet.nb_metab)

addstoichMat_columns.append(mnet.reactions_id.index('R_EX_photon_h_s'))
addstoichMat_columns.append(mnet.reactions_id.index('R_EX_photon_h_l'))

# input ratio

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

#                               STEM EFFICIENCY

#############################################################################

eff_reduction = [i for i in np.arange(1, 18, 0.2)]

photons_eff_list = []
eff_reduction_list = []
carbon_prod_stem_list = []
sucrose_uptake_stem_list = []
PSIHL = []
PSIHS = []
PSIIHS = []
PSIIHL = []
TKT1HS = []
TAHS = []
RPEHR = []
PYRDCS = []
ARODHS = []
TPIHS = []
XYLI2S = []
FBAS = []
PFKS = []
ENOS = []
FUMS = []
FUMMS = []
TPIS = []
ENOHS = []
HEX7S = []
RPEHS = []
RPEHL = []
PRUKL = []
PRUKS = []
PRUKR = []
TKT1HR = []
GLUKBRS = []
GLUKS = []
SUCSS = []
RPIHR = []
G6PBDHHR = []
PGIBHR = []
PFKHR = []
SUCR_ROOT = []
SUCR_LEAF = []
SUCR_PRODS = []

solution_status = []

for eff in eff_reduction:
    LEAF_STEM_PHOTOSYNTHESIS = eff

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

        addstoichMat_values[num_ratio_eff] = LEAF_STEM_PHOTOSYNTHESIS

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
        photons_eff_list.append(sum_photons_min)
        solution_status.append(
            probFBA.solution.status[probFBA.solution.get_status()])

        # PRINT RESULTS

        print ("Solution status = ",
               probFBA.solution.get_status(), ":",)
        print (probFBA.solution.status[probFBA.solution.get_status()])
        print ("Objective value  = ",
               probFBA.solution.get_objective_value())
        print ("Iteration count = ",
               probFBA.solution.progress.get_num_iterations())

        print(sum_photons_min)

        print('PS stem', x[
            mnet.reactions_id.index('R_EX_photon_h_s')])
        print('PS leaf', x[
            mnet.reactions_id.index('R_EX_photon_h_l')])

        # new constraints for flux minimization

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

        addstoichMat_values_2.append(mnet.STEM_WEIGHT)
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

            # PRINT RESULTS

            flux_sum = probFBA.solution.get_objective_value()

            print("flux sum", flux_sum)

            x2 = probFBA.solution.get_values()

            eff_reduction_list.append(eff)

            PSIHL.append(x2[mnet.reactions_id.index('R_PSIh_l')])
            PSIHS.append(x2[mnet.reactions_id.index('R_PSIh_s')])
            PSIIHS.append(x2[mnet.reactions_id.index('R_PSII_h_s')])
            PSIIHL.append(x2[mnet.reactions_id.index('R_PSII_h_l')])
            TKT1HS.append(x2[mnet.reactions_id.index('R_TKT1h_s')])
            TAHS.append(x2[mnet.reactions_id.index('R_TAh_s')])
            RPEHR.append(x2[mnet.reactions_id.index('R_RPEh_r')])
            PYRDCS.append(x2[mnet.reactions_id.index('R_PYRDC_s')])
            ARODHS.append(x2[mnet.reactions_id.index('R_ARODH_c_s')])
            TPIHS.append(x2[mnet.reactions_id.index('R_TPIh_s')])
            XYLI2S.append(x2[mnet.reactions_id.index('R_XYLI2_s')])
            FBAS.append(x2[mnet.reactions_id.index('R_FBA_s')])
            PFKS.append(x2[mnet.reactions_id.index('R_PFK_s')])
            ENOS.append(x2[mnet.reactions_id.index('R_ENO_s')])
            FUMS.append(x2[mnet.reactions_id.index('R_FUM_s')])
            FUMMS.append(x2[mnet.reactions_id.index('R_FUMm_s')])
            TPIS.append(x2[mnet.reactions_id.index('R_TPI_s')])
            ENOHS.append(x2[mnet.reactions_id.index('R_ENO_h_s')])
            HEX7S.append(x2[mnet.reactions_id.index('R_HEX7_s')])
            RPEHS.append(x2[mnet.reactions_id.index('R_RPEh_s')])
            RPEHL.append(x2[mnet.reactions_id.index('R_RPEh_l')])
            PRUKL.append(x2[mnet.reactions_id.index('R_PRUK_1_l')])
            PRUKS.append(x2[mnet.reactions_id.index('R_PRUK_1_s')])
            PRUKR.append(x2[mnet.reactions_id.index('R_PRUK_1_r')])
            TKT1HR.append(x2[mnet.reactions_id.index('R_TKT1h_r')])
            GLUKBRS.append(x2[mnet.reactions_id.index('R_GLUKBh_s')])
            GLUKS.append(x2[mnet.reactions_id.index('R_GLUK_s')])
            SUCSS.append(x2[mnet.reactions_id.index('R_SUCS_s')])
            RPIHR.append(x2[mnet.reactions_id.index('R_RPIh_r')])
            G6PBDHHR.append(x2[mnet.reactions_id.index('R_G6PBDHh_r')])
            PGIBHR.append(x2[mnet.reactions_id.index('R_PGIBh_r')])
            PFKHR.append(x2[mnet.reactions_id.index('R_PFKh_r')])
            SUCR_ROOT.append(
                x2[mnet.reactions_id.index('Exch_phl_r_M_sucr_c')])
            SUCR_PRODS.append(
                x2[mnet.reactions_id.index('Exch_s_phl_M_sucr_c')])
            carbon_prod_stem_list.append(
                x2[mnet.reactions_id.index('R_RBPCh_s')])
            sucrose_uptake_stem_list.append(
                x2[mnet.reactions_id.index('Exch_phl_s_M_sucr_c')])
            SUCR_LEAF.append(
                x2[mnet.reactions_id.index('Exch_l_phl_M_sucr_c')])
            print(eff)

        except CplexError as exc:
            print(exc)

    except CplexError as exc:
        print(exc)

# save results

df = pd.DataFrame({'Ratio efficiency leaf:stem': eff_reduction_list,
                   'Photon uptake': photons_eff_list,
                   'Carbon production stem': carbon_prod_stem_list,
                   'Sucrose uptake stem': sucrose_uptake_stem_list,
                   'R_PSIh_leaf': PSIHL,
                   'R_PSIh_stem': PSIHS,
                   'R_PSII_h_stem': PSIIHS,
                   'R_PSII_h_leaf': PSIIHL,
                   'R_TKT1h_stem': TKT1HS,
                   'R_TAh_stem': TAHS,
                   'R_RPEh_root': RPEHR,
                   'R_PYRDC_stem': PYRDCS,
                   'R_ARODH_c_stem': ARODHS,
                   'R_TPIh_stem': TPIHS,
                   'R_XYLI2_stem': XYLI2S,
                   'R_FBA_stem': FBAS,
                   'R_PFK_stem': PFKS,
                   'R_ENO_stem': ENOS,
                   'R_FUM_stem': FUMS,
                   'R_FUMm_stem': FUMMS,
                   'R_TPI_stem': TPIS,
                   'R_ENO_h_stem': ENOHS,
                   'R_HEX7_stem': HEX7S,
                   'R_RPEh_stem': RPEHS,
                   'R_RPEh_leaf': RPEHL,
                   'R_PRUK_1_leaf': PRUKL,
                   'R_PRUK_1_stem': PRUKS,
                   'R_PRUK_1_root': PRUKR,
                   'R_TKT1h_root': TKT1HR,
                   'R_GLUKBh_stem': GLUKBRS,
                   'R_GLUK_stem': GLUKS,
                   'R_SUCS_stem': SUCSS,
                   'R_RPIh_root': RPIHR,
                   'R_G6PBDHh_root': G6PBDHHR,
                   'R_PGIBh_root': PGIBHR,
                   'R_PFKh_root': PFKHR,
                   'Exchange_phloem_root_M_sucr_c': SUCR_ROOT,
                   'Exchange_leaf_phloem_M_sucr_c': SUCR_LEAF,
                  'Exchange_stem_phloem_M_sucr_c': SUCR_PRODS
                   })

df = df.reindex(columns=['Ratio efficiency leaf:stem',
                         'Photon uptake',
                         'Carbon production stem',
                         'Sucrose uptake stem',
                         'R_PSIh_leaf',
                         'R_PSIh_stem',
                         'R_PSII_h_stem',
                         'R_PSII_h_leaf',
                         'R_TKT1h_stem',
                         'R_TAh_stem',
                         'R_RPEh_root',
                         'R_PYRDC_stem',
                         'R_ARODH_c_stem',
                         'R_TPIh_stem',
                         'R_XYLI2_stem',
                         'R_FBA_stem',
                         'R_PFK_stem',
                         'R_ENO_stem',
                         'R_FUM_stem',
                         'R_FUMm_stem',
                         'R_TPI_stem',
                         'R_ENO_h_stem',
                         'R_HEX7_stem',
                         'R_RPEh_stem',
                         'R_RPEh_leaf',
                         'R_PRUK_1_leaf',
                         'R_PRUK_1_stem',
                         'R_PRUK_1_root',
                         'R_TKT1h_root',
                         'R_GLUKBh_stem',
                         'R_GLUK_stem',
                         'R_SUCS_stem',
                         'R_RPIh_root',
                         'R_G6PBDHh_root',
                         'R_PGIBh_root',
                         'R_PFKh_root',
                         'Exchange_phloem_root_M_sucr_c',
                         'Exchange_leaf_phloem_M_sucr_c',
                         'Exchange_stem_phloem_M_sucr_c'
                         ])

df.to_excel('output/surface_limit_stem.xlsx',
            sheet_name='sheet1', index=False)
