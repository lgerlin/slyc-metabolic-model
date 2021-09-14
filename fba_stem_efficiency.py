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
ub = [cplex.infinity] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))


# Lower bounds. Set to 0 for exchange reactions to avoid infinite uptake of substrates other than those wished
lb = [-cplex.infinity if rev == 1 else 0 for rev in mnet.reactions_reversible] \
     + [0] * mnet.nb_exchangereac * 1


# Additional constraints such as ATP maintenance ou leaf biomass growth
#leaf biomass growth
lb[mnet.reactions_id.index('R_BIOMASS_LEAF_l')] = mnet.leaf_biomass_growth

#ATP maintenance
lb[mnet.reactions_id.index('R_ATPS_l')] = mnet.atpm_l
lb[mnet.reactions_id.index('R_ATPS_s')] = mnet.atpm_s
lb[mnet.reactions_id.index('R_ATPS_r')] = mnet.atpm_r


# Change of bounds for boundary metabolites allowed to be exchanged
infinite_exchange = ['M_ca2_b_r',
                     'M_cl_b_r',
                     'M_so4_b_r',
                     'M_pi_b_r',
                     'M_fe2_b_r',
                     'M_k_b_r',
                     'M_mg2_b_r',
                     'M_na_b_r',
                     'M_nh4_b_r',
                     'M_no3_b_r',
                     'M_h2o_b_r',
                     'M_co2_b_r', 'M_co2_b_s', 'M_co2_b_l',
                     'M_o2_b_r', 'M_o2_b_s', 'M_o2_b_l',
                     'M_h_b_r', 'M_h_b_s', 'M_h_b_l',
                     'M_photon_b_s', 'M_photon_b_l',
                     'M_nadh_work_b_r', 'M_nadh_work_b_s', 'M_nadh_work_b_l',
                     'M_nadph_work_b_r', 'M_nadph_work_b_s', 'M_nadph_work_b_l',
                     'M_atp_work_b_r', 'M_atp_work_b_s', 'M_atp_work_b_l',
                     'M_biomass_leaf_b_l',
                     'M_biomass_root_b_r',
                     'M_biomass_stem_b_s']

for id in infinite_exchange:
    lb[mnet.nb_reac + mnet.exchangereactions_id.index('Exch_'+id)] = -cplex.infinity


# [OPTIONAL] Constraint on xylem composition (taken from xylem.csv file)
if mnet.xyl == 1:
    xylem = pd.read_csv("input/xylem.csv", sep=";")
    for num, id in enumerate(xylem.id):
        lb[mnet.reactions_id.index(id)] = xylem.lb[num]
        ub[mnet.reactions_id.index(id)] = xylem.ub[num]


# Right-hand side of the LP problem
rhs = [0] * mnet.nb_metab


# Senses of the LP problem
sense = "E" * mnet.nb_metab


# Objective function of the LP problem (minimizations of photons import, stem and leaf)
obj = [0.0] * (len(mnet.reactions_id) + len(mnet.exchangereactions_id))


#leaf
num_R_photons_leaf = len(mnet.reactions_id) + mnet.exchangereactions_id.index('Exch_M_photon_b_l')
lb[num_R_photons_leaf] = -cplex.infinity
ub[num_R_photons_leaf] = cplex.infinity

#stem
num_R_photons_stem =len(mnet.reactions_id) + mnet.exchangereactions_id.index('Exch_M_photon_b_s')
lb[num_R_photons_stem] = -cplex.infinity
ub[num_R_photons_stem] = cplex.infinity

obj[num_R_photons_leaf] = float(-1.0 * mnet.LEAF_WEIGHT)
obj[num_R_photons_stem] = float(-1.0 * mnet.STEM_WEIGHT)

# Variables' initialisation for additional constraints of the LP problem, not on lb or ub
addnames = []
addstoichMat_rows = []
addstoichMat_values = []
addstoichMat_columns = []
addrhs = []
addsense = ""


# Constraint for relative growth rates of each organs (based on experimental data)
addnames.append("Root to Leaf relative growth rate")

addstoichMat_values.append(0.17/0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab + len(addsense))
addstoichMat_rows.append(mnet.nb_metab + len(addsense))

addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_LEAF_l'))
addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_ROOT_r'))

addrhs.append(0)

addsense = addsense + "E"


addnames.append("Stem to Leaf relative growth rate")
addstoichMat_values.append(0.26/0.21)
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab + len(addsense))
addstoichMat_rows.append(mnet.nb_metab + len(addsense))

addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_LEAF_l'))
addstoichMat_columns.append(mnet.reactions_id.index('R_BIOMASS_STEM_s'))

addrhs.append(0)

addsense = addsense + "E"

# leaf stem photosynthesis constraint
addnames.append('leaf stem photosynthesis constraint')

addstoichMat_values.append(mnet.leaf_stem_photosynthesis)
num_constraint_ratio_eff = len(addstoichMat_values) - 1
addstoichMat_values.append(-1)


addstoichMat_rows.append(mnet.nb_metab + len(addsense))
addstoichMat_rows.append(mnet.nb_metab + len(addsense))

addstoichMat_columns.append(mnet.reactions_id.index('R_EX_photon_h_s'))
addstoichMat_columns.append(mnet.reactions_id.index('R_EX_photon_h_l'))

addrhs.append(0)

addsense = addsense + "L"


# [OPTIONAL] nitrogen ratio constraint
if not np.isnan(mnet.nh4_no3):
    addnames.append("NH4 to NO3 uptake ratio")

    addstoichMat_values.append(mnet.nh4_no3)
    addstoichMat_values.append(-1)

    addstoichMat_rows.append(mnet.nb_metab + len(addsense))
    addstoichMat_rows.append(mnet.nb_metab + len(addsense))

    addstoichMat_columns.append(mnet.reactions_id.index('R_EX_no3_c_r'))
    addstoichMat_columns.append(mnet.reactions_id.index('R_EX_nh4_c_r'))

    addrhs.append(0)

    addsense = addsense + "E"

# [OPTIONAL] carboxylase oxygenase Rubisco activity ratio
if not np.isnan(mnet.ratio_Rubisco_carboxylase_oxygenase):
    addnames.append("Rubisco CO2 vs O2 ratio activity leaf")

    addstoichMat_values.append(1)
    addstoichMat_values.append(-mnet.ratio_Rubisco_carboxylase_oxygenase)

    addstoichMat_rows.append(mnet.nb_metab + len(addsense))
    addstoichMat_rows.append(mnet.nb_metab + len(addsense))

    addstoichMat_columns.append(mnet.reactions_id.index('R_RBPCh_l'))
    addstoichMat_columns.append(mnet.reactions_id.index('R_RBCh_1_l'))

    addrhs.append(0)

    addsense = addsense + "E"

    addnames.append("Rubisco CO2 vs O2 ratio activity stem")

    addstoichMat_values.append(1)
    addstoichMat_values.append(-mnet.ratio_Rubisco_carboxylase_oxygenase)

    addstoichMat_rows.append(mnet.nb_metab + len(addsense))
    addstoichMat_rows.append(mnet.nb_metab + len(addsense))

    addstoichMat_columns.append(mnet.reactions_id.index('R_RBPCh_s'))
    addstoichMat_columns.append(mnet.reactions_id.index('R_RBCh_1_s'))

    addrhs.append(0)

    addsense = addsense + "E"

#############################################################################

#                               STEM EFFICIENCY

#############################################################################

eff_reduction = [i for i in np.arange(1, 18, 0.2)]


reac_ids = ['R_RBPCh_s', 'Exch_phl_s_M_sucr_c','R_PSIh_l',
            'R_PSIh_s','R_PSII_h_s', 'R_PSII_h_l','R_TKT1h_s',
            'R_TAh_s', 'R_RPEh_r','R_PYRDC_s', 'R_ARODH_c_s',
            'R_TPIh_s', 'R_XYLI2_s', 'R_FBA_s', 'R_PFK_s',
            'R_ENO_s', 'R_FUM_s', 'R_FUMm_s', 'R_TPI_s',
            'R_ENO_h_s', 'R_HEX7_s', 'R_RPEh_s', 'R_RPEh_l',
            'R_PRUK_1_l', 'R_PRUK_1_s', 'R_PRUK_1_r', 'R_TKT1h_r',
            'R_GLUKBh_s', 'R_GLUK_s', 'R_SUCS_s', 'R_RPIh_r',
            'R_G6PBDHh_r', 'R_PGIBh_r', 'R_PFKh_r', 'Exch_phl_r_M_sucr_c',
            'Exch_l_phl_M_sucr_c', 'Exch_s_phl_M_sucr_c']

reac_fluxes = {'Ratio efficiency leaf:stem': eff_reduction, 'Photon uptake': []}

for reac_id in reac_ids:
    reac_fluxes[reac_id] = []

reac_fluxes['Sum Fluxes'] = []
reac_fluxes['Status Min Photon'] = []
reac_fluxes['Status Min Sum Fluxes'] = []

print(reac_fluxes)



for eff in eff_reduction:
    leaf_stem_photosynthesis = eff

    #############################################################################

    #                      FLUX BALANCE ANALYSIS
    #                     1) Photon uptake minimization
    #                     2) ABS(Flux) sum minimization

    #############################################################################

    addstoichMat_values[num_constraint_ratio_eff] = leaf_stem_photosynthesis

    try:
        #####################################
        # FBA 1) Photon uptake minimization
        #####################################

        probFBA = cplex.Cplex()

        probFBA.set_problem_name("FBATomatoMinPhotons")

        probFBA.objective.set_sense(probFBA.objective.sense.minimize)

        probFBA.linear_constraints.add(rhs=rhs + addrhs, senses=sense + addsense,
                                       names=["QSSA_" + id for id in mnet.metabolites_id] + addnames)

        probFBA.variables.add(obj=obj, lb=lb, ub=ub, names=mnet.reactions_id + mnet.exchangereactions_id)

        probFBA.linear_constraints.set_coefficients(
            zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows + addstoichMat_rows,
                mnet.stoichMat_columns + mnet.exchangestoichMat_columns + addstoichMat_columns,
                mnet.stoichMat_values + mnet.exchangestoichMat_values + addstoichMat_values))

        #probFBA.write("output/LP/" + "FBATomatoMinPhotons" + "_probFBA.lp")

        probFBA.solve()

        x = probFBA.solution.get_values()


        # Print results
        print("Solution status = ", probFBA.solution.get_status(), ":",)
        print(probFBA.solution.status[probFBA.solution.get_status()])
        print("Objective value  = ", probFBA.solution.get_objective_value())
        #print("Iteration count = ", probFBA.solution.progress.get_num_iterations())

        #print("NH4", x[mnet.reactions_id.index('R_EX_nh4_c_r')])
        #print("NO3", x[mnet.reactions_id.index('R_EX_no3_c_r')])
        #print("co2_l", x[mnet.reactions_id.index('R_EX_co2_e_l')])
        #print("co2_s", x[mnet.reactions_id.index('R_EX_co2_e_s')])
        #print("co2_r", x[mnet.reactions_id.index('R_EX_co2_e_r')])
        #print('PS stem', x[mnet.reactions_id.index('R_EX_photon_h_s')])
        #print('PS leaf', x[mnet.reactions_id.index('R_EX_photon_h_l')])

        #store results
        reac_fluxes['Photon uptake'].append(probFBA.solution.get_objective_value())
        reac_fluxes['Status Min Photon'].append(probFBA.solution.get_status())


        # Create new constraints for FBA 2)

        # define new bounds for FBA 2)
        lb2 = lb + [-cplex.infinity] * len(mnet.reactions_id)

        ub2 = ub + [cplex.infinity] * len(mnet.reactions_id)

        #add constraints for photon uptake from previous FBA
        lb2[mnet.reactions_id.index('R_EX_photon_h_s')] = x[mnet.reactions_id.index('R_EX_photon_h_s')]
        ub2[mnet.reactions_id.index('R_EX_photon_h_s')] = x[mnet.reactions_id.index('R_EX_photon_h_s')]

        lb2[mnet.reactions_id.index('R_EX_photon_h_l')] = x[mnet.reactions_id.index('R_EX_photon_h_l')]
        ub2[mnet.reactions_id.index('R_EX_photon_h_l')] = x[mnet.reactions_id.index('R_EX_photon_h_l')]

        # define new rhs for FBA 2)
        rhs2 = rhs + addrhs + [0] * len(mnet.reactions_id) * 2

        # define new sense for FBA 2)
        sense2 = sense + addsense + "L" * len(mnet.reactions_id) * 2

        # define new obj FBA 2)
        obj2 = [0]*(len(mnet.reactions_id) + len(mnet.exchangereactions_id)) + [1] * len(mnet.reactions_id)

        # Variables' initialisation for additional constraints of the LP problem, not on lb or ub
        minnames = ["MinFluxConstraint1_" + id for id in mnet.reactions_id] + \
                   ["MinFluxConstraint2_" + id for id in mnet.reactions_id]

        minstoichMat_rows = [index for index in range(mnet.nb_metab + len(addsense), mnet.nb_metab + len(addsense) + len(mnet.reactions_id))] + \
                            [index for index in range(mnet.nb_metab + len(addsense), mnet.nb_metab + len(addsense) + len(mnet.reactions_id))] + \
                            [index for index in range(mnet.nb_metab + len(addsense) + len(mnet.reactions_id), mnet.nb_metab + len(addsense) + len(mnet.reactions_id) * 2)] +\
                            [index for index in range(mnet.nb_metab + len(addsense) + len(mnet.reactions_id), mnet.nb_metab + len(addsense) + len(mnet.reactions_id) * 2)]
        minstoichMat_columns = [index for index in range(len(mnet.reactions_id))] + \
                               [index for index in range(len(mnet.exchangereactions_id) + len(mnet.reactions_id), len(mnet.exchangereactions_id) + len(mnet.reactions_id)*2)] + \
                               [index for index in range(len(mnet.reactions_id))] + \
                               [index for index in range(len(mnet.exchangereactions_id) + len(mnet.reactions_id), len(mnet.exchangereactions_id) + len(mnet.reactions_id)*2)]

        minstoichMat_values = [1] * len(mnet.reactions_id) + [-1] * len(mnet.reactions_id) * 3

        try:
            #####################################
            # FBA 2) ABS(Flux) sum minimization
            #####################################

            probFBA2 = cplex.Cplex()

            probFBA2.set_problem_name("FBATomatoMinFluxSum")

            probFBA2.objective.set_sense(probFBA2.objective.sense.minimize)

            probFBA2.linear_constraints.add(rhs=rhs2, senses=sense2,
                                            names=["QSSA_" + id for id in mnet.metabolites_id] + addnames + minnames)

            probFBA2.variables.add(obj=obj2, lb=lb2, ub=ub2,
                                   names=mnet.reactions_id + mnet.exchangereactions_id
                                         + ["MinVar_" + id for id in mnet.reactions_id])
            probFBA2.linear_constraints.set_coefficients(
                zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows + addstoichMat_rows + minstoichMat_rows,
                    mnet.stoichMat_columns + mnet.exchangestoichMat_columns + addstoichMat_columns + minstoichMat_columns,
                    mnet.stoichMat_values + mnet.exchangestoichMat_values + addstoichMat_values + minstoichMat_values))

            #probFBA2.write("output/LP/" + "FBATomatoMinFluxSum" + "_probFBA.lp")

            probFBA2.solve()

            x2 = probFBA2.solution.get_values()

            # Print results
            print("Solution status = ", probFBA2.solution.get_status(), ":", )
            #print(probFBA2.solution.status[probFBA2.solution.get_status()])
            print("Objective value  = ", probFBA2.solution.get_objective_value())
            #print("Iteration count = ", probFBA2.solution.progress.get_num_iterations())

            #print("NH4", x2[mnet.reactions_id.index('R_EX_nh4_c_r')])
            #print("NO3", x2[mnet.reactions_id.index('R_EX_no3_c_r')])
            #print("co2_l", x2[mnet.reactions_id.index('R_EX_co2_e_l')])
            #print("co2_s", x2[mnet.reactions_id.index('R_EX_co2_e_s')])
            #print("co2_r", x2[mnet.reactions_id.index('R_EX_co2_e_r')])
            #print('PS stem', x2[mnet.reactions_id.index('R_EX_photon_h_s')])
            #print('PS leaf', x2[mnet.reactions_id.index('R_EX_photon_h_l')])


            #store results
            reac_fluxes['Sum Fluxes'].append(probFBA2.solution.get_objective_value())
            reac_fluxes['Status Min Sum Fluxes'].append(probFBA2.solution.get_status())

            for reac_id in reac_ids:
                reac_fluxes[reac_id].append(x2[mnet.reactions_id.index(reac_id)])

        except CplexError as exc:
            print(exc)

    except CplexError as exc:
        print(exc)

# save results
df = pd.DataFrame(reac_fluxes)


df.to_excel('output/FBA_3comp_stem_efficiency.xlsx',
            sheet_name='results_stem_eff', index=False)
