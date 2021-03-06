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
addstoichMat_values.append(-1)

addstoichMat_rows.append(mnet.nb_metab + len(addsense))
addstoichMat_rows.append(mnet.nb_metab + len(addsense))

addstoichMat_columns.append(mnet.reactions_id.index('R_EX_photon_h_s'))
addstoichMat_columns.append(mnet.reactions_id.index('R_EX_photon_h_l'))

addrhs.append(0)

addsense = addsense + "L"


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


#constraint on nh4/no3 ratio
nh4_no3_values = ['no constraints', 'no NH4'] \
                 + [0.01, 0.02, 0.1, 0.25, 0.33, 1, 3, 4, 10, 50, 100]

# create variable to store results of the different simulation
reac_ids = ['R_EX_photon_h_l', 'R_EX_photon_h_s', 'R_EX_nh4_c_r', 'R_EX_no3_c_r', 'Exch_r_xyl_M_nh4_c',
            'Exch_r_xyl_M_no2_c', 'Exch_r_xyl_M_no3_c','Exch_r_xyl_M_phe_L_c', 'Exch_r_xyl_M_tyr_L_c',
            'Exch_r_xyl_M_asn_L_c','Exch_r_xyl_M_asp_L_c', 'Exch_r_xyl_M_gln_L_c', 'Exch_r_xyl_M_arg_L_c',
            'Exch_r_xyl_M_sucr_c', 'Exch_r_xyl_M_fum_c', 'Exch_r_xyl_M_thr_L_c', 'Exch_r_xyl_M_pro_L_c',
            'Exch_r_xyl_M_val_L_c', 'Exch_r_xyl_M_ile_L_c', 'Exch_r_xyl_M_leu_L_c', 'Exch_r_xyl_M_lys_L_c',
            'Exch_r_xyl_M_etoh_c', 'Exch_r_xyl_M_ala_L_c', 'Exch_r_xyl_M_glc_D_c']

reac_fluxes = {'NH4/NO3 uptake ratio': nh4_no3_values, 'Photon uptake': []}

for reac_id in reac_ids:
    reac_fluxes[reac_id] = []

reac_fluxes['Sum Fluxes'] = []
reac_fluxes['Status Min Photon'] = []
reac_fluxes['Status Min Sum Fluxes'] = []


for nh4_no3_value in nh4_no3_values:

    print()
    print(nh4_no3_value)
    print()

    if nh4_no3_value =='no constraints':
        lb[mnet.reactions_id.index('R_EX_nh4_c_r')] = -cplex.infinity
        ub[mnet.reactions_id.index('R_EX_nh4_c_r')] = cplex.infinity

        lb[mnet.reactions_id.index('R_EX_no3_c_r')] = -cplex.infinity
        ub[mnet.reactions_id.index('R_EX_no3_c_r')] = cplex.infinity

        nh4_no3_stoichMat_values = []
        nh4_no3_stoichMat_rows = []
        nh4_no3_stoichMat_columns = []
        nh4_no3_sense = ""
        nh4_no3_rhs = []
        nh4_no3_name = []


    elif nh4_no3_value == 'no NH4':
        lb[mnet.reactions_id.index('R_EX_nh4_c_r')] = 0
        ub[mnet.reactions_id.index('R_EX_nh4_c_r')] = 0

        lb[mnet.reactions_id.index('R_EX_no3_c_r')] = -cplex.infinity
        ub[mnet.reactions_id.index('R_EX_no3_c_r')] = cplex.infinity

        nh4_no3_stoichMat_values = []
        nh4_no3_stoichMat_rows = []
        nh4_no3_stoichMat_columns = []
        nh4_no3_sense = ""
        nh4_no3_rhs = []
        nh4_no3_name = []


    else:
        lb[mnet.reactions_id.index('R_EX_nh4_c_r')] = -cplex.infinity
        ub[mnet.reactions_id.index('R_EX_nh4_c_r')] = cplex.infinity

        lb[mnet.reactions_id.index('R_EX_no3_c_r')] = -cplex.infinity
        ub[mnet.reactions_id.index('R_EX_no3_c_r')] = cplex.infinity

        nh4_no3_stoichMat_values = [nh4_no3_value, -1]
        nh4_no3_stoichMat_rows = [mnet.nb_metab + len(addsense), mnet.nb_metab + len(addsense)]
        nh4_no3_stoichMat_columns = [mnet.reactions_id.index('R_EX_no3_c_r'), mnet.reactions_id.index('R_EX_nh4_c_r')]
        nh4_no3_sense = "E"
        nh4_no3_rhs = [0]
        nh4_no3_name = ['NH4 to NO3 uptake ratio']



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

        probFBA.linear_constraints.add(rhs=rhs + addrhs + nh4_no3_rhs, senses=sense + addsense + nh4_no3_sense,
                                       names=["QSSA_" + id for id in mnet.metabolites_id] + addnames
                                             + nh4_no3_name)

        probFBA.variables.add(obj=obj, lb=lb, ub=ub, names=mnet.reactions_id + mnet.exchangereactions_id)

        probFBA.linear_constraints.set_coefficients(
            zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows + addstoichMat_rows + nh4_no3_stoichMat_rows,
                mnet.stoichMat_columns + mnet.exchangestoichMat_columns + addstoichMat_columns + nh4_no3_stoichMat_columns,
                mnet.stoichMat_values + mnet.exchangestoichMat_values + addstoichMat_values + nh4_no3_stoichMat_values))

        probFBA.write("output/LP/" + "FBATomatoMinPhotons" + "_probFBA.lp")

        probFBA.solve()

        x = probFBA.solution.get_values()


        # Print results
        print("Solution status = ", probFBA.solution.get_status(), ":",)
        print(probFBA.solution.status[probFBA.solution.get_status()])
        print("Objective value  = ", probFBA.solution.get_objective_value())
        #print("Iteration count = ", probFBA.solution.progress.get_num_iterations())

        print("NH4", x[mnet.reactions_id.index('R_EX_nh4_c_r')])
        print("NO3", x[mnet.reactions_id.index('R_EX_no3_c_r')])
        #print("co2_l", x[mnet.reactions_id.index('R_EX_co2_e_l')])
        #print("co2_s", x[mnet.reactions_id.index('R_EX_co2_e_s')])
        #print("co2_r", x[mnet.reactions_id.index('R_EX_co2_e_r')])
        #print('PS stem', x[mnet.reactions_id.index('R_EX_photon_h_s')])
        #print('PS leaf', x[mnet.reactions_id.index('R_EX_photon_h_l')])

        # store results
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
        rhs2 = rhs + addrhs + nh4_no3_rhs + [0] * len(mnet.reactions_id) * 2

        # define new sense for FBA 2)
        sense2 = sense + addsense + nh4_no3_sense + "L" * len(mnet.reactions_id) * 2

        # define new obj FBA 2)
        obj2 = [0]*(len(mnet.reactions_id) + len(mnet.exchangereactions_id)) + [1] * len(mnet.reactions_id)

        # Variables' initialisation for additional constraints of the LP problem, not on lb or ub
        minnames = ["MinFluxConstraint1_" + id for id in mnet.reactions_id] + \
                   ["MinFluxConstraint2_" + id for id in mnet.reactions_id]

        minstoichMat_rows = [index for index in range(mnet.nb_metab + len(addsense) + len(nh4_no3_sense), mnet.nb_metab + len(addsense) + len(nh4_no3_sense) + len(mnet.reactions_id))] + \
                            [index for index in range(mnet.nb_metab + len(addsense) + len(nh4_no3_sense), mnet.nb_metab + len(addsense) + len(nh4_no3_sense) + len(mnet.reactions_id))] + \
                            [index for index in range(mnet.nb_metab + len(addsense) + len(nh4_no3_sense) + len(mnet.reactions_id), mnet.nb_metab + len(addsense) + len(nh4_no3_sense) + len(mnet.reactions_id) * 2)] +\
                            [index for index in range(mnet.nb_metab + len(addsense) + len(nh4_no3_sense) + len(mnet.reactions_id), mnet.nb_metab + len(addsense) + len(nh4_no3_sense) + len(mnet.reactions_id) * 2)]
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
                                            names=["QSSA_" + id for id in mnet.metabolites_id]
                                                  + addnames + nh4_no3_name + minnames)

            probFBA2.variables.add(obj=obj2, lb=lb2, ub=ub2,
                                   names=mnet.reactions_id + mnet.exchangereactions_id
                                         + ["MinVar_" + id for id in mnet.reactions_id])
            probFBA2.linear_constraints.set_coefficients(
                zip(mnet.stoichMat_rows + mnet.exchangestoichMat_rows + addstoichMat_rows
                    + nh4_no3_stoichMat_rows + minstoichMat_rows,
                    mnet.stoichMat_columns + mnet.exchangestoichMat_columns + addstoichMat_columns
                    + nh4_no3_stoichMat_columns + minstoichMat_columns,
                    mnet.stoichMat_values + mnet.exchangestoichMat_values + addstoichMat_values
                    + nh4_no3_stoichMat_values + minstoichMat_values))

            probFBA2.write("output/LP/" + "FBATomatoMinFluxSum" + "_probFBA.lp")

            probFBA2.solve()

            x2 = probFBA2.solution.get_values()

            # Print results
            print("Solution status = ", probFBA2.solution.get_status(), ":", )
            print(probFBA2.solution.status[probFBA2.solution.get_status()])
            print("Objective value  = ", probFBA2.solution.get_objective_value())
            #print("Iteration count = ", probFBA2.solution.progress.get_num_iterations())

            print("NH4", x2[mnet.reactions_id.index('R_EX_nh4_c_r')])
            print("NO3", x2[mnet.reactions_id.index('R_EX_no3_c_r')])
            #print("co2_l", x2[mnet.reactions_id.index('R_EX_co2_e_l')])
            #print("co2_s", x2[mnet.reactions_id.index('R_EX_co2_e_s')])
            #print("co2_r", x2[mnet.reactions_id.index('R_EX_co2_e_r')])
            #print('PS stem', x2[mnet.reactions_id.index('R_EX_photon_h_s')])
            #print('PS leaf', x2[mnet.reactions_id.index('R_EX_photon_h_l')])

            # store results
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


df.to_excel('output/FBA_3comp_nh4_no3_ratio.xlsx',
            sheet_name='results_stem_proportion', index=False)