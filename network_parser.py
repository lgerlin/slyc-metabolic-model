# -*- coding:Utf8 -*-

#############################################################################
# Program Python type
# authors: Gerlin et al., 2021
#############################################################################

#############################################################################
# External functions
import pandas as pd
from lxml import etree


#############################################################################
# Local functions

def save_data(data, name):
    with open('{0}.txt'.format(name), 'w') as filehandle:
        for listitem in data:
            filehandle.write('%s\n' % listitem)


#############################################################################

# SBML file
filename = 'input/Sl2184.xml'


# open calib file and assign values
calib = pd.read_csv("input/calibration.csv", sep=";", na_values="na", index_col=0)
print(calib)
print()

nh4_no3 = float(calib.ix['nh4_to_no3_ratio'].value)
leaf_stem_photosynthesis = float(calib.ix['photosynthesis'])
cost = float(calib.ix['cost'].value)
xyl = int(calib.ix['xylem_constraint'].value)
atpm_l = float(calib.ix['ATP_maintenance_leaf'].value)
atpm_r = float(calib.ix['ATP_maintenance_root'].value)
atpm_s = float(calib.ix['ATP_maintenance_stem'].value)
leaf_biomass_growth = float(calib.ix['Leaf_biomass_growth'].value)
unit = calib.ix['unit'].value
ratio_Rubisco_carboxylase_oxygenase = float(calib.ix['Rubisco_carboxylase_oxygenase_ratio'].value)

# organs relative weight from calib file
LEAF_WEIGHT = float(calib.ix['leaf_relative_weight'].value)
STEM_WEIGHT = float(calib.ix['stem_relative_weight'].value)
ROOT_WEIGHT = float(calib.ix['root_relative_weight'].value)


#############################################################################

#   PARSE THE SBML FILE AND GENERATE 3 ORGAN MODEL

#############################################################################


# Exchange compartments
organs = ['root', 'leaf', 'stem']
organs_id = ['r', 'l', 's']
nb_organs = len(organs_id)
exchange_compartments_names = ['xylem', 'phloem']
exchange_compartments_id = ['xyl', 'phl']
exchange_compartments_metabolites = {
    'xyl': ['M_phe_L_c', 'M_tyr_L_c', 'M_asn_L_c', 'M_asp_L_c',
            'M_gln_L_c', 'M_arg_L_c', 'M_sucr_c', 'M_glc_D_c',
            'M_thr_L_c', 'M_pro_L_c', 'M_val_L_c', 'M_so4_c',
            'M_ile_L_c', 'M_leu_L_c', 'M_lys_L_c', 'M_etoh_c',
            'M_no3_c', 'M_nh4_c', 'M_pi_c', 'M_ca2_c', 'M_mg2_c',
            'M_ala_L_c', 'M_fum_c', 'M_h2o_c', 'M_k_c', 'M_no2_c',
            'M_na_c', 'M_cl_c'
            ],
    'phl': ['M_glc_D_c', 'M_sucr_c', 'M_pro_L_c', 'M_mal_L_c',
            'M_4abut_c', 'M_cit_c', 'M_quin_c']}

exchange_compartments_direction_RL = {'xyl': ['r', 'l'],
                                      'phl': ['l', 'r']
                                      }

exchange_compartments_direction_S = {'xyl': ['s'],
                                     'phl': ['s']
                                     }

organs_weights = {'r': ROOT_WEIGHT, 'l': LEAF_WEIGHT}

PLANT_WEIGHT = LEAF_WEIGHT + STEM_WEIGHT + ROOT_WEIGHT


# create a tree by parsing the xml file
tree = etree.parse(filename)

# get the root of the tree
root = tree.getroot()

# get and parse compartments
compartments_name = []
compartments_id = []

listOfCompartments = root[0][4]
print('number of compartments in genome scale network: ', len(listOfCompartments))

for compartment in listOfCompartments:
    compartments_name.append(compartment.get('name'))
    compartments_id.append(compartment.get('id'))

# get genes nodes and print statistics
listOfGenes = root[0][2]
print('number of genes in single organ genome scale network: ', len(listOfGenes))

# get metabolites nodes and print statistics
listOfSpecies = root[0][5]
print('number of metabolites in single organ genome scale network: ', len(listOfSpecies))

# get reactions nodes and print statistics
listOfReactions = root[0][7]
print('number of reactions in single organ genome scale network: ', len(listOfReactions))

# get the metabolites and reactions of the 1 organ network
metab_id_1comp = []
react_id_1comp = []

for metabolite in listOfSpecies:
    metab_id_1comp.append(metabolite.get('id'))

for reaction in listOfReactions:
    react_id_1comp.append(reaction.get('id'))



# generation of the 3 comp network
# generation of metabolites
metabolites_name = []
metabolites_id = []
metabolites_compartment = []
metabolites_boundaryCondition = []
list_boundary = []

for metabolite in listOfSpecies:
    for i, organ in enumerate(organs):
        metab_name = metabolite.get('name') + '_' + organ
        metab_id = metabolite.get('id') + '_' + organs_id[i]
        metabolites_name.append(metab_name)
        metabolites_id.append(metab_id)

        if metabolite.get('boundaryCondition') == 'true':
            metabolites_boundaryCondition.append(1)
            list_boundary.append(metab_id)
        else:
            metabolites_boundaryCondition.append(0)


# metabolites in exchange compartment (xylem, phloeme)
for i, comp in enumerate(exchange_compartments_id):
    for metab in exchange_compartments_metabolites[comp]:
        metabolites_id.append(metab[0:-1] + comp)
        metabolites_name.append(metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                + exchange_compartments_names[i])



# generation of reactions
reactions_name = []
reactions_id = []
reactions_reversible = []
reactions_formulas = []
reactions_formulasNames = []

stoichMat_rows = []
stoichMat_columns = []
stoichMat_values = []


for ind, reaction in enumerate(listOfReactions):

    for i, organ in enumerate(organs):
        reaction_formula = ''
        reaction_formulaName = ''
        reactions_name.append(reaction.get('name') + '_' + organ)
        reactions_id.append(reaction.get('id') + '_' + organs_id[i])

        listOfReactants = reaction[3]

        for substrate in listOfReactants:
            stoichMat_rows.append(metabolites_id.index(substrate.get('species') + '_' + organs_id[i]))
            stoichMat_columns.append(nb_organs * ind + i)
            if 'R_BIOMASS' in reaction.get('id') and unit=='mM':
                stoichMat_values.append(- float(substrate.get('stoichiometry'))/1000)
            else:
                stoichMat_values.append(- float(substrate.get('stoichiometry')))

            reaction_formula = reaction_formula + \
                               substrate.get('stoichiometry') + \
                               ' ' + substrate.get('species') + '_' + \
                               organs_id[i] + ' + '

            reaction_formulaName = reaction_formulaName + \
                                       substrate.get('stoichiometry') + \
                                       ' ' + metabolites_name[
                                           metabolites_id.index(
                                               substrate.get('species') +
                                               '_' + organs_id[i])] + ' + '

        reaction_formula = reaction_formula[0:-2]
        reaction_formulaName = reaction_formulaName[0:-2]

        if reaction.get('reversible') == 'true':

            reactions_reversible.append(1)
            reaction_formula = reaction_formula + '<--> '
            reaction_formulaName = reaction_formulaName + '<--> '

        else:
            reactions_reversible.append(0)
            reaction_formula = reaction_formula + '--> '
            reaction_formulaName = reaction_formulaName + '--> '

        listOfProducts = reaction[4]

        for product in listOfProducts:
            stoichMat_rows.append(metabolites_id.index(product.get('species') + '_' + organs_id[i]))

            stoichMat_columns.append(nb_organs * ind + i)

            stoichMat_values.append(float(product.get('stoichiometry')))

            reaction_formula = reaction_formula + ' ' + \
                               product.get('stoichiometry') + ' ' + \
                               product.get('species') + '_' + \
                               organs_id[i] + ' + '

            reaction_formulaName = reaction_formulaName + ' ' + \
                                       product.get(
                                           'stoichiometry') + ' ' + \
                                       metabolites_name[
                                           metabolites_id.index(
                                               product.get(
                                                   'species'
                                               ) + '_' + organs_id[i]
                                           )
                                       ] + ' + '

        reaction_formula = reaction_formula[0:-3]
        reaction_formulaName = reaction_formulaName[0:-3]

        reactions_formulas.append(reaction_formula)
        reactions_formulasNames.append(reaction_formulaName)


# Exchanges between compartments
nb_reac = len(reactions_id)

atp_cost_list = []
indexes_constraint_stem_weight = []
indexes_constraint_leaf_root_weight = []

for i, id in enumerate(exchange_compartments_id):

    for j, metab in enumerate(exchange_compartments_metabolites[id]):

        # Reaction names and reversibility
        # organs -> exchange compartement (xylem or phloem)
        reactions_name.append('Exchange_' + organs[organs_id.index(exchange_compartments_direction_RL[id][0])]
                              + '_' + exchange_compartments_names[i] + '_' + metab)

        reactions_id.append('Exch_' + exchange_compartments_direction_RL[id][0]
                            + '_' + id + '_' + metab)
        reactions_reversible.append(0)

        # exchange compartement (xylem or phloem) -> organs
        reactions_name.append('Exchange_' + exchange_compartments_names[i]
                              + '_' + organs[organs_id.index(exchange_compartments_direction_RL[id][1])]
                              + '_' + metab)
        reactions_id.append('Exch_' + id + '_' + exchange_compartments_direction_RL[id][1]
                            + '_' + metab)
        reactions_reversible.append(0)

        # Stoichiometric matrix
        # organs -> exchange compartement (xylem or phloem)
        if metab !='M_h2o_c':
            stoichMat_rows.append(metabolites_id.index('M_atp_cost_b_' + exchange_compartments_direction_RL[id][0]))
            stoichMat_columns.append(nb_reac)
            stoichMat_values.append(-cost)
            atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + id))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(organs_weights[exchange_compartments_direction_RL[id][0]] / PLANT_WEIGHT)
        indexes_constraint_leaf_root_weight.append(len(stoichMat_values) - 1)


        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + 'c_' + exchange_compartments_direction_RL[id][0]))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(-1)

        # exchange compartement (xylem or phloem) -> organs
        if metab != 'M_h2o_c':
            stoichMat_rows.append(metabolites_id.index('M_atp_cost_b_' + exchange_compartments_direction_RL[id][1]))
            stoichMat_columns.append(nb_reac + 1)
            stoichMat_values.append(-cost)
            atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + id))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(-organs_weights[exchange_compartments_direction_RL[id][1]] / PLANT_WEIGHT)
        indexes_constraint_leaf_root_weight.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + 'c_' + exchange_compartments_direction_RL[id][1]))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(1)

        nb_reac = nb_reac + 2

        # Reaction Formula
        # organs -> exchange compartement (xylem or phloem)
        if metab !='M_h2o_c':
            reactions_formulas.append('1 ' + metab + '_' + exchange_compartments_direction_RL[id][0] + ' + '
                                  + str(cost) + ' M_atp_cost_b_' + exchange_compartments_direction_RL[id][0]
                                  + ' --> '
                                  + str(organs_weights[exchange_compartments_direction_RL[id][0]] /PLANT_WEIGHT)
                                  + ' ' + metab[0:-1] + id)

            reactions_formulasNames.append('1 ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                       + organs[organs_id.index(exchange_compartments_direction_RL[id][0])]
                                       + ' + ' + str(cost)
                                       + ' ATPCost_Cyto_'
                                       + organs[organs_id.index(exchange_compartments_direction_RL[id][0])]
                                       + ' --> '
                                       + str(organs_weights[exchange_compartments_direction_RL[id][0]] /PLANT_WEIGHT)
                                       + ' ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                       + exchange_compartments_names[exchange_compartments_id.index(id)])
        else:
            reactions_formulas.append('1 ' + metab + '_' + exchange_compartments_direction_RL[id][0]
                                      + ' --> '
                                      + str(organs_weights[exchange_compartments_direction_RL[id][0]] / PLANT_WEIGHT)
                                      + ' ' + metab[0:-1] + id)

            reactions_formulasNames.append('1 ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + organs[organs_id.index(exchange_compartments_direction_RL[id][0])]
                                           + ' --> '
                                           + str(
                organs_weights[exchange_compartments_direction_RL[id][0]] / PLANT_WEIGHT)
                                           + ' ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + exchange_compartments_names[exchange_compartments_id.index(id)])

        # exchange compartement (xylem or phloem) -> organs
        if metab != 'M_h2o_c':
            reactions_formulas.append(str(organs_weights[exchange_compartments_direction_RL[id][1]] /PLANT_WEIGHT)
                                      + ' ' + metab[0:-1] + id + ' + '
                                      + str(cost) + ' M_atp_cost_b_' + exchange_compartments_direction_RL[id][1] +
                                      ' --> 1 ' + metab + '_' + exchange_compartments_direction_RL[id][1])

            reactions_formulasNames.append(str(organs_weights[exchange_compartments_direction_RL[id][1]] /PLANT_WEIGHT)
                                           + ' ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + exchange_compartments_names[exchange_compartments_id.index(id)]
                                           + ' + ' + str(cost)
                                           + ' ATPCost_Cyto_'
                                           + organs[organs_id.index(exchange_compartments_direction_RL[id][1])]
                                           + ' --> 1 ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + organs[organs_id.index(exchange_compartments_direction_RL[id][1])])
        else:
            reactions_formulas.append(str(organs_weights[exchange_compartments_direction_RL[id][1]] / PLANT_WEIGHT)
                                      + ' ' + metab[0:-1] + id +
                                      ' --> 1 ' + metab + '_' + exchange_compartments_direction_RL[id][1])

            reactions_formulasNames.append(str(organs_weights[exchange_compartments_direction_RL[id][1]] / PLANT_WEIGHT)
                                           + ' ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + exchange_compartments_names[exchange_compartments_id.index(id)]
                                           + ' --> 1 ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + organs[organs_id.index(exchange_compartments_direction_RL[id][1])])


# Additional exchanges with stem (uptake and export can both happen for both phloem and xylem)
for i, id in enumerate(exchange_compartments_id):
    for j, metab in enumerate(exchange_compartments_metabolites[id]):

        # Reaction names and reversibility
        # stem -> exchange compartement (xylem or phloem)

        reactions_name.append('Exchange_' + organs[organs_id.index(exchange_compartments_direction_S[id][0])]
                              + '_' + exchange_compartments_names[i] + '_' + metab)
        reactions_id.append('Exch_' + exchange_compartments_direction_S[id][0] + '_' + id + '_' + metab)
        reactions_reversible.append(0)

        # exchange compartement (xylem or phloem) -> stem
        reactions_name.append('Exchange_' + exchange_compartments_names[i] + '_'
                              + organs[organs_id.index(exchange_compartments_direction_S[id][0])] + '_' + metab)
        reactions_id.append('Exch_' + id + '_' + exchange_compartments_direction_S[id][0] + '_' + metab)
        reactions_reversible.append(0)

        # Stoichiometric matrix
        # stem -> exchange compartement (xylem or phloem)

        if metab != 'M_h2o_c':
            stoichMat_rows.append(metabolites_id.index('M_atp_cost_b_' + exchange_compartments_direction_S[id][0]))
            stoichMat_columns.append(nb_reac)
            stoichMat_values.append(-cost)
            atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + id))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(STEM_WEIGHT/PLANT_WEIGHT)
        indexes_constraint_stem_weight.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(metabolites_id.index(metab + '_' + exchange_compartments_direction_S[id][0]))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(-1.0)


        # exchange compartement (xylem or phloem) -> stem
        if metab != 'M_h2o_c':
            stoichMat_rows.append(metabolites_id.index('M_atp_cost_b_' + exchange_compartments_direction_S[id][0]))
            stoichMat_columns.append(nb_reac + 1)
            stoichMat_values.append(-cost)
            atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + id))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(-STEM_WEIGHT/PLANT_WEIGHT)
        indexes_constraint_stem_weight.append(len(stoichMat_values) - 1)


        stoichMat_rows.append(metabolites_id.index(metab + '_' + exchange_compartments_direction_S[id][0]))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(1.0)

        nb_reac = nb_reac + 2

        # Reaction Formula
        # stem -> exchange compartement (xylem or phloem)
        if metab != 'M_h2o_c':
            reactions_formulas.append('1 ' + metab + '_' +
                exchange_compartments_direction_S[id][0] + ' + ' +
                str(cost) + ' M_atp_cost_b_' + exchange_compartments_direction_S[id][0] +
                ' --> ' + str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab[0:-1] + id)

            reactions_formulasNames.append('1 ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + organs[organs_id.index(exchange_compartments_direction_S[id][0])]
                                           + ' + ' + str(cost)
                                           + ' ATPCost_Cyto_'
                                           + organs[organs_id.index(exchange_compartments_direction_S[id][0])]
                                           + ' --> ' + str(STEM_WEIGHT / PLANT_WEIGHT) + ' '
                                           + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + exchange_compartments_names[
                    exchange_compartments_id.index(id)])
        else:
            reactions_formulas.append('1 ' + metab + '_' + exchange_compartments_direction_S[id][0] +
                                      ' --> ' + str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab[0:-1] + id)

            reactions_formulasNames.append('1 ' + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + organs[organs_id.index(exchange_compartments_direction_S[id][0])]
                                           + ' --> ' + str(STEM_WEIGHT / PLANT_WEIGHT) + ' '
                                           + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + exchange_compartments_names[
                                               exchange_compartments_id.index(id)])

        # exchange compartement (xylem or phloem) -> stem
        if metab != 'M_h2o_c':
            reactions_formulas.append(str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab[0:-1] + id + ' + ' +
                                      str(cost) + ' M_atp_cost_b_' + exchange_compartments_direction_S[id][0] +
                                      ' --> 1 ' +
                                      metab + '_' + exchange_compartments_direction_S[id][0])

            reactions_formulasNames.append(str(STEM_WEIGHT / PLANT_WEIGHT) + ' '
                                           + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + exchange_compartments_names[exchange_compartments_id.index(id)] +
                                           ' + ' + str(cost) + ' ATPCost_Cyto_'
                                           + organs[organs_id.index(exchange_compartments_direction_S[id][0])]
                                           + ' --> 1 '
                                           + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + organs[organs_id.index(exchange_compartments_direction_S[id][0])])
        else:
            reactions_formulas.append(str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab[0:-1] + id + ' --> 1 ' +
                                      metab + '_' + exchange_compartments_direction_S[id][0])

            reactions_formulasNames.append(str(STEM_WEIGHT / PLANT_WEIGHT) + ' '
                                           + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + exchange_compartments_names[exchange_compartments_id.index(id)]
                                           + ' --> 1 '
                                           + metabolites_name[metabolites_id.index(metab + '_l')][0:-4]
                                           + organs[organs_id.index(exchange_compartments_direction_S[id][0])])


# reactions for boundary metabolites
exchangereactions_id = []
exchangereactions_name = []
exchangereactions_formulas = []
exchangereactions_formulasNames = []
exchangestoichMat_rows = []
exchangestoichMat_columns = []
exchangestoichMat_values = []

boundary_metabolites = [i for i, e in enumerate(metabolites_boundaryCondition) if e == 1]
nb_exchangereac = len(boundary_metabolites)

for index, metab in enumerate(boundary_metabolites):
    exchangestoichMat_columns.append(nb_reac + index)
    exchangestoichMat_values.append(-1)
    exchangereactions_id.append('Exch_' + metabolites_id[metab])
    exchangereactions_name.append('Exchange of ' + metabolites_name[metab])
    exchangereactions_formulas.append(metabolites_id[metab] + ' <--> ')
    exchangereactions_formulasNames.append(
        metabolites_name[metab] + ' <--> ')
    exchangestoichMat_rows.append(metab)

nb_exchangereac = len(exchangereactions_id)

#get some statistics
nb_metab = len(metabolites_id)
nb_reac = len(reactions_id)

print('number of reactions in three compartment network: ', nb_reac)
print('number of metabolites in three compartment network: ', nb_metab)

# saved the parsed network file
save_data(reactions_name, "output/parser/reactions_name")
save_data(exchangereactions_name, "output/parser/exchangereactions_name")
save_data(reactions_reversible, "output/parser/reactions_reversible")
save_data(boundary_metabolites, "output/parser/boundary_metabolites")
save_data(stoichMat_values, "output/parser/coefficients")
save_data(stoichMat_rows, "output/parser/rows")
save_data(stoichMat_columns, "output/parser/columns")
save_data(exchangestoichMat_columns, "output/parser/exchangestoichMat_columns")
save_data(exchangestoichMat_rows, "output/parser/exchangestoichMat_rows")
save_data(exchangereactions_id, "output/parser/exchangereactions_id")
save_data(metabolites_id, "output/parser/metabolites_id")
save_data(metabolites_name, "output/parser/metabolites_name")
save_data(exchangestoichMat_values, "output/parser/exchangestoichMat_values")
save_data(reactions_formulas, "output/parser/reactions_formulas")
save_data(reactions_formulasNames, "output/parser/reactions_formulasNames")
save_data(reactions_id, "output/parser/reactions_id")
save_data(exchangereactions_formulasNames, "output/parser/exchangereactions_formulasNames")
save_data(exchangereactions_formulas, "output/parser/exchangereactions_formulas")
save_data(list_boundary, 'output/parser/boundary')
save_data(react_id_1comp, "output/parser/react_id_1comp")
save_data(metab_id_1comp, "output/parser/metab_id_1comp")
