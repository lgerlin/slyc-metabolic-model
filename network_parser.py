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

# GEM file
filename = 'input/Sl2186.xml'

# open calib file and assign values

calib = pd.read_csv("input/calibration.csv", sep=";")
calib_names = calib['calibration'].values.tolist()
NH4_NO3 = calib.value[calib_names.index('nitrogen')]
LEAF_STEM_PHOTOSYNTHESIS = calib.value[calib_names.index('photosynthesis')]

cost = float(calib.value[calib_names.index('cost')])
xyl = int(calib.value[calib_names.index('xylem_constraint')])

#############################################################################

#   PARSE THE GEM FILE AND GENERATE 3 ORGAN MODEL

#############################################################################


# EXCHANGE COMPARTMENTS
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

exchange_compartments_direction_RL = {
    'xyl': ['r', 'l'] + [0] * len(
        exchange_compartments_metabolites["xyl"]),
    'phl': ['l', 'r'] + [0] * len(
        exchange_compartments_metabolites["phl"])}

exchange_compartments_direction_S = {
    'xyl': ['s'] + [0] * len(
        exchange_compartments_metabolites["xyl"]),
    'phl': ['s'] + [0] * len(
        exchange_compartments_metabolites["phl"])}

# calib

LEAF_WEIGHT = calib.value[calib_names.index('leaf_relative_weight')]
STEM_WEIGHT = calib.value[calib_names.index('stem_relative_weight')]
ROOT_WEIGHT = calib.value[calib_names.index('root_relative_weight')]

organs_weights = {'r': ROOT_WEIGHT, 'l': LEAF_WEIGHT}

PLANT_WEIGHT = LEAF_WEIGHT + STEM_WEIGHT + ROOT_WEIGHT
organs_weights = {'r': ROOT_WEIGHT, 'l': LEAF_WEIGHT}
PLANT_WEIGHT = LEAF_WEIGHT + STEM_WEIGHT + ROOT_WEIGHT

# create a tree by parsing the xml file
tree = etree.parse(filename)

# get the root of the tree
root = tree.getroot()

# get compartments
compartments_name = []
compartments_id = []

listOfGenes = root[0][2]
print('number of genes in genome scale network: ', len(listOfGenes))

listOfCompartments = root[0][4]
print('number of compartments in genome scale network: ', len(listOfCompartments))

for compartment in listOfCompartments:
    compartments_name.append(compartment.get('name'))
    compartments_id.append(compartment.get('id'))

# get metabolites
listOfSpecies = root[0][5]
print('number of metabolites in genome scale network: ', len(listOfSpecies))

# get reactions
listOfReactions = root[0][7]
print('number of reactions in genome scale network: ', len(listOfReactions))

metabolites_name = []
metabolites_id = []
metabolites_compartment = []
metabolites_boundaryCondition = []
species = []
nb_metab = len(listOfSpecies)
atp_cost_list = []
list_boundary = []

# data in the network

metab_id_1comp = []
react_id_1comp = []

for metabolite in listOfSpecies:
    metab_id_1comp.append(metabolite.get('id'))

for reaction in listOfReactions:
    react_id_1comp.append(reaction.get('id'))



# generation of the 3 comp network

for metabolite in listOfSpecies:
    for i, organ in enumerate(organs):
        metab_name = metabolite.get('name') + '_' + organ
        metab_id = metabolite.get('id') + '_' + organs_id[
            organs.index(organ)]
        metabolites_name.append(metab_name)
        metabolites_id.append(metab_id)
        if metabolite.get('boundaryCondition') == 'true':
            metabolites_boundaryCondition.append(1)
            list_boundary.append(metab_id)
        else:
            metabolites_boundaryCondition.append(0)

# metabolites in exchange compartment

for i, comp in enumerate(exchange_compartments_id):
    for metab in exchange_compartments_metabolites[comp]:
        metabolites_name.append(
            metab[0:-1] + exchange_compartments_names[i])
        metabolites_id.append(metab[0:-1] + comp)

nb_metab = len(metabolites_id)

nb_reac = len(listOfReactions)

reactions_name = []
reactions_id = []
reactions_reversible = []
reactions_formulas = []
reactions_formulasNames = []

stoichMat_rows = []
stoichMat_columns = []
stoichMat_values = []

reactions_exch_comp = []

index_plant_stem = []
index_organ_stem = []

for ind, reaction in enumerate(listOfReactions):

    for i, organ in enumerate(organs):

        reaction_formula = ''
        reaction_formulaName = ''
        reactions_name.append(reaction.get('name') + '_' + organ)
        reactions_id.append(reaction.get('id') + '_' + organs_id[i])

        listOfReactants = reaction[3]

        for substrate in listOfReactants:

            try:
                stoichMat_rows.append(metabolites_id.index
                                      (substrate.get('species') +
                                       '_' + organs_id[i]))

            except:
                stoichMat_rows.append(
                    metabolites_id.index(substrate.get('species')))

            stoichMat_columns.append(nb_organs * ind + i)
            stoichMat_values.append(
                - float(substrate.get('stoichiometry')))

            reaction_formula = reaction_formula + \
                               substrate.get('stoichiometry') + \
                               ' ' + substrate.get('species') + '_' + \
                               organs_id[i] + ' + '

            try:
                reaction_formulaName = reaction_formulaName + \
                                       substrate.get('stoichiometry') + \
                                       ' ' + metabolites_name[
                                           metabolites_id.index(
                                               substrate.get('species') +
                                               '_' + organs_id[i])] + ' + '
            except:
                reaction_formulaName = reaction_formulaName + \
                                       substrate.get('stoichiometry') + \
                                       ' ' + metabolites_name[
                                           metabolites_id.index(
                                               substrate.get(
                                                   'species'))] + \
                                       ' + '

        reaction_formula = reaction_formula[0:-2]
        reaction_formulaName = reaction_formulaName[0:-2]

        if reaction.get('reversible') == 'true':

            reactions_reversible.append(1)
            reaction_formula = reaction_formula + '<--> '
            reaction_formulaName = reaction_formulaName + '<--> '

        else:
            # for j in enumerate(organs):
            reactions_reversible.append(0)
            reaction_formula = reaction_formula + '-->'
            reaction_formulaName = reaction_formulaName + '-->'

        listOfProducts = reaction[4]

        for product in listOfProducts:
            try:
                stoichMat_rows.append(metabolites_id.
                                      index(product.get('species') +
                                            '_' + organs_id[i]))

            except:
                stoichMat_rows.append(metabolites_id.
                                      index(product.get('species')))

            stoichMat_columns.append(nb_organs * ind + i)
            stoichMat_values.append(float(product.get('stoichiometry')))

            reaction_formula = reaction_formula + ' ' + \
                               product.get('stoichiometry') + ' ' + \
                               product.get('species') + '_' + \
                               organs_id[i] + ' + '

            try:
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
            except:
                reaction_formulaName = reaction_formulaName + ' ' + \
                                       product.get(
                                           'stoichiometry'
                                       ) + ' ' + \
                                       metabolites_name[
                                           metabolites_id.index(
                                               product.get(
                                                   'species'
                                               ))
                                       ] + ' + '

        reaction_formula = reaction_formula[0:-3]
        reaction_formulaName = reaction_formulaName[0:-3]

        reactions_formulas.append(reaction_formula)
        reactions_formulasNames.append(reaction_formulaName)

# EXCHANGES BETWEEN COMPARTMENTS
nb_reac = len(reactions_id)

for i, id in enumerate(exchange_compartments_id):

    for j, metab in enumerate(exchange_compartments_metabolites[id]):

        # REAC NAMES

        # ORGANS -> EXCH COMP

        reactions_name.append(
            'Exchange_' + organs[organs_id.index(
                exchange_compartments_direction_RL[id][0])] + '_' +
            exchange_compartments_names[i] + '_' + metab)

        reactions_id.append(
            'Exch_' + exchange_compartments_direction_RL[id][0] +
            '_' + id + '_' + metab)
        reactions_reversible.append(
            exchange_compartments_direction_RL[id][j + 2])

        # EXCH COMP -> ORGANS

        reactions_name.append(
            'Exchange_' + exchange_compartments_names[i] + '_' + organs[
                organs_id.index(
                    exchange_compartments_direction_RL[id][1])
            ] +
            '_' + metab)
        reactions_id.append(
            'Exch_' + id + '_' + exchange_compartments_direction_RL[
                id][1] +
            '_' + metab)
        reactions_reversible.append(
            exchange_compartments_direction_RL[id][j + 2])

        # ORGANS -> EXCH COMP

        stoichMat_rows.append(metabolites_id.index(
            'M_atp_cost_b_' + exchange_compartments_direction_RL[id][0]))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(-cost)
        atp_cost_list.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + id))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(
            organs_weights[
                exchange_compartments_direction_RL[id][0]] / PLANT_WEIGHT
        )

        index_organ_stem.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(
            metabolites_id.index(
                metab[0:-1] + 'c_' +
                exchange_compartments_direction_RL[id][0]))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(-1)

        # EXCH COMP -> ORGANS

        stoichMat_rows.append(metabolites_id.index(
            'M_atp_cost_b_' + exchange_compartments_direction_RL[id][1]))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(-cost)
        atp_cost_list.append(len(stoichMat_values) - 1)
        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + id))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(
            -organs_weights[
                exchange_compartments_direction_RL[id][1]] / PLANT_WEIGHT)
        index_organ_stem.append(len(stoichMat_values) - 1)

        stoichMat_rows.append(
            metabolites_id.index(metab[0:-1] + 'c_' +
                                 exchange_compartments_direction_RL[
                                     id][1]))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(1)

        nb_reac = nb_reac + 2

        # REACTION FORMULA

        if exchange_compartments_direction_RL[id][j + 2] == 0:  # irrev
            # ORGANS -> EXCH COMP
            reactions_formulas.append(
                str(1) + ' ' + metab + '_' +
                exchange_compartments_direction_RL[id][0] + ' --> ' +
                str(organs_weights[
                        exchange_compartments_direction_RL[id][0]] /
                    PLANT_WEIGHT) + ' ' + metab[0:-1] + id)
            reactions_formulasNames.append(
                str(1) + ' ' + metab + '_' +
                exchange_compartments_direction_RL[id][0] + ' --> ' +
                str(organs_weights[
                        exchange_compartments_direction_RL[id][0]] /
                    PLANT_WEIGHT) + ' ' + metab[0:-1] +
                exchange_compartments_names[
                    exchange_compartments_id.index(id)])

            # EXCH COMP -> ORGANS
            reactions_formulas.append(
                str(organs_weights[
                        exchange_compartments_direction_RL[id][1]] /
                    PLANT_WEIGHT) + ' ' + metab[0:-1] + id + ' --> 1 ' +
                metab + '_' + exchange_compartments_direction_RL[id][1])
            reactions_formulasNames.append(
                str(organs_weights[
                        exchange_compartments_direction_RL[id][1]] /
                    PLANT_WEIGHT) + ' ' + metab[0:-1] +
                exchange_compartments_names[
                    exchange_compartments_id.index(id)] +
                ' --> 1 ' + metab + '_' + organs[
                    organs_id.index(
                        exchange_compartments_direction_RL[id][1])])

        else:  # reversibility
            # ORGANS -> EXCH COMP
            reactions_formulas.append(
                str(1) + ' ' + metab + '_' +
                exchange_compartments_direction_RL[id][0] + ' <--> 1 ' +
                str(organs_weights[
                        exchange_compartments_direction_RL[id][0]] /
                    PLANT_WEIGHT) + ' ' + metab[0:-1] + id)
            reactions_formulasNames.append(
                str(1) + ' ' + metab + '_' +
                exchange_compartments_direction_RL[id][0] + ' <--> ' +
                str(
                    organs_weights[
                        exchange_compartments_direction_RL[id][0]] /
                    PLANT_WEIGHT) + ' ' + metab[0:-1] +
                exchange_compartments_names[
                    exchange_compartments_id.index(id)])

            # EXCH COMP -> ORGANS
            reactions_formulas.append(
                str(organs_weights[
                        exchange_compartments_direction_RL[id][1]] /
                    PLANT_WEIGHT) + ' ' + metab[0:-1] + id +
                ' <--> 1 ' + metab + '_' +
                exchange_compartments_direction_RL[id][1])
            reactions_formulasNames.append(
                str(organs_weights[
                        exchange_compartments_direction_RL[id][1]] /
                    PLANT_WEIGHT) + ' ' + metab[0:-1] +
                exchange_compartments_names[
                    exchange_compartments_id.index(id)] + ' <--> 1 ' +
                ' ' + metab + '_' + organs[
                    organs_id.index(
                        exchange_compartments_direction_RL[id][1])])

# EXCHANGES WITH STEM

for i, id in enumerate(exchange_compartments_id):
    for j, metab in enumerate(exchange_compartments_metabolites[id]):

        # REAC NAMES

        # ORGANS -> EXCH

        reactions_name.append(
            'Exchange_' + organs[
                organs_id.index(
                    exchange_compartments_direction_S[id][0])] + '_' +
            exchange_compartments_names[i] + '_' + metab)
        reactions_id.append(
            'Exch_' + exchange_compartments_direction_S[id][0] +
            '_' + id + '_' + metab)
        reactions_reversible.append(
            exchange_compartments_direction_S[id][j + 1])

        # EXCH -> ORGANS

        reactions_name.append(
            'Exchange_' + exchange_compartments_names[i] + '_' +
            organs[
                organs_id.index(
                    exchange_compartments_direction_S[id][0])] +
            '_' + metab)
        reactions_id.append(
            'Exch_' + id + '_' + exchange_compartments_direction_S[id][0] +
            '_' + metab)
        reactions_reversible.append(
            exchange_compartments_direction_S[id][j + 1])

        # STOICH MAT

        # ORGANS -> EXCH COMP

        stoichMat_rows.append(
            metabolites_id.index(
                'M_atp_cost_b_' + exchange_compartments_direction_S[
                    id][0]))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(-cost)
        atp_cost_list.append(len(stoichMat_values) - 1)
        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + id))
        stoichMat_columns.append(nb_reac)
        stoichMat_values.append(1)

        stoichMat_values.append(-1.0 * PLANT_WEIGHT / STEM_WEIGHT)
        index_plant_stem.append(len(stoichMat_values) - 1)
        stoichMat_rows.append(
            metabolites_id.index(
                metab + '_' + exchange_compartments_direction_S[id][0]))
        stoichMat_columns.append(nb_reac)

        # EXCH COMP -> ORGANS

        stoichMat_rows.append(
            metabolites_id.index(
                'M_atp_cost_b_' + exchange_compartments_direction_S[
                    id][0]))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(-cost)
        atp_cost_list.append(len(stoichMat_values) - 1)
        stoichMat_rows.append(metabolites_id.index(metab[0:-1] + id))
        stoichMat_columns.append(nb_reac + 1)
        stoichMat_values.append(-1)

        stoichMat_values.append(1.0 * PLANT_WEIGHT / STEM_WEIGHT)
        index_plant_stem.append(len(stoichMat_values) - 1)
        stoichMat_rows.append(
            metabolites_id.index(
                metab + '_' + exchange_compartments_direction_S[id][0]))
        stoichMat_columns.append(nb_reac + 1)

        nb_reac = nb_reac + 2

        # REACTION FORMULA

        if exchange_compartments_direction_S[id][j + 1] == 0:
            reactions_formulas.append(
                str(PLANT_WEIGHT / STEM_WEIGHT) + ' ' + metab + '_' +
                exchange_compartments_direction_S[id][0] + ' --> 1 ' +
                metab[0:-1] + id)
            reactions_formulasNames.append(
                str(
                    PLANT_WEIGHT / STEM_WEIGHT
                ) + ' ' + metab + '_' + organs[
                    organs_id.index(
                        exchange_compartments_direction_S[id][0])] +
                ' --> 1 ' + metab[0:-1] + exchange_compartments_names[
                    exchange_compartments_id.index(id)])

            reactions_formulas.append(
                '1 ' + metab[0:-1] + id + ' --> ' +
                str(PLANT_WEIGHT / STEM_WEIGHT) + ' ' +  # modif
                metab + '_' + exchange_compartments_direction_S[id][0])

            reactions_formulasNames.append(
                metab[0:-1] +
                exchange_compartments_names[
                    exchange_compartments_id.index(id)] + ' -->' +
                str(PLANT_WEIGHT / STEM_WEIGHT) + ' ' + metab + '_' +
                organs[organs_id.index(
                    exchange_compartments_direction_S[id][0])])

        else:  # reversibility
            reactions_formulas.append(
                str(PLANT_WEIGHT / STEM_WEIGHT) + ' ' + metab + '_' +
                exchange_compartments_direction_S[id][0] + ' <--> 1 ' +
                metab[0:-1] + id)
            reactions_formulasNames.append(
                str(PLANT_WEIGHT / STEM_WEIGHT) + ' ' + metab + '_' +
                organs[organs_id.index(
                    exchange_compartments_direction_S[id][0])] +
                ' <--> 1 ' + metab[0:-1] + exchange_compartments_names[
                    exchange_compartments_id.index(id)])

            reactions_formulas.append(
                '1 ' + metab[0:-1] + id + ' <--> ' +
                str(
                    STEM_WEIGHT / PLANT_WEIGHT
                ) +
                ' ' + metab + '_' + exchange_compartments_direction_S[
                    id][0])
            reactions_formulasNames.append(
                metab[0:-1] +
                exchange_compartments_names[
                    exchange_compartments_id.index(id)]) + ' <--> 1' + \
            str(STEM_WEIGHT / PLANT_WEIGHT) + ' ' + metab + '_' + \
            organs[
                organs_id.index(
                    exchange_compartments_direction_S[id][0])]

exchangereactions_id = []
exchangereactions_formulas = []
exchangereactions_formulasNames = []
exchangestoichMat_rows = []
exchangestoichMat_columns = []
exchangestoichMat_values = []

# reactions for boundary metabolites

boundary_metabolites = [i for i, e in
                        enumerate(
                            metabolites_boundaryCondition) if e == 1]
nb_exchangereac = len(boundary_metabolites)

for index, metab in enumerate(boundary_metabolites):
    exchangestoichMat_columns.append(nb_reac + index)
    exchangestoichMat_values.append(-1)
    exchangereactions_id.append('exch_' + metabolites_id[metab])
    exchangereactions_formulas.append(metabolites_id[metab] + ' <--> ')
    exchangereactions_formulasNames.append(
        metabolites_name[metab] + ' <--> ')
    exchangestoichMat_rows.append(metab)

nb_exchangereac = len(exchangereactions_id)

print('number of reactions in three compartment network: ', len(reactions_id))
print('number of metabolites in three compartment network: ', len(metabolites_id))

# saved the parsed network file

save_data(reactions_name, "output/parser/reactions_name")
save_data(reactions_reversible, "output/parser/reactions_reversible")
save_data(boundary_metabolites, "output/parser/boundary_metabolites")
save_data(stoichMat_values, "output/parser/coefficients")
save_data(stoichMat_rows, "output/parser/rows")
save_data(stoichMat_columns, "output/parser/columns")
save_data(exchangestoichMat_columns,
          "output/parser/exchangestoichMat_columns")
save_data(exchangestoichMat_rows, "output/parser/exchangestoichMat_rows")
save_data(exchangereactions_id, "output/parser/exchangereactions_id")
save_data(metabolites_id, "output/parser/metabolites_id")
save_data(exchangestoichMat_values, "output/parser/exchangestoichMat_values")
save_data(reactions_formulas, "output/parser/reactions_formulas")
save_data(reactions_formulasNames, "output/parser/reactions_formulasNames")
save_data(reactions_id, "output/parser/reactions_id")
save_data(exchangereactions_formulasNames,
          "output/parser/exchangereactions_formulasNames")
save_data(exchangereactions_formulas,
          "output/parser/exchangereactions_formulas")
save_data(list_boundary, 'output/parser/boundary')
save_data(react_id_1comp,
          "output/parser/react_id_1comp")
save_data(metab_id_1comp,
          "output/parser/metab_id_1comp")
