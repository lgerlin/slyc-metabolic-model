# slyc-metabolic-model / input

This repository contains the input files for simulations on VYTOP (Virtual Young TOmato Plant).


Sl2184.xml: cell metabolic model of tomato plant, with biomass equation for leaf, stem, root


calibration.csv: file to calibrate different parameters of the model
-nh4_to_no3_ratio: ammonium to nitrate ratio constraint written as NO3_root_uptake = value * NH4_root_uptake. Optional (na if not desired)

-photosynthesis: maximum photosynthesis contribution of the stem compared to the leaf. Constraint written as value * photons_stems <= photons_leaf

-leaf_relative_weight: leaf relative weight compared to the root. Leaf_weight = value * root_weight. Used for exchange reactions between leaf and xylem and phloem

-stem_relative_weight: stem relative weight compared to the root. Stem_weight = value * root_weight. Used for exchange reactions between stem and xylem and phloem

-root_relative_weight: root relative weight compared to the root (supposed to be 1). Used for exchange reactions between root and xylem and phloem

-xylem_constraint: option to constraint xylem fluxes to experimental data measured (if xylem is 1) or not (if xylem is 0). Xylem experimental values are in xylem.csv file

-cost: transport cost in terms of mol of atp required per mol exchanged between organ and xylem/phloem

- ATP_maintenance_leaf: ATP maintenance term (Non Growth Associated Maintenance NGAM) for leaf. Has to be in chosen unit (mM or uM).h-1.gBleaf-1

- ATP_maintenance_stem: ATP maintenance term (Non Growth Associated Maintenance NGAM) for stem. Has to be in chosen unit (mM or uM).h-1.gBleaf-1

- ATP_maintenance_root: ATP maintenance term (Non Growth Associated Maintenance NGAM) for root. Has to be in chosen unit (mM or uM).h-1.gBleaf-1

- Leaf_biomass_growth: biomass growth of leaf.

- Rubisco_carboxylase_oxygenase_ratio: ratio of carboxylase to oxygenase activity of the Rubisco, same constraint used for stem and leaf. Constraints written as R_RBPCh_l = value * R_RBCh_1_l and R_RBPCh_s = value * R_RBCh_1_s

-unit: unit in which similations should be performed (uM or mM).

FVAreactions.csv: list of reaction IDs for which we want to perform FVA (Flux Variability Analysis)

xylem.csv: lower and upper bound definition for xylem fluxes according to experimental data
they constrain the model if xylem value is 1 in calibration.csv. Has to be in in chosen unit (mM or uM).h-1.gBroot-1 

reaction_inactivation.csv: list of reaction IDs for which we want to perform knock out or flux reduction
- reaction column: the reaction id on which you want to simulate mutant behavior
- status: the effect you want to impose on the associated reaction
'off' will assign the flux to zero
'red' will reduce the flux to between 50 percent of its value and zero
other status assigned will not affect the reaction flux