# slyc-metabolic-model / input

This repository contains the input files for simulations on VYTOP (Virtual Young TOmato Plant).

Sl2186.xml: cell metabolic model of tomato plant, with biomass equation for leaf, stem, root

calibration.csv: file to calibrate different parameters of the model
-nitrogen ratio
-photosynthesis contribution (leaf/stem)
-organ relative weights
-constraining xylem fluxes (if xylem is 1) or not (if xylem is 0)
-transport cost: mol of atp required per mol exchanged
the user can modify these files depending on his/her biological question

constraints.csv: specification of upper and/or lower bound of fluxes when it is required (e.g no photon uptake in root)

exchanges.csv: specification of upper and/or lower bound of uptake fluxes for boundary metabolites

objective.csv: definition of the fluxes to minimize

variable.csv: list of reaction IDs for which we want to perform FVA (Flux Variability Analysis)

xylem.csv: lower and upper bound definition for xylem fluxes according to experimental data
they constrain the model if xylem value is 1 in calibration.csv




