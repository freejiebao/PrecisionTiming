#!/bin/bash

crab status -d crab_DY_noPU_TuneCUETP8M1
crab status -d crab_DY_noPU_0J
crab status -d crab_DY_noPU_1J
crab status -d crab_DY_noPU_2J
crab status -d crab_DY_noPU_3J
crab status -d crab_DY_PU200_TuneCUETP8M1_ext
crab status -d crab_DY_PU200_TuneCUETP8M1
crab status -d crab_DY_PU200_0J_ext
crab status -d crab_DY_PU200_0J
crab status -d crab_DY_PU200_1J_ext
crab status -d crab_DY_PU200_1J
crab status -d crab_DY_PU200_2J_ext
crab status -d crab_DY_PU200_2J
crab status -d crab_DY_PU200_3J
crab status -d crab_QCD_noPU
crab status -d crab_QCD_PU200_ext
crab status -d crab_QCD_PU200


crab status -d crab_DY_PU200_GEN-SIM-RECO_ext
crab status -d crab_QCD_PU200_GEN-SIM-RECO
crab status -d crab_QCD_noPU_GEN-SIM-RECO
crab status -d crab_DY_1J_PU200_GEN-SIM-RECO
crab status -d crab_DY_1J_PU200_GEN-SIM-RECO_ext
crab status -d crab_DY_1J_noPU_GEN-SIM-RECO
crab status -d crab_DY_2J_PU200_GEN-SIM-RECO
crab status -d crab_DY_2J_noPU_GEN-SIM-RECO
crab status -d crab_DY_3J_PU200_GEN-SIM-RECO

finished:

crab status -d crab_DY_noPU_GEN-SIM-RECO
crab status -d crab_DY_PU200_GEN-SIM-RECO
crab status -d crab_QCD_PU200_GEN-SIM-RECO_ext
crab status -d crab_DY_2J_PU200_GEN-SIM-RECO_ext
crab status -d crab_DY_3J_noPU_GEN-SIM-RECO

on tape
crab_DY_2J_noPU_GEN-SIM-RECO
crab_DY_3J_PU200_GEN-SIM-RECO
crab_DY_1J_noPU_GEN-SIM-RECO

delete crab_DY_3J_noPU_GEN-SIM-RECO
delete crab_DY_2J_PU200_GEN-SIM-RECO_ext
crab_QCD_PU200_GEN-SIM-RECO_ext
delete crab_DY_PU200_GEN-SIM-RECO
delete crab_DY_noPU_GEN-SIM-RECO
