"""
Specification of the Ros-dependent Wnt/beta-Catenin pathway.

Staehlke, S., Haack, F., Waldner, A.C., Koczan, D., Moerke, C., Mueller, P.,
Uhrmacher, A.M., Nebe, J.B.: Ros dependent wnt/β-catenin pathway and its reg-
ulation on defined micro-pillars—a combined in vitro and in silico study. Cells 9(8)
(2020). https://doi.org/10.3390/cells9081784

Original author: Fiete Haack
Adapted for use with evolib by: Justin Kreikemeyer
"""

from evolib.reaction import Reaction
from evolib.reaction_library import ReactionLibrary

Axin_u = 0
Axin_p = 1
Bcat_a = 2
Bcat_i = 3
Dvl_i = 4
Dvl_a = 5
Nrx_nO = 6
Nrx_O = 7
DvlNrx = 8
#Ros_i = 9  # not in the original (historic reasons)
Ros_a = 9
DvlAxin_u = 10
DvlAxin_p = 11
ICAT = 12

Nuc_Bcat_a = 13 # "free and *a*ctive"
Nuc_Bcat_i = 14 # BcatIcat ("*i*nactive")
Nuc_Bcat_c = 15 # BcatTCF ("*c*omplex")
Nuc_TCF = 16
#Nuc_Sox17 = 17
Nuc_ICAT = 17
Counter = 18

species_names = [key for key, val in locals().items() if isinstance(val, int)]
num_species = len(species_names)
# separate species by cytosol and nucleus
species_subgroups = [[0,1,2,3,4,5,6,7,8,9,10,11,12,18], [13,14,15,16,17]]


# ** beta-catenin signalling **
nbetacyt = 12989
nbetanuc = 5282
nAxin = 252
nAxinP = 219
# ** Ros-Dvl signalling **
nRos = 0
nDvl = 855
nNrx = 18
nDvlNrx = 166200
# ** regulatory factors **
nICAT = 1200
#nSox = 100
nTCF = 7714

# from SESSL experiment
nRos = 0
nICAT = 1200
#nSox = 100
nTCF =  7714

_init = {
  Counter: 0,
  Dvl_i: nDvl,
  Dvl_a: 0,
  Nrx_nO: nNrx,
  Nrx_O: 0,
  DvlNrx: nDvlNrx,
  DvlAxin_u: 0,
  DvlAxin_p: 0,
  #Ros_i: nRos,
  Ros_a: 0,
  Bcat_a: nbetacyt,
  Bcat_i: 0,
  ICAT: 0,
  Axin_u: nAxin,
  Axin_p: nAxinP,
  # Nuc
  Nuc_Bcat_a: nbetanuc,
  Nuc_Bcat_i: 0,
  Nuc_Bcat_c: 0,
  Nuc_ICAT: nICAT,
  #Nuc_Sox17: nSox,
  Nuc_TCF: nTCF,
}
wnt_init = list(x[1] for x in sorted(_init.items(), key=lambda x: x[0]))

kApA = 0.03; # k7r Basal dephosphorylation of AxinP
kAAp = 0.03; # k7 Phosphorylation of Axin
kApdeg = 4.48E-3; # k6 Degradation of phosphorylated Axin
kAdeg = 4.48E-3; # k6 Degradation of unphosphorylated Axin
kAsyn = 4E-4; # k16 Axin synthesis (beta-catenin mediated)

kbetasyn = 600.0; # k9 beta-catenin synthesis
kbetadeg_act = 2.1E-4; # k8 axin-induced degradation of beta-catenin
kbetadeg = 1.13E-4; # k9r basal degradation of beta-catenin
kbetain = 0.0549; # k10 beta-catenin shuttling into nucleus
kbetaout = 0.135; # k10r beta-catenin shuttling out of nucleus

# ** Ros-Dvl Signalling **
kRosSyn = 100.0; # k1 ROS synthesis

# Nrx
kNrxRos = 5E2; # k2_h Oxidation of Nrx by Ros
kNrxNo = 2E-2; # k3r_h Reduction of Nrx

# Dvl
kDvlSponAgg = 5E-04; # k4 spontaneous aggregation of Dvl
kDvldisAgg = 0.5; # k4r basal disscoiation of Dvl aggregates

# Dvl-Nrx
kDvlNrxBind = 22.5; # k3r Binding of Nrx and Dvl
kDvlNrxUnbind = 2.3E-2; # k3 Basal unbinding of Nrx and Dvl
kDvlNrxRos = 3.2E2; # k2 Unbinding of Nrx and Dvl complex by ROS 

# Dvl-Axin
kDvlAxinBind = 0.075; # k5 Binding of (activated) DVL and Axin
kDvlAxinUnbind = 6.8E-2; # k5r Unbinding of Dvl/Axin complex

# regulatory factors
kICATsyn = 250.0; # k_h2 Increasing ICAT concentration
ka_IcatBcat = 0.1; # k11 Binding of ICAT and beta-catenin
kd_IcatBcat = 0.032; # k11r Unbinding of ICAT and beta-catenin complex
ka_TcfBcat = 0.00196; # k13 Binding of TCF and beta-catenin
kd_TcfBcat = 0.0141; # k13r Unbinding of TCF and beta-catenin
kTcfSyn = 0.029; # k15 Synthesis of TCF
kSox = 2.3e-4; # k14 degradation of TCF/beta-catenin complex

# from SESSL experiment
kTcfSyn = 100.0
kICATsyn = 250.0
kRosSyn = 200.0

Reac = lambda x, y, r: Reaction(x, y, r, num_species, species_names=species_names)
wnt_model = ReactionLibrary([

  Reac([], [Counter], 1.0),

  # **** Ros-Dvl Signalling ****

  # (R1) Ros Synthesis
  Reac([], [Ros_a], kRosSyn),
  Reac([Counter], [Counter, Ros_a], 0.2),
  # (R2) Forced Unbinding of Dvl from Nrx by Ros
  Reac([DvlNrx, Ros_a], [Dvl_i, Nrx_O], kDvlNrxRos),
  # (R2_h) Oxidation of Nrx by Ros
  Reac([Nrx_nO, Ros_a], [Nrx_O], kNrxRos),
  # (R3) Basal unbinding of Dvl from Nrx
  Reac([DvlNrx], [Dvl_i, Nrx_nO], kDvlNrxUnbind),
  # (R3r) Binding of Dvl by Nrx
  Reac([Dvl_i, Nrx_nO], [DvlNrx], kDvlNrxBind),
  # (R3r_h) Reduction of Nrx
  Reac([Nrx_O], [Nrx_nO], kNrxNo),
  # (R3r/R4) Forced Disaggregation of Dvl by un-oxidized Nrx
  Reac([Dvl_a, Nrx_nO], [DvlNrx], kDvlNrxBind),
  # (R4) Activation (by e.g. aggregation) of Dvl
  Reac([Dvl_i], [Dvl_a], kDvlSponAgg),
  # (R4r) Dynamic deactivation (e.g. by disaggregation) of Dvl
  Reac([Dvl_a], [Dvl_i], kDvldisAgg),

  # **** Axin Dvl signalling ****

  # (R5) Axin binding by activated Dvl
  Reac([Dvl_a, Axin_p], [DvlAxin_p], kDvlAxinBind),
  Reac([Dvl_a, Axin_u], [DvlAxin_u], kDvlAxinBind),
  # (R5r) Axin Dvl unbinding
  Reac([DvlAxin_p], [Dvl_a, Axin_p], kDvlAxinUnbind),
  Reac([DvlAxin_u], [Dvl_a, Axin_u], kDvlAxinUnbind),
  # (R6) Axin degradation
  Reac([Axin_p], [], kApdeg),
  Reac([Axin_u], [], kApdeg),
  # (R7r) Basal AxinP dephosphorylation
  Reac([Axin_p], [Axin_u], kApA),
  # (R7) Axin phosphorylation
  Reac([Axin_u], [Axin_p], kAAp),

  # **** Beta-catenin signalling ****

  # (R8) Activated beta-catenin degradation
  Reac([Axin_p, Bcat_a], [Axin_p], kbetadeg_act),
  Reac([Axin_p, Bcat_i], [Axin_p], kbetadeg_act),
  #Reac([Axin_p, Bcat_c], [Axin_p], kbetadeg_act), # (no Bcat_c can be in cell, only Nuc)
  # (R9) Beta-catenin synthesis
  Reac([], [Bcat_a], kbetasyn),
  # (R9r) Basal beta-catenin degradation
  Reac([Bcat_a], [], kbetadeg),
  Reac([Bcat_i], [], kbetadeg),
  #Reac([Bcat_c], [], kbetadeg),
  Reac([Nuc_Bcat_a], [], kbetadeg),
  Reac([Nuc_Bcat_i], [], kbetadeg),
  Reac([Nuc_Bcat_c], [], kbetadeg),
  # (R10) Beta-catenin shuttling into the nucleus
  Reac([Bcat_a], [Nuc_Bcat_a], kbetain),
  Reac([Bcat_i], [Nuc_Bcat_i], kbetain),
  # (R10r) Beta-catenin shuttling out of the nucleus
  Reac([Nuc_Bcat_a], [Bcat_a], kbetaout),
  Reac([Nuc_Bcat_i], [Bcat_i], kbetaout),

  # **** ICAT and SOX17 signaling ****

  # (H2) increasing ICAT concentration
  Reac([], [ICAT], kICATsyn),
  # (R11) ICAT binding beta-catenin
  Reac([ICAT, Bcat_a], [Bcat_i], ka_IcatBcat),
  Reac([Nuc_ICAT, Nuc_Bcat_a], [Nuc_Bcat_i], ka_IcatBcat),
  # (R11r) unbinding of ICAT and beta-catenin
  Reac([Bcat_i], [ICAT, Bcat_a], kd_IcatBcat),
  Reac([Nuc_Bcat_i], [Nuc_ICAT, Nuc_Bcat_a], kd_IcatBcat),
  # (R12) ICAT shuttling into the nucleus
  Reac([ICAT], [Nuc_ICAT], kbetain),
  # (R12r) ICAT shuttling out of the nucleus
  Reac([Nuc_ICAT], [ICAT], kbetaout),
  # (R13) TCF binding beta-catenin
  Reac([Nuc_Bcat_a, Nuc_TCF], [Nuc_Bcat_c], ka_TcfBcat),
  # (R13r) unbinding of TCF and beta-catenin
  Reac([Nuc_Bcat_c], [Nuc_Bcat_a, Nuc_TCF], kd_TcfBcat),
  # (R14) Degradation of TCF/beta-catenin complex
  #Reac([Nuc_Sox17, Nuc_Bcat_c], [Nuc_Sox17], kSox),
  Reac([Nuc_Bcat_c], [], kSox * 100.0), # Sox does not change from 100
  # (R15) TCF synthesis
  Reac([], [Nuc_TCF], kTcfSyn),
  # (R16) Axin synthesis
  Reac([Nuc_Bcat_c], [Nuc_Bcat_c, Axin_u], kAsyn),

], num_species, species_names=species_names)

fiete_fitted_reactions = [
  Reac([], [ICAT], kICATsyn),
  Reac([], [Ros_a], kRosSyn),
  Reac([ICAT, Bcat_a], [Bcat_i], ka_IcatBcat),
  Reac([Nuc_ICAT, Nuc_Bcat_a], [Nuc_Bcat_i], ka_IcatBcat),
  Reac([Bcat_i], [ICAT, Bcat_a], kd_IcatBcat),
  Reac([Nuc_Bcat_i], [Nuc_ICAT, Nuc_Bcat_a], kd_IcatBcat),
  Reac([], [Nuc_TCF], kTcfSyn),
  #Reac([Nuc_Sox17, Nuc_Bcat_c], [Nuc_Sox17], kSox),
  Reac([Nuc_Bcat_c], [], kSox * 100.0),
]
# fixate all reactions except those Fiete fitted
fixated_reactions = [reac for reac in wnt_model.reactions if reac not in fiete_fitted_reactions]

def wnt_penalty(member):
  """Penalty function for the genetic fitness to include some background knowledge."""
  penalty = 0
  # penalize use of Counter species
  for idx, r in enumerate(member.non_fixated_reactions):
    if 18 in r.reactands or 18 in r.products:
      penalty += 0.01
  shuttle_in = False
  shuttle_out = False
  # penalize absence of shuttling reactions in/out of nucleus with 1 (each)
  for r in member.reactions:
    if any(reac in species_subgroups[0] for reac in r.reactands) and any(prod in species_subgroups[1] for prod in r.products):
      shuttle_in = True
    if any(reac in species_subgroups[1] for reac in r.reactands) and any(prod in species_subgroups[0] for prod in r.products):
      shuttle_out = True
  penalty += 1 - int(shuttle_in) + 1 - int(shuttle_out)
  return penalty

if __name__ == "__main__":
  from evolib.integrate import gen_data_for, plot_sim_trace_of
  print(wnt_model, len(wnt_model.reactions))
  data = gen_data_for(wnt_model, wnt_init, t_end=1440, n_points=0)
  data.to_csv("test.csv", index=False)
  wnt_model.ref_data_path = "test.csv"
  plot_sim_trace_of(wnt_model, t_end=1440, n_points=0)



