##################################################
############# MOFA - Pareto ######################
##################################################

require(tibble)
require(reticulate)

devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"), force = TRUE)
reticulate::use_condaenv("~/miniconda3/envs/mofa/")
library(MOFA2)

path.MOFA = "/path/to/preprocessed/matrices/"

## Load matrices
load(paste0(path.MOFA,"D_expr_MOFA.RData"))
load(paste0(path.MOFA,"D_met.proB_MOFA.RData"))
load(paste0(path.MOFA,"D_met.bodB_MOFA.RData"))
load(paste0(path.MOFA,"D_met.enhB_MOFA.RData"))
load(paste0(path.MOFA,"D_cnv_MOFA.RData"))
load(file=paste0(path.MOFA,"D_loh_MOFA.RData"))
load(file=paste0(path.MOFA,"D_alt_MOFA.RData"))

## Prepare MOFA
# MOFA object creation
MOFAobject = create_mofa(list("RNA" = D_expr_MOFA,"MethPro" = D_met.proB_MOFA,"MethBod" = D_met.bodB_MOFA,"MethEnh" = D_met.enhB_MOFA, 
                              "Total" = D_cnv_MOFA, "Minor" = D_loh_MOFA, "Alt" = D_alt_MOFA)) 

data_opts = get_default_data_options(MOFAobject)
model_opts = get_default_model_options(MOFAobject)
model_opts$num_factors = 10
train_opts = get_default_training_options(MOFAobject)
train_opts$convergence_mode = "slow"
stochastic_opts = get_default_stochastic_options(MOFAobject)
model_opts$likelihoods["Alt"] = "bernoulli"
train_opts$maxiter = 10000

# Prepare
MOFAobject = prepare_mofa(object=MOFAobject, data_options=data_opts, model_options=model_opts, training_options=train_opts )

## Run MOFA
MOFAobject.trained = run_mofa(MOFAobject, save_data = T, outfile = "/path/to/MOFAobject.hdf5")

# Extract LFs coordinates
LFs = as.data.frame(get_factors(MOFAobject.trained)$group1)

## Arc proportions
BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo"))
reticulate::use_python("~/.local/share/r-miniconda/envs/reticulate_PCHA/bin/python")
library(ParetoTI)

path.Pareto = "/path/to/phenotypic/map/data/"

LFs = read.table("/path/to/supplementary/table2/LFs_73.Cordin.txt",header = TRUE, sep = "\t")
colnames(LFs)[1:4] = paste0(c("Ploidy", "Morphology", "Adaptive-response", "CIMP-index"), "_factor")
LFs = LFs[,c(2,3,4,1)] # RNA contribution order

arc_ks.LFs = lapply(2:ncol(LFs), function(k) k_fit_pch( t(LFs[,1:k]), ks = 2:6, check_installed = T,
                                                        bootstrap = T, bootstrap_N = 200, maxiter = 1000,
                                                        bootstrap_type = "m", seed = 2543, 
                                                        volume_ratio = "t_ratio", 
                                                        delta=0, conv_crit = 1e-04, order_type = "align",
                                                        sample_prop = 0.75))

arc_ks.LFs.noboot = lapply(2:ncol(LFs), function(k) k_fit_pch( t(LFs[,1:k]), ks = 2:6, check_installed = T,
                                                               bootstrap = F, 
                                                               volume_ratio = "t_ratio" ) )

tern_arc = as_tibble( t( arc_ks.LFs.noboot[[1]]$pch_fits$S[[2]] ) )
tern_arc$sample = rownames(t( arc_ks.LFs.noboot[[1]]$pch_fits$S[[2]] ))
write.table(tern_arc, paste0(path.Pareto,"tern_arc.txt"), row.names = T, col.names = T, quote = F, sep="\t")

arc_boot_best.boot = fit_pch_bootstrap(t(LFs[,1:2]), n = 200, sample_prop = 0.75, noc = 3, delta = 0, conv_crit = 1e-04, type = "s")
arc_boot_best = average_pch_fits(arc_boot_best.boot) 
arc = as.data.frame(arc_boot_best$XC)
colnames(arc) = c("Arc1","Arc2","Arc3")
write.table(arc, paste0(path.Pareto,"arc_pos.txt"),row.names = T,col.names = T,quote = F, sep="\t")
