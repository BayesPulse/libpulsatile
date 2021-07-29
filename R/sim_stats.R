# #-------------------------------------------------------------------------------
# # Summary stats functions
# #
# #   Main functions: bias(es) and coverage(s)
# #     Results are by dataset, plural versions calculate for a provided vector of
# #     variables.
# #-------------------------------------------------------------------------------
# 
# #---------------------------------------
# # calculate mode
# #---------------------------------------
# Mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# } 
# 
# 
# #---------------------------------------
# # match_simparm_name()
# #   match provided variable name from common 
# #   mcmc chain with population parameter from 
# #   simulation (simparms.Rds)
# #---------------------------------------
# match_simparm_name <- 
#   function(mcmc.var.name) {
#     switch(mcmc.var.name,
#            num.pulses       = "ipimean",
#            baseline         = "constant.baseline",
#            half.life        = "constant.halflife",
#            mean.pulse.mass  = "mua",
#            mean.pulse.width = "muw",
#            model.error      = "vare",
#            var.mass         = "vara",
#            var.width        = "varw",
#            total.num.pulses = "ipimean")
#   }
# 
# #-----------------------------------------------------------
# # bias() Function
# #   Calculates bias and 95% credible interval of bias
# #   Depends on match_simparms_name() in helperfns.R
# #-----------------------------------------------------------
# bias <- 
#   function(var, sim.parms, common.chain) {
#     require(dplyr)
#     require(lazyeval)
# 
#     # get matching variable/parameter name for sim data
#     sim.var <- match_simparm_name(var)
# 
#     # get name of sim and mcmc case from common.chain
#     curr.sim.case  <- unique(common.chain$sim.case)
#     curr.mcmc.case <- unique(common.chain$mcmc.case)
# 
#     # grab true value and sim name
#     sim.parms %<>%
#       filter(case %in% curr.sim.case) %>%
#       select_("case", sim.var) %>%
#       mutate(mcmc.case = curr.mcmc.case) %>%
#       rename_("true.value" = sim.var)
#     true.value <- sim.parms$true.value
# 
#     # calculate bias for chains from each dataset
#     dataset.mean <- 
#       common.chain %>%
#         group_by(dataset) %>%
#         select_("dataset", var) %>%
#         summarise_(post.mean = interp(~mean(var), var = as.name(var))) 
# 
#     dataset.bias <- 
#       dataset.mean %>%
#         group_by(dataset) %>% # this doesn't look necessary
#         mutate(bias = post.mean - true.value, 
#                relative.bias = (post.mean - true.value)/true.value) %>%
#         ungroup 
# 
#     # Combine results with useful detail
#     results <-
#       cbind(dataset.bias, 
#             true.value, 
#             "sim.case" = curr.sim.case, 
#             "mcmc.case" = curr.mcmc.case) %>% tbl_df
# 
#     return(results)
# 
#   }
# 
# #-----------------------------------------------------------
# # biases() Function
# #   wrapper for bias for calculating all variables
# #-----------------------------------------------------------
# biases <- 
#   function(vars, sim.parms, common.chain) {
#     result <- 
#       lapply(vars, function(x) {
#           df1 <- bias(x, sim.parms, common.chain)
#           df1 %<>% mutate(variable = x)
#         }) %>% do.call(rbind, .)
#     return(result)
#   }
# 
# 
# 
# #-----------------------------------------------------------
# # coverage() Function
# #   Calculates 95% coverage of provided variable
# #   Depends on match_simparms_name() in helperfns.R
# #-----------------------------------------------------------
# coverage <- 
#   function(var, sim.parms, common.chain) {
#     require(dplyr)
#     require(lazyeval)
# 
#     # get matching variable/parameter name for sim data
#     sim.var <- match_simparm_name(var)
# 
#     # get name of sim and mcmc case from common.chain
#     curr.sim.case  <- unique(common.chain$sim.case)
#     curr.mcmc.case <- unique(common.chain$mcmc.case)
# 
#     # grab true value and sim name
#     sim.parms %<>%
#       filter(case %in% curr.sim.case) %>%
#       select_("case", sim.var) %>%
#       mutate(mcmc.case = curr.mcmc.case) %>%
#       rename_("true.value" = sim.var)
#     true.value <- sim.parms$true.value
# 
#     # calculate coverage for chains from each dataset
#     coverage.df <- 
#       common.chain %>% 
#         group_by(dataset) %>%
#           select_("dataset", var) %>%
#           summarise_(lower = interp(~quantile(var, 0.025), var = as.name(var)),
#                      upper = interp(~quantile(var, 0.975), var = as.name(var))) %>%
#           # this gives credible intervals for each dataset (may need later)
#           mutate(in.ci = ifelse(true.value >= lower & true.value <= upper, 1, 0)) %>%
#           ungroup()
# 
#     # Combine results with useful detail
#     results <-
#       cbind(coverage.df, 
#             true.value, 
#             "sim.case" = curr.sim.case, 
#             "mcmc.case" = curr.mcmc.case) %>% tbl_df
# 
#     return(results)
# 
#   }
# 
# #-----------------------------------------------------------
# # coverages() Function
# #   wrapper for coverage() for calculating all variables
# #-----------------------------------------------------------
# coverages <- 
#   function(vars, sim.parms, common.chain) {
#     result <- 
#       lapply(vars, function(x) {
#           df1 <- coverage(x, sim.parms, common.chain)
#           df1 %<>% mutate(variable = x)
#         }) %>% do.call(rbind, .)
#     return(result)
#   }
# 
