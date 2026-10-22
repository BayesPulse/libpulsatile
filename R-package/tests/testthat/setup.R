# Test fits use short chains (burnin far below pv_adjust_max_iter), so silence
# the burn-in/adaptation warning globally; tests of the warning re-enable it.
.bayespulse_old_opts <- options(bayespulse.warn_burnin_adaptation = FALSE)
