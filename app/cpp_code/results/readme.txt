export PYTHONPATH="/mnt/e/1.PhD/factor_graphs/g2o/ideas/unconstrained/mpc/mpc_space/../lib:/mnt/e/1.PhD/factor_graphs/g2o/ideas/unconstrained/mpc/mpc_space:$PYTHONPATH"

python3 -m results.results_harvesting 260217_1357_AMPL_0_3.bin
python3 -m results.results_harvesting_both 260220_1610_AMPL_0_G2O_ISPD_gn_GN_13_3.bin
