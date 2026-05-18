:: Run SIMMAP predictions. ----
:: DR.LH
Rscript r_script_05_DR.LH_predict_simmap_1_save_2_dep.r
TIMEOUT /T 10
Rscript r_script_05_DR.LH_predict_simmap_2_load_2_dep.r
TIMEOUT /T 10
:: DEP2.L
Rscript r_script_05_DEP2.L_predict_simmap_1_save_2_dep.r
TIMEOUT /T 10
Rscript r_script_05_DEP2.L_predict_simmap_2_load_2_dep.r
TIMEOUT /T 10
:: DEP2.M
Rscript r_script_05_DEP2.M_predict_simmap_1_save_2_dep.r
TIMEOUT /T 10
Rscript r_script_05_DEP2.M_predict_simmap_2_load_2_dep.r
TIMEOUT /T 10
:: DEP2.H
Rscript r_script_05_DEP2.H_predict_simmap_1_save_2_dep.r
TIMEOUT /T 10
Rscript r_script_05_DEP2.H_predict_simmap_2_load_2_dep.r
TIMEOUT /T 10
