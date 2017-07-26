#ifndef PARAMETERS_H
#define PARAMETERS_H

void get_expermient_vs_N_parameters(int N, bool WITH_COMB, double &Bcst_loc,
                                    double &Bcst_est, double &Comb_cst,
                                    int &loc_loops, int &est_loops,
                                    int &threshold_loops, int &comb_loops,
                                    double &tolerance_loc,
                                    double &tolerance_est);

void get_expermient_vs_K_parameters(int K, bool WITH_COMB, double &Bcst_loc,
                                    double &Bcst_est, double &Comb_cst,
                                    int &loc_loops, int &est_loops,
                                    int &threshold_loops, int &comb_loops,
                                    double &tolerance_loc,
                                    double &tolerance_est);

#endif
