struct rxSolveF {
  //
  // std::string estStr;
  // std::string gradStr;
  // std::string obfStr;
  //
  t_dydt dydt = NULL;
  t_calc_jac calc_jac = NULL;
  t_calc_lhs calc_lhs = NULL;
  t_update_inis update_inis = NULL;
  t_dydt_lsoda_dum dydt_lsoda_dum = NULL;
  t_dydt_liblsoda dydt_liblsoda = NULL;
  t_jdum_lsoda jdum_lsoda = NULL;
  t_set_solve set_solve = NULL;
  t_get_solve get_solve = NULL;
  int global_jt = 2;
  int global_mf = 22;
  int global_debug = 0;
  int neq = NA_INTEGER;
};

typedef int (*iniSubjectI_t)(int solveid, int inLhs, rx_solving_options_ind *ind, rx_solving_options *op, rx_solve *rx,
                             t_update_inis u_inis);

extern "C" {
  typedef void (*ind_solve_t)(rx_solve *rx, unsigned int cid, t_dydt_liblsoda dydt_lls,
                              t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
                              t_dydt c_dydt, t_update_inis u_inis, int jt);
  typedef rx_solve* (*getRxSolve_t)();
  typedef double (*getTime_t)(int idx, rx_solving_options_ind *ind);
  typedef int (*isRstudio_t)();

}
