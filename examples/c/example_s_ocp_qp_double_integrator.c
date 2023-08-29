/**************************************************************************************************
 *                                                                                                 *
 * This file is part of HPIPM. *
 *                                                                                                 *
 * HPIPM -- High-Performance Interior Point Method. * Copyright (C) 2019 by
 *Gianluca Frison.                                                          *
 * Developed at IMTEK (University of Freiburg) under the supervision of Moritz
 *Diehl.              * All rights reserved. *
 *                                                                                                 *
 * The 2-Clause BSD License *
 *                                                                                                 *
 * Redistribution and use in source and binary forms, with or without *
 * modification, are permitted provided that the following conditions are met: *
 *                                                                                                 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *this                  * list of conditions and the following disclaimer. *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 ** this list of conditions and the following disclaimer in the documentation *
 *    and/or other materials provided with the distribution. *
 *                                                                                                 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND                 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *LIMITED TO, THE IMPLIED                   * WARRANTIES OF MERCHANTABILITY AND
 *FITNESS FOR A PARTICULAR PURPOSE ARE                          * DISCLAIMED. IN
 *NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR * ANY DIRECT,
 *INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 ** LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ** ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 ** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. *
 *                                                                                                 *
 * Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de *
 *                                                                                                 *
 **************************************************************************************************/

/**
 *
 * Basic example of solving an OCP QP with data provided in the format generated
 *by the routine d_ocp_qp_codegen(...);
 *
 **/

#include <stdio.h>
#include <stdlib.h>

#include <blasfeo_s_aux_ext_dep.h>

#include <hpipm_s_ocp_qp.h>
#include <hpipm_s_ocp_qp_dim.h>
#include <hpipm_s_ocp_qp_ipm.h>
#include <hpipm_s_ocp_qp_sol.h>
#include <hpipm_s_ocp_qp_utils.h>
#include <hpipm_timing.h>

// qp data as global data
extern int N;
extern int *nx;
extern int *nu;
extern int *nbu;
extern int *nbx;
extern int *ng;
extern int *nsbx;
extern int *nsbu;
extern int *nsg;
extern int *nbue;
extern int *nbxe;
extern int *nge;
extern float **hA;
extern float **hB;
extern float **hb;
extern float **hQ;
extern float **hR;
extern float **hS;
extern float **hq;
extern float **hr;
extern int **hidxbx;
extern float **hlbx;
extern float **hubx;
extern int **hidxbu;
extern float **hlbu;
extern float **hubu;
extern float **hC;
extern float **hD;
extern float **hlg;
extern float **hug;
extern float **hZl;
extern float **hZu;
extern float **hzl;
extern float **hzu;
extern int **hidxs;
extern float **hlls;
extern float **hlus;
extern int **hidxe;
// arg
extern int mode;
extern int iter_max;
extern float alpha_min;
extern float mu0;
extern float tol_stat;
extern float tol_eq;
extern float tol_ineq;
extern float tol_comp;
extern float reg_prim;
extern int warm_start;
extern int pred_corr;
extern int ric_alg;
extern int split_step;

// main
int main() {

  int ii, jj;

  int hpipm_status;

  int rep, nrep = 1000;

  hpipm_timer timer;

  /************************************************
   * ocp qp dim
   ************************************************/

  hpipm_size_t dim_size = s_ocp_qp_dim_memsize(N);
  void *dim_mem = malloc(dim_size);

  struct s_ocp_qp_dim dim;
  s_ocp_qp_dim_create(N, &dim, dim_mem);

  s_ocp_qp_dim_set_all(nx, nu, nbx, nbu, ng, nsbx, nsbu, nsg, &dim);

  //	d_ocp_qp_dim_codegen("examples/c/data/test_d_ocp_data.c", "w", &dim);

  /************************************************
   * ocp qp
   ************************************************/

  hpipm_size_t qp_size = s_ocp_qp_memsize(&dim);
  void *qp_mem = malloc(qp_size);

  struct s_ocp_qp qp;
  s_ocp_qp_create(&dim, &qp, qp_mem);

  s_ocp_qp_set_all(hA, hB, hb, hQ, hS, hR, hq, hr, hidxbx, hlbx, hubx, hidxbu,
                   hlbu, hubu, hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hidxs,
                   hlls, hlus, &qp);

  //	d_ocp_qp_codegen("examples/c/data/test_d_ocp_data.c", "a", &dim, &qp);

  /************************************************
   * ocp qp sol
   ************************************************/

  hpipm_size_t qp_sol_size = s_ocp_qp_sol_memsize(&dim);
  void *qp_sol_mem = malloc(qp_sol_size);

  struct s_ocp_qp_sol qp_sol;
  s_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

  /************************************************
   * ipm arg
   ************************************************/

  hpipm_size_t ipm_arg_size = s_ocp_qp_ipm_arg_memsize(&dim);
  void *ipm_arg_mem = malloc(ipm_arg_size);

  struct s_ocp_qp_ipm_arg arg;
  s_ocp_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);

  s_ocp_qp_ipm_arg_set_default(mode, &arg);

  s_ocp_qp_ipm_arg_set_mu0(&mu0, &arg);
  s_ocp_qp_ipm_arg_set_iter_max(&iter_max, &arg);
  s_ocp_qp_ipm_arg_set_alpha_min(&alpha_min, &arg);
  s_ocp_qp_ipm_arg_set_tol_stat(&tol_stat, &arg);
  s_ocp_qp_ipm_arg_set_tol_eq(&tol_eq, &arg);
  s_ocp_qp_ipm_arg_set_tol_ineq(&tol_ineq, &arg);
  s_ocp_qp_ipm_arg_set_tol_comp(&tol_comp, &arg);
  s_ocp_qp_ipm_arg_set_reg_prim(&reg_prim, &arg);
  s_ocp_qp_ipm_arg_set_warm_start(&warm_start, &arg);
  s_ocp_qp_ipm_arg_set_pred_corr(&pred_corr, &arg);
  s_ocp_qp_ipm_arg_set_ric_alg(&ric_alg, &arg);
  s_ocp_qp_ipm_arg_set_split_step(&split_step, &arg);

  //	d_ocp_qp_ipm_arg_codegen("examples/c/data/test_d_ocp_data.c", "a", &dim,
  //&arg);

  /************************************************
   * ipm workspace
   ************************************************/

  hpipm_size_t ipm_size = s_ocp_qp_ipm_ws_memsize(&dim, &arg);
  void *ipm_mem = malloc(ipm_size);

  struct s_ocp_qp_ipm_ws workspace;
  s_ocp_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

  /************************************************
   * ipm solver
   ************************************************/

  hpipm_tic(&timer);

  for (rep = 0; rep < nrep; rep++) {
    // call solver
    s_ocp_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);
    s_ocp_qp_ipm_get_status(&workspace, &hpipm_status);
  }

  double time_ipm = hpipm_toc(&timer) / nrep;

  /************************************************
   * print solution info
   ************************************************/

  printf("\nHPIPM returned with flag %i.\n", hpipm_status);
  if (hpipm_status == 0) {
    printf("\n -> QP solved!\n");
  } else if (hpipm_status == 1) {
    printf("\n -> Solver failed! Maximum number of iterations reached\n");
  } else if (hpipm_status == 2) {
    printf("\n -> Solver failed! Minimum step lenght reached\n");
  } else if (hpipm_status == 3) {
    printf("\n -> Solver failed! NaN in computations\n");
  } else {
    printf("\n -> Solver failed! Unknown return flag\n");
  }
  printf("\nAverage solution time over %i runs: %e [s]\n", nrep, time_ipm);
  printf("\n\n");

  /************************************************
   * extract and print solution
   ************************************************/

  // u

  int nu_max = nu[0];
  for (ii = 1; ii <= N; ii++)
    if (nu[ii] > nu_max)
      nu_max = nu[ii];

  float *u = malloc(nu_max * sizeof(float));

  printf("\nu = \n");
  for (ii = 0; ii <= N; ii++) {
    s_ocp_qp_sol_get_u(ii, &qp_sol, u);
    s_print_mat(1, nu[ii], u, 1);
  }

  // x

  int nx_max = nx[0];
  for (ii = 1; ii <= N; ii++)
    if (nx[ii] > nx_max)
      nx_max = nx[ii];

  float *x = malloc(nx_max * sizeof(float));

  printf("\nx = \n");
  for (ii = 0; ii <= N; ii++) {
    s_ocp_qp_sol_get_x(ii, &qp_sol, x);
    s_print_mat(1, nx[ii], x, 1);
  }

  // pi

  float *pi = malloc(nx_max * sizeof(float));

  printf("\npi = \n");
  for (ii = 0; ii < N; ii++) {
    s_ocp_qp_sol_get_pi(ii, &qp_sol, pi);
    s_print_mat(1, nx[ii + 1], pi, 1);
  }

  // all solution components at once

  float **u1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(u1 + ii, nu[ii], 1);
  float **x1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(x1 + ii, nx[ii], 1);
  float **ls1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(ls1 + ii, nsbu[ii] + nsbx[ii] + nsg[ii], 1);
  float **us1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(us1 + ii, nsbu[ii] + nsbx[ii] + nsg[ii], 1);
  float **pi1 = malloc((N) * sizeof(float *));
  for (ii = 0; ii < N; ii++)
    s_zeros(pi1 + ii, nx[ii + 1], 1);
  float **lam_lb1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(lam_lb1 + ii, nbu[ii] + nbx[ii], 1);
  float **lam_ub1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(lam_ub1 + ii, nbu[ii] + nbx[ii], 1);
  float **lam_lg1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(lam_lg1 + ii, ng[ii], 1);
  float **lam_ug1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(lam_ug1 + ii, ng[ii], 1);
  float **lam_ls1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(lam_ls1 + ii, nsbu[ii] + nsbx[ii] + nsg[ii], 1);
  float **lam_us1 = malloc((N + 1) * sizeof(float *));
  for (ii = 0; ii <= N; ii++)
    s_zeros(lam_us1 + ii, nsbu[ii] + nsbx[ii] + nsg[ii], 1);

  s_ocp_qp_sol_get_all(&qp_sol, u1, x1, ls1, us1, pi1, lam_lb1, lam_ub1,
                       lam_lg1, lam_ug1, lam_ls1, lam_us1);

  /************************************************
   * print ipm statistics
   ************************************************/

  int iter;
  s_ocp_qp_ipm_get_iter(&workspace, &iter);
  float res_stat;
  s_ocp_qp_ipm_get_max_res_stat(&workspace, &res_stat);
  float res_eq;
  s_ocp_qp_ipm_get_max_res_eq(&workspace, &res_eq);
  float res_ineq;
  s_ocp_qp_ipm_get_max_res_ineq(&workspace, &res_ineq);
  float res_comp;
  s_ocp_qp_ipm_get_max_res_comp(&workspace, &res_comp);
  float *stat;
  s_ocp_qp_ipm_get_stat(&workspace, &stat);
  int stat_m;
  s_ocp_qp_ipm_get_stat_m(&workspace, &stat_m);

  printf("\nipm return = %d\n", hpipm_status);
  printf(
      "\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n",
      res_stat, res_eq, res_ineq, res_comp);

  printf("\nipm iter = %d\n", iter);
  printf(
      "\nalpha_aff\tmu_aff\t\tsigma\t\talpha_prim\talpha_dual\tmu\t\tres_"
      "stat\tres_eq\t\tres_ineq\tres_comp\tobj\t\tlq fact\t\titref pred\titref "
      "corr\tlin res stat\tlin res eq\tlin res ineq\tlin res comp\n");
  s_print_exp_tran_mat(stat_m, iter + 1, stat, stat_m);

  printf("\nocp ipm time = %e [s]\n\n", time_ipm);

  /************************************************
   * get riccati matrices and vectors
   ************************************************/

#if 1
  printf("\nget Riccati recursion matrices and vectors\n");

  float *Lr0 = malloc(nu[0] * nu[0] * sizeof(float));
  float *Ls0 = malloc(nx[0] * nu[0] * sizeof(float));
  float *P0 = malloc(nx[0] * nx[0] * sizeof(float));
  float *lr0 = malloc(nu[0] * sizeof(float));
  float *p0 = malloc(nx[0] * sizeof(float));
  float *K0 = malloc(nu[0] * nx[0] * sizeof(float));
  float *k0 = malloc(nu[0] * sizeof(float));

  float *Lr1 = malloc(nu[1] * nu[1] * sizeof(float));
  float *Ls1 = malloc(nx[1] * nu[1] * sizeof(float));
  float *P1 = malloc(nx[1] * nx[1] * sizeof(float));
  float *lr1 = malloc(nu[1] * sizeof(float));
  float *p1 = malloc(nx[1] * sizeof(float));
  float *K1 = malloc(nu[1] * nx[1] * sizeof(float));
  float *k1 = malloc(nu[1] * sizeof(float));

  //
  s_ocp_qp_ipm_get_ric_Lr(&qp, &arg, &workspace, 0, Lr0);
  printf("\nLr0\n");
  s_print_exp_mat(nu[0], nu[0], Lr0, nu[0]);
  //
  s_ocp_qp_ipm_get_ric_Ls(&qp, &arg, &workspace, 0, Ls0);
  printf("\nLs0\n");
  s_print_exp_mat(nx[0], nu[0], Ls0, nx[0]);
  //
  s_ocp_qp_ipm_get_ric_P(&qp, &arg, &workspace, 0, P0);
  printf("\nP0\n");
  s_print_exp_mat(nx[0], nx[0], P0, nx[0]);
  //
  s_ocp_qp_ipm_get_ric_lr(&qp, &arg, &workspace, 0, lr0);
  printf("\nlr0\n");
  s_print_exp_mat(1, nu[0], lr0, 1);
  //
  s_ocp_qp_ipm_get_ric_p(&qp, &arg, &workspace, 0, p0);
  printf("\np0\n");
  s_print_exp_mat(1, nx[0], p0, 1);
  //
  s_ocp_qp_ipm_get_ric_K(&qp, &arg, &workspace, 0, K0);
  printf("\nK0\n");
  s_print_exp_mat(nu[0], nx[0], K0, nu[0]);
  //
  s_ocp_qp_ipm_get_ric_k(&qp, &arg, &workspace, 0, k0);
  printf("\nk0\n");
  s_print_exp_mat(1, nu[0], k0, 1);

  //
  s_ocp_qp_ipm_get_ric_Lr(&qp, &arg, &workspace, 1, Lr1);
  printf("\nLr1\n");
  s_print_exp_mat(nu[1], nu[1], Lr1, nu[1]);
  //
  s_ocp_qp_ipm_get_ric_Ls(&qp, &arg, &workspace, 1, Ls1);
  printf("\nLs1\n");
  s_print_exp_mat(nx[1], nu[1], Ls1, nx[1]);
  //
  s_ocp_qp_ipm_get_ric_P(&qp, &arg, &workspace, 1, P1);
  printf("\nP1\n");
  s_print_exp_mat(nx[1], nx[1], P1, nx[1]);
  //
  s_ocp_qp_ipm_get_ric_lr(&qp, &arg, &workspace, 1, lr1);
  printf("\nlr1\n");
  s_print_exp_mat(1, nu[1], lr1, 1);
  //
  s_ocp_qp_ipm_get_ric_p(&qp, &arg, &workspace, 1, p1);
  printf("\np1\n");
  s_print_exp_mat(1, nx[1], p1, 1);
  //
  s_ocp_qp_ipm_get_ric_K(&qp, &arg, &workspace, 1, K1);
  printf("\nK1\n");
  s_print_exp_mat(nu[1], nx[1], K1, nu[1]);
  //
  s_ocp_qp_ipm_get_ric_k(&qp, &arg, &workspace, 1, k1);
  printf("\nk1\n");
  s_print_exp_mat(1, nu[1], k1, 1);

  free(Lr0);
  free(Ls0);
  free(P0);
  free(lr0);
  free(p0);
  free(K0);
  free(k0);

  free(Lr1);
  free(Ls1);
  free(P1);
  free(lr1);
  free(p1);
  free(K1);
  free(k1);
#endif

  /************************************************
   * free memory and return
   ************************************************/

  free(dim_mem);
  free(qp_mem);
  free(qp_sol_mem);
  free(ipm_arg_mem);
  free(ipm_mem);

  free(u);
  free(x);
  free(pi);

  for (ii = 0; ii < N; ii++) {
    s_free(u1[ii]);
    s_free(x1[ii]);
    s_free(ls1[ii]);
    s_free(us1[ii]);
    s_free(pi1[ii]);
    s_free(lam_lb1[ii]);
    s_free(lam_ub1[ii]);
    s_free(lam_lg1[ii]);
    s_free(lam_ug1[ii]);
    s_free(lam_ls1[ii]);
    s_free(lam_us1[ii]);
  }
  s_free(u1[ii]);
  s_free(x1[ii]);
  s_free(ls1[ii]);
  s_free(us1[ii]);
  s_free(lam_lb1[ii]);
  s_free(lam_ub1[ii]);
  s_free(lam_lg1[ii]);
  s_free(lam_ug1[ii]);
  s_free(lam_ls1[ii]);
  s_free(lam_us1[ii]);

  free(u1);
  free(x1);
  free(ls1);
  free(us1);
  free(lam_lb1);
  free(lam_ub1);
  free(lam_lg1);
  free(lam_ug1);
  free(lam_ls1);
  free(lam_us1);

  return 0;
}
