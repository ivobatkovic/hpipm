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

#include <stdlib.h>

// QP size

// horizon lenght
int N = 5;
// number of input
static int nnu[6] = {1, 1, 1, 1, 1, 0};
// number of states
static int nnx[6] = {2, 2, 2, 2, 2, 2};
// number of input box constraints
static int nnbu[6] = {0, 0, 0, 0, 0, 0};
// number of states box constraints
static int nnbx[6] = {2, 0, 0, 0, 0, 0};
// number of general constraints
static int nng[6] = {0, 0, 0, 0, 0, 0};
// number of softed constraints on state box constraints
static int nnsbx[6] = {0, 0, 0, 0, 0, 0};
// number of softed constraints on input box constraints
static int nnsbu[6] = {0, 0, 0, 0, 0, 0};
// number of softed constraints on general constraints
static int nnsg[6] = {0, 0, 0, 0, 0, 0};
// number of input box constraints considered as equality
static int nnbue[6] = {0, 0, 0, 0, 0, 0};
// number of states box constraints considered as equality
static int nnbxe[6] = {0, 0, 0, 0, 0, 0};
// number of general constraints considered as equality
static int nnge[6] = {0, 0, 0, 0, 0, 0};

// QP data

// float integrator
const float ts = 1;
static float A[] = {1, 0, ts, 1};
//
static float B[] = {ts * ts / 2.0, ts};
//
static float b[] = {0, 0};

//
static float Q[] = {0, 0, 0, 1};
//
static float R[] = {.1};
//
static float S[] = {0, 0};
//
static float q[] = {0, 0};
//
static float r[] = {0};

//
static float lbx0[] = {1, 1};
//
static float ubx0[] = {1, 1};
//
static int idxbx0[] = {0, 1};

//
static float u_guess[] = {0};
//
static float x_guess[] = {0, 0};
//
static float sl_guess[] = {};
//
static float su_guess[] = {};

// array of pointers

//
static float *AA[5] = {A, A, A, A, A};
//
static float *BB[5] = {B, B, B, B, B};
//
static float *bb[5] = {b, b, b, b, b};
//
static float *QQ[6] = {Q, Q, Q, Q, Q, Q};
//
static float *RR[6] = {R, R, R, R, R, R};
//
static float *SS[6] = {S, S, S, S, S, S};
//
static float *qq[6] = {q, q, q, q, q, q};
//
static float *rr[6] = {r, r, r, r, r, r};
//
static int *iidxbx[6] = {idxbx0, NULL, NULL, NULL, NULL, NULL};
//
static float *llbx[6] = {lbx0, NULL, NULL, NULL, NULL, NULL};
//
static float *uubx[6] = {ubx0, NULL, NULL, NULL, NULL, NULL};
//
static int *iidxbu[6] = {};
//
static float *llbu[6] = {};
//
static float *uubu[6] = {};
//
static float *CC[6] = {};
//
static float *DD[6] = {};
//
static float *llg[6] = {};
//
static float *uug[6] = {};
//
static float *ZZl[6] = {};
//
static float *ZZu[6] = {};
//
static float *zzl[6] = {};
//
static float *zzu[6] = {};
//
static int *iidxs[6] = {};
//
static float *llls[6] = {};
//
static float *llus[6] = {};
//
static float *iidxe[6] = {};

//
static float *uu_guess[6] = {u_guess, u_guess, u_guess,
                             u_guess, u_guess, u_guess};
//
static float *xx_guess[6] = {x_guess, x_guess, x_guess,
                             x_guess, x_guess, x_guess};
//
static float *ssl_guess[6] = {sl_guess, sl_guess, sl_guess,
                              sl_guess, sl_guess, sl_guess};
//
static float *ssu_guess[6] = {su_guess, su_guess, su_guess,
                              su_guess, su_guess, su_guess};

// export as global data

int *nu = nnu;
int *nx = nnx;
int *nbu = nnbu;
int *nbx = nnbx;
int *ng = nng;
int *nsbx = nnsbx;
int *nsbu = nnsbu;
int *nsg = nnsg;
int *nbue = nnbue;
int *nbxe = nnbxe;
int *nge = nnge;

float **hA = AA;
float **hB = BB;
float **hb = bb;
float **hQ = QQ;
float **hR = RR;
float **hS = SS;
float **hq = qq;
float **hr = rr;
int **hidxbx = iidxbx;
float **hlbx = llbx;
float **hubx = uubx;
int **hidxbu = iidxbu;
float **hlbu = llbu;
float **hubu = uubu;
float **hC = CC;
float **hD = DD;
float **hlg = llg;
float **hug = uug;
float **hZl = ZZl;
float **hZu = ZZu;
float **hzl = zzl;
float **hzu = zzu;
int **hidxs = iidxs;
float **hlls = llls;
float **hlus = llus;
float **hidxe = iidxe;

float **hu_guess = uu_guess;
float **hx_guess = xx_guess;
float **hsl_guess = ssl_guess;
float **hsu_guess = ssu_guess;

// arg
int mode = 0;
int iter_max = 30;
float alpha_min = 1e-8;
float mu0 = 1e2;
float tol_stat = 1e-4;
float tol_eq = 1e-5;
float tol_ineq = 1e-5;
float tol_comp = 1e-5;
float reg_prim = 1e-12;
int warm_start = 0;
int pred_corr = 1;
int ric_alg = 0;
int split_step = 1;
