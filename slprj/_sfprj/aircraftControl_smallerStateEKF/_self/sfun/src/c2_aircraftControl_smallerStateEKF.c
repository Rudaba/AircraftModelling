/* Include files */

#include <stddef.h>
#include "blas.h"
#include "aircraftControl_smallerStateEKF_sfun.h"
#include "c2_aircraftControl_smallerStateEKF.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "aircraftControl_smallerStateEKF_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c2_debug_family_names[89] = { "g", "vs", "R2D", "rho", "S",
  "m", "b", "cbar", "Ix", "Iy", "Iz", "Ixz", "xCG_ref", "xCG", "N", "E", "D",
  "vN", "vE", "vD", "phi", "theta", "psi", "p", "q", "r", "H", "throttle",
  "deltaA", "deltaE", "deltaR", "Cbn", "Cnb", "wind", "vrelN", "vrelE", "vrelD",
  "Vt", "v_body", "windOmega", "omegaAero", "alpha", "beta", "dE", "dA", "dR",
  "phidot", "pdot", "thetadot", "qdot", "psidot", "rdot", "qbar", "hT", "Tmax",
  "T", "CX", "CY", "CZ", "Cl", "Cm", "Cn", "anav", "vNdot", "vEdot", "vDdot",
  "navAccel", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "heng",
  "Ndot", "Edot", "Ddot", "nargin", "nargout", "y", "u", "loopFlag", "yDots",
  "abody", "ALP", "BET" };

static const char * c2_b_debug_family_names[10] = { "cosph", "sinph", "costh",
  "sinth", "cosps", "sinps", "nargin", "nargout", "euler", "Cbn" };

/* Function Declarations */
static void initialize_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void initialize_params_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void enable_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void disable_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void c2_update_debugger_state_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void set_sim_state_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_st);
static void finalize_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void sf_gateway_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void c2_chartstep_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void initSimStructsc2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber);
static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData);
static real_T c2_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_BET, const char_T *c2_identifier);
static real_T c2_b_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_c_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_abody, const char_T *c2_identifier, real_T c2_y[3]);
static void c2_d_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[3]);
static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_e_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_yDots, const char_T *c2_identifier, real_T c2_y[12]);
static void c2_f_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[12]);
static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_g_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[9]);
static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_info_helper(const mxArray **c2_info);
static const mxArray *c2_emlrt_marshallOut(const char * c2_u);
static const mxArray *c2_b_emlrt_marshallOut(const uint32_T c2_u);
static void c2_eul2cbn(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_euler[3], real_T c2_Cbn[9]);
static real_T c2_mpower(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_a);
static void c2_eml_scalar_eg(SFc2_aircraftControl_smallerStateEKFInstanceStruct *
  chartInstance);
static real_T c2_sqrt(SFc2_aircraftControl_smallerStateEKFInstanceStruct
                      *chartInstance, real_T c2_x);
static void c2_eml_error(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance);
static void c2_b_eml_scalar_eg
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);
static void c2_eml_xgemm(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_A[9], real_T c2_B[3], real_T c2_C[3], real_T c2_b_C
  [3]);
static real_T c2_asin(SFc2_aircraftControl_smallerStateEKFInstanceStruct
                      *chartInstance, real_T c2_x);
static void c2_b_eml_error(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance);
static real_T c2_b_mpower(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_a);
static real_T c2_c_mpower(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_a);
static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static int32_T c2_h_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static uint8_T c2_i_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_b_is_active_c2_aircraftControl_smallerStateEKF, const char_T
   *c2_identifier);
static uint8_T c2_j_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_b_sqrt(SFc2_aircraftControl_smallerStateEKFInstanceStruct
                      *chartInstance, real_T *c2_x);
static void c2_b_eml_xgemm(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_A[9], real_T c2_B[3], real_T c2_C[3]);
static void c2_b_asin(SFc2_aircraftControl_smallerStateEKFInstanceStruct
                      *chartInstance, real_T *c2_x);
static void init_dsm_address_info
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  chartInstance->c2_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c2_is_active_c2_aircraftControl_smallerStateEKF = 0U;
}

static void initialize_params_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c2_update_debugger_state_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_y = NULL;
  real_T c2_hoistedGlobal;
  real_T c2_u;
  const mxArray *c2_b_y = NULL;
  real_T c2_b_hoistedGlobal;
  real_T c2_b_u;
  const mxArray *c2_c_y = NULL;
  int32_T c2_i0;
  real_T c2_c_u[3];
  const mxArray *c2_d_y = NULL;
  int32_T c2_i1;
  real_T c2_d_u[12];
  const mxArray *c2_e_y = NULL;
  uint8_T c2_c_hoistedGlobal;
  uint8_T c2_e_u;
  const mxArray *c2_f_y = NULL;
  real_T *c2_ALP;
  real_T *c2_BET;
  real_T (*c2_yDots)[12];
  real_T (*c2_abody)[3];
  c2_BET = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_ALP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_abody = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c2_yDots = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
  c2_st = NULL;
  c2_st = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createcellmatrix(5, 1), false);
  c2_hoistedGlobal = *c2_ALP;
  c2_u = c2_hoistedGlobal;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 0, c2_b_y);
  c2_b_hoistedGlobal = *c2_BET;
  c2_b_u = c2_b_hoistedGlobal;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_b_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 1, c2_c_y);
  for (c2_i0 = 0; c2_i0 < 3; c2_i0++) {
    c2_c_u[c2_i0] = (*c2_abody)[c2_i0];
  }

  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", c2_c_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c2_y, 2, c2_d_y);
  for (c2_i1 = 0; c2_i1 < 12; c2_i1++) {
    c2_d_u[c2_i1] = (*c2_yDots)[c2_i1];
  }

  c2_e_y = NULL;
  sf_mex_assign(&c2_e_y, sf_mex_create("y", c2_d_u, 0, 0U, 1U, 0U, 1, 12), false);
  sf_mex_setcell(c2_y, 3, c2_e_y);
  c2_c_hoistedGlobal =
    chartInstance->c2_is_active_c2_aircraftControl_smallerStateEKF;
  c2_e_u = c2_c_hoistedGlobal;
  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", &c2_e_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_y, 4, c2_f_y);
  sf_mex_assign(&c2_st, c2_y, false);
  return c2_st;
}

static void set_sim_state_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_st)
{
  const mxArray *c2_u;
  real_T c2_dv0[3];
  int32_T c2_i2;
  real_T c2_dv1[12];
  int32_T c2_i3;
  real_T *c2_ALP;
  real_T *c2_BET;
  real_T (*c2_abody)[3];
  real_T (*c2_yDots)[12];
  c2_BET = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_ALP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_abody = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c2_yDots = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c2_doneDoubleBufferReInit = true;
  c2_u = sf_mex_dup(c2_st);
  *c2_ALP = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 0)),
    "ALP");
  *c2_BET = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 1)),
    "BET");
  c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 2)),
                        "abody", c2_dv0);
  for (c2_i2 = 0; c2_i2 < 3; c2_i2++) {
    (*c2_abody)[c2_i2] = c2_dv0[c2_i2];
  }

  c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 3)),
                        "yDots", c2_dv1);
  for (c2_i3 = 0; c2_i3 < 12; c2_i3++) {
    (*c2_yDots)[c2_i3] = c2_dv1[c2_i3];
  }

  chartInstance->c2_is_active_c2_aircraftControl_smallerStateEKF =
    c2_i_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 4)),
    "is_active_c2_aircraftControl_smallerStateEKF");
  sf_mex_destroy(&c2_u);
  c2_update_debugger_state_c2_aircraftControl_smallerStateEKF(chartInstance);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  int32_T c2_i4;
  int32_T c2_i5;
  int32_T c2_i6;
  int32_T c2_i7;
  int32_T c2_i8;
  real_T *c2_ALP;
  real_T *c2_BET;
  real_T (*c2_loopFlag)[3];
  real_T (*c2_abody)[3];
  real_T (*c2_yDots)[12];
  real_T (*c2_u)[4];
  real_T (*c2_y)[12];
  c2_BET = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_ALP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_loopFlag = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
  c2_abody = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c2_yDots = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
  c2_u = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 1);
  c2_y = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  for (c2_i4 = 0; c2_i4 < 12; c2_i4++) {
    _SFD_DATA_RANGE_CHECK((*c2_y)[c2_i4], 0U);
  }

  for (c2_i5 = 0; c2_i5 < 4; c2_i5++) {
    _SFD_DATA_RANGE_CHECK((*c2_u)[c2_i5], 1U);
  }

  chartInstance->c2_sfEvent = CALL_EVENT;
  c2_chartstep_c2_aircraftControl_smallerStateEKF(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY
    (_aircraftControl_smallerStateEKFMachineNumber_, chartInstance->chartNumber,
     chartInstance->instanceNumber);
  for (c2_i6 = 0; c2_i6 < 12; c2_i6++) {
    _SFD_DATA_RANGE_CHECK((*c2_yDots)[c2_i6], 2U);
  }

  for (c2_i7 = 0; c2_i7 < 3; c2_i7++) {
    _SFD_DATA_RANGE_CHECK((*c2_abody)[c2_i7], 3U);
  }

  for (c2_i8 = 0; c2_i8 < 3; c2_i8++) {
    _SFD_DATA_RANGE_CHECK((*c2_loopFlag)[c2_i8], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*c2_ALP, 5U);
  _SFD_DATA_RANGE_CHECK(*c2_BET, 6U);
}

static void c2_chartstep_c2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  int32_T c2_i9;
  real_T c2_y[12];
  int32_T c2_i10;
  real_T c2_u[4];
  int32_T c2_i11;
  real_T c2_loopFlag[3];
  uint32_T c2_debug_family_var_map[89];
  real_T c2_g;
  real_T c2_vs;
  real_T c2_R2D;
  real_T c2_rho;
  real_T c2_S;
  real_T c2_m;
  real_T c2_b;
  real_T c2_cbar;
  real_T c2_Ix;
  real_T c2_Iy;
  real_T c2_Iz;
  real_T c2_Ixz;
  real_T c2_xCG_ref;
  real_T c2_xCG;
  real_T c2_N;
  real_T c2_E;
  real_T c2_D;
  real_T c2_vN;
  real_T c2_vE;
  real_T c2_vD;
  real_T c2_phi;
  real_T c2_theta;
  real_T c2_psi;
  real_T c2_p;
  real_T c2_q;
  real_T c2_r;
  real_T c2_H;
  real_T c2_throttle;
  real_T c2_deltaA;
  real_T c2_deltaE;
  real_T c2_deltaR;
  real_T c2_Cbn[9];
  real_T c2_Cnb[9];
  real_T c2_wind[3];
  real_T c2_vrelN;
  real_T c2_vrelE;
  real_T c2_vrelD;
  real_T c2_Vt;
  real_T c2_v_body[3];
  real_T c2_windOmega[3];
  real_T c2_omegaAero[3];
  real_T c2_alpha;
  real_T c2_beta;
  real_T c2_dE;
  real_T c2_dA;
  real_T c2_dR;
  real_T c2_phidot;
  real_T c2_pdot;
  real_T c2_thetadot;
  real_T c2_qdot;
  real_T c2_psidot;
  real_T c2_rdot;
  real_T c2_qbar;
  real_T c2_hT;
  real_T c2_Tmax;
  real_T c2_T;
  real_T c2_CX;
  real_T c2_CY;
  real_T c2_CZ;
  real_T c2_Cl;
  real_T c2_Cm;
  real_T c2_Cn;
  real_T c2_anav[3];
  real_T c2_vNdot;
  real_T c2_vEdot;
  real_T c2_vDdot;
  real_T c2_navAccel[3];
  real_T c2_c1;
  real_T c2_c2;
  real_T c2_c3;
  real_T c2_c4;
  real_T c2_c5;
  real_T c2_c6;
  real_T c2_c7;
  real_T c2_c8;
  real_T c2_c9;
  real_T c2_heng;
  real_T c2_Ndot;
  real_T c2_Edot;
  real_T c2_Ddot;
  real_T c2_nargin = 3.0;
  real_T c2_nargout = 4.0;
  real_T c2_yDots[12];
  real_T c2_abody[3];
  real_T c2_ALP;
  real_T c2_BET;
  int32_T c2_i12;
  real_T c2_b_phi[3];
  real_T c2_dv2[9];
  int32_T c2_i13;
  int32_T c2_i14;
  int32_T c2_i15;
  int32_T c2_i16;
  int32_T c2_i17;
  int32_T c2_i18;
  int32_T c2_i19;
  real_T c2_a[9];
  real_T c2_b_b[3];
  int32_T c2_i20;
  int32_T c2_i21;
  int32_T c2_i22;
  real_T c2_dv3[9];
  int32_T c2_i23;
  real_T c2_dv4[3];
  int32_T c2_i24;
  real_T c2_dv5[9];
  int32_T c2_i25;
  real_T c2_dv6[3];
  int32_T c2_i26;
  real_T c2_A;
  real_T c2_B;
  real_T c2_x;
  real_T c2_b_y;
  real_T c2_b_x;
  real_T c2_c_y;
  real_T c2_c_x;
  real_T c2_d_y;
  real_T c2_e_y;
  real_T c2_d_x;
  real_T c2_e_x;
  real_T c2_b_A;
  real_T c2_b_B;
  real_T c2_f_x;
  real_T c2_f_y;
  real_T c2_g_x;
  real_T c2_g_y;
  real_T c2_h_x;
  real_T c2_h_y;
  real_T c2_i_y;
  real_T c2_c_A;
  real_T c2_i_x;
  real_T c2_j_x;
  real_T c2_k_x;
  real_T c2_d_A;
  real_T c2_l_x;
  real_T c2_m_x;
  real_T c2_n_x;
  real_T c2_j_y;
  real_T c2_e_A;
  real_T c2_o_x;
  real_T c2_p_x;
  real_T c2_q_x;
  real_T c2_k_y;
  real_T c2_f_A;
  real_T c2_r_x;
  real_T c2_s_x;
  real_T c2_t_x;
  real_T c2_l_y;
  real_T c2_g_A;
  real_T c2_u_x;
  real_T c2_v_x;
  real_T c2_w_x;
  real_T c2_m_y;
  real_T c2_h_A;
  real_T c2_x_x;
  real_T c2_y_x;
  real_T c2_ab_x;
  real_T c2_i_A;
  real_T c2_bb_x;
  real_T c2_cb_x;
  real_T c2_db_x;
  real_T c2_n_y;
  real_T c2_j_A;
  real_T c2_eb_x;
  real_T c2_fb_x;
  real_T c2_gb_x;
  real_T c2_o_y;
  real_T c2_k_A;
  real_T c2_c_B;
  real_T c2_hb_x;
  real_T c2_p_y;
  real_T c2_ib_x;
  real_T c2_q_y;
  real_T c2_jb_x;
  real_T c2_r_y;
  real_T c2_s_y;
  real_T c2_d_B;
  real_T c2_t_y;
  real_T c2_u_y;
  real_T c2_v_y;
  real_T c2_w_y;
  real_T c2_l_A;
  real_T c2_kb_x;
  real_T c2_lb_x;
  real_T c2_mb_x;
  real_T c2_x_y;
  real_T c2_m_A;
  real_T c2_nb_x;
  real_T c2_ob_x;
  real_T c2_pb_x;
  real_T c2_y_y;
  real_T c2_n_A;
  real_T c2_e_B;
  real_T c2_qb_x;
  real_T c2_ab_y;
  real_T c2_rb_x;
  real_T c2_bb_y;
  real_T c2_sb_x;
  real_T c2_cb_y;
  real_T c2_db_y;
  real_T c2_f_B;
  real_T c2_eb_y;
  real_T c2_fb_y;
  real_T c2_gb_y;
  real_T c2_hb_y;
  real_T c2_o_A;
  real_T c2_tb_x;
  real_T c2_ub_x;
  real_T c2_vb_x;
  real_T c2_ib_y;
  real_T c2_p_A;
  real_T c2_wb_x;
  real_T c2_xb_x;
  real_T c2_yb_x;
  real_T c2_jb_y;
  real_T c2_q_A;
  real_T c2_g_B;
  real_T c2_ac_x;
  real_T c2_kb_y;
  real_T c2_bc_x;
  real_T c2_lb_y;
  real_T c2_cc_x;
  real_T c2_mb_y;
  real_T c2_nb_y;
  real_T c2_h_B;
  real_T c2_ob_y;
  real_T c2_pb_y;
  real_T c2_qb_y;
  real_T c2_rb_y;
  int32_T c2_i27;
  int32_T c2_i28;
  int32_T c2_i29;
  int32_T c2_i30;
  int32_T c2_i31;
  int32_T c2_i32;
  real_T c2_dv7[9];
  int32_T c2_i33;
  real_T c2_dv8[3];
  int32_T c2_i34;
  real_T c2_dv9[9];
  int32_T c2_i35;
  real_T c2_dv10[3];
  real_T c2_dv11[3];
  int32_T c2_i36;
  real_T c2_dc_x;
  real_T c2_ec_x;
  real_T c2_fc_x;
  real_T c2_gc_x;
  real_T c2_hc_x;
  real_T c2_ic_x;
  real_T c2_jc_x;
  real_T c2_kc_x;
  real_T c2_lc_x;
  real_T c2_mc_x;
  real_T c2_nc_x;
  real_T c2_oc_x;
  real_T c2_pc_x;
  real_T c2_qc_x;
  real_T c2_rc_x;
  real_T c2_sc_x;
  real_T c2_r_A;
  real_T c2_i_B;
  real_T c2_tc_x;
  real_T c2_sb_y;
  real_T c2_uc_x;
  real_T c2_tb_y;
  real_T c2_vc_x;
  real_T c2_ub_y;
  int32_T c2_i37;
  int32_T c2_i38;
  real_T *c2_b_BET;
  real_T *c2_b_ALP;
  real_T (*c2_b_yDots)[12];
  real_T (*c2_b_abody)[3];
  real_T (*c2_b_loopFlag)[3];
  real_T (*c2_b_u)[4];
  real_T (*c2_vb_y)[12];
  c2_b_BET = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c2_b_ALP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c2_b_loopFlag = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
  c2_b_abody = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c2_b_yDots = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
  c2_b_u = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 1);
  c2_vb_y = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  for (c2_i9 = 0; c2_i9 < 12; c2_i9++) {
    c2_y[c2_i9] = (*c2_vb_y)[c2_i9];
  }

  for (c2_i10 = 0; c2_i10 < 4; c2_i10++) {
    c2_u[c2_i10] = (*c2_b_u)[c2_i10];
  }

  for (c2_i11 = 0; c2_i11 < 3; c2_i11++) {
    c2_loopFlag[c2_i11] = (*c2_b_loopFlag)[c2_i11];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 89U, 89U, c2_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_g, 0U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_vs, 1U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_R2D, 2U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_rho, 3U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_S, 4U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_m, 5U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_b, 6U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_cbar, 7U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_Ix, 8U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_Iy, 9U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_Iz, 10U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_Ixz, 11U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_xCG_ref, 12U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_xCG, 13U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_N, 14U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_E, 15U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_D, 16U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vN, 17U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vE, 18U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vD, 19U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_phi, 20U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_theta, 21U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_psi, 22U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_p, 23U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_q, 24U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_r, 25U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_H, 26U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_throttle, 27U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_deltaA, 28U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_deltaE, 29U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_deltaR, 30U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Cbn, 31U, c2_e_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Cnb, 32U, c2_e_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_wind, 33U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vrelN, 34U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vrelE, 35U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vrelD, 36U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Vt, 37U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_v_body, 38U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_windOmega, 39U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_omegaAero, 40U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_alpha, 41U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_beta, 42U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_dE, 43U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_dA, 44U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_dR, 45U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_phidot, 46U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_pdot, 47U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_thetadot, 48U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_qdot, 49U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_psidot, 50U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_rdot, 51U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_qbar, 52U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_hT, 53U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Tmax, 54U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_T, 55U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_CX, 56U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_CY, 57U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_CZ, 58U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Cl, 59U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Cm, 60U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Cn, 61U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_anav, 62U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vNdot, 63U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vEdot, 64U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_vDdot, 65U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_navAccel, 66U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c1, 67U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c2, 68U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c3, 69U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c4, 70U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c5, 71U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c6, 72U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c7, 73U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c8, 74U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_c9, 75U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_heng, 76U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Ndot, 77U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Edot, 78U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_Ddot, 79U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 80U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 81U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_y, 82U, c2_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_u, 83U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_loopFlag, 84U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_yDots, 85U, c2_c_sf_marshallOut,
    c2_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_abody, 86U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_ALP, 87U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_BET, 88U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 4);
  for (c2_i12 = 0; c2_i12 < 12; c2_i12++) {
    c2_yDots[c2_i12] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 7);
  c2_g = 9.81;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 8);
  c2_vs = 340.3;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 9);
  c2_R2D = 57.295779513082323;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 10);
  c2_rho = 1.225;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 13);
  c2_S = 20.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 14);
  c2_m = 1177.0410005333333;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 15);
  c2_b = 6.666666666666667;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 16);
  c2_cbar = 3.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 17);
  c2_Ix = 2256.9839826075154;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 18);
  c2_Iy = 11044.488299351713;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 19);
  c2_Iz = 12636.21789221188;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 20);
  c2_Ixz = 106.20569401537166;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 21);
  c2_xCG_ref = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 22);
  c2_xCG = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 25);
  c2_N = c2_y[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 26);
  c2_E = c2_y[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 27);
  c2_D = c2_y[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 28);
  c2_vN = c2_y[3];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 29);
  c2_vE = c2_y[4];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 30);
  c2_vD = c2_y[5];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 31);
  c2_phi = c2_y[6];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 32);
  c2_theta = c2_y[7];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 33);
  c2_psi = c2_y[8];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 34);
  c2_p = c2_y[9];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 35);
  c2_q = c2_y[10];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 36);
  c2_r = c2_y[11];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 38);
  c2_H = -c2_D;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 40);
  c2_throttle = c2_u[3];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 41);
  c2_deltaA = c2_u[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 42);
  c2_deltaE = c2_u[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 43);
  c2_deltaR = c2_u[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 47);
  c2_b_phi[0] = c2_phi;
  c2_b_phi[1] = c2_theta;
  c2_b_phi[2] = c2_psi;
  c2_eul2cbn(chartInstance, c2_b_phi, c2_dv2);
  for (c2_i13 = 0; c2_i13 < 9; c2_i13++) {
    c2_Cbn[c2_i13] = c2_dv2[c2_i13];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 48);
  c2_i14 = 0;
  for (c2_i15 = 0; c2_i15 < 3; c2_i15++) {
    c2_i16 = 0;
    for (c2_i17 = 0; c2_i17 < 3; c2_i17++) {
      c2_Cnb[c2_i17 + c2_i14] = c2_Cbn[c2_i16 + c2_i15];
      c2_i16 += 3;
    }

    c2_i14 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 50);
  for (c2_i18 = 0; c2_i18 < 3; c2_i18++) {
    c2_wind[c2_i18] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 51);
  c2_vrelN = c2_vN - c2_wind[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 52);
  c2_vrelE = c2_vE - c2_wind[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 53);
  c2_vrelD = c2_vD - c2_wind[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 55);
  c2_Vt = (c2_mpower(chartInstance, c2_vrelN) + c2_mpower(chartInstance,
            c2_vrelE)) + c2_mpower(chartInstance, c2_vrelD);
  c2_b_sqrt(chartInstance, &c2_Vt);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 57);
  for (c2_i19 = 0; c2_i19 < 9; c2_i19++) {
    c2_a[c2_i19] = c2_Cnb[c2_i19];
  }

  c2_b_b[0] = c2_vrelN;
  c2_b_b[1] = c2_vrelE;
  c2_b_b[2] = c2_vrelD;
  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i20 = 0; c2_i20 < 3; c2_i20++) {
    c2_v_body[c2_i20] = 0.0;
  }

  for (c2_i21 = 0; c2_i21 < 3; c2_i21++) {
    c2_v_body[c2_i21] = 0.0;
  }

  for (c2_i22 = 0; c2_i22 < 9; c2_i22++) {
    c2_dv3[c2_i22] = c2_a[c2_i22];
  }

  for (c2_i23 = 0; c2_i23 < 3; c2_i23++) {
    c2_dv4[c2_i23] = c2_b_b[c2_i23];
  }

  for (c2_i24 = 0; c2_i24 < 9; c2_i24++) {
    c2_dv5[c2_i24] = c2_dv3[c2_i24];
  }

  for (c2_i25 = 0; c2_i25 < 3; c2_i25++) {
    c2_dv6[c2_i25] = c2_dv4[c2_i25];
  }

  c2_b_eml_xgemm(chartInstance, c2_dv5, c2_dv6, c2_v_body);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 59);
  for (c2_i26 = 0; c2_i26 < 3; c2_i26++) {
    c2_windOmega[c2_i26] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 60);
  c2_omegaAero[0] = c2_p + c2_windOmega[0];
  c2_omegaAero[1] = c2_q + c2_windOmega[1];
  c2_omegaAero[2] = c2_r + c2_windOmega[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 63);
  c2_A = c2_v_body[2];
  c2_B = c2_v_body[0];
  c2_x = c2_A;
  c2_b_y = c2_B;
  c2_b_x = c2_x;
  c2_c_y = c2_b_y;
  c2_c_x = c2_b_x;
  c2_d_y = c2_c_y;
  c2_e_y = c2_c_x / c2_d_y;
  c2_d_x = c2_e_y;
  c2_alpha = c2_d_x;
  c2_e_x = c2_alpha;
  c2_alpha = c2_e_x;
  c2_alpha = muDoubleScalarAtan(c2_alpha);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 64);
  c2_b_A = c2_v_body[1];
  c2_b_B = c2_Vt;
  c2_f_x = c2_b_A;
  c2_f_y = c2_b_B;
  c2_g_x = c2_f_x;
  c2_g_y = c2_f_y;
  c2_h_x = c2_g_x;
  c2_h_y = c2_g_y;
  c2_i_y = c2_h_x / c2_h_y;
  c2_beta = c2_i_y;
  c2_b_asin(chartInstance, &c2_beta);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 66);
  c2_ALP = c2_alpha * 57.295779513082323;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 67);
  c2_BET = c2_beta * 57.295779513082323;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 68);
  c2_dE = c2_deltaE * 57.295779513082323;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 69);
  c2_dA = c2_deltaA * 57.295779513082323;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 70);
  c2_dR = c2_deltaR * 57.295779513082323;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 74);
  if (CV_EML_IF(0, 1, 0, c2_loopFlag[0] == 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 75);
    c2_dA = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 76);
    c2_phidot = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 77);
    c2_pdot = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 80);
  if (CV_EML_IF(0, 1, 1, c2_loopFlag[1] == 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 81);
    c2_dE = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 82);
    c2_thetadot = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 83);
    c2_qdot = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 87);
  if (CV_EML_IF(0, 1, 2, c2_loopFlag[2] == 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 88);
    c2_dR = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 89);
    c2_psidot = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 90);
    c2_rdot = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 95);
  c2_qbar = 0.6125 * c2_mpower(chartInstance, c2_Vt);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 98);
  c2_c_A = c2_H;
  c2_i_x = c2_c_A;
  c2_j_x = c2_i_x;
  c2_k_x = c2_j_x;
  c2_hT = c2_k_x / 3048.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 100);
  c2_d_A = c2_Vt;
  c2_l_x = c2_d_A;
  c2_m_x = c2_l_x;
  c2_n_x = c2_m_x;
  c2_j_y = c2_n_x / 340.3;
  c2_e_A = c2_Vt;
  c2_o_x = c2_e_A;
  c2_p_x = c2_o_x;
  c2_q_x = c2_p_x;
  c2_k_y = c2_q_x / 340.3;
  c2_f_A = c2_Vt;
  c2_r_x = c2_f_A;
  c2_s_x = c2_r_x;
  c2_t_x = c2_s_x;
  c2_l_y = c2_t_x / 340.3;
  c2_g_A = c2_Vt;
  c2_u_x = c2_g_A;
  c2_v_x = c2_u_x;
  c2_w_x = c2_v_x;
  c2_m_y = c2_w_x / 340.3;
  c2_h_A = ((((((((30.21 - 0.668 * c2_hT) - 6.877 * c2_mpower(chartInstance,
    c2_hT)) + 1.951 * c2_b_mpower(chartInstance, c2_hT)) - 0.1512 * c2_c_mpower
                (chartInstance, c2_hT)) + c2_j_y * ((((-33.8 + 3.347 * c2_hT) +
    18.13 * c2_mpower(chartInstance, c2_hT)) - 5.865 * c2_b_mpower(chartInstance,
    c2_hT)) + 0.4757 * c2_c_mpower(chartInstance, c2_hT))) + c2_mpower
              (chartInstance, c2_k_y) * ((((100.8 - 77.56 * c2_hT) + 5.441 *
    c2_mpower(chartInstance, c2_hT)) + 2.864 * c2_b_mpower(chartInstance, c2_hT))
    - 0.3355 * c2_c_mpower(chartInstance, c2_hT))) + c2_b_mpower(chartInstance,
              c2_l_y) * ((((-78.99 + 101.4 * c2_hT) - 30.28 * c2_mpower
    (chartInstance, c2_hT)) + 3.236 * c2_b_mpower(chartInstance, c2_hT)) -
              0.1089 * c2_c_mpower(chartInstance, c2_hT))) + c2_c_mpower
            (chartInstance, c2_m_y) * ((((18.74 - 31.6 * c2_hT) + 12.04 *
    c2_mpower(chartInstance, c2_hT)) - 1.785 * c2_b_mpower(chartInstance, c2_hT))
             + 0.09417 * c2_c_mpower(chartInstance, c2_hT))) * 4448.22;
  c2_x_x = c2_h_A;
  c2_y_x = c2_x_x;
  c2_ab_x = c2_y_x;
  c2_Tmax = c2_ab_x / 20.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 106);
  c2_T = c2_Tmax * c2_throttle;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 109);
  c2_i_A = 180.0 * c2_omegaAero[1] * 3.0;
  c2_bb_x = c2_i_A;
  c2_cb_x = c2_bb_x;
  c2_db_x = c2_cb_x;
  c2_n_y = c2_db_x / 3.1415926535897931;
  c2_j_A = c2_n_y;
  c2_eb_x = c2_j_A;
  c2_fb_x = c2_eb_x;
  c2_gb_x = c2_fb_x;
  c2_o_y = c2_gb_x / 2.0;
  c2_k_A = c2_o_y;
  c2_c_B = c2_Vt;
  c2_hb_x = c2_k_A;
  c2_p_y = c2_c_B;
  c2_ib_x = c2_hb_x;
  c2_q_y = c2_p_y;
  c2_jb_x = c2_ib_x;
  c2_r_y = c2_q_y;
  c2_s_y = c2_jb_x / c2_r_y;
  c2_CX = (((((-0.0434 + 0.00239 * c2_ALP) + 2.53E-5 * c2_mpower(chartInstance,
    c2_BET)) - 1.07E-6 * c2_ALP * c2_mpower(chartInstance, c2_BET)) + 0.00095 *
            c2_dE) - 8.5E-7 * c2_dE * c2_mpower(chartInstance, c2_BET)) + c2_s_y
    * ((0.00873 + 0.001 * c2_ALP) - 0.000175 * c2_mpower(chartInstance, c2_ALP));
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 110);
  c2_d_B = c2_Vt;
  c2_t_y = c2_d_B;
  c2_u_y = c2_t_y;
  c2_v_y = c2_u_y;
  c2_w_y = 190.9859317102744 / c2_v_y;
  c2_CY = ((-0.012 * c2_BET + 0.00155 * c2_dR) - 8.0E-6 * c2_dR * c2_ALP) +
    c2_w_y * (((0.00225 * c2_omegaAero[0] + 0.0117 * c2_omegaAero[2]) - 0.000367
               * c2_omegaAero[2] * c2_ALP) + 0.000175 * c2_omegaAero[2] * c2_dE);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 111);
  c2_l_A = 180.0 * c2_omegaAero[1] * 3.0;
  c2_kb_x = c2_l_A;
  c2_lb_x = c2_kb_x;
  c2_mb_x = c2_lb_x;
  c2_x_y = c2_mb_x / 3.1415926535897931;
  c2_m_A = c2_x_y;
  c2_nb_x = c2_m_A;
  c2_ob_x = c2_nb_x;
  c2_pb_x = c2_ob_x;
  c2_y_y = c2_pb_x / 2.0;
  c2_n_A = c2_y_y;
  c2_e_B = c2_Vt;
  c2_qb_x = c2_n_A;
  c2_ab_y = c2_e_B;
  c2_rb_x = c2_qb_x;
  c2_bb_y = c2_ab_y;
  c2_sb_x = c2_rb_x;
  c2_cb_y = c2_bb_y;
  c2_db_y = c2_sb_x / c2_cb_y;
  c2_CZ = ((((-0.131 - 0.0538 * c2_ALP) - 0.00476 * c2_dE) - 3.3E-5 * c2_dE *
            c2_ALP) - 7.5E-5 * c2_mpower(chartInstance, c2_dA)) + c2_db_y *
    ((-0.111 + 0.00517 * c2_ALP) - 0.0011 * c2_mpower(chartInstance, c2_ALP));
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 112);
  c2_f_B = c2_Vt;
  c2_eb_y = c2_f_B;
  c2_fb_y = c2_eb_y;
  c2_gb_y = c2_fb_y;
  c2_hb_y = 190.9859317102744 / c2_gb_y;
  c2_Cl = ((((-0.000598 * c2_BET - 0.000283 * c2_ALP * c2_BET) + 1.51E-5 *
             c2_mpower(chartInstance, c2_ALP) * c2_BET) - c2_dA * ((0.00061 +
              2.5E-5 * c2_ALP) - 2.6E-6 * c2_mpower(chartInstance, c2_ALP))) -
           c2_dR * (-0.00023 + 4.5E-6 * c2_ALP)) + c2_hb_y * (((((-0.00412 *
    c2_omegaAero[0] - 0.000524 * c2_omegaAero[0] * c2_ALP) + 4.36E-5 *
    c2_omegaAero[0] * c2_mpower(chartInstance, c2_ALP)) + 0.000436 *
    c2_omegaAero[2]) + 0.000105 * c2_omegaAero[2] * c2_ALP) + 5.24E-5 *
    c2_omegaAero[2] * c2_dE);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 113);
  c2_o_A = 180.0 * c2_omegaAero[1] * 3.0;
  c2_tb_x = c2_o_A;
  c2_ub_x = c2_tb_x;
  c2_vb_x = c2_ub_x;
  c2_ib_y = c2_vb_x / 2.0;
  c2_p_A = c2_ib_y;
  c2_wb_x = c2_p_A;
  c2_xb_x = c2_wb_x;
  c2_yb_x = c2_xb_x;
  c2_jb_y = c2_yb_x / 3.1415926535897931;
  c2_q_A = c2_jb_y;
  c2_g_B = c2_Vt;
  c2_ac_x = c2_q_A;
  c2_kb_y = c2_g_B;
  c2_bc_x = c2_ac_x;
  c2_lb_y = c2_kb_y;
  c2_cc_x = c2_bc_x;
  c2_mb_y = c2_lb_y;
  c2_nb_y = c2_cc_x / c2_mb_y;
  c2_Cm = ((((((((-0.00661 - 0.00267 * c2_ALP) - 6.48E-5 * c2_mpower
                 (chartInstance, c2_BET)) - 2.65E-6 * c2_ALP * c2_mpower
                (chartInstance, c2_BET)) - 0.00654 * c2_dE) - 8.49E-5 * c2_dE *
              c2_ALP) + 3.74E-6 * c2_dE * c2_mpower(chartInstance, c2_BET)) -
            3.5E-5 * c2_mpower(chartInstance, c2_dA)) + c2_nb_y * (-0.0473 -
            0.00157 * c2_ALP)) + 0.0 * c2_CZ;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 114);
  c2_h_B = c2_Vt;
  c2_ob_y = c2_h_B;
  c2_pb_y = c2_ob_y;
  c2_qb_y = c2_pb_y;
  c2_rb_y = 190.9859317102744 / c2_qb_y;
  c2_Cn = ((((((0.00228 * c2_BET + 1.79E-6 * c2_b_mpower(chartInstance, c2_BET))
               + 1.4E-5 * c2_dA) + 7.0E-6 * c2_dA * c2_ALP) - 0.0009 * c2_dR) +
            4.0E-6 * c2_dR * c2_ALP) + c2_rb_y * (((((-6.63E-5 * c2_omegaAero[0]
    - 1.92E-5 * c2_omegaAero[0] * c2_ALP) + 5.06E-6 * c2_omegaAero[0] *
    c2_mpower(chartInstance, c2_ALP)) - 0.00606 * c2_omegaAero[2]) - 8.73E-5 *
             c2_omegaAero[2] * c2_dE) + 8.7E-6 * c2_omegaAero[2] * c2_dE *
            c2_ALP)) - 0.0 * c2_CY;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 121);
  c2_b_b[0] = c2_qbar * 20.0 * c2_CX + c2_T;
  c2_b_b[1] = c2_qbar * 20.0 * c2_CY;
  c2_b_b[2] = c2_qbar * 20.0 * c2_CZ;
  for (c2_i27 = 0; c2_i27 < 3; c2_i27++) {
    c2_abody[c2_i27] = c2_b_b[c2_i27] / 1177.0410005333333;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 122);
  for (c2_i28 = 0; c2_i28 < 9; c2_i28++) {
    c2_a[c2_i28] = c2_Cbn[c2_i28];
  }

  for (c2_i29 = 0; c2_i29 < 3; c2_i29++) {
    c2_b_b[c2_i29] = c2_abody[c2_i29];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i30 = 0; c2_i30 < 3; c2_i30++) {
    c2_anav[c2_i30] = 0.0;
  }

  for (c2_i31 = 0; c2_i31 < 3; c2_i31++) {
    c2_anav[c2_i31] = 0.0;
  }

  for (c2_i32 = 0; c2_i32 < 9; c2_i32++) {
    c2_dv7[c2_i32] = c2_a[c2_i32];
  }

  for (c2_i33 = 0; c2_i33 < 3; c2_i33++) {
    c2_dv8[c2_i33] = c2_b_b[c2_i33];
  }

  for (c2_i34 = 0; c2_i34 < 9; c2_i34++) {
    c2_dv9[c2_i34] = c2_dv7[c2_i34];
  }

  for (c2_i35 = 0; c2_i35 < 3; c2_i35++) {
    c2_dv10[c2_i35] = c2_dv8[c2_i35];
  }

  c2_b_eml_xgemm(chartInstance, c2_dv9, c2_dv10, c2_anav);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 124);
  c2_vNdot = c2_anav[0];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 125);
  c2_vEdot = c2_anav[1];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 126);
  c2_vDdot = c2_g + c2_anav[2];
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 128U);
  c2_dv11[0] = 0.0;
  c2_dv11[1] = 0.0;
  c2_dv11[2] = c2_g;
  for (c2_i36 = 0; c2_i36 < 3; c2_i36++) {
    c2_navAccel[c2_i36] = c2_anav[c2_i36] + c2_dv11[c2_i36];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 138U);
  c2_c1 = -0.70592099279383547;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 139U);
  c2_c2 = 0.014338034095370216;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 140U);
  c2_c3 = 0.000443244465804795;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 141U);
  c2_c4 = 3.7254094944251367E-6;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 142U);
  c2_c5 = 0.93976593829282273;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 143U);
  c2_c6 = 0.0096161715361322529;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 144U);
  c2_c7 = 9.054290003265229E-5;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 145U);
  c2_c8 = -0.69609284164731566;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 146U);
  c2_c9 = 7.916891495812398E-5;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 148U);
  c2_heng = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 149U);
  c2_pdot = (-0.70592099279383547 * c2_r + 0.014338034095370216 * c2_p) * c2_q +
    c2_qbar * 20.0 * 6.666666666666667 * (0.000443244465804795 * c2_Cl +
    3.7254094944251367E-6 * c2_Cn);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 150U);
  c2_qdot = (0.93976593829282273 * c2_p * c2_r - 0.0096161715361322529 *
             (c2_mpower(chartInstance, c2_p) - c2_mpower(chartInstance, c2_r)))
    + c2_qbar * 20.0 * 3.0 * 9.054290003265229E-5 * c2_Cm;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 151U);
  c2_rdot = (-0.69609284164731566 * c2_p - 0.014338034095370216 * c2_r) * c2_q +
    c2_qbar * 20.0 * 6.666666666666667 * (3.7254094944251367E-6 * c2_Cl +
    7.916891495812398E-5 * c2_Cn);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 157U);
  c2_Ndot = c2_vN;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 158U);
  c2_Edot = c2_vE;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 159U);
  c2_Ddot = c2_vD;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 161U);
  c2_dc_x = c2_theta;
  c2_ec_x = c2_dc_x;
  c2_ec_x = muDoubleScalarTan(c2_ec_x);
  c2_fc_x = c2_phi;
  c2_gc_x = c2_fc_x;
  c2_gc_x = muDoubleScalarSin(c2_gc_x);
  c2_hc_x = c2_phi;
  c2_ic_x = c2_hc_x;
  c2_ic_x = muDoubleScalarCos(c2_ic_x);
  c2_phidot = c2_p + c2_ec_x * (c2_q * c2_gc_x + c2_r * c2_ic_x);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 162U);
  c2_jc_x = c2_phi;
  c2_kc_x = c2_jc_x;
  c2_kc_x = muDoubleScalarCos(c2_kc_x);
  c2_lc_x = c2_phi;
  c2_mc_x = c2_lc_x;
  c2_mc_x = muDoubleScalarSin(c2_mc_x);
  c2_thetadot = c2_q * c2_kc_x - c2_r * c2_mc_x;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 163U);
  c2_nc_x = c2_phi;
  c2_oc_x = c2_nc_x;
  c2_oc_x = muDoubleScalarSin(c2_oc_x);
  c2_pc_x = c2_phi;
  c2_qc_x = c2_pc_x;
  c2_qc_x = muDoubleScalarCos(c2_qc_x);
  c2_rc_x = c2_theta;
  c2_sc_x = c2_rc_x;
  c2_sc_x = muDoubleScalarCos(c2_sc_x);
  c2_r_A = c2_q * c2_oc_x + c2_r * c2_qc_x;
  c2_i_B = c2_sc_x;
  c2_tc_x = c2_r_A;
  c2_sb_y = c2_i_B;
  c2_uc_x = c2_tc_x;
  c2_tb_y = c2_sb_y;
  c2_vc_x = c2_uc_x;
  c2_ub_y = c2_tb_y;
  c2_psidot = c2_vc_x / c2_ub_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 166U);
  if (CV_EML_IF(0, 1, 3, c2_loopFlag[0] == 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 167U);
    c2_dA = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 168U);
    c2_phidot = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 169U);
    c2_pdot = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 172U);
  if (CV_EML_IF(0, 1, 4, c2_loopFlag[1] == 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 173U);
    c2_dE = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 174U);
    c2_thetadot = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 175U);
    c2_qdot = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 179U);
  if (CV_EML_IF(0, 1, 5, c2_loopFlag[2] == 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 180U);
    c2_dR = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 181U);
    c2_psidot = 0.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 182U);
    c2_rdot = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 187U);
  c2_yDots[0] = c2_Ndot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 188U);
  c2_yDots[1] = c2_Edot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 189U);
  c2_yDots[2] = c2_Ddot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 190U);
  c2_yDots[3] = c2_vNdot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 191U);
  c2_yDots[4] = c2_vEdot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 192U);
  c2_yDots[5] = c2_vDdot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 193U);
  c2_yDots[6] = c2_phidot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 194U);
  c2_yDots[7] = c2_thetadot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 195U);
  c2_yDots[8] = c2_psidot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 196U);
  c2_yDots[9] = c2_pdot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 197U);
  c2_yDots[10] = c2_qdot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 198U);
  c2_yDots[11] = c2_rdot;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, -198);
  _SFD_SYMBOL_SCOPE_POP();
  for (c2_i37 = 0; c2_i37 < 12; c2_i37++) {
    (*c2_b_yDots)[c2_i37] = c2_yDots[c2_i37];
  }

  for (c2_i38 = 0; c2_i38 < 3; c2_i38++) {
    (*c2_b_abody)[c2_i38] = c2_abody[c2_i38];
  }

  *c2_b_ALP = c2_ALP;
  *c2_b_BET = c2_BET;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
}

static void initSimStructsc2_aircraftControl_smallerStateEKF
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber)
{
  (void)c2_machineNumber;
  (void)c2_chartNumber;
  (void)c2_instanceNumber;
}

static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  real_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(real_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static real_T c2_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_BET, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_BET), &c2_thisId);
  sf_mex_destroy(&c2_BET);
  return c2_y;
}

static real_T c2_b_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d0, 1, 0, 0U, 0, 0U, 0);
  c2_y = c2_d0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_BET;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_BET = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_BET), &c2_thisId);
  sf_mex_destroy(&c2_BET);
  *(real_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i39;
  real_T c2_b_inData[3];
  int32_T c2_i40;
  real_T c2_u[3];
  const mxArray *c2_y = NULL;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i39 = 0; c2_i39 < 3; c2_i39++) {
    c2_b_inData[c2_i39] = (*(real_T (*)[3])c2_inData)[c2_i39];
  }

  for (c2_i40 = 0; c2_i40 < 3; c2_i40++) {
    c2_u[c2_i40] = c2_b_inData[c2_i40];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_c_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_abody, const char_T *c2_identifier, real_T c2_y[3])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_abody), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_abody);
}

static void c2_d_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[3])
{
  real_T c2_dv12[3];
  int32_T c2_i41;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv12, 1, 0, 0U, 1, 0U, 1, 3);
  for (c2_i41 = 0; c2_i41 < 3; c2_i41++) {
    c2_y[c2_i41] = c2_dv12[c2_i41];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_abody;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[3];
  int32_T c2_i42;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_abody = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_abody), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_abody);
  for (c2_i42 = 0; c2_i42 < 3; c2_i42++) {
    (*(real_T (*)[3])c2_outData)[c2_i42] = c2_y[c2_i42];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i43;
  real_T c2_b_inData[12];
  int32_T c2_i44;
  real_T c2_u[12];
  const mxArray *c2_y = NULL;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i43 = 0; c2_i43 < 12; c2_i43++) {
    c2_b_inData[c2_i43] = (*(real_T (*)[12])c2_inData)[c2_i43];
  }

  for (c2_i44 = 0; c2_i44 < 12; c2_i44++) {
    c2_u[c2_i44] = c2_b_inData[c2_i44];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 12), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_e_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_yDots, const char_T *c2_identifier, real_T c2_y[12])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_yDots), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_yDots);
}

static void c2_f_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[12])
{
  real_T c2_dv13[12];
  int32_T c2_i45;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv13, 1, 0, 0U, 1, 0U, 1, 12);
  for (c2_i45 = 0; c2_i45 < 12; c2_i45++) {
    c2_y[c2_i45] = c2_dv13[c2_i45];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_yDots;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[12];
  int32_T c2_i46;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_yDots = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_yDots), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_yDots);
  for (c2_i46 = 0; c2_i46 < 12; c2_i46++) {
    (*(real_T (*)[12])c2_outData)[c2_i46] = c2_y[c2_i46];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i47;
  real_T c2_b_inData[4];
  int32_T c2_i48;
  real_T c2_u[4];
  const mxArray *c2_y = NULL;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i47 = 0; c2_i47 < 4; c2_i47++) {
    c2_b_inData[c2_i47] = (*(real_T (*)[4])c2_inData)[c2_i47];
  }

  for (c2_i48 = 0; c2_i48 < 4; c2_i48++) {
    c2_u[c2_i48] = c2_b_inData[c2_i48];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i49;
  int32_T c2_i50;
  int32_T c2_i51;
  real_T c2_b_inData[9];
  int32_T c2_i52;
  int32_T c2_i53;
  int32_T c2_i54;
  real_T c2_u[9];
  const mxArray *c2_y = NULL;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_i49 = 0;
  for (c2_i50 = 0; c2_i50 < 3; c2_i50++) {
    for (c2_i51 = 0; c2_i51 < 3; c2_i51++) {
      c2_b_inData[c2_i51 + c2_i49] = (*(real_T (*)[9])c2_inData)[c2_i51 + c2_i49];
    }

    c2_i49 += 3;
  }

  c2_i52 = 0;
  for (c2_i53 = 0; c2_i53 < 3; c2_i53++) {
    for (c2_i54 = 0; c2_i54 < 3; c2_i54++) {
      c2_u[c2_i54 + c2_i52] = c2_b_inData[c2_i54 + c2_i52];
    }

    c2_i52 += 3;
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static void c2_g_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[9])
{
  real_T c2_dv14[9];
  int32_T c2_i55;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv14, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c2_i55 = 0; c2_i55 < 9; c2_i55++) {
    c2_y[c2_i55] = c2_dv14[c2_i55];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_Cnb;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[9];
  int32_T c2_i56;
  int32_T c2_i57;
  int32_T c2_i58;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_Cnb = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_Cnb), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_Cnb);
  c2_i56 = 0;
  for (c2_i57 = 0; c2_i57 < 3; c2_i57++) {
    for (c2_i58 = 0; c2_i58 < 3; c2_i58++) {
      (*(real_T (*)[9])c2_outData)[c2_i58 + c2_i56] = c2_y[c2_i58 + c2_i56];
    }

    c2_i56 += 3;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

const mxArray
  *sf_c2_aircraftControl_smallerStateEKF_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  sf_mex_assign(&c2_nameCaptureInfo, sf_mex_createstruct("structure", 2, 46, 1),
                false);
  c2_info_helper(&c2_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c2_nameCaptureInfo);
  return c2_nameCaptureInfo;
}

static void c2_info_helper(const mxArray **c2_info)
{
  const mxArray *c2_rhs0 = NULL;
  const mxArray *c2_lhs0 = NULL;
  const mxArray *c2_rhs1 = NULL;
  const mxArray *c2_lhs1 = NULL;
  const mxArray *c2_rhs2 = NULL;
  const mxArray *c2_lhs2 = NULL;
  const mxArray *c2_rhs3 = NULL;
  const mxArray *c2_lhs3 = NULL;
  const mxArray *c2_rhs4 = NULL;
  const mxArray *c2_lhs4 = NULL;
  const mxArray *c2_rhs5 = NULL;
  const mxArray *c2_lhs5 = NULL;
  const mxArray *c2_rhs6 = NULL;
  const mxArray *c2_lhs6 = NULL;
  const mxArray *c2_rhs7 = NULL;
  const mxArray *c2_lhs7 = NULL;
  const mxArray *c2_rhs8 = NULL;
  const mxArray *c2_lhs8 = NULL;
  const mxArray *c2_rhs9 = NULL;
  const mxArray *c2_lhs9 = NULL;
  const mxArray *c2_rhs10 = NULL;
  const mxArray *c2_lhs10 = NULL;
  const mxArray *c2_rhs11 = NULL;
  const mxArray *c2_lhs11 = NULL;
  const mxArray *c2_rhs12 = NULL;
  const mxArray *c2_lhs12 = NULL;
  const mxArray *c2_rhs13 = NULL;
  const mxArray *c2_lhs13 = NULL;
  const mxArray *c2_rhs14 = NULL;
  const mxArray *c2_lhs14 = NULL;
  const mxArray *c2_rhs15 = NULL;
  const mxArray *c2_lhs15 = NULL;
  const mxArray *c2_rhs16 = NULL;
  const mxArray *c2_lhs16 = NULL;
  const mxArray *c2_rhs17 = NULL;
  const mxArray *c2_lhs17 = NULL;
  const mxArray *c2_rhs18 = NULL;
  const mxArray *c2_lhs18 = NULL;
  const mxArray *c2_rhs19 = NULL;
  const mxArray *c2_lhs19 = NULL;
  const mxArray *c2_rhs20 = NULL;
  const mxArray *c2_lhs20 = NULL;
  const mxArray *c2_rhs21 = NULL;
  const mxArray *c2_lhs21 = NULL;
  const mxArray *c2_rhs22 = NULL;
  const mxArray *c2_lhs22 = NULL;
  const mxArray *c2_rhs23 = NULL;
  const mxArray *c2_lhs23 = NULL;
  const mxArray *c2_rhs24 = NULL;
  const mxArray *c2_lhs24 = NULL;
  const mxArray *c2_rhs25 = NULL;
  const mxArray *c2_lhs25 = NULL;
  const mxArray *c2_rhs26 = NULL;
  const mxArray *c2_lhs26 = NULL;
  const mxArray *c2_rhs27 = NULL;
  const mxArray *c2_lhs27 = NULL;
  const mxArray *c2_rhs28 = NULL;
  const mxArray *c2_lhs28 = NULL;
  const mxArray *c2_rhs29 = NULL;
  const mxArray *c2_lhs29 = NULL;
  const mxArray *c2_rhs30 = NULL;
  const mxArray *c2_lhs30 = NULL;
  const mxArray *c2_rhs31 = NULL;
  const mxArray *c2_lhs31 = NULL;
  const mxArray *c2_rhs32 = NULL;
  const mxArray *c2_lhs32 = NULL;
  const mxArray *c2_rhs33 = NULL;
  const mxArray *c2_lhs33 = NULL;
  const mxArray *c2_rhs34 = NULL;
  const mxArray *c2_lhs34 = NULL;
  const mxArray *c2_rhs35 = NULL;
  const mxArray *c2_lhs35 = NULL;
  const mxArray *c2_rhs36 = NULL;
  const mxArray *c2_lhs36 = NULL;
  const mxArray *c2_rhs37 = NULL;
  const mxArray *c2_lhs37 = NULL;
  const mxArray *c2_rhs38 = NULL;
  const mxArray *c2_lhs38 = NULL;
  const mxArray *c2_rhs39 = NULL;
  const mxArray *c2_lhs39 = NULL;
  const mxArray *c2_rhs40 = NULL;
  const mxArray *c2_lhs40 = NULL;
  const mxArray *c2_rhs41 = NULL;
  const mxArray *c2_lhs41 = NULL;
  const mxArray *c2_rhs42 = NULL;
  const mxArray *c2_lhs42 = NULL;
  const mxArray *c2_rhs43 = NULL;
  const mxArray *c2_lhs43 = NULL;
  const mxArray *c2_rhs44 = NULL;
  const mxArray *c2_lhs44 = NULL;
  const mxArray *c2_rhs45 = NULL;
  const mxArray *c2_lhs45 = NULL;
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("mrdivide"), "name", "name", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1388424096U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1369981086U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c2_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c2_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("rdivide"), "name", "name", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363677880U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c2_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c2_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286786396U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c2_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_div"), "name", "name", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c2_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c2_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("cos"), "name", "name", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801572U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c2_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286786322U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c2_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("sin"), "name", "name", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801586U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c2_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286786336U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c2_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("mpower"), "name", "name", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363677878U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c2_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c2_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("ismatrix"), "name", "name", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1331268858U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c2_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("power"), "name", "name", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363677880U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c2_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c2_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c2_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c2_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c2_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c2_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("floor"), "name", "name", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363677854U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c2_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c2_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286786326U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c2_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 23);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 23);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c2_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 24);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("sqrt"), "name", "name", 24);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801586U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c2_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_error"), "name", "name",
                  25);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801558U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c2_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 26);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286786338U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c2_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 27);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 27);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1383841294U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c2_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 28);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 28);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c2_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 29);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 29);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c2_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 30);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 30);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c2_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 31);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  31);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375951890U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c2_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 32);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c2_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 33);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 33);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c2_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 34);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 34);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c2_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 35);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 35);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c2_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 36);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 36);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1381817900U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c2_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 37);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 37);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c2_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 38);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 38);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c2_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 39);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("atan"), "name", "name", 39);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801572U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c2_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_atan"), "name",
                  "name", 40);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286786318U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c2_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 41);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("asin"), "name", "name", 41);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/asin.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801570U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c2_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/asin.m"), "context",
                  "context", 42);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_error"), "name", "name",
                  42);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801558U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c2_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/asin.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_asin"), "name",
                  "name", 43);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_asin.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801576U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c2_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 44);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("tan"), "name", "name", 44);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tan.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1343801586U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c2_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tan.m"), "context",
                  "context", 45);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_tan"), "name",
                  "name", 45);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_tan.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286786338U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c2_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs45), "lhs", "lhs",
                  45);
  sf_mex_destroy(&c2_rhs0);
  sf_mex_destroy(&c2_lhs0);
  sf_mex_destroy(&c2_rhs1);
  sf_mex_destroy(&c2_lhs1);
  sf_mex_destroy(&c2_rhs2);
  sf_mex_destroy(&c2_lhs2);
  sf_mex_destroy(&c2_rhs3);
  sf_mex_destroy(&c2_lhs3);
  sf_mex_destroy(&c2_rhs4);
  sf_mex_destroy(&c2_lhs4);
  sf_mex_destroy(&c2_rhs5);
  sf_mex_destroy(&c2_lhs5);
  sf_mex_destroy(&c2_rhs6);
  sf_mex_destroy(&c2_lhs6);
  sf_mex_destroy(&c2_rhs7);
  sf_mex_destroy(&c2_lhs7);
  sf_mex_destroy(&c2_rhs8);
  sf_mex_destroy(&c2_lhs8);
  sf_mex_destroy(&c2_rhs9);
  sf_mex_destroy(&c2_lhs9);
  sf_mex_destroy(&c2_rhs10);
  sf_mex_destroy(&c2_lhs10);
  sf_mex_destroy(&c2_rhs11);
  sf_mex_destroy(&c2_lhs11);
  sf_mex_destroy(&c2_rhs12);
  sf_mex_destroy(&c2_lhs12);
  sf_mex_destroy(&c2_rhs13);
  sf_mex_destroy(&c2_lhs13);
  sf_mex_destroy(&c2_rhs14);
  sf_mex_destroy(&c2_lhs14);
  sf_mex_destroy(&c2_rhs15);
  sf_mex_destroy(&c2_lhs15);
  sf_mex_destroy(&c2_rhs16);
  sf_mex_destroy(&c2_lhs16);
  sf_mex_destroy(&c2_rhs17);
  sf_mex_destroy(&c2_lhs17);
  sf_mex_destroy(&c2_rhs18);
  sf_mex_destroy(&c2_lhs18);
  sf_mex_destroy(&c2_rhs19);
  sf_mex_destroy(&c2_lhs19);
  sf_mex_destroy(&c2_rhs20);
  sf_mex_destroy(&c2_lhs20);
  sf_mex_destroy(&c2_rhs21);
  sf_mex_destroy(&c2_lhs21);
  sf_mex_destroy(&c2_rhs22);
  sf_mex_destroy(&c2_lhs22);
  sf_mex_destroy(&c2_rhs23);
  sf_mex_destroy(&c2_lhs23);
  sf_mex_destroy(&c2_rhs24);
  sf_mex_destroy(&c2_lhs24);
  sf_mex_destroy(&c2_rhs25);
  sf_mex_destroy(&c2_lhs25);
  sf_mex_destroy(&c2_rhs26);
  sf_mex_destroy(&c2_lhs26);
  sf_mex_destroy(&c2_rhs27);
  sf_mex_destroy(&c2_lhs27);
  sf_mex_destroy(&c2_rhs28);
  sf_mex_destroy(&c2_lhs28);
  sf_mex_destroy(&c2_rhs29);
  sf_mex_destroy(&c2_lhs29);
  sf_mex_destroy(&c2_rhs30);
  sf_mex_destroy(&c2_lhs30);
  sf_mex_destroy(&c2_rhs31);
  sf_mex_destroy(&c2_lhs31);
  sf_mex_destroy(&c2_rhs32);
  sf_mex_destroy(&c2_lhs32);
  sf_mex_destroy(&c2_rhs33);
  sf_mex_destroy(&c2_lhs33);
  sf_mex_destroy(&c2_rhs34);
  sf_mex_destroy(&c2_lhs34);
  sf_mex_destroy(&c2_rhs35);
  sf_mex_destroy(&c2_lhs35);
  sf_mex_destroy(&c2_rhs36);
  sf_mex_destroy(&c2_lhs36);
  sf_mex_destroy(&c2_rhs37);
  sf_mex_destroy(&c2_lhs37);
  sf_mex_destroy(&c2_rhs38);
  sf_mex_destroy(&c2_lhs38);
  sf_mex_destroy(&c2_rhs39);
  sf_mex_destroy(&c2_lhs39);
  sf_mex_destroy(&c2_rhs40);
  sf_mex_destroy(&c2_lhs40);
  sf_mex_destroy(&c2_rhs41);
  sf_mex_destroy(&c2_lhs41);
  sf_mex_destroy(&c2_rhs42);
  sf_mex_destroy(&c2_lhs42);
  sf_mex_destroy(&c2_rhs43);
  sf_mex_destroy(&c2_lhs43);
  sf_mex_destroy(&c2_rhs44);
  sf_mex_destroy(&c2_lhs44);
  sf_mex_destroy(&c2_rhs45);
  sf_mex_destroy(&c2_lhs45);
}

static const mxArray *c2_emlrt_marshallOut(const char * c2_u)
{
  const mxArray *c2_y = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c2_u)), false);
  return c2_y;
}

static const mxArray *c2_b_emlrt_marshallOut(const uint32_T c2_u)
{
  const mxArray *c2_y = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 7, 0U, 0U, 0U, 0), false);
  return c2_y;
}

static void c2_eul2cbn(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_euler[3], real_T c2_Cbn[9])
{
  uint32_T c2_debug_family_var_map[10];
  real_T c2_cosph;
  real_T c2_sinph;
  real_T c2_costh;
  real_T c2_sinth;
  real_T c2_cosps;
  real_T c2_sinps;
  real_T c2_nargin = 1.0;
  real_T c2_nargout = 1.0;
  real_T c2_x;
  real_T c2_b_x;
  real_T c2_c_x;
  real_T c2_d_x;
  real_T c2_e_x;
  real_T c2_f_x;
  real_T c2_g_x;
  real_T c2_h_x;
  real_T c2_i_x;
  real_T c2_j_x;
  real_T c2_k_x;
  real_T c2_l_x;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 10U, 10U, c2_b_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_cosph, 0U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_sinph, 1U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_costh, 2U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_sinth, 3U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_cosps, 4U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_sinps, 5U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 6U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 7U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_euler, 8U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Cbn, 9U, c2_e_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_EML_FCN(0, 1);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 205U);
  c2_x = c2_euler[0];
  c2_cosph = c2_x;
  c2_b_x = c2_cosph;
  c2_cosph = c2_b_x;
  c2_cosph = muDoubleScalarCos(c2_cosph);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 206U);
  c2_c_x = c2_euler[0];
  c2_sinph = c2_c_x;
  c2_d_x = c2_sinph;
  c2_sinph = c2_d_x;
  c2_sinph = muDoubleScalarSin(c2_sinph);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 207U);
  c2_e_x = c2_euler[1];
  c2_costh = c2_e_x;
  c2_f_x = c2_costh;
  c2_costh = c2_f_x;
  c2_costh = muDoubleScalarCos(c2_costh);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 208U);
  c2_g_x = c2_euler[1];
  c2_sinth = c2_g_x;
  c2_h_x = c2_sinth;
  c2_sinth = c2_h_x;
  c2_sinth = muDoubleScalarSin(c2_sinth);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 209U);
  c2_i_x = c2_euler[2];
  c2_cosps = c2_i_x;
  c2_j_x = c2_cosps;
  c2_cosps = c2_j_x;
  c2_cosps = muDoubleScalarCos(c2_cosps);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 210U);
  c2_k_x = c2_euler[2];
  c2_sinps = c2_k_x;
  c2_l_x = c2_sinps;
  c2_sinps = c2_l_x;
  c2_sinps = muDoubleScalarSin(c2_sinps);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 212U);
  c2_Cbn[0] = c2_costh * c2_cosps;
  c2_Cbn[3] = -c2_cosph * c2_sinps + c2_sinph * c2_sinth * c2_cosps;
  c2_Cbn[6] = c2_sinph * c2_sinps + c2_cosph * c2_sinth * c2_cosps;
  c2_Cbn[1] = c2_costh * c2_sinps;
  c2_Cbn[4] = c2_cosph * c2_cosps + c2_sinph * c2_sinth * c2_sinps;
  c2_Cbn[7] = -c2_sinph * c2_cosps + c2_cosph * c2_sinth * c2_sinps;
  c2_Cbn[2] = -c2_sinth;
  c2_Cbn[5] = c2_sinph * c2_costh;
  c2_Cbn[8] = c2_cosph * c2_costh;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, -212);
  _SFD_SYMBOL_SCOPE_POP();
}

static real_T c2_mpower(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_a)
{
  real_T c2_b_a;
  real_T c2_c_a;
  real_T c2_ak;
  real_T c2_d_a;
  c2_b_a = c2_a;
  c2_c_a = c2_b_a;
  c2_eml_scalar_eg(chartInstance);
  c2_ak = c2_c_a;
  c2_d_a = c2_ak;
  c2_eml_scalar_eg(chartInstance);
  return c2_d_a * c2_d_a;
}

static void c2_eml_scalar_eg(SFc2_aircraftControl_smallerStateEKFInstanceStruct *
  chartInstance)
{
  (void)chartInstance;
}

static real_T c2_sqrt(SFc2_aircraftControl_smallerStateEKFInstanceStruct
                      *chartInstance, real_T c2_x)
{
  real_T c2_b_x;
  c2_b_x = c2_x;
  c2_b_sqrt(chartInstance, &c2_b_x);
  return c2_b_x;
}

static void c2_eml_error(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance)
{
  int32_T c2_i59;
  static char_T c2_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c2_u[30];
  const mxArray *c2_y = NULL;
  int32_T c2_i60;
  static char_T c2_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c2_b_u[4];
  const mxArray *c2_b_y = NULL;
  (void)chartInstance;
  for (c2_i59 = 0; c2_i59 < 30; c2_i59++) {
    c2_u[c2_i59] = c2_cv0[c2_i59];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c2_i60 = 0; c2_i60 < 4; c2_i60++) {
    c2_b_u[c2_i60] = c2_cv1[c2_i60];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c2_y, 14, c2_b_y));
}

static void c2_b_eml_scalar_eg
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c2_eml_xgemm(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_A[9], real_T c2_B[3], real_T c2_C[3], real_T c2_b_C
  [3])
{
  int32_T c2_i61;
  int32_T c2_i62;
  real_T c2_b_A[9];
  int32_T c2_i63;
  real_T c2_b_B[3];
  for (c2_i61 = 0; c2_i61 < 3; c2_i61++) {
    c2_b_C[c2_i61] = c2_C[c2_i61];
  }

  for (c2_i62 = 0; c2_i62 < 9; c2_i62++) {
    c2_b_A[c2_i62] = c2_A[c2_i62];
  }

  for (c2_i63 = 0; c2_i63 < 3; c2_i63++) {
    c2_b_B[c2_i63] = c2_B[c2_i63];
  }

  c2_b_eml_xgemm(chartInstance, c2_b_A, c2_b_B, c2_b_C);
}

static real_T c2_asin(SFc2_aircraftControl_smallerStateEKFInstanceStruct
                      *chartInstance, real_T c2_x)
{
  real_T c2_b_x;
  c2_b_x = c2_x;
  c2_b_asin(chartInstance, &c2_b_x);
  return c2_b_x;
}

static void c2_b_eml_error(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance)
{
  int32_T c2_i64;
  static char_T c2_cv2[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c2_u[30];
  const mxArray *c2_y = NULL;
  int32_T c2_i65;
  static char_T c2_cv3[4] = { 'a', 's', 'i', 'n' };

  char_T c2_b_u[4];
  const mxArray *c2_b_y = NULL;
  (void)chartInstance;
  for (c2_i64 = 0; c2_i64 < 30; c2_i64++) {
    c2_u[c2_i64] = c2_cv2[c2_i64];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  for (c2_i65 = 0; c2_i65 < 4; c2_i65++) {
    c2_b_u[c2_i65] = c2_cv3[c2_i65];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c2_y, 14, c2_b_y));
}

static real_T c2_b_mpower(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_a)
{
  real_T c2_b_a;
  real_T c2_c_a;
  real_T c2_ak;
  real_T c2_d_a;
  real_T c2_ar;
  c2_b_a = c2_a;
  c2_c_a = c2_b_a;
  c2_eml_scalar_eg(chartInstance);
  c2_ak = c2_c_a;
  c2_d_a = c2_ak;
  c2_eml_scalar_eg(chartInstance);
  c2_ar = c2_d_a;
  return muDoubleScalarPower(c2_ar, 3.0);
}

static real_T c2_c_mpower(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_a)
{
  real_T c2_b_a;
  real_T c2_c_a;
  real_T c2_ak;
  real_T c2_d_a;
  real_T c2_ar;
  c2_b_a = c2_a;
  c2_c_a = c2_b_a;
  c2_eml_scalar_eg(chartInstance);
  c2_ak = c2_c_a;
  c2_d_a = c2_ak;
  c2_eml_scalar_eg(chartInstance);
  c2_ar = c2_d_a;
  return muDoubleScalarPower(c2_ar, 4.0);
}

static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(int32_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, false);
  return c2_mxArrayOutData;
}

static int32_T c2_h_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  int32_T c2_y;
  int32_T c2_i66;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_i66, 1, 6, 0U, 0, 0U, 0);
  c2_y = c2_i66;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_sfEvent;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  int32_T c2_y;
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
    chartInstanceVoid;
  c2_b_sfEvent = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_sfEvent),
    &c2_thisId);
  sf_mex_destroy(&c2_b_sfEvent);
  *(int32_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static uint8_T c2_i_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_b_is_active_c2_aircraftControl_smallerStateEKF, const char_T
   *c2_identifier)
{
  uint8_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_j_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c2_b_is_active_c2_aircraftControl_smallerStateEKF), &c2_thisId);
  sf_mex_destroy(&c2_b_is_active_c2_aircraftControl_smallerStateEKF);
  return c2_y;
}

static uint8_T c2_j_emlrt_marshallIn
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance, const
   mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_y;
  uint8_T c2_u0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_y = c2_u0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_b_sqrt(SFc2_aircraftControl_smallerStateEKFInstanceStruct
                      *chartInstance, real_T *c2_x)
{
  if (*c2_x < 0.0) {
    c2_eml_error(chartInstance);
  }

  *c2_x = muDoubleScalarSqrt(*c2_x);
}

static void c2_b_eml_xgemm(SFc2_aircraftControl_smallerStateEKFInstanceStruct
  *chartInstance, real_T c2_A[9], real_T c2_B[3], real_T c2_C[3])
{
  int32_T c2_i67;
  int32_T c2_i68;
  int32_T c2_i69;
  (void)chartInstance;
  for (c2_i67 = 0; c2_i67 < 3; c2_i67++) {
    c2_C[c2_i67] = 0.0;
    c2_i68 = 0;
    for (c2_i69 = 0; c2_i69 < 3; c2_i69++) {
      c2_C[c2_i67] += c2_A[c2_i68 + c2_i67] * c2_B[c2_i69];
      c2_i68 += 3;
    }
  }
}

static void c2_b_asin(SFc2_aircraftControl_smallerStateEKFInstanceStruct
                      *chartInstance, real_T *c2_x)
{
  boolean_T guard1 = false;
  guard1 = false;
  if (*c2_x < -1.0) {
    guard1 = true;
  } else {
    if (1.0 < *c2_x) {
      guard1 = true;
    }
  }

  if (guard1 == true) {
    c2_b_eml_error(chartInstance);
  }

  *c2_x = muDoubleScalarAsin(*c2_x);
}

static void init_dsm_address_info
  (SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c2_aircraftControl_smallerStateEKF_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2957806761U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1737119201U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2985843587U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4129876773U);
}

mxArray *sf_c2_aircraftControl_smallerStateEKF_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("CID6rmT4eMQ7PIXwQUuvDD");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c2_aircraftControl_smallerStateEKF_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c2_aircraftControl_smallerStateEKF_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c2_aircraftControl_smallerStateEKF
  (void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x5'type','srcId','name','auxInfo'{{M[1],M[13],T\"ALP\",},{M[1],M[14],T\"BET\",},{M[1],M[10],T\"abody\",},{M[1],M[5],T\"yDots\",},{M[8],M[0],T\"is_active_c2_aircraftControl_smallerStateEKF\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 5, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_aircraftControl_smallerStateEKF_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _aircraftControl_smallerStateEKFMachineNumber_,
           2,
           1,
           1,
           0,
           7,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation
          (_aircraftControl_smallerStateEKFMachineNumber_,
           chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,
             _aircraftControl_smallerStateEKFMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _aircraftControl_smallerStateEKFMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"y");
          _SFD_SET_DATA_PROPS(1,1,1,0,"u");
          _SFD_SET_DATA_PROPS(2,2,0,1,"yDots");
          _SFD_SET_DATA_PROPS(3,2,0,1,"abody");
          _SFD_SET_DATA_PROPS(4,1,1,0,"loopFlag");
          _SFD_SET_DATA_PROPS(5,2,0,1,"ALP");
          _SFD_SET_DATA_PROPS(6,2,0,1,"BET");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,2,6,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,5696);
        _SFD_CV_INIT_EML_FCN(0,1,"eul2cbn",5696,-1,6329);
        _SFD_CV_INIT_EML_IF(0,1,0,1628,1647,-1,1702);
        _SFD_CV_INIT_EML_IF(0,1,1,1704,1723,-1,1790);
        _SFD_CV_INIT_EML_IF(0,1,2,1793,1812,-1,1861);
        _SFD_CV_INIT_EML_IF(0,1,3,5159,5178,-1,5233);
        _SFD_CV_INIT_EML_IF(0,1,4,5235,5254,-1,5321);
        _SFD_CV_INIT_EML_IF(0,1,5,5324,5343,-1,5392);

        {
          unsigned int dimVector[1];
          dimVector[0]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)
            c2_c_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)
            c2_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);

        {
          real_T *c2_ALP;
          real_T *c2_BET;
          real_T (*c2_y)[12];
          real_T (*c2_u)[4];
          real_T (*c2_yDots)[12];
          real_T (*c2_abody)[3];
          real_T (*c2_loopFlag)[3];
          c2_BET = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
          c2_ALP = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
          c2_loopFlag = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
          c2_abody = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
          c2_yDots = (real_T (*)[12])ssGetOutputPortSignal(chartInstance->S, 1);
          c2_u = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 1);
          c2_y = (real_T (*)[12])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c2_y);
          _SFD_SET_DATA_VALUE_PTR(1U, *c2_u);
          _SFD_SET_DATA_VALUE_PTR(2U, *c2_yDots);
          _SFD_SET_DATA_VALUE_PTR(3U, *c2_abody);
          _SFD_SET_DATA_VALUE_PTR(4U, *c2_loopFlag);
          _SFD_SET_DATA_VALUE_PTR(5U, c2_ALP);
          _SFD_SET_DATA_VALUE_PTR(6U, c2_BET);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _aircraftControl_smallerStateEKFMachineNumber_,
        chartInstance->chartNumber,chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "bGrY11S7Y0iWG6vIOmiGtC";
}

static void sf_opaque_initialize_c2_aircraftControl_smallerStateEKF(void
  *chartInstanceVar)
{
  chart_debug_initialization
    (((SFc2_aircraftControl_smallerStateEKFInstanceStruct*) chartInstanceVar)->S,
     0);
  initialize_params_c2_aircraftControl_smallerStateEKF
    ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*) chartInstanceVar);
  initialize_c2_aircraftControl_smallerStateEKF
    ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c2_aircraftControl_smallerStateEKF(void
  *chartInstanceVar)
{
  enable_c2_aircraftControl_smallerStateEKF
    ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c2_aircraftControl_smallerStateEKF(void
  *chartInstanceVar)
{
  disable_c2_aircraftControl_smallerStateEKF
    ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c2_aircraftControl_smallerStateEKF(void
  *chartInstanceVar)
{
  sf_gateway_c2_aircraftControl_smallerStateEKF
    ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*) chartInstanceVar);
}

extern const mxArray*
  sf_internal_get_sim_state_c2_aircraftControl_smallerStateEKF(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c2_aircraftControl_smallerStateEKF
    ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*)
     chartInfo->chartInstance);        /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_aircraftControl_smallerStateEKF();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c2_aircraftControl_smallerStateEKF
  (SimStruct* S, const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c2_aircraftControl_smallerStateEKF();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c2_aircraftControl_smallerStateEKF
    ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*)
     chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c2_aircraftControl_smallerStateEKF
  (SimStruct* S)
{
  return sf_internal_get_sim_state_c2_aircraftControl_smallerStateEKF(S);
}

static void sf_opaque_set_sim_state_c2_aircraftControl_smallerStateEKF(SimStruct*
  S, const mxArray *st)
{
  sf_internal_set_sim_state_c2_aircraftControl_smallerStateEKF(S, st);
}

static void sf_opaque_terminate_c2_aircraftControl_smallerStateEKF(void
  *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*)
                    chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_aircraftControl_smallerStateEKF_optimization_info();
    }

    finalize_c2_aircraftControl_smallerStateEKF
      ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc2_aircraftControl_smallerStateEKF
    ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_aircraftControl_smallerStateEKF(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c2_aircraftControl_smallerStateEKF
      ((SFc2_aircraftControl_smallerStateEKFInstanceStruct*)
       (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c2_aircraftControl_smallerStateEKF(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_aircraftControl_smallerStateEKF_optimization_info
      ();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,2,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,2);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,2,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,2,4);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=4; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1230447615U));
  ssSetChecksum1(S,(2814533331U));
  ssSetChecksum2(S,(2929743696U));
  ssSetChecksum3(S,(81889865U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c2_aircraftControl_smallerStateEKF(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_aircraftControl_smallerStateEKF(SimStruct *S)
{
  SFc2_aircraftControl_smallerStateEKFInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc2_aircraftControl_smallerStateEKFInstanceStruct *)utMalloc
    (sizeof(SFc2_aircraftControl_smallerStateEKFInstanceStruct));
  memset(chartInstance, 0, sizeof
         (SFc2_aircraftControl_smallerStateEKFInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.mdlStart =
    mdlStart_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c2_aircraftControl_smallerStateEKF;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c2_aircraftControl_smallerStateEKF_method_dispatcher(SimStruct *S, int_T
  method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_aircraftControl_smallerStateEKF(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_aircraftControl_smallerStateEKF(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_aircraftControl_smallerStateEKF(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_aircraftControl_smallerStateEKF_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
