/* Include files */

#include <stddef.h>
#include "blas.h"
#include "aircraftControl_FullStateFilters_sfun.h"
#include "c3_aircraftControl_FullStateFilters.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "aircraftControl_FullStateFilters_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c3_debug_family_names[13] = { "nargin", "nargout",
  "accelerations", "omega", "ALP", "BET", "u", "Vt", "H", "states_OUT",
  "cov_OUT", "innovation", "S" };

static const char * c3_b_debug_family_names[130] = { "Qcov", "Rcov", "R2D", "g",
  "vs", "rho", "dt", "m", "Ix", "Iy", "Iz", "Ixz", "b", "cbar", "S", "xCG_ref",
  "xCG", "dA", "dE", "dR", "throttle", "measurement", "qbar", "hT", "Tmax", "T",
  "kappa", "na", "N", "sigma", "sigma_meas", "W", "matrixSQRT", "i", "p", "q",
  "r", "prev_pdot", "prev_qdot", "prev_rdot", "CX_dE1_mp", "CX_dE2_mp",
  "CY_dE1_mp", "CY_dR1_mp", "CY_dR2_mp", "CZ_dE1_mp", "CZ_dE2_mp", "CZ_dA1_mp",
  "Cl_dA1_mp", "Cl_dA2_mp", "Cl_dA3_mp", "Cl_dR1_mp", "Cl_dR2_mp", "Cl_dE1_mp",
  "Cm_dE1_mp", "Cm_dE2_mp", "Cm_dE3_mp", "Cm_dA1_mp", "Cn_dA1_mp", "Cn_dA2_mp",
  "Cn_dE1_mp", "Cn_dE2_mp", "Cn_dR1_mp", "Cn_dR2_mp", "CX_dE1", "CX_dE2",
  "CY_dE1", "CY_dR1", "CY_dR2", "CZ_dE1", "CZ_dE2", "CZ_dA1", "Cl_dA1", "Cl_dA2",
  "Cl_dA3", "Cl_dR1", "Cl_dR2", "Cl_dE1", "Cm_dE1", "Cm_dE2", "Cm_dE3", "Cm_dA1",
  "Cn_dA1", "Cn_dA2", "Cn_dE1", "Cn_dE2", "Cn_dR1", "Cn_dR2", "CX", "CY", "CZ",
  "Cl", "Cm", "Cn", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "heng",
  "pdot", "qdot", "rdot", "ymeas", "Pxz", "Pzz", "Scov", "k", "correction",
  "nargin", "nargout", "u", "Vt", "ALP", "BET", "accelerations", "omega",
  "Height", "pStatesOUT", "pCovarianceOUT", "innovation", "SCovOUT", "pStates",
  "pCovariance", "covDiag", "prevOmegaDots" };

/* Function Declarations */
static void initialize_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void initialize_params_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void enable_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void disable_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_update_debugger_state_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void set_sim_state_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_st);
static void finalize_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void sf_gateway_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void initSimStructsc3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_simulateUKF(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_u[4], real_T c3_Vt, real_T c3_ALP, real_T c3_BET,
  real_T c3_accelerations[3], real_T c3_omega[3], real_T c3_Height, real_T
  c3_pStatesOUT[30], real_T c3_pCovarianceOUT[30], real_T c3_innovation[6],
  real_T c3_SCovOUT[6]);
static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber);
static void c3_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_sprintf, const char_T *c3_identifier, char_T c3_y[14]);
static void c3_b_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, char_T c3_y[14]);
static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData);
static void c3_c_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_S, const char_T *c3_identifier, real_T c3_y[6]);
static void c3_d_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[6]);
static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_e_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_cov_OUT, const char_T *c3_identifier, real_T c3_y[30]);
static void c3_f_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[30]);
static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static real_T c3_g_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_f_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_h_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_prevOmegaDots, const char_T *c3_identifier, real_T c3_y[183]);
static void c3_i_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[183]);
static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_g_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_j_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_covDiag, const char_T *c3_identifier, real_T c3_y[30]);
static void c3_k_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[30]);
static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_h_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_l_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_pCovariance, const char_T *c3_identifier, real_T c3_y[900]);
static void c3_m_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[900]);
static void c3_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_i_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_n_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_pStates, const char_T *c3_identifier, real_T c3_y[30]);
static void c3_o_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[30]);
static void c3_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_p_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[3]);
static void c3_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_q_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[4]);
static void c3_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_j_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_r_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[180]);
static void c3_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_k_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_s_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[36]);
static void c3_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_l_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_t_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[900]);
static void c3_l_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_m_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_u_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[61]);
static void c3_m_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_n_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_v_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[366]);
static void c3_n_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_o_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_w_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[1830]);
static void c3_o_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_info_helper(const mxArray **c3_info);
static const mxArray *c3_emlrt_marshallOut(const char * c3_u);
static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u);
static void c3_b_info_helper(const mxArray **c3_info);
static void c3_c_info_helper(const mxArray **c3_info);
static void c3_d_info_helper(const mxArray **c3_info);
static real_T c3_mpower(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a);
static real_T c3_power(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a);
static void c3_eml_scalar_eg(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance);
static void c3_diag(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                    *chartInstance, real_T c3_v[30], real_T c3_d[900]);
static void c3_eml_switch_helper
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_b_eml_switch_helper
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_b_power(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a[900], real_T c3_y[900]);
static real_T c3_b_mpower(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a);
static real_T c3_c_mpower(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a);
static void c3_eml_error(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance);
static void c3_eml_matlab_zpotrf
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, real_T
   c3_A[900], real_T c3_b_A[900], int32_T *c3_info);
static real_T c3_eml_xdotc(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_n, real_T c3_x[900], int32_T c3_ix0, real_T c3_y
  [900], int32_T c3_iy0);
static void c3_check_forloop_overflow_error
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, boolean_T
   c3_overflow);
static void c3_eml_xgemv(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_m, int32_T c3_n, int32_T c3_ia0, int32_T c3_ix0,
  real_T c3_y[900], int32_T c3_iy0, real_T c3_b_y[900]);
static void c3_below_threshold
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_b_below_threshold
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_b_eml_error(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance);
static void c3_inv(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                   *chartInstance, real_T c3_x[36], real_T c3_y[36]);
static void c3_eps(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                   *chartInstance);
static void c3_eml_matlab_zgetrf
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, real_T
   c3_A[36], real_T c3_b_A[36], int32_T c3_ipiv[6], int32_T *c3_info);
static int32_T c3_eml_ixamax(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_n, real_T c3_x[36], int32_T c3_ix0);
static void c3_threshold(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance);
static void c3_eml_xgeru(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_m, int32_T c3_n, real_T c3_alpha1, int32_T c3_ix0,
  int32_T c3_iy0, real_T c3_A[36], int32_T c3_ia0, real_T c3_b_A[36]);
static void c3_eml_ipiv2perm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_ipiv[6], int32_T c3_perm[6]);
static void c3_eml_xtrsm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[36], real_T c3_B[36], real_T c3_b_B[36]);
static void c3_b_threshold(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance);
static real_T c3_norm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                      *chartInstance, real_T c3_x[36]);
static void c3_eml_warning(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance);
static void c3_b_eml_warning(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, char_T c3_varargin_2[14]);
static void c3_b_eml_scalar_eg
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_eml_xgemm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[180], real_T c3_B[36], real_T c3_C[180], real_T
  c3_b_C[180]);
static void c3_c_threshold(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance);
static void c3_c_eml_scalar_eg
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_d_eml_scalar_eg
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);
static void c3_b_eml_xgemm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[180], real_T c3_B[180], real_T c3_C[900], real_T
  c3_b_C[900]);
static void c3_b_diag(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                      *chartInstance, real_T c3_v[36], real_T c3_d[6]);
static void c3_c_diag(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                      *chartInstance, real_T c3_v[900], real_T c3_d[30]);
static const mxArray *c3_p_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static int32_T c3_x_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_p_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static uint8_T c3_y_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_is_active_c3_aircraftControl_FullStateFilters, const char_T
   *c3_identifier);
static uint8_T c3_ab_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static int32_T c3_b_eml_matlab_zpotrf
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, real_T
   c3_A[900]);
static void c3_b_eml_xgemv(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_m, int32_T c3_n, int32_T c3_ia0, int32_T c3_ix0,
  real_T c3_y[900], int32_T c3_iy0);
static void c3_b_eml_matlab_zgetrf
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, real_T
   c3_A[36], int32_T c3_ipiv[6], int32_T *c3_info);
static void c3_b_eml_xgeru(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_m, int32_T c3_n, real_T c3_alpha1, int32_T c3_ix0,
  int32_T c3_iy0, real_T c3_A[36], int32_T c3_ia0);
static void c3_b_eml_xtrsm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[36], real_T c3_B[36]);
static void c3_c_eml_xgemm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[180], real_T c3_B[36], real_T c3_C[180]);
static void c3_d_eml_xgemm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[180], real_T c3_B[180], real_T c3_C[900]);
static void init_dsm_address_info
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  chartInstance->c3_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c3_pStates_not_empty = false;
  chartInstance->c3_pCovariance_not_empty = false;
  chartInstance->c3_covDiag_not_empty = false;
  chartInstance->c3_prevOmegaDots_not_empty = false;
  chartInstance->c3_is_active_c3_aircraftControl_FullStateFilters = 0U;
}

static void initialize_params_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c3_update_debugger_state_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  const mxArray *c3_st;
  const mxArray *c3_y = NULL;
  int32_T c3_i0;
  real_T c3_u[6];
  const mxArray *c3_b_y = NULL;
  int32_T c3_i1;
  real_T c3_b_u[30];
  const mxArray *c3_c_y = NULL;
  int32_T c3_i2;
  real_T c3_c_u[6];
  const mxArray *c3_d_y = NULL;
  int32_T c3_i3;
  real_T c3_d_u[30];
  const mxArray *c3_e_y = NULL;
  int32_T c3_i4;
  real_T c3_e_u[30];
  const mxArray *c3_f_y = NULL;
  int32_T c3_i5;
  real_T c3_f_u[900];
  const mxArray *c3_g_y = NULL;
  int32_T c3_i6;
  real_T c3_g_u[30];
  const mxArray *c3_h_y = NULL;
  int32_T c3_i7;
  real_T c3_h_u[183];
  const mxArray *c3_i_y = NULL;
  uint8_T c3_hoistedGlobal;
  uint8_T c3_i_u;
  const mxArray *c3_j_y = NULL;
  real_T (*c3_states_OUT)[30];
  real_T (*c3_innovation)[6];
  real_T (*c3_cov_OUT)[30];
  real_T (*c3_S)[6];
  c3_S = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 4);
  c3_innovation = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 3);
  c3_cov_OUT = (real_T (*)[30])ssGetOutputPortSignal(chartInstance->S, 2);
  c3_states_OUT = (real_T (*)[30])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_st = NULL;
  c3_st = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createcellmatrix(9, 1), false);
  for (c3_i0 = 0; c3_i0 < 6; c3_i0++) {
    c3_u[c3_i0] = (*c3_S)[c3_i0];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_setcell(c3_y, 0, c3_b_y);
  for (c3_i1 = 0; c3_i1 < 30; c3_i1++) {
    c3_b_u[c3_i1] = (*c3_cov_OUT)[c3_i1];
  }

  c3_c_y = NULL;
  sf_mex_assign(&c3_c_y, sf_mex_create("y", c3_b_u, 0, 0U, 1U, 0U, 1, 30), false);
  sf_mex_setcell(c3_y, 1, c3_c_y);
  for (c3_i2 = 0; c3_i2 < 6; c3_i2++) {
    c3_c_u[c3_i2] = (*c3_innovation)[c3_i2];
  }

  c3_d_y = NULL;
  sf_mex_assign(&c3_d_y, sf_mex_create("y", c3_c_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_setcell(c3_y, 2, c3_d_y);
  for (c3_i3 = 0; c3_i3 < 30; c3_i3++) {
    c3_d_u[c3_i3] = (*c3_states_OUT)[c3_i3];
  }

  c3_e_y = NULL;
  sf_mex_assign(&c3_e_y, sf_mex_create("y", c3_d_u, 0, 0U, 1U, 0U, 1, 30), false);
  sf_mex_setcell(c3_y, 3, c3_e_y);
  for (c3_i4 = 0; c3_i4 < 30; c3_i4++) {
    c3_e_u[c3_i4] = chartInstance->c3_covDiag[c3_i4];
  }

  c3_f_y = NULL;
  if (!chartInstance->c3_covDiag_not_empty) {
    sf_mex_assign(&c3_f_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c3_f_y, sf_mex_create("y", c3_e_u, 0, 0U, 1U, 0U, 1, 30),
                  false);
  }

  sf_mex_setcell(c3_y, 4, c3_f_y);
  for (c3_i5 = 0; c3_i5 < 900; c3_i5++) {
    c3_f_u[c3_i5] = chartInstance->c3_pCovariance[c3_i5];
  }

  c3_g_y = NULL;
  if (!chartInstance->c3_pCovariance_not_empty) {
    sf_mex_assign(&c3_g_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c3_g_y, sf_mex_create("y", c3_f_u, 0, 0U, 1U, 0U, 2, 30, 30),
                  false);
  }

  sf_mex_setcell(c3_y, 5, c3_g_y);
  for (c3_i6 = 0; c3_i6 < 30; c3_i6++) {
    c3_g_u[c3_i6] = chartInstance->c3_pStates[c3_i6];
  }

  c3_h_y = NULL;
  if (!chartInstance->c3_pStates_not_empty) {
    sf_mex_assign(&c3_h_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c3_h_y, sf_mex_create("y", c3_g_u, 0, 0U, 1U, 0U, 1, 30),
                  false);
  }

  sf_mex_setcell(c3_y, 6, c3_h_y);
  for (c3_i7 = 0; c3_i7 < 183; c3_i7++) {
    c3_h_u[c3_i7] = chartInstance->c3_prevOmegaDots[c3_i7];
  }

  c3_i_y = NULL;
  if (!chartInstance->c3_prevOmegaDots_not_empty) {
    sf_mex_assign(&c3_i_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c3_i_y, sf_mex_create("y", c3_h_u, 0, 0U, 1U, 0U, 2, 3, 61),
                  false);
  }

  sf_mex_setcell(c3_y, 7, c3_i_y);
  c3_hoistedGlobal =
    chartInstance->c3_is_active_c3_aircraftControl_FullStateFilters;
  c3_i_u = c3_hoistedGlobal;
  c3_j_y = NULL;
  sf_mex_assign(&c3_j_y, sf_mex_create("y", &c3_i_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c3_y, 8, c3_j_y);
  sf_mex_assign(&c3_st, c3_y, false);
  return c3_st;
}

static void set_sim_state_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_st)
{
  const mxArray *c3_u;
  real_T c3_dv0[6];
  int32_T c3_i8;
  real_T c3_dv1[30];
  int32_T c3_i9;
  real_T c3_dv2[6];
  int32_T c3_i10;
  real_T c3_dv3[30];
  int32_T c3_i11;
  real_T c3_dv4[30];
  int32_T c3_i12;
  real_T c3_dv5[900];
  int32_T c3_i13;
  real_T c3_dv6[30];
  int32_T c3_i14;
  real_T c3_dv7[183];
  int32_T c3_i15;
  real_T (*c3_S)[6];
  real_T (*c3_cov_OUT)[30];
  real_T (*c3_innovation)[6];
  real_T (*c3_states_OUT)[30];
  c3_S = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 4);
  c3_innovation = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 3);
  c3_cov_OUT = (real_T (*)[30])ssGetOutputPortSignal(chartInstance->S, 2);
  c3_states_OUT = (real_T (*)[30])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c3_doneDoubleBufferReInit = true;
  c3_u = sf_mex_dup(c3_st);
  c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 0)), "S",
                        c3_dv0);
  for (c3_i8 = 0; c3_i8 < 6; c3_i8++) {
    (*c3_S)[c3_i8] = c3_dv0[c3_i8];
  }

  c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 1)),
                        "cov_OUT", c3_dv1);
  for (c3_i9 = 0; c3_i9 < 30; c3_i9++) {
    (*c3_cov_OUT)[c3_i9] = c3_dv1[c3_i9];
  }

  c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 2)),
                        "innovation", c3_dv2);
  for (c3_i10 = 0; c3_i10 < 6; c3_i10++) {
    (*c3_innovation)[c3_i10] = c3_dv2[c3_i10];
  }

  c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 3)),
                        "states_OUT", c3_dv3);
  for (c3_i11 = 0; c3_i11 < 30; c3_i11++) {
    (*c3_states_OUT)[c3_i11] = c3_dv3[c3_i11];
  }

  c3_j_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 4)),
                        "covDiag", c3_dv4);
  for (c3_i12 = 0; c3_i12 < 30; c3_i12++) {
    chartInstance->c3_covDiag[c3_i12] = c3_dv4[c3_i12];
  }

  c3_l_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 5)),
                        "pCovariance", c3_dv5);
  for (c3_i13 = 0; c3_i13 < 900; c3_i13++) {
    chartInstance->c3_pCovariance[c3_i13] = c3_dv5[c3_i13];
  }

  c3_n_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 6)),
                        "pStates", c3_dv6);
  for (c3_i14 = 0; c3_i14 < 30; c3_i14++) {
    chartInstance->c3_pStates[c3_i14] = c3_dv6[c3_i14];
  }

  c3_h_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 7)),
                        "prevOmegaDots", c3_dv7);
  for (c3_i15 = 0; c3_i15 < 183; c3_i15++) {
    chartInstance->c3_prevOmegaDots[c3_i15] = c3_dv7[c3_i15];
  }

  chartInstance->c3_is_active_c3_aircraftControl_FullStateFilters =
    c3_y_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 8)),
    "is_active_c3_aircraftControl_FullStateFilters");
  sf_mex_destroy(&c3_u);
  c3_update_debugger_state_c3_aircraftControl_FullStateFilters(chartInstance);
  sf_mex_destroy(&c3_st);
}

static void finalize_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  int32_T c3_i16;
  int32_T c3_i17;
  int32_T c3_i18;
  real_T c3_hoistedGlobal;
  real_T c3_b_hoistedGlobal;
  real_T c3_c_hoistedGlobal;
  real_T c3_d_hoistedGlobal;
  int32_T c3_i19;
  real_T c3_accelerations[3];
  int32_T c3_i20;
  real_T c3_omega[3];
  real_T c3_ALP;
  real_T c3_BET;
  int32_T c3_i21;
  real_T c3_u[4];
  real_T c3_Vt;
  real_T c3_H;
  uint32_T c3_debug_family_var_map[13];
  real_T c3_nargin = 7.0;
  real_T c3_nargout = 4.0;
  real_T c3_states_OUT[30];
  real_T c3_cov_OUT[30];
  real_T c3_innovation[6];
  real_T c3_S[6];
  int32_T c3_i22;
  real_T c3_b_u[4];
  int32_T c3_i23;
  real_T c3_b_accelerations[3];
  int32_T c3_i24;
  real_T c3_b_omega[3];
  real_T c3_b_S[6];
  real_T c3_b_innovation[6];
  real_T c3_b_cov_OUT[30];
  real_T c3_b_states_OUT[30];
  int32_T c3_i25;
  int32_T c3_i26;
  int32_T c3_i27;
  int32_T c3_i28;
  int32_T c3_i29;
  int32_T c3_i30;
  int32_T c3_i31;
  int32_T c3_i32;
  int32_T c3_i33;
  int32_T c3_i34;
  int32_T c3_i35;
  int32_T c3_i36;
  real_T *c3_b_ALP;
  real_T *c3_b_BET;
  real_T *c3_b_Vt;
  real_T *c3_b_H;
  real_T (*c3_c_states_OUT)[30];
  real_T (*c3_c_cov_OUT)[30];
  real_T (*c3_c_innovation)[6];
  real_T (*c3_c_S)[6];
  real_T (*c3_c_u)[4];
  real_T (*c3_c_omega)[3];
  real_T (*c3_c_accelerations)[3];
  c3_c_S = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 4);
  c3_c_innovation = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 3);
  c3_b_H = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
  c3_c_cov_OUT = (real_T (*)[30])ssGetOutputPortSignal(chartInstance->S, 2);
  c3_b_Vt = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c3_c_states_OUT = (real_T (*)[30])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_c_u = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 4);
  c3_b_BET = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c3_b_ALP = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c3_c_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c3_c_accelerations = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  for (c3_i16 = 0; c3_i16 < 3; c3_i16++) {
    _SFD_DATA_RANGE_CHECK((*c3_c_accelerations)[c3_i16], 0U);
  }

  for (c3_i17 = 0; c3_i17 < 3; c3_i17++) {
    _SFD_DATA_RANGE_CHECK((*c3_c_omega)[c3_i17], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c3_b_ALP, 2U);
  _SFD_DATA_RANGE_CHECK(*c3_b_BET, 3U);
  for (c3_i18 = 0; c3_i18 < 4; c3_i18++) {
    _SFD_DATA_RANGE_CHECK((*c3_c_u)[c3_i18], 4U);
  }

  chartInstance->c3_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  c3_hoistedGlobal = *c3_b_ALP;
  c3_b_hoistedGlobal = *c3_b_BET;
  c3_c_hoistedGlobal = *c3_b_Vt;
  c3_d_hoistedGlobal = *c3_b_H;
  for (c3_i19 = 0; c3_i19 < 3; c3_i19++) {
    c3_accelerations[c3_i19] = (*c3_c_accelerations)[c3_i19];
  }

  for (c3_i20 = 0; c3_i20 < 3; c3_i20++) {
    c3_omega[c3_i20] = (*c3_c_omega)[c3_i20];
  }

  c3_ALP = c3_hoistedGlobal;
  c3_BET = c3_b_hoistedGlobal;
  for (c3_i21 = 0; c3_i21 < 4; c3_i21++) {
    c3_u[c3_i21] = (*c3_c_u)[c3_i21];
  }

  c3_Vt = c3_c_hoistedGlobal;
  c3_H = c3_d_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 13U, 13U, c3_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 0U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 1U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_accelerations, 2U, c3_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_omega, 3U, c3_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_ALP, 4U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_BET, 5U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_u, 6U, c3_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_Vt, 7U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_H, 8U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_states_OUT, 9U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_cov_OUT, 10U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_innovation, 11U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_S, 12U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 4);
  for (c3_i22 = 0; c3_i22 < 4; c3_i22++) {
    c3_b_u[c3_i22] = c3_u[c3_i22];
  }

  for (c3_i23 = 0; c3_i23 < 3; c3_i23++) {
    c3_b_accelerations[c3_i23] = c3_accelerations[c3_i23];
  }

  for (c3_i24 = 0; c3_i24 < 3; c3_i24++) {
    c3_b_omega[c3_i24] = c3_omega[c3_i24];
  }

  c3_simulateUKF(chartInstance, c3_b_u, c3_Vt, c3_ALP, c3_BET,
                 c3_b_accelerations, c3_b_omega, c3_H, c3_b_states_OUT,
                 c3_b_cov_OUT, c3_b_innovation, c3_b_S);
  for (c3_i25 = 0; c3_i25 < 30; c3_i25++) {
    c3_states_OUT[c3_i25] = c3_b_states_OUT[c3_i25];
  }

  for (c3_i26 = 0; c3_i26 < 30; c3_i26++) {
    c3_cov_OUT[c3_i26] = c3_b_cov_OUT[c3_i26];
  }

  for (c3_i27 = 0; c3_i27 < 6; c3_i27++) {
    c3_innovation[c3_i27] = c3_b_innovation[c3_i27];
  }

  for (c3_i28 = 0; c3_i28 < 6; c3_i28++) {
    c3_S[c3_i28] = c3_b_S[c3_i28];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, -4);
  _SFD_SYMBOL_SCOPE_POP();
  for (c3_i29 = 0; c3_i29 < 30; c3_i29++) {
    (*c3_c_states_OUT)[c3_i29] = c3_states_OUT[c3_i29];
  }

  for (c3_i30 = 0; c3_i30 < 30; c3_i30++) {
    (*c3_c_cov_OUT)[c3_i30] = c3_cov_OUT[c3_i30];
  }

  for (c3_i31 = 0; c3_i31 < 6; c3_i31++) {
    (*c3_c_innovation)[c3_i31] = c3_innovation[c3_i31];
  }

  for (c3_i32 = 0; c3_i32 < 6; c3_i32++) {
    (*c3_c_S)[c3_i32] = c3_S[c3_i32];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY
    (_aircraftControl_FullStateFiltersMachineNumber_, chartInstance->chartNumber,
     chartInstance->instanceNumber);
  for (c3_i33 = 0; c3_i33 < 30; c3_i33++) {
    _SFD_DATA_RANGE_CHECK((*c3_c_states_OUT)[c3_i33], 5U);
  }

  _SFD_DATA_RANGE_CHECK(*c3_b_Vt, 6U);
  for (c3_i34 = 0; c3_i34 < 30; c3_i34++) {
    _SFD_DATA_RANGE_CHECK((*c3_c_cov_OUT)[c3_i34], 7U);
  }

  _SFD_DATA_RANGE_CHECK(*c3_b_H, 8U);
  for (c3_i35 = 0; c3_i35 < 6; c3_i35++) {
    _SFD_DATA_RANGE_CHECK((*c3_c_innovation)[c3_i35], 9U);
  }

  for (c3_i36 = 0; c3_i36 < 6; c3_i36++) {
    _SFD_DATA_RANGE_CHECK((*c3_c_S)[c3_i36], 10U);
  }
}

static void initSimStructsc3_aircraftControl_FullStateFilters
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_simulateUKF(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_u[4], real_T c3_Vt, real_T c3_ALP, real_T c3_BET,
  real_T c3_accelerations[3], real_T c3_omega[3], real_T c3_Height, real_T
  c3_pStatesOUT[30], real_T c3_pCovarianceOUT[30], real_T c3_innovation[6],
  real_T c3_SCovOUT[6])
{
  uint32_T c3_debug_family_var_map[130];
  real_T c3_Qcov[900];
  real_T c3_Rcov[36];
  real_T c3_R2D;
  real_T c3_g;
  real_T c3_vs;
  real_T c3_rho;
  real_T c3_dt;
  real_T c3_m;
  real_T c3_Ix;
  real_T c3_Iy;
  real_T c3_Iz;
  real_T c3_Ixz;
  real_T c3_b;
  real_T c3_cbar;
  real_T c3_S;
  real_T c3_xCG_ref;
  real_T c3_xCG;
  real_T c3_dA;
  real_T c3_dE;
  real_T c3_dR;
  real_T c3_throttle;
  real_T c3_measurement[6];
  real_T c3_qbar;
  real_T c3_hT;
  real_T c3_Tmax;
  real_T c3_T;
  real_T c3_kappa;
  real_T c3_na;
  real_T c3_N;
  real_T c3_sigma[1830];
  real_T c3_sigma_meas[366];
  real_T c3_W[61];
  real_T c3_matrixSQRT[900];
  real_T c3_i;
  real_T c3_p;
  real_T c3_q;
  real_T c3_r;
  real_T c3_prev_pdot;
  real_T c3_prev_qdot;
  real_T c3_prev_rdot;
  real_T c3_CX_dE1_mp;
  real_T c3_CX_dE2_mp;
  real_T c3_CY_dE1_mp;
  real_T c3_CY_dR1_mp;
  real_T c3_CY_dR2_mp;
  real_T c3_CZ_dE1_mp;
  real_T c3_CZ_dE2_mp;
  real_T c3_CZ_dA1_mp;
  real_T c3_Cl_dA1_mp;
  real_T c3_Cl_dA2_mp;
  real_T c3_Cl_dA3_mp;
  real_T c3_Cl_dR1_mp;
  real_T c3_Cl_dR2_mp;
  real_T c3_Cl_dE1_mp;
  real_T c3_Cm_dE1_mp;
  real_T c3_Cm_dE2_mp;
  real_T c3_Cm_dE3_mp;
  real_T c3_Cm_dA1_mp;
  real_T c3_Cn_dA1_mp;
  real_T c3_Cn_dA2_mp;
  real_T c3_Cn_dE1_mp;
  real_T c3_Cn_dE2_mp;
  real_T c3_Cn_dR1_mp;
  real_T c3_Cn_dR2_mp;
  real_T c3_CX_dE1;
  real_T c3_CX_dE2;
  real_T c3_CY_dE1;
  real_T c3_CY_dR1;
  real_T c3_CY_dR2;
  real_T c3_CZ_dE1;
  real_T c3_CZ_dE2;
  real_T c3_CZ_dA1;
  real_T c3_Cl_dA1;
  real_T c3_Cl_dA2;
  real_T c3_Cl_dA3;
  real_T c3_Cl_dR1;
  real_T c3_Cl_dR2;
  real_T c3_Cl_dE1;
  real_T c3_Cm_dE1;
  real_T c3_Cm_dE2;
  real_T c3_Cm_dE3;
  real_T c3_Cm_dA1;
  real_T c3_Cn_dA1;
  real_T c3_Cn_dA2;
  real_T c3_Cn_dE1;
  real_T c3_Cn_dE2;
  real_T c3_Cn_dR1;
  real_T c3_Cn_dR2;
  real_T c3_CX;
  real_T c3_CY;
  real_T c3_CZ;
  real_T c3_Cl;
  real_T c3_Cm;
  real_T c3_Cn;
  real_T c3_c1;
  real_T c3_c2;
  real_T c3_c3;
  real_T c3_c4;
  real_T c3_c5;
  real_T c3_c6;
  real_T c3_c7;
  real_T c3_c8;
  real_T c3_c9;
  real_T c3_heng;
  real_T c3_pdot;
  real_T c3_qdot;
  real_T c3_rdot;
  real_T c3_ymeas[6];
  real_T c3_Pxz[180];
  real_T c3_Pzz[36];
  real_T c3_Scov[36];
  real_T c3_k[180];
  real_T c3_correction[30];
  real_T c3_nargin = 7.0;
  real_T c3_nargout = 4.0;
  int32_T c3_i37;
  int32_T c3_i38;
  int32_T c3_i39;
  real_T c3_d0;
  int32_T c3_i40;
  int32_T c3_i41;
  real_T c3_hoistedGlobal[30];
  int32_T c3_i42;
  int32_T c3_i43;
  real_T c3_b_hoistedGlobal[30];
  real_T c3_dv8[900];
  int32_T c3_i44;
  int32_T c3_i45;
  int32_T c3_i46;
  static real_T c3_dv9[900] = { 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.01 };

  real_T c3_dv10[900];
  real_T c3_dv11[900];
  int32_T c3_i47;
  int32_T c3_i48;
  static real_T c3_dv12[36] = { 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001 };

  real_T c3_A;
  real_T c3_x;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_b_A;
  real_T c3_d_x;
  real_T c3_e_x;
  real_T c3_f_x;
  real_T c3_c_A;
  real_T c3_g_x;
  real_T c3_h_x;
  real_T c3_i_x;
  int32_T c3_i49;
  int32_T c3_i50;
  real_T c3_d_A;
  real_T c3_j_x;
  real_T c3_k_x;
  real_T c3_l_x;
  real_T c3_e_A;
  real_T c3_m_x;
  real_T c3_n_x;
  real_T c3_o_x;
  real_T c3_y;
  real_T c3_f_A;
  real_T c3_p_x;
  real_T c3_q_x;
  real_T c3_r_x;
  real_T c3_b_y;
  real_T c3_g_A;
  real_T c3_s_x;
  real_T c3_t_x;
  real_T c3_u_x;
  real_T c3_c_y;
  real_T c3_h_A;
  real_T c3_v_x;
  real_T c3_w_x;
  real_T c3_x_x;
  real_T c3_d_y;
  real_T c3_i_A;
  real_T c3_y_x;
  real_T c3_ab_x;
  real_T c3_bb_x;
  int32_T c3_i51;
  int32_T c3_i52;
  int32_T c3_i53;
  int32_T c3_i54;
  int32_T c3_i55;
  real_T c3_c_hoistedGlobal[900];
  int32_T c3_i56;
  int32_T c3_j;
  int32_T c3_b_j;
  int32_T c3_info;
  int32_T c3_b_info;
  int32_T c3_c_info;
  int32_T c3_d_info;
  int32_T c3_jmax;
  int32_T c3_a;
  int32_T c3_b_a;
  int32_T c3_b_jmax;
  int32_T c3_b_b;
  int32_T c3_c_b;
  boolean_T c3_overflow;
  int32_T c3_c_j;
  int32_T c3_c_a;
  int32_T c3_d_a;
  int32_T c3_i57;
  int32_T c3_c_jmax;
  int32_T c3_e_a;
  int32_T c3_d_b;
  int32_T c3_f_a;
  int32_T c3_e_b;
  boolean_T c3_b_overflow;
  int32_T c3_b_i;
  int32_T c3_c_i;
  int32_T c3_i58;
  int32_T c3_d_i;
  int32_T c3_e_i;
  int32_T c3_f_i;
  int32_T c3_i59;
  int32_T c3_g_i;
  int32_T c3_h_i;
  int32_T c3_i60;
  int32_T c3_i_i;
  real_T c3_j_A;
  real_T c3_cb_x;
  real_T c3_db_x;
  real_T c3_eb_x;
  real_T c3_e_y;
  real_T c3_k_A;
  real_T c3_fb_x;
  real_T c3_gb_x;
  real_T c3_hb_x;
  real_T c3_f_y;
  real_T c3_l_A;
  real_T c3_B;
  real_T c3_ib_x;
  real_T c3_g_y;
  real_T c3_jb_x;
  real_T c3_h_y;
  real_T c3_kb_x;
  real_T c3_i_y;
  real_T c3_j_y;
  real_T c3_b_B;
  real_T c3_k_y;
  real_T c3_l_y;
  real_T c3_m_y;
  real_T c3_n_y;
  real_T c3_m_A;
  real_T c3_lb_x;
  real_T c3_mb_x;
  real_T c3_nb_x;
  real_T c3_o_y;
  real_T c3_n_A;
  real_T c3_ob_x;
  real_T c3_pb_x;
  real_T c3_qb_x;
  real_T c3_p_y;
  real_T c3_o_A;
  real_T c3_c_B;
  real_T c3_rb_x;
  real_T c3_q_y;
  real_T c3_sb_x;
  real_T c3_r_y;
  real_T c3_tb_x;
  real_T c3_s_y;
  real_T c3_t_y;
  real_T c3_d_B;
  real_T c3_u_y;
  real_T c3_v_y;
  real_T c3_w_y;
  real_T c3_x_y;
  real_T c3_p_A;
  real_T c3_ub_x;
  real_T c3_vb_x;
  real_T c3_wb_x;
  real_T c3_y_y;
  real_T c3_q_A;
  real_T c3_xb_x;
  real_T c3_yb_x;
  real_T c3_ac_x;
  real_T c3_ab_y;
  real_T c3_r_A;
  real_T c3_e_B;
  real_T c3_bc_x;
  real_T c3_bb_y;
  real_T c3_cc_x;
  real_T c3_cb_y;
  real_T c3_dc_x;
  real_T c3_db_y;
  real_T c3_eb_y;
  real_T c3_f_B;
  real_T c3_fb_y;
  real_T c3_gb_y;
  real_T c3_hb_y;
  real_T c3_ib_y;
  real_T c3_s_A;
  real_T c3_ec_x;
  real_T c3_fc_x;
  real_T c3_gc_x;
  real_T c3_jb_y;
  real_T c3_t_A;
  real_T c3_hc_x;
  real_T c3_ic_x;
  real_T c3_jc_x;
  real_T c3_kb_y;
  real_T c3_u_A;
  real_T c3_kc_x;
  real_T c3_lc_x;
  real_T c3_mc_x;
  real_T c3_lb_y;
  real_T c3_v_A;
  real_T c3_nc_x;
  real_T c3_oc_x;
  real_T c3_pc_x;
  real_T c3_mb_y;
  real_T c3_w_A;
  real_T c3_qc_x;
  real_T c3_rc_x;
  real_T c3_sc_x;
  real_T c3_nb_y;
  real_T c3_x_A;
  real_T c3_tc_x;
  real_T c3_uc_x;
  real_T c3_vc_x;
  real_T c3_ob_y;
  int32_T c3_j_i;
  int32_T c3_i61;
  int32_T c3_i62;
  int32_T c3_i63;
  int32_T c3_k_i;
  int32_T c3_i64;
  real_T c3_g_a;
  int32_T c3_l_i;
  int32_T c3_i65;
  real_T c3_C[30];
  int32_T c3_i66;
  int32_T c3_i67;
  int32_T c3_m_i;
  int32_T c3_i68;
  int32_T c3_i69;
  real_T c3_h_a;
  int32_T c3_n_i;
  int32_T c3_i70;
  int32_T c3_i71;
  int32_T c3_i72;
  int32_T c3_o_i;
  int32_T c3_i73;
  real_T c3_f_b[30];
  int32_T c3_i74;
  int32_T c3_i75;
  int32_T c3_i76;
  real_T c3_pb_y[900];
  int32_T c3_i77;
  int32_T c3_i78;
  int32_T c3_p_i;
  real_T c3_i_a;
  int32_T c3_q_i;
  int32_T c3_i79;
  real_T c3_g_b[6];
  int32_T c3_i80;
  int32_T c3_i81;
  int32_T c3_i82;
  int32_T c3_i83;
  int32_T c3_r_i;
  int32_T c3_i84;
  real_T c3_j_a;
  int32_T c3_s_i;
  int32_T c3_i85;
  int32_T c3_i86;
  int32_T c3_t_i;
  int32_T c3_i87;
  real_T c3_h_b[6];
  int32_T c3_i88;
  int32_T c3_i89;
  int32_T c3_i90;
  real_T c3_qb_y[180];
  int32_T c3_i91;
  real_T c3_k_a;
  int32_T c3_u_i;
  int32_T c3_i92;
  int32_T c3_i93;
  int32_T c3_v_i;
  int32_T c3_i94;
  int32_T c3_i95;
  int32_T c3_i96;
  int32_T c3_i97;
  real_T c3_i_b[36];
  int32_T c3_i98;
  int32_T c3_i99;
  int32_T c3_i100;
  real_T c3_l_a[180];
  int32_T c3_i101;
  real_T c3_b_Scov[36];
  int32_T c3_i102;
  int32_T c3_i103;
  int32_T c3_i104;
  real_T c3_dv13[180];
  int32_T c3_i105;
  real_T c3_dv14[36];
  int32_T c3_i106;
  real_T c3_dv15[180];
  int32_T c3_i107;
  real_T c3_dv16[36];
  int32_T c3_i108;
  int32_T c3_i109;
  int32_T c3_i110;
  int32_T c3_i111;
  int32_T c3_i112;
  int32_T c3_i113;
  int32_T c3_i114;
  int32_T c3_i115;
  int32_T c3_i116;
  int32_T c3_i117;
  int32_T c3_i118;
  int32_T c3_i119;
  int32_T c3_i120;
  int32_T c3_i121;
  int32_T c3_i122;
  int32_T c3_i123;
  int32_T c3_i124;
  int32_T c3_i125;
  real_T c3_m_a[180];
  int32_T c3_i126;
  real_T c3_j_b[36];
  int32_T c3_i127;
  int32_T c3_i128;
  int32_T c3_i129;
  int32_T c3_i130;
  real_T c3_k_b[180];
  int32_T c3_i131;
  int32_T c3_i132;
  real_T c3_rb_y[180];
  int32_T c3_i133;
  real_T c3_l_b[180];
  int32_T c3_i134;
  int32_T c3_i135;
  real_T c3_c_Scov[36];
  real_T c3_dv17[6];
  int32_T c3_i136;
  int32_T c3_i137;
  int32_T c3_i138;
  real_T c3_dv18[900];
  real_T c3_dv19[30];
  int32_T c3_i139;
  boolean_T guard1 = false;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 130U, 130U, c3_b_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Qcov, 0U, c3_l_sf_marshallOut,
    c3_l_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Rcov, 1U, c3_k_sf_marshallOut,
    c3_k_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_R2D, 2U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g, 3U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_vs, 4U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_rho, 5U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_dt, 6U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_m, 7U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_Ix, 8U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_Iy, 9U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_Iz, 10U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_Ixz, 11U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_b, 12U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_cbar, 13U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_S, 14U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_xCG_ref, 15U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_xCG, 16U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_dA, 17U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_dE, 18U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_dR, 19U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_throttle, 20U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_measurement, 21U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_qbar, 22U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_hT, 23U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Tmax, 24U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_T, 25U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_kappa, 26U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_na, 27U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_N, 28U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_sigma, 29U, c3_o_sf_marshallOut,
    c3_o_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_sigma_meas, 30U, c3_n_sf_marshallOut,
    c3_n_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_W, 31U, c3_m_sf_marshallOut,
    c3_m_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_matrixSQRT, 32U, c3_l_sf_marshallOut,
    c3_l_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i, 33U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p, 34U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q, 35U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r, 36U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_prev_pdot, 37U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_prev_qdot, 38U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_prev_rdot, 39U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CX_dE1_mp, 40U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CX_dE2_mp, 41U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CY_dE1_mp, 42U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CY_dR1_mp, 43U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CY_dR2_mp, 44U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CZ_dE1_mp, 45U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CZ_dE2_mp, 46U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CZ_dA1_mp, 47U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dA1_mp, 48U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dA2_mp, 49U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dA3_mp, 50U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dR1_mp, 51U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dR2_mp, 52U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dE1_mp, 53U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm_dE1_mp, 54U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm_dE2_mp, 55U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm_dE3_mp, 56U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm_dA1_mp, 57U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dA1_mp, 58U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dA2_mp, 59U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dE1_mp, 60U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dE2_mp, 61U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dR1_mp, 62U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dR2_mp, 63U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CX_dE1, 64U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CX_dE2, 65U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CY_dE1, 66U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CY_dR1, 67U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CY_dR2, 68U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CZ_dE1, 69U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CZ_dE2, 70U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CZ_dA1, 71U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dA1, 72U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dA2, 73U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dA3, 74U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dR1, 75U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dR2, 76U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl_dE1, 77U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm_dE1, 78U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm_dE2, 79U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm_dE3, 80U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm_dA1, 81U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dA1, 82U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dA2, 83U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dE1, 84U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dE2, 85U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dR1, 86U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn_dR2, 87U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CX, 88U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CY, 89U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_CZ, 90U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cl, 91U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cm, 92U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Cn, 93U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c1, 94U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c2, 95U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c3, 96U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c4, 97U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c5, 98U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c6, 99U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c7, 100U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c8, 101U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_c9, 102U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c3_heng, 103U, c3_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_pdot, 104U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_qdot, 105U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_rdot, 106U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_ymeas, 107U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Pxz, 108U, c3_j_sf_marshallOut,
    c3_j_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Pzz, 109U, c3_k_sf_marshallOut,
    c3_k_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_Scov, 110U, c3_k_sf_marshallOut,
    c3_k_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_k, 111U, c3_j_sf_marshallOut,
    c3_j_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_correction, 112U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 113U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 114U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_u, 115U, c3_d_sf_marshallOut,
    c3_i_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Vt, 116U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ALP, 117U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_BET, 118U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_accelerations, 119U,
    c3_e_sf_marshallOut, c3_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_omega, 120U, c3_e_sf_marshallOut,
    c3_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_Height, 121U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_pStatesOUT, 122U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_pCovarianceOUT, 123U,
    c3_b_sf_marshallOut, c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_innovation, 124U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_SCovOUT, 125U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(chartInstance->c3_pStates, 126U,
    c3_i_sf_marshallOut, c3_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(chartInstance->c3_pCovariance, 127U,
    c3_h_sf_marshallOut, c3_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(chartInstance->c3_covDiag, 128U,
    c3_g_sf_marshallOut, c3_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(chartInstance->c3_prevOmegaDots, 129U,
    c3_f_sf_marshallOut, c3_d_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 4);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 5);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 6);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 7);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 9);
  if (CV_SCRIPT_IF(0, 0, !chartInstance->c3_pStates_not_empty)) {
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 11);
    for (c3_i37 = 0; c3_i37 < 30; c3_i37++) {
      chartInstance->c3_pStates[c3_i37] = 0.0;
    }

    chartInstance->c3_pStates_not_empty = true;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 12);
    for (c3_i38 = 0; c3_i38 < 30; c3_i38++) {
      chartInstance->c3_covDiag[c3_i38] = 0.0;
    }

    chartInstance->c3_covDiag_not_empty = true;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 13);
    chartInstance->c3_pStates[0] = c3_accelerations[0];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 13);
    chartInstance->c3_covDiag[0] = c3_mpower(chartInstance, 0.01);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 14);
    chartInstance->c3_pStates[1] = c3_accelerations[1];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 14);
    chartInstance->c3_covDiag[1] = c3_mpower(chartInstance, 0.01);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 15);
    chartInstance->c3_pStates[2] = c3_accelerations[2];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 15);
    chartInstance->c3_covDiag[2] = c3_mpower(chartInstance, 0.01);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 17);
    chartInstance->c3_pStates[3] = c3_omega[0];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 17);
    chartInstance->c3_covDiag[3] = c3_mpower(chartInstance, 0.01);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 18);
    chartInstance->c3_pStates[4] = c3_omega[1];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 18);
    chartInstance->c3_covDiag[4] = c3_mpower(chartInstance, 0.01);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 19);
    chartInstance->c3_pStates[5] = c3_omega[2];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 19);
    chartInstance->c3_covDiag[5] = c3_mpower(chartInstance, 0.01);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 22);
    chartInstance->c3_pStates[6] = 0.00095;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 22);
    chartInstance->c3_covDiag[6] = c3_mpower(chartInstance, 1.0E-5);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 23);
    chartInstance->c3_pStates[7] = 8.5E-7;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 23);
    chartInstance->c3_covDiag[7] = c3_mpower(chartInstance, 1.0E-8);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 24);
    chartInstance->c3_pStates[8] = 0.000175;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 24);
    chartInstance->c3_covDiag[8] = c3_mpower(chartInstance, 1.0E-5);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 25);
    chartInstance->c3_pStates[9] = 0.00155;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 25);
    chartInstance->c3_covDiag[9] = c3_mpower(chartInstance, 0.0001);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 26);
    chartInstance->c3_pStates[10] = 8.0E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 26);
    chartInstance->c3_covDiag[10] = c3_mpower(chartInstance, 1.0E-7);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 27);
    chartInstance->c3_pStates[11] = 0.00476;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 27);
    chartInstance->c3_covDiag[11] = c3_mpower(chartInstance, 0.0001);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 28);
    chartInstance->c3_pStates[12] = 3.3E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 28);
    chartInstance->c3_covDiag[12] = c3_mpower(chartInstance, 1.0E-6);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 29);
    chartInstance->c3_pStates[13] = 7.5E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 29);
    chartInstance->c3_covDiag[13] = c3_mpower(chartInstance, 1.0E-6);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 30);
    chartInstance->c3_pStates[14] = 0.00061;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 30);
    chartInstance->c3_covDiag[14] = c3_mpower(chartInstance, 1.0E-5);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 31);
    chartInstance->c3_pStates[15] = 2.5E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 31);
    chartInstance->c3_covDiag[15] = c3_mpower(chartInstance, 1.0E-6);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 32);
    chartInstance->c3_pStates[16] = 2.6E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 32);
    chartInstance->c3_covDiag[16] = c3_mpower(chartInstance, 1.0E-7);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 33);
    chartInstance->c3_pStates[17] = -0.00023;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 33);
    chartInstance->c3_covDiag[17] = c3_mpower(chartInstance, 1.0E-5);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 34);
    chartInstance->c3_pStates[18] = 4.5E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 34);
    chartInstance->c3_covDiag[18] = c3_mpower(chartInstance, 1.0E-7);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 35);
    chartInstance->c3_pStates[19] = 5.24E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 35);
    chartInstance->c3_covDiag[19] = c3_mpower(chartInstance, 1.0E-6);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 36);
    chartInstance->c3_pStates[20] = 0.00654;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 36);
    chartInstance->c3_covDiag[20] = c3_mpower(chartInstance, 0.0001);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 37);
    chartInstance->c3_pStates[21] = 8.49E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 37);
    chartInstance->c3_covDiag[21] = c3_mpower(chartInstance, 1.0E-6);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 38);
    chartInstance->c3_pStates[22] = 3.74E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 38);
    chartInstance->c3_covDiag[22] = c3_mpower(chartInstance, 1.0E-7);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 39);
    chartInstance->c3_pStates[23] = 3.5E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 39);
    chartInstance->c3_covDiag[23] = c3_mpower(chartInstance, 1.0E-6);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 40);
    chartInstance->c3_pStates[24] = 1.4E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 40);
    chartInstance->c3_covDiag[24] = c3_mpower(chartInstance, 1.0E-6);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 41);
    chartInstance->c3_pStates[25] = 7.0E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 41);
    chartInstance->c3_covDiag[25] = c3_mpower(chartInstance, 1.0E-7);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 42);
    chartInstance->c3_pStates[26] = 8.73E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 42);
    chartInstance->c3_covDiag[26] = c3_mpower(chartInstance, 1.0E-6);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 43);
    chartInstance->c3_pStates[27] = 8.7E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 43);
    chartInstance->c3_covDiag[27] = c3_mpower(chartInstance, 1.0E-7);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 44);
    chartInstance->c3_pStates[28] = 0.0009;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 44);
    chartInstance->c3_covDiag[28] = c3_mpower(chartInstance, 1.0E-5);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 45);
    chartInstance->c3_pStates[29] = 4.0E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 45);
    chartInstance->c3_covDiag[29] = c3_mpower(chartInstance, 1.0E-7);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 47);
    for (c3_i39 = 0; c3_i39 < 24; c3_i39++) {
      chartInstance->c3_pStates[c3_i39 + 6] = 1.0;
    }

    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 48);
    c3_d0 = c3_mpower(chartInstance, 0.01);
    for (c3_i40 = 0; c3_i40 < 24; c3_i40++) {
      chartInstance->c3_covDiag[c3_i40 + 6] = c3_d0;
    }

    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 51);
    for (c3_i41 = 0; c3_i41 < 30; c3_i41++) {
      c3_hoistedGlobal[c3_i41] = chartInstance->c3_covDiag[c3_i41];
    }

    for (c3_i42 = 0; c3_i42 < 30; c3_i42++) {
      c3_hoistedGlobal[c3_i42] *= 0.1;
    }

    for (c3_i43 = 0; c3_i43 < 30; c3_i43++) {
      c3_b_hoistedGlobal[c3_i43] = c3_hoistedGlobal[c3_i43];
    }

    c3_diag(chartInstance, c3_b_hoistedGlobal, c3_dv8);
    for (c3_i44 = 0; c3_i44 < 900; c3_i44++) {
      chartInstance->c3_pCovariance[c3_i44] = c3_dv8[c3_i44];
    }

    chartInstance->c3_pCovariance_not_empty = true;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 53);
    for (c3_i45 = 0; c3_i45 < 183; c3_i45++) {
      chartInstance->c3_prevOmegaDots[c3_i45] = 0.0;
    }

    chartInstance->c3_prevOmegaDots_not_empty = true;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 59);
  for (c3_i46 = 0; c3_i46 < 900; c3_i46++) {
    c3_dv10[c3_i46] = c3_dv9[c3_i46];
  }

  c3_b_power(chartInstance, c3_dv10, c3_dv11);
  for (c3_i47 = 0; c3_i47 < 900; c3_i47++) {
    c3_Qcov[c3_i47] = c3_dv11[c3_i47];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 60);
  c3_Qcov[0] = c3_mpower(chartInstance, 0.01);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 61);
  c3_Qcov[31] = c3_mpower(chartInstance, 0.01);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 62);
  c3_Qcov[62] = c3_mpower(chartInstance, 0.01);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 63);
  c3_Qcov[93] = c3_mpower(chartInstance, 0.01);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 64);
  c3_Qcov[124] = c3_mpower(chartInstance, 0.01);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 65);
  c3_Qcov[155] = c3_mpower(chartInstance, 0.01);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 67);
  for (c3_i48 = 0; c3_i48 < 36; c3_i48++) {
    c3_Rcov[c3_i48] = c3_dv12[c3_i48];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 70);
  c3_R2D = 57.295779513082323;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 71);
  c3_g = 9.81;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 72);
  c3_vs = 340.3;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 73);
  c3_R2D = 57.295779513082323;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 74);
  c3_rho = 1.225;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 75);
  c3_dt = 0.005;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 76);
  c3_m = 1177.0410005333333;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 78);
  c3_rho = 1.225;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 79);
  c3_Ix = 2256.9839826075154;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 80);
  c3_Iy = 11044.488299351713;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 81);
  c3_Iz = 12636.21789221188;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 82);
  c3_Ixz = 106.20569401537166;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 83);
  c3_b = 6.666666666666667;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 84);
  c3_cbar = 3.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 85);
  c3_S = 20.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 87);
  c3_xCG_ref = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 88);
  c3_xCG = 0.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 90);
  c3_A = c3_u[0] * 180.0;
  c3_x = c3_A;
  c3_b_x = c3_x;
  c3_c_x = c3_b_x;
  c3_dA = c3_c_x / 3.1415926535897931;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 91);
  c3_b_A = c3_u[1] * 180.0;
  c3_d_x = c3_b_A;
  c3_e_x = c3_d_x;
  c3_f_x = c3_e_x;
  c3_dE = c3_f_x / 3.1415926535897931;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 92);
  c3_c_A = c3_u[2] * 180.0;
  c3_g_x = c3_c_A;
  c3_h_x = c3_g_x;
  c3_i_x = c3_h_x;
  c3_dR = c3_i_x / 3.1415926535897931;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 94);
  c3_throttle = c3_u[3];
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 96);
  for (c3_i49 = 0; c3_i49 < 3; c3_i49++) {
    c3_measurement[c3_i49] = c3_accelerations[c3_i49];
  }

  for (c3_i50 = 0; c3_i50 < 3; c3_i50++) {
    c3_measurement[c3_i50 + 3] = c3_omega[c3_i50];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 99);
  c3_qbar = 0.6125 * c3_mpower(chartInstance, c3_Vt);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 102);
  c3_d_A = c3_Height;
  c3_j_x = c3_d_A;
  c3_k_x = c3_j_x;
  c3_l_x = c3_k_x;
  c3_hT = c3_l_x / 3048.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 104);
  c3_e_A = c3_Vt;
  c3_m_x = c3_e_A;
  c3_n_x = c3_m_x;
  c3_o_x = c3_n_x;
  c3_y = c3_o_x / 340.3;
  c3_f_A = c3_Vt;
  c3_p_x = c3_f_A;
  c3_q_x = c3_p_x;
  c3_r_x = c3_q_x;
  c3_b_y = c3_r_x / 340.3;
  c3_g_A = c3_Vt;
  c3_s_x = c3_g_A;
  c3_t_x = c3_s_x;
  c3_u_x = c3_t_x;
  c3_c_y = c3_u_x / 340.3;
  c3_h_A = c3_Vt;
  c3_v_x = c3_h_A;
  c3_w_x = c3_v_x;
  c3_x_x = c3_w_x;
  c3_d_y = c3_x_x / 340.3;
  c3_i_A = ((((((((30.21 - 0.668 * c3_hT) - 6.877 * c3_mpower(chartInstance,
    c3_hT)) + 1.951 * c3_b_mpower(chartInstance, c3_hT)) - 0.1512 * c3_c_mpower
                (chartInstance, c3_hT)) + c3_y * ((((-33.8 + 3.347 * c3_hT) +
    18.13 * c3_mpower(chartInstance, c3_hT)) - 5.865 * c3_b_mpower(chartInstance,
    c3_hT)) + 0.4757 * c3_c_mpower(chartInstance, c3_hT))) + c3_mpower
              (chartInstance, c3_b_y) * ((((100.8 - 77.56 * c3_hT) + 5.441 *
    c3_mpower(chartInstance, c3_hT)) + 2.864 * c3_b_mpower(chartInstance, c3_hT))
    - 0.3355 * c3_c_mpower(chartInstance, c3_hT))) + c3_b_mpower(chartInstance,
              c3_c_y) * ((((-78.99 + 101.4 * c3_hT) - 30.28 * c3_mpower
    (chartInstance, c3_hT)) + 3.236 * c3_b_mpower(chartInstance, c3_hT)) -
              0.1089 * c3_c_mpower(chartInstance, c3_hT))) + c3_c_mpower
            (chartInstance, c3_d_y) * ((((18.74 - 31.6 * c3_hT) + 12.04 *
    c3_mpower(chartInstance, c3_hT)) - 1.785 * c3_b_mpower(chartInstance, c3_hT))
             + 0.09417 * c3_c_mpower(chartInstance, c3_hT))) * 4448.22;
  c3_y_x = c3_i_A;
  c3_ab_x = c3_y_x;
  c3_bb_x = c3_ab_x;
  c3_Tmax = c3_bb_x / 20.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 110);
  c3_T = c3_Tmax * c3_throttle;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 115);
  c3_kappa = 0.0001;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 116);
  c3_na = 30.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 117);
  c3_N = 6.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 119);
  for (c3_i51 = 0; c3_i51 < 1830; c3_i51++) {
    c3_sigma[c3_i51] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 120);
  for (c3_i52 = 0; c3_i52 < 366; c3_i52++) {
    c3_sigma_meas[c3_i52] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 121);
  for (c3_i53 = 0; c3_i53 < 61; c3_i53++) {
    c3_W[c3_i53] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 123);
  for (c3_i54 = 0; c3_i54 < 30; c3_i54++) {
    c3_sigma[c3_i54] = chartInstance->c3_pStates[c3_i54];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 124);
  c3_W[0] = 3.3333222222592593E-6;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 126);
  for (c3_i55 = 0; c3_i55 < 900; c3_i55++) {
    c3_c_hoistedGlobal[c3_i55] = chartInstance->c3_pCovariance[c3_i55];
  }

  for (c3_i56 = 0; c3_i56 < 900; c3_i56++) {
    c3_c_hoistedGlobal[c3_i56] *= 30.0001;
  }

  for (c3_j = 1; c3_j < 31; c3_j++) {
    c3_b_j = c3_j;
  }

  c3_info = c3_b_eml_matlab_zpotrf(chartInstance, c3_c_hoistedGlobal);
  c3_b_info = c3_info;
  c3_c_info = c3_b_info;
  c3_d_info = c3_c_info;
  if (c3_d_info == 0) {
    c3_jmax = 30;
  } else {
    c3_b_eml_error(chartInstance);
    c3_a = c3_d_info;
    c3_b_a = c3_a - 1;
    c3_jmax = c3_b_a;
  }

  c3_b_jmax = c3_jmax;
  c3_b_b = c3_b_jmax;
  c3_c_b = c3_b_b;
  if (1 > c3_c_b) {
    c3_overflow = false;
  } else {
    c3_eml_switch_helper(chartInstance);
    c3_overflow = (c3_c_b > 2147483646);
  }

  if (c3_overflow) {
    c3_check_forloop_overflow_error(chartInstance, c3_overflow);
  }

  for (c3_c_j = 1; c3_c_j <= c3_b_jmax; c3_c_j++) {
    c3_b_j = c3_c_j;
    c3_c_a = c3_b_j;
    c3_d_a = c3_c_a + 1;
    c3_i57 = c3_d_a;
    c3_c_jmax = c3_jmax;
    c3_e_a = c3_i57;
    c3_d_b = c3_c_jmax;
    c3_f_a = c3_e_a;
    c3_e_b = c3_d_b;
    if (c3_f_a > c3_e_b) {
      c3_b_overflow = false;
    } else {
      c3_eml_switch_helper(chartInstance);
      c3_b_overflow = (c3_e_b > 2147483646);
    }

    if (c3_b_overflow) {
      c3_check_forloop_overflow_error(chartInstance, c3_b_overflow);
    }

    for (c3_b_i = c3_i57; c3_b_i <= c3_c_jmax; c3_b_i++) {
      c3_c_i = c3_b_i;
      c3_c_hoistedGlobal[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c3_c_i), 1, 30, 1, 0) + 30 *
                          (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c3_b_j), 1, 30, 2, 0) - 1)) - 1] = 0.0;
    }
  }

  for (c3_i58 = 0; c3_i58 < 900; c3_i58++) {
    c3_matrixSQRT[c3_i58] = c3_c_hoistedGlobal[c3_i58];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 128U);
  c3_i = 1.0;
  c3_d_i = 0;
  while (c3_d_i < 60) {
    c3_i = 1.0 + (real_T)c3_d_i;
    CV_SCRIPT_FOR(0, 0, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 130U);
    guard1 = false;
    if (CV_SCRIPT_COND(0, 0, c3_i >= 1.0)) {
      if (CV_SCRIPT_COND(0, 1, c3_i <= c3_na)) {
        CV_SCRIPT_MCDC(0, 0, true);
        CV_SCRIPT_IF(0, 1, true);
        _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 131U);
        c3_e_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
          _SFD_INTEGER_CHECK("i+1", c3_i + 1.0), 1, 61, 2, 0) - 1;
        c3_f_i = _SFD_EML_ARRAY_BOUNDS_CHECK("matrixSQRT", (int32_T)
          _SFD_INTEGER_CHECK("i", c3_i), 1, 30, 1, 0) - 1;
        for (c3_i59 = 0; c3_i59 < 30; c3_i59++) {
          c3_sigma[c3_i59 + 30 * c3_e_i] = chartInstance->c3_pStates[c3_i59] +
            c3_matrixSQRT[c3_f_i + 30 * c3_i59];
        }
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      CV_SCRIPT_MCDC(0, 0, false);
      CV_SCRIPT_IF(0, 1, false);
      _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 133U);
      c3_g_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK(
        "i+1", c3_i + 1.0), 1, 61, 2, 0) - 1;
      c3_h_i = _SFD_EML_ARRAY_BOUNDS_CHECK("matrixSQRT", (int32_T)
        _SFD_INTEGER_CHECK("i-na", c3_i - c3_na), 1, 30, 1, 0) - 1;
      for (c3_i60 = 0; c3_i60 < 30; c3_i60++) {
        c3_sigma[c3_i60 + 30 * c3_g_i] = chartInstance->c3_pStates[c3_i60] -
          c3_matrixSQRT[c3_h_i + 30 * c3_i60];
      }
    }

    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 136U);
    c3_W[_SFD_EML_ARRAY_BOUNDS_CHECK("W", (int32_T)_SFD_INTEGER_CHECK("i+1",
      c3_i + 1.0), 1, 61, 1, 0) - 1] = 0.016666611111296296;
    c3_d_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_SCRIPT_FOR(0, 0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 141U);
  c3_i = 1.0;
  c3_i_i = 0;
  while (c3_i_i < 61) {
    c3_i = 1.0 + (real_T)c3_i_i;
    CV_SCRIPT_FOR(0, 1, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 143U);
    c3_p = c3_sigma[3 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 144U);
    c3_q = c3_sigma[4 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 145U);
    c3_r = c3_sigma[5 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 147U);
    c3_prev_pdot = chartInstance->c3_prevOmegaDots[3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("prevOmegaDots", (int32_T)_SFD_INTEGER_CHECK(
         "i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 148U);
    c3_prev_qdot = chartInstance->c3_prevOmegaDots[1 + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("prevOmegaDots", (int32_T)_SFD_INTEGER_CHECK(
         "i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 149U);
    c3_prev_rdot = chartInstance->c3_prevOmegaDots[2 + 3 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("prevOmegaDots", (int32_T)_SFD_INTEGER_CHECK(
         "i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 151U);
    c3_CX_dE1_mp = c3_sigma[6 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 152U);
    c3_CX_dE2_mp = c3_sigma[7 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 153U);
    c3_CY_dE1_mp = c3_sigma[8 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 154U);
    c3_CY_dR1_mp = c3_sigma[9 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 155U);
    c3_CY_dR2_mp = c3_sigma[10 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 156U);
    c3_CZ_dE1_mp = c3_sigma[11 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 157U);
    c3_CZ_dE2_mp = c3_sigma[12 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 158U);
    c3_CZ_dA1_mp = c3_sigma[13 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 159U);
    c3_Cl_dA1_mp = c3_sigma[14 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 160U);
    c3_Cl_dA2_mp = c3_sigma[15 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 161U);
    c3_Cl_dA3_mp = c3_sigma[16 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 162U);
    c3_Cl_dR1_mp = c3_sigma[17 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 163U);
    c3_Cl_dR2_mp = c3_sigma[18 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 164U);
    c3_Cl_dE1_mp = c3_sigma[19 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 165U);
    c3_Cm_dE1_mp = c3_sigma[20 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 166U);
    c3_Cm_dE2_mp = c3_sigma[21 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 167U);
    c3_Cm_dE3_mp = c3_sigma[22 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 168U);
    c3_Cm_dA1_mp = c3_sigma[23 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 169U);
    c3_Cn_dA1_mp = c3_sigma[24 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 170U);
    c3_Cn_dA2_mp = c3_sigma[25 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 171U);
    c3_Cn_dE1_mp = c3_sigma[26 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 172U);
    c3_Cn_dE2_mp = c3_sigma[27 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 173U);
    c3_Cn_dR1_mp = c3_sigma[28 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 174U);
    c3_Cn_dR2_mp = c3_sigma[29 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma",
      (int32_T)_SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 176U);
    c3_CX_dE1 = c3_CX_dE1_mp * 0.00095;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 177U);
    c3_CX_dE2 = c3_CX_dE2_mp * 8.5E-7;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 178U);
    c3_CY_dE1 = c3_CY_dE1_mp * 0.000175;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 179U);
    c3_CY_dR1 = c3_CY_dR1_mp * 0.00155;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 180U);
    c3_CY_dR2 = c3_CY_dR2_mp * 8.0E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 181U);
    c3_CZ_dE1 = c3_CZ_dE1_mp * 0.00476;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 182U);
    c3_CZ_dE2 = c3_CZ_dE2_mp * 3.3E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 183U);
    c3_CZ_dA1 = c3_CZ_dA1_mp * 7.5E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 184U);
    c3_Cl_dA1 = c3_Cl_dA1_mp * 0.00061;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 185U);
    c3_Cl_dA2 = c3_Cl_dA2_mp * 2.5E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 186U);
    c3_Cl_dA3 = c3_Cl_dA3_mp * 2.6E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 187U);
    c3_Cl_dR1 = c3_Cl_dR1_mp * -0.00023;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 188U);
    c3_Cl_dR2 = c3_Cl_dR2_mp * 4.5E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 189U);
    c3_Cl_dE1 = c3_Cl_dE1_mp * 5.24E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 190U);
    c3_Cm_dE1 = c3_Cm_dE1_mp * 0.00654;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 191U);
    c3_Cm_dE2 = c3_Cm_dE2_mp * 8.49E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 192U);
    c3_Cm_dE3 = c3_Cm_dE3_mp * 3.74E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 193U);
    c3_Cm_dA1 = c3_Cm_dA1_mp * 3.5E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 194U);
    c3_Cn_dA1 = c3_Cn_dA1_mp * 1.4E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 195U);
    c3_Cn_dA2 = c3_Cn_dA2_mp * 7.0E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 196U);
    c3_Cn_dE1 = c3_Cn_dE1_mp * 8.73E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 197U);
    c3_Cn_dE2 = c3_Cn_dE2_mp * 8.7E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 198U);
    c3_Cn_dR1 = c3_Cn_dR1_mp * 0.0009;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 199U);
    c3_Cn_dR2 = c3_Cn_dR2_mp * 4.0E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 208U);
    c3_j_A = 180.0 * c3_q * 3.0;
    c3_cb_x = c3_j_A;
    c3_db_x = c3_cb_x;
    c3_eb_x = c3_db_x;
    c3_e_y = c3_eb_x / 3.1415926535897931;
    c3_k_A = c3_e_y;
    c3_fb_x = c3_k_A;
    c3_gb_x = c3_fb_x;
    c3_hb_x = c3_gb_x;
    c3_f_y = c3_hb_x / 2.0;
    c3_l_A = c3_f_y;
    c3_B = c3_Vt;
    c3_ib_x = c3_l_A;
    c3_g_y = c3_B;
    c3_jb_x = c3_ib_x;
    c3_h_y = c3_g_y;
    c3_kb_x = c3_jb_x;
    c3_i_y = c3_h_y;
    c3_j_y = c3_kb_x / c3_i_y;
    c3_CX = (((((-0.0434 + 0.00239 * c3_ALP) + 2.53E-5 * c3_mpower(chartInstance,
      c3_BET)) - 1.07E-6 * c3_ALP * c3_mpower(chartInstance, c3_BET)) +
              c3_CX_dE1 * c3_dE) - c3_CX_dE2 * c3_dE * c3_mpower(chartInstance,
              c3_BET)) + c3_j_y * ((0.00873 + 0.001 * c3_ALP) - 0.000175 *
      c3_mpower(chartInstance, c3_ALP));
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 209U);
    c3_b_B = c3_Vt;
    c3_k_y = c3_b_B;
    c3_l_y = c3_k_y;
    c3_m_y = c3_l_y;
    c3_n_y = 190.9859317102744 / c3_m_y;
    c3_CY = ((-0.012 * c3_BET + c3_CY_dR1 * c3_dR) - c3_CY_dR2 * c3_dR * c3_ALP)
      + c3_n_y * (((0.00225 * c3_p + 0.0117 * c3_r) - 0.000367 * c3_r * c3_ALP)
                  + c3_CY_dE1 * c3_r * c3_dE);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 210U);
    c3_m_A = 180.0 * c3_q * 3.0;
    c3_lb_x = c3_m_A;
    c3_mb_x = c3_lb_x;
    c3_nb_x = c3_mb_x;
    c3_o_y = c3_nb_x / 3.1415926535897931;
    c3_n_A = c3_o_y;
    c3_ob_x = c3_n_A;
    c3_pb_x = c3_ob_x;
    c3_qb_x = c3_pb_x;
    c3_p_y = c3_qb_x / 2.0;
    c3_o_A = c3_p_y;
    c3_c_B = c3_Vt;
    c3_rb_x = c3_o_A;
    c3_q_y = c3_c_B;
    c3_sb_x = c3_rb_x;
    c3_r_y = c3_q_y;
    c3_tb_x = c3_sb_x;
    c3_s_y = c3_r_y;
    c3_t_y = c3_tb_x / c3_s_y;
    c3_CZ = ((((-0.131 - 0.0538 * c3_ALP) - c3_CZ_dE1 * c3_dE) - c3_CZ_dE2 *
              c3_dE * c3_ALP) - c3_CZ_dA1 * c3_power(chartInstance, c3_dA)) +
      c3_t_y * ((-0.111 + 0.00517 * c3_ALP) - 0.0011 * c3_mpower(chartInstance,
      c3_ALP));
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 211U);
    c3_d_B = c3_Vt;
    c3_u_y = c3_d_B;
    c3_v_y = c3_u_y;
    c3_w_y = c3_v_y;
    c3_x_y = 190.9859317102744 / c3_w_y;
    c3_Cl = ((((-0.000598 * c3_BET - 0.000283 * c3_ALP * c3_BET) + 1.51E-5 *
               c3_mpower(chartInstance, c3_ALP) * c3_BET) - c3_dA * ((c3_Cl_dA1
                + c3_Cl_dA2 * c3_ALP) - c3_Cl_dA3 * c3_mpower(chartInstance,
                c3_ALP))) - c3_dR * (c3_Cl_dR1 + c3_Cl_dR2 * c3_ALP)) + c3_x_y *
      (((((-0.00412 * c3_p - 0.000524 * c3_p * c3_ALP) + 4.36E-5 * c3_p *
          c3_mpower(chartInstance, c3_ALP)) + 0.000436 * c3_r) + 0.000105 * c3_r
        * c3_ALP) + c3_Cl_dE1 * c3_r * c3_dE);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 212U);
    c3_p_A = 180.0 * c3_q * 3.0;
    c3_ub_x = c3_p_A;
    c3_vb_x = c3_ub_x;
    c3_wb_x = c3_vb_x;
    c3_y_y = c3_wb_x / 2.0;
    c3_q_A = c3_y_y;
    c3_xb_x = c3_q_A;
    c3_yb_x = c3_xb_x;
    c3_ac_x = c3_yb_x;
    c3_ab_y = c3_ac_x / 3.1415926535897931;
    c3_r_A = c3_ab_y;
    c3_e_B = c3_Vt;
    c3_bc_x = c3_r_A;
    c3_bb_y = c3_e_B;
    c3_cc_x = c3_bc_x;
    c3_cb_y = c3_bb_y;
    c3_dc_x = c3_cc_x;
    c3_db_y = c3_cb_y;
    c3_eb_y = c3_dc_x / c3_db_y;
    c3_Cm = ((((((((-0.00661 - 0.00267 * c3_ALP) - 6.48E-5 * c3_mpower
                   (chartInstance, c3_BET)) - 2.65E-6 * c3_ALP * c3_mpower
                  (chartInstance, c3_BET)) - c3_Cm_dE1 * c3_dE) - c3_Cm_dE2 *
                c3_dE * c3_ALP) + c3_Cm_dE3 * c3_dE * c3_mpower(chartInstance,
                c3_BET)) - c3_Cm_dA1 * c3_power(chartInstance, c3_dA)) + c3_eb_y
             * (-0.0473 - 0.00157 * c3_ALP)) + 0.0 * c3_CZ;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 213U);
    c3_f_B = c3_Vt;
    c3_fb_y = c3_f_B;
    c3_gb_y = c3_fb_y;
    c3_hb_y = c3_gb_y;
    c3_ib_y = 190.9859317102744 / c3_hb_y;
    c3_Cn = ((((((0.00228 * c3_BET + 1.79E-6 * c3_b_mpower(chartInstance, c3_BET))
                 + c3_Cn_dA1 * c3_dA) + c3_Cn_dA2 * c3_dA * c3_ALP) - c3_Cn_dR1 *
               c3_dR) + c3_Cn_dR2 * c3_dR * c3_ALP) + c3_ib_y * (((((-6.63E-5 *
      c3_p - 1.92E-5 * c3_p * c3_ALP) + 5.06E-6 * c3_p * c3_mpower(chartInstance,
      c3_ALP)) - 0.00606 * c3_r) - c3_Cn_dE1 * c3_r * c3_dE) + c3_Cn_dE2 * c3_r *
              c3_dE * c3_ALP)) - 0.0 * c3_CY;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 215U);
    c3_c1 = -0.70592099279383547;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 216U);
    c3_c2 = 0.014338034095370216;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 217U);
    c3_c3 = 0.000443244465804795;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 218U);
    c3_c4 = 3.7254094944251367E-6;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 219U);
    c3_c5 = 0.93976593829282273;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 220U);
    c3_c6 = 0.0096161715361322529;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 221U);
    c3_c7 = 9.054290003265229E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 222U);
    c3_c8 = -0.69609284164731566;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 223U);
    c3_c9 = 7.916891495812398E-5;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 225U);
    c3_heng = 0.0;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 226U);
    c3_pdot = (-0.70592099279383547 * c3_r + 0.014338034095370216 * c3_p) * c3_q
      + c3_qbar * 20.0 * 6.666666666666667 * (0.000443244465804795 * c3_Cl +
      3.7254094944251367E-6 * c3_Cn);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 227U);
    c3_qdot = (0.93976593829282273 * c3_p * c3_r - 0.0096161715361322529 *
               (c3_power(chartInstance, c3_p) - c3_power(chartInstance, c3_r)))
      + c3_qbar * 20.0 * 3.0 * 9.054290003265229E-5 * c3_Cm;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 228U);
    c3_rdot = (-0.69609284164731566 * c3_p - 0.014338034095370216 * c3_r) * c3_q
      + c3_qbar * 20.0 * 6.666666666666667 * (3.7254094944251367E-6 * c3_Cl +
      7.916891495812398E-5 * c3_Cn);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 232U);
    c3_s_A = c3_qbar * 20.0 * c3_CX + c3_T;
    c3_ec_x = c3_s_A;
    c3_fc_x = c3_ec_x;
    c3_gc_x = c3_fc_x;
    c3_jb_y = c3_gc_x / 1177.0410005333333;
    c3_sigma[30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_jb_y;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 233U);
    c3_t_A = c3_qbar * 20.0 * c3_CY;
    c3_hc_x = c3_t_A;
    c3_ic_x = c3_hc_x;
    c3_jc_x = c3_ic_x;
    c3_kb_y = c3_jc_x / 1177.0410005333333;
    c3_sigma[1 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_kb_y;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 234U);
    c3_u_A = c3_qbar * 20.0 * c3_CZ;
    c3_kc_x = c3_u_A;
    c3_lc_x = c3_kc_x;
    c3_mc_x = c3_lc_x;
    c3_lb_y = c3_mc_x / 1177.0410005333333;
    c3_sigma[2 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_lb_y;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 238U);
    c3_v_A = 0.005 * (c3_prev_pdot + c3_pdot);
    c3_nc_x = c3_v_A;
    c3_oc_x = c3_nc_x;
    c3_pc_x = c3_oc_x;
    c3_mb_y = c3_pc_x / 2.0;
    c3_sigma[3 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_p + c3_mb_y;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 239U);
    c3_w_A = 0.005 * (c3_prev_qdot + c3_qdot);
    c3_qc_x = c3_w_A;
    c3_rc_x = c3_qc_x;
    c3_sc_x = c3_rc_x;
    c3_nb_y = c3_sc_x / 2.0;
    c3_sigma[4 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_q + c3_nb_y;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 240U);
    c3_x_A = 0.005 * (c3_prev_rdot + c3_rdot);
    c3_tc_x = c3_x_A;
    c3_uc_x = c3_tc_x;
    c3_vc_x = c3_uc_x;
    c3_ob_y = c3_vc_x / 2.0;
    c3_sigma[5 + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_r + c3_ob_y;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 242U);
    c3_j_i = _SFD_EML_ARRAY_BOUNDS_CHECK("prevOmegaDots", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1;
    chartInstance->c3_prevOmegaDots[3 * c3_j_i] = c3_pdot;
    chartInstance->c3_prevOmegaDots[1 + 3 * c3_j_i] = c3_qdot;
    chartInstance->c3_prevOmegaDots[2 + 3 * c3_j_i] = c3_rdot;
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 245U);
    c3_sigma_meas[6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_sigma[30 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 246U);
    c3_sigma_meas[1 + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_sigma[1 + 30 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 247U);
    c3_sigma_meas[2 + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_sigma[2 + 30 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 249U);
    c3_sigma_meas[3 + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_sigma[3 + 30 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 250U);
    c3_sigma_meas[4 + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_sigma[4 + 30 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 61, 2, 0) - 1)];
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 251U);
    c3_sigma_meas[5 + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1)] = c3_sigma[5 + 30 *
      (_SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK("i",
         c3_i), 1, 61, 2, 0) - 1)];
    c3_i_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_SCRIPT_FOR(0, 1, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, MAX_uint8_T);
  for (c3_i61 = 0; c3_i61 < 30; c3_i61++) {
    chartInstance->c3_pStates[c3_i61] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 256);
  for (c3_i62 = 0; c3_i62 < 900; c3_i62++) {
    chartInstance->c3_pCovariance[c3_i62] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 257);
  for (c3_i63 = 0; c3_i63 < 6; c3_i63++) {
    c3_ymeas[c3_i63] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 259);
  c3_i = 1.0;
  c3_k_i = 0;
  while (c3_k_i < 61) {
    c3_i = 1.0 + (real_T)c3_k_i;
    CV_SCRIPT_FOR(0, 2, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 260);
    for (c3_i64 = 0; c3_i64 < 30; c3_i64++) {
      c3_hoistedGlobal[c3_i64] = chartInstance->c3_pStates[c3_i64];
    }

    c3_g_a = c3_W[_SFD_EML_ARRAY_BOUNDS_CHECK("W", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 1, 0) - 1];
    c3_l_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 2, 0) - 1;
    for (c3_i65 = 0; c3_i65 < 30; c3_i65++) {
      c3_C[c3_i65] = c3_sigma[c3_i65 + 30 * c3_l_i];
    }

    for (c3_i66 = 0; c3_i66 < 30; c3_i66++) {
      c3_C[c3_i66] *= c3_g_a;
    }

    for (c3_i67 = 0; c3_i67 < 30; c3_i67++) {
      chartInstance->c3_pStates[c3_i67] = c3_hoistedGlobal[c3_i67] + c3_C[c3_i67];
    }

    c3_k_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_SCRIPT_FOR(0, 2, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 263);
  c3_i = 1.0;
  c3_m_i = 0;
  while (c3_m_i < 61) {
    c3_i = 1.0 + (real_T)c3_m_i;
    CV_SCRIPT_FOR(0, 3, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 264);
    for (c3_i68 = 0; c3_i68 < 900; c3_i68++) {
      c3_c_hoistedGlobal[c3_i68] = chartInstance->c3_pCovariance[c3_i68];
    }

    for (c3_i69 = 0; c3_i69 < 30; c3_i69++) {
      c3_hoistedGlobal[c3_i69] = chartInstance->c3_pStates[c3_i69];
    }

    c3_h_a = c3_W[_SFD_EML_ARRAY_BOUNDS_CHECK("W", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 1, 0) - 1];
    c3_n_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 2, 0) - 1;
    for (c3_i70 = 0; c3_i70 < 30; c3_i70++) {
      c3_hoistedGlobal[c3_i70] = c3_sigma[c3_i70 + 30 * c3_n_i] -
        c3_hoistedGlobal[c3_i70];
    }

    for (c3_i71 = 0; c3_i71 < 30; c3_i71++) {
      c3_hoistedGlobal[c3_i71] *= c3_h_a;
    }

    for (c3_i72 = 0; c3_i72 < 30; c3_i72++) {
      c3_C[c3_i72] = chartInstance->c3_pStates[c3_i72];
    }

    c3_o_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 2, 0) - 1;
    for (c3_i73 = 0; c3_i73 < 30; c3_i73++) {
      c3_f_b[c3_i73] = c3_sigma[c3_i73 + 30 * c3_o_i] - c3_C[c3_i73];
    }

    for (c3_i74 = 0; c3_i74 < 30; c3_i74++) {
      c3_i75 = 0;
      for (c3_i76 = 0; c3_i76 < 30; c3_i76++) {
        c3_pb_y[c3_i75 + c3_i74] = c3_hoistedGlobal[c3_i74] * c3_f_b[c3_i76];
        c3_i75 += 30;
      }
    }

    for (c3_i77 = 0; c3_i77 < 900; c3_i77++) {
      chartInstance->c3_pCovariance[c3_i77] = c3_c_hoistedGlobal[c3_i77] +
        c3_pb_y[c3_i77];
    }

    c3_m_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_SCRIPT_FOR(0, 3, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 267);
  for (c3_i78 = 0; c3_i78 < 900; c3_i78++) {
    chartInstance->c3_pCovariance[c3_i78] += c3_Qcov[c3_i78];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 269);
  c3_i = 1.0;
  c3_p_i = 0;
  while (c3_p_i < 61) {
    c3_i = 1.0 + (real_T)c3_p_i;
    CV_SCRIPT_FOR(0, 4, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 270);
    c3_i_a = c3_W[_SFD_EML_ARRAY_BOUNDS_CHECK("W", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 1, 0) - 1];
    c3_q_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1;
    for (c3_i79 = 0; c3_i79 < 6; c3_i79++) {
      c3_g_b[c3_i79] = c3_sigma_meas[c3_i79 + 6 * c3_q_i];
    }

    for (c3_i80 = 0; c3_i80 < 6; c3_i80++) {
      c3_g_b[c3_i80] *= c3_i_a;
    }

    for (c3_i81 = 0; c3_i81 < 6; c3_i81++) {
      c3_ymeas[c3_i81] += c3_g_b[c3_i81];
    }

    c3_p_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_SCRIPT_FOR(0, 4, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 274);
  for (c3_i82 = 0; c3_i82 < 180; c3_i82++) {
    c3_Pxz[c3_i82] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 275);
  for (c3_i83 = 0; c3_i83 < 36; c3_i83++) {
    c3_Pzz[c3_i83] = 0.0;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 277);
  c3_i = 1.0;
  c3_r_i = 0;
  while (c3_r_i < 61) {
    c3_i = 1.0 + (real_T)c3_r_i;
    CV_SCRIPT_FOR(0, 5, 1);
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 278);
    for (c3_i84 = 0; c3_i84 < 30; c3_i84++) {
      c3_hoistedGlobal[c3_i84] = chartInstance->c3_pStates[c3_i84];
    }

    c3_j_a = c3_W[_SFD_EML_ARRAY_BOUNDS_CHECK("W", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 1, 0) - 1];
    c3_s_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 2, 0) - 1;
    for (c3_i85 = 0; c3_i85 < 30; c3_i85++) {
      c3_hoistedGlobal[c3_i85] = c3_sigma[c3_i85 + 30 * c3_s_i] -
        c3_hoistedGlobal[c3_i85];
    }

    for (c3_i86 = 0; c3_i86 < 30; c3_i86++) {
      c3_hoistedGlobal[c3_i86] *= c3_j_a;
    }

    c3_t_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1;
    for (c3_i87 = 0; c3_i87 < 6; c3_i87++) {
      c3_h_b[c3_i87] = c3_sigma_meas[c3_i87 + 6 * c3_t_i] - c3_ymeas[c3_i87];
    }

    for (c3_i88 = 0; c3_i88 < 30; c3_i88++) {
      c3_i89 = 0;
      for (c3_i90 = 0; c3_i90 < 6; c3_i90++) {
        c3_qb_y[c3_i89 + c3_i88] = c3_hoistedGlobal[c3_i88] * c3_h_b[c3_i90];
        c3_i89 += 30;
      }
    }

    for (c3_i91 = 0; c3_i91 < 180; c3_i91++) {
      c3_Pxz[c3_i91] += c3_qb_y[c3_i91];
    }

    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 279);
    c3_k_a = c3_W[_SFD_EML_ARRAY_BOUNDS_CHECK("W", (int32_T)_SFD_INTEGER_CHECK(
      "i", c3_i), 1, 61, 1, 0) - 1];
    c3_u_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1;
    for (c3_i92 = 0; c3_i92 < 6; c3_i92++) {
      c3_g_b[c3_i92] = c3_sigma_meas[c3_i92 + 6 * c3_u_i] - c3_ymeas[c3_i92];
    }

    for (c3_i93 = 0; c3_i93 < 6; c3_i93++) {
      c3_g_b[c3_i93] *= c3_k_a;
    }

    c3_v_i = _SFD_EML_ARRAY_BOUNDS_CHECK("sigma_meas", (int32_T)
      _SFD_INTEGER_CHECK("i", c3_i), 1, 61, 2, 0) - 1;
    for (c3_i94 = 0; c3_i94 < 6; c3_i94++) {
      c3_h_b[c3_i94] = c3_sigma_meas[c3_i94 + 6 * c3_v_i] - c3_ymeas[c3_i94];
    }

    for (c3_i95 = 0; c3_i95 < 6; c3_i95++) {
      c3_i96 = 0;
      for (c3_i97 = 0; c3_i97 < 6; c3_i97++) {
        c3_i_b[c3_i96 + c3_i95] = c3_g_b[c3_i95] * c3_h_b[c3_i97];
        c3_i96 += 6;
      }
    }

    for (c3_i98 = 0; c3_i98 < 36; c3_i98++) {
      c3_Pzz[c3_i98] += c3_i_b[c3_i98];
    }

    c3_r_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_SCRIPT_FOR(0, 5, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 282);
  for (c3_i99 = 0; c3_i99 < 36; c3_i99++) {
    c3_Scov[c3_i99] = c3_Rcov[c3_i99] + c3_Pzz[c3_i99];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 283);
  for (c3_i100 = 0; c3_i100 < 180; c3_i100++) {
    c3_l_a[c3_i100] = c3_Pxz[c3_i100];
  }

  for (c3_i101 = 0; c3_i101 < 36; c3_i101++) {
    c3_b_Scov[c3_i101] = c3_Scov[c3_i101];
  }

  c3_inv(chartInstance, c3_b_Scov, c3_i_b);
  c3_b_eml_scalar_eg(chartInstance);
  c3_b_eml_scalar_eg(chartInstance);
  for (c3_i102 = 0; c3_i102 < 180; c3_i102++) {
    c3_k[c3_i102] = 0.0;
  }

  for (c3_i103 = 0; c3_i103 < 180; c3_i103++) {
    c3_k[c3_i103] = 0.0;
  }

  for (c3_i104 = 0; c3_i104 < 180; c3_i104++) {
    c3_dv13[c3_i104] = c3_l_a[c3_i104];
  }

  for (c3_i105 = 0; c3_i105 < 36; c3_i105++) {
    c3_dv14[c3_i105] = c3_i_b[c3_i105];
  }

  for (c3_i106 = 0; c3_i106 < 180; c3_i106++) {
    c3_dv15[c3_i106] = c3_dv13[c3_i106];
  }

  for (c3_i107 = 0; c3_i107 < 36; c3_i107++) {
    c3_dv16[c3_i107] = c3_dv14[c3_i107];
  }

  c3_c_eml_xgemm(chartInstance, c3_dv15, c3_dv16, c3_k);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 284);
  for (c3_i108 = 0; c3_i108 < 6; c3_i108++) {
    c3_innovation[c3_i108] = c3_measurement[c3_i108] - c3_ymeas[c3_i108];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 285);
  for (c3_i109 = 0; c3_i109 < 180; c3_i109++) {
    c3_l_a[c3_i109] = c3_k[c3_i109];
  }

  for (c3_i110 = 0; c3_i110 < 6; c3_i110++) {
    c3_g_b[c3_i110] = c3_innovation[c3_i110];
  }

  c3_c_eml_scalar_eg(chartInstance);
  c3_c_eml_scalar_eg(chartInstance);
  for (c3_i111 = 0; c3_i111 < 30; c3_i111++) {
    c3_correction[c3_i111] = 0.0;
  }

  for (c3_i112 = 0; c3_i112 < 30; c3_i112++) {
    c3_correction[c3_i112] = 0.0;
  }

  for (c3_i113 = 0; c3_i113 < 30; c3_i113++) {
    c3_C[c3_i113] = c3_correction[c3_i113];
  }

  for (c3_i114 = 0; c3_i114 < 30; c3_i114++) {
    c3_correction[c3_i114] = c3_C[c3_i114];
  }

  c3_c_threshold(chartInstance);
  for (c3_i115 = 0; c3_i115 < 30; c3_i115++) {
    c3_C[c3_i115] = c3_correction[c3_i115];
  }

  for (c3_i116 = 0; c3_i116 < 30; c3_i116++) {
    c3_correction[c3_i116] = c3_C[c3_i116];
  }

  for (c3_i117 = 0; c3_i117 < 30; c3_i117++) {
    c3_correction[c3_i117] = 0.0;
    c3_i118 = 0;
    for (c3_i119 = 0; c3_i119 < 6; c3_i119++) {
      c3_correction[c3_i117] += c3_l_a[c3_i118 + c3_i117] * c3_g_b[c3_i119];
      c3_i118 += 30;
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 286);
  for (c3_i120 = 0; c3_i120 < 30; c3_i120++) {
    chartInstance->c3_pStates[c3_i120] += c3_correction[c3_i120];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 287);
  for (c3_i121 = 0; c3_i121 < 900; c3_i121++) {
    c3_c_hoistedGlobal[c3_i121] = chartInstance->c3_pCovariance[c3_i121];
  }

  for (c3_i122 = 0; c3_i122 < 180; c3_i122++) {
    c3_l_a[c3_i122] = c3_k[c3_i122];
  }

  for (c3_i123 = 0; c3_i123 < 36; c3_i123++) {
    c3_i_b[c3_i123] = c3_Scov[c3_i123];
  }

  c3_b_eml_scalar_eg(chartInstance);
  c3_b_eml_scalar_eg(chartInstance);
  for (c3_i124 = 0; c3_i124 < 180; c3_i124++) {
    c3_qb_y[c3_i124] = 0.0;
  }

  for (c3_i125 = 0; c3_i125 < 180; c3_i125++) {
    c3_m_a[c3_i125] = c3_l_a[c3_i125];
  }

  for (c3_i126 = 0; c3_i126 < 36; c3_i126++) {
    c3_j_b[c3_i126] = c3_i_b[c3_i126];
  }

  c3_c_eml_xgemm(chartInstance, c3_m_a, c3_j_b, c3_qb_y);
  c3_i127 = 0;
  for (c3_i128 = 0; c3_i128 < 30; c3_i128++) {
    c3_i129 = 0;
    for (c3_i130 = 0; c3_i130 < 6; c3_i130++) {
      c3_k_b[c3_i130 + c3_i127] = c3_k[c3_i129 + c3_i128];
      c3_i129 += 30;
    }

    c3_i127 += 6;
  }

  c3_d_eml_scalar_eg(chartInstance);
  c3_d_eml_scalar_eg(chartInstance);
  for (c3_i131 = 0; c3_i131 < 900; c3_i131++) {
    c3_pb_y[c3_i131] = 0.0;
  }

  for (c3_i132 = 0; c3_i132 < 180; c3_i132++) {
    c3_rb_y[c3_i132] = c3_qb_y[c3_i132];
  }

  for (c3_i133 = 0; c3_i133 < 180; c3_i133++) {
    c3_l_b[c3_i133] = c3_k_b[c3_i133];
  }

  c3_d_eml_xgemm(chartInstance, c3_rb_y, c3_l_b, c3_pb_y);
  for (c3_i134 = 0; c3_i134 < 900; c3_i134++) {
    chartInstance->c3_pCovariance[c3_i134] = c3_c_hoistedGlobal[c3_i134] -
      c3_pb_y[c3_i134];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 289);
  for (c3_i135 = 0; c3_i135 < 36; c3_i135++) {
    c3_c_Scov[c3_i135] = c3_Scov[c3_i135];
  }

  c3_b_diag(chartInstance, c3_c_Scov, c3_dv17);
  for (c3_i136 = 0; c3_i136 < 6; c3_i136++) {
    c3_SCovOUT[c3_i136] = c3_dv17[c3_i136];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 290);
  for (c3_i137 = 0; c3_i137 < 30; c3_i137++) {
    c3_pStatesOUT[c3_i137] = chartInstance->c3_pStates[c3_i137];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 291);
  for (c3_i138 = 0; c3_i138 < 900; c3_i138++) {
    c3_dv18[c3_i138] = chartInstance->c3_pCovariance[c3_i138];
  }

  c3_c_diag(chartInstance, c3_dv18, c3_dv19);
  for (c3_i139 = 0; c3_i139 < 30; c3_i139++) {
    c3_pCovarianceOUT[c3_i139] = c3_dv19[c3_i139];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, -291);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber)
{
  (void)c3_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c3_chartNumber, c3_instanceNumber, 0U,
    sf_debug_get_script_id(
    "C:\\Users\\Rudaba\\Documents\\PhD Take 2\\Simulations\\AircraftModelling\\simulateUKF.m"));
}

static void c3_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_sprintf, const char_T *c3_identifier, char_T c3_y[14])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_sprintf), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_sprintf);
}

static void c3_b_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, char_T c3_y[14])
{
  char_T c3_cv0[14];
  int32_T c3_i140;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_cv0, 1, 10, 0U, 1, 0U, 2, 1,
                14);
  for (c3_i140 = 0; c3_i140 < 14; c3_i140++) {
    c3_y[c3_i140] = c3_cv0[c3_i140];
  }

  sf_mex_destroy(&c3_u);
}

static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i141;
  real_T c3_b_inData[6];
  int32_T c3_i142;
  real_T c3_u[6];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i141 = 0; c3_i141 < 6; c3_i141++) {
    c3_b_inData[c3_i141] = (*(real_T (*)[6])c3_inData)[c3_i141];
  }

  for (c3_i142 = 0; c3_i142 < 6; c3_i142++) {
    c3_u[c3_i142] = c3_b_inData[c3_i142];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_c_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_S, const char_T *c3_identifier, real_T c3_y[6])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_S), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_S);
}

static void c3_d_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[6])
{
  real_T c3_dv20[6];
  int32_T c3_i143;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv20, 1, 0, 0U, 1, 0U, 1, 6);
  for (c3_i143 = 0; c3_i143 < 6; c3_i143++) {
    c3_y[c3_i143] = c3_dv20[c3_i143];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_S;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[6];
  int32_T c3_i144;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_S = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_S), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_S);
  for (c3_i144 = 0; c3_i144 < 6; c3_i144++) {
    (*(real_T (*)[6])c3_outData)[c3_i144] = c3_y[c3_i144];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i145;
  real_T c3_b_inData[30];
  int32_T c3_i146;
  real_T c3_u[30];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i145 = 0; c3_i145 < 30; c3_i145++) {
    c3_b_inData[c3_i145] = (*(real_T (*)[30])c3_inData)[c3_i145];
  }

  for (c3_i146 = 0; c3_i146 < 30; c3_i146++) {
    c3_u[c3_i146] = c3_b_inData[c3_i146];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 30), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_e_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_cov_OUT, const char_T *c3_identifier, real_T c3_y[30])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_cov_OUT), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_cov_OUT);
}

static void c3_f_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[30])
{
  real_T c3_dv21[30];
  int32_T c3_i147;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv21, 1, 0, 0U, 1, 0U, 1, 30);
  for (c3_i147 = 0; c3_i147 < 30; c3_i147++) {
    c3_y[c3_i147] = c3_dv21[c3_i147];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_cov_OUT;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[30];
  int32_T c3_i148;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_cov_OUT = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_cov_OUT), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_cov_OUT);
  for (c3_i148 = 0; c3_i148 < 30; c3_i148++) {
    (*(real_T (*)[30])c3_outData)[c3_i148] = c3_y[c3_i148];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  real_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(real_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i149;
  real_T c3_b_inData[4];
  int32_T c3_i150;
  real_T c3_u[4];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i149 = 0; c3_i149 < 4; c3_i149++) {
    c3_b_inData[c3_i149] = (*(real_T (*)[4])c3_inData)[c3_i149];
  }

  for (c3_i150 = 0; c3_i150 < 4; c3_i150++) {
    c3_u[c3_i150] = c3_b_inData[c3_i150];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 4), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i151;
  real_T c3_b_inData[3];
  int32_T c3_i152;
  real_T c3_u[3];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i151 = 0; c3_i151 < 3; c3_i151++) {
    c3_b_inData[c3_i151] = (*(real_T (*)[3])c3_inData)[c3_i151];
  }

  for (c3_i152 = 0; c3_i152 < 3; c3_i152++) {
    c3_u[c3_i152] = c3_b_inData[c3_i152];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static real_T c3_g_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  real_T c3_y;
  real_T c3_d1;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_d1, 1, 0, 0U, 0, 0U, 0);
  c3_y = c3_d1;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_nargout;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_nargout = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_nargout), &c3_thisId);
  sf_mex_destroy(&c3_nargout);
  *(real_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_f_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i153;
  int32_T c3_i154;
  int32_T c3_i155;
  real_T c3_b_inData[183];
  int32_T c3_i156;
  int32_T c3_i157;
  int32_T c3_i158;
  real_T c3_u[183];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i153 = 0;
  for (c3_i154 = 0; c3_i154 < 61; c3_i154++) {
    for (c3_i155 = 0; c3_i155 < 3; c3_i155++) {
      c3_b_inData[c3_i155 + c3_i153] = (*(real_T (*)[183])c3_inData)[c3_i155 +
        c3_i153];
    }

    c3_i153 += 3;
  }

  c3_i156 = 0;
  for (c3_i157 = 0; c3_i157 < 61; c3_i157++) {
    for (c3_i158 = 0; c3_i158 < 3; c3_i158++) {
      c3_u[c3_i158 + c3_i156] = c3_b_inData[c3_i158 + c3_i156];
    }

    c3_i156 += 3;
  }

  c3_y = NULL;
  if (!chartInstance->c3_prevOmegaDots_not_empty) {
    sf_mex_assign(&c3_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0), false);
  } else {
    sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 3, 61),
                  false);
  }

  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_h_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_prevOmegaDots, const char_T *c3_identifier, real_T c3_y[183])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_prevOmegaDots),
                        &c3_thisId, c3_y);
  sf_mex_destroy(&c3_b_prevOmegaDots);
}

static void c3_i_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[183])
{
  real_T c3_dv22[183];
  int32_T c3_i159;
  if (mxIsEmpty(c3_u)) {
    chartInstance->c3_prevOmegaDots_not_empty = false;
  } else {
    chartInstance->c3_prevOmegaDots_not_empty = true;
    sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv22, 1, 0, 0U, 1, 0U, 2, 3,
                  61);
    for (c3_i159 = 0; c3_i159 < 183; c3_i159++) {
      c3_y[c3_i159] = c3_dv22[c3_i159];
    }
  }

  sf_mex_destroy(&c3_u);
}

static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_prevOmegaDots;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[183];
  int32_T c3_i160;
  int32_T c3_i161;
  int32_T c3_i162;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_b_prevOmegaDots = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_prevOmegaDots),
                        &c3_thisId, c3_y);
  sf_mex_destroy(&c3_b_prevOmegaDots);
  c3_i160 = 0;
  for (c3_i161 = 0; c3_i161 < 61; c3_i161++) {
    for (c3_i162 = 0; c3_i162 < 3; c3_i162++) {
      (*(real_T (*)[183])c3_outData)[c3_i162 + c3_i160] = c3_y[c3_i162 + c3_i160];
    }

    c3_i160 += 3;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_g_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i163;
  real_T c3_b_inData[30];
  int32_T c3_i164;
  real_T c3_u[30];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i163 = 0; c3_i163 < 30; c3_i163++) {
    c3_b_inData[c3_i163] = (*(real_T (*)[30])c3_inData)[c3_i163];
  }

  for (c3_i164 = 0; c3_i164 < 30; c3_i164++) {
    c3_u[c3_i164] = c3_b_inData[c3_i164];
  }

  c3_y = NULL;
  if (!chartInstance->c3_covDiag_not_empty) {
    sf_mex_assign(&c3_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0), false);
  } else {
    sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 30), false);
  }

  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_j_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_covDiag, const char_T *c3_identifier, real_T c3_y[30])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_covDiag), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_b_covDiag);
}

static void c3_k_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[30])
{
  real_T c3_dv23[30];
  int32_T c3_i165;
  if (mxIsEmpty(c3_u)) {
    chartInstance->c3_covDiag_not_empty = false;
  } else {
    chartInstance->c3_covDiag_not_empty = true;
    sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv23, 1, 0, 0U, 1, 0U, 1, 30);
    for (c3_i165 = 0; c3_i165 < 30; c3_i165++) {
      c3_y[c3_i165] = c3_dv23[c3_i165];
    }
  }

  sf_mex_destroy(&c3_u);
}

static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_covDiag;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[30];
  int32_T c3_i166;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_b_covDiag = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_covDiag), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_b_covDiag);
  for (c3_i166 = 0; c3_i166 < 30; c3_i166++) {
    (*(real_T (*)[30])c3_outData)[c3_i166] = c3_y[c3_i166];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_h_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i167;
  int32_T c3_i168;
  int32_T c3_i169;
  real_T c3_b_inData[900];
  int32_T c3_i170;
  int32_T c3_i171;
  int32_T c3_i172;
  real_T c3_u[900];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i167 = 0;
  for (c3_i168 = 0; c3_i168 < 30; c3_i168++) {
    for (c3_i169 = 0; c3_i169 < 30; c3_i169++) {
      c3_b_inData[c3_i169 + c3_i167] = (*(real_T (*)[900])c3_inData)[c3_i169 +
        c3_i167];
    }

    c3_i167 += 30;
  }

  c3_i170 = 0;
  for (c3_i171 = 0; c3_i171 < 30; c3_i171++) {
    for (c3_i172 = 0; c3_i172 < 30; c3_i172++) {
      c3_u[c3_i172 + c3_i170] = c3_b_inData[c3_i172 + c3_i170];
    }

    c3_i170 += 30;
  }

  c3_y = NULL;
  if (!chartInstance->c3_pCovariance_not_empty) {
    sf_mex_assign(&c3_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0), false);
  } else {
    sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 30, 30),
                  false);
  }

  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_l_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_pCovariance, const char_T *c3_identifier, real_T c3_y[900])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_pCovariance), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_b_pCovariance);
}

static void c3_m_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[900])
{
  real_T c3_dv24[900];
  int32_T c3_i173;
  if (mxIsEmpty(c3_u)) {
    chartInstance->c3_pCovariance_not_empty = false;
  } else {
    chartInstance->c3_pCovariance_not_empty = true;
    sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv24, 1, 0, 0U, 1, 0U, 2, 30,
                  30);
    for (c3_i173 = 0; c3_i173 < 900; c3_i173++) {
      c3_y[c3_i173] = c3_dv24[c3_i173];
    }
  }

  sf_mex_destroy(&c3_u);
}

static void c3_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_pCovariance;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[900];
  int32_T c3_i174;
  int32_T c3_i175;
  int32_T c3_i176;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_b_pCovariance = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_pCovariance), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_b_pCovariance);
  c3_i174 = 0;
  for (c3_i175 = 0; c3_i175 < 30; c3_i175++) {
    for (c3_i176 = 0; c3_i176 < 30; c3_i176++) {
      (*(real_T (*)[900])c3_outData)[c3_i176 + c3_i174] = c3_y[c3_i176 + c3_i174];
    }

    c3_i174 += 30;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_i_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i177;
  real_T c3_b_inData[30];
  int32_T c3_i178;
  real_T c3_u[30];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i177 = 0; c3_i177 < 30; c3_i177++) {
    c3_b_inData[c3_i177] = (*(real_T (*)[30])c3_inData)[c3_i177];
  }

  for (c3_i178 = 0; c3_i178 < 30; c3_i178++) {
    c3_u[c3_i178] = c3_b_inData[c3_i178];
  }

  c3_y = NULL;
  if (!chartInstance->c3_pStates_not_empty) {
    sf_mex_assign(&c3_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0), false);
  } else {
    sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 30), false);
  }

  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_n_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_pStates, const char_T *c3_identifier, real_T c3_y[30])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_pStates), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_b_pStates);
}

static void c3_o_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[30])
{
  real_T c3_dv25[30];
  int32_T c3_i179;
  if (mxIsEmpty(c3_u)) {
    chartInstance->c3_pStates_not_empty = false;
  } else {
    chartInstance->c3_pStates_not_empty = true;
    sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv25, 1, 0, 0U, 1, 0U, 1, 30);
    for (c3_i179 = 0; c3_i179 < 30; c3_i179++) {
      c3_y[c3_i179] = c3_dv25[c3_i179];
    }
  }

  sf_mex_destroy(&c3_u);
}

static void c3_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_pStates;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[30];
  int32_T c3_i180;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_b_pStates = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_pStates), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_b_pStates);
  for (c3_i180 = 0; c3_i180 < 30; c3_i180++) {
    (*(real_T (*)[30])c3_outData)[c3_i180] = c3_y[c3_i180];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static void c3_p_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[3])
{
  real_T c3_dv26[3];
  int32_T c3_i181;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv26, 1, 0, 0U, 1, 0U, 1, 3);
  for (c3_i181 = 0; c3_i181 < 3; c3_i181++) {
    c3_y[c3_i181] = c3_dv26[c3_i181];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_omega;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[3];
  int32_T c3_i182;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_omega = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_p_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_omega), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_omega);
  for (c3_i182 = 0; c3_i182 < 3; c3_i182++) {
    (*(real_T (*)[3])c3_outData)[c3_i182] = c3_y[c3_i182];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static void c3_q_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[4])
{
  real_T c3_dv27[4];
  int32_T c3_i183;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv27, 1, 0, 0U, 1, 0U, 1, 4);
  for (c3_i183 = 0; c3_i183 < 4; c3_i183++) {
    c3_y[c3_i183] = c3_dv27[c3_i183];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_u;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[4];
  int32_T c3_i184;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_u = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_q_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_u), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_u);
  for (c3_i184 = 0; c3_i184 < 4; c3_i184++) {
    (*(real_T (*)[4])c3_outData)[c3_i184] = c3_y[c3_i184];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_j_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i185;
  int32_T c3_i186;
  int32_T c3_i187;
  real_T c3_b_inData[180];
  int32_T c3_i188;
  int32_T c3_i189;
  int32_T c3_i190;
  real_T c3_u[180];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i185 = 0;
  for (c3_i186 = 0; c3_i186 < 6; c3_i186++) {
    for (c3_i187 = 0; c3_i187 < 30; c3_i187++) {
      c3_b_inData[c3_i187 + c3_i185] = (*(real_T (*)[180])c3_inData)[c3_i187 +
        c3_i185];
    }

    c3_i185 += 30;
  }

  c3_i188 = 0;
  for (c3_i189 = 0; c3_i189 < 6; c3_i189++) {
    for (c3_i190 = 0; c3_i190 < 30; c3_i190++) {
      c3_u[c3_i190 + c3_i188] = c3_b_inData[c3_i190 + c3_i188];
    }

    c3_i188 += 30;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 30, 6), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_r_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[180])
{
  real_T c3_dv28[180];
  int32_T c3_i191;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv28, 1, 0, 0U, 1, 0U, 2, 30,
                6);
  for (c3_i191 = 0; c3_i191 < 180; c3_i191++) {
    c3_y[c3_i191] = c3_dv28[c3_i191];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_k;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[180];
  int32_T c3_i192;
  int32_T c3_i193;
  int32_T c3_i194;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_k = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_r_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_k), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_k);
  c3_i192 = 0;
  for (c3_i193 = 0; c3_i193 < 6; c3_i193++) {
    for (c3_i194 = 0; c3_i194 < 30; c3_i194++) {
      (*(real_T (*)[180])c3_outData)[c3_i194 + c3_i192] = c3_y[c3_i194 + c3_i192];
    }

    c3_i192 += 30;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_k_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i195;
  int32_T c3_i196;
  int32_T c3_i197;
  real_T c3_b_inData[36];
  int32_T c3_i198;
  int32_T c3_i199;
  int32_T c3_i200;
  real_T c3_u[36];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i195 = 0;
  for (c3_i196 = 0; c3_i196 < 6; c3_i196++) {
    for (c3_i197 = 0; c3_i197 < 6; c3_i197++) {
      c3_b_inData[c3_i197 + c3_i195] = (*(real_T (*)[36])c3_inData)[c3_i197 +
        c3_i195];
    }

    c3_i195 += 6;
  }

  c3_i198 = 0;
  for (c3_i199 = 0; c3_i199 < 6; c3_i199++) {
    for (c3_i200 = 0; c3_i200 < 6; c3_i200++) {
      c3_u[c3_i200 + c3_i198] = c3_b_inData[c3_i200 + c3_i198];
    }

    c3_i198 += 6;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 6, 6), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_s_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[36])
{
  real_T c3_dv29[36];
  int32_T c3_i201;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv29, 1, 0, 0U, 1, 0U, 2, 6, 6);
  for (c3_i201 = 0; c3_i201 < 36; c3_i201++) {
    c3_y[c3_i201] = c3_dv29[c3_i201];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_Scov;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[36];
  int32_T c3_i202;
  int32_T c3_i203;
  int32_T c3_i204;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_Scov = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_s_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_Scov), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_Scov);
  c3_i202 = 0;
  for (c3_i203 = 0; c3_i203 < 6; c3_i203++) {
    for (c3_i204 = 0; c3_i204 < 6; c3_i204++) {
      (*(real_T (*)[36])c3_outData)[c3_i204 + c3_i202] = c3_y[c3_i204 + c3_i202];
    }

    c3_i202 += 6;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_l_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i205;
  int32_T c3_i206;
  int32_T c3_i207;
  real_T c3_b_inData[900];
  int32_T c3_i208;
  int32_T c3_i209;
  int32_T c3_i210;
  real_T c3_u[900];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i205 = 0;
  for (c3_i206 = 0; c3_i206 < 30; c3_i206++) {
    for (c3_i207 = 0; c3_i207 < 30; c3_i207++) {
      c3_b_inData[c3_i207 + c3_i205] = (*(real_T (*)[900])c3_inData)[c3_i207 +
        c3_i205];
    }

    c3_i205 += 30;
  }

  c3_i208 = 0;
  for (c3_i209 = 0; c3_i209 < 30; c3_i209++) {
    for (c3_i210 = 0; c3_i210 < 30; c3_i210++) {
      c3_u[c3_i210 + c3_i208] = c3_b_inData[c3_i210 + c3_i208];
    }

    c3_i208 += 30;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 30, 30), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_t_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[900])
{
  real_T c3_dv30[900];
  int32_T c3_i211;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv30, 1, 0, 0U, 1, 0U, 2, 30,
                30);
  for (c3_i211 = 0; c3_i211 < 900; c3_i211++) {
    c3_y[c3_i211] = c3_dv30[c3_i211];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_l_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_matrixSQRT;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[900];
  int32_T c3_i212;
  int32_T c3_i213;
  int32_T c3_i214;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_matrixSQRT = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_t_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_matrixSQRT), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_matrixSQRT);
  c3_i212 = 0;
  for (c3_i213 = 0; c3_i213 < 30; c3_i213++) {
    for (c3_i214 = 0; c3_i214 < 30; c3_i214++) {
      (*(real_T (*)[900])c3_outData)[c3_i214 + c3_i212] = c3_y[c3_i214 + c3_i212];
    }

    c3_i212 += 30;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_m_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i215;
  real_T c3_b_inData[61];
  int32_T c3_i216;
  real_T c3_u[61];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i215 = 0; c3_i215 < 61; c3_i215++) {
    c3_b_inData[c3_i215] = (*(real_T (*)[61])c3_inData)[c3_i215];
  }

  for (c3_i216 = 0; c3_i216 < 61; c3_i216++) {
    c3_u[c3_i216] = c3_b_inData[c3_i216];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 61), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_u_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[61])
{
  real_T c3_dv31[61];
  int32_T c3_i217;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv31, 1, 0, 0U, 1, 0U, 1, 61);
  for (c3_i217 = 0; c3_i217 < 61; c3_i217++) {
    c3_y[c3_i217] = c3_dv31[c3_i217];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_m_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_W;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[61];
  int32_T c3_i218;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_W = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_u_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_W), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_W);
  for (c3_i218 = 0; c3_i218 < 61; c3_i218++) {
    (*(real_T (*)[61])c3_outData)[c3_i218] = c3_y[c3_i218];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_n_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i219;
  int32_T c3_i220;
  int32_T c3_i221;
  real_T c3_b_inData[366];
  int32_T c3_i222;
  int32_T c3_i223;
  int32_T c3_i224;
  real_T c3_u[366];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i219 = 0;
  for (c3_i220 = 0; c3_i220 < 61; c3_i220++) {
    for (c3_i221 = 0; c3_i221 < 6; c3_i221++) {
      c3_b_inData[c3_i221 + c3_i219] = (*(real_T (*)[366])c3_inData)[c3_i221 +
        c3_i219];
    }

    c3_i219 += 6;
  }

  c3_i222 = 0;
  for (c3_i223 = 0; c3_i223 < 61; c3_i223++) {
    for (c3_i224 = 0; c3_i224 < 6; c3_i224++) {
      c3_u[c3_i224 + c3_i222] = c3_b_inData[c3_i224 + c3_i222];
    }

    c3_i222 += 6;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 6, 61), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_v_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[366])
{
  real_T c3_dv32[366];
  int32_T c3_i225;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv32, 1, 0, 0U, 1, 0U, 2, 6,
                61);
  for (c3_i225 = 0; c3_i225 < 366; c3_i225++) {
    c3_y[c3_i225] = c3_dv32[c3_i225];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_n_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_sigma_meas;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[366];
  int32_T c3_i226;
  int32_T c3_i227;
  int32_T c3_i228;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_sigma_meas = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_v_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_sigma_meas), &c3_thisId,
                        c3_y);
  sf_mex_destroy(&c3_sigma_meas);
  c3_i226 = 0;
  for (c3_i227 = 0; c3_i227 < 61; c3_i227++) {
    for (c3_i228 = 0; c3_i228 < 6; c3_i228++) {
      (*(real_T (*)[366])c3_outData)[c3_i228 + c3_i226] = c3_y[c3_i228 + c3_i226];
    }

    c3_i226 += 6;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_o_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i229;
  int32_T c3_i230;
  int32_T c3_i231;
  real_T c3_b_inData[1830];
  int32_T c3_i232;
  int32_T c3_i233;
  int32_T c3_i234;
  real_T c3_u[1830];
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i229 = 0;
  for (c3_i230 = 0; c3_i230 < 61; c3_i230++) {
    for (c3_i231 = 0; c3_i231 < 30; c3_i231++) {
      c3_b_inData[c3_i231 + c3_i229] = (*(real_T (*)[1830])c3_inData)[c3_i231 +
        c3_i229];
    }

    c3_i229 += 30;
  }

  c3_i232 = 0;
  for (c3_i233 = 0; c3_i233 < 61; c3_i233++) {
    for (c3_i234 = 0; c3_i234 < 30; c3_i234++) {
      c3_u[c3_i234 + c3_i232] = c3_b_inData[c3_i234 + c3_i232];
    }

    c3_i232 += 30;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 30, 61), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static void c3_w_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[1830])
{
  real_T c3_dv33[1830];
  int32_T c3_i235;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv33, 1, 0, 0U, 1, 0U, 2, 30,
                61);
  for (c3_i235 = 0; c3_i235 < 1830; c3_i235++) {
    c3_y[c3_i235] = c3_dv33[c3_i235];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_o_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_sigma;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[1830];
  int32_T c3_i236;
  int32_T c3_i237;
  int32_T c3_i238;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_sigma = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_w_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_sigma), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_sigma);
  c3_i236 = 0;
  for (c3_i237 = 0; c3_i237 < 61; c3_i237++) {
    for (c3_i238 = 0; c3_i238 < 30; c3_i238++) {
      (*(real_T (*)[1830])c3_outData)[c3_i238 + c3_i236] = c3_y[c3_i238 +
        c3_i236];
    }

    c3_i236 += 30;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

const mxArray
  *sf_c3_aircraftControl_FullStateFilters_get_eml_resolved_functions_info(void)
{
  const mxArray *c3_nameCaptureInfo = NULL;
  c3_nameCaptureInfo = NULL;
  sf_mex_assign(&c3_nameCaptureInfo, sf_mex_createstruct("structure", 2, 255, 1),
                false);
  c3_info_helper(&c3_nameCaptureInfo);
  c3_b_info_helper(&c3_nameCaptureInfo);
  c3_c_info_helper(&c3_nameCaptureInfo);
  c3_d_info_helper(&c3_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c3_nameCaptureInfo);
  return c3_nameCaptureInfo;
}

static void c3_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs0 = NULL;
  const mxArray *c3_lhs0 = NULL;
  const mxArray *c3_rhs1 = NULL;
  const mxArray *c3_lhs1 = NULL;
  const mxArray *c3_rhs2 = NULL;
  const mxArray *c3_lhs2 = NULL;
  const mxArray *c3_rhs3 = NULL;
  const mxArray *c3_lhs3 = NULL;
  const mxArray *c3_rhs4 = NULL;
  const mxArray *c3_lhs4 = NULL;
  const mxArray *c3_rhs5 = NULL;
  const mxArray *c3_lhs5 = NULL;
  const mxArray *c3_rhs6 = NULL;
  const mxArray *c3_lhs6 = NULL;
  const mxArray *c3_rhs7 = NULL;
  const mxArray *c3_lhs7 = NULL;
  const mxArray *c3_rhs8 = NULL;
  const mxArray *c3_lhs8 = NULL;
  const mxArray *c3_rhs9 = NULL;
  const mxArray *c3_lhs9 = NULL;
  const mxArray *c3_rhs10 = NULL;
  const mxArray *c3_lhs10 = NULL;
  const mxArray *c3_rhs11 = NULL;
  const mxArray *c3_lhs11 = NULL;
  const mxArray *c3_rhs12 = NULL;
  const mxArray *c3_lhs12 = NULL;
  const mxArray *c3_rhs13 = NULL;
  const mxArray *c3_lhs13 = NULL;
  const mxArray *c3_rhs14 = NULL;
  const mxArray *c3_lhs14 = NULL;
  const mxArray *c3_rhs15 = NULL;
  const mxArray *c3_lhs15 = NULL;
  const mxArray *c3_rhs16 = NULL;
  const mxArray *c3_lhs16 = NULL;
  const mxArray *c3_rhs17 = NULL;
  const mxArray *c3_lhs17 = NULL;
  const mxArray *c3_rhs18 = NULL;
  const mxArray *c3_lhs18 = NULL;
  const mxArray *c3_rhs19 = NULL;
  const mxArray *c3_lhs19 = NULL;
  const mxArray *c3_rhs20 = NULL;
  const mxArray *c3_lhs20 = NULL;
  const mxArray *c3_rhs21 = NULL;
  const mxArray *c3_lhs21 = NULL;
  const mxArray *c3_rhs22 = NULL;
  const mxArray *c3_lhs22 = NULL;
  const mxArray *c3_rhs23 = NULL;
  const mxArray *c3_lhs23 = NULL;
  const mxArray *c3_rhs24 = NULL;
  const mxArray *c3_lhs24 = NULL;
  const mxArray *c3_rhs25 = NULL;
  const mxArray *c3_lhs25 = NULL;
  const mxArray *c3_rhs26 = NULL;
  const mxArray *c3_lhs26 = NULL;
  const mxArray *c3_rhs27 = NULL;
  const mxArray *c3_lhs27 = NULL;
  const mxArray *c3_rhs28 = NULL;
  const mxArray *c3_lhs28 = NULL;
  const mxArray *c3_rhs29 = NULL;
  const mxArray *c3_lhs29 = NULL;
  const mxArray *c3_rhs30 = NULL;
  const mxArray *c3_lhs30 = NULL;
  const mxArray *c3_rhs31 = NULL;
  const mxArray *c3_lhs31 = NULL;
  const mxArray *c3_rhs32 = NULL;
  const mxArray *c3_lhs32 = NULL;
  const mxArray *c3_rhs33 = NULL;
  const mxArray *c3_lhs33 = NULL;
  const mxArray *c3_rhs34 = NULL;
  const mxArray *c3_lhs34 = NULL;
  const mxArray *c3_rhs35 = NULL;
  const mxArray *c3_lhs35 = NULL;
  const mxArray *c3_rhs36 = NULL;
  const mxArray *c3_lhs36 = NULL;
  const mxArray *c3_rhs37 = NULL;
  const mxArray *c3_lhs37 = NULL;
  const mxArray *c3_rhs38 = NULL;
  const mxArray *c3_lhs38 = NULL;
  const mxArray *c3_rhs39 = NULL;
  const mxArray *c3_lhs39 = NULL;
  const mxArray *c3_rhs40 = NULL;
  const mxArray *c3_lhs40 = NULL;
  const mxArray *c3_rhs41 = NULL;
  const mxArray *c3_lhs41 = NULL;
  const mxArray *c3_rhs42 = NULL;
  const mxArray *c3_lhs42 = NULL;
  const mxArray *c3_rhs43 = NULL;
  const mxArray *c3_lhs43 = NULL;
  const mxArray *c3_rhs44 = NULL;
  const mxArray *c3_lhs44 = NULL;
  const mxArray *c3_rhs45 = NULL;
  const mxArray *c3_lhs45 = NULL;
  const mxArray *c3_rhs46 = NULL;
  const mxArray *c3_lhs46 = NULL;
  const mxArray *c3_rhs47 = NULL;
  const mxArray *c3_lhs47 = NULL;
  const mxArray *c3_rhs48 = NULL;
  const mxArray *c3_lhs48 = NULL;
  const mxArray *c3_rhs49 = NULL;
  const mxArray *c3_lhs49 = NULL;
  const mxArray *c3_rhs50 = NULL;
  const mxArray *c3_lhs50 = NULL;
  const mxArray *c3_rhs51 = NULL;
  const mxArray *c3_lhs51 = NULL;
  const mxArray *c3_rhs52 = NULL;
  const mxArray *c3_lhs52 = NULL;
  const mxArray *c3_rhs53 = NULL;
  const mxArray *c3_lhs53 = NULL;
  const mxArray *c3_rhs54 = NULL;
  const mxArray *c3_lhs54 = NULL;
  const mxArray *c3_rhs55 = NULL;
  const mxArray *c3_lhs55 = NULL;
  const mxArray *c3_rhs56 = NULL;
  const mxArray *c3_lhs56 = NULL;
  const mxArray *c3_rhs57 = NULL;
  const mxArray *c3_lhs57 = NULL;
  const mxArray *c3_rhs58 = NULL;
  const mxArray *c3_lhs58 = NULL;
  const mxArray *c3_rhs59 = NULL;
  const mxArray *c3_lhs59 = NULL;
  const mxArray *c3_rhs60 = NULL;
  const mxArray *c3_lhs60 = NULL;
  const mxArray *c3_rhs61 = NULL;
  const mxArray *c3_lhs61 = NULL;
  const mxArray *c3_rhs62 = NULL;
  const mxArray *c3_lhs62 = NULL;
  const mxArray *c3_rhs63 = NULL;
  const mxArray *c3_lhs63 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("simulateUKF"), "name", "name",
                  0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1435884894U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c3_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mpower"), "name", "name", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677878U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c3_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c3_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("ismatrix"), "name", "name", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1331268858U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c3_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("power"), "name", "name", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677880U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c3_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c3_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c3_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c3_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c3_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c3_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("floor"), "name", "name", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677854U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c3_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c3_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786326U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c3_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c3_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1383841294U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c3_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c3_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("diag"), "name", "name", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677854U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c3_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("ismatrix"), "name", "name", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1331268858U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c3_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c3_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c3_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c3_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c3_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381817900U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c3_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eye"), "name", "name", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381817898U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c3_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1368154230U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c3_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c3_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("isinf"), "name", "name", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677856U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c3_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c3_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786382U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c3_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c3_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmin"), "name", "name", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c3_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381817900U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c3_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326692322U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c3_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381817900U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c3_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326691996U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c3_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmin"), "name", "name", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c3_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c3_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c3_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c3_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("power"), "name", "name", 39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677880U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c3_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("mrdivide"), "name", "name", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1388424096U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1369981086U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c3_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c3_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("rdivide"), "name", "name", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677880U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c3_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c3_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786396U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c3_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c3_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c3_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("chol"), "name", "name", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1344443234U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c3_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c3_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("ismatrix"), "name", "name", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1331268858U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c3_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c3_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_error"), "name", "name",
                  51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1343801558U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c3_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xpotrf"), "name", "name",
                  52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xpotrf.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786408U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c3_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xpotrf.m"),
                  "context", "context", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_lapack_xpotrf"), "name",
                  "name", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xpotrf.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786412U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c3_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xpotrf.m"),
                  "context", "context", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_matlab_zpotrf"), "name",
                  "name", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786424U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c3_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c3_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c3_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c3_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c3_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c3_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c3_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c3_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xdotc"), "name", "name",
                  62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951890U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c3_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c3_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c3_rhs0);
  sf_mex_destroy(&c3_lhs0);
  sf_mex_destroy(&c3_rhs1);
  sf_mex_destroy(&c3_lhs1);
  sf_mex_destroy(&c3_rhs2);
  sf_mex_destroy(&c3_lhs2);
  sf_mex_destroy(&c3_rhs3);
  sf_mex_destroy(&c3_lhs3);
  sf_mex_destroy(&c3_rhs4);
  sf_mex_destroy(&c3_lhs4);
  sf_mex_destroy(&c3_rhs5);
  sf_mex_destroy(&c3_lhs5);
  sf_mex_destroy(&c3_rhs6);
  sf_mex_destroy(&c3_lhs6);
  sf_mex_destroy(&c3_rhs7);
  sf_mex_destroy(&c3_lhs7);
  sf_mex_destroy(&c3_rhs8);
  sf_mex_destroy(&c3_lhs8);
  sf_mex_destroy(&c3_rhs9);
  sf_mex_destroy(&c3_lhs9);
  sf_mex_destroy(&c3_rhs10);
  sf_mex_destroy(&c3_lhs10);
  sf_mex_destroy(&c3_rhs11);
  sf_mex_destroy(&c3_lhs11);
  sf_mex_destroy(&c3_rhs12);
  sf_mex_destroy(&c3_lhs12);
  sf_mex_destroy(&c3_rhs13);
  sf_mex_destroy(&c3_lhs13);
  sf_mex_destroy(&c3_rhs14);
  sf_mex_destroy(&c3_lhs14);
  sf_mex_destroy(&c3_rhs15);
  sf_mex_destroy(&c3_lhs15);
  sf_mex_destroy(&c3_rhs16);
  sf_mex_destroy(&c3_lhs16);
  sf_mex_destroy(&c3_rhs17);
  sf_mex_destroy(&c3_lhs17);
  sf_mex_destroy(&c3_rhs18);
  sf_mex_destroy(&c3_lhs18);
  sf_mex_destroy(&c3_rhs19);
  sf_mex_destroy(&c3_lhs19);
  sf_mex_destroy(&c3_rhs20);
  sf_mex_destroy(&c3_lhs20);
  sf_mex_destroy(&c3_rhs21);
  sf_mex_destroy(&c3_lhs21);
  sf_mex_destroy(&c3_rhs22);
  sf_mex_destroy(&c3_lhs22);
  sf_mex_destroy(&c3_rhs23);
  sf_mex_destroy(&c3_lhs23);
  sf_mex_destroy(&c3_rhs24);
  sf_mex_destroy(&c3_lhs24);
  sf_mex_destroy(&c3_rhs25);
  sf_mex_destroy(&c3_lhs25);
  sf_mex_destroy(&c3_rhs26);
  sf_mex_destroy(&c3_lhs26);
  sf_mex_destroy(&c3_rhs27);
  sf_mex_destroy(&c3_lhs27);
  sf_mex_destroy(&c3_rhs28);
  sf_mex_destroy(&c3_lhs28);
  sf_mex_destroy(&c3_rhs29);
  sf_mex_destroy(&c3_lhs29);
  sf_mex_destroy(&c3_rhs30);
  sf_mex_destroy(&c3_lhs30);
  sf_mex_destroy(&c3_rhs31);
  sf_mex_destroy(&c3_lhs31);
  sf_mex_destroy(&c3_rhs32);
  sf_mex_destroy(&c3_lhs32);
  sf_mex_destroy(&c3_rhs33);
  sf_mex_destroy(&c3_lhs33);
  sf_mex_destroy(&c3_rhs34);
  sf_mex_destroy(&c3_lhs34);
  sf_mex_destroy(&c3_rhs35);
  sf_mex_destroy(&c3_lhs35);
  sf_mex_destroy(&c3_rhs36);
  sf_mex_destroy(&c3_lhs36);
  sf_mex_destroy(&c3_rhs37);
  sf_mex_destroy(&c3_lhs37);
  sf_mex_destroy(&c3_rhs38);
  sf_mex_destroy(&c3_lhs38);
  sf_mex_destroy(&c3_rhs39);
  sf_mex_destroy(&c3_lhs39);
  sf_mex_destroy(&c3_rhs40);
  sf_mex_destroy(&c3_lhs40);
  sf_mex_destroy(&c3_rhs41);
  sf_mex_destroy(&c3_lhs41);
  sf_mex_destroy(&c3_rhs42);
  sf_mex_destroy(&c3_lhs42);
  sf_mex_destroy(&c3_rhs43);
  sf_mex_destroy(&c3_lhs43);
  sf_mex_destroy(&c3_rhs44);
  sf_mex_destroy(&c3_lhs44);
  sf_mex_destroy(&c3_rhs45);
  sf_mex_destroy(&c3_lhs45);
  sf_mex_destroy(&c3_rhs46);
  sf_mex_destroy(&c3_lhs46);
  sf_mex_destroy(&c3_rhs47);
  sf_mex_destroy(&c3_lhs47);
  sf_mex_destroy(&c3_rhs48);
  sf_mex_destroy(&c3_lhs48);
  sf_mex_destroy(&c3_rhs49);
  sf_mex_destroy(&c3_lhs49);
  sf_mex_destroy(&c3_rhs50);
  sf_mex_destroy(&c3_lhs50);
  sf_mex_destroy(&c3_rhs51);
  sf_mex_destroy(&c3_lhs51);
  sf_mex_destroy(&c3_rhs52);
  sf_mex_destroy(&c3_lhs52);
  sf_mex_destroy(&c3_rhs53);
  sf_mex_destroy(&c3_lhs53);
  sf_mex_destroy(&c3_rhs54);
  sf_mex_destroy(&c3_lhs54);
  sf_mex_destroy(&c3_rhs55);
  sf_mex_destroy(&c3_lhs55);
  sf_mex_destroy(&c3_rhs56);
  sf_mex_destroy(&c3_lhs56);
  sf_mex_destroy(&c3_rhs57);
  sf_mex_destroy(&c3_lhs57);
  sf_mex_destroy(&c3_rhs58);
  sf_mex_destroy(&c3_lhs58);
  sf_mex_destroy(&c3_rhs59);
  sf_mex_destroy(&c3_lhs59);
  sf_mex_destroy(&c3_rhs60);
  sf_mex_destroy(&c3_lhs60);
  sf_mex_destroy(&c3_rhs61);
  sf_mex_destroy(&c3_lhs61);
  sf_mex_destroy(&c3_rhs62);
  sf_mex_destroy(&c3_lhs62);
  sf_mex_destroy(&c3_rhs63);
  sf_mex_destroy(&c3_lhs63);
}

static const mxArray *c3_emlrt_marshallOut(const char * c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c3_u)), false);
  return c3_y;
}

static const mxArray *c3_b_emlrt_marshallOut(const uint32_T c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 7, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static void c3_b_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs64 = NULL;
  const mxArray *c3_lhs64 = NULL;
  const mxArray *c3_rhs65 = NULL;
  const mxArray *c3_lhs65 = NULL;
  const mxArray *c3_rhs66 = NULL;
  const mxArray *c3_lhs66 = NULL;
  const mxArray *c3_rhs67 = NULL;
  const mxArray *c3_lhs67 = NULL;
  const mxArray *c3_rhs68 = NULL;
  const mxArray *c3_lhs68 = NULL;
  const mxArray *c3_rhs69 = NULL;
  const mxArray *c3_lhs69 = NULL;
  const mxArray *c3_rhs70 = NULL;
  const mxArray *c3_lhs70 = NULL;
  const mxArray *c3_rhs71 = NULL;
  const mxArray *c3_lhs71 = NULL;
  const mxArray *c3_rhs72 = NULL;
  const mxArray *c3_lhs72 = NULL;
  const mxArray *c3_rhs73 = NULL;
  const mxArray *c3_lhs73 = NULL;
  const mxArray *c3_rhs74 = NULL;
  const mxArray *c3_lhs74 = NULL;
  const mxArray *c3_rhs75 = NULL;
  const mxArray *c3_lhs75 = NULL;
  const mxArray *c3_rhs76 = NULL;
  const mxArray *c3_lhs76 = NULL;
  const mxArray *c3_rhs77 = NULL;
  const mxArray *c3_lhs77 = NULL;
  const mxArray *c3_rhs78 = NULL;
  const mxArray *c3_lhs78 = NULL;
  const mxArray *c3_rhs79 = NULL;
  const mxArray *c3_lhs79 = NULL;
  const mxArray *c3_rhs80 = NULL;
  const mxArray *c3_lhs80 = NULL;
  const mxArray *c3_rhs81 = NULL;
  const mxArray *c3_lhs81 = NULL;
  const mxArray *c3_rhs82 = NULL;
  const mxArray *c3_lhs82 = NULL;
  const mxArray *c3_rhs83 = NULL;
  const mxArray *c3_lhs83 = NULL;
  const mxArray *c3_rhs84 = NULL;
  const mxArray *c3_lhs84 = NULL;
  const mxArray *c3_rhs85 = NULL;
  const mxArray *c3_lhs85 = NULL;
  const mxArray *c3_rhs86 = NULL;
  const mxArray *c3_lhs86 = NULL;
  const mxArray *c3_rhs87 = NULL;
  const mxArray *c3_lhs87 = NULL;
  const mxArray *c3_rhs88 = NULL;
  const mxArray *c3_lhs88 = NULL;
  const mxArray *c3_rhs89 = NULL;
  const mxArray *c3_lhs89 = NULL;
  const mxArray *c3_rhs90 = NULL;
  const mxArray *c3_lhs90 = NULL;
  const mxArray *c3_rhs91 = NULL;
  const mxArray *c3_lhs91 = NULL;
  const mxArray *c3_rhs92 = NULL;
  const mxArray *c3_lhs92 = NULL;
  const mxArray *c3_rhs93 = NULL;
  const mxArray *c3_lhs93 = NULL;
  const mxArray *c3_rhs94 = NULL;
  const mxArray *c3_lhs94 = NULL;
  const mxArray *c3_rhs95 = NULL;
  const mxArray *c3_lhs95 = NULL;
  const mxArray *c3_rhs96 = NULL;
  const mxArray *c3_lhs96 = NULL;
  const mxArray *c3_rhs97 = NULL;
  const mxArray *c3_lhs97 = NULL;
  const mxArray *c3_rhs98 = NULL;
  const mxArray *c3_lhs98 = NULL;
  const mxArray *c3_rhs99 = NULL;
  const mxArray *c3_lhs99 = NULL;
  const mxArray *c3_rhs100 = NULL;
  const mxArray *c3_lhs100 = NULL;
  const mxArray *c3_rhs101 = NULL;
  const mxArray *c3_lhs101 = NULL;
  const mxArray *c3_rhs102 = NULL;
  const mxArray *c3_lhs102 = NULL;
  const mxArray *c3_rhs103 = NULL;
  const mxArray *c3_lhs103 = NULL;
  const mxArray *c3_rhs104 = NULL;
  const mxArray *c3_lhs104 = NULL;
  const mxArray *c3_rhs105 = NULL;
  const mxArray *c3_lhs105 = NULL;
  const mxArray *c3_rhs106 = NULL;
  const mxArray *c3_lhs106 = NULL;
  const mxArray *c3_rhs107 = NULL;
  const mxArray *c3_lhs107 = NULL;
  const mxArray *c3_rhs108 = NULL;
  const mxArray *c3_lhs108 = NULL;
  const mxArray *c3_rhs109 = NULL;
  const mxArray *c3_lhs109 = NULL;
  const mxArray *c3_rhs110 = NULL;
  const mxArray *c3_lhs110 = NULL;
  const mxArray *c3_rhs111 = NULL;
  const mxArray *c3_lhs111 = NULL;
  const mxArray *c3_rhs112 = NULL;
  const mxArray *c3_lhs112 = NULL;
  const mxArray *c3_rhs113 = NULL;
  const mxArray *c3_lhs113 = NULL;
  const mxArray *c3_rhs114 = NULL;
  const mxArray *c3_lhs114 = NULL;
  const mxArray *c3_rhs115 = NULL;
  const mxArray *c3_lhs115 = NULL;
  const mxArray *c3_rhs116 = NULL;
  const mxArray *c3_lhs116 = NULL;
  const mxArray *c3_rhs117 = NULL;
  const mxArray *c3_lhs117 = NULL;
  const mxArray *c3_rhs118 = NULL;
  const mxArray *c3_lhs118 = NULL;
  const mxArray *c3_rhs119 = NULL;
  const mxArray *c3_lhs119 = NULL;
  const mxArray *c3_rhs120 = NULL;
  const mxArray *c3_lhs120 = NULL;
  const mxArray *c3_rhs121 = NULL;
  const mxArray *c3_lhs121 = NULL;
  const mxArray *c3_rhs122 = NULL;
  const mxArray *c3_lhs122 = NULL;
  const mxArray *c3_rhs123 = NULL;
  const mxArray *c3_lhs123 = NULL;
  const mxArray *c3_rhs124 = NULL;
  const mxArray *c3_lhs124 = NULL;
  const mxArray *c3_rhs125 = NULL;
  const mxArray *c3_lhs125 = NULL;
  const mxArray *c3_rhs126 = NULL;
  const mxArray *c3_lhs126 = NULL;
  const mxArray *c3_rhs127 = NULL;
  const mxArray *c3_lhs127 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotc.m"), "context",
                  "context", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xdotc"),
                  "name", "name", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdotc.p"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c3_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdotc.p"),
                  "context", "context", 65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xdot"),
                  "name", "name", 65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "resolved", "resolved", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c3_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "context", "context", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c3_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p!below_threshold"),
                  "context", "context", 67);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 67);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 67);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c3_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 68);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 68);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 68);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381817900U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c3_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p!below_threshold"),
                  "context", "context", 69);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("length"), "name", "name", 69);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 69);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1303117406U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c3_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 70);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 70);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 70);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c3_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xdot.p"),
                  "context", "context", 71);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xdot"),
                  "name", "name", 71);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdot.p"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c3_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdot.p"),
                  "context", "context", 72);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xdotx"),
                  "name", "name", 72);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c3_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 73);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 73);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c3_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 74);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 74);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 74);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c3_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xdotx.p"),
                  "context", "context", 75);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 75);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 75);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c3_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 76);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 76);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 76);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 76);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c3_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 77);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 77);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 77);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 77);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c3_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 78);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgemv"), "name", "name",
                  78);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 78);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951890U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c3_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"), "context",
                  "context", 79);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 79);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 79);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c3_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"), "context",
                  "context", 80);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xgemv"),
                  "name", "name", 80);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c3_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "context", "context", 81);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 81);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 81);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c3_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!below_threshold"),
                  "context", "context", 82);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 82);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 82);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c3_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!below_threshold"),
                  "context", "context", 83);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("length"), "name", "name", 83);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 83);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1303117406U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c3_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!below_threshold"),
                  "context", "context", 84);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 84);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 84);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c3_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!below_threshold"),
                  "context", "context", 85);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 85);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 85);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c3_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "context", "context", 86);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 86);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 86);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c3_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "context", "context", 87);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xgemv"),
                  "name", "name", 87);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c3_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 88);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 88);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 88);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c3_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 89);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 89);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 89);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c3_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 90);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 90);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 90);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 90);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c3_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 91);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 91);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c3_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemv.p"),
                  "context", "context", 92);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 92);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 92);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c3_rhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p"),
                  "context", "context", 93);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 93);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 93);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c3_rhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!ceval_xgemv"),
                  "context", "context", 94);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.size_ptr"),
                  "name", "name", 94);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/size_ptr.p"),
                  "resolved", "resolved", 94);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c3_rhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemv.p!ceval_xgemv"),
                  "context", "context", 95);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.c_cast"),
                  "name", "name", 95);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("int32"), "dominantType",
                  "dominantType", 95);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/c_cast.p"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c3_rhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 96);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 96);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 96);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 96);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c3_rhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zpotrf.m"),
                  "context", "context", 97);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xscal"), "name", "name",
                  97);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 97);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"),
                  "resolved", "resolved", 97);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951892U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c3_rhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 98);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 98);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 98);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c3_rhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xscal.m"), "context",
                  "context", 99);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xscal"),
                  "name", "name", 99);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "resolved", "resolved", 99);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c3_rhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 100);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 100);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c3_rhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs100), "rhs", "rhs",
                  100);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs100), "lhs", "lhs",
                  100);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p!below_threshold"),
                  "context", "context", 101);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 101);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 101);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 101);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c3_rhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs101), "rhs", "rhs",
                  101);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs101), "lhs", "lhs",
                  101);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p!below_threshold"),
                  "context", "context", 102);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("length"), "name", "name", 102);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 102);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1303117406U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c3_rhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs102), "rhs", "rhs",
                  102);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs102), "lhs", "lhs",
                  102);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 103);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 103);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 103);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 103);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c3_rhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs103), "rhs", "rhs",
                  103);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs103), "lhs", "lhs",
                  103);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xscal.p"),
                  "context", "context", 104);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xscal"),
                  "name", "name", 104);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "resolved", "resolved", 104);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c3_rhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs104), "rhs", "rhs",
                  104);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs104), "lhs", "lhs",
                  104);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 105);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 105);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 105);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c3_rhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs105), "rhs", "rhs",
                  105);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs105), "lhs", "lhs",
                  105);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 106);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 106);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 106);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 106);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c3_rhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs106), "rhs", "rhs",
                  106);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs106), "lhs", "lhs",
                  106);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 107);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 107);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 107);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c3_rhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs107), "rhs", "rhs",
                  107);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs107), "lhs", "lhs",
                  107);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xscal.p"),
                  "context", "context", 108);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 108);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c3_rhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs108), "rhs", "rhs",
                  108);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs108), "lhs", "lhs",
                  108);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 109);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 109);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 109);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 109);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c3_rhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs109), "rhs", "rhs",
                  109);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs109), "lhs", "lhs",
                  109);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/chol.m!cholesky"),
                  "context", "context", 110);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 110);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 110);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c3_rhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs110), "rhs", "rhs",
                  110);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs110), "lhs", "lhs",
                  110);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 111);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 111);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 111);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 111);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c3_rhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs111), "rhs", "rhs",
                  111);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs111), "lhs", "lhs",
                  111);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 112);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("length"), "name", "name", 112);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 112);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1303117406U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c3_rhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs112), "rhs", "rhs",
                  112);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs112), "lhs", "lhs",
                  112);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[E]C:/Users/Rudaba/Documents/PhD Take 2/Simulations/AircraftModelling/simulateUKF.m"),
                  "context", "context", 113);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("inv"), "name", "name", 113);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 113);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1305289200U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c3_rhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs113), "rhs", "rhs",
                  113);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs113), "lhs", "lhs",
                  113);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 114);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 114);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 114);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c3_rhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs114), "rhs", "rhs",
                  114);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs114), "lhs", "lhs",
                  114);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 115);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgetrf"), "name", "name",
                  115);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786406U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c3_rhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs115), "rhs", "rhs",
                  115);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs115), "lhs", "lhs",
                  115);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "context", "context", 116);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_lapack_xgetrf"), "name",
                  "name", 116);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "resolved", "resolved", 116);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786410U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c3_rhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs116), "rhs", "rhs",
                  116);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs116), "lhs", "lhs",
                  116);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "context", "context", 117);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_matlab_zgetrf"), "name",
                  "name", 117);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "resolved", "resolved", 117);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1302660194U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c3_rhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs117), "rhs", "rhs",
                  117);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs117), "lhs", "lhs",
                  117);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 118);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("realmin"), "name", "name", 118);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 118);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 118);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1307622442U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c3_rhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs118), "rhs", "rhs",
                  118);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs118), "lhs", "lhs",
                  118);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 119);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_realmin"), "name", "name",
                  119);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 119);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1307622444U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c3_rhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs119), "rhs", "rhs",
                  119);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs119), "lhs", "lhs",
                  119);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 120);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 120);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 120);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326691996U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c3_rhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs120), "rhs", "rhs",
                  120);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs120), "lhs", "lhs",
                  120);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 121);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eps"), "name", "name", 121);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 121);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 121);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326691996U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c3_rhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs121), "rhs", "rhs",
                  121);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs121), "lhs", "lhs",
                  121);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 122);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 122);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 122);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786382U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c3_rhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs122), "rhs", "rhs",
                  122);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs122), "lhs", "lhs",
                  122);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 123);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_eps"), "name", "name", 123);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 123);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 123);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326691996U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c3_rhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs123), "rhs", "rhs",
                  123);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs123), "lhs", "lhs",
                  123);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 124);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 124);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 124);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 124);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326691996U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c3_rhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs124), "rhs", "rhs",
                  124);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs124), "lhs", "lhs",
                  124);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 125);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("min"), "name", "name", 125);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 125);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 125);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1311226518U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c3_rhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs125), "rhs", "rhs",
                  125);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs125), "lhs", "lhs",
                  125);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 126);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 126);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 126);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1378267184U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c3_rhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs126), "rhs", "rhs",
                  126);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs126), "lhs", "lhs",
                  126);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 127);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 127);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 127);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 127);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c3_rhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs127), "rhs", "rhs",
                  127);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs127), "lhs", "lhs",
                  127);
  sf_mex_destroy(&c3_rhs64);
  sf_mex_destroy(&c3_lhs64);
  sf_mex_destroy(&c3_rhs65);
  sf_mex_destroy(&c3_lhs65);
  sf_mex_destroy(&c3_rhs66);
  sf_mex_destroy(&c3_lhs66);
  sf_mex_destroy(&c3_rhs67);
  sf_mex_destroy(&c3_lhs67);
  sf_mex_destroy(&c3_rhs68);
  sf_mex_destroy(&c3_lhs68);
  sf_mex_destroy(&c3_rhs69);
  sf_mex_destroy(&c3_lhs69);
  sf_mex_destroy(&c3_rhs70);
  sf_mex_destroy(&c3_lhs70);
  sf_mex_destroy(&c3_rhs71);
  sf_mex_destroy(&c3_lhs71);
  sf_mex_destroy(&c3_rhs72);
  sf_mex_destroy(&c3_lhs72);
  sf_mex_destroy(&c3_rhs73);
  sf_mex_destroy(&c3_lhs73);
  sf_mex_destroy(&c3_rhs74);
  sf_mex_destroy(&c3_lhs74);
  sf_mex_destroy(&c3_rhs75);
  sf_mex_destroy(&c3_lhs75);
  sf_mex_destroy(&c3_rhs76);
  sf_mex_destroy(&c3_lhs76);
  sf_mex_destroy(&c3_rhs77);
  sf_mex_destroy(&c3_lhs77);
  sf_mex_destroy(&c3_rhs78);
  sf_mex_destroy(&c3_lhs78);
  sf_mex_destroy(&c3_rhs79);
  sf_mex_destroy(&c3_lhs79);
  sf_mex_destroy(&c3_rhs80);
  sf_mex_destroy(&c3_lhs80);
  sf_mex_destroy(&c3_rhs81);
  sf_mex_destroy(&c3_lhs81);
  sf_mex_destroy(&c3_rhs82);
  sf_mex_destroy(&c3_lhs82);
  sf_mex_destroy(&c3_rhs83);
  sf_mex_destroy(&c3_lhs83);
  sf_mex_destroy(&c3_rhs84);
  sf_mex_destroy(&c3_lhs84);
  sf_mex_destroy(&c3_rhs85);
  sf_mex_destroy(&c3_lhs85);
  sf_mex_destroy(&c3_rhs86);
  sf_mex_destroy(&c3_lhs86);
  sf_mex_destroy(&c3_rhs87);
  sf_mex_destroy(&c3_lhs87);
  sf_mex_destroy(&c3_rhs88);
  sf_mex_destroy(&c3_lhs88);
  sf_mex_destroy(&c3_rhs89);
  sf_mex_destroy(&c3_lhs89);
  sf_mex_destroy(&c3_rhs90);
  sf_mex_destroy(&c3_lhs90);
  sf_mex_destroy(&c3_rhs91);
  sf_mex_destroy(&c3_lhs91);
  sf_mex_destroy(&c3_rhs92);
  sf_mex_destroy(&c3_lhs92);
  sf_mex_destroy(&c3_rhs93);
  sf_mex_destroy(&c3_lhs93);
  sf_mex_destroy(&c3_rhs94);
  sf_mex_destroy(&c3_lhs94);
  sf_mex_destroy(&c3_rhs95);
  sf_mex_destroy(&c3_lhs95);
  sf_mex_destroy(&c3_rhs96);
  sf_mex_destroy(&c3_lhs96);
  sf_mex_destroy(&c3_rhs97);
  sf_mex_destroy(&c3_lhs97);
  sf_mex_destroy(&c3_rhs98);
  sf_mex_destroy(&c3_lhs98);
  sf_mex_destroy(&c3_rhs99);
  sf_mex_destroy(&c3_lhs99);
  sf_mex_destroy(&c3_rhs100);
  sf_mex_destroy(&c3_lhs100);
  sf_mex_destroy(&c3_rhs101);
  sf_mex_destroy(&c3_lhs101);
  sf_mex_destroy(&c3_rhs102);
  sf_mex_destroy(&c3_lhs102);
  sf_mex_destroy(&c3_rhs103);
  sf_mex_destroy(&c3_lhs103);
  sf_mex_destroy(&c3_rhs104);
  sf_mex_destroy(&c3_lhs104);
  sf_mex_destroy(&c3_rhs105);
  sf_mex_destroy(&c3_lhs105);
  sf_mex_destroy(&c3_rhs106);
  sf_mex_destroy(&c3_lhs106);
  sf_mex_destroy(&c3_rhs107);
  sf_mex_destroy(&c3_lhs107);
  sf_mex_destroy(&c3_rhs108);
  sf_mex_destroy(&c3_lhs108);
  sf_mex_destroy(&c3_rhs109);
  sf_mex_destroy(&c3_lhs109);
  sf_mex_destroy(&c3_rhs110);
  sf_mex_destroy(&c3_lhs110);
  sf_mex_destroy(&c3_rhs111);
  sf_mex_destroy(&c3_lhs111);
  sf_mex_destroy(&c3_rhs112);
  sf_mex_destroy(&c3_lhs112);
  sf_mex_destroy(&c3_rhs113);
  sf_mex_destroy(&c3_lhs113);
  sf_mex_destroy(&c3_rhs114);
  sf_mex_destroy(&c3_lhs114);
  sf_mex_destroy(&c3_rhs115);
  sf_mex_destroy(&c3_lhs115);
  sf_mex_destroy(&c3_rhs116);
  sf_mex_destroy(&c3_lhs116);
  sf_mex_destroy(&c3_rhs117);
  sf_mex_destroy(&c3_lhs117);
  sf_mex_destroy(&c3_rhs118);
  sf_mex_destroy(&c3_lhs118);
  sf_mex_destroy(&c3_rhs119);
  sf_mex_destroy(&c3_lhs119);
  sf_mex_destroy(&c3_rhs120);
  sf_mex_destroy(&c3_lhs120);
  sf_mex_destroy(&c3_rhs121);
  sf_mex_destroy(&c3_lhs121);
  sf_mex_destroy(&c3_rhs122);
  sf_mex_destroy(&c3_lhs122);
  sf_mex_destroy(&c3_rhs123);
  sf_mex_destroy(&c3_lhs123);
  sf_mex_destroy(&c3_rhs124);
  sf_mex_destroy(&c3_lhs124);
  sf_mex_destroy(&c3_rhs125);
  sf_mex_destroy(&c3_lhs125);
  sf_mex_destroy(&c3_rhs126);
  sf_mex_destroy(&c3_lhs126);
  sf_mex_destroy(&c3_rhs127);
  sf_mex_destroy(&c3_lhs127);
}

static void c3_c_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs128 = NULL;
  const mxArray *c3_lhs128 = NULL;
  const mxArray *c3_rhs129 = NULL;
  const mxArray *c3_lhs129 = NULL;
  const mxArray *c3_rhs130 = NULL;
  const mxArray *c3_lhs130 = NULL;
  const mxArray *c3_rhs131 = NULL;
  const mxArray *c3_lhs131 = NULL;
  const mxArray *c3_rhs132 = NULL;
  const mxArray *c3_lhs132 = NULL;
  const mxArray *c3_rhs133 = NULL;
  const mxArray *c3_lhs133 = NULL;
  const mxArray *c3_rhs134 = NULL;
  const mxArray *c3_lhs134 = NULL;
  const mxArray *c3_rhs135 = NULL;
  const mxArray *c3_lhs135 = NULL;
  const mxArray *c3_rhs136 = NULL;
  const mxArray *c3_lhs136 = NULL;
  const mxArray *c3_rhs137 = NULL;
  const mxArray *c3_lhs137 = NULL;
  const mxArray *c3_rhs138 = NULL;
  const mxArray *c3_lhs138 = NULL;
  const mxArray *c3_rhs139 = NULL;
  const mxArray *c3_lhs139 = NULL;
  const mxArray *c3_rhs140 = NULL;
  const mxArray *c3_lhs140 = NULL;
  const mxArray *c3_rhs141 = NULL;
  const mxArray *c3_lhs141 = NULL;
  const mxArray *c3_rhs142 = NULL;
  const mxArray *c3_lhs142 = NULL;
  const mxArray *c3_rhs143 = NULL;
  const mxArray *c3_lhs143 = NULL;
  const mxArray *c3_rhs144 = NULL;
  const mxArray *c3_lhs144 = NULL;
  const mxArray *c3_rhs145 = NULL;
  const mxArray *c3_lhs145 = NULL;
  const mxArray *c3_rhs146 = NULL;
  const mxArray *c3_lhs146 = NULL;
  const mxArray *c3_rhs147 = NULL;
  const mxArray *c3_lhs147 = NULL;
  const mxArray *c3_rhs148 = NULL;
  const mxArray *c3_lhs148 = NULL;
  const mxArray *c3_rhs149 = NULL;
  const mxArray *c3_lhs149 = NULL;
  const mxArray *c3_rhs150 = NULL;
  const mxArray *c3_lhs150 = NULL;
  const mxArray *c3_rhs151 = NULL;
  const mxArray *c3_lhs151 = NULL;
  const mxArray *c3_rhs152 = NULL;
  const mxArray *c3_lhs152 = NULL;
  const mxArray *c3_rhs153 = NULL;
  const mxArray *c3_lhs153 = NULL;
  const mxArray *c3_rhs154 = NULL;
  const mxArray *c3_lhs154 = NULL;
  const mxArray *c3_rhs155 = NULL;
  const mxArray *c3_lhs155 = NULL;
  const mxArray *c3_rhs156 = NULL;
  const mxArray *c3_lhs156 = NULL;
  const mxArray *c3_rhs157 = NULL;
  const mxArray *c3_lhs157 = NULL;
  const mxArray *c3_rhs158 = NULL;
  const mxArray *c3_lhs158 = NULL;
  const mxArray *c3_rhs159 = NULL;
  const mxArray *c3_lhs159 = NULL;
  const mxArray *c3_rhs160 = NULL;
  const mxArray *c3_lhs160 = NULL;
  const mxArray *c3_rhs161 = NULL;
  const mxArray *c3_lhs161 = NULL;
  const mxArray *c3_rhs162 = NULL;
  const mxArray *c3_lhs162 = NULL;
  const mxArray *c3_rhs163 = NULL;
  const mxArray *c3_lhs163 = NULL;
  const mxArray *c3_rhs164 = NULL;
  const mxArray *c3_lhs164 = NULL;
  const mxArray *c3_rhs165 = NULL;
  const mxArray *c3_lhs165 = NULL;
  const mxArray *c3_rhs166 = NULL;
  const mxArray *c3_lhs166 = NULL;
  const mxArray *c3_rhs167 = NULL;
  const mxArray *c3_lhs167 = NULL;
  const mxArray *c3_rhs168 = NULL;
  const mxArray *c3_lhs168 = NULL;
  const mxArray *c3_rhs169 = NULL;
  const mxArray *c3_lhs169 = NULL;
  const mxArray *c3_rhs170 = NULL;
  const mxArray *c3_lhs170 = NULL;
  const mxArray *c3_rhs171 = NULL;
  const mxArray *c3_lhs171 = NULL;
  const mxArray *c3_rhs172 = NULL;
  const mxArray *c3_lhs172 = NULL;
  const mxArray *c3_rhs173 = NULL;
  const mxArray *c3_lhs173 = NULL;
  const mxArray *c3_rhs174 = NULL;
  const mxArray *c3_lhs174 = NULL;
  const mxArray *c3_rhs175 = NULL;
  const mxArray *c3_lhs175 = NULL;
  const mxArray *c3_rhs176 = NULL;
  const mxArray *c3_lhs176 = NULL;
  const mxArray *c3_rhs177 = NULL;
  const mxArray *c3_lhs177 = NULL;
  const mxArray *c3_rhs178 = NULL;
  const mxArray *c3_lhs178 = NULL;
  const mxArray *c3_rhs179 = NULL;
  const mxArray *c3_lhs179 = NULL;
  const mxArray *c3_rhs180 = NULL;
  const mxArray *c3_lhs180 = NULL;
  const mxArray *c3_rhs181 = NULL;
  const mxArray *c3_lhs181 = NULL;
  const mxArray *c3_rhs182 = NULL;
  const mxArray *c3_lhs182 = NULL;
  const mxArray *c3_rhs183 = NULL;
  const mxArray *c3_lhs183 = NULL;
  const mxArray *c3_rhs184 = NULL;
  const mxArray *c3_lhs184 = NULL;
  const mxArray *c3_rhs185 = NULL;
  const mxArray *c3_lhs185 = NULL;
  const mxArray *c3_rhs186 = NULL;
  const mxArray *c3_lhs186 = NULL;
  const mxArray *c3_rhs187 = NULL;
  const mxArray *c3_lhs187 = NULL;
  const mxArray *c3_rhs188 = NULL;
  const mxArray *c3_lhs188 = NULL;
  const mxArray *c3_rhs189 = NULL;
  const mxArray *c3_lhs189 = NULL;
  const mxArray *c3_rhs190 = NULL;
  const mxArray *c3_lhs190 = NULL;
  const mxArray *c3_rhs191 = NULL;
  const mxArray *c3_lhs191 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 128);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 128);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 128);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 128);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c3_rhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs128), "rhs", "rhs",
                  128);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs128), "lhs", "lhs",
                  128);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 129);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 129);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 129);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 129);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c3_rhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs129), "rhs", "rhs",
                  129);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs129), "lhs", "lhs",
                  129);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 130);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 130);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 130);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 130);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c3_rhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs130), "rhs", "rhs",
                  130);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs130), "lhs", "lhs",
                  130);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 131);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 131);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 131);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c3_rhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs131), "rhs", "rhs",
                  131);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs131), "lhs", "lhs",
                  131);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 132);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 132);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 132);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 132);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c3_rhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs132), "rhs", "rhs",
                  132);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs132), "lhs", "lhs",
                  132);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 133);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 133);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 133);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 133);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c3_rhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs133), "rhs", "rhs",
                  133);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs133), "lhs", "lhs",
                  133);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 134);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("colon"), "name", "name", 134);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 134);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1378267188U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c3_rhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs134), "rhs", "rhs",
                  134);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs134), "lhs", "lhs",
                  134);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 135);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("colon"), "name", "name", 135);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 135);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1378267188U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c3_rhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs135), "rhs", "rhs",
                  135);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs135), "lhs", "lhs",
                  135);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 136);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 136);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 136);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c3_rhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs136), "rhs", "rhs",
                  136);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs136), "lhs", "lhs",
                  136);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 137);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 137);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 137);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 137);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c3_rhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs137), "rhs", "rhs",
                  137);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs137), "lhs", "lhs",
                  137);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 138);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("floor"), "name", "name", 138);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 138);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677854U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c3_rhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs138), "rhs", "rhs",
                  138);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs138), "lhs", "lhs",
                  138);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 139);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmin"), "name", "name", 139);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 139);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c3_rhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs139), "rhs", "rhs",
                  139);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs139), "lhs", "lhs",
                  139);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 140);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 140);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 140);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c3_rhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs140), "rhs", "rhs",
                  140);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs140), "lhs", "lhs",
                  140);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 141);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmin"), "name", "name", 141);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 141);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 141);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c3_rhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs141), "rhs", "rhs",
                  141);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs141), "lhs", "lhs",
                  141);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 142);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 142);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 142);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 142);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c3_rhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs142), "rhs", "rhs",
                  142);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs142), "lhs", "lhs",
                  142);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 143);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_isa_uint"), "name", "name",
                  143);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 143);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 143);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c3_rhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs143), "rhs", "rhs",
                  143);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs143), "lhs", "lhs",
                  143);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "context",
                  "context", 144);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.isaUint"),
                  "name", "name", 144);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 144);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/isaUint.p"),
                  "resolved", "resolved", 144);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c3_rhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs144), "rhs", "rhs",
                  144);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs144), "lhs", "lhs",
                  144);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 145);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_unsigned_class"), "name",
                  "name", 145);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "resolved", "resolved", 145);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c3_rhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs145), "rhs", "rhs",
                  145);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs145), "lhs", "lhs",
                  145);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "context", "context", 146);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.unsignedClass"),
                  "name", "name", 146);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "resolved", "resolved", 146);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c3_rhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs146), "rhs", "rhs",
                  146);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs146), "lhs", "lhs",
                  146);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "context", "context", 147);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 147);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 147);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1381817900U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c3_rhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs147), "rhs", "rhs",
                  147);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs147), "lhs", "lhs",
                  147);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "context", "context", 148);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 148);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 148);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c3_rhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs148), "rhs", "rhs",
                  148);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs148), "lhs", "lhs",
                  148);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 149);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 149);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 149);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 149);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c3_rhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs149), "rhs", "rhs",
                  149);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs149), "lhs", "lhs",
                  149);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 150);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 150);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 150);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c3_rhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs150), "rhs", "rhs",
                  150);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs150), "lhs", "lhs",
                  150);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 151);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_isa_uint"), "name", "name",
                  151);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 151);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 151);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c3_rhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs151), "rhs", "rhs",
                  151);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs151), "lhs", "lhs",
                  151);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 152);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 152);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 152);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c3_rhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs152), "rhs", "rhs",
                  152);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs152), "lhs", "lhs",
                  152);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon"),
                  "context", "context", 153);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 153);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 153);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c3_rhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs153), "rhs", "rhs",
                  153);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs153), "lhs", "lhs",
                  153);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 154);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 154);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 154);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c3_rhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs154), "rhs", "rhs",
                  154);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs154), "lhs", "lhs",
                  154);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 155);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 155);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 155);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c3_rhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs155), "rhs", "rhs",
                  155);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs155), "lhs", "lhs",
                  155);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 156);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 156);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 156);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c3_rhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs156), "rhs", "rhs",
                  156);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs156), "lhs", "lhs",
                  156);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 157);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 157);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 157);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 157);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c3_rhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs157), "rhs", "rhs",
                  157);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs157), "lhs", "lhs",
                  157);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 158);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 158);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 158);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c3_rhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs158), "rhs", "rhs",
                  158);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs158), "lhs", "lhs",
                  158);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 159);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 159);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 159);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 159);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c3_rhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs159), "rhs", "rhs",
                  159);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs159), "lhs", "lhs",
                  159);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 160);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 160);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 160);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 160);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c3_rhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs160), "rhs", "rhs",
                  160);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs160), "lhs", "lhs",
                  160);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 161);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 161);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 161);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 161);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c3_rhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs161), "rhs", "rhs",
                  161);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs161), "lhs", "lhs",
                  161);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 162);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_ixamax"), "name", "name",
                  162);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 162);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "resolved", "resolved", 162);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c3_rhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs162), "rhs", "rhs",
                  162);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs162), "lhs", "lhs",
                  162);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 163);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 163);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c3_rhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs163), "rhs", "rhs",
                  163);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs163), "lhs", "lhs",
                  163);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 164);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.ixamax"),
                  "name", "name", 164);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 164);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "resolved", "resolved", 164);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c3_rhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs164), "rhs", "rhs",
                  164);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs164), "lhs", "lhs",
                  164);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "context", "context", 165);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 165);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 165);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 165);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c3_rhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs165), "rhs", "rhs",
                  165);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs165), "lhs", "lhs",
                  165);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p!below_threshold"),
                  "context", "context", 166);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 166);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 166);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 166);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c3_rhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs166), "rhs", "rhs",
                  166);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs166), "lhs", "lhs",
                  166);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p!below_threshold"),
                  "context", "context", 167);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("length"), "name", "name", 167);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 167);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1303117406U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c3_rhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs167), "rhs", "rhs",
                  167);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs167), "lhs", "lhs",
                  167);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "context", "context", 168);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.ixamax"),
                  "name", "name", 168);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "resolved", "resolved", 168);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c3_rhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs168), "rhs", "rhs",
                  168);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs168), "lhs", "lhs",
                  168);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 169);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xcabs1"),
                  "name", "name", 169);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "resolved", "resolved", 169);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c3_rhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs169), "rhs", "rhs",
                  169);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs169), "lhs", "lhs",
                  169);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "context", "context", 170);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("abs"), "name", "name", 170);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 170);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677852U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c3_rhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs170), "rhs", "rhs",
                  170);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs170), "lhs", "lhs",
                  170);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 171);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 171);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 171);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c3_rhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs171), "rhs", "rhs",
                  171);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs171), "lhs", "lhs",
                  171);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 172);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 172);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 172);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786312U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c3_rhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs172), "rhs", "rhs",
                  172);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs172), "lhs", "lhs",
                  172);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 173);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 173);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 173);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 173);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c3_rhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs173), "rhs", "rhs",
                  173);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs173), "lhs", "lhs",
                  173);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 174);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 174);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 174);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 174);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c3_rhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs174), "rhs", "rhs",
                  174);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs174), "lhs", "lhs",
                  174);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 175);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xswap"), "name", "name",
                  175);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 175);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"),
                  "resolved", "resolved", 175);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951892U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c3_rhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs175), "rhs", "rhs",
                  175);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs175), "lhs", "lhs",
                  175);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 176);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 176);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 176);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c3_rhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs176), "rhs", "rhs",
                  176);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs176), "lhs", "lhs",
                  176);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 177);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xswap"),
                  "name", "name", 177);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "resolved", "resolved", 177);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c3_rhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs177), "rhs", "rhs",
                  177);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs177), "lhs", "lhs",
                  177);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 178);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 178);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 178);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 178);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c3_rhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs178), "rhs", "rhs",
                  178);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs178), "lhs", "lhs",
                  178);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p!below_threshold"),
                  "context", "context", 179);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 179);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 179);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 179);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c3_rhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs179), "rhs", "rhs",
                  179);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs179), "lhs", "lhs",
                  179);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 180);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xswap"),
                  "name", "name", 180);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 180);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "resolved", "resolved", 180);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c3_rhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs180), "rhs", "rhs",
                  180);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs180), "lhs", "lhs",
                  180);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 181);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("abs"), "name", "name", 181);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 181);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 181);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677852U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c3_rhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs181), "rhs", "rhs",
                  181);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs181), "lhs", "lhs",
                  181);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 182);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 182);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 182);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 182);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c3_rhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs182), "rhs", "rhs",
                  182);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs182), "lhs", "lhs",
                  182);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 183);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 183);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 183);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 183);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786312U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c3_rhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs183), "rhs", "rhs",
                  183);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs183), "lhs", "lhs",
                  183);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 184);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 184);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 184);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 184);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 184);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 184);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 184);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 184);
  sf_mex_assign(&c3_rhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs184), "rhs", "rhs",
                  184);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs184), "lhs", "lhs",
                  184);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 185);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 185);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 185);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 185);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 185);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 185);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 185);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 185);
  sf_mex_assign(&c3_rhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs185), "rhs", "rhs",
                  185);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs185), "lhs", "lhs",
                  185);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 186);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_div"), "name", "name", 186);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 186);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 186);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 186);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 186);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 186);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 186);
  sf_mex_assign(&c3_rhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs186), "rhs", "rhs",
                  186);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs186), "lhs", "lhs",
                  186);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 187);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgeru"), "name", "name",
                  187);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 187);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"),
                  "resolved", "resolved", 187);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951890U), "fileTimeLo",
                  "fileTimeLo", 187);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 187);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 187);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 187);
  sf_mex_assign(&c3_rhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs187), "rhs", "rhs",
                  187);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs187), "lhs", "lhs",
                  187);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 188);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 188);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 188);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 188);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 188);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 188);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 188);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 188);
  sf_mex_assign(&c3_rhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs188), "rhs", "rhs",
                  188);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs188), "lhs", "lhs",
                  188);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 189);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xgeru"),
                  "name", "name", 189);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 189);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgeru.p"),
                  "resolved", "resolved", 189);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 189);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 189);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 189);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 189);
  sf_mex_assign(&c3_rhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs189), "rhs", "rhs",
                  189);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs189), "lhs", "lhs",
                  189);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgeru.p"),
                  "context", "context", 190);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xger"),
                  "name", "name", 190);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 190);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "resolved", "resolved", 190);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 190);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 190);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 190);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 190);
  sf_mex_assign(&c3_rhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs190), "rhs", "rhs",
                  190);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs190), "lhs", "lhs",
                  190);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "context", "context", 191);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 191);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 191);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 191);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 191);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 191);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 191);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 191);
  sf_mex_assign(&c3_rhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs191), "rhs", "rhs",
                  191);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs191), "lhs", "lhs",
                  191);
  sf_mex_destroy(&c3_rhs128);
  sf_mex_destroy(&c3_lhs128);
  sf_mex_destroy(&c3_rhs129);
  sf_mex_destroy(&c3_lhs129);
  sf_mex_destroy(&c3_rhs130);
  sf_mex_destroy(&c3_lhs130);
  sf_mex_destroy(&c3_rhs131);
  sf_mex_destroy(&c3_lhs131);
  sf_mex_destroy(&c3_rhs132);
  sf_mex_destroy(&c3_lhs132);
  sf_mex_destroy(&c3_rhs133);
  sf_mex_destroy(&c3_lhs133);
  sf_mex_destroy(&c3_rhs134);
  sf_mex_destroy(&c3_lhs134);
  sf_mex_destroy(&c3_rhs135);
  sf_mex_destroy(&c3_lhs135);
  sf_mex_destroy(&c3_rhs136);
  sf_mex_destroy(&c3_lhs136);
  sf_mex_destroy(&c3_rhs137);
  sf_mex_destroy(&c3_lhs137);
  sf_mex_destroy(&c3_rhs138);
  sf_mex_destroy(&c3_lhs138);
  sf_mex_destroy(&c3_rhs139);
  sf_mex_destroy(&c3_lhs139);
  sf_mex_destroy(&c3_rhs140);
  sf_mex_destroy(&c3_lhs140);
  sf_mex_destroy(&c3_rhs141);
  sf_mex_destroy(&c3_lhs141);
  sf_mex_destroy(&c3_rhs142);
  sf_mex_destroy(&c3_lhs142);
  sf_mex_destroy(&c3_rhs143);
  sf_mex_destroy(&c3_lhs143);
  sf_mex_destroy(&c3_rhs144);
  sf_mex_destroy(&c3_lhs144);
  sf_mex_destroy(&c3_rhs145);
  sf_mex_destroy(&c3_lhs145);
  sf_mex_destroy(&c3_rhs146);
  sf_mex_destroy(&c3_lhs146);
  sf_mex_destroy(&c3_rhs147);
  sf_mex_destroy(&c3_lhs147);
  sf_mex_destroy(&c3_rhs148);
  sf_mex_destroy(&c3_lhs148);
  sf_mex_destroy(&c3_rhs149);
  sf_mex_destroy(&c3_lhs149);
  sf_mex_destroy(&c3_rhs150);
  sf_mex_destroy(&c3_lhs150);
  sf_mex_destroy(&c3_rhs151);
  sf_mex_destroy(&c3_lhs151);
  sf_mex_destroy(&c3_rhs152);
  sf_mex_destroy(&c3_lhs152);
  sf_mex_destroy(&c3_rhs153);
  sf_mex_destroy(&c3_lhs153);
  sf_mex_destroy(&c3_rhs154);
  sf_mex_destroy(&c3_lhs154);
  sf_mex_destroy(&c3_rhs155);
  sf_mex_destroy(&c3_lhs155);
  sf_mex_destroy(&c3_rhs156);
  sf_mex_destroy(&c3_lhs156);
  sf_mex_destroy(&c3_rhs157);
  sf_mex_destroy(&c3_lhs157);
  sf_mex_destroy(&c3_rhs158);
  sf_mex_destroy(&c3_lhs158);
  sf_mex_destroy(&c3_rhs159);
  sf_mex_destroy(&c3_lhs159);
  sf_mex_destroy(&c3_rhs160);
  sf_mex_destroy(&c3_lhs160);
  sf_mex_destroy(&c3_rhs161);
  sf_mex_destroy(&c3_lhs161);
  sf_mex_destroy(&c3_rhs162);
  sf_mex_destroy(&c3_lhs162);
  sf_mex_destroy(&c3_rhs163);
  sf_mex_destroy(&c3_lhs163);
  sf_mex_destroy(&c3_rhs164);
  sf_mex_destroy(&c3_lhs164);
  sf_mex_destroy(&c3_rhs165);
  sf_mex_destroy(&c3_lhs165);
  sf_mex_destroy(&c3_rhs166);
  sf_mex_destroy(&c3_lhs166);
  sf_mex_destroy(&c3_rhs167);
  sf_mex_destroy(&c3_lhs167);
  sf_mex_destroy(&c3_rhs168);
  sf_mex_destroy(&c3_lhs168);
  sf_mex_destroy(&c3_rhs169);
  sf_mex_destroy(&c3_lhs169);
  sf_mex_destroy(&c3_rhs170);
  sf_mex_destroy(&c3_lhs170);
  sf_mex_destroy(&c3_rhs171);
  sf_mex_destroy(&c3_lhs171);
  sf_mex_destroy(&c3_rhs172);
  sf_mex_destroy(&c3_lhs172);
  sf_mex_destroy(&c3_rhs173);
  sf_mex_destroy(&c3_lhs173);
  sf_mex_destroy(&c3_rhs174);
  sf_mex_destroy(&c3_lhs174);
  sf_mex_destroy(&c3_rhs175);
  sf_mex_destroy(&c3_lhs175);
  sf_mex_destroy(&c3_rhs176);
  sf_mex_destroy(&c3_lhs176);
  sf_mex_destroy(&c3_rhs177);
  sf_mex_destroy(&c3_lhs177);
  sf_mex_destroy(&c3_rhs178);
  sf_mex_destroy(&c3_lhs178);
  sf_mex_destroy(&c3_rhs179);
  sf_mex_destroy(&c3_lhs179);
  sf_mex_destroy(&c3_rhs180);
  sf_mex_destroy(&c3_lhs180);
  sf_mex_destroy(&c3_rhs181);
  sf_mex_destroy(&c3_lhs181);
  sf_mex_destroy(&c3_rhs182);
  sf_mex_destroy(&c3_lhs182);
  sf_mex_destroy(&c3_rhs183);
  sf_mex_destroy(&c3_lhs183);
  sf_mex_destroy(&c3_rhs184);
  sf_mex_destroy(&c3_lhs184);
  sf_mex_destroy(&c3_rhs185);
  sf_mex_destroy(&c3_lhs185);
  sf_mex_destroy(&c3_rhs186);
  sf_mex_destroy(&c3_lhs186);
  sf_mex_destroy(&c3_rhs187);
  sf_mex_destroy(&c3_lhs187);
  sf_mex_destroy(&c3_rhs188);
  sf_mex_destroy(&c3_lhs188);
  sf_mex_destroy(&c3_rhs189);
  sf_mex_destroy(&c3_lhs189);
  sf_mex_destroy(&c3_rhs190);
  sf_mex_destroy(&c3_lhs190);
  sf_mex_destroy(&c3_rhs191);
  sf_mex_destroy(&c3_lhs191);
}

static void c3_d_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs192 = NULL;
  const mxArray *c3_lhs192 = NULL;
  const mxArray *c3_rhs193 = NULL;
  const mxArray *c3_lhs193 = NULL;
  const mxArray *c3_rhs194 = NULL;
  const mxArray *c3_lhs194 = NULL;
  const mxArray *c3_rhs195 = NULL;
  const mxArray *c3_lhs195 = NULL;
  const mxArray *c3_rhs196 = NULL;
  const mxArray *c3_lhs196 = NULL;
  const mxArray *c3_rhs197 = NULL;
  const mxArray *c3_lhs197 = NULL;
  const mxArray *c3_rhs198 = NULL;
  const mxArray *c3_lhs198 = NULL;
  const mxArray *c3_rhs199 = NULL;
  const mxArray *c3_lhs199 = NULL;
  const mxArray *c3_rhs200 = NULL;
  const mxArray *c3_lhs200 = NULL;
  const mxArray *c3_rhs201 = NULL;
  const mxArray *c3_lhs201 = NULL;
  const mxArray *c3_rhs202 = NULL;
  const mxArray *c3_lhs202 = NULL;
  const mxArray *c3_rhs203 = NULL;
  const mxArray *c3_lhs203 = NULL;
  const mxArray *c3_rhs204 = NULL;
  const mxArray *c3_lhs204 = NULL;
  const mxArray *c3_rhs205 = NULL;
  const mxArray *c3_lhs205 = NULL;
  const mxArray *c3_rhs206 = NULL;
  const mxArray *c3_lhs206 = NULL;
  const mxArray *c3_rhs207 = NULL;
  const mxArray *c3_lhs207 = NULL;
  const mxArray *c3_rhs208 = NULL;
  const mxArray *c3_lhs208 = NULL;
  const mxArray *c3_rhs209 = NULL;
  const mxArray *c3_lhs209 = NULL;
  const mxArray *c3_rhs210 = NULL;
  const mxArray *c3_lhs210 = NULL;
  const mxArray *c3_rhs211 = NULL;
  const mxArray *c3_lhs211 = NULL;
  const mxArray *c3_rhs212 = NULL;
  const mxArray *c3_lhs212 = NULL;
  const mxArray *c3_rhs213 = NULL;
  const mxArray *c3_lhs213 = NULL;
  const mxArray *c3_rhs214 = NULL;
  const mxArray *c3_lhs214 = NULL;
  const mxArray *c3_rhs215 = NULL;
  const mxArray *c3_lhs215 = NULL;
  const mxArray *c3_rhs216 = NULL;
  const mxArray *c3_lhs216 = NULL;
  const mxArray *c3_rhs217 = NULL;
  const mxArray *c3_lhs217 = NULL;
  const mxArray *c3_rhs218 = NULL;
  const mxArray *c3_lhs218 = NULL;
  const mxArray *c3_rhs219 = NULL;
  const mxArray *c3_lhs219 = NULL;
  const mxArray *c3_rhs220 = NULL;
  const mxArray *c3_lhs220 = NULL;
  const mxArray *c3_rhs221 = NULL;
  const mxArray *c3_lhs221 = NULL;
  const mxArray *c3_rhs222 = NULL;
  const mxArray *c3_lhs222 = NULL;
  const mxArray *c3_rhs223 = NULL;
  const mxArray *c3_lhs223 = NULL;
  const mxArray *c3_rhs224 = NULL;
  const mxArray *c3_lhs224 = NULL;
  const mxArray *c3_rhs225 = NULL;
  const mxArray *c3_lhs225 = NULL;
  const mxArray *c3_rhs226 = NULL;
  const mxArray *c3_lhs226 = NULL;
  const mxArray *c3_rhs227 = NULL;
  const mxArray *c3_lhs227 = NULL;
  const mxArray *c3_rhs228 = NULL;
  const mxArray *c3_lhs228 = NULL;
  const mxArray *c3_rhs229 = NULL;
  const mxArray *c3_lhs229 = NULL;
  const mxArray *c3_rhs230 = NULL;
  const mxArray *c3_lhs230 = NULL;
  const mxArray *c3_rhs231 = NULL;
  const mxArray *c3_lhs231 = NULL;
  const mxArray *c3_rhs232 = NULL;
  const mxArray *c3_lhs232 = NULL;
  const mxArray *c3_rhs233 = NULL;
  const mxArray *c3_lhs233 = NULL;
  const mxArray *c3_rhs234 = NULL;
  const mxArray *c3_lhs234 = NULL;
  const mxArray *c3_rhs235 = NULL;
  const mxArray *c3_lhs235 = NULL;
  const mxArray *c3_rhs236 = NULL;
  const mxArray *c3_lhs236 = NULL;
  const mxArray *c3_rhs237 = NULL;
  const mxArray *c3_lhs237 = NULL;
  const mxArray *c3_rhs238 = NULL;
  const mxArray *c3_lhs238 = NULL;
  const mxArray *c3_rhs239 = NULL;
  const mxArray *c3_lhs239 = NULL;
  const mxArray *c3_rhs240 = NULL;
  const mxArray *c3_lhs240 = NULL;
  const mxArray *c3_rhs241 = NULL;
  const mxArray *c3_lhs241 = NULL;
  const mxArray *c3_rhs242 = NULL;
  const mxArray *c3_lhs242 = NULL;
  const mxArray *c3_rhs243 = NULL;
  const mxArray *c3_lhs243 = NULL;
  const mxArray *c3_rhs244 = NULL;
  const mxArray *c3_lhs244 = NULL;
  const mxArray *c3_rhs245 = NULL;
  const mxArray *c3_lhs245 = NULL;
  const mxArray *c3_rhs246 = NULL;
  const mxArray *c3_lhs246 = NULL;
  const mxArray *c3_rhs247 = NULL;
  const mxArray *c3_lhs247 = NULL;
  const mxArray *c3_rhs248 = NULL;
  const mxArray *c3_lhs248 = NULL;
  const mxArray *c3_rhs249 = NULL;
  const mxArray *c3_lhs249 = NULL;
  const mxArray *c3_rhs250 = NULL;
  const mxArray *c3_lhs250 = NULL;
  const mxArray *c3_rhs251 = NULL;
  const mxArray *c3_lhs251 = NULL;
  const mxArray *c3_rhs252 = NULL;
  const mxArray *c3_lhs252 = NULL;
  const mxArray *c3_rhs253 = NULL;
  const mxArray *c3_lhs253 = NULL;
  const mxArray *c3_rhs254 = NULL;
  const mxArray *c3_lhs254 = NULL;
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 192);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 192);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 192);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 192);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 192);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 192);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 192);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 192);
  sf_mex_assign(&c3_rhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs192), "rhs", "rhs",
                  192);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs192), "lhs", "lhs",
                  192);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 193);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 193);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 193);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 193);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 193);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 193);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 193);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 193);
  sf_mex_assign(&c3_rhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs193), "rhs", "rhs",
                  193);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs193), "lhs", "lhs",
                  193);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 194);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmax"), "name", "name", 194);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 194);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 194);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 194);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 194);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 194);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 194);
  sf_mex_assign(&c3_rhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs194), "rhs", "rhs",
                  194);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs194), "lhs", "lhs",
                  194);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 195);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("min"), "name", "name", 195);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 195);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 195);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1311226518U), "fileTimeLo",
                  "fileTimeLo", 195);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 195);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 195);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 195);
  sf_mex_assign(&c3_rhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs195), "rhs", "rhs",
                  195);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs195), "lhs", "lhs",
                  195);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 196);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 196);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 196);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 196);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 196);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 196);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 196);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 196);
  sf_mex_assign(&c3_rhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs196), "rhs", "rhs",
                  196);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs196), "lhs", "lhs",
                  196);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 197);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 197);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 197);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 197);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 197);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 197);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 197);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 197);
  sf_mex_assign(&c3_rhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs197), "rhs", "rhs",
                  197);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs197), "lhs", "lhs",
                  197);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 198);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 198);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 198);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 198);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 198);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 198);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 198);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 198);
  sf_mex_assign(&c3_rhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs198), "rhs", "rhs",
                  198);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs198), "lhs", "lhs",
                  198);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 199);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 199);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 199);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 199);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 199);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 199);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 199);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 199);
  sf_mex_assign(&c3_rhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs199), "rhs", "rhs",
                  199);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs199), "lhs", "lhs",
                  199);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "context", "context", 200);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xger"),
                  "name", "name", 200);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 200);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xger.p"),
                  "resolved", "resolved", 200);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 200);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 200);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 200);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 200);
  sf_mex_assign(&c3_rhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs200), "rhs", "rhs",
                  200);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs200), "lhs", "lhs",
                  200);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xger.p"),
                  "context", "context", 201);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xgerx"),
                  "name", "name", 201);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 201);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "resolved", "resolved", 201);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 201);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 201);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 201);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 201);
  sf_mex_assign(&c3_rhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs201), "rhs", "rhs",
                  201);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs201), "lhs", "lhs",
                  201);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 202);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("abs"), "name", "name", 202);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 202);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 202);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677852U), "fileTimeLo",
                  "fileTimeLo", 202);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 202);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 202);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 202);
  sf_mex_assign(&c3_rhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs202), "rhs", "rhs",
                  202);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs202), "lhs", "lhs",
                  202);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 203);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 203);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 203);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 203);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 203);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 203);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 203);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 203);
  sf_mex_assign(&c3_rhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs203), "rhs", "rhs",
                  203);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs203), "lhs", "lhs",
                  203);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 204);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 204);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 204);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 204);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 204);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 204);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 204);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 204);
  sf_mex_assign(&c3_rhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs204), "rhs", "rhs",
                  204);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs204), "lhs", "lhs",
                  204);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 205);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 205);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 205);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 205);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 205);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 205);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 205);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 205);
  sf_mex_assign(&c3_rhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs205), "rhs", "rhs",
                  205);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs205), "lhs", "lhs",
                  205);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 206);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 206);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 206);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 206);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372554360U), "fileTimeLo",
                  "fileTimeLo", 206);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 206);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 206);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 206);
  sf_mex_assign(&c3_rhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs206), "rhs", "rhs",
                  206);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs206), "lhs", "lhs",
                  206);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 207);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_ipiv2perm"), "name",
                  "name", 207);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 207);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_ipiv2perm.m"), "resolved",
                  "resolved", 207);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786382U), "fileTimeLo",
                  "fileTimeLo", 207);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 207);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 207);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 207);
  sf_mex_assign(&c3_rhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs207), "rhs", "rhs",
                  207);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs207), "lhs", "lhs",
                  207);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_ipiv2perm.m"), "context",
                  "context", 208);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("colon"), "name", "name", 208);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 208);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 208);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1378267188U), "fileTimeLo",
                  "fileTimeLo", 208);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 208);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 208);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 208);
  sf_mex_assign(&c3_rhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs208), "rhs", "rhs",
                  208);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs208), "lhs", "lhs",
                  208);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_ipiv2perm.m"), "context",
                  "context", 209);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 209);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 209);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 209);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 209);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 209);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 209);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 209);
  sf_mex_assign(&c3_rhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs209), "rhs", "rhs",
                  209);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs209), "lhs", "lhs",
                  209);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_ipiv2perm.m"), "context",
                  "context", 210);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 210);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 210);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 210);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326692322U), "fileTimeLo",
                  "fileTimeLo", 210);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 210);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 210);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 210);
  sf_mex_assign(&c3_rhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs210), "rhs", "rhs",
                  210);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs210), "lhs", "lhs",
                  210);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 211);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 211);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 211);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 211);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 211);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 211);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 211);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 211);
  sf_mex_assign(&c3_rhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs211), "rhs", "rhs",
                  211);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs211), "lhs", "lhs",
                  211);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 212);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 212);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 212);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 212);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1372553616U), "fileTimeLo",
                  "fileTimeLo", 212);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 212);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 212);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 212);
  sf_mex_assign(&c3_rhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs212), "rhs", "rhs",
                  212);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs212), "lhs", "lhs",
                  212);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 213);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 213);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 213);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 213);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 213);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 213);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 213);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 213);
  sf_mex_assign(&c3_rhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs213), "rhs", "rhs",
                  213);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs213), "lhs", "lhs",
                  213);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 214);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xtrsm"), "name", "name",
                  214);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 214);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"),
                  "resolved", "resolved", 214);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951892U), "fileTimeLo",
                  "fileTimeLo", 214);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 214);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 214);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 214);
  sf_mex_assign(&c3_rhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs214), "rhs", "rhs",
                  214);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs214), "lhs", "lhs",
                  214);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 215);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 215);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 215);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 215);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 215);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 215);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 215);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 215);
  sf_mex_assign(&c3_rhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs215), "rhs", "rhs",
                  215);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs215), "lhs", "lhs",
                  215);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 216);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xtrsm"),
                  "name", "name", 216);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 216);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "resolved", "resolved", 216);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 216);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 216);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 216);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 216);
  sf_mex_assign(&c3_rhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs216), "rhs", "rhs",
                  216);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs216), "lhs", "lhs",
                  216);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 217);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 217);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 217);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 217);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 217);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 217);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 217);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 217);
  sf_mex_assign(&c3_rhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs217), "rhs", "rhs",
                  217);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs217), "lhs", "lhs",
                  217);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p!below_threshold"),
                  "context", "context", 218);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 218);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 218);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 218);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 218);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 218);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 218);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 218);
  sf_mex_assign(&c3_rhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs218), "rhs", "rhs",
                  218);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs218), "lhs", "lhs",
                  218);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 219);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 219);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 219);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 219);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 219);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 219);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 219);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 219);
  sf_mex_assign(&c3_rhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs219), "rhs", "rhs",
                  219);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs219), "lhs", "lhs",
                  219);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 220);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xtrsm"),
                  "name", "name", 220);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 220);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "resolved", "resolved", 220);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 220);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 220);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 220);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 220);
  sf_mex_assign(&c3_rhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs220), "rhs", "rhs",
                  220);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs220), "lhs", "lhs",
                  220);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 221);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 221);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 221);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 221);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 221);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 221);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 221);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 221);
  sf_mex_assign(&c3_rhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs221), "rhs", "rhs",
                  221);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs221), "lhs", "lhs",
                  221);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 222);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 222);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 222);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 222);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 222);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 222);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 222);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 222);
  sf_mex_assign(&c3_rhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs222), "rhs", "rhs",
                  222);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs222), "lhs", "lhs",
                  222);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 223);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 223);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 223);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 223);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 223);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 223);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 223);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 223);
  sf_mex_assign(&c3_rhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs223), "rhs", "rhs",
                  223);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs223), "lhs", "lhs",
                  223);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 224);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("intmin"), "name", "name", 224);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 224);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 224);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1362225882U), "fileTimeLo",
                  "fileTimeLo", 224);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 224);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 224);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 224);
  sf_mex_assign(&c3_rhs224, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs224, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs224), "rhs", "rhs",
                  224);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs224), "lhs", "lhs",
                  224);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 225);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("rdivide"), "name", "name", 225);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 225);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 225);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677880U), "fileTimeLo",
                  "fileTimeLo", 225);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 225);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 225);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 225);
  sf_mex_assign(&c3_rhs225, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs225, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs225), "rhs", "rhs",
                  225);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs225), "lhs", "lhs",
                  225);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 226);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("norm"), "name", "name", 226);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 226);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 226);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677868U), "fileTimeLo",
                  "fileTimeLo", 226);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 226);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 226);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 226);
  sf_mex_assign(&c3_rhs226, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs226, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs226), "rhs", "rhs",
                  226);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs226), "lhs", "lhs",
                  226);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 227);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 227);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 227);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 227);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 227);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 227);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 227);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 227);
  sf_mex_assign(&c3_rhs227, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs227, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs227), "rhs", "rhs",
                  227);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs227), "lhs", "lhs",
                  227);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 228);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("abs"), "name", "name", 228);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 228);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 228);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677852U), "fileTimeLo",
                  "fileTimeLo", 228);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 228);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 228);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 228);
  sf_mex_assign(&c3_rhs228, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs228, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs228), "rhs", "rhs",
                  228);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs228), "lhs", "lhs",
                  228);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 229);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("isnan"), "name", "name", 229);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 229);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 229);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677858U), "fileTimeLo",
                  "fileTimeLo", 229);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 229);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 229);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 229);
  sf_mex_assign(&c3_rhs229, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs229, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs229), "rhs", "rhs",
                  229);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs229), "lhs", "lhs",
                  229);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 230);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 230);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 230);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 230);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363678556U), "fileTimeLo",
                  "fileTimeLo", 230);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 230);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 230);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 230);
  sf_mex_assign(&c3_rhs230, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs230, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs230), "rhs", "rhs",
                  230);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs230), "lhs", "lhs",
                  230);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 231);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 231);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 231);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 231);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786376U), "fileTimeLo",
                  "fileTimeLo", 231);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 231);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 231);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 231);
  sf_mex_assign(&c3_rhs231, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs231, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs231), "rhs", "rhs",
                  231);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs231), "lhs", "lhs",
                  231);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 232);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 232);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 232);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 232);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786382U), "fileTimeLo",
                  "fileTimeLo", 232);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 232);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 232);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 232);
  sf_mex_assign(&c3_rhs232, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs232, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs232), "rhs", "rhs",
                  232);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs232), "lhs", "lhs",
                  232);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 233);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_warning"), "name", "name",
                  233);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 233);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 233);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786402U), "fileTimeLo",
                  "fileTimeLo", 233);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 233);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 233);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 233);
  sf_mex_assign(&c3_rhs233, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs233, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs233), "rhs", "rhs",
                  233);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs233), "lhs", "lhs",
                  233);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 234);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("isnan"), "name", "name", 234);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 234);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 234);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677858U), "fileTimeLo",
                  "fileTimeLo", 234);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 234);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 234);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 234);
  sf_mex_assign(&c3_rhs234, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs234, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs234), "rhs", "rhs",
                  234);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs234), "lhs", "lhs",
                  234);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 235);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eps"), "name", "name", 235);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 235);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 235);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1326691996U), "fileTimeLo",
                  "fileTimeLo", 235);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 235);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 235);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 235);
  sf_mex_assign(&c3_rhs235, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs235, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs235), "rhs", "rhs",
                  235);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs235), "lhs", "lhs",
                  235);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 236);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_flt2str"), "name", "name",
                  236);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 236);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "resolved",
                  "resolved", 236);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1360246350U), "fileTimeLo",
                  "fileTimeLo", 236);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 236);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 236);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 236);
  sf_mex_assign(&c3_rhs236, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs236, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs236), "rhs", "rhs",
                  236);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs236), "lhs", "lhs",
                  236);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "context",
                  "context", 237);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "name", "name", 237);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 237);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m"), "resolved",
                  "resolved", 237);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1319697568U), "fileTimeLo",
                  "fileTimeLo", 237);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 237);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 237);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 237);
  sf_mex_assign(&c3_rhs237, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs237, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs237), "rhs", "rhs",
                  237);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs237), "lhs", "lhs",
                  237);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 238);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 238);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 238);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 238);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 238);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 238);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 238);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 238);
  sf_mex_assign(&c3_rhs238, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs238, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs238), "rhs", "rhs",
                  238);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs238), "lhs", "lhs",
                  238);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 239);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 239);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 239);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 239);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951888U), "fileTimeLo",
                  "fileTimeLo", 239);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 239);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 239);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 239);
  sf_mex_assign(&c3_rhs239, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs239, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs239), "rhs", "rhs",
                  239);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs239), "lhs", "lhs",
                  239);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 240);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  240);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 240);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 240);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1375951890U), "fileTimeLo",
                  "fileTimeLo", 240);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 240);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 240);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 240);
  sf_mex_assign(&c3_rhs240, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs240, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs240), "rhs", "rhs",
                  240);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs240), "lhs", "lhs",
                  240);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 241);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 241);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 241);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 241);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 241);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 241);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 241);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 241);
  sf_mex_assign(&c3_rhs241, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs241, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs241), "rhs", "rhs",
                  241);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs241), "lhs", "lhs",
                  241);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 242);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 242);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 242);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 242);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 242);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 242);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 242);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 242);
  sf_mex_assign(&c3_rhs242, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs242, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs242), "rhs", "rhs",
                  242);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs242), "lhs", "lhs",
                  242);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 243);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 243);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 243);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 243);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 243);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 243);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 243);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 243);
  sf_mex_assign(&c3_rhs243, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs243, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs243), "rhs", "rhs",
                  243);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs243), "lhs", "lhs",
                  243);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 244);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 244);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 244);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 244);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 244);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 244);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 244);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 244);
  sf_mex_assign(&c3_rhs244, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs244, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs244), "rhs", "rhs",
                  244);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs244), "lhs", "lhs",
                  244);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 245);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 245);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 245);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 245);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 245);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 245);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 245);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 245);
  sf_mex_assign(&c3_rhs245, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs245, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs245), "rhs", "rhs",
                  245);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs245), "lhs", "lhs",
                  245);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 246);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 246);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 246);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 246);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 246);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 246);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 246);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 246);
  sf_mex_assign(&c3_rhs246, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs246, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs246), "rhs", "rhs",
                  246);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs246), "lhs", "lhs",
                  246);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 247);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 247);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 247);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 247);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 247);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 247);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 247);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 247);
  sf_mex_assign(&c3_rhs247, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs247, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs247), "rhs", "rhs",
                  247);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs247), "lhs", "lhs",
                  247);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!ceval_xgemm"),
                  "context", "context", 248);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.blas.size_ptr"),
                  "name", "name", 248);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 248);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/size_ptr.p"),
                  "resolved", "resolved", 248);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271922U), "fileTimeLo",
                  "fileTimeLo", 248);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 248);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 248);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 248);
  sf_mex_assign(&c3_rhs248, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs248, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs248), "rhs", "rhs",
                  248);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs248), "lhs", "lhs",
                  248);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!ceval_xgemm"),
                  "context", "context", 249);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.c_cast"),
                  "name", "name", 249);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("int32"), "dominantType",
                  "dominantType", 249);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/c_cast.p"),
                  "resolved", "resolved", 249);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1389271920U), "fileTimeLo",
                  "fileTimeLo", 249);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 249);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 249);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 249);
  sf_mex_assign(&c3_rhs249, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs249, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs249), "rhs", "rhs",
                  249);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs249), "lhs", "lhs",
                  249);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 250);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("isequal"), "name", "name", 250);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 250);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 250);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786358U), "fileTimeLo",
                  "fileTimeLo", 250);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 250);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 250);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 250);
  sf_mex_assign(&c3_rhs250, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs250, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs250), "rhs", "rhs",
                  250);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs250), "lhs", "lhs",
                  250);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 251);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 251);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 251);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 251);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1286786386U), "fileTimeLo",
                  "fileTimeLo", 251);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 251);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 251);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 251);
  sf_mex_assign(&c3_rhs251, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs251, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs251), "rhs", "rhs",
                  251);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs251), "lhs", "lhs",
                  251);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m!isequal_scalar"),
                  "context", "context", 252);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("isnan"), "name", "name", 252);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 252);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 252);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1363677858U), "fileTimeLo",
                  "fileTimeLo", 252);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 252);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 252);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 252);
  sf_mex_assign(&c3_rhs252, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs252, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs252), "rhs", "rhs",
                  252);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs252), "lhs", "lhs",
                  252);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m!calclen"), "context",
                  "context", 253);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 253);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 253);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 253);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1323134578U), "fileTimeLo",
                  "fileTimeLo", 253);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 253);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 253);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 253);
  sf_mex_assign(&c3_rhs253, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs253, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs253), "rhs", "rhs",
                  253);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs253), "lhs", "lhs",
                  253);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m!calclen"), "context",
                  "context", 254);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("min"), "name", "name", 254);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 254);
  sf_mex_addfield(*c3_info, c3_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 254);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(1311226518U), "fileTimeLo",
                  "fileTimeLo", 254);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 254);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 254);
  sf_mex_addfield(*c3_info, c3_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 254);
  sf_mex_assign(&c3_rhs254, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs254, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs254), "rhs", "rhs",
                  254);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs254), "lhs", "lhs",
                  254);
  sf_mex_destroy(&c3_rhs192);
  sf_mex_destroy(&c3_lhs192);
  sf_mex_destroy(&c3_rhs193);
  sf_mex_destroy(&c3_lhs193);
  sf_mex_destroy(&c3_rhs194);
  sf_mex_destroy(&c3_lhs194);
  sf_mex_destroy(&c3_rhs195);
  sf_mex_destroy(&c3_lhs195);
  sf_mex_destroy(&c3_rhs196);
  sf_mex_destroy(&c3_lhs196);
  sf_mex_destroy(&c3_rhs197);
  sf_mex_destroy(&c3_lhs197);
  sf_mex_destroy(&c3_rhs198);
  sf_mex_destroy(&c3_lhs198);
  sf_mex_destroy(&c3_rhs199);
  sf_mex_destroy(&c3_lhs199);
  sf_mex_destroy(&c3_rhs200);
  sf_mex_destroy(&c3_lhs200);
  sf_mex_destroy(&c3_rhs201);
  sf_mex_destroy(&c3_lhs201);
  sf_mex_destroy(&c3_rhs202);
  sf_mex_destroy(&c3_lhs202);
  sf_mex_destroy(&c3_rhs203);
  sf_mex_destroy(&c3_lhs203);
  sf_mex_destroy(&c3_rhs204);
  sf_mex_destroy(&c3_lhs204);
  sf_mex_destroy(&c3_rhs205);
  sf_mex_destroy(&c3_lhs205);
  sf_mex_destroy(&c3_rhs206);
  sf_mex_destroy(&c3_lhs206);
  sf_mex_destroy(&c3_rhs207);
  sf_mex_destroy(&c3_lhs207);
  sf_mex_destroy(&c3_rhs208);
  sf_mex_destroy(&c3_lhs208);
  sf_mex_destroy(&c3_rhs209);
  sf_mex_destroy(&c3_lhs209);
  sf_mex_destroy(&c3_rhs210);
  sf_mex_destroy(&c3_lhs210);
  sf_mex_destroy(&c3_rhs211);
  sf_mex_destroy(&c3_lhs211);
  sf_mex_destroy(&c3_rhs212);
  sf_mex_destroy(&c3_lhs212);
  sf_mex_destroy(&c3_rhs213);
  sf_mex_destroy(&c3_lhs213);
  sf_mex_destroy(&c3_rhs214);
  sf_mex_destroy(&c3_lhs214);
  sf_mex_destroy(&c3_rhs215);
  sf_mex_destroy(&c3_lhs215);
  sf_mex_destroy(&c3_rhs216);
  sf_mex_destroy(&c3_lhs216);
  sf_mex_destroy(&c3_rhs217);
  sf_mex_destroy(&c3_lhs217);
  sf_mex_destroy(&c3_rhs218);
  sf_mex_destroy(&c3_lhs218);
  sf_mex_destroy(&c3_rhs219);
  sf_mex_destroy(&c3_lhs219);
  sf_mex_destroy(&c3_rhs220);
  sf_mex_destroy(&c3_lhs220);
  sf_mex_destroy(&c3_rhs221);
  sf_mex_destroy(&c3_lhs221);
  sf_mex_destroy(&c3_rhs222);
  sf_mex_destroy(&c3_lhs222);
  sf_mex_destroy(&c3_rhs223);
  sf_mex_destroy(&c3_lhs223);
  sf_mex_destroy(&c3_rhs224);
  sf_mex_destroy(&c3_lhs224);
  sf_mex_destroy(&c3_rhs225);
  sf_mex_destroy(&c3_lhs225);
  sf_mex_destroy(&c3_rhs226);
  sf_mex_destroy(&c3_lhs226);
  sf_mex_destroy(&c3_rhs227);
  sf_mex_destroy(&c3_lhs227);
  sf_mex_destroy(&c3_rhs228);
  sf_mex_destroy(&c3_lhs228);
  sf_mex_destroy(&c3_rhs229);
  sf_mex_destroy(&c3_lhs229);
  sf_mex_destroy(&c3_rhs230);
  sf_mex_destroy(&c3_lhs230);
  sf_mex_destroy(&c3_rhs231);
  sf_mex_destroy(&c3_lhs231);
  sf_mex_destroy(&c3_rhs232);
  sf_mex_destroy(&c3_lhs232);
  sf_mex_destroy(&c3_rhs233);
  sf_mex_destroy(&c3_lhs233);
  sf_mex_destroy(&c3_rhs234);
  sf_mex_destroy(&c3_lhs234);
  sf_mex_destroy(&c3_rhs235);
  sf_mex_destroy(&c3_lhs235);
  sf_mex_destroy(&c3_rhs236);
  sf_mex_destroy(&c3_lhs236);
  sf_mex_destroy(&c3_rhs237);
  sf_mex_destroy(&c3_lhs237);
  sf_mex_destroy(&c3_rhs238);
  sf_mex_destroy(&c3_lhs238);
  sf_mex_destroy(&c3_rhs239);
  sf_mex_destroy(&c3_lhs239);
  sf_mex_destroy(&c3_rhs240);
  sf_mex_destroy(&c3_lhs240);
  sf_mex_destroy(&c3_rhs241);
  sf_mex_destroy(&c3_lhs241);
  sf_mex_destroy(&c3_rhs242);
  sf_mex_destroy(&c3_lhs242);
  sf_mex_destroy(&c3_rhs243);
  sf_mex_destroy(&c3_lhs243);
  sf_mex_destroy(&c3_rhs244);
  sf_mex_destroy(&c3_lhs244);
  sf_mex_destroy(&c3_rhs245);
  sf_mex_destroy(&c3_lhs245);
  sf_mex_destroy(&c3_rhs246);
  sf_mex_destroy(&c3_lhs246);
  sf_mex_destroy(&c3_rhs247);
  sf_mex_destroy(&c3_lhs247);
  sf_mex_destroy(&c3_rhs248);
  sf_mex_destroy(&c3_lhs248);
  sf_mex_destroy(&c3_rhs249);
  sf_mex_destroy(&c3_lhs249);
  sf_mex_destroy(&c3_rhs250);
  sf_mex_destroy(&c3_lhs250);
  sf_mex_destroy(&c3_rhs251);
  sf_mex_destroy(&c3_lhs251);
  sf_mex_destroy(&c3_rhs252);
  sf_mex_destroy(&c3_lhs252);
  sf_mex_destroy(&c3_rhs253);
  sf_mex_destroy(&c3_lhs253);
  sf_mex_destroy(&c3_rhs254);
  sf_mex_destroy(&c3_lhs254);
}

static real_T c3_mpower(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a)
{
  return c3_power(chartInstance, c3_a);
}

static real_T c3_power(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a)
{
  real_T c3_b_a;
  real_T c3_ak;
  real_T c3_c_a;
  c3_b_a = c3_a;
  c3_eml_scalar_eg(chartInstance);
  c3_ak = c3_b_a;
  c3_c_a = c3_ak;
  c3_eml_scalar_eg(chartInstance);
  return c3_c_a * c3_c_a;
}

static void c3_eml_scalar_eg(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c3_diag(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                    *chartInstance, real_T c3_v[30], real_T c3_d[900])
{
  int32_T c3_i239;
  int32_T c3_j;
  int32_T c3_b_j;
  (void)chartInstance;
  for (c3_i239 = 0; c3_i239 < 900; c3_i239++) {
    c3_d[c3_i239] = 0.0;
  }

  for (c3_j = 1; c3_j < 31; c3_j++) {
    c3_b_j = c3_j;
    c3_d[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_b_j), 1, 30, 1, 0) + 30 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_b_j), 1, 30, 2, 0) - 1))
      - 1] = c3_v[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c3_b_j), 1, 30, 1, 0) - 1];
  }
}

static void c3_eml_switch_helper
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_b_eml_switch_helper
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_b_power(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a[900], real_T c3_y[900])
{
  int32_T c3_k;
  real_T c3_b_k;
  real_T c3_ak;
  real_T c3_b_a;
  real_T c3_b_y;
  for (c3_k = 0; c3_k < 900; c3_k++) {
    c3_b_k = 1.0 + (real_T)c3_k;
    c3_ak = c3_a[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      c3_b_k), 1, 900, 1, 0) - 1];
    c3_b_a = c3_ak;
    c3_eml_scalar_eg(chartInstance);
    c3_b_y = c3_b_a * c3_b_a;
    c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", c3_b_k),
      1, 900, 1, 0) - 1] = c3_b_y;
  }
}

static real_T c3_b_mpower(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a)
{
  real_T c3_b_a;
  real_T c3_c_a;
  real_T c3_ak;
  real_T c3_d_a;
  real_T c3_ar;
  c3_b_a = c3_a;
  c3_c_a = c3_b_a;
  c3_eml_scalar_eg(chartInstance);
  c3_ak = c3_c_a;
  c3_d_a = c3_ak;
  c3_eml_scalar_eg(chartInstance);
  c3_ar = c3_d_a;
  return muDoubleScalarPower(c3_ar, 3.0);
}

static real_T c3_c_mpower(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_a)
{
  real_T c3_b_a;
  real_T c3_c_a;
  real_T c3_ak;
  real_T c3_d_a;
  real_T c3_ar;
  c3_b_a = c3_a;
  c3_c_a = c3_b_a;
  c3_eml_scalar_eg(chartInstance);
  c3_ak = c3_c_a;
  c3_d_a = c3_ak;
  c3_eml_scalar_eg(chartInstance);
  c3_ar = c3_d_a;
  return muDoubleScalarPower(c3_ar, 4.0);
}

static void c3_eml_error(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance)
{
  int32_T c3_i240;
  static char_T c3_cv1[48] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'c', 'h', 'o', 'l', '_', 'm', 'a', 't', 'r', 'i', 'x', 'M',
    'u', 's', 't', 'B', 'e', 'P', 'o', 's', 'D', 'e', 'f', 'W', 'i', 't', 'h',
    'R', 'e', 'a', 'l', 'D', 'i', 'a', 'g' };

  char_T c3_u[48];
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  for (c3_i240 = 0; c3_i240 < 48; c3_i240++) {
    c3_u[c3_i240] = c3_cv1[c3_i240];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 48), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c3_y));
}

static void c3_eml_matlab_zpotrf
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, real_T
   c3_A[900], real_T c3_b_A[900], int32_T *c3_info)
{
  int32_T c3_i241;
  for (c3_i241 = 0; c3_i241 < 900; c3_i241++) {
    c3_b_A[c3_i241] = c3_A[c3_i241];
  }

  *c3_info = c3_b_eml_matlab_zpotrf(chartInstance, c3_b_A);
}

static real_T c3_eml_xdotc(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_n, real_T c3_x[900], int32_T c3_ix0, real_T c3_y
  [900], int32_T c3_iy0)
{
  real_T c3_d;
  int32_T c3_b_n;
  int32_T c3_b_ix0;
  int32_T c3_b_iy0;
  int32_T c3_c_n;
  int32_T c3_c_ix0;
  int32_T c3_c_iy0;
  int32_T c3_d_n;
  int32_T c3_d_ix0;
  int32_T c3_d_iy0;
  int32_T c3_e_n;
  int32_T c3_e_ix0;
  int32_T c3_e_iy0;
  int32_T c3_ix;
  int32_T c3_iy;
  int32_T c3_f_n;
  int32_T c3_b;
  int32_T c3_b_b;
  boolean_T c3_overflow;
  int32_T c3_k;
  int32_T c3_a;
  int32_T c3_b_a;
  c3_b_n = c3_n;
  c3_b_ix0 = c3_ix0;
  c3_b_iy0 = c3_iy0;
  c3_c_n = c3_b_n;
  c3_c_ix0 = c3_b_ix0;
  c3_c_iy0 = c3_b_iy0;
  c3_d_n = c3_c_n;
  c3_d_ix0 = c3_c_ix0;
  c3_d_iy0 = c3_c_iy0;
  c3_e_n = c3_d_n;
  c3_e_ix0 = c3_d_ix0;
  c3_e_iy0 = c3_d_iy0;
  c3_d = 0.0;
  if (c3_e_n < 1) {
  } else {
    c3_ix = c3_e_ix0;
    c3_iy = c3_e_iy0;
    c3_f_n = c3_e_n;
    c3_b = c3_f_n;
    c3_b_b = c3_b;
    if (1 > c3_b_b) {
      c3_overflow = false;
    } else {
      c3_eml_switch_helper(chartInstance);
      c3_overflow = (c3_b_b > 2147483646);
    }

    if (c3_overflow) {
      c3_check_forloop_overflow_error(chartInstance, c3_overflow);
    }

    for (c3_k = 1; c3_k <= c3_f_n; c3_k++) {
      c3_d += c3_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c3_ix), 1, 900, 1, 0) - 1] *
        c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c3_iy), 1, 900, 1, 0) - 1];
      c3_a = c3_ix + 1;
      c3_ix = c3_a;
      c3_b_a = c3_iy + 1;
      c3_iy = c3_b_a;
    }
  }

  return c3_d;
}

static void c3_check_forloop_overflow_error
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, boolean_T
   c3_overflow)
{
  int32_T c3_i242;
  static char_T c3_cv2[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c3_u[34];
  const mxArray *c3_y = NULL;
  int32_T c3_i243;
  static char_T c3_cv3[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c3_b_u[23];
  const mxArray *c3_b_y = NULL;
  (void)chartInstance;
  if (!c3_overflow) {
  } else {
    for (c3_i242 = 0; c3_i242 < 34; c3_i242++) {
      c3_u[c3_i242] = c3_cv2[c3_i242];
    }

    c3_y = NULL;
    sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  false);
    for (c3_i243 = 0; c3_i243 < 23; c3_i243++) {
      c3_b_u[c3_i243] = c3_cv3[c3_i243];
    }

    c3_b_y = NULL;
    sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 2U, 14, c3_y, 14, c3_b_y));
  }
}

static void c3_eml_xgemv(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_m, int32_T c3_n, int32_T c3_ia0, int32_T c3_ix0,
  real_T c3_y[900], int32_T c3_iy0, real_T c3_b_y[900])
{
  int32_T c3_i244;
  for (c3_i244 = 0; c3_i244 < 900; c3_i244++) {
    c3_b_y[c3_i244] = c3_y[c3_i244];
  }

  c3_b_eml_xgemv(chartInstance, c3_m, c3_n, c3_ia0, c3_ix0, c3_b_y, c3_iy0);
}

static void c3_below_threshold
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_b_below_threshold
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_b_eml_error(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance)
{
  int32_T c3_i245;
  static char_T c3_cv4[19] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T', 'L',
    'A', 'B', ':', 'p', 'o', 's', 'd', 'e', 'f' };

  char_T c3_u[19];
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  for (c3_i245 = 0; c3_i245 < 19; c3_i245++) {
    c3_u[c3_i245] = c3_cv4[c3_i245];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 19), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c3_y));
}

static void c3_inv(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                   *chartInstance, real_T c3_x[36], real_T c3_y[36])
{
  int32_T c3_i246;
  real_T c3_b_x[36];
  int32_T c3_i247;
  int32_T c3_info;
  int32_T c3_ipiv[6];
  int32_T c3_i248;
  int32_T c3_b_ipiv[6];
  int32_T c3_k;
  int32_T c3_b_k;
  int32_T c3_c;
  int32_T c3_c_k;
  int32_T c3_a;
  int32_T c3_b_a;
  boolean_T c3_overflow;
  int32_T c3_j;
  int32_T c3_b_j;
  int32_T c3_c_a;
  int32_T c3_d_a;
  int32_T c3_i249;
  int32_T c3_e_a;
  int32_T c3_f_a;
  boolean_T c3_b_overflow;
  int32_T c3_i;
  int32_T c3_b_i;
  int32_T c3_i250;
  real_T c3_c_x[36];
  int32_T c3_i251;
  real_T c3_d_x[36];
  real_T c3_n1x;
  int32_T c3_i252;
  real_T c3_b_y[36];
  real_T c3_n1xinv;
  real_T c3_rc;
  real_T c3_e_x;
  boolean_T c3_b;
  real_T c3_f_x;
  int32_T c3_i253;
  static char_T c3_cv5[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c3_u[8];
  const mxArray *c3_c_y = NULL;
  real_T c3_b_u;
  const mxArray *c3_d_y = NULL;
  real_T c3_c_u;
  const mxArray *c3_e_y = NULL;
  real_T c3_d_u;
  const mxArray *c3_f_y = NULL;
  char_T c3_str[14];
  int32_T c3_i254;
  char_T c3_b_str[14];
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c3_i246 = 0; c3_i246 < 36; c3_i246++) {
    c3_b_x[c3_i246] = c3_x[c3_i246];
  }

  for (c3_i247 = 0; c3_i247 < 36; c3_i247++) {
    c3_y[c3_i247] = 0.0;
  }

  c3_b_eml_matlab_zgetrf(chartInstance, c3_b_x, c3_ipiv, &c3_info);
  for (c3_i248 = 0; c3_i248 < 6; c3_i248++) {
    c3_b_ipiv[c3_i248] = c3_ipiv[c3_i248];
  }

  c3_eml_ipiv2perm(chartInstance, c3_b_ipiv, c3_ipiv);
  for (c3_k = 1; c3_k < 7; c3_k++) {
    c3_b_k = c3_k;
    c3_c = c3_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
      "", (real_T)c3_b_k), 1, 6, 1, 0) - 1];
    c3_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_b_k), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_c), 1, 6, 2, 0) - 1)) - 1]
      = 1.0;
    c3_c_k = c3_b_k;
    c3_a = c3_c_k;
    c3_b_a = c3_a;
    if (c3_b_a > 6) {
      c3_overflow = false;
    } else {
      c3_eml_switch_helper(chartInstance);
      c3_overflow = false;
    }

    if (c3_overflow) {
      c3_check_forloop_overflow_error(chartInstance, c3_overflow);
    }

    for (c3_j = c3_c_k; c3_j < 7; c3_j++) {
      c3_b_j = c3_j;
      if (c3_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
             (real_T)c3_b_j), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
             (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_c), 1, 6, 2, 0) - 1)) -
          1] != 0.0) {
        c3_c_a = c3_b_j;
        c3_d_a = c3_c_a + 1;
        c3_i249 = c3_d_a;
        c3_e_a = c3_i249;
        c3_f_a = c3_e_a;
        if (c3_f_a > 6) {
          c3_b_overflow = false;
        } else {
          c3_eml_switch_helper(chartInstance);
          c3_b_overflow = false;
        }

        if (c3_b_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_b_overflow);
        }

        for (c3_i = c3_i249; c3_i < 7; c3_i++) {
          c3_b_i = c3_i;
          c3_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c3_b_i), 1, 6, 1, 0) + 6 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c3_c), 1, 6, 2, 0) - 1)) - 1] = c3_y
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c3_b_i), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_c), 1, 6, 2, 0) -
               1)) - 1] - c3_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c3_b_j), 1, 6, 1, 0) + 6 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c3_c), 1, 6, 2, 0) - 1)) - 1] *
            c3_b_x[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c3_b_i), 1, 6, 1, 0) + 6 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c3_b_j), 1, 6, 2, 0) - 1)) - 1];
        }
      }
    }
  }

  for (c3_i250 = 0; c3_i250 < 36; c3_i250++) {
    c3_c_x[c3_i250] = c3_b_x[c3_i250];
  }

  c3_b_eml_xtrsm(chartInstance, c3_c_x, c3_y);
  for (c3_i251 = 0; c3_i251 < 36; c3_i251++) {
    c3_d_x[c3_i251] = c3_x[c3_i251];
  }

  c3_n1x = c3_norm(chartInstance, c3_d_x);
  for (c3_i252 = 0; c3_i252 < 36; c3_i252++) {
    c3_b_y[c3_i252] = c3_y[c3_i252];
  }

  c3_n1xinv = c3_norm(chartInstance, c3_b_y);
  c3_rc = 1.0 / (c3_n1x * c3_n1xinv);
  guard1 = false;
  guard2 = false;
  if (c3_n1x == 0.0) {
    guard2 = true;
  } else if (c3_n1xinv == 0.0) {
    guard2 = true;
  } else if (c3_rc == 0.0) {
    guard1 = true;
  } else {
    c3_e_x = c3_rc;
    c3_b = muDoubleScalarIsNaN(c3_e_x);
    guard3 = false;
    if (c3_b) {
      guard3 = true;
    } else {
      c3_eps(chartInstance);
      if (c3_rc < 2.2204460492503131E-16) {
        guard3 = true;
      }
    }

    if (guard3 == true) {
      c3_f_x = c3_rc;
      for (c3_i253 = 0; c3_i253 < 8; c3_i253++) {
        c3_u[c3_i253] = c3_cv5[c3_i253];
      }

      c3_c_y = NULL;
      sf_mex_assign(&c3_c_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 8),
                    false);
      c3_b_u = 14.0;
      c3_d_y = NULL;
      sf_mex_assign(&c3_d_y, sf_mex_create("y", &c3_b_u, 0, 0U, 0U, 0U, 0),
                    false);
      c3_c_u = 6.0;
      c3_e_y = NULL;
      sf_mex_assign(&c3_e_y, sf_mex_create("y", &c3_c_u, 0, 0U, 0U, 0U, 0),
                    false);
      c3_d_u = c3_f_x;
      c3_f_y = NULL;
      sf_mex_assign(&c3_f_y, sf_mex_create("y", &c3_d_u, 0, 0U, 0U, 0U, 0),
                    false);
      c3_emlrt_marshallIn(chartInstance, sf_mex_call_debug
                          (sfGlobalDebugInstanceStruct, "sprintf", 1U, 2U, 14,
                           sf_mex_call_debug(sfGlobalDebugInstanceStruct,
        "sprintf", 1U, 3U, 14, c3_c_y, 14, c3_d_y, 14, c3_e_y), 14, c3_f_y),
                          "sprintf", c3_str);
      for (c3_i254 = 0; c3_i254 < 14; c3_i254++) {
        c3_b_str[c3_i254] = c3_str[c3_i254];
      }

      c3_b_eml_warning(chartInstance, c3_b_str);
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c3_eml_warning(chartInstance);
  }
}

static void c3_eps(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                   *chartInstance)
{
  (void)chartInstance;
}

static void c3_eml_matlab_zgetrf
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, real_T
   c3_A[36], real_T c3_b_A[36], int32_T c3_ipiv[6], int32_T *c3_info)
{
  int32_T c3_i255;
  for (c3_i255 = 0; c3_i255 < 36; c3_i255++) {
    c3_b_A[c3_i255] = c3_A[c3_i255];
  }

  c3_b_eml_matlab_zgetrf(chartInstance, c3_b_A, c3_ipiv, c3_info);
}

static int32_T c3_eml_ixamax(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_n, real_T c3_x[36], int32_T c3_ix0)
{
  int32_T c3_idxmax;
  int32_T c3_b_n;
  int32_T c3_b_ix0;
  int32_T c3_c_n;
  int32_T c3_c_ix0;
  int32_T c3_ix;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_d_x;
  real_T c3_y;
  real_T c3_e_x;
  real_T c3_f_x;
  real_T c3_b_y;
  real_T c3_smax;
  int32_T c3_d_n;
  int32_T c3_b;
  int32_T c3_b_b;
  boolean_T c3_overflow;
  int32_T c3_k;
  int32_T c3_b_k;
  int32_T c3_a;
  real_T c3_g_x;
  real_T c3_h_x;
  real_T c3_i_x;
  real_T c3_c_y;
  real_T c3_j_x;
  real_T c3_k_x;
  real_T c3_d_y;
  real_T c3_s;
  c3_b_n = c3_n;
  c3_b_ix0 = c3_ix0;
  c3_c_n = c3_b_n;
  c3_c_ix0 = c3_b_ix0;
  if (c3_c_n < 1) {
    c3_idxmax = 0;
  } else {
    c3_idxmax = 1;
    if (c3_c_n > 1) {
      c3_ix = c3_c_ix0;
      c3_b_x = c3_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c3_ix), 1, 36, 1, 0) - 1];
      c3_c_x = c3_b_x;
      c3_d_x = c3_c_x;
      c3_y = muDoubleScalarAbs(c3_d_x);
      c3_e_x = 0.0;
      c3_f_x = c3_e_x;
      c3_b_y = muDoubleScalarAbs(c3_f_x);
      c3_smax = c3_y + c3_b_y;
      c3_d_n = c3_c_n;
      c3_b = c3_d_n;
      c3_b_b = c3_b;
      if (2 > c3_b_b) {
        c3_overflow = false;
      } else {
        c3_eml_switch_helper(chartInstance);
        c3_overflow = (c3_b_b > 2147483646);
      }

      if (c3_overflow) {
        c3_check_forloop_overflow_error(chartInstance, c3_overflow);
      }

      for (c3_k = 2; c3_k <= c3_d_n; c3_k++) {
        c3_b_k = c3_k;
        c3_a = c3_ix + 1;
        c3_ix = c3_a;
        c3_g_x = c3_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c3_ix), 1, 36, 1, 0) - 1];
        c3_h_x = c3_g_x;
        c3_i_x = c3_h_x;
        c3_c_y = muDoubleScalarAbs(c3_i_x);
        c3_j_x = 0.0;
        c3_k_x = c3_j_x;
        c3_d_y = muDoubleScalarAbs(c3_k_x);
        c3_s = c3_c_y + c3_d_y;
        if (c3_s > c3_smax) {
          c3_idxmax = c3_b_k;
          c3_smax = c3_s;
        }
      }
    }
  }

  return c3_idxmax;
}

static void c3_threshold(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c3_eml_xgeru(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_m, int32_T c3_n, real_T c3_alpha1, int32_T c3_ix0,
  int32_T c3_iy0, real_T c3_A[36], int32_T c3_ia0, real_T c3_b_A[36])
{
  int32_T c3_i256;
  for (c3_i256 = 0; c3_i256 < 36; c3_i256++) {
    c3_b_A[c3_i256] = c3_A[c3_i256];
  }

  c3_b_eml_xgeru(chartInstance, c3_m, c3_n, c3_alpha1, c3_ix0, c3_iy0, c3_b_A,
                 c3_ia0);
}

static void c3_eml_ipiv2perm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_ipiv[6], int32_T c3_perm[6])
{
  int32_T c3_i257;
  int32_T c3_k;
  real_T c3_b_k;
  int32_T c3_ipk;
  int32_T c3_a;
  real_T c3_b;
  int32_T c3_b_a;
  real_T c3_b_b;
  int32_T c3_idx;
  real_T c3_flt;
  boolean_T c3_p;
  int32_T c3_pipk;
  for (c3_i257 = 0; c3_i257 < 6; c3_i257++) {
    c3_perm[c3_i257] = 1 + c3_i257;
  }

  for (c3_k = 0; c3_k < 5; c3_k++) {
    c3_b_k = 1.0 + (real_T)c3_k;
    c3_ipk = c3_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", c3_b_k), 1, 6, 1, 0) - 1];
    c3_a = c3_ipk;
    c3_b = c3_b_k;
    c3_b_a = c3_a;
    c3_b_b = c3_b;
    c3_b_eml_switch_helper(chartInstance);
    c3_idx = c3_b_a;
    c3_flt = c3_b_b;
    c3_p = ((real_T)c3_idx > c3_flt);
    if (c3_p) {
      c3_pipk = c3_perm[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c3_ipk), 1, 6, 1, 0) - 1];
      c3_perm[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c3_ipk), 1, 6, 1, 0) - 1] = c3_perm[_SFD_EML_ARRAY_BOUNDS_CHECK(
        "", (int32_T)_SFD_INTEGER_CHECK("", c3_b_k), 1, 6, 1, 0) - 1];
      c3_perm[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c3_b_k), 1, 6, 1, 0) - 1] = c3_pipk;
    }
  }
}

static void c3_eml_xtrsm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[36], real_T c3_B[36], real_T c3_b_B[36])
{
  int32_T c3_i258;
  int32_T c3_i259;
  real_T c3_b_A[36];
  for (c3_i258 = 0; c3_i258 < 36; c3_i258++) {
    c3_b_B[c3_i258] = c3_B[c3_i258];
  }

  for (c3_i259 = 0; c3_i259 < 36; c3_i259++) {
    c3_b_A[c3_i259] = c3_A[c3_i259];
  }

  c3_b_eml_xtrsm(chartInstance, c3_b_A, c3_b_B);
}

static void c3_b_threshold(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static real_T c3_norm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                      *chartInstance, real_T c3_x[36])
{
  real_T c3_y;
  int32_T c3_j;
  real_T c3_b_j;
  real_T c3_s;
  int32_T c3_i;
  real_T c3_b_i;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_b_y;
  real_T c3_d_x;
  boolean_T c3_b;
  boolean_T exitg1;
  (void)chartInstance;
  c3_y = 0.0;
  c3_j = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c3_j < 6)) {
    c3_b_j = 1.0 + (real_T)c3_j;
    c3_s = 0.0;
    for (c3_i = 0; c3_i < 6; c3_i++) {
      c3_b_i = 1.0 + (real_T)c3_i;
      c3_b_x = c3_x[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c3_b_i), 1, 6, 1, 0) + 6 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c3_b_j), 1, 6, 2, 0) - 1)) - 1];
      c3_c_x = c3_b_x;
      c3_b_y = muDoubleScalarAbs(c3_c_x);
      c3_s += c3_b_y;
    }

    c3_d_x = c3_s;
    c3_b = muDoubleScalarIsNaN(c3_d_x);
    if (c3_b) {
      c3_y = rtNaN;
      exitg1 = true;
    } else {
      if (c3_s > c3_y) {
        c3_y = c3_s;
      }

      c3_j++;
    }
  }

  return c3_y;
}

static void c3_eml_warning(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance)
{
  int32_T c3_i260;
  static char_T c3_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c3_u[27];
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  for (c3_i260 = 0; c3_i260 < 27; c3_i260++) {
    c3_u[c3_i260] = c3_varargin_1[c3_i260];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 27), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c3_y));
}

static void c3_b_eml_warning(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, char_T c3_varargin_2[14])
{
  int32_T c3_i261;
  static char_T c3_varargin_1[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 'i', 'l', 'l', 'C', 'o', 'n', 'd', 'i', 't', 'i',
    'o', 'n', 'e', 'd', 'M', 'a', 't', 'r', 'i', 'x' };

  char_T c3_u[33];
  const mxArray *c3_y = NULL;
  int32_T c3_i262;
  char_T c3_b_u[14];
  const mxArray *c3_b_y = NULL;
  (void)chartInstance;
  for (c3_i261 = 0; c3_i261 < 33; c3_i261++) {
    c3_u[c3_i261] = c3_varargin_1[c3_i261];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 33), false);
  for (c3_i262 = 0; c3_i262 < 14; c3_i262++) {
    c3_b_u[c3_i262] = c3_varargin_2[c3_i262];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_b_u, 10, 0U, 1U, 0U, 2, 1, 14),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c3_y, 14, c3_b_y));
}

static void c3_b_eml_scalar_eg
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_eml_xgemm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[180], real_T c3_B[36], real_T c3_C[180], real_T
  c3_b_C[180])
{
  int32_T c3_i263;
  int32_T c3_i264;
  real_T c3_b_A[180];
  int32_T c3_i265;
  real_T c3_b_B[36];
  for (c3_i263 = 0; c3_i263 < 180; c3_i263++) {
    c3_b_C[c3_i263] = c3_C[c3_i263];
  }

  for (c3_i264 = 0; c3_i264 < 180; c3_i264++) {
    c3_b_A[c3_i264] = c3_A[c3_i264];
  }

  for (c3_i265 = 0; c3_i265 < 36; c3_i265++) {
    c3_b_B[c3_i265] = c3_B[c3_i265];
  }

  c3_c_eml_xgemm(chartInstance, c3_b_A, c3_b_B, c3_b_C);
}

static void c3_c_threshold(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c3_c_eml_scalar_eg
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_d_eml_scalar_eg
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c3_b_eml_xgemm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[180], real_T c3_B[180], real_T c3_C[900], real_T
  c3_b_C[900])
{
  int32_T c3_i266;
  int32_T c3_i267;
  real_T c3_b_A[180];
  int32_T c3_i268;
  real_T c3_b_B[180];
  for (c3_i266 = 0; c3_i266 < 900; c3_i266++) {
    c3_b_C[c3_i266] = c3_C[c3_i266];
  }

  for (c3_i267 = 0; c3_i267 < 180; c3_i267++) {
    c3_b_A[c3_i267] = c3_A[c3_i267];
  }

  for (c3_i268 = 0; c3_i268 < 180; c3_i268++) {
    c3_b_B[c3_i268] = c3_B[c3_i268];
  }

  c3_d_eml_xgemm(chartInstance, c3_b_A, c3_b_B, c3_b_C);
}

static void c3_b_diag(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                      *chartInstance, real_T c3_v[36], real_T c3_d[6])
{
  int32_T c3_j;
  int32_T c3_b_j;
  (void)chartInstance;
  for (c3_j = 1; c3_j < 7; c3_j++) {
    c3_b_j = c3_j;
    c3_d[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c3_b_j), 1, 6, 1, 0) - 1] = c3_v[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)(1 + (c3_b_j - 1) * 7)), 1, 36, 1, 0) - 1];
  }
}

static void c3_c_diag(SFc3_aircraftControl_FullStateFiltersInstanceStruct
                      *chartInstance, real_T c3_v[900], real_T c3_d[30])
{
  int32_T c3_j;
  int32_T c3_b_j;
  (void)chartInstance;
  for (c3_j = 1; c3_j < 31; c3_j++) {
    c3_b_j = c3_j;
    c3_d[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
      c3_b_j), 1, 30, 1, 0) - 1] = c3_v[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)(1 + (c3_b_j - 1) * 31)), 1, 900, 1, 0) - 1];
  }
}

static const mxArray *c3_p_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(int32_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, false);
  return c3_mxArrayOutData;
}

static int32_T c3_x_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  int32_T c3_y;
  int32_T c3_i269;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_i269, 1, 6, 0U, 0, 0U, 0);
  c3_y = c3_i269;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_p_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_sfEvent;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  int32_T c3_y;
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    chartInstanceVoid;
  c3_b_sfEvent = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_x_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_sfEvent),
    &c3_thisId);
  sf_mex_destroy(&c3_b_sfEvent);
  *(int32_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static uint8_T c3_y_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_b_is_active_c3_aircraftControl_FullStateFilters, const char_T
   *c3_identifier)
{
  uint8_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_ab_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_b_is_active_c3_aircraftControl_FullStateFilters), &c3_thisId);
  sf_mex_destroy(&c3_b_is_active_c3_aircraftControl_FullStateFilters);
  return c3_y;
}

static uint8_T c3_ab_emlrt_marshallIn
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, const
   mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  uint8_T c3_y;
  uint8_T c3_u0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_u0, 1, 3, 0U, 0, 0U, 0);
  c3_y = c3_u0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static int32_T c3_b_eml_matlab_zpotrf
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, real_T
   c3_A[900])
{
  int32_T c3_info;
  int32_T c3_colj;
  int32_T c3_j;
  int32_T c3_b_j;
  int32_T c3_a;
  int32_T c3_b_a;
  int32_T c3_jm1;
  int32_T c3_c_a;
  int32_T c3_b;
  int32_T c3_d_a;
  int32_T c3_b_b;
  int32_T c3_jj;
  int32_T c3_i270;
  int32_T c3_i271;
  int32_T c3_i272;
  real_T c3_b_A[900];
  int32_T c3_i273;
  int32_T c3_i274;
  int32_T c3_i275;
  real_T c3_c_A[900];
  real_T c3_ajj;
  int32_T c3_c_b;
  int32_T c3_d_b;
  int32_T c3_nmj;
  int32_T c3_e_a;
  int32_T c3_f_a;
  int32_T c3_jjp1;
  int32_T c3_g_a;
  int32_T c3_h_a;
  int32_T c3_coljp1;
  int32_T c3_b_jm1;
  int32_T c3_e_b;
  int32_T c3_f_b;
  boolean_T c3_overflow;
  int32_T c3_k;
  int32_T c3_b_k;
  int32_T c3_c_jm1;
  int32_T c3_g_b;
  int32_T c3_h_b;
  boolean_T c3_b_overflow;
  int32_T c3_c_k;
  real_T c3_y;
  real_T c3_b_y;
  real_T c3_z;
  int32_T c3_n;
  real_T c3_i_a;
  int32_T c3_ix0;
  int32_T c3_b_n;
  real_T c3_j_a;
  int32_T c3_b_ix0;
  int32_T c3_c_n;
  real_T c3_k_a;
  int32_T c3_c_ix0;
  int32_T c3_d_ix0;
  int32_T c3_l_a;
  int32_T c3_c;
  int32_T c3_i_b;
  int32_T c3_b_c;
  int32_T c3_m_a;
  int32_T c3_j_b;
  int32_T c3_i276;
  int32_T c3_n_a;
  int32_T c3_k_b;
  int32_T c3_o_a;
  int32_T c3_l_b;
  boolean_T c3_c_overflow;
  int32_T c3_d_k;
  int32_T c3_e_k;
  boolean_T exitg1;
  c3_info = 0;
  c3_colj = 1;
  c3_j = 1;
  exitg1 = false;
  while ((exitg1 == false) && (c3_j < 31)) {
    c3_b_j = c3_j;
    c3_a = c3_b_j;
    c3_b_a = c3_a - 1;
    c3_jm1 = c3_b_a;
    c3_c_a = c3_colj;
    c3_b = c3_jm1;
    c3_d_a = c3_c_a;
    c3_b_b = c3_b;
    c3_jj = c3_d_a + c3_b_b;
    c3_i270 = 0;
    for (c3_i271 = 0; c3_i271 < 30; c3_i271++) {
      for (c3_i272 = 0; c3_i272 < 30; c3_i272++) {
        c3_b_A[c3_i272 + c3_i270] = c3_A[c3_i272 + c3_i270];
      }

      c3_i270 += 30;
    }

    c3_i273 = 0;
    for (c3_i274 = 0; c3_i274 < 30; c3_i274++) {
      for (c3_i275 = 0; c3_i275 < 30; c3_i275++) {
        c3_c_A[c3_i275 + c3_i273] = c3_A[c3_i275 + c3_i273];
      }

      c3_i273 += 30;
    }

    c3_ajj = c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c3_jj), 1, 900, 1, 0) - 1] - c3_eml_xdotc(chartInstance, c3_jm1,
      c3_b_A, c3_colj, c3_c_A, c3_colj);
    if (c3_ajj > 0.0) {
      c3_ajj = muDoubleScalarSqrt(c3_ajj);
      c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c3_jj), 1, 900, 1, 0) - 1] = c3_ajj;
      if (c3_b_j < 30) {
        c3_c_b = c3_b_j;
        c3_d_b = c3_c_b;
        c3_nmj = 30 - c3_d_b;
        c3_e_a = c3_jj;
        c3_f_a = c3_e_a + 30;
        c3_jjp1 = c3_f_a;
        c3_g_a = c3_colj;
        c3_h_a = c3_g_a + 30;
        c3_coljp1 = c3_h_a;
        c3_b_jm1 = c3_jm1;
        c3_e_b = c3_b_jm1;
        c3_f_b = c3_e_b;
        if (1 > c3_f_b) {
          c3_overflow = false;
        } else {
          c3_eml_switch_helper(chartInstance);
          c3_overflow = (c3_f_b > 2147483646);
        }

        if (c3_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_overflow);
        }

        for (c3_k = 1; c3_k <= c3_b_jm1; c3_k++) {
          c3_b_k = c3_k;
          c3_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c3_b_k), 1, 30, 1, 0) + 30 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c3_b_j), 1, 30, 2, 0) - 1)) - 1] = c3_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c3_b_k), 1, 30, 1, 0) + 30 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c3_b_j), 1, 30, 2, 0) - 1)) - 1];
        }

        c3_b_eml_xgemv(chartInstance, c3_jm1, c3_nmj, c3_coljp1, c3_colj, c3_A,
                       c3_jjp1);
        c3_c_jm1 = c3_jm1;
        c3_g_b = c3_c_jm1;
        c3_h_b = c3_g_b;
        if (1 > c3_h_b) {
          c3_b_overflow = false;
        } else {
          c3_eml_switch_helper(chartInstance);
          c3_b_overflow = (c3_h_b > 2147483646);
        }

        if (c3_b_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_b_overflow);
        }

        for (c3_c_k = 1; c3_c_k <= c3_c_jm1; c3_c_k++) {
          c3_b_k = c3_c_k;
          c3_A[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c3_b_k), 1, 30, 1, 0) + 30 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c3_b_j), 1, 30, 2, 0) - 1)) - 1] = c3_A
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c3_b_k), 1, 30, 1, 0) + 30 *
              (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                 (real_T)c3_b_j), 1, 30, 2, 0) - 1)) - 1];
        }

        c3_y = c3_ajj;
        c3_b_y = c3_y;
        c3_z = 1.0 / c3_b_y;
        c3_n = c3_nmj;
        c3_i_a = c3_z;
        c3_ix0 = c3_jjp1;
        c3_b_n = c3_n;
        c3_j_a = c3_i_a;
        c3_b_ix0 = c3_ix0;
        c3_b_below_threshold(chartInstance);
        c3_c_n = c3_b_n;
        c3_k_a = c3_j_a;
        c3_c_ix0 = c3_b_ix0;
        c3_d_ix0 = c3_c_ix0;
        c3_l_a = c3_c_n;
        c3_c = c3_l_a;
        c3_i_b = c3_c - 1;
        c3_b_c = 30 * c3_i_b;
        c3_m_a = c3_c_ix0;
        c3_j_b = c3_b_c;
        c3_i276 = c3_m_a + c3_j_b;
        c3_n_a = c3_d_ix0;
        c3_k_b = c3_i276;
        c3_o_a = c3_n_a;
        c3_l_b = c3_k_b;
        if (c3_o_a > c3_l_b) {
          c3_c_overflow = false;
        } else {
          c3_eml_switch_helper(chartInstance);
          c3_c_overflow = (c3_l_b > 2147483617);
        }

        if (c3_c_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_c_overflow);
        }

        for (c3_d_k = c3_d_ix0; c3_d_k <= c3_i276; c3_d_k += 30) {
          c3_e_k = c3_d_k;
          c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_e_k), 1, 900, 1, 0) - 1] = c3_k_a *
            c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_e_k), 1, 900, 1, 0) - 1];
        }

        c3_colj = c3_coljp1;
      }

      c3_j++;
    } else {
      c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c3_jj), 1, 900, 1, 0) - 1] = c3_ajj;
      c3_info = c3_b_j;
      exitg1 = true;
    }
  }

  return c3_info;
}

static void c3_b_eml_xgemv(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_m, int32_T c3_n, int32_T c3_ia0, int32_T c3_ix0,
  real_T c3_y[900], int32_T c3_iy0)
{
  int32_T c3_b_m;
  int32_T c3_b_n;
  int32_T c3_b_ia0;
  int32_T c3_b_ix0;
  int32_T c3_b_iy0;
  int32_T c3_c_m;
  int32_T c3_c_n;
  real_T c3_alpha1;
  int32_T c3_c_ia0;
  int32_T c3_c_ix0;
  real_T c3_beta1;
  int32_T c3_c_iy0;
  char_T c3_TRANSA;
  int32_T c3_var;
  ptrdiff_t c3_m_t;
  int32_T c3_b_var;
  ptrdiff_t c3_n_t;
  ptrdiff_t c3_lda_t;
  ptrdiff_t c3_incx_t;
  ptrdiff_t c3_incy_t;
  double * c3_alpha1_t;
  double * c3_beta1_t;
  double * c3_yiy0_t;
  double * c3_yix0_t;
  double * c3_yia0_t;
  c3_b_m = c3_m;
  c3_b_n = c3_n;
  c3_b_ia0 = c3_ia0;
  c3_b_ix0 = c3_ix0;
  c3_b_iy0 = c3_iy0;
  c3_below_threshold(chartInstance);
  if (c3_b_m < 1) {
  } else if (c3_b_n < 1) {
  } else {
    c3_c_m = c3_b_m;
    c3_c_n = c3_b_n;
    c3_alpha1 = -1.0;
    c3_c_ia0 = c3_b_ia0;
    c3_c_ix0 = c3_b_ix0;
    c3_beta1 = 1.0;
    c3_c_iy0 = c3_b_iy0;
    c3_TRANSA = 'T';
    c3_var = c3_c_m;
    c3_m_t = (ptrdiff_t)(c3_var);
    c3_b_var = c3_c_n;
    c3_n_t = (ptrdiff_t)(c3_b_var);
    c3_lda_t = (ptrdiff_t)(30);
    c3_incx_t = (ptrdiff_t)(1);
    c3_incy_t = (ptrdiff_t)(30);
    c3_alpha1_t = (double *)(&c3_alpha1);
    c3_beta1_t = (double *)(&c3_beta1);
    c3_yiy0_t = (double *)(&c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c3_c_iy0), 1, 900, 1, 0) - 1]);
    c3_yix0_t = (double *)(&c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c3_c_ix0), 1, 900, 1, 0) - 1]);
    c3_yia0_t = (double *)(&c3_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c3_c_ia0), 1, 900, 1, 0) - 1]);
    dgemv(&c3_TRANSA, &c3_m_t, &c3_n_t, c3_alpha1_t, c3_yia0_t, &c3_lda_t,
          c3_yix0_t, &c3_incx_t, c3_beta1_t, c3_yiy0_t, &c3_incy_t);
  }
}

static void c3_b_eml_matlab_zgetrf
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance, real_T
   c3_A[36], int32_T c3_ipiv[6], int32_T *c3_info)
{
  int32_T c3_i277;
  int32_T c3_j;
  int32_T c3_b_j;
  int32_T c3_a;
  int32_T c3_b_a;
  int32_T c3_jm1;
  int32_T c3_b;
  int32_T c3_b_b;
  int32_T c3_mmj;
  int32_T c3_c_a;
  int32_T c3_d_a;
  int32_T c3_c;
  int32_T c3_c_b;
  int32_T c3_d_b;
  int32_T c3_jj;
  int32_T c3_e_a;
  int32_T c3_f_a;
  int32_T c3_jp1j;
  int32_T c3_g_a;
  int32_T c3_h_a;
  int32_T c3_b_c;
  int32_T c3_i278;
  int32_T c3_i279;
  int32_T c3_i280;
  real_T c3_b_A[36];
  int32_T c3_i_a;
  int32_T c3_j_a;
  int32_T c3_jpiv_offset;
  int32_T c3_k_a;
  int32_T c3_e_b;
  int32_T c3_l_a;
  int32_T c3_f_b;
  int32_T c3_jpiv;
  int32_T c3_m_a;
  int32_T c3_g_b;
  int32_T c3_n_a;
  int32_T c3_h_b;
  int32_T c3_c_c;
  int32_T c3_i_b;
  int32_T c3_j_b;
  int32_T c3_jrow;
  int32_T c3_o_a;
  int32_T c3_k_b;
  int32_T c3_p_a;
  int32_T c3_l_b;
  int32_T c3_jprow;
  int32_T c3_ix0;
  int32_T c3_iy0;
  int32_T c3_b_ix0;
  int32_T c3_b_iy0;
  int32_T c3_c_ix0;
  int32_T c3_c_iy0;
  int32_T c3_ix;
  int32_T c3_iy;
  int32_T c3_k;
  real_T c3_temp;
  int32_T c3_q_a;
  int32_T c3_r_a;
  int32_T c3_b_jp1j;
  int32_T c3_s_a;
  int32_T c3_t_a;
  int32_T c3_d_c;
  int32_T c3_u_a;
  int32_T c3_m_b;
  int32_T c3_v_a;
  int32_T c3_n_b;
  int32_T c3_i281;
  int32_T c3_w_a;
  int32_T c3_o_b;
  int32_T c3_x_a;
  int32_T c3_p_b;
  boolean_T c3_overflow;
  int32_T c3_i;
  int32_T c3_b_i;
  real_T c3_x;
  real_T c3_y;
  real_T c3_b_x;
  real_T c3_b_y;
  real_T c3_z;
  int32_T c3_q_b;
  int32_T c3_r_b;
  int32_T c3_e_c;
  int32_T c3_y_a;
  int32_T c3_ab_a;
  int32_T c3_f_c;
  int32_T c3_bb_a;
  int32_T c3_cb_a;
  int32_T c3_g_c;
  real_T c3_d2;
  c3_eps(chartInstance);
  for (c3_i277 = 0; c3_i277 < 6; c3_i277++) {
    c3_ipiv[c3_i277] = 1 + c3_i277;
  }

  *c3_info = 0;
  for (c3_j = 1; c3_j < 6; c3_j++) {
    c3_b_j = c3_j;
    c3_a = c3_b_j;
    c3_b_a = c3_a - 1;
    c3_jm1 = c3_b_a;
    c3_b = c3_b_j;
    c3_b_b = c3_b;
    c3_mmj = 6 - c3_b_b;
    c3_c_a = c3_jm1;
    c3_d_a = c3_c_a;
    c3_c = c3_d_a * 7;
    c3_c_b = c3_c;
    c3_d_b = c3_c_b + 1;
    c3_jj = c3_d_b;
    c3_e_a = c3_jj;
    c3_f_a = c3_e_a + 1;
    c3_jp1j = c3_f_a;
    c3_g_a = c3_mmj;
    c3_h_a = c3_g_a;
    c3_b_c = c3_h_a;
    c3_i278 = 0;
    for (c3_i279 = 0; c3_i279 < 6; c3_i279++) {
      for (c3_i280 = 0; c3_i280 < 6; c3_i280++) {
        c3_b_A[c3_i280 + c3_i278] = c3_A[c3_i280 + c3_i278];
      }

      c3_i278 += 6;
    }

    c3_i_a = c3_eml_ixamax(chartInstance, c3_b_c + 1, c3_b_A, c3_jj);
    c3_j_a = c3_i_a - 1;
    c3_jpiv_offset = c3_j_a;
    c3_k_a = c3_jj;
    c3_e_b = c3_jpiv_offset;
    c3_l_a = c3_k_a;
    c3_f_b = c3_e_b;
    c3_jpiv = c3_l_a + c3_f_b;
    if (c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c3_jpiv), 1, 36, 1, 0) - 1] != 0.0) {
      if (c3_jpiv_offset != 0) {
        c3_m_a = c3_b_j;
        c3_g_b = c3_jpiv_offset;
        c3_n_a = c3_m_a;
        c3_h_b = c3_g_b;
        c3_c_c = c3_n_a + c3_h_b;
        c3_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c3_b_j), 1, 6, 1, 0) - 1] = c3_c_c;
        c3_i_b = c3_jm1;
        c3_j_b = c3_i_b + 1;
        c3_jrow = c3_j_b;
        c3_o_a = c3_jrow;
        c3_k_b = c3_jpiv_offset;
        c3_p_a = c3_o_a;
        c3_l_b = c3_k_b;
        c3_jprow = c3_p_a + c3_l_b;
        c3_ix0 = c3_jrow;
        c3_iy0 = c3_jprow;
        c3_b_ix0 = c3_ix0;
        c3_b_iy0 = c3_iy0;
        c3_threshold(chartInstance);
        c3_c_ix0 = c3_b_ix0;
        c3_c_iy0 = c3_b_iy0;
        c3_ix = c3_c_ix0;
        c3_iy = c3_c_iy0;
        for (c3_k = 1; c3_k < 7; c3_k++) {
          c3_temp = c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c3_ix), 1, 36, 1, 0) - 1];
          c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_ix), 1, 36, 1, 0) - 1] = c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK
            ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c3_iy), 1, 36, 1, 0) -
            1];
          c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_iy), 1, 36, 1, 0) - 1] = c3_temp;
          c3_q_a = c3_ix + 6;
          c3_ix = c3_q_a;
          c3_r_a = c3_iy + 6;
          c3_iy = c3_r_a;
        }
      }

      c3_b_jp1j = c3_jp1j;
      c3_s_a = c3_mmj;
      c3_t_a = c3_s_a;
      c3_d_c = c3_t_a;
      c3_u_a = c3_jp1j;
      c3_m_b = c3_d_c - 1;
      c3_v_a = c3_u_a;
      c3_n_b = c3_m_b;
      c3_i281 = c3_v_a + c3_n_b;
      c3_w_a = c3_b_jp1j;
      c3_o_b = c3_i281;
      c3_x_a = c3_w_a;
      c3_p_b = c3_o_b;
      if (c3_x_a > c3_p_b) {
        c3_overflow = false;
      } else {
        c3_eml_switch_helper(chartInstance);
        c3_overflow = (c3_p_b > 2147483646);
      }

      if (c3_overflow) {
        c3_check_forloop_overflow_error(chartInstance, c3_overflow);
      }

      for (c3_i = c3_b_jp1j; c3_i <= c3_i281; c3_i++) {
        c3_b_i = c3_i;
        c3_x = c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c3_b_i), 1, 36, 1, 0) - 1];
        c3_y = c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c3_jj), 1, 36, 1, 0) - 1];
        c3_b_x = c3_x;
        c3_b_y = c3_y;
        c3_z = c3_b_x / c3_b_y;
        c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c3_b_i), 1, 36, 1, 0) - 1] = c3_z;
      }
    } else {
      *c3_info = c3_b_j;
    }

    c3_q_b = c3_b_j;
    c3_r_b = c3_q_b;
    c3_e_c = 6 - c3_r_b;
    c3_y_a = c3_jj;
    c3_ab_a = c3_y_a;
    c3_f_c = c3_ab_a;
    c3_bb_a = c3_jj;
    c3_cb_a = c3_bb_a;
    c3_g_c = c3_cb_a;
    c3_d2 = -1.0;
    c3_b_eml_xgeru(chartInstance, c3_mmj, c3_e_c, c3_d2, c3_jp1j, c3_f_c + 6,
                   c3_A, c3_g_c + 7);
  }

  if (*c3_info == 0) {
    if (!(c3_A[35] != 0.0)) {
      *c3_info = 6;
    }
  }
}

static void c3_b_eml_xgeru(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, int32_T c3_m, int32_T c3_n, real_T c3_alpha1, int32_T c3_ix0,
  int32_T c3_iy0, real_T c3_A[36], int32_T c3_ia0)
{
  int32_T c3_b_m;
  int32_T c3_b_n;
  real_T c3_b_alpha1;
  int32_T c3_b_ix0;
  int32_T c3_b_iy0;
  int32_T c3_b_ia0;
  int32_T c3_c_m;
  int32_T c3_c_n;
  real_T c3_c_alpha1;
  int32_T c3_c_ix0;
  int32_T c3_c_iy0;
  int32_T c3_c_ia0;
  int32_T c3_d_m;
  int32_T c3_d_n;
  real_T c3_d_alpha1;
  int32_T c3_d_ix0;
  int32_T c3_d_iy0;
  int32_T c3_d_ia0;
  int32_T c3_e_m;
  int32_T c3_e_n;
  real_T c3_e_alpha1;
  int32_T c3_e_ix0;
  int32_T c3_e_iy0;
  int32_T c3_e_ia0;
  int32_T c3_ixstart;
  int32_T c3_a;
  int32_T c3_jA;
  int32_T c3_jy;
  int32_T c3_f_n;
  int32_T c3_b;
  int32_T c3_b_b;
  boolean_T c3_overflow;
  int32_T c3_j;
  real_T c3_yjy;
  real_T c3_temp;
  int32_T c3_ix;
  int32_T c3_c_b;
  int32_T c3_i282;
  int32_T c3_b_a;
  int32_T c3_d_b;
  int32_T c3_i283;
  int32_T c3_c_a;
  int32_T c3_e_b;
  int32_T c3_d_a;
  int32_T c3_f_b;
  boolean_T c3_b_overflow;
  int32_T c3_ijA;
  int32_T c3_b_ijA;
  int32_T c3_e_a;
  int32_T c3_f_a;
  int32_T c3_g_a;
  c3_b_m = c3_m;
  c3_b_n = c3_n;
  c3_b_alpha1 = c3_alpha1;
  c3_b_ix0 = c3_ix0;
  c3_b_iy0 = c3_iy0;
  c3_b_ia0 = c3_ia0;
  c3_c_m = c3_b_m;
  c3_c_n = c3_b_n;
  c3_c_alpha1 = c3_b_alpha1;
  c3_c_ix0 = c3_b_ix0;
  c3_c_iy0 = c3_b_iy0;
  c3_c_ia0 = c3_b_ia0;
  c3_d_m = c3_c_m;
  c3_d_n = c3_c_n;
  c3_d_alpha1 = c3_c_alpha1;
  c3_d_ix0 = c3_c_ix0;
  c3_d_iy0 = c3_c_iy0;
  c3_d_ia0 = c3_c_ia0;
  c3_e_m = c3_d_m;
  c3_e_n = c3_d_n;
  c3_e_alpha1 = c3_d_alpha1;
  c3_e_ix0 = c3_d_ix0;
  c3_e_iy0 = c3_d_iy0;
  c3_e_ia0 = c3_d_ia0;
  if (c3_e_alpha1 == 0.0) {
  } else {
    c3_ixstart = c3_e_ix0;
    c3_a = c3_e_ia0 - 1;
    c3_jA = c3_a;
    c3_jy = c3_e_iy0;
    c3_f_n = c3_e_n;
    c3_b = c3_f_n;
    c3_b_b = c3_b;
    if (1 > c3_b_b) {
      c3_overflow = false;
    } else {
      c3_eml_switch_helper(chartInstance);
      c3_overflow = (c3_b_b > 2147483646);
    }

    if (c3_overflow) {
      c3_check_forloop_overflow_error(chartInstance, c3_overflow);
    }

    for (c3_j = 1; c3_j <= c3_f_n; c3_j++) {
      c3_yjy = c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c3_jy), 1, 36, 1, 0) - 1];
      if (c3_yjy != 0.0) {
        c3_temp = c3_yjy * c3_e_alpha1;
        c3_ix = c3_ixstart;
        c3_c_b = c3_jA + 1;
        c3_i282 = c3_c_b;
        c3_b_a = c3_e_m;
        c3_d_b = c3_jA;
        c3_i283 = c3_b_a + c3_d_b;
        c3_c_a = c3_i282;
        c3_e_b = c3_i283;
        c3_d_a = c3_c_a;
        c3_f_b = c3_e_b;
        if (c3_d_a > c3_f_b) {
          c3_b_overflow = false;
        } else {
          c3_eml_switch_helper(chartInstance);
          c3_b_overflow = (c3_f_b > 2147483646);
        }

        if (c3_b_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_b_overflow);
        }

        for (c3_ijA = c3_i282; c3_ijA <= c3_i283; c3_ijA++) {
          c3_b_ijA = c3_ijA;
          c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_b_ijA), 1, 36, 1, 0) - 1] =
            c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_b_ijA), 1, 36, 1, 0) - 1] +
            c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c3_ix), 1, 36, 1, 0) - 1] * c3_temp;
          c3_e_a = c3_ix + 1;
          c3_ix = c3_e_a;
        }
      }

      c3_f_a = c3_jy + 6;
      c3_jy = c3_f_a;
      c3_g_a = c3_jA + 6;
      c3_jA = c3_g_a;
    }
  }
}

static void c3_b_eml_xtrsm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[36], real_T c3_B[36])
{
  int32_T c3_j;
  int32_T c3_b_j;
  int32_T c3_jBcol;
  int32_T c3_k;
  int32_T c3_b_k;
  int32_T c3_kAcol;
  real_T c3_x;
  real_T c3_y;
  real_T c3_b_x;
  real_T c3_b_y;
  real_T c3_c_x;
  real_T c3_c_y;
  real_T c3_z;
  int32_T c3_i284;
  int32_T c3_b;
  int32_T c3_b_b;
  boolean_T c3_overflow;
  int32_T c3_i;
  int32_T c3_b_i;
  c3_b_threshold(chartInstance);
  for (c3_j = 1; c3_j < 7; c3_j++) {
    c3_b_j = c3_j - 1;
    c3_jBcol = 6 * c3_b_j;
    for (c3_k = 6; c3_k > 0; c3_k--) {
      c3_b_k = c3_k;
      c3_kAcol = 6 * (c3_b_k - 1);
      if (c3_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c3_b_k + c3_jBcol)), 1, 36, 1, 0) - 1] != 0.0) {
        c3_x = c3_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)(c3_b_k + c3_jBcol)), 1, 36, 1, 0) - 1];
        c3_y = c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)(c3_b_k + c3_kAcol)), 1, 36, 1, 0) - 1];
        c3_b_x = c3_x;
        c3_b_y = c3_y;
        c3_c_x = c3_b_x;
        c3_c_y = c3_b_y;
        c3_z = c3_c_x / c3_c_y;
        c3_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)(c3_b_k + c3_jBcol)), 1, 36, 1, 0) - 1] = c3_z;
        c3_i284 = c3_b_k - 1;
        c3_b = c3_i284;
        c3_b_b = c3_b;
        if (1 > c3_b_b) {
          c3_overflow = false;
        } else {
          c3_eml_switch_helper(chartInstance);
          c3_overflow = (c3_b_b > 2147483646);
        }

        if (c3_overflow) {
          c3_check_forloop_overflow_error(chartInstance, c3_overflow);
        }

        for (c3_i = 1; c3_i <= c3_i284; c3_i++) {
          c3_b_i = c3_i;
          c3_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c3_b_i + c3_jBcol)), 1, 36, 1, 0) - 1] =
            c3_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c3_b_i + c3_jBcol)), 1, 36, 1, 0) - 1] -
            c3_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c3_b_k + c3_jBcol)), 1, 36, 1, 0) - 1] *
            c3_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c3_b_i + c3_kAcol)), 1, 36, 1, 0) - 1];
        }
      }
    }
  }
}

static void c3_c_eml_xgemm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[180], real_T c3_B[36], real_T c3_C[180])
{
  real_T c3_alpha1;
  real_T c3_beta1;
  char_T c3_TRANSB;
  char_T c3_TRANSA;
  ptrdiff_t c3_m_t;
  ptrdiff_t c3_n_t;
  ptrdiff_t c3_k_t;
  ptrdiff_t c3_lda_t;
  ptrdiff_t c3_ldb_t;
  ptrdiff_t c3_ldc_t;
  double * c3_alpha1_t;
  double * c3_Aia0_t;
  double * c3_Bib0_t;
  double * c3_beta1_t;
  double * c3_Cic0_t;
  c3_c_threshold(chartInstance);
  c3_alpha1 = 1.0;
  c3_beta1 = 0.0;
  c3_TRANSB = 'N';
  c3_TRANSA = 'N';
  c3_m_t = (ptrdiff_t)(30);
  c3_n_t = (ptrdiff_t)(6);
  c3_k_t = (ptrdiff_t)(6);
  c3_lda_t = (ptrdiff_t)(30);
  c3_ldb_t = (ptrdiff_t)(6);
  c3_ldc_t = (ptrdiff_t)(30);
  c3_alpha1_t = (double *)(&c3_alpha1);
  c3_Aia0_t = (double *)(&c3_A[0]);
  c3_Bib0_t = (double *)(&c3_B[0]);
  c3_beta1_t = (double *)(&c3_beta1);
  c3_Cic0_t = (double *)(&c3_C[0]);
  dgemm(&c3_TRANSA, &c3_TRANSB, &c3_m_t, &c3_n_t, &c3_k_t, c3_alpha1_t,
        c3_Aia0_t, &c3_lda_t, c3_Bib0_t, &c3_ldb_t, c3_beta1_t, c3_Cic0_t,
        &c3_ldc_t);
}

static void c3_d_eml_xgemm(SFc3_aircraftControl_FullStateFiltersInstanceStruct
  *chartInstance, real_T c3_A[180], real_T c3_B[180], real_T c3_C[900])
{
  real_T c3_alpha1;
  real_T c3_beta1;
  char_T c3_TRANSB;
  char_T c3_TRANSA;
  ptrdiff_t c3_m_t;
  ptrdiff_t c3_n_t;
  ptrdiff_t c3_k_t;
  ptrdiff_t c3_lda_t;
  ptrdiff_t c3_ldb_t;
  ptrdiff_t c3_ldc_t;
  double * c3_alpha1_t;
  double * c3_Aia0_t;
  double * c3_Bib0_t;
  double * c3_beta1_t;
  double * c3_Cic0_t;
  c3_c_threshold(chartInstance);
  c3_alpha1 = 1.0;
  c3_beta1 = 0.0;
  c3_TRANSB = 'N';
  c3_TRANSA = 'N';
  c3_m_t = (ptrdiff_t)(30);
  c3_n_t = (ptrdiff_t)(30);
  c3_k_t = (ptrdiff_t)(6);
  c3_lda_t = (ptrdiff_t)(30);
  c3_ldb_t = (ptrdiff_t)(6);
  c3_ldc_t = (ptrdiff_t)(30);
  c3_alpha1_t = (double *)(&c3_alpha1);
  c3_Aia0_t = (double *)(&c3_A[0]);
  c3_Bib0_t = (double *)(&c3_B[0]);
  c3_beta1_t = (double *)(&c3_beta1);
  c3_Cic0_t = (double *)(&c3_C[0]);
  dgemm(&c3_TRANSA, &c3_TRANSB, &c3_m_t, &c3_n_t, &c3_k_t, c3_alpha1_t,
        c3_Aia0_t, &c3_lda_t, c3_Bib0_t, &c3_ldb_t, c3_beta1_t, c3_Cic0_t,
        &c3_ldc_t);
}

static void init_dsm_address_info
  (SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance)
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

void sf_c3_aircraftControl_FullStateFilters_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1635658176U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2647582931U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(239511733U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3467139067U);
}

mxArray *sf_c3_aircraftControl_FullStateFilters_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("706fMF9jq2bMlpieGmPLYG");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,7,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,6,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,6,"type",mxType);
    }

    mxSetField(mxData,6,"complexity",mxCreateDoubleScalar(0));
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
      pr[0] = (double)(30);
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
      pr[0] = (double)(30);
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
      pr[0] = (double)(6);
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
      pr[0] = (double)(6);
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

mxArray *sf_c3_aircraftControl_FullStateFilters_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c3_aircraftControl_FullStateFilters_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c3_aircraftControl_FullStateFilters
  (void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x9'type','srcId','name','auxInfo'{{M[1],M[14],T\"S\",},{M[1],M[10],T\"cov_OUT\",},{M[1],M[13],T\"innovation\",},{M[1],M[5],T\"states_OUT\",},{M[4],M[0],T\"covDiag\",S'l','i','p'{{M1x2[174 181],M[1],T\"C:\\Users\\Rudaba\\Documents\\PhD Take 2\\Simulations\\AircraftModelling\\simulateUKF.m\"}}},{M[4],M[0],T\"pCovariance\",S'l','i','p'{{M1x2[150 161],M[1],T\"C:\\Users\\Rudaba\\Documents\\PhD Take 2\\Simulations\\AircraftModelling\\simulateUKF.m\"}}},{M[4],M[0],T\"pStates\",S'l','i','p'{{M1x2[130 137],M[1],T\"C:\\Users\\Rudaba\\Documents\\PhD Take 2\\Simulations\\AircraftModelling\\simulateUKF.m\"}}},{M[4],M[0],T\"prevOmegaDots\",S'l','i','p'{{M1x2[194 207],M[1],T\"C:\\Users\\Rudaba\\Documents\\PhD Take 2\\Simulations\\AircraftModelling\\simulateUKF.m\"}}},{M[8],M[0],T\"is_active_c3_aircraftControl_FullStateFilters\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 9, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c3_aircraftControl_FullStateFilters_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _aircraftControl_FullStateFiltersMachineNumber_,
           3,
           1,
           1,
           0,
           11,
           0,
           0,
           0,
           0,
           1,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation
          (_aircraftControl_FullStateFiltersMachineNumber_,
           chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,
             _aircraftControl_FullStateFiltersMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _aircraftControl_FullStateFiltersMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"accelerations");
          _SFD_SET_DATA_PROPS(1,1,1,0,"omega");
          _SFD_SET_DATA_PROPS(2,1,1,0,"ALP");
          _SFD_SET_DATA_PROPS(3,1,1,0,"BET");
          _SFD_SET_DATA_PROPS(4,1,1,0,"u");
          _SFD_SET_DATA_PROPS(5,2,0,1,"states_OUT");
          _SFD_SET_DATA_PROPS(6,1,1,0,"Vt");
          _SFD_SET_DATA_PROPS(7,2,0,1,"cov_OUT");
          _SFD_SET_DATA_PROPS(8,1,1,0,"H");
          _SFD_SET_DATA_PROPS(9,2,0,1,"innovation");
          _SFD_SET_DATA_PROPS(10,2,0,1,"S");
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
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,198);
        _SFD_CV_INIT_SCRIPT(0,1,2,0,0,0,6,0,2,1);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"simulateUKF",0,-1,11123);
        _SFD_CV_INIT_SCRIPT_IF(0,0,210,229,-1,2492);
        _SFD_CV_INIT_SCRIPT_IF(0,1,4517,4537,4593,4659);
        _SFD_CV_INIT_SCRIPT_FOR(0,0,4484,4512,4712);
        _SFD_CV_INIT_SCRIPT_FOR(0,1,4748,4765,10238);
        _SFD_CV_INIT_SCRIPT_FOR(0,2,10362,10379,10423);
        _SFD_CV_INIT_SCRIPT_FOR(0,3,10425,10442,10526);
        _SFD_CV_INIT_SCRIPT_FOR(0,4,10563,10580,10625);
        _SFD_CV_INIT_SCRIPT_FOR(0,5,10673,10690,10832);

        {
          static int condStart[] = { 4520, 4530 };

          static int condEnd[] = { 4526, 4537 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_SCRIPT_MCDC(0,0,4520,4537,2,0,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_c_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 30;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)
            c3_b_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 30;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_b_sf_marshallOut,(MexInFcnForType)
            c3_b_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        {
          real_T *c3_ALP;
          real_T *c3_BET;
          real_T *c3_Vt;
          real_T *c3_H;
          real_T (*c3_accelerations)[3];
          real_T (*c3_omega)[3];
          real_T (*c3_u)[4];
          real_T (*c3_states_OUT)[30];
          real_T (*c3_cov_OUT)[30];
          real_T (*c3_innovation)[6];
          real_T (*c3_S)[6];
          c3_S = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 4);
          c3_innovation = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S,
            3);
          c3_H = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
          c3_cov_OUT = (real_T (*)[30])ssGetOutputPortSignal(chartInstance->S, 2);
          c3_Vt = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
          c3_states_OUT = (real_T (*)[30])ssGetOutputPortSignal(chartInstance->S,
            1);
          c3_u = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 4);
          c3_BET = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c3_ALP = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c3_omega = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
          c3_accelerations = (real_T (*)[3])ssGetInputPortSignal
            (chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c3_accelerations);
          _SFD_SET_DATA_VALUE_PTR(1U, *c3_omega);
          _SFD_SET_DATA_VALUE_PTR(2U, c3_ALP);
          _SFD_SET_DATA_VALUE_PTR(3U, c3_BET);
          _SFD_SET_DATA_VALUE_PTR(4U, *c3_u);
          _SFD_SET_DATA_VALUE_PTR(5U, *c3_states_OUT);
          _SFD_SET_DATA_VALUE_PTR(6U, c3_Vt);
          _SFD_SET_DATA_VALUE_PTR(7U, *c3_cov_OUT);
          _SFD_SET_DATA_VALUE_PTR(8U, c3_H);
          _SFD_SET_DATA_VALUE_PTR(9U, *c3_innovation);
          _SFD_SET_DATA_VALUE_PTR(10U, *c3_S);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _aircraftControl_FullStateFiltersMachineNumber_,
        chartInstance->chartNumber,chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "2wHByxgJNPJRlDmCjplk0G";
}

static void sf_opaque_initialize_c3_aircraftControl_FullStateFilters(void
  *chartInstanceVar)
{
  chart_debug_initialization
    (((SFc3_aircraftControl_FullStateFiltersInstanceStruct*) chartInstanceVar)
     ->S,0);
  initialize_params_c3_aircraftControl_FullStateFilters
    ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*) chartInstanceVar);
  initialize_c3_aircraftControl_FullStateFilters
    ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c3_aircraftControl_FullStateFilters(void
  *chartInstanceVar)
{
  enable_c3_aircraftControl_FullStateFilters
    ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c3_aircraftControl_FullStateFilters(void
  *chartInstanceVar)
{
  disable_c3_aircraftControl_FullStateFilters
    ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c3_aircraftControl_FullStateFilters(void
  *chartInstanceVar)
{
  sf_gateway_c3_aircraftControl_FullStateFilters
    ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*) chartInstanceVar);
}

extern const mxArray*
  sf_internal_get_sim_state_c3_aircraftControl_FullStateFilters(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c3_aircraftControl_FullStateFilters
    ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*)
     chartInfo->chartInstance);        /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c3_aircraftControl_FullStateFilters
    ();                                /* state var info */
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

extern void sf_internal_set_sim_state_c3_aircraftControl_FullStateFilters
  (SimStruct* S, const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c3_aircraftControl_FullStateFilters
    ();                                /* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c3_aircraftControl_FullStateFilters
    ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*)
     chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray*
  sf_opaque_get_sim_state_c3_aircraftControl_FullStateFilters(SimStruct* S)
{
  return sf_internal_get_sim_state_c3_aircraftControl_FullStateFilters(S);
}

static void sf_opaque_set_sim_state_c3_aircraftControl_FullStateFilters
  (SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c3_aircraftControl_FullStateFilters(S, st);
}

static void sf_opaque_terminate_c3_aircraftControl_FullStateFilters(void
  *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*)
                    chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_aircraftControl_FullStateFilters_optimization_info();
    }

    finalize_c3_aircraftControl_FullStateFilters
      ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc3_aircraftControl_FullStateFilters
    ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c3_aircraftControl_FullStateFilters(SimStruct
  *S)
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
    initialize_params_c3_aircraftControl_FullStateFilters
      ((SFc3_aircraftControl_FullStateFiltersInstanceStruct*)
       (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c3_aircraftControl_FullStateFilters(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct =
      load_aircraftControl_FullStateFilters_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,3);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,3,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,3,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,3);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,3,7);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,3,4);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=4; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 7; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,3);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2137117815U));
  ssSetChecksum1(S,(3553118480U));
  ssSetChecksum2(S,(1816332075U));
  ssSetChecksum3(S,(1558164765U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c3_aircraftControl_FullStateFilters(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c3_aircraftControl_FullStateFilters(SimStruct *S)
{
  SFc3_aircraftControl_FullStateFiltersInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc3_aircraftControl_FullStateFiltersInstanceStruct *)
    utMalloc(sizeof(SFc3_aircraftControl_FullStateFiltersInstanceStruct));
  memset(chartInstance, 0, sizeof
         (SFc3_aircraftControl_FullStateFiltersInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.mdlStart =
    mdlStart_c3_aircraftControl_FullStateFilters;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c3_aircraftControl_FullStateFilters;
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

void c3_aircraftControl_FullStateFilters_method_dispatcher(SimStruct *S, int_T
  method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c3_aircraftControl_FullStateFilters(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c3_aircraftControl_FullStateFilters(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c3_aircraftControl_FullStateFilters(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c3_aircraftControl_FullStateFilters_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
