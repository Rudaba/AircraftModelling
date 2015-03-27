#ifndef __c1_aircraftControl_smallerStateEKF_h__
#define __c1_aircraftControl_smallerStateEKF_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc1_aircraftControl_smallerStateEKFInstanceStruct
#define typedef_SFc1_aircraftControl_smallerStateEKFInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_aircraftControl_smallerStateEKF;
  real_T c1_states[19];
  boolean_T c1_states_not_empty;
  real_T c1_cov[361];
  boolean_T c1_cov_not_empty;
  real_T c1_prev_p;
  boolean_T c1_prev_p_not_empty;
  real_T c1_prev_q;
  boolean_T c1_prev_q_not_empty;
  real_T c1_prev_r;
  boolean_T c1_prev_r_not_empty;
} SFc1_aircraftControl_smallerStateEKFInstanceStruct;

#endif                                 /*typedef_SFc1_aircraftControl_smallerStateEKFInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray
  *sf_c1_aircraftControl_smallerStateEKF_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_aircraftControl_smallerStateEKF_get_check_sum(mxArray *plhs[]);
extern void c1_aircraftControl_smallerStateEKF_method_dispatcher(SimStruct *S,
  int_T method, void *data);

#endif
