#ifndef __c1_aircraftControl_FullEKF_h__
#define __c1_aircraftControl_FullEKF_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc1_aircraftControl_FullEKFInstanceStruct
#define typedef_SFc1_aircraftControl_FullEKFInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_aircraftControl_FullEKF;
  real_T c1_pStates[30];
  boolean_T c1_pStates_not_empty;
  real_T c1_pCovariance[900];
  boolean_T c1_pCovariance_not_empty;
  real_T c1_covDiag[30];
  boolean_T c1_covDiag_not_empty;
  real_T c1_prevOmegaDots[3];
  boolean_T c1_prevOmegaDots_not_empty;
} SFc1_aircraftControl_FullEKFInstanceStruct;

#endif                                 /*typedef_SFc1_aircraftControl_FullEKFInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray
  *sf_c1_aircraftControl_FullEKF_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_aircraftControl_FullEKF_get_check_sum(mxArray *plhs[]);
extern void c1_aircraftControl_FullEKF_method_dispatcher(SimStruct *S, int_T
  method, void *data);

#endif
