%
% Copyright 2013-2021 INRIA
%

function options = test_lexlse_define()
%
% defines tests
%

load_labels

i = 0;

% =========================================================================================
% =========================================================================================
% REGULARIZATION_TIKHONOV
% =========================================================================================
% =========================================================================================

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV;
options(i).regularization_factors  = [1,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 3;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 3;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV;
options(i).regularization_factors  = [0,2,3,4];

% =========================================================================================
% =========================================================================================
% REGULARIZATION_TIKHONOV_1
% =========================================================================================
% =========================================================================================

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_1;
options(i).regularization_factors  = [1,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_1;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_1;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_1;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_1;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_1;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 3;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_1;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 3;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_1;
options(i).regularization_factors  = [0,2,3,4];

% =========================================================================================
% =========================================================================================
% REGULARIZATION_TIKHONOV_2
% =========================================================================================
% =========================================================================================

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_2;
options(i).regularization_factors  = [1,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_2;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_2;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_2;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_2;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_2;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 3;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_2;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 3;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_TIKHONOV_2;
options(i).regularization_factors  = [0,2,3,4];


% =========================================================================================
% =========================================================================================
% REGULARIZATION_R
% =========================================================================================
% =========================================================================================


% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_R;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_R;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_R;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_R;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_R;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_R;
options(i).regularization_factors  = [0,2,3,4];

% =========================================================================================
% =========================================================================================
% REGULARIZATION_R_NO_Z
% =========================================================================================
% =========================================================================================

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_R_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_R_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_R_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_R_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_R_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_R_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% =========================================================================================
% =========================================================================================
% REGULARIZATION_RT_NO_Z
% =========================================================================================
% =========================================================================================

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_RT_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 0;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_RT_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_RT_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 1;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_RT_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 1;
options(i).regularization_type     = REGULARIZATION_RT_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

% ------------------------------------------------------------------
i = i + 1;
% ------------------------------------------------------------------
options(i).get_least_norm_solution = 2;
options(i).enable_fixed_variables  = 0;
options(i).regularization_type     = REGULARIZATION_RT_NO_Z;
options(i).regularization_factors  = [0,2,3,4];

%%%EOF
