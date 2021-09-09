%
% Copyright 2013-2021 INRIA
%

function wset = wset_modify_status(wset, obj_index, ctr_index, ctr_type)
%%%
%
% obj_index and ctr_index start from zero
%

    wset{obj_index+1}(ctr_index+1) = ctr_type;

%%%EOF
