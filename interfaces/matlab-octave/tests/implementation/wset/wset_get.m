%
% Copyright 2013-2021 INRIA
%

function wset = wset_get(wset, wlog, iter_num)
%%%
%
% get the working set at iteration number iter_num
%

    wset = wset_deactivate_all_ineq(wset);

    for i=1:iter_num
	wset = wset_modify_status(wset, wlog(i).obj_index, wlog(i).ctr_index, wlog(i).ctr_type);
    end

%%%EOF
