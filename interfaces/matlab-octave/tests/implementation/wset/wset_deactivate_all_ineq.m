function wset = wset_deactivate_all_ineq(wset)
  %%%
  %
  %
  %

  nObj = length(wset);

  for i=1:nObj
    nCtr = length(wset{i});
    for j=1:nCtr
      if (wset{i}(j) == 1) || (wset{i}(j) == 2)
        wset{i}(j) == 0;
      end
    end
  end

  %%%EOF
