function [lexqr_struct] = lexqr_lambda(lexqr_struct)
  %%%
  %
  % calls lexqr_lambda_obj.m to form the matrix of Lagrange multipliers
  %

  %% -------------------------------------
  %% input handling
  %% -------------------------------------
  nObj = lexqr_struct.nObj;
  nCtr = lexqr_struct.nCtr;
  nVar = lexqr_struct.nVar;
  dim  = [lexqr_struct.obj_info(:).dim];

  %% -------------------------------------
  %% determine the right-hand-side
  %% -------------------------------------

  if isfield(lexqr_struct, 'X_mu')
    RT  = lexqr_struct.Rr(:,1:end-1)';
    P   = lexqr_struct.P;
    iRT = pinv(RT);
  else
    rhs_all = zeros(nVar,1);
  end

  %% -------------------------------------

  L = zeros(nCtr,nObj);
  for k = 1:nObj
    ind = 1:sum(dim(1:k));

    if isfield(lexqr_struct, 'X_mu')
      mu = lexqr_struct.options.regularization_factors(k);

      rhs_all = -mu^2*iRT*P'*lexqr_struct.Xd(:,k);
      %%tmp = -mu^2*iRT*P'*lexqr_struct.X_mu(:,k);

      %{
      n = size(RT,2);
      y = -mu^2*P'*lexqr_struct.Xd(:,k);
      rr = RT(1:n,1:n)\y(1:n);
      keyboard
      %}

      L(ind,k) = lexqr_lambda_obj(lexqr_struct,k,rhs_all,1);
    else
      L(ind,k) = lexqr_lambda_obj(lexqr_struct,k,rhs_all,0);
    end
  end

  lexqr_struct.lambda = L;

  %%%EOF
