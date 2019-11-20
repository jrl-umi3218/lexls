function L = lexqr_lambda_compact(lexqr_struct)
  %%%
  %
  % form the Lagrange multipliers from the compact decomposition (naive implementation)
  % (each objective is assumed to be minimize 0.5*r'*r)
  %
  % requires the residuals lexqr_struct.v
  %

  %% -------------------------------------
  %% input handling
  %% -------------------------------------
  lexqr = lexqr_struct.lexqr;
  nObj  = lexqr_struct.nObj;
  nCtr  = lexqr_struct.nCtr;
  r0    = [lexqr_struct.obj_info(:).rank];
  m0    = [lexqr_struct.obj_info(:).dim];
  Qr    = lexqr_struct.Qr;
  v     = lexqr_struct.v;

  %% -------------------------------------
  %% initializations
  %% -------------------------------------

  QrT   = Qr';
  L     = zeros(nCtr,nObj);

  %% -------------------------------------
  %% index handling
  %% -------------------------------------

  dim = 0; % dimension
  rnk = 0; % rank

  for i=1:nObj
    ind(i).d = dim+1:dim+m0(i);
    ind(i).r = rnk+1:rnk+r0(i);

    dim = dim + m0(i);
    rnk = rnk + r0(i);
  end

  %% -------------------------------------
  %% compute the multipliers
  %% -------------------------------------

  for i = nObj:-1:1 % sensitivity of the i-th objective

    rhs = zeros(sum(r0(1:i)),1);

    for j = i:-1:1 % constraints in the j-th objective

      d = ind(j).d;
      r = ind(j).r;
      I = 1:sum(r0(1:i-1));

      if j == i
        L(d,i) = v(d);
      else
        L(d,i) = QrT(r,d)'*rhs(r);
      end

      rhs(I) = rhs(I) - QrT(I,d)*L(d,i);
    end

  end

  %%%EOF
