function vec = lexqr_residual_compact(lexqr_struct)
  %%%
  %
  % form the residuals from the compact decomposition (naive implementation)
  %

  %% -------------------------------------
  %% input handling
  %% -------------------------------------
  lexqr = lexqr_struct.lexqr;
  nObj  = lexqr_struct.nObj;
  r0    = [lexqr_struct.obj_info(:).rank];
  m0    = [lexqr_struct.obj_info(:).dim];
  Qr    = lexqr_struct.Qr;
  Qs    = lexqr_struct.Qs;

  Ab    = obj2array(lexqr_struct.obj);
  rhs0  = Ab(:,end);

  %% -------------------------------------
  %% index handling
  %% -------------------------------------
  dim = 0; % dimension
  rnk = 0; % rank
  sng = 0; % singular indexes

  for i=1:nObj
    ind(i).d = dim+1:dim+m0(i);
    ind(i).r = rnk+1:rnk+r0(i);
    ind(i).s = sng+1:sng+m0(i)-r0(i);

    dim = dim + m0(i);
    rnk = rnk + r0(i);
    sng = sng + m0(i) - r0(i);
  end

  %% -------------------------------------
  %% compute the residuals
  %% -------------------------------------

  for i = 1:nObj

    d = ind(i).d;
    r = ind(i).r;
    s = ind(i).s;

    rhs = rhs0(d);
    for j = 1:i-1
      rj = ind(j).r;

      rhs = rhs - Qr(d,rj)*y{j}(1:length(rj));
    end

    Q = [Qr(d,r), Qs(d,s)];

    y{i} = Q'*rhs;
    v{i} = -Qs(d,s)*y{i}(length(r)+1:end);
  end

  vec = [];
  for i=1:nObj
    vec = [vec;v{i}];
  end

  %%%EOF
