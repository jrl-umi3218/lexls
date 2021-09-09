%
% Copyright 2013-2021 INRIA
%

function lexqr_struct = options_handler(lexqr_struct)
%%%
%
%
%

    %% -----------------------------------------------------------------------------------
    %% dafault options
    %% -----------------------------------------------------------------------------------

    def.tol_linear_dependence  = 1e-12;
    def.regularization_type    = 0;
    def.regularization_factors = [];

    out = def;

    %% -----------------------------------------------------------------------------------
    %% modify the dafault options
    %% -----------------------------------------------------------------------------------

    if isfield(lexqr_struct,'options')
	options = lexqr_struct.options;

	names = fieldnames(options);
	for i = 1:length(names)

	    value = getfield(options, char(names(i)));

	    out = setfield(out, char(names(i)), value);
	    if ~isfield(def,char(names(i)))
		fprintf('WARNING: set non-default option %s \n', char(names(i)));
	    end
	end
    end

    lexqr_struct.options = out;

%%%EOF
