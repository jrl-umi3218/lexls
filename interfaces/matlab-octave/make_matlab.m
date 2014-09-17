% on MACOS
% copy mexopts.sh from ~/.matlab/R2011b to lexls/interface and set MACOSX_DEPLOYMENT_TARGET='10.7' 
%
% http://stackoverflow.com/questions/15568925/mex-files-on-mountain-lion-explicit-instantiation-error

system('rm -rf *.mexmaci*');

% ------------------------------------------------------------
% compilation settings
% ------------------------------------------------------------
EIGEN_INCLUDE = '/usr/local/include/eigen3';
LEXLS_INCLUDE = '../../include';

FILE{1}       = 'lexlse.cpp';
FILE{2}       = 'lexlsi.cpp';

% ------------------------------------------------------------
% compilation 
% ------------------------------------------------------------
INCLUDES = ['-I', EIGEN_INCLUDE, ' -I', LEXLS_INCLUDE];

for i=1:length(FILE)
  cc = ['mex ', ' -v ', INCLUDES, ' ', FILE{i}];

  disp(cc);

  eval(cc);
end

%%%EOF
