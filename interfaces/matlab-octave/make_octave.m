%
% Copyright 2013-2021 INRIA
%

% on MACOS
% copy mexopts.sh from ~/.matlab/R2011b to lexls/interface and set MACOSX_DEPLOYMENT_TARGET='10.7'
%
% http://stackoverflow.com/questions/15568925/mex-files-on-mountain-lion-explicit-instantiation-error


% ------------------------------------------------------------
% compilation settings
% ------------------------------------------------------------
EIGEN_INCLUDE = '/usr/local/include/eigen3';
LEXLS_INCLUDE = '../../include';

FILE{1}       = 'lexlse.cpp';
FILE{2}       = 'lexlsi.cpp';

%setenv('CXXFLAGS', cstrcat(octave_config_info.CXXFLAGS, ' -pedantic -Wall -std=c++98 '))
%setenv('CXXFLAGS', ' -g -pedantic -Wall -std=c++98 -O0 -fno-strict-aliasing  -D_THREAD_SAFE -pthread ')
%setenv('DL_LDFLAGS', ' -g -shared -Wl,-x ');
%setenv('CXXFLAGS', cstrcat(octave_config_info.CXXFLAGS, ' -pedantic -Wall -O3 -DNDEBUG '))
setenv('CXXFLAGS', '-march=native -fstack-protector -fno-strict-aliasing -D_THREAD_SAFE -pthread -pedantic -Wall -O3 -DNDEBUG ')
% ------------------------------------------------------------
% compilation
% ------------------------------------------------------------
INCLUDES = ['-I', EIGEN_INCLUDE, ' -I', LEXLS_INCLUDE];

for i=1:length(FILE)
    %%cc = ['mex -v -g ', ' ', INCLUDES, ' ', FILE{i}];
    cc = ['mex -v', ' ', INCLUDES, ' ', FILE{i}];

  disp(cc);

  eval(cc);
end

%%%EOF
