%
% Copyright 2013-2021 INRIA
%

FILE{1}       = '@CMAKE_CURRENT_SOURCE_DIR@/lexlse.cpp';
FILE{2}       = '@CMAKE_CURRENT_SOURCE_DIR@/lexlsi.cpp';

function retval = flatten_cmake_list(l, joinc = ' ')
  retval = strsplit(l, ';');
  retval = retval(~cellfun('isempty',retval));
  retval = strjoin(retval, joinc);
endfunction

LEXLS_COMPILE_DEFINITIONS = flatten_cmake_list('@LEXLS_COMPILE_DEFINITIONS@');
LEXLS_COMPILE_FLAGS = flatten_cmake_list('@LEXLS_COMPILE_FLAGS@');
LEXLS_LINK_FLAGS = flatten_cmake_list('@LEXLS_LINK_FLAGS@');
LEXLS_INCLUDE_DIRECTORIES = flatten_cmake_list('@LEXLS_INCLUDE_DIRECTORIES@', ' -I');
if length(LEXLS_INCLUDE_DIRECTORIES) ~= 0
  LEXLS_INCLUDE_DIRECTORIES = ['-I', LEXLS_INCLUDE_DIRECTORIES];
end

CMAKE_CXX_FLAGS = {
  '@CMAKE_CXX_FLAGS@',
  '$<$<CONFIG:DEBUG>:@CMAKE_CXX_FLAGS_DEBUG@>'
  '$<$<CONFIG:MINSIZEREL>:@CMAKE_CXX_FLAGS_MINSIZEREL@>'
  '$<$<CONFIG:RELWITHDEBINFO>:@CMAKE_CXX_FLAGS_RELWITHDEBINFO@>'
  '$<$<CONFIG:RELEASE>:@CMAKE_CXX_FLAGS_RELEASE@>'
};
CMAKE_CXX_FLAGS = CMAKE_CXX_FLAGS(~cellfun('isempty',CMAKE_CXX_FLAGS));
CMAKE_CXX_FLAGS = strjoin(CMAKE_CXX_FLAGS, ' ');

setenv('CXXFLAGS', [LEXLS_COMPILE_DEFINITIONS, ' ', LEXLS_COMPILE_FLAGS', ' ', '@CMAKE_CXX_FLAGS@ $<$<CONFIG:RelWithDebInfo>:@CMAKE_CXX_FLAGS_RELWITHDEBINFO@>'])

% ------------------------------------------------------------
% compilation
% ------------------------------------------------------------
for i=1:length(FILE)
  cc = ['mex -v', ' ', LEXLS_INCLUDE_DIRECTORIES, ' ', FILE{i}];
  disp(cc);
  ret = eval(cc);
  if ret ~= 0
    exit(ret);
  end
end

%%%EOF
