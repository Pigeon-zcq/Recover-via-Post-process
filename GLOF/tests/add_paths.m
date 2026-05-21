%ADD_PATHS  Add the MATLAB experiment folder and tests folder to the path.

this_file = mfilename('fullpath');
tests_dir = fileparts(this_file);
matlab_dir = fileparts(tests_dir);

addpath(matlab_dir);
addpath(tests_dir);
