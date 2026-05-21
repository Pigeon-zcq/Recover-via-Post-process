%RUN_REQUESTED_TESTS  Run the requested Fourier and GLOF tests.

clear; clc; close all;
add_paths;

fprintf('Running requested MATLAB tests in %s\n', fileparts(mfilename('fullpath')));

test_fourier_convergence;
test_glof_projection_convergence;
test_glof_from_fourier_recovery;
test_glof_from_fourier_sin_multiplier;

fprintf('\nAll requested tests completed.\n');
