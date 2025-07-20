addpath("Function");

processDate = "2025-01-01";
dir_ephemeris_txt = "Ephemeris";
dir_ephemeris_mat = "Ephemeris_mat";
dir_deltaE_mat = "DeltaE_mat";
% Step 1
% Compile the necessary sofa_c code
Make_src;
% Step 2
% Transform the Starlink Ephemeris from txt to mat data file.
Ephemeris_txt2mat(processDate, dir_ephemeris_txt, dir_ephemeris_mat);
% Step 3
% Calculate the Mechanical Energy Loss from Starlink Ephemeris.
Ephemeris2DeltaE(processDate, dir_ephemeris_mat, dir_deltaE_mat);


