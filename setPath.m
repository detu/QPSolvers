%SETPATH add the subfolder FEM_library and RB_library to the path 

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath(genpath(strcat(pwd,sslash,'LinearSolver')));
addpath(genpath(strcat(pwd,sslash,'FilterQP')));
addpath(genpath(strcat(pwd,sslash,'QPA')));
addpath(genpath(strcat(pwd,sslash,'QPB')));
addpath(genpath(strcat(pwd,sslash,'QPC')));
rmpath(genpath(strcat(pwd,sslash,'LinearSolver/Mumps/Libraries')));

