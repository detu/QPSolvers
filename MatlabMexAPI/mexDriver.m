%Tested with Matlab R2019b and Visual Studio 2015
clc; clear all
inScalar =  153;
inVector =  [1 2 3];
inMatrix = [1 2 3; 11 22 33; 111 222 333; 1111 2222 3333];
inCharArray = 'c:/temp/a/';
inString = "c:/temp/a/"; %Note double quotes
inStruct.x = [1 2 3];
inStruct.y = [11 22; 111 222];
mex mexWorkspaceDemo.cpp
mexWorkspaceDemo();