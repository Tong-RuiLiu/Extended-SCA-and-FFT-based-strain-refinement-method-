% Source code of 3-stage extended full-field self-consistent clustering analysis
% for 3D anisotrpic woven composites 
% The code is distributed under BSD 3-Clause License
% Copyright (c) 2023, Tong-Rui Liu
% email: tongrui.liu18@imperial.ac.uk,tongruiliu1996@outlook.com    
% Imperial College London

% If using this code repository for research (Only!), please cite:
% Tong-Rui Liu, Yang Yang, Omar R. Bacarreza, Shaoqiang Tang and M.H. Aliabadi
% An extended full field self-consistent cluster analysis framework for woven composite
% International Journal of Solids and Structures 281: 112407 (2023)
% https://doi.org/10.1016/j.ijsolstr.2023.112407
%%  This subroutine is used for setpath 
clc;
clear variables;
close all
warning off
home=pwd;
addpath(genpath(home)); 