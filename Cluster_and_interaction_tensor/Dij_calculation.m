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
%%
% This subroutime concerns about computing the interaction tensor in
% offline computing stage, using the subroutine ParallelcomputeDij.m
clc
clear variables
clc
% size of domain/UC parameter length/width/height
m=128;n=128;l=32; % Number of discritizations for voxel 
Lx =128;Ly=128;Lz =32;% size of RVE
% Number of clusters
load('ClusterData-64-16-4.mat')
%% Predefine the Number of the pattern (Defined by the users)
Ncluster_m=64;
Ncluster_mU=16;
Ncluster_mA=4;
TotNc=Ncluster_m+4*(Ncluster_mA*Ncluster_mU);% Total number of clusters

list = ClusterData(:,1);
idx  = ClusterData(:,8);

%% Continious Green's operator based on Moulinec-Suquet discretization 
% tic
% [D1,D2] = ParallelcomputeDij(TotNc,m,n,l,Lx,Ly,Lz,list,idx,1);
% toc
% save('D1Suquet.mat','D1')
% save('D2Suquet.mat','D2')
% clear D1 D2
% clear Ncluster Ncluster_m Ncluster_f list idx
%% Discrete Green's operator based on Willot's discretization 
tic
[D1,D2] = ParallelcomputeDij(TotNc,m,n,l,Lx,Ly,Lz,list,idx,2);
toc
save('D1Willot.mat','D1')
save('D2Willot.mat','D2')
clear D1 D2
clear Ncluster Ncluster_m Ncluster_f list idx

