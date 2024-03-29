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
%% This is the strain refinement stage in 3 stage full-field self-consistent clustering analysis
clc
clear variables
close all
disp('input the material parameter')
%%  initialize the the cluster number, strain and stress 
% Load cluster ID 
load('ClusterData-64-16-4.mat')
% % create the stress tensor to store the stress vector of all elements
load('SCAstran.mat')
load('SCAstress.mat')
Npm = 64;
Npy = 16;
Npmeso = Npm+4*4*Npy;
%%  Material Elasticity parameter (Yarn in local coordinate)
E11=5.312e+10; % out of plane Young's modulus
E22=1.4660e+10; % in plane Young's modulus
E33=E22;
v12=0.266;v13=v12;v23=0.268;
G12=4.240e+09;G13=G12;G23=5.780e+09;
%Material Elasticity parameter (Matrix)
Em=3.170e+9;num=0.35;
%  input your material property in each voxel element size of the mesh
disp('input the mesh size')
m=128;n=128;l=32;
Lx =128;Ly = 128;Lz = 32;
cons=zeros(m*n*l,36); % Define the constitutive law of each voxel in UC
%% initialize the Raw material constitutive
disp('initialize the Raw material constitutive') 
% read the Undulation angle file
load('ORI.mat')
% read the YarnX ID
load('YarnX.mat')
% read the YarnY ID 
load('YarnY.mat')
% another are matrix ID 
load ('Matrix.mat')
% update the constitutive matrix (put the constitutive matrix in each row)
disp('update the constitutive matrix')
for p = 1: length(Matrix)
     cons(Matrix(p),:)=reshape(Elasm3D(Em,num),[1,36]);
end
for q = 1: length(YarnX)
     cons(YarnX(q),:)=reshape(Elasm3DYarnX(E11,E22,E33,v12,v13,v23,G12,G13,G23,ORI(YarnX(q),1),ORI(YarnX(q),3)),[1,36]);
end
for r = 1: length(YarnY)
     cons(YarnY(r),:)=reshape(Elasm3DYarnY(E11,E22,E33,v12,v13,v23,G12,G13,G23,ORI(YarnY(r),2),ORI(YarnY(r),3)),[1,36]);
end
disp('Calculating the constitutive matrix end')
%% Calculate the eigenvalue of anisotropic tensor
disp('input the reference material')
% Here we set the Green's Function Lamada = 0 and G0 dependent on eigenvalue of stiffness matrix 
lamada0=0;
Evlist=zeros(m*n*l,6);
for i = 1: m*n*l 
    Evlist(i,:) = (eig(reshape(cons(i,:),[6,6])))';
end 
Evmin = min(min(Evlist));
Evmax = max(max(Evlist));
G0=(Evmin+Evmax)/4; 
% Constants in Green's function
c1 = 1/4/G0; c2 = (lamada0+G0)/G0/(lamada0+2*G0);
%% Create the stran, stress, effective property, iteration number to save 
Evec = zeros(m*n*l,36);
Svec = zeros(m*n*l,36);
CH = zeros(6,6);
iter = zeros(6,1);
%% Controlling the tolerance of strain refinement stage (fixed point iteration solver)
tolerance_picard = 1e-6;
%% input the disired strain BC
%% case 1 Strain=[1,0,0,0,0,0]'
disp('input the BC1')
Strain=[1,0,0,0,0,0]';
tic % First case start
[s11,s22,s33,s12,s13,s23,e11,e22,e33,e12,e13,e23,CH(1,1),CH(2,1),...
              CH(3,1),CH(4,1),CH(5,1),CH(6,1),iter(1)] = ReDGO(Npmeso,m,n,l,Lx,Ly,Lz,Strain,cons,lamada0,...
              G0,c1,c2,SCAstran(:,1),SCAstress(:,1),ClusterData,tolerance_picard);
toc %First case end 
Svec(:,1)=reshape (s11,[],1);Svec(:,2)=reshape (s22,[],1);Svec(:,3)=reshape (s33,[],1);
Svec(:,4)=reshape (s12,[],1);Svec(:,5)=reshape (s13,[],1);Svec(:,6)=reshape (s23,[],1);
Evec(:,1)=reshape (e11,[],1);Evec(:,2)=reshape (e22,[],1);Evec(:,3)=reshape (e33,[],1);
Evec(:,4)=reshape (e12,[],1);Evec(:,5)=reshape (e13,[],1);Evec(:,6)=reshape (e23,[],1);
clear s11 s22 s33 s12 s13 s23 e11 e22 e33 e12 e13 e23
%% case 2 Strain=[0,1,0,0,0,0]'
disp('input the BC2')
Strain=[0,1,0,0,0,0]';
tic % Second case start 
[s11,s22,s33,s12,s13,s23,e11,e22,e33,e12,e13,e23,CH(1,2),CH(2,2),...
              CH(3,2),CH(4,2),CH(5,2),CH(6,2),iter(2)] = ReDGO(Npmeso,m,n,l,Lx,Ly,Lz,Strain,cons,lamada0,...
              G0,c1,c2,SCAstran(:,2),SCAstress(:,2),ClusterData,tolerance_picard);
toc % second case end 
Svec(:,7)=reshape (s11,[],1);Svec(:,8)=reshape (s22,[],1);Svec(:,9)=reshape (s33,[],1);
Svec(:,10)=reshape (s12,[],1);Svec(:,11)=reshape (s13,[],1);Svec(:,12)=reshape (s23,[],1);
Evec(:,7)=reshape (e11,[],1);Evec(:,8)=reshape (e22,[],1);Evec(:,9)=reshape (e33,[],1);
Evec(:,10)=reshape (e12,[],1);Evec(:,11)=reshape (e13,[],1);Evec(:,12)=reshape (e23,[],1);
clear s11 s22 s33 s12 s13 s23 e11 e22 e33 e12 e13 e23
%% case 3 Strain=[0,0,1,0,0,0]'
disp('input the BC3')
Strain=[0,0,1,0,0,0]';
tic % third case start 
[s11,s22,s33,s12,s13,s23,e11,e22,e33,e12,e13,e23,CH(1,3),CH(2,3),...
              CH(3,3),CH(4,3),CH(5,3),CH(6,3),iter(3)] = ReDGO(Npmeso,m,n,l,Lx,Ly,Lz,Strain,cons,lamada0,...
              G0,c1,c2,SCAstran(:,3),SCAstress(:,3),ClusterData,tolerance_picard);
toc % third case end 
Svec(:,13)=reshape (s11,[],1);Svec(:,14)=reshape (s22,[],1);Svec(:,15)=reshape (s33,[],1);
Svec(:,16)=reshape (s12,[],1);Svec(:,17)=reshape (s13,[],1);Svec(:,18)=reshape (s23,[],1);
Evec(:,13)=reshape (e11,[],1);Evec(:,14)=reshape (e22,[],1);Evec(:,15)=reshape (e33,[],1);
Evec(:,16)=reshape (e12,[],1);Evec(:,17)=reshape (e13,[],1);Evec(:,18)=reshape (e23,[],1);
clear s11 s22 s33 s12 s13 s23 e11 e22 e33 e12 e13 e23
%% case 4 Strain=[0,0,0,1,0,0]'
disp('input the BC4')
Strain=[0,0,0,1,0,0]';
tic % fourth case start 
[s11,s22,s33,s12,s13,s23,e11,e22,e33,e12,e13,e23,CH(1,4),CH(2,4),...
              CH(3,4),CH(4,4),CH(5,4),CH(6,4),iter(4)] = ReDGO(Npmeso,m,n,l,Lx,Ly,Lz,Strain,cons,lamada0,...
              G0,c1,c2,SCAstran(:,4),SCAstress(:,4),ClusterData,tolerance_picard);
toc % fourth case end
Svec(:,19)=reshape (s11,[],1);Svec(:,20)=reshape (s22,[],1);Svec(:,21)=reshape (s33,[],1);
Svec(:,22)=reshape (s12,[],1);Svec(:,23)=reshape (s13,[],1);Svec(:,24)=reshape (s23,[],1);
Evec(:,19)=reshape (e11,[],1);Evec(:,20)=reshape (e22,[],1);Evec(:,21)=reshape (e33,[],1);
Evec(:,22)=reshape (e12,[],1);Evec(:,23)=reshape (e13,[],1);Evec(:,24)=reshape (e23,[],1);
clear s11 s22 s33 s12 s13 s23 e11 e22 e33 e12 e13 e23
%% case 5 Strain=[0,0,0,0,1,0]'
disp('input the BC5')
Strain=[0,0,0,0,1,0]';
tic % 5-th case start 
[s11,s22,s33,s12,s13,s23,e11,e22,e33,e12,e13,e23,CH(1,5),CH(2,5),...
              CH(3,5),CH(4,5),CH(5,5),CH(6,5),iter(5)] = ReDGO(Npmeso,m,n,l,Lx,Ly,Lz,Strain,cons,lamada0,...
              G0,c1,c2,SCAstran(:,5),SCAstress(:,5),ClusterData,tolerance_picard);
toc %5-th case end
Svec(:,25)=reshape (s11,[],1);Svec(:,26)=reshape (s22,[],1);Svec(:,27)=reshape (s33,[],1);
Svec(:,28)=reshape (s12,[],1);Svec(:,29)=reshape (s13,[],1);Svec(:,30)=reshape (s23,[],1);
Evec(:,25)=reshape (e11,[],1);Evec(:,26)=reshape (e22,[],1);Evec(:,27)=reshape (e33,[],1);
Evec(:,28)=reshape (e12,[],1);Evec(:,29)=reshape (e13,[],1);Evec(:,30)=reshape (e23,[],1);
clear s11 s22 s33 s12 s13 s23 e11 e22 e33 e12 e13 e23
%% case 6 Strain=[0,0,0,0,0,1]'
disp('input the BC6')
Strain=[0,0,0,0,0,1]';
tic %6-th case start 
[s11,s22,s33,s12,s13,s23,e11,e22,e33,e12,e13,e23,CH(1,6),CH(2,6),...
              CH(3,6),CH(4,6),CH(5,6),CH(6,6),iter(6)] = ReDGO(Npmeso,m,n,l,Lx,Ly,Lz,Strain,cons,lamada0,...
              G0,c1,c2,SCAstran(:,6),SCAstress(:,6),ClusterData,tolerance_picard);
toc % 6-th case end
Svec(:,31)=reshape (s11,[],1);Svec(:,32)=reshape (s22,[],1);Svec(:,33)=reshape (s33,[],1);
Svec(:,34)=reshape (s12,[],1);Svec(:,35)=reshape (s13,[],1);Svec(:,36)=reshape (s23,[],1);
Evec(:,31)=reshape (e11,[],1);Evec(:,32)=reshape (e22,[],1);Evec(:,33)=reshape (e33,[],1);
Evec(:,34)=reshape (e12,[],1);Evec(:,35)=reshape (e13,[],1);Evec(:,36)=reshape (e23,[],1);
clear s11 s22 s33 s12 s13 s23 e11 e22 e33 e12 e13 e23
save ('Svec.mat','Svec')
save ('Evec.mat','Evec')
save ('iter.mat','iter')
save ('CH.mat','CH')
