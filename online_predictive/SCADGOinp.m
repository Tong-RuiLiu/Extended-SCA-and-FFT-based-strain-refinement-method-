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
%% This is a SCA online prediction stage for DGO (Willot's Green Operator)
format long
clc
clear variables
%% You must load the cluster-wised inputfile
load('vfc-64-16-4.mat')% load volume fraction of each cluster
load('A&U-64-16-4.mat')% Strain concentration tensor and undulation angle 
load('D1willot.mat') % first part of interaction tensor 
load('D2willot.mat') % second part of interaction tensor 
%% Initialize the pattern ID, number, voulme...
Np_m=64;
Np_yarnxpos=64;Np_yarnxneg=64;Np_yarnypos=64;Np_yarnyneg=64;    
Npmeso=Np_m+Np_yarnxpos+Np_yarnxneg+Np_yarnypos+Np_yarnyneg;
mesoinc = 1 ;           % mesostep                  
ndofd  = 6*Npmeso;           % dimension of Dij and MIJ
ntens = 6;                   % stress and strain component
toler = 0.00001;          % tolerance
%% Tensor define
J=[1/3, 1/3, 1/3, 0, 0, 0;
   1/3, 1/3, 1/3, 0, 0, 0;
   1/3, 1/3, 1/3, 0, 0, 0;
   0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0;];% Deviatoric Tensor
II=[eye(3) zeros(3); zeros(3) diag([0.5 0.5 0.5])];% Fourth-Order symmetric identity matrix
K=II-J;
%% Material Elasticity parameter (Matrix)
Em=3.170e+9;
num=0.35;
%% Material Elasticity parameter (Yarn in local coordinate)
E11=5.312e+10;E22=1.4660e+10;E33=E22;
v12=0.266;v13=0.266;v23=0.268;
G12=4.240e+09;G13=G12;G23=5.780e+09;
%% Copy the uyarnpos into each pattern
p=length(Ayarnxpos)/length(uyarnxpos);
temp_uyarnxpos=zeros(4*p,1);
temp_uyarnxneg=zeros(4*p,1);
temp_uyarnypos=zeros(4*p,1);
temp_uyarnyneg=zeros(4*p,1);
for i=1:length(uyarnxpos)
    temp_uyarnxpos(4*i-3:4*i,:)=uyarnxpos(i);
    temp_uyarnxneg(4*i-3:4*i,:)=uyarnxneg(i);
    temp_uyarnypos(4*i-3:4*i,:)=uyarnypos(i);
    temp_uyarnyneg(4*i-3:4*i,:)=uyarnyneg(i);
end
uyarnxpos=temp_uyarnxpos;
uyarnxneg=temp_uyarnxneg;
uyarnypos=temp_uyarnypos;
uyarnyneg=temp_uyarnyneg;
clear temp_uyarnxpos
clear temp_uyarnxneg
clear temp_uyarnypos
clear temp_uyarnyneg
ulist=[ones(Np_m,1);uyarnxpos;uyarnxneg;uyarnypos;uyarnyneg];
%% Define the stiffness matrix for each cluster
C_i=zeros(6,6*Npmeso);
for i=1:Npmeso
    if (i<=Np_m)
        C_i(1:6,6*i-5:6*i)=Elasm3D(Em,num);
    elseif (i>Np_m&&i<=Np_m+Np_yarnxpos+Np_yarnxneg)
        C_i(1:6,6*i-5:6*i)=Elasm3DYarnX(E11,E22,E33,v12,v13,v23,G12,G13,G23,sqrt(1-ulist(i)^2),ulist(i));
    else
        C_i(1:6,6*i-5:6*i)=Elasm3DYarnY(E11,E22,E33,v12,v13,v23,G12,G13,G23,sqrt(1-ulist(i)^2),ulist(i));
    end
end
%% Define the C_0, C_ave matrix, and the interaction tensor 
C_ave=zeros(6,6);
A=[Amatrix;Ayarnxpos;Ayarnxneg;Ayarnypos;Ayarnyneg];
for i=1:Npmeso
    C_ave=C_ave+vfc(i,1)*C_i(1:6,6*i-5:6*i)*reshape(A(i,:),6,[]);
end
mu0=1/10*(sum(sum(K.*C_ave)));% Shear Modulus
lamada0=(sum(sum(J.*C_ave))-2*mu0)/3;% Lame Constant(1)
k0=lamada0+2/3*mu0;% Bulk   Modulus
C_0=sum(sum(J.*C_ave))*J+1/5*(sum(sum(K.*C_ave)))*K;% Reference material Stiffness matrix
D=1/4/mu0*D1+(lamada0+mu0)/mu0/(lamada0+2*mu0)*D2;
%% Generate the initial stiffness matrix
MIJmeso=zeros(ntens*Npmeso,ntens*Npmeso);% System Jacobian Matrix
for k = 1: Npmeso
    dCal_diag(1:6,6*k-5:6*k)=C_i(1:6,6*k-5:6*k)-C_0(1:6,1:6);
end 
for i=1:Npmeso
    MIJmeso(6*i-5:6*i,6*i-5:6*i)=MIJmeso(6*i-5:6*i,6*i-5:6*i)+eye(6,6);
    for kk=1:Npmeso
        MIJmeso(6*i-5:6*i,6*kk-5:6*kk)=MIJmeso(6*i-5:6*i,6*kk-5:6*kk)+D(6*i-5:6*i,6*kk-5:6*kk)*dCal_diag(1:6,6*kk-5:6*kk);
    end
end
%% Create a array to save all the computational results
SCAstran = zeros (ndofd,ntens);
SCAstress = zeros (ndofd,ntens);
SCAsave  = zeros(ntens,ntens);
SCAeave  = zeros(ntens,ntens);
%% Apply the loading, which contains 6 loading conditions 
for k  = 1 : ntens
    strinc = zeros(ntens,1);
    strinc(k) = 1;
    stranold =  zeros(ndofd+ntens,1); % Initial condition
    tic 
    [stranold,sigma]=SCAmeso(stranold,strinc,MIJmeso,Npmeso,ndofd,ntens,mesoinc,C_i,C_0,D,toler);
    toc
    SCAsave(:,k)=sum(reshape(sigma,6,Npmeso)*diag(vfc(:,1)),2);
    SCAeave(:,k)=sum(reshape(stranold(7:end),6,Npmeso)*diag(vfc(:,1)),2);
    SCAstran(:,k) = stranold(7:end);
    SCAstress(:,k) = sigma;
end
save('SCAstran.mat','SCAstran')
save('SCAstress.mat','SCAstress')
save('SCAsave.mat','SCAsave')
save('SCAeave.mat','SCAeave')