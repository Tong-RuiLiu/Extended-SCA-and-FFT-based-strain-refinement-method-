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
%% This is the constitutive laws for yarn materials with centerline stretches in the X direction. 
function[Cglobal]=Elasm3DYarnX(E11,E22,E33,v12,v13,v23,G12,G13,G23,csi,si)
Clocal=zeros(6);
rmx=[csi,0,si;0,cos(pi),0;si,0,-csi];% Yarn0 Yarn1 Yarn5 (YarnX)
%rmx=[csi,0,-si;0,1,0;si,0,csi];% Yarn0 Yarn1 Yarn5 (YarnX)
%% Local coordinate stiffness matrix (Transverse Isotropic material)
Clocal(1,1)=1/E11;
Clocal(2,2)=1/E22;
Clocal(3,3)=1/E33;
Clocal(2,1)=-v12/E11;
Clocal(3,1)=-v13/E11;
Clocal(3,2)=-v23/E22;
Clocal(1,2)=Clocal(2,1);
Clocal(1,3)=Clocal(3,1);
Clocal(2,3)=Clocal(3,2);
Clocal(4,4)=1/G12;
Clocal(5,5)=1/G13;
Clocal(6,6)=1/G23;
Clocal=inv(Clocal);
%% Stress transformation matrix(From Local to global)
T=[rmx(1,1)^2,rmx(1,2)^2,rmx(1,3)^2,2*rmx(1,1)*rmx(1,2),2*rmx(1,1)*rmx(1,3),2*rmx(1,2)*rmx(1,3);
   rmx(2,1)^2,rmx(2,2)^2,rmx(2,3)^2,2*rmx(2,1)*rmx(2,2),2*rmx(2,1)*rmx(2,3),2*rmx(2,2)*rmx(2,3);
   rmx(3,1)^2,rmx(3,2)^2,rmx(3,3)^2,2*rmx(3,1)*rmx(3,2),2*rmx(3,1)*rmx(3,3),2*rmx(3,2)*rmx(3,3);
   rmx(1,1)*rmx(2,1),rmx(1,2)*rmx(2,2),rmx(1,3)*rmx(2,3),rmx(1,1)*rmx(2,2)+rmx(1,2)*rmx(2,1),rmx(1,1)*rmx(2,3)+rmx(1,3)*rmx(2,1),rmx(1,2)*rmx(2,3)+rmx(1,3)*rmx(2,2);
   rmx(3,1)*rmx(1,1),rmx(1,2)*rmx(3,2),rmx(3,3)*rmx(1,3),rmx(3,1)*rmx(1,2)+rmx(1,1)*rmx(3,2),rmx(3,1)*rmx(1,3)+rmx(3,3)*rmx(1,1),rmx(3,2)*rmx(1,3)+rmx(3,3)*rmx(1,2);
   rmx(2,1)*rmx(3,1),rmx(2,2)*rmx(3,2),rmx(3,3)*rmx(2,3),rmx(2,1)*rmx(3,2)+rmx(2,2)*rmx(3,1),rmx(3,1)*rmx(2,3)+rmx(3,3)*rmx(2,1),rmx(2,2)*rmx(3,3)+rmx(2,3)*rmx(3,2);];
%% Calcuate elasticity matrix in global coordinate system
Cglobal=(T*Clocal)*transpose(T);