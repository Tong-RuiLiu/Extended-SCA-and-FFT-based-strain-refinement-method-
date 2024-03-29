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
function [s11,s22,s33,s12,s13,s23,e11,e22,e33,e12,e13,e23,S11ave,S22ave,...
              S33ave,S12ave,S13ave,S23ave,iter] = ReDGO(Np,m,n,l,Lx,Ly,Lz,Strain,cons,lamada0,...
              G0,c1,c2,SCAstran,SCAstress,ClusterData,tol_picard)
% Initialize the stress strain and polarization stress
e11= zeros(m*n*l,1);e22= zeros(m*n*l,1);e33= zeros(m*n*l,1);e12= zeros(m*n*l,1);e13= zeros(m*n*l,1);e23= zeros(m*n*l,1);
s11= zeros(m*n*l,1);s22= zeros(m*n*l,1);s33= zeros(m*n*l,1);s12= zeros(m*n*l,1);s13= zeros(m*n*l,1);s23= zeros(m*n*l,1);
p11= zeros(m*n*l,1);p22= zeros(m*n*l,1);p33= zeros(m*n*l,1);p12= zeros(m*n*l,1);p13= zeros(m*n*l,1);p23= zeros(m*n*l,1);
for i = 1:Np
    % Distribute the stress strain ...
    id = find(ClusterData(:,8)==i);
    s11(id) = SCAstress(6*i-5);s22(id) = SCAstress(6*i-4);s33(id) =SCAstress(6*i-3);
    s12(id) = SCAstress(6*i-2);s13(id) = SCAstress(6*i-1);s23(id) = SCAstress(6*i);
    
    e11(id) = SCAstran(6*i-5);e22(id) = SCAstran(6*i-4);e33(id) = SCAstran(6*i-3);
    e12(id) = SCAstran(6*i-2);e13(id) = SCAstran(6*i-1);e23(id) = SCAstran(6*i);  
end
% Modify the macro strain 
Emacro11 = Strain(1)*ones(m,n,l);
Emacro22 = Strain(2)*ones(m,n,l);
Emacro33 = Strain(3)*ones(m,n,l);
Emacro12 = Strain(4)*ones(m,n,l);
Emacro13 = Strain(5)*ones(m,n,l);
Emacro23 = Strain(6)*ones(m,n,l);
%% Compute the Green function in the frequency domain
%Define the Green's Operator 
T11 = zeros(m,n,l);T12 = zeros(m,n,l);T13 = zeros(m,n,l);T14 = zeros(m,n,l);T15 = zeros(m,n,l);T16 = zeros(m,n,l);
T22 = zeros(m,n,l);T23 = zeros(m,n,l);T24 = zeros(m,n,l);T25 = zeros(m,n,l);T26 = zeros(m,n,l);
T33 = zeros(m,n,l);T34 = zeros(m,n,l);T35 = zeros(m,n,l);T36 = zeros(m,n,l);
T44 = zeros(m,n,l);T45 = zeros(m,n,l);T46 = zeros(m,n,l);
T55 = zeros(m,n,l);T56 = zeros(m,n,l);
T66 = zeros(m,n,l);
% Discrete frequency point, assuming the length of each element is 1
%%   No-Need FFTshift and Ifftshift (Nyquist = 0)
em =2*pi*[0:m/2-1,-m/2:-1]/Lx;
en =2*pi*[0:n/2-1,-n/2:-1]/Ly;
el = 2*pi*[0:l/2-1,-l/2:-1]/Lz;
E0= G0*(3*lamada0+2*G0)/(lamada0+G0);
v0=lamada0/(2*(lamada0+G0));
C0 = E0/((1+v0)*(1-2*v0))* [1-v0,v0,v0,0,0,0;
                                    v0,1-v0,v0,0,0,0;
                                    v0,v0,1-v0,0,0,0;
                                    0,0,0,(1-2*v0)/2,0,0;
                                    0,0,0,0,(1-2*v0)/2,0;
                                    0,0,0,0,0,(1-2*v0)/2];
Diff = cons- kron (reshape(C0,1,[]),ones(m*n*l,1));
%% Distribute the polarization stress,strain, stress
for i = 1:m*n*l
    p11(i) = s11(i) - (C0(1,1)*e11(i)+C0(1,2)*e22(i)+C0(1,3)*e33(i)+C0(1,4)*e12(i)+C0(1,5)*e13(i)+C0(1,6)*e23(i));
    p22(i) = s22(i) - (C0(2,1)*e11(i)+C0(2,2)*e22(i)+C0(2,3)*e33(i)+C0(2,4)*e12(i)+C0(2,5)*e13(i)+C0(2,6)*e23(i));
    p33(i) = s33(i) - (C0(3,1)*e11(i)+C0(3,2)*e22(i)+C0(3,3)*e33(i)+C0(3,4)*e12(i)+C0(3,5)*e13(i)+C0(3,6)*e23(i));
    p12(i) = s12(i) - (C0(4,1)*e11(i)+C0(4,2)*e22(i)+C0(4,3)*e33(i)+C0(4,4)*e12(i)+C0(4,5)*e13(i)+C0(4,6)*e23(i));
    p13(i) = s13(i) - (C0(5,1)*e11(i)+C0(5,2)*e22(i)+C0(5,3)*e33(i)+C0(5,4)*e12(i)+C0(5,5)*e13(i)+C0(5,6)*e23(i));
    p23(i) = s23(i) - (C0(6,1)*e11(i)+C0(6,2)*e22(i)+C0(6,3)*e33(i)+C0(6,4)*e12(i)+C0(6,5)*e13(i)+C0(6,6)*e23(i));
end
disp('initialize the (polarization) stress at each element/Real space')
s11 = reshape(s11,[m,n,l]);s22= reshape(s22,[m,n,l]);s33 = reshape(s33,[m,n,l]);
s12 = reshape(s12,[m,n,l]);s13= reshape(s13,[m,n,l]);s23 = reshape(s23,[m,n,l]);
p11 = reshape(p11,[m,n,l]);p22= reshape(p22,[m,n,l]);p33 = reshape(p33,[m,n,l]);
p12 = reshape(p12,[m,n,l]);p13= reshape(p13,[m,n,l]);p23 = reshape(p23,[m,n,l]);
disp('initialize the strain at each element/Real space')
e11 = reshape(e11,[m,n,l]);e22= reshape(e22,[m,n,l]);e33 = reshape(e33,[m,n,l]);
e12 = reshape(e12,[m,n,l]);e13= reshape(e13,[m,n,l]);e23 = reshape(e23,[m,n,l]);
%% Reference materal and Green's Function construction
for i = 1:m
    for j = 1:n
        for k = 1:l
            %% Willot's discrete Green's operator
            e1 = 2*m/Lx*sin(em(i)/2)*cos(en(j)/2)*cos(el(k)/2);
            e2 = 2*n/Ly*cos(em(i)/2)*sin(en(j)/2)*cos(el(k)/2);
            e3 = 2*l/Lz*cos(em(i)/2)*cos(en(j)/2)*sin(el(k)/2);
            %% Moulinec-Suquet continuous Green's operator
            %e1 = em(i); e2 = en(j); e3 = el(k);
            %% Compute the Green's operators
            E = norm([e1, e2, e3]);
            T11(i,j,k) = c1/E^2*(4*e1^2)-c2*e1^4/E^4;
            T12(i,j,k) = -c2*e1^2*e2^2/E^4;
            T13(i,j,k) = -c2*e1^2*e3^2/E^4;
            T14(i,j,k) = c1/E^2*(2*e1*e2)-c2*e1^3*e2/E^4;
            T15(i,j,k) = c1/E^2*(2*e1*e3)-c2*e1^3*e3/E^4;
            T16(i,j,k) = -c2*e1^2*e2*e3/E^4;
            
            T22(i,j,k) = c1/E^2*(4*e2^2)-c2*e2^4/E^4;
            T23(i,j,k) = -c2*e2^2*e3^2/E^4;
            T24(i,j,k) = c1/E^2*(2*e1*e2)-c2*e2^3*e1/E^4;
            T25(i,j,k) = -c2*e2^2*e1*e3/E^4;
            T26(i,j,k) = c1/E^2*(2*e2*e3)-c2*e2^3*e3/E^4;
            
            T33(i,j,k) = c1/E^2*(4*e3^2)-c2*e3^4/E^4;
            T34(i,j,k) =-c2*e3^2*e1*e2/E^4;
            T35(i,j,k) = c1/E^2*(2*e3*e1)-c2*e3^3*e1/E^4;
            T36(i,j,k) = c1/E^2*(2*e3*e2)-c2*e3^3*e2/E^4;
            
            T44(i,j,k) = c1/E^2*(e1^2+e2^2)-c2*e1^2*e2^2/E^4;
            T45(i,j,k) = c1/E^2*(e2*e3)-c2*e1^2*e2*e3/E^4;
            T46(i,j,k) = c1/E^2*(e1*e3)-c2*e2^2*e1*e3/E^4;
            
            T55(i,j,k) = c1/E^2*(e1^2+e3^2)-c2*e1^2*e3^2/E^4;
            T56(i,j,k) = c1/E^2*(e1*e2)-c2*e3^2*e1*e2/E^4;
            
            T66(i,j,k) = c1/E^2*(e2^2+e3^2)-c2*e2^2*e3^2/E^4;
        end
    end
end
% modify the Green function tensor according to voigt notation
T14 = T14*2; T15 = T15*2; T16 = T16*2; T24 = T24*2; T25 = T25*2; T26 = T26*2;
T34 = T34*2; T35 = T35*2; T36 = T36*2; 
T44 = T44*4; T45 = T45*4; T46 = T46*4; T55 = T55*4; T56 = T56*4; T66 = T66*4;
% Complete the symmetric part
T21 = T12; T31 = T13; T41 = T14; T51 = T15; T61=T16;
                  T32 = T23; T42 = T24; T52 = T25; T62=T26;
                                    T43 = T34; T53 = T35; T63=T36;
                                                      T54 = T45; T64=T46;
                                                                        T65= T56;
% Modify the Green tensor at the zero freqency
disp('Modify the Green tensor at the zero freqency')
cpx = 1;cpy = 1;cpz= 1; % Does not need FFT shift at all!!!
T11(cpx,cpy,cpz) = 0; T12(cpx,cpy,cpz) = 0; T13(cpx,cpy,cpz) = 0; T14(cpx,cpy,cpz) = 0; T15(cpx,cpy,cpz) = 0; T16(cpx,cpy,cpz) = 0;
T21(cpx,cpy,cpz) = 0; T22(cpx,cpy,cpz) = 0; T23(cpx,cpy,cpz) = 0; T24(cpx,cpy,cpz) = 0; T25(cpx,cpy,cpz) = 0; T26(cpx,cpy,cpz) = 0;
T31(cpx,cpy,cpz) = 0; T32(cpx,cpy,cpz) = 0; T33(cpx,cpy,cpz) = 0; T34(cpx,cpy,cpz) = 0; T35(cpx,cpy,cpz) = 0; T36(cpx,cpy,cpz) = 0;
T41(cpx,cpy,cpz) = 0; T42(cpx,cpy,cpz) = 0; T43(cpx,cpy,cpz) = 0; T44(cpx,cpy,cpz) = 0; T45(cpx,cpy,cpz) = 0; T46(cpx,cpy,cpz) = 0;
T51(cpx,cpy,cpz) = 0; T52(cpx,cpy,cpz) = 0; T53(cpx,cpy,cpz) = 0; T54(cpx,cpy,cpz) = 0; T55(cpx,cpy,cpz) = 0; T56(cpx,cpy,cpz) = 0;
T61(cpx,cpy,cpz) = 0; T62(cpx,cpy,cpz) = 0; T63(cpx,cpy,cpz) = 0; T64(cpx,cpy,cpz) = 0; T65(cpx,cpy,cpz) = 0; T66(cpx,cpy,cpz) = 0;
%%  FFT on the strain
disp('FFT on initialize the strain at each element/Fourier space') 
Snew= sum(sum(sum(sqrt(s11.^2+s22.^2+s33.^2+s12.^2+s13.^2+s23.^2))))/m/n/l;
Sold = 0;iter = 0;

%% Go to the iteration for FFT main program
disp('Start the iteration')
while abs(Snew-Sold)/Snew >tol_picard

    iter = iter+1;
    disp(['iter=',num2str(iter)]);
    % FFT on the stress   
    p11fft = (fftn(p11)); p22fft = (fftn(p22)); p33fft = (fftn(p33));  
    p12fft = (fftn(p12)); p13fft = (fftn(p13)); p23fft = (fftn(p23));
    % Update the strain
    e11 = Emacro11 - real(ifftn((T11.*p11fft + T12.*p22fft + T13.*p33fft) + (T14.*p12fft + T15.*p13fft + T16.*p23fft)));
    e22 = Emacro22 - real(ifftn((T21.*p11fft + T22.*p22fft + T23.*p33fft) + (T24.*p12fft + T25.*p13fft + T26.*p23fft)));
    e33 = Emacro33 - real(ifftn((T31.*p11fft + T32.*p22fft + T33.*p33fft) + (T34.*p12fft + T35.*p13fft + T36.*p23fft)));
    e12 = Emacro12 - real(ifftn((T41.*p11fft + T42.*p22fft + T43.*p33fft) + (T44.*p12fft + T45.*p13fft + T46.*p23fft)));
    e13 = Emacro13 - real(ifftn((T51.*p11fft + T52.*p22fft + T53.*p33fft) + (T54.*p12fft + T55.*p13fft + T56.*p23fft)));
    e23 = Emacro23 - real(ifftn((T61.*p11fft + T62.*p22fft + T63.*p33fft) + (T64.*p12fft + T65.*p13fft + T66.*p23fft)));
    %% For Parallel Computing
    e11 = reshape(e11,[],1);s11 = reshape(s11,[],1);p11 = reshape(p11,[],1);
    e22 = reshape(e22,[],1);s22 = reshape(s22,[],1);p22 = reshape(p22,[],1);
    e33 = reshape(e33,[],1);s33 = reshape(s33,[],1);p33 = reshape(p33,[],1);
    e12 = reshape(e12,[],1);s12 = reshape(s12,[],1);p12 = reshape(p12,[],1);
    e13 = reshape(e13,[],1);s13 = reshape(s13,[],1);p13 = reshape(p13,[],1);
    e23 = reshape(e23,[],1);s23 = reshape(s23,[],1);p23 = reshape(p23,[],1);
    parfor mm  = 1: m*n*l
        Con = reshape(cons(mm,:),[6,6]);
        Dif = reshape(Diff(mm,:),[6,6]);
        s11(mm) = Con(1,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        s22(mm) = Con(2,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        s33(mm) = Con(3,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        s12(mm) = Con(4,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        s13(mm) = Con(5,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        s23(mm) = Con(6,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        p11(mm) = Dif(1,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        p22(mm) = Dif(2,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        p33(mm) = Dif(3,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        p12(mm) = Dif(4,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        p13(mm) = Dif(5,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
        p23(mm) = Dif(6,:)*[e11(mm);e22(mm);e33(mm);e12(mm);e13(mm);e23(mm)];
    end
    e11 = reshape(e11,[m,n,l]);s11 = reshape(s11,[m,n,l]);p11 = reshape(p11,[m,n,l]);
    e22 = reshape(e22,[m,n,l]);s22 = reshape(s22,[m,n,l]);p22 = reshape(p22,[m,n,l]);
    e33 = reshape(e33,[m,n,l]);s33 = reshape(s33,[m,n,l]);p33 = reshape(p33,[m,n,l]);
    e12 = reshape(e12,[m,n,l]);s12 = reshape(s12,[m,n,l]);p12 = reshape(p12,[m,n,l]);
    e13 = reshape(e13,[m,n,l]);s13 = reshape(s13,[m,n,l]);p13 = reshape(p13,[m,n,l]);
    e23 = reshape(e23,[m,n,l]);s23 = reshape(s23,[m,n,l]);p23 = reshape(p23,[m,n,l]);
    Sold = Snew;
    Snew = sum(sum(sum(sqrt(s11.^2+s22.^2+s33.^2+s12.^2+s13.^2+s23.^2))))/m/n/l;
end
%% compute the average stress of RVE
disp('Compute the average/homogenized stress')
S11ave=mean(reshape(s11,[m*n*l,1]));
S22ave=mean(reshape(s22,[m*n*l,1]));
S33ave=mean(reshape(s33,[m*n*l,1]));
S12ave=mean(reshape(s12,[m*n*l,1]));
S13ave=mean(reshape(s13,[m*n*l,1]));
S23ave=mean(reshape(s23,[m*n*l,1]));

end