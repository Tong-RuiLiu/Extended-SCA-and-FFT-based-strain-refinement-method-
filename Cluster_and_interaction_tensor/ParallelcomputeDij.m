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
function [D1,D2] = ParallelcomputeDij(Ncluster,m,n,l,Lx,Ly,Lz,list,idx,indicator)
%% Wave vector (continious) define to decide need FFTshift or not. Herein, need FFTshift
ks1 =2*pi*(ceil(-m/2):floor((m-1)/2))/Lx;
ks2 =2*pi*(ceil(-n/2):floor((n-1)/2))/Ly;
ks3= 2*pi*(ceil(-l/2):floor((l-1)/2))/Lz;
%% Calculate the interaction tensors between the patterns using fft method
disp('Interaction tensor Dij calculation in Fourier space via DFT...');
%%  First part of the interaction tensor ---------------------------------------
disp('Start Calculation of the interaction tensor  the 1st and 2nd part!');
% Update the Green's function
Ghat11 = zeros(m,n,l);Ghat12 = zeros(m,n,l);Ghat13 = zeros(m,n,l);Ghat14 = zeros(m,n,l);Ghat15 = zeros(m,n,l);Ghat16 = zeros(m,n,l);
Ghat22 = zeros(m,n,l);Ghat23 = zeros(m,n,l);Ghat24 = zeros(m,n,l);Ghat25 = zeros(m,n,l);Ghat26 = zeros(m,n,l);
Ghat33 = zeros(m,n,l);Ghat34 = zeros(m,n,l);Ghat35 = zeros(m,n,l);Ghat36 = zeros(m,n,l);
Ghat44 = zeros(m,n,l);Ghat45 = zeros(m,n,l);Ghat46 = zeros(m,n,l);
Ghat55 = zeros(m,n,l);Ghat56 = zeros(m,n,l);
Ghat66 = zeros(m,n,l);
G2hat11 = zeros(m,n,l);G2hat12 = zeros(m,n,l);G2hat13 = zeros(m,n,l);G2hat14 = zeros(m,n,l);G2hat15 = zeros(m,n,l);G2hat16 = zeros(m,n,l);
G2hat22 = zeros(m,n,l);G2hat23 = zeros(m,n,l);G2hat24 = zeros(m,n,l);G2hat25 = zeros(m,n,l);G2hat26 = zeros(m,n,l);
G2hat33 = zeros(m,n,l);G2hat34 = zeros(m,n,l);G2hat35 = zeros(m,n,l);G2hat36 = zeros(m,n,l);
G2hat44 = zeros(m,n,l);G2hat45 = zeros(m,n,l);G2hat46 = zeros(m,n,l);
G2hat55 = zeros(m,n,l);G2hat56 = zeros(m,n,l);
G2hat66 = zeros(m,n,l);
for i = 1:m
    for j = 1:n
        for k = 1:l
            if indicator == 1 % CGO
                e1 = ks1(i); e2 = ks2(j); e3 = ks3(k);
            elseif indicator == 2 %DGO
                e1 = 2*m/Lx*sin(ks1(i)/2)*cos(ks2(j)/2)*cos(ks3(k)/2); % Willot's Scheme based on rotated Grid
                e2 = 2*n/Ly*cos(ks1(i)/2)*sin(ks2(j)/2)*cos(ks3(k)/2);
                e3 = 2*l/Lz*cos(ks1(i)/2)*cos(ks2(j)/2)*sin(ks3(k)/2);
            end
            E = norm([e1, e2, e3]);%(Norm of kesi Vector)
            Ghat11(i,j,k) = 1/E^2*(4*e1^2);
            Ghat12(i,j,k) = 0;
            Ghat13(i,j,k) = 0;
            Ghat14(i,j,k) = 1/E^2*(2*e1*e2);
            Ghat15(i,j,k) = 1/E^2*(2*e1*e3);
            Ghat16(i,j,k) = 0;
            Ghat22(i,j,k) = 1/E^2*(4*e2^2);
            Ghat23(i,j,k) = 0;
            Ghat24(i,j,k) = 1/E^2*(2*e1*e2);
            Ghat25(i,j,k) = 0;
            Ghat26(i,j,k) = 1/E^2*(2*e2*e3);
            Ghat33(i,j,k) = 1/E^2*(4*e3^2);
            Ghat34(i,j,k) = 0;
            Ghat35(i,j,k) = 1/E^2*(2*e3*e1);
            Ghat36(i,j,k) = 1/E^2*(2*e3*e2);
            Ghat44(i,j,k) = 1/E^2*(e1^2+e2^2);
            Ghat45(i,j,k) = 1/E^2*(e2*e3);
            Ghat46(i,j,k) = 1/E^2*(e1*e3);
            Ghat55(i,j,k) = 1/E^2*(e1^2+e3^2);
            Ghat56(i,j,k) = 1/E^2*(e1*e2);
            Ghat66(i,j,k) = 1/E^2*(e2^2+e3^2);
            %% Second Part Green's Operator
            G2hat11(i,j,k) = -e1^4/E^4;
            G2hat12(i,j,k) = -e1^2*e2^2/E^4;
            G2hat13(i,j,k) = -e1^2*e3^2/E^4;
            G2hat14(i,j,k) = -e1^3*e2/E^4;
            G2hat15(i,j,k) = -e1^3*e3/E^4;
            G2hat16(i,j,k) = -e1^2*e2*e3/E^4;
            G2hat22(i,j,k) = -e2^4/E^4;
            G2hat23(i,j,k) = -e2^2*e3^2/E^4;
            G2hat24(i,j,k) = -e2^3*e1/E^4;
            G2hat25(i,j,k) = -e2^2*e1*e3/E^4;
            G2hat26(i,j,k) = -e2^3*e3/E^4;
            G2hat33(i,j,k) = -e3^4/E^4;
            G2hat34(i,j,k) = -e3^2*e1*e2/E^4;
            G2hat35(i,j,k) = -e3^3*e1/E^4;
            G2hat36(i,j,k) = -e3^3*e2/E^4;
            G2hat44(i,j,k) = -e1^2*e2^2/E^4;
            G2hat45(i,j,k) = -e1^2*e2*e3/E^4;
            G2hat46(i,j,k) = -e2^2*e1*e3/E^4;
            G2hat55(i,j,k) = -e1^2*e3^2/E^4;
            G2hat56(i,j,k) = -e3^2*e1*e2/E^4;
            G2hat66(i,j,k) = -e2^2*e3^2/E^4;
        end
    end
end
%% modify the Green function tensor according to voigt notation
% first part
Ghat14 = Ghat14*2; Ghat15 = Ghat15*2; Ghat16 = Ghat16*2; Ghat24 = Ghat24*2; Ghat25 = Ghat25*2; Ghat26 = Ghat26*2;
Ghat34 = Ghat34*2; Ghat35 = Ghat35*2; Ghat36 = Ghat36*2; 
Ghat44 = Ghat44*4; Ghat45 = Ghat45*4; Ghat46 = Ghat46*4; Ghat55 = Ghat55*4; Ghat56 = Ghat56*4; Ghat66 = Ghat66*4;
% second part
G2hat14 = G2hat14*2; G2hat15 = G2hat15*2; G2hat16 = G2hat16*2; G2hat24 = G2hat24*2; G2hat25 = G2hat25*2; G2hat26 = G2hat26*2;
G2hat34 = G2hat34*2; G2hat35 = G2hat35*2; G2hat36 = G2hat36*2; 
G2hat44 = G2hat44*4; G2hat45 = G2hat45*4; G2hat46 = G2hat46*4; G2hat55 = G2hat55*4; G2hat56 = G2hat56*4; G2hat66 = G2hat66*4;
%% Modify the Green tensor at the 0 freqency
cpx = ceil((m+1)/2); cpy = ceil((n+1)/2); cpz = ceil((l+1)/2);
%first part
Ghat11(cpx,cpy,cpz) = 0; Ghat12(cpx,cpy,cpz) = 0; Ghat13(cpx,cpy,cpz) = 0; Ghat14(cpx,cpy,cpz) = 0; Ghat15(cpx,cpy,cpz) = 0; Ghat16(cpx,cpy,cpz) = 0;
Ghat22(cpx,cpy,cpz) = 0; Ghat23(cpx,cpy,cpz) = 0; Ghat24(cpx,cpy,cpz) = 0; Ghat25(cpx,cpy,cpz) = 0; Ghat26(cpx,cpy,cpz) = 0;
Ghat33(cpx,cpy,cpz) = 0; Ghat34(cpx,cpy,cpz) = 0; Ghat35(cpx,cpy,cpz) = 0; Ghat36(cpx,cpy,cpz) = 0;
Ghat44(cpx,cpy,cpz) = 0; Ghat45(cpx,cpy,cpz) = 0; Ghat46(cpx,cpy,cpz) = 0;
Ghat55(cpx,cpy,cpz) = 0; Ghat56(cpx,cpy,cpz) = 0;
Ghat66(cpx,cpy,cpz) = 0;
%second part
G2hat11(cpx,cpy,cpz) = 0; G2hat12(cpx,cpy,cpz) = 0; G2hat13(cpx,cpy,cpz) = 0; G2hat14(cpx,cpy,cpz) = 0; G2hat15(cpx,cpy,cpz) = 0; G2hat16(cpx,cpy,cpz) = 0;
G2hat22(cpx,cpy,cpz) = 0; G2hat23(cpx,cpy,cpz) = 0; G2hat24(cpx,cpy,cpz) = 0; G2hat25(cpx,cpy,cpz) = 0; G2hat26(cpx,cpy,cpz) = 0;
G2hat33(cpx,cpy,cpz) = 0; G2hat34(cpx,cpy,cpz) = 0; G2hat35(cpx,cpy,cpz) = 0; G2hat36(cpx,cpy,cpz) = 0;
G2hat44(cpx,cpy,cpz) = 0; G2hat45(cpx,cpy,cpz) = 0; G2hat46(cpx,cpy,cpz) = 0;
G2hat55(cpx,cpy,cpz) = 0; G2hat56(cpx,cpy,cpz) = 0;
G2hat66(cpx,cpy,cpz) = 0;
%% Define the interaction tensor Database (Initialization)
D1= zeros(6*Ncluster,6*Ncluster);
D2= zeros(6*Ncluster,6*Ncluster);
%% First part of Dij
Done11= zeros(Ncluster,Ncluster);Done12= zeros(Ncluster,Ncluster);Done13= zeros(Ncluster,Ncluster);Done14= zeros(Ncluster,Ncluster);Done15= zeros(Ncluster,Ncluster);Done16= zeros(Ncluster,Ncluster);
Done22= zeros(Ncluster,Ncluster);Done23= zeros(Ncluster,Ncluster);Done24= zeros(Ncluster,Ncluster);Done25= zeros(Ncluster,Ncluster);Done26= zeros(Ncluster,Ncluster);
Done33= zeros(Ncluster,Ncluster);Done34= zeros(Ncluster,Ncluster);Done35= zeros(Ncluster,Ncluster);Done36= zeros(Ncluster,Ncluster);
Done44= zeros(Ncluster,Ncluster);Done45= zeros(Ncluster,Ncluster);Done46= zeros(Ncluster,Ncluster);
Done55= zeros(Ncluster,Ncluster);Done56= zeros(Ncluster,Ncluster);
Done66= zeros(Ncluster,Ncluster);
%% Second part of Dij
Dtwo11= zeros(Ncluster,Ncluster);Dtwo12= zeros(Ncluster,Ncluster);Dtwo13= zeros(Ncluster,Ncluster);Dtwo14= zeros(Ncluster,Ncluster);Dtwo15= zeros(Ncluster,Ncluster);Dtwo16= zeros(Ncluster,Ncluster);
Dtwo22= zeros(Ncluster,Ncluster);Dtwo23= zeros(Ncluster,Ncluster);Dtwo24= zeros(Ncluster,Ncluster);Dtwo25= zeros(Ncluster,Ncluster);Dtwo26= zeros(Ncluster,Ncluster);
Dtwo33= zeros(Ncluster,Ncluster);Dtwo34= zeros(Ncluster,Ncluster);Dtwo35= zeros(Ncluster,Ncluster);Dtwo36= zeros(Ncluster,Ncluster);
Dtwo44= zeros(Ncluster,Ncluster);Dtwo45= zeros(Ncluster,Ncluster);Dtwo46= zeros(Ncluster,Ncluster);
Dtwo55= zeros(Ncluster,Ncluster);Dtwo56= zeros(Ncluster,Ncluster);
Dtwo66= zeros(Ncluster,Ncluster);
parfor i = 1:Ncluster
    disp([num2str(i),'-th clusters in UD fiber']);
    Chivector = zeros(m*n*l,1);
    Chivector(list) = (idx==i);
    Chi = reshape(Chivector,m,n,l);
    fftChi = fftn(Chi);
    fftChi = fftshift(fftChi);
    %%  FFT transformation for the first part 
    fftD11 = Ghat11.*fftChi; fftD12 = Ghat12.*fftChi; fftD13 = Ghat13.*fftChi;fftD14 = Ghat14.*fftChi; fftD15 = Ghat15.*fftChi; fftD16 = Ghat16.*fftChi;
    fftD22 = Ghat22.*fftChi; fftD23 = Ghat23.*fftChi;fftD24 = Ghat24.*fftChi; fftD25 = Ghat25.*fftChi; fftD26 = Ghat26.*fftChi;
    fftD33 = Ghat33.*fftChi;fftD34 = Ghat34.*fftChi; fftD35 = Ghat35.*fftChi; fftD36 = Ghat36.*fftChi;
    fftD44 = Ghat44.*fftChi; fftD45 = Ghat45.*fftChi; fftD46 = Ghat46.*fftChi;
    fftD55 = Ghat55.*fftChi; fftD56 = Ghat56.*fftChi;
    fftD66 = Ghat66.*fftChi;
    % Reverse back
    D11 = ifftn(ifftshift(fftD11)); D12 = ifftn(ifftshift(fftD12)); D13 = ifftn(ifftshift(fftD13)); D14 = ifftn(ifftshift(fftD14)); D15 = ifftn(ifftshift(fftD15)); D16 = ifftn(ifftshift(fftD16));
    D22 = ifftn(ifftshift(fftD22)); D23 = ifftn(ifftshift(fftD23)); D24 = ifftn(ifftshift(fftD24)); D25 = ifftn(ifftshift(fftD25)); D26 = ifftn(ifftshift(fftD26));
    D33 = ifftn(ifftshift(fftD33)); D34 = ifftn(ifftshift(fftD34)); D35 = ifftn(ifftshift(fftD35)); D36 = ifftn(ifftshift(fftD36));
    D44 = ifftn(ifftshift(fftD44)); D45 = ifftn(ifftshift(fftD45)); D46 = ifftn(ifftshift(fftD46));
    D55 = ifftn(ifftshift(fftD55)); D56 = ifftn(ifftshift(fftD56));
    D66 = ifftn(ifftshift(fftD66));
    % If m, n, l are even, D_ij (i ~= j) are complex numbers, choose the real part
    D11 = real(reshape(D11,[],1)); D12 = real(reshape(D12,[],1)); D13 = real(reshape(D13,[],1)); D14 = real(reshape(D14,[],1)); D15 = real(reshape(D15,[],1)); D16 = real(reshape(D16,[],1));
    D22 = real(reshape(D22,[],1)); D23 = real(reshape(D23,[],1)); D24 = real(reshape(D24,[],1)); D25 = real(reshape(D25,[],1)); D26 = real(reshape(D26,[],1));
    D33 = real(reshape(D33,[],1)); D34 = real(reshape(D34,[],1)); D35 = real(reshape(D35,[],1)); D36 = real(reshape(D36,[],1));
    D44 = real(reshape(D44,[],1)); D45 = real(reshape(D45,[],1)); D46 = real(reshape(D46,[],1));
    D55 = real(reshape(D55,[],1)); D56 = real(reshape(D56,[],1));
    D66 = real(reshape(D66,[],1));
    %%  FFT transformation for second part
    fftDT11 = G2hat11.*fftChi; fftDT12 = G2hat12.*fftChi; fftDT13 = G2hat13.*fftChi; fftDT14 = G2hat14.*fftChi; fftDT15 = G2hat15.*fftChi; fftDT16 = G2hat16.*fftChi;
    fftDT22 = G2hat22.*fftChi; fftDT23 = G2hat23.*fftChi; fftDT24 = G2hat24.*fftChi; fftDT25 = G2hat25.*fftChi; fftDT26 = G2hat26.*fftChi;
    fftDT33 = G2hat33.*fftChi; fftDT34 = G2hat34.*fftChi; fftDT35 = G2hat35.*fftChi; fftDT36 = G2hat36.*fftChi;
    fftDT44 = G2hat44.*fftChi; fftDT45 = G2hat45.*fftChi; fftDT46 = G2hat46.*fftChi;
    fftDT55 = G2hat55.*fftChi; fftDT56 = G2hat56.*fftChi;
    fftDT66 = G2hat66.*fftChi;
    % Reverse back
    DT11 = ifftn(ifftshift(fftDT11)); DT12 = ifftn(ifftshift(fftDT12)); DT13 = ifftn(ifftshift(fftDT13)); DT14 = ifftn(ifftshift(fftDT14)); DT15 = ifftn(ifftshift(fftDT15)); DT16 = ifftn(ifftshift(fftDT16));
    DT22 = ifftn(ifftshift(fftDT22)); DT23 = ifftn(ifftshift(fftDT23)); DT24 = ifftn(ifftshift(fftDT24)); DT25 = ifftn(ifftshift(fftDT25)); DT26 = ifftn(ifftshift(fftDT26));
    DT33 = ifftn(ifftshift(fftDT33)); DT34 = ifftn(ifftshift(fftDT34)); DT35 = ifftn(ifftshift(fftDT35)); DT36 = ifftn(ifftshift(fftDT36));
    DT44 = ifftn(ifftshift(fftDT44)); DT45 = ifftn(ifftshift(fftDT45)); DT46 = ifftn(ifftshift(fftDT46));
    DT55 = ifftn(ifftshift(fftDT55)); DT56 = ifftn(ifftshift(fftDT56));
    DT66 = ifftn(ifftshift(fftDT66));
    % If m, n, l are even, D_ij (i ~= j) are complex numbers, choose the real part
    DT11 = real(reshape(DT11,[],1)); DT12 = real(reshape(DT12,[],1)); DT13 = real(reshape(DT13,[],1)); DT14 = real(reshape(DT14,[],1)); DT15 = real(reshape(DT15,[],1)); DT16 = real(reshape(DT16,[],1));
    DT22 = real(reshape(DT22,[],1)); DT23 = real(reshape(DT23,[],1)); DT24 = real(reshape(DT24,[],1)); DT25 = real(reshape(DT25,[],1)); DT26 = real(reshape(DT26,[],1));
    DT33 = real(reshape(DT33,[],1)); DT34 = real(reshape(DT34,[],1)); DT35 = real(reshape(DT35,[],1)); DT36 = real(reshape(DT36,[],1));
    DT44 = real(reshape(DT44,[],1)); DT45 = real(reshape(DT45,[],1)); DT46 = real(reshape(DT46,[],1));
    DT55 = real(reshape(DT55,[],1)); DT56 = real(reshape(DT56,[],1));
    DT66 = real(reshape(DT66,[],1));
    for j = 1:Ncluster
        %% index of cluster
        Cindex =list&(idx==j);
        %% First part
        Done11(j,i)= sum(D11(Cindex))/sum(Cindex);Done12(j,i)= sum(D12(Cindex))/sum(Cindex);Done13(j,i)= sum(D13(Cindex))/sum(Cindex);Done14(j,i)= sum(D14(Cindex))/sum(Cindex);Done15(j,i)= sum(D15(Cindex))/sum(Cindex);Done16(j,i)= sum(D16(Cindex))/sum(Cindex);
        Done22(j,i)= sum(D22(Cindex))/sum(Cindex);Done23(j,i)= sum(D23(Cindex))/sum(Cindex);Done24(j,i)= sum(D24(Cindex))/sum(Cindex);Done25(j,i)= sum(D25(Cindex))/sum(Cindex);Done26(j,i)= sum(D26(Cindex))/sum(Cindex);
        Done33(j,i)= sum(D33(Cindex))/sum(Cindex);Done34(j,i)= sum(D34(Cindex))/sum(Cindex);Done35(j,i)= sum(D35(Cindex))/sum(Cindex);Done36(j,i)= sum(D36(Cindex))/sum(Cindex);
        Done44(j,i)= sum(D44(Cindex))/sum(Cindex);Done45(j,i)= sum(D45(Cindex))/sum(Cindex);Done46(j,i)= sum(D46(Cindex))/sum(Cindex);
        Done55(j,i)= sum(D55(Cindex))/sum(Cindex);Done56(j,i)= sum(D56(Cindex))/sum(Cindex);
        Done66(j,i)= sum(D66(Cindex))/sum(Cindex);
        %% Second part
        Dtwo11(j,i)= sum(DT11(Cindex))/sum(Cindex);Dtwo12(j,i)= sum(DT12(Cindex))/sum(Cindex);Dtwo13(j,i)= sum(DT13(Cindex))/sum(Cindex);Dtwo14(j,i)= sum(DT14(Cindex))/sum(Cindex);Dtwo15(j,i)= sum(DT15(Cindex))/sum(Cindex);Dtwo16(j,i)= sum(DT16(Cindex))/sum(Cindex);
        Dtwo22(j,i)= sum(DT22(Cindex))/sum(Cindex);Dtwo23(j,i)= sum(DT23(Cindex))/sum(Cindex);Dtwo24(j,i)= sum(DT24(Cindex))/sum(Cindex);Dtwo25(j,i)= sum(DT25(Cindex))/sum(Cindex);Dtwo26(j,i)= sum(DT26(Cindex))/sum(Cindex);
        Dtwo33(j,i)= sum(DT33(Cindex))/sum(Cindex);Dtwo34(j,i)= sum(DT34(Cindex))/sum(Cindex);Dtwo35(j,i)= sum(DT35(Cindex))/sum(Cindex);Dtwo36(j,i)= sum(DT36(Cindex))/sum(Cindex);
        Dtwo44(j,i)= sum(DT44(Cindex))/sum(Cindex);Dtwo45(j,i)= sum(DT45(Cindex))/sum(Cindex);Dtwo46(j,i)= sum(DT46(Cindex))/sum(Cindex);
        Dtwo55(j,i)= sum(DT55(Cindex))/sum(Cindex);Dtwo56(j,i)= sum(DT56(Cindex))/sum(Cindex);
        Dtwo66(j,i)= sum(DT66(Cindex))/sum(Cindex);
    end
end
toc
%% Postprocessing for 1st and 2nd Part of interaction tensor (store into D1 and D2)
for  i = 1:Ncluster
    for j = 1:Ncluster
        %% The first part
        D1(6*i-5,6*j-5) = Done11(i,j);D1(6*i-5,6*j-4) = Done12(i,j);D1(6*i-5,6*j-3) = Done13(i,j);D1(6*i-5,6*j-2) = Done14(i,j);D1(6*i-5,6*j-1) = Done15(i,j);D1(6*i-5,6*j) = Done16(i,j);
        D1(6*i-4,6*j-5) = Done12(i,j);D1(6*i-4,6*j-4) = Done22(i,j);D1(6*i-4,6*j-3) = Done23(i,j);D1(6*i-4,6*j-2) = Done24(i,j);D1(6*i-4,6*j-1) = Done25(i,j);D1(6*i-4,6*j) = Done26(i,j);
        D1(6*i-3,6*j-5) = Done13(i,j);D1(6*i-3,6*j-4) = Done23(i,j);D1(6*i-3,6*j-3) = Done33(i,j);D1(6*i-3,6*j-2) = Done34(i,j);D1(6*i-3,6*j-1) = Done35(i,j);D1(6*i-3,6*j) = Done36(i,j);
        D1(6*i-2,6*j-5) = Done14(i,j);D1(6*i-2,6*j-4) = Done24(i,j);D1(6*i-2,6*j-3) = Done34(i,j);D1(6*i-2,6*j-2) = Done44(i,j);D1(6*i-2,6*j-1) = Done45(i,j);D1(6*i-2,6*j) = Done46(i,j);
        D1(6*i-1,6*j-5) = Done15(i,j);D1(6*i-1,6*j-4) = Done25(i,j);D1(6*i-1,6*j-3) = Done35(i,j);D1(6*i-1,6*j-2) = Done45(i,j);D1(6*i-1,6*j-1) = Done55(i,j);D1(6*i-1,6*j) = Done56(i,j);
        D1(6*i,6*j-5)    = Done16(i,j);D1(6*i,6*j-4)    = Done26(i,j);D1(6*i,6*j-3)    = Done36(i,j);D1(6*i,6*j-2)    = Done46(i,j);D1(6*i,6*j-1)    = Done56(i,j);D1(6*i,6*j)    = Done66(i,j);
        %% The second part
        D2(6*i-5,6*j-5) = Dtwo11(i,j);D2(6*i-5,6*j-4) = Dtwo12(i,j);D2(6*i-5,6*j-3) = Dtwo13(i,j);D2(6*i-5,6*j-2) = Dtwo14(i,j);D2(6*i-5,6*j-1) = Dtwo15(i,j);D2(6*i-5,6*j) = Dtwo16(i,j);
        D2(6*i-4,6*j-5) = Dtwo12(i,j);D2(6*i-4,6*j-4) = Dtwo22(i,j);D2(6*i-4,6*j-3) = Dtwo23(i,j);D2(6*i-4,6*j-2) = Dtwo24(i,j);D2(6*i-4,6*j-1) = Dtwo25(i,j);D2(6*i-4,6*j) = Dtwo26(i,j);
        D2(6*i-3,6*j-5) = Dtwo13(i,j);D2(6*i-3,6*j-4) = Dtwo23(i,j);D2(6*i-3,6*j-3) = Dtwo33(i,j);D2(6*i-3,6*j-2) = Dtwo34(i,j);D2(6*i-3,6*j-1) = Dtwo35(i,j);D2(6*i-3,6*j) = Dtwo36(i,j);
        D2(6*i-2,6*j-5) = Dtwo14(i,j);D2(6*i-2,6*j-4) = Dtwo24(i,j);D2(6*i-2,6*j-3) = Dtwo34(i,j);D2(6*i-2,6*j-2) = Dtwo44(i,j);D2(6*i-2,6*j-1) = Dtwo45(i,j);D2(6*i-2,6*j) = Dtwo46(i,j);
        D2(6*i-1,6*j-5) = Dtwo15(i,j);D2(6*i-1,6*j-4) = Dtwo25(i,j);D2(6*i-1,6*j-3) = Dtwo35(i,j);D2(6*i-1,6*j-2) = Dtwo45(i,j);D2(6*i-1,6*j-1) = Dtwo55(i,j);D2(6*i-1,6*j) = Dtwo56(i,j);
        D2(6*i,6*j-5)    = Dtwo16(i,j);D2(6*i,6*j-4)    = Dtwo26(i,j);D2(6*i,6*j-3)    = Dtwo36(i,j);D2(6*i,6*j-2)    = Dtwo46(i,j);D2(6*i,6*j-1)    = Dtwo56(i,j);D2(6*i,6*j)    = Dtwo66(i,j);
    end
end