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
%% This is SCA online predictive solver/ fixed point iteration solver
% For elasticity, we can use just one step and abondon the newton iteration
% For generalization, the code is written in non-linear manner
% This can be easily extended to nonlinear material stage
function [stranold,sigma] = SCAmeso(stranold,strinc,MIJmeso,npmeso,ndofd,ntens,mesoinc,C_i,C_0,D,toler)
%% Define the matrix and vector 
% input stranold
% output strannew
strannew= zeros(ndofd+ntens,1);
de0 = zeros(ntens,1);
de= zeros(ndofd,1);
sigmaC0 = zeros(ndofd,1);
Resmeso = zeros(ndofd,1);
sigma =  zeros(ndofd,1);
%%  
for  j = 1:mesoinc
            strannew(1:ntens)=stranold(1:ntens)+1/mesoinc*strinc(1:ntens);
            de0 (1:ntens) = 1/mesoinc*strinc(1:ntens);
            for k = 1: npmeso
                  de (6*k-5:6*k) = de0 (1:ntens);
                  strannew(6*k+1:6*k+6)=stranold(6*k+1:6*k+6)+de(6*k-5:6*k);
            end
            mesoiter = 1;
            while true %(fixed number of iteration)
                  disp(['iter=',num2str(mesoiter)]);
                  mesoiter = mesoiter + 1;
                  for  k =1:npmeso
                        sigmaC0(6*k-5:6*k) =(C_i(1:6,6*k-5:6*k)*de(6*k-5:6*k))-(C_0(1:ntens,1:ntens)*de(6*k-5:6*k));
                  end

                  for  k = 1:npmeso
                      Resmeso(6*k-5:6*k) = de(6*k-5:6*k) - de0(1:ntens);
                      for kk = 1:npmeso
                          Resmeso(6*k-5:6*k) =Resmeso(6*k-5:6*k) +D(6*k-5:6*k,6*kk-5:6*kk)*sigmaC0(6*kk-5:6*kk);
                      end
                  end
                   Resmeso(1:ndofd) = -(MIJmeso)\Resmeso(1:ndofd);
%                   Resmeso(1:ndofd) = gmres(MIJmeso,Resmeso(1:ndofd),[],0.0001,4);
                 
                  de(1:ndofd) = de(1:ndofd) + Resmeso(1:ndofd);
                  for k = 1: npmeso
                        strannew(6*k+1:6*k+6)=stranold(6*k+1:6*k+6)+de (6*k-5:6*k);
                  end
                   Ndeltaemeso=norm(Resmeso(1:ndofd));
                   Ndemeso = norm(de(1:ndofd));
                   CVmeso =Ndeltaemeso/Ndemeso;
                   disp (num2str(CVmeso))
                   if (CVmeso<toler) % Convergence check
                        break;
                   end
            end
            
            
            stranold(1:ntens)=strannew(1:ntens);
            stranold(1+ntens:ndofd+ntens)=strannew(1+ntens:ndofd+ntens);
            for k=1:npmeso % calculate and save the stress component at current step in each cluster
                sigma(6*k-5:6*k) = C_i(1:6,6*k-5:6*k)* stranold(6*k+1:6*k+6);
            end

end
end