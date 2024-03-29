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
clc
close all
clear variables
%% Load File 
ucp=zeros(1,3);
ucp(1)=128;% Length 
ucp(2)=128;% Width 
ucp(3)=32;% Height 
load('YarnX.mat'); % The element ID of  YarnX (YarnX represents yarn centerline is streching along X Axis)
load('YarnY.mat'); %  The element ID of  YarnY (YarnY represents yarn centerline is streching along Y Axis)
load('Matrix.mat'); % The element (Voxel) ID of matrix material 
load('ElementConnective.mat'); % Element connective matrix 
load('NodeCoordinate.mat'); % Nodal coordinate of C3D8R voxel/element
load('ORI.mat') % Orientation file of yarn material 
% Load strain concentration tensor of each C3D8R voxel/element
load('A1.mat'); load('A2.mat'); load('A3.mat'); load('A4.mat'); load('A5.mat'); load('A6.mat');
gvector = [A1 A2 A3 A4 A5 A6];
% Later the YarnX is decomposed into Yarnxpos and Yarnxneg
% Later the YarnY is decomposed into Yarnypos and Yarnyneg
%% Define the vectors:Upper and Lower Coordinate,center of voxel mesh
flag=zeros(length(ElementConnective),1);
flag(YarnX)=1;
flag(YarnY)=3;
flag(Matrix)=0;
ULcoordinate=zeros(length(flag),6);%% Upper and lower coordinate
Middlecoordinate=zeros(length(flag),3);% center coordinate of each voxel
for i=1:length(flag)
    N=zeros(1,8);
    for j=1:8
        N(j)=ElementConnective(i,j+1);
    end
    ULcoordinate(i,1)=max([NodeCoordinate(N(1),2) NodeCoordinate(N(2),2) NodeCoordinate(N(3),2) NodeCoordinate(N(4),2) NodeCoordinate(N(5),2) NodeCoordinate(N(6),2) NodeCoordinate(N(7),2) NodeCoordinate(N(8),2)]);%Upper X
    ULcoordinate(i,2)=min([NodeCoordinate(N(1),2) NodeCoordinate(N(2),2) NodeCoordinate(N(3),2) NodeCoordinate(N(4),2) NodeCoordinate(N(5),2) NodeCoordinate(N(6),2) NodeCoordinate(N(7),2) NodeCoordinate(N(8),2)]);%Lower X
    ULcoordinate(i,3)=max([NodeCoordinate(N(1),3) NodeCoordinate(N(2),3) NodeCoordinate(N(3),3) NodeCoordinate(N(4),3) NodeCoordinate(N(5),3) NodeCoordinate(N(6),3) NodeCoordinate(N(7),3) NodeCoordinate(N(8),3)]);%Upper Y
    ULcoordinate(i,4)=min([NodeCoordinate(N(1),3) NodeCoordinate(N(2),3) NodeCoordinate(N(3),3) NodeCoordinate(N(4),3) NodeCoordinate(N(5),3) NodeCoordinate(N(6),3) NodeCoordinate(N(7),3) NodeCoordinate(N(8),3)]);%Lower Y
    ULcoordinate(i,5)=max([NodeCoordinate(N(1),4) NodeCoordinate(N(2),4) NodeCoordinate(N(3),4) NodeCoordinate(N(4),4) NodeCoordinate(N(5),4) NodeCoordinate(N(6),4) NodeCoordinate(N(7),4) NodeCoordinate(N(8),4)]);%Upper Z
    ULcoordinate(i,6)=min([NodeCoordinate(N(1),4) NodeCoordinate(N(2),4) NodeCoordinate(N(3),4) NodeCoordinate(N(4),4) NodeCoordinate(N(5),4) NodeCoordinate(N(6),4) NodeCoordinate(N(7),4) NodeCoordinate(N(8),4)]);%Lower Z
    Middlecoordinate(i,1)=0.5*(ULcoordinate(i,1)+ULcoordinate(i,2));%X
    Middlecoordinate(i,2)=0.5*(ULcoordinate(i,3)+ULcoordinate(i,4));%Y
    Middlecoordinate(i,3)=0.5*(ULcoordinate(i,5)+ULcoordinate(i,6));%Z
end
%% Upper and Lower Coordinate
diagram=[flag,ULcoordinate,ORI];
negx=0;
negy=0;
[len,width]=size(diagram);
for j=1:length(flag)
    if((diagram(j,1)==1)&&(diagram(j,width)<0))
        diagram(j,1)=2;
        negx=negx+1;
    end
    if((diagram(j,1)==3)&&(diagram(j,width)<0))
        diagram(j,1)=4;
        negy=negy+1;
    end
end
Diagram=[diagram(:,1),Middlecoordinate,ORI(:,3)];
%% About the clustering
%% Identify the ID Number of Gauss Point 
% Fill Yarn and Warp yarn
Idmat=find(Diagram(:,1)==0);
Idyarnxpos=find(Diagram(:,1)==1);
Idyarnxneg=find(Diagram(:,1)==2);
Idyarnypos=find(Diagram(:,1)==3);
Idyarnyneg=find(Diagram(:,1)==4);
Nmat=length(Idmat);
Nyarnxpos=length(Idyarnxpos);
Nyarnxneg=length(Idyarnxneg);
Nyarnypos=length(Idyarnypos);
Nyarnyneg=length(Idyarnyneg);
%% Predefine the Number of the pattern (Defined by the users)
Ncluster_m=64; % Purely based on strain concentration tensor 
Ncluster_mU=16; % Based on undulation angle 
Ncluster_mA=4; % Based on strain concentration tensor 
% Define the database which is need to be plot in the next section 
% eg. Strain concentration tensor and undulation angle
Ayarnxpos=zeros(Ncluster_mU*Ncluster_mA,36);
Ayarnxneg=zeros(Ncluster_mU*Ncluster_mA,36);
Ayarnypos=zeros(Ncluster_mU*Ncluster_mA,36);
Ayarnyneg=zeros(Ncluster_mU*Ncluster_mA,36);
%equal number of the pattern between fill and warp
Ncluster_iyarnxpos=Ncluster_mA*Ncluster_mU;
Ncluster_iyarnxneg=Ncluster_mA*Ncluster_mU;
Ncluster_iyarnypos=Ncluster_mA*Ncluster_mU;
Ncluster_iyarnyneg=Ncluster_mA*Ncluster_mU;
TotNc=Ncluster_m+4*(Ncluster_mA*Ncluster_mU);
idx = zeros(length(gvector),1); %ID number of the Gauss point in total in case of clarifying the cluster ID
%% Extract the Undulation Angle in Warp and Fill
Uangle_yarnxpos=Diagram(Idyarnxpos,:);
Uangle_yarnxneg=Diagram(Idyarnxneg,:);
Uangle_yarnypos=Diagram(Idyarnypos,:);
Uangle_yarnyneg=Diagram(Idyarnyneg,:);
%% k-mean clustering OF Matrix based on the strain conectration tensor (gvector)
disp('Start kmeans clustering of Matrix ...');
gvector_m = gvector(Idmat,:);
[idx_m,Amatrix] = kmeans(gvector_m,Ncluster_m,'MaxIter',10000);
% The total index of the patterns
idx(Idmat) = idx_m;
ResultMat=[Diagram(Idmat,:),idx(Idmat)];
figure(1)
scatter3(ResultMat(:,2),ResultMat(:,3),ResultMat(:,4),16,ResultMat(:,6),'filled','square')
title('Clustering Based on Strain Concentration Tensor (Matrix)','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis([0 ucp(1) 0 ucp(2) 0 ucp(3)])
xlabel('x','FontWeight','bold','FontName','Times New Roman','FontSize',12);
ylabel('y','FontWeight','bold','FontName','Times New Roman','FontSize',12);
zlabel('z','FontWeight','bold','FontName','Times New Roman','FontSize',12);
box on;grid on;colorbar;h=colorbar;
ylabel(h,'Cluster ID','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis equal
%% k-mean clustering OF YarnXpos
disp('Start kmeans clustering of YarnXpos...');
%Clustering Based on Undulation Angle
Ncluster_mUyarnxpos=Ncluster_mU;
%In each Undulation Cluster, cluster with A
Ncluster_mAyarnxpos=Ncluster_mA;
%Clustering with Undulation Angle and return the cluster label for each%Gauss point
gvector_iyarnxUpos = Uangle_yarnxpos; % Flag,position,Undulation angle,and A tensor,
[idx_yarnxUpos,uyarnxpos] = kmeans(gvector_iyarnxUpos(:,5),Ncluster_mUyarnxpos,'MaxIter',10000);
%Preparation for the next clustering step based on A tensor,
gvector_yarnxApos = [gvector_iyarnxUpos,gvector(Idyarnxpos,:),idx_yarnxUpos,Idyarnxpos,zeros(Nyarnxpos,1)];
for i=1:Ncluster_mUyarnxpos
    %Find the A tensor from cluster ID from 1:Ncluster_mUw1
    temgvecxApos=gvector_yarnxApos(gvector_yarnxApos(:,42)==i,:);
    [idx_yarnxApos,Ayarnxpos((i-1)*Ncluster_mAyarnxpos+1:i*Ncluster_mAyarnxpos,:)]=kmeans(temgvecxApos(:,6:41), Ncluster_mAyarnxpos,'MaxIter',5000);
    gvector_yarnxApos(gvector_yarnxApos(:,42)==i,44)=Ncluster_m+idx_yarnxApos+Ncluster_mAyarnxpos*(i-1);
end
figure(2)
scatter3(gvector_yarnxApos(:,2),gvector_yarnxApos(:,3),gvector_yarnxApos(:,4),16,gvector_yarnxApos(:,44),'filled')
title('Clustering Based on Strain Concentration Tensor (YarnXpos)','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis([0 ucp(1) 0 ucp(2) 0 ucp(3)])
xlabel('x','FontWeight','bold','FontName','Times New Roman','FontSize',12);
ylabel('y','FontWeight','bold','FontName','Times New Roman','FontSize',12);
zlabel('z','FontWeight','bold','FontName','Times New Roman','FontSize',12);
box on;grid on;colorbar;h=colorbar;
ylabel(h,'Cluster ID','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis equal
%% k-mean clustering OF YarnXneg
disp('Start kmeans clustering of YarnXneg...');
%Clustering Based on Undulation Angle
Ncluster_mUyarnxneg=Ncluster_mU;
%In each Undulation Cluster, cluster with A
Ncluster_mAyarnxneg=Ncluster_mA;
%Clustering with Undulation Angle and return the cluster label for each%Gauss point
gvector_iyarnxUneg = Uangle_yarnxneg; % Flag,position,Undulation angle,and A tensor,
[idx_yarnxUneg,uyarnxneg] = kmeans(gvector_iyarnxUneg(:,5),Ncluster_mUyarnxneg,'MaxIter',10000);
%Preparation for the next clustering step based on A tensor,
gvector_yarnxAneg = [gvector_iyarnxUneg,gvector(Idyarnxneg,:),idx_yarnxUneg,Idyarnxneg,zeros(Nyarnxneg,1)];
for i=1:Ncluster_mUyarnxneg
    %Find the A tensor from cluster ID from 1:Ncluster_mUw1
    temgvecxAneg=gvector_yarnxAneg(gvector_yarnxAneg(:,42)==i,:);
    [idx_yarnxAneg,Ayarnxneg((i-1)*Ncluster_mAyarnxneg+1:i*Ncluster_mAyarnxneg,:)]=kmeans(temgvecxAneg(:,6:41), Ncluster_mAyarnxneg,'MaxIter',5000);
    gvector_yarnxAneg(gvector_yarnxAneg(:,42)==i,44)=Ncluster_m+Ncluster_iyarnxpos+idx_yarnxAneg+Ncluster_mAyarnxneg*(i-1);
end
figure(3)
scatter3(gvector_yarnxAneg(:,2),gvector_yarnxAneg(:,3),gvector_yarnxAneg(:,4),16,gvector_yarnxAneg(:,44),'filled')
title('Clustering Based on Strain Concentration Tensor (YarnXneg)','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis([0 ucp(1) 0 ucp(2) 0 ucp(3)])
xlabel('x','FontWeight','bold','FontName','Times New Roman','FontSize',12);
ylabel('y','FontWeight','bold','FontName','Times New Roman','FontSize',12);
zlabel('z','FontWeight','bold','FontName','Times New Roman','FontSize',12);
box on;grid on;colorbar;h=colorbar;
ylabel(h,'Cluster ID','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis equal
%% k-mean clustering OF YarnYpos
disp('Start kmeans clustering of YarnYpos...');
%Clustering Based on Undulation Angle
Ncluster_mUyarnypos=Ncluster_mU;
%In each Undulation Cluster, cluster with A
Ncluster_mAyarnypos=Ncluster_mA;
%Clustering with Undulation Angle and return the cluster label for each%Gauss point
gvector_iyarnyUpos = Uangle_yarnypos; % Flag,position,Undulation angle,and A tensor,
[idx_yarnyUpos,uyarnypos] = kmeans(gvector_iyarnyUpos(:,5),Ncluster_mUyarnypos,'MaxIter',10000);
%Preparation for the next clustering step based on A tensor,
gvector_yarnyApos = [gvector_iyarnyUpos,gvector(Idyarnypos,:),idx_yarnyUpos,Idyarnypos,zeros(Nyarnypos,1)];
for i=1:Ncluster_mUyarnypos
    %Find the A tensor from cluster ID from 1:Ncluster_mUw1
    temgvecyApos=gvector_yarnyApos(gvector_yarnyApos(:,42)==i,:);
    [idx_yarnyApos,Ayarnypos((i-1)*Ncluster_mAyarnypos+1:i*Ncluster_mAyarnypos,:)]=kmeans(temgvecyApos(:,6:41), Ncluster_mAyarnypos,'MaxIter',5000);
    gvector_yarnyApos(gvector_yarnyApos(:,42)==i,44)=Ncluster_m+Ncluster_iyarnxpos+Ncluster_iyarnxneg+idx_yarnyApos+Ncluster_mAyarnypos*(i-1);
end
figure(4)
scatter3(gvector_yarnyApos(:,2),gvector_yarnyApos(:,3),gvector_yarnyApos(:,4),16,gvector_yarnyApos(:,44),'filled')
title('Clustering Based on Strain Concentration Tensor (YarnYpos)','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis([0 ucp(1) 0 ucp(2) 0 ucp(3)])
xlabel('x','FontWeight','bold','FontName','Times New Roman','FontSize',12);
ylabel('y','FontWeight','bold','FontName','Times New Roman','FontSize',12);
zlabel('z','FontWeight','bold','FontName','Times New Roman','FontSize',12);
box on;grid on;colorbar;h=colorbar;
ylabel(h,'Cluster ID','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis equal
%% k-mean clustering OF YarnYneg
disp('Start kmeans clustering of YarnYneg...');
%Clustering Based on Undulation Angle
Ncluster_mUyarnyneg=Ncluster_mU;
%In each Undulation Cluster, cluster with A
Ncluster_mAyarnyneg=Ncluster_mA;
%Clustering with Undulation Angle and return the cluster label for each%Gauss point
gvector_iyarnyUneg = Uangle_yarnyneg; % Flag,position,Undulation angle,and A tensor,
[idx_yarnyUneg,uyarnyneg] = kmeans(gvector_iyarnyUneg(:,5),Ncluster_mUyarnyneg,'MaxIter',10000);
%Preparation for the next clustering step based on A tensor,
gvector_yarnyAneg = [gvector_iyarnyUneg,gvector(Idyarnyneg,:),idx_yarnyUneg,Idyarnyneg,zeros(Nyarnyneg,1)];
for i=1:Ncluster_mUyarnyneg
    %Find the A tensor from cluster ID from 1:Ncluster_mUw1
    temgvecyAneg=gvector_yarnyAneg(gvector_yarnyAneg(:,42)==i,:);
    [idx_yarnyAneg,Ayarnyneg((i-1)*Ncluster_mAyarnyneg+1:i*Ncluster_mAyarnyneg,:)]=kmeans(temgvecyAneg(:,6:41), Ncluster_mAyarnyneg,'MaxIter',5000);
    gvector_yarnyAneg(gvector_yarnyAneg(:,42)==i,44)=Ncluster_m+Ncluster_iyarnxpos+Ncluster_iyarnxneg+Ncluster_iyarnypos+idx_yarnyAneg+Ncluster_mAyarnyneg*(i-1);
end
figure(5)
scatter3(gvector_yarnyAneg(:,2),gvector_yarnyAneg(:,3),gvector_yarnyAneg(:,4),16,gvector_yarnyAneg(:,44),'filled')
title('Clustering Based on Strain Concentration Tensor (YarnYneg)','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis([0 ucp(1) 0 ucp(2) 0 ucp(3)])
xlabel('x','FontWeight','bold','FontName','Times New Roman','FontSize',12);
ylabel('y','FontWeight','bold','FontName','Times New Roman','FontSize',12);
zlabel('z','FontWeight','bold','FontName','Times New Roman','FontSize',12);
box on;grid on;colorbar;h=colorbar;
ylabel(h,'Cluster ID','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis equal
%% Postprocessing and give the table 
clusterIDmatrix=[Idmat ResultMat(:,2:4) ResultMat(:,6)];
clusterIDYarnxpos=[gvector_yarnxApos(:,43) gvector_yarnxApos(:,2:4) gvector_yarnxApos(:,44)];
clusterIDYarnxneg=[gvector_yarnxAneg(:,43) gvector_yarnxAneg(:,2:4) gvector_yarnxAneg(:,44)];
clusterIDYarnypos=[gvector_yarnyApos(:,43) gvector_yarnyApos(:,2:4) gvector_yarnyApos(:,44)];
clusterIDYarnyneg=[gvector_yarnyAneg(:,43) gvector_yarnyAneg(:,2:4) gvector_yarnyAneg(:,44)];
ClusterID=[clusterIDmatrix;clusterIDYarnxpos;clusterIDYarnxneg;clusterIDYarnypos;clusterIDYarnyneg];
ClusterID=sortrows(ClusterID,1);
ClusterData=[ClusterID(:,1) ULcoordinate ClusterID(:,5)];
vfc=zeros(TotNc,2);
for i=1:TotNc
    for j=1:length(ClusterID)
        if ClusterID(j,5)==i
            vfc(i,1)=vfc(i,1)+1;
        end
    end
end
vfc(:,1)=vfc(:,1)/length(ClusterID); % Volume fraction 
vfc(:,2)=vfc(:,1)*ucp(1)*ucp(2)*ucp(3);

scatter3(ClusterID(:,2),ClusterID(:,3),ClusterID(:,4),16,ClusterID(:,5),'filled')
title('Clustering Based on Strain Concentration Tensor','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis([0 ucp(1) 0 ucp(2) 0 ucp(3)])
xlabel('x','FontWeight','bold','FontName','Times New Roman','FontSize',12);
ylabel('y','FontWeight','bold','FontName','Times New Roman','FontSize',12);
zlabel('z','FontWeight','bold','FontName','Times New Roman','FontSize',12);
box on;
grid on;
colorbar;
h=colorbar;
ylabel(h,'Cluster ID','FontWeight','bold','FontName','Times New Roman','FontSize',12)
axis equal
%% Save data for SCA online predictive stage
save('A&U-64-16-4.mat','Amatrix','uyarnxpos','uyarnxneg','uyarnypos','uyarnyneg','Ayarnxpos','Ayarnxneg','Ayarnypos','Ayarnyneg')
save('ClusterData-64-16-4.mat','ClusterData')
save('vfc-64-16-4.mat','vfc') % Volume fraction of each cluster 
