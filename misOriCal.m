clear
load('boundaryAnalysis.mat')
clear pi

% graindoundary (phase1, phase2, st_x,st_y,en_x, en_y, boundary
% length(7),orientation of gb(8)
% grainorientation: grainori
%misOri;
%grainori=grainorientation;
for n=1:size(grainBoundary,1)
    
    % skip the very tine segment
    if grainBoundary(n,7)<2
        continue
    end
    boundaryOri=grainBoundary(n,8);
    ori_1=rem((rem(-grainori(grainBoundary(n,1))*180/pi,120)-60),120);
    ori_2=rem((rem(-grainori(grainBoundary(n,2))*180/pi,120)-60),120);
    
    % because the three-fold symmetry, the misorientation of two grain must
    % be smaller than 60 degree.
    ori=[ori_1,ori_2;
        ori_1+120,ori_2+120;
        ori_1-120,ori_2-120;];
    % the orientation of grainboundary is in the range of (0,180), boundaryOri
    % and boundaryOri-180 should be considered,  find the range which the boundaryOri 
    % locate in.
    
    for y=1:3
        if boundaryOri<=max(ori(y,2),ori(y,1))&&boundaryOri>=min(ori(y,2),ori(y,1))
                inclination_1=ori(y,1)-boundaryOri;  % beta_1 
                inclination_2=ori(y,2)-boundaryOri; % beta_2
        elseif (boundaryOri-180)<=max(ori(y,2),ori(y,1))&&(boundaryOri-180)>=min(ori(y,2),ori(y,1));
                inclination_1=ori(y,1)-boundaryOri+180;  % beta_1 
                inclination_2=ori(y,2)-boundaryOri+180; % beta_2
        else
            continue
        end
    end
    
    misOri(n,1:2)=grainBoundary(n,1:2); % grain_1,grain_2
    misOri(n,3)=inclination_1;
    misOri(n,4)=inclination_2;
    % symmetry parameter eta
    misOri(n,5)= abs(abs(misOri(n,3))-abs(misOri(n,4)))/(abs(misOri(n,3))+abs(misOri(n,4)));
    
    % for length distribution 
    misOri(n,6)=grainBoundary(n,7);
    
%     grainBoundary(n,9)=vpa((ori_1-grainBoundary(n,8)),4);
%     grainBoundary(n,10)=vpa((ori_2-grainBoundary(n,8)),4);
%     grainBoundary(n,11)=rem(180/pi.*abs(grainorientation(grainBoundary(n,1))-grainorientation(grainBoundary(n,2))), 60)-30;
    n=n+1;
end

delv=[];
for n=1:size(misOri,1)
    if misOri(n,3)==0&&misOri(n,4)==0
        delv=[delv,n];
    end
end
misOri(delv,:)=[];

x=1;
for n=1:size(misOri,1)
    for m=1:floor(misOri(n,6))
        distDistri(x)=misOri(n,5);
        x=x+1;
    end
end
misOri
%hist(distDistri,linspace(0.1,0.9,9))
% hist(misOri(:,5), 10)
% figure
% imagesc(mImage)