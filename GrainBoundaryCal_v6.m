% Modified by 5/5/2016
% Add inclination calculation

% Modified by 5/12/2016
% Change the result, each row represents a straight segment
% Add a restrain to the inclination calculatation, exp. 11222 to 11 222
% not 112 22

% Modified by 5/17/2016
% Add the start and end points into the cc matrix

% Problems remain: 
% 1. The orientation of straight sgements should be weighted average value?
% 2. The symmetry has not been applied.  

% Modified by 9/9/2016
% To track the boundary, mlinePoints need to be passed to the calculation part. 
% Add the code to rebuild the boundary after applying the threhold.

% Modified by 9/23/2016
% Add the misorientation calculation, to the grainBoundary matrix as the angle
% between two adjacent phases. 
% grainBoundary(phase_1,phase_2,startPoint_x,startPoint_y,endPoint_x,endPoint_y,boundary length, boundary orientation, misoriantation)

% input file grainconfigend.mat and grainorientation

clear
clc

load('TranglarCrystalGrowth_1.mat')

%load('partHex.mat')
%image=grainconfigend;
image=grainconfig(:,:,end);
%grainorientation=grainori;
% load('HexP_2.mat')
% image=pImage;

mImage=[];
% scan image to find grain number, triple jucntion pixels and double junction pixels

imageSize=size(image);

% mimage is the image with addtional border(-1)
mImage=ones(size(image,2)+2,size(image,1)+2).*(-1);

bImage=mImage;% Highlight the triple jucntion pixels and double junction pixels
for n=1:size(image,1)
    mImage(n+1,:)=[-1,image(n,:),-1];
end
xmax=size(mImage,2);
ymax=size(mImage,1);
surroundPixels=[-1,-1,-1,-1];
grainNumber=[];
triplePoints=[];
doublePoints=[];
unik=0;
ti=1; %triple point index
di=1; %double point index
% eight direction around target pixel
edirections=[0,1;1,1;1,0;1,-1;0,-1;-1,-1;-1,0;-1,1]; 
% four direction around target pixel
fdirections=[1,0;0,1;0,-1;-1,0];
% diagonal direction
ddirections=[1,1;1,-1;-1,-1;-1,1];

%% scan image to find grain number, triple jucntion pixels and double junction pixels.
% The scan process is made of two steps, bulk scan and border scan.
% The corner pixels are not counted.

% scan bulk pixels, triple jucntion pixels are boundary start/end points 
for i=2:(ymax-1)
    for j=2:(xmax-1)
        % find triple jucntion pixels and double junction pixels
        % image values of surrounding pixels in four direction
        for m=1:4
            surroundPixels(m)=mImage(i+fdirections(m,1),j+fdirections(m,2));
        end
        %surroundPixels(m+1)=mImage(i,j);
        unik=unique(surroundPixels);
        % find triple junction pixels triplePoint(x,y,phase1,phase2,phase3,imagevalue)
        if size(unik,2)==3
            triplePoints(ti,1)=i;
            triplePoints(ti,2)=j;
            triplePoints(ti,3)=unik(1);
            triplePoints(ti,4)=unik(2);
            triplePoints(ti,5)=unik(3);
            ti=ti+1;
            bImage(i,j)=150;%mImage(i,j)+150; % highlight triple jucntion pixels 
        % find double junction pixels doublePoint(x,y,phase1,phase2)   
        elseif size(unik,2)==2
            doublePoints(di,1)=i;
            doublePoints(di,2)=j;
            doublePoints(di,3)=unik(1);
            doublePoints(di,4)=unik(2);
            di=di+1;
            bImage(i,j)=50;%mImage(i,j)+50; % highlight double jucntion pixels
        else
            bImage(i,j)=-1; 
        end
    end
end

% get adjacent phases and grain number based on the triple junction pixels group
adjacentPhases=[];
api=1;
for n=1:size(triplePoints,1)
    p_1=triplePoints(n,3);p_2=triplePoints(n,4);p_3=triplePoints(n,5); %p_1<p_2<p_3
    if p_1*p_2*p_3>0
    adjacentPhases(api,:)=[p_1,p_2];
    adjacentPhases(api+1,:)=[p_1,p_3]; 
    adjacentPhases(api+2,:)=[p_2,p_3]; 
    api=api+3;
    grainNumber=[grainNumber,p_1,p_2,p_3]; 
    else
        adjacentPhases(api,:)=[p_2,p_3];
        grainNumber=[grainNumber,p_2,p_3];
        api=api+1;
    end
end
grainNumber=unique(grainNumber);
% simplify the adjacent phases matrix
ap=[];
api=1;
for n=1:length(adjacentPhases)
    x_phase1=adjacentPhases(n,1);
    y_phase1=adjacentPhases(n,2);
    for m=(n+1):length(adjacentPhases)
        x_phase2=adjacentPhases(m,1);
        y_phase2=adjacentPhases(m,2);
        if [x_phase1,y_phase1]==[x_phase2,y_phase2]|[y_phase1,x_phase1]==[x_phase2,y_phase2];
            ap(api)=m;
            api=api+1;
        end
    end
end
adjacentPhases(unique(ap),:)=[];

% start/end pixels group
sePixels=triplePoints;
sePixels=sortrows(sePixels,1);
sePoint=[];
boundaryLength=[];
szap=size(adjacentPhases); %size of adjacentPhases

%record chaincode cc(phase1,phase2,chaincode)
cc=zeros(xmax,ymax);
ccr=1;
%record boundary points in the same order of cc boundaryPoint(boudary_1(x),
%boudary_1(y),...)
boundaryPoint=zeros(xmax,ymax);
bpi=1;

%% get the chaincode
for x=1:szap(1,1)
    mlinePoint=[];
    phase_1=adjacentPhases(x,1);
    phase_2=adjacentPhases(x,2);
    phase_1f=find(sePixels'==phase_1)/5;
    pffi=1;
    for n=1:length(phase_1f)
        if (rem(phase_1f(n)*0.99999,1)-0.5)>=0
        phase_1k(pffi)=ceil(phase_1f(n));
        pffi=pffi+1;
        end
    end
    pfsi=1;
    phase_2f=find(sePixels'==phase_2)/5;
    for n=1:length(phase_2f)
        if (rem(phase_2f(n)*0.99999,1)-0.5)>=0
            phase_2k(pfsi)=ceil(phase_2f(n));
            pfsi=pfsi+1;
        end
    end
% start/end pixels group
        sePoint=intersect(phase_1k,phase_2k);
        
% neglect the phases only have one pixel        
        if size(sePoint)<2
            continue
        end
        phase_1k=[];
        phase_2k=[];
% modify the start/end point group based on the image value
        mse=1;
        msePoint=[];
        for n=1:length(sePoint)
            aa=sePixels(sePoint(n),:);
            if mImage(aa(1,1),aa(1,2))==phase_1||mImage(aa(1,1),aa(1,2))==phase_2
                msePoint(mse)=sePoint(n);
                mse=mse+1;
            end
        end 
        if length(msePoint)>2
            rr=1;
            rse=[];
            for n=1:length(msePoint)
                ta=sePixels(msePoint(n),:);
                for m=(n+1):length(msePoint)
                    tb=sePixels(msePoint(m),:);
                    if norm(ta(1,1:2)-tb(1,1:2))<2
                        rse(rr)=m;
                        rr=rr+1;
                    end
                end
            end
            msePoint(unique(rse))=[];     
        end
        
        if size(msePoint,2)<2
             continue   
        end
        
    startPoint=sePixels(msePoint(1),:);
    endPoint=sePixels(msePoint(2),:);
    %get the pixels on boundary
    
    szd=size(doublePoints);
    linePoint=[];
    lpi=1; %linpoint index
    for n=1:szd
        if doublePoints(n,3)==phase_1&&doublePoints(n,4)==phase_2
        linePoint(lpi,:)=doublePoints(n,:);
        lpi=lpi+1;
        end
    end
     if size(linePoint,1)<1
         
         continue
     end
    
    % mlinePoint is the boundary, mlinePoint(x,y), start with the start point
    % surroundsp is the pixels around the start point, the intersect of
    % surroundsp and slinepoint gives the second pixels in the chain
    % mlinePoint is the boundary
   
    %pixels around start point
    for m=1:8
            surroundsp(m,1)=startPoint(1,1)+edirections(m,1);
            surroundsp(m,2)=startPoint(1,2)+edirections(m,2);
    end
    ssp=[];
    sspi=1;
    % find intercept between surroundsp and linePoint
    for n=1:8
        for m=1:size(linePoint,1)
            if surroundsp(n,1)==linePoint(m,1)&&surroundsp(n,2)==linePoint(m,2)
                ssp(sspi,:)=linePoint(m,:);
                sspi=sspi+1;
            end
        end
    end
   boundaryPhase=-2; % the image value of boundary pixels
   for n=1:size(ssp,1)
       if mImage(ssp(n,1),ssp(n,2))==mImage(startPoint(1,1),startPoint(1,2))
           boundaryPhase=mImage(startPoint(1,1),startPoint(1,2));
       end
   end
   if boundaryPhase==-2;
       boundaryPhase=mImage(ssp(n,1),ssp(n,2));
   end
   
   %slinePoint(x,y,phase1,phase2)
   slp=1;
   slinePoint=[]; 

   for n=1:size(linePoint,1)
        % the image values of boundary pixels are the same of start point 
        if mImage(linePoint(n,1),linePoint(n,2))==boundaryPhase 
            slinePoint(slp,:)=linePoint(n,:);
            slp=slp+1;
        end
   end
   % if the end point value is the same with boundary phase, then add end
   % point to the slinePoint
   if boundaryPhase==mImage(endPoint(1,1),endPoint(1,2))
       slinePoint(slp,:)=endPoint(1,1:4);
   end
   
   szslp=size(slinePoint);
   chaincode=[];
   cci=1; % index of chaincode
   mlinePoint(1,:)=startPoint(1,1:2);
   secondPoint=[];
   scp=1;
   imagevalue=[];
   del=[];
   
    %the second counted pixels
    for n=1:szslp(1,1)
        for m=1:8
            if slinePoint(n,1)==surroundsp(m,1)&&slinePoint(n,2)==surroundsp(m,2)
                secondPoint(scp,1:2)=slinePoint(n,1:2);
                secondPoint(scp,3:4)=[n,m];
                scp=scp+1;
            end
        end
    end
    % the counted pixels are resigned coordinate (-1,-1)
    
    szsp=size(secondPoint);
    if szsp(1,1)==1
        mlinePoint(2,:)=secondPoint(1,1:2);
        chaincode(cci)=secondPoint(1,4);
        cci=cci+1;
    else
        for m=1:szsp(1,1)
            for l=1:8
                srundsp(l,1)=secondPoint(m,1)+edirections(l,1);
                srundsp(l,2)=secondPoint(m,2)+edirections(l,2);
                imagevalue(l)=mImage(srundsp(1,1),srundsp(1,1));
            end
            if unique(imagevalue)~=1
                   del=m;
            end
        end
        secondPoint(del,:)=[];
        mlinePoint(2,:)=secondPoint(1,1:2);
        chaincode(cci)=secondPoint(1,4);
        cci=cci+1;
    end
    slinePoint(secondPoint(:,3),1:2)=-1;
    dir=mlinePoint(2,:);
    mlp=3;
    for n=1:(szslp(1,1)-1)
        for m=1:8
            surroundPoints(m,1)=dir(1,1)+edirections(m,1);
            surroundPoints(m,2)=dir(1,2)+edirections(m,2);
        end
        nextPoint=[];
        npi=1;
        for o=1:8
            for l=1:szslp(1,1)
                if slinePoint(l,1)==surroundPoints(o,1)&&slinePoint(l,2)==surroundPoints(o,2)
                    nextPoint(npi,1:2)=slinePoint(l,1:2);
                    nextPoint(npi,3)=l; % location in slinePoint
                    nextPoint(npi,4)=o; % location in the surround dirction
                    npi=npi+1;
                end
            end
        end
        if npi>1 % if there are two pixels surround the point
            del=[];
            snp=[];snpi=1;
            for k=1:size(nextPoint,1)
                for kk=1:8
                    snp(snpi,:)=nextPoint(k,1:2)+edirections(kk,:);
                    snpi=snpi+1;
                end
                for kk=1:8
                    kki=0;
                    for kkk=1:size(slinePoint,1)
                        if snp(kk,1:2)==slinePoint(kkk,1:2)
                            kki=kki+1;
                        end
                    end
                end
                if kki>1
                    del=[del,k];
                end
            end
            nextPoint(del)=[];
        elseif size(nextPoint,1)==0
            
                continue
        end
        mlinePoint(mlp,:)=nextPoint(1,1:2);
        dir(1,1)=nextPoint(1,1);
        dir(1,2)=nextPoint(1,2);
        chaincode(cci)=nextPoint(1,4);
        cci=cci+1; mlp=mlp+1;
        slinePoint(nextPoint(1,3),1:2)=[-1,-1];
    end
    
    % record chaincode 
    cc(ccr,1)=phase_1;
    cc(ccr,2)=phase_2;
    cc(ccr,3:4)=mlinePoint(1,:);
    cc(ccr,5:6)=mlinePoint(end,:);
    cc(ccr,7:size(chaincode,2)+6)=chaincode;
    ccr=ccr+1;
    % record boundary points, mlinePoint
    boundaryPoint(1:size(mlinePoint,1),bpi)=mlinePoint(:,1);
    boundaryPoint(1:size(mlinePoint,1),bpi+1)=mlinePoint(:,2);
    bpi=bpi+2;
end

% simplify the cc matrix
delv=[];
for n=1:size(cc,1)
    if cc(n,1)==0
        delv=n;
    end
end
cc(n:ymax,:)=[];


grainBoundary=[]; 
gbi=1;% index of grainBoundary
% a sample for single Hexgonal test
% grainorientation=[0.7,1,1.1];
for y=1:size(cc,1)
    temp=cc(y,7:end);
    chaincode=temp(find(temp));
    cc_2=chaincode;
    szcc=size(chaincode,2);
    %paste the chaincode pattern, the end number of each row is the repeat
    %time. 
    pattern=zeros(size(chaincode,2)+3); 
    pi=1; %pattern index
    ccIndex=1; % corner coordinate
    nPoint=1; % number of elements in the pattern
    cornerPoint=[];
    cpi=1;

% In case of straight line
    if size(unique(chaincode),2)==1
        %%%%%%%%
        grainBoundary(y,1:7)=[cc(y,1:6),(size(cc_2,2)+1)];
        if rem(chaincode(1),4)==0
            grainBoundary(y,8)=45;
        elseif rem(chaincode(1),4)==1
            grainBoundary(y,8)=0;
        elseif rem(chaincode(1),4)==2
            grainBoundary(y,8)=135;
        else
            grainBoundary(y,8)=90;
        end
        continue
    else
    % find pattern, based on the first different number, pattern(phase1,phase2,corner point, pattern,repeat times)
    % pattern is the repeated element in a segment
    for n=1:szcc
        largeAngle=0;
        if chaincode(n)==0
            continue
        end
        % in the case that only straigth line left
        if chaincode(1)==0&&size(unique(chaincode),2)==2
            %%%%%
            break
        end
        for m=(n+1):szcc
          % in case of the presentation of large angle, special case:'1''8'
            if abs(chaincode(m)-chaincode(n))>1&&abs(chaincode(m)-chaincode(n))~=7
                cornerPoint(cpi)=m-1;
                pattern(pi,1:2)=[cc(y,1),cc(y,2)];
                %pattern(pi,3)=m-1;
                pattern(pi,3:(m-n)+2)=chaincode(n:(m-1));
                pattern(pi,end)=1;
                chaincode(n:(m-1))=0;
                pi=pi+1;
                cpi=cpi+1;
                largeAngle=1;
                break
            end
            % pattern recorded and delete it 
            if chaincode(m)~=chaincode(n)
                if (m+1)<szcc&&chaincode(m)==chaincode(m+1)
                    pattern(pi,1:2)=[cc(y,1),cc(y,2)];
                    ccIndex=m-1;
                    nPoint=m-n;
                    pattern(pi,3:nPoint+2)=chaincode(n:(m-1));
                    pattern(pi,szcc+3)=1; %pattern repeat times
                    pi=pi+1;
                    chaincode(n:m-1)=0;
                else
                    ccIndex=m;
                    nPoint=m-n+1;
                    pattern(pi,1:2)=[cc(y,1),cc(y,2)];
                    pattern(pi,3:nPoint+2)=chaincode(n:m);
                    pattern(pi,szcc+3)=1; %pattern repeat times
                    pi=pi+1;
                    chaincode(n:m)=0;
                end
                break
            else
                ccIndex=m;
            end
        end
        %in case of reaching the end 
        if ccIndex==szcc
            break
        end
        if largeAngle==1
            continue
        end
        % compare following chaincode with pattern just recorded
        for m=(ccIndex+1):nPoint:szcc
            pa=pattern(pi-1,3:szcc);
            pt=pa(find(pa));
            if (m+nPoint-1)>szcc
                cornerPoint(cpi)=m-1;
                cpi=cpi+1;
                break
            end
            % The chaincode has same pattern, change the repeat time and
            % delete it from chaincode
            if chaincode(m:(m+nPoint)-1)==pt
                % record repeat times 
                pattern((pi-1),szcc+3)=pattern((pi-1),szcc+3)+1;
                ccIndex=m+nPoint-1;
                chaincode(m:(m+nPoint-1))=0;
            else %chaincode(m:(m+nPoint-1))~=pt
                cornerPoint(cpi)=m-1;
                cpi=cpi+1;
                break
            end
        end
    end
    end
    
% record the straight line left
    if ccIndex~=szcc
        pattern(pi,1:2)=[cc(y,1),cc(y,2)];
        pattern(pi,3:(szcc-ccIndex)+2)=chaincode((ccIndex+1):szcc);
        pattern(pi,end)=1;
    end 
    
% Simpliy the pattern matrix
    delv=[];
    delh=[];
    for n=1:size(pattern,2)
        if pattern(n,1)==0
            delv=[delv,n];
        end
        if unique(pattern(:,n))==0
            delh=[delh,n];
        end
    end
    pattern(delv,:)=[];
    pattern(:,delh)=[];
    %pattern
    
% grainBoundary(phase1, phase2, boundary length, misorientation, orientation of gb, inclination_phase1, inclination_phase2)
% add startpoint and end point to the graindoundary (phase1, phase2, st_x,st_y,en_x, en_y, boundary length,orientation of gb)

     for n=1:size(pattern,1)
        grainBoundary(gbi,1:2)=pattern(n,1:2);
        
        % in case of being out of cornerPoint matrix
        if n==1
            grainBoundary(gbi,3:4)=boundaryPoint(1,(2*y-1):2*y);
        else
            grainBoundary(gbi,3:4)=boundaryPoint(cornerPoint(n-1)+1,(2*y-1):2*y);
        end
        if n~=size(pattern,1)
            grainBoundary(gbi,5:6)=boundaryPoint(cornerPoint(n)+1,(2*y-1):2*y);
        else
            grainBoundary(gbi,5:6)=boundaryPoint(size(cc_2,2)+1,(2*y-1):2*y);
        end
        
        % grain boundary length should be the distance between start point
        % and end point
        grainBoundary(gbi,7)=sqrt((grainBoundary(n,5)-(grainBoundary(n,3)))^2+(grainBoundary(n,6)-(grainBoundary(n,4)))^2);
        %grainBoundary(gbi,8)=atan((-grainBoundary(gbi,6)+grainBoundary(gbi,4))/(grainBoundary(gbi,5)-grainBoundary(gbi,3)))*180/3.141592654;
        grainBoundary(n,8)=atan((-grainBoundary(n,5)+grainBoundary(n,3))/(grainBoundary(n,6)-grainBoundary(n,4)))*180/3.141592654;
        grainBoundary(n,8)=rem(grainBoundary(n,8)+180,180);
    gbi=gbi+1;
        
    end
end

% set up threhold 
% and record which segments are combined
%combinedCorner=[];
%ccpi=1; % combined cornerPoints index

for n=1:size(grainBoundary,1)
    if grainBoundary(n,1)==0
        continue
    end
    for m=n+1:size(grainBoundary,1)
        if grainBoundary(m,1)==grainBoundary(n,1)&&grainBoundary(n,2)==grainBoundary(m,2)
            if abs(grainBoundary(n,8)-grainBoundary(m,8))<18
                grainBoundary(n,5:6)=grainBoundary(m,5:6);
              % re-calculate the boundary length or just add them
              % grainBoundary(n,3)=grainBoundary(n,3)+grainBoundary(m,3);
                grainBoundary(n,7)=sqrt((grainBoundary(n,5)-(grainBoundary(n,3)))^2+(grainBoundary(n,6)-(grainBoundary(n,4)))^2);
                %grainBoundary(n,8)=atan((-grainBoundary(n,6)+grainBoundary(n,4))/(grainBoundary(n,5)-grainBoundary(n,3)))*180/3.141592654;
                grainBoundary(n,8)=atan((-grainBoundary(n,5)+grainBoundary(n,3))/(grainBoundary(n,6)-grainBoundary(n,4)))*180/3.14159265;
                grainBoundary(n,8)=rem(grainBoundary(n,8)+180,180);
                %combinedCorner(ccpi)=cornerPoint(m);
                %ccpi=ccpi+1;
                grainBoundary(m,:)=[0,0,0,0,0,0,0,0];
            else
                break
            end
        else
            break
        end
    end
end

% delete the 0 rows

delv=[];
for n=1:size(grainBoundary,1)
    if grainBoundary(n,1)==0
        delv=[delv,n];
    end
end
grainBoundary(delv,:)=[];


% % Misorientation Calculation 
% for n=1:size(grainBoundary,1)
%     ori_1=rem(grainorientation(grainBoundary(n,1))*180/3.14159265,60);
%     ori_2=rem(grainorientation(grainBoundary(n,2))*180/3.14159265,60);
%     grainBoundary(n,9)=vpa((ori_1-grainBoundary(n,8)),4);
%     grainBoundary(n,10)=vpa((ori_2-grainBoundary(n,8)),4);
%     grainBoundary(n,11)=rem(180/pi.*abs(grainorientation(grainBoundary(n,1))-grainorientation(grainBoundary(n,2))), 60)-30;
%     n=n+1;
% end

%grainBoundary
% hist(grainBoundary(:,11), 10);




save boundaryAnalysis.mat grainBoundary grainori boundaryPoint mImage 