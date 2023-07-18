function [radiusMin, angMax, distMax,angle, rApical] =  bronchial_architecture_case(LIDCcase)
%% read obj

folder  = ['C:\Users\TBD\Documents\Lung_CT_cases\LIDC-IDRI\' LIDCcase '\'];
type = 'Centerline';
[vertexCL, lineCL] = bronchial_readNetwork(folder,type);
rigidTip = 10;
minRadius = 6;
%% Plot Branches
lineNum = length(lineCL);
scrsz = get(0,'ScreenSize');
h=figure('Position',[100 10 scrsz(3)/3 scrsz(4)/1.5]);
title(folder(max(strfind(folder(1:end-1),'\'))+1:end-1))
xlabel('Right -> Left (L) [mm]'),ylabel('(P) Posterior <- Anterior [mm]'),zlabel('Feet -> Head (S) [mm]')
caz = 28;
cel = 8; %tilted view from Left lung
view(caz, cel),hold on
for ln=[1:lineNum]
    plot3(vertexCL(lineCL{ln},1),vertexCL(lineCL{ln},2),vertexCL(lineCL{ln},3),'.')
    ind = nonzeros(lineCL{ln});
    hText(ln) = text(vertexCL(lineCL{ln}(end),1),vertexCL(lineCL{ln}(end),2),vertexCL(lineCL{ln}(end),3),num2str(ln));
    %pause(1)
end
hLabel = get(gca,'ylabel');
 set(hLabel,'rotation',caz+cel,'HorizontalAlignment','right','VerticalAlignment','bottom')
 hLabel = get(gca,'xlabel');
 set(hLabel,'rotation',cel/2,'HorizontalAlignment','center','VerticalAlignment','bottom')
grid minor
axis equal
%% Display RUL branch
m = 1; % branch in right lobes
str = input('Right Upper Lobe Apical Segment Branch#','s'); %RUL branch
while isempty(str2num(str)) || str2num(str)<1 || str2num(str)>lineNum
    str = input('Right Upper Lobe Apical Segment Branch#','s');
end
rApical = str2num(str);

shortLn = find(cellfun('length',lineCL)<3); % to detect branch that is not part of main tree, i.e. usually just two points
lineCL(shortLn)=[];
display(['Removed Branch#' num2str(shortLn)]) 
lineNum = lineNum - length(shortLn);

for ln=[1:lineNum]
    minpts = min(length(lineCL{ln}),length(lineCL{rApical}));
    numPts(ln) = length(find(vecnorm(vertexCL(lineCL{ln}(1:minpts),:)'-vertexCL(lineCL{rApical}(1:minpts),:)')<1));
end

[comPts,splitLn,~] = uniquetol(numPts,10/max(numPts)); %if within 10 points, merge them
comPts = sort(comPts);

carinaPt = vertexCL(lineCL{rApical}(comPts(1)),:);
tracheaLength = sum(vecnorm(diff(vertexCL(lineCL{rApical}(1:min(comPts)),:))'));
bifurPt = vertexCL(lineCL{rApical}(comPts(2:end-1)),:);
beyondTracheaPts = vertexCL(lineCL{rApical}(min(comPts)+1:end),:);
%figure(h), hold on, 
plot3(carinaPt(1),carinaPt(2),carinaPt(3),'.','MarkerSize',28)
plot3(bifurPt(:,1),bifurPt(:,2),bifurPt(:,3),'.','MarkerSize',28)
%annotation('textarrow',[0.3,0.5],[0.6,0.5],'String','Carina');
%plot3(beyondTracheaPts(:,1),beyondTracheaPts(:,2),beyondTracheaPts(:,3),'.','MarkerSize',1)

%% Calculate Curvature every 10-mm length along Centerline
curvepts = [chopSeg(rigidTip,'e',tracheaPts); unique(beyondTracheaPts,'rows','stable')];
lastSeg = chopSeg(rigidTip,'e',curvepts);
%figure,hold on,view(caz, cel)
h2 = figure; view(caz, cel),hold on, axis equal
title({[folder(max(strfind(folder(1:end-1),'\'))+1:end-1) ' RUL Normal Vectors'],...
    'Star  Sign: Angle btw Neighboring Normal Vectors > 30deg','Red Arrow: Radius < 6mm'})
% title([folder(max(strfind(folder(1:end-1),'\'))+1:end-1) ' RUL Normal Vectors'])
plot3(curvepts(:,1),curvepts(:,2),curvepts(:,3),'.-')
gen = 1; %airway generation count
for p =1:length(curvepts)-length(lastSeg)+1
    tofit = chopSeg(rigidTip,'b',curvepts(p:end,:));
    %if sum(ismember(tofit,bifurPt,'rows'))
    endpt(p) = p+length(tofit)-1;
    [radius(p),normVec(:,p)] = circfit3(tofit,rigidTip,minRadius); %minRadius is 6mm
    if radius(p) < minRadius
        %mark the position
    end
    [azimuth, elevation] = cart2sph(normVec(1,:),normVec(2,:),normVec(3,:));
    azimuth = azimuth * 180 / pi;
    elevation = elevation * 180 / pi;
end

%[B,I] =maxk(acosd(dot(normVec(:,1:end-1),normVec(:,2:end))),5);
angNeighVec = acosd(dot(normVec(:,1:end-1),normVec(:,2:end)));
[I] = find(angNeighVec>30);
angMax = max(angNeighVec);
radiusMin = min(radius);
plot3(curvepts(endpt(I),1),curvepts(endpt(I),2),curvepts(endpt(I),3),'y*','LineWidth',3,'MarkerSize',17)

savefig(['C:\Users\C20479\Documents\SNAKE\Lung_CT_cases\LIDC-Results\' LIDCcase 'RUL.fig'])
saveas(h2,['C:\Users\C20479\Documents\SNAKE\Lung_CT_cases\LIDC-Results\' LIDCcase 'RUL.png'])
figure,plot(radius,'.','MarkerSize',10),hold on, plot(1:length(radius),ones(size(radius))*minRadius)



%% Read in Lung model
fileName = dir([folder 'Lung.obj']).name;
vertexLung= (importdata([folder fileName], ' ', 4).data);

%% Stretch endpoint to pleura
toplot=1;
figure(h)
title(folder(max(strfind(folder(1:end-1),'\'))+1:end-1))
for ln = 1:lineNum
    %Ltotal = arclength(vertexCL(lineCL{ln} ,1),vertexCL(lineCL{ln} ,2),vertexCL(lineCL{ln} ,3),'spline');
    pts = vertexCL(lineCL{ln},:);
    ptDist = vecnorm(diff(vertexCL(lineCL{ln},:),1,1)');
    totalDist = floor(sum(ptDist));
    lastSeg = pts(min(find(cumsum(ptDist)>totalDist-rigidTip)):end,:);
    
    
    lastSegDist = sum(vecnorm(diff(lastSeg,1,1)'));
    %sum(ptDist(end-length(lastSeg)+2:end))
    %ptsInterpt = fnplt(cscvn(ptsAll'));
    lastSegLine = fitline3(lastSeg')';
    if dot(lastSegLine(2,:)-lastSegLine(1,:),lastSeg(end,:)-lastSeg(1,:))<0
        lastSegLine = flipud(lastSegLine);
    end
    if mean(lastSegLine(:,3))<mean(vertexLung(:,3))
        adjustDist = mean(vertexLung(:,1))+0.1*(carinaPt(1)-mean(lastSegLine(:,1)));
    else
        adjustDist = mean(vertexLung(:,1))+0.5*(carinaPt(1)-mean(lastSegLine(:,1)));
    end
    vertexQualify = vertexLung(vecnorm(cross(vertexLung-lastSegLine(1,:),vertexLung-lastSegLine(2,:),2)')/vecnorm(diff(lastSegLine,1,1))<0.5...
                                      & (vertexLung(:,1)'-(adjustDist))*(lastSegLine(1,1)-(adjustDist))>0 ,:);
                                        % &  (vertexLung(:,1)'-carinaPt(1))*(lastSegLine(1,1)-carinaPt(1))>0,:);
    
    [minDist,indMin] = max(dot(vertexQualify-lastSegLine(2,:),repmat(lastSegLine(2,:)-lastSegLine(1,:),size(vertexQualify,1),1),2)./...
         vecnorm(vertexQualify'-lastSegLine(2,:)')'/vecnorm(lastSegLine(2,:)'-lastSegLine(1,:)'));
    carinaToPleura(ln)= sum(vecnorm(diff([pts;vertexQualify(indMin,:)],1,1)'))-tracheaLength;
    if 0% ln==37
        mArrow3(carinaPt, vertexQualify(indMin,:),'color','k')
        sprintf('Carina straight to Diaphragm is %.f [mm]\n',vecnorm((carinaPt-vertexQualify(indMin,:))'))
        projXZ = [vertexQualify(indMin,1),carinaPt(2),vertexQualify(indMin,3)];
        mArrow3(carinaPt, projXZ,'color','k')
        sprintf('Carina straight to Diaphragm projected to X-Z plane is %.f [mm]\n',vecnorm((carinaPt-projXZ)'))
    end
    if toplot
        delete(hText);    
        text(vertexQualify(indMin,1),vertexQualify(indMin,2),vertexQualify(indMin,3),num2str(ln))

        if carinaToPleura(ln)<150
            mArrow3(lastSegLine(2,:), vertexQualify(indMin,:),'color','b');
        elseif carinaToPleura(ln)<200
            mArrow3(lastSegLine(2,:), vertexQualify(indMin,:),'color','g');
        else
            mArrow3(lastSegLine(2,:), vertexQualify(indMin,:),'color','r');
        end        
    end
    %pause(0.1)
end
savefig(['C:\Users\C20479\Documents\SNAKE\Lung_CT_cases\LIDC-Results\' LIDCcase 'Tree.fig'])
saveas(h,['C:\Users\C20479\Documents\SNAKE\Lung_CT_cases\LIDC-Results\' LIDCcase 'Tree.png'])

[distMax, I] = max(carinaToPleura);
text
%% Straighten segment for Angle Measurement
lineNum=plotBranches(lineCL,vertexCL);
title(folder(max(strfind(folder(1:end-1),'\'))+1:end-1))
method =2;% 1 = linear fit subbranch; 2 = linear fit branching point neighbor 10mm segment; 3 = connect branching points


plot3(carinaPt(1),carinaPt(2),carinaPt(3),'.','MarkerSize',28)

ctr = 1;
comPts = chopSeg(rigidTip,'e',vertexCL(lineCL{splitLn(ctr)}(1:min(comPts)),:));
segLine{ctr} = fitline3(comPts');
segLineE{ctr} = fitline3(comPts');
plot3(comPts(:,1),comPts(:,2),comPts(:,3),'.')
mArrow3(segLine{ctr}(:,1),segLine{ctr}(:,2), 'stemWidth',0.3, 'facealpha',0.3);

for ln = splitLn(2:end)' 
    ctr = ctr+1;
    comPts = vertexCL(lineCL{ln}(numPts(splitLn(ctr-1)):numPts(ln)),:);
    size(comPts)
    switch method
        case 1 % linear fit subbranch
            segLine{ctr}= fitline3(comPts');
            segVec1 = diff(segLine{ctr-1},1,2);
            segVec2 = diff(segLine{ctr},1,2);
            mArrow3(segLine{ctr}(:,1),segLine{ctr}(:,2), 'stemWidth',0.3, 'facealpha',0.3);
        case 2 %linear fit branching point neighbor 10mm segment
            segLineB{ctr} = fitline3(chopSeg(rigidTip,'b',comPts)');%chopSeg(rigidTip,location,pts,dim)
            segLineE{ctr} = fitline3(chopSeg(rigidTip,'e',comPts)');%chopSeg(rigidTip,location,pts,dim)
            segVec1 = diff(segLineE{ctr-1},1,2);
            segVec2 = diff(segLineB{ctr},1,2);
            mArrow3(segLineB{ctr}(:,1),segLineB{ctr}(:,2), 'stemWidth',0.3, 'facealpha',0.3);
            mArrow3(segLineE{ctr}(:,1),segLineE{ctr}(:,2), 'stemWidth',0.3, 'facealpha',0.3);
            %mArrow3(segLineE{ctr}(:,2),segLineE{ctr}(:,2)*2-segLineE{ctr}(:,1), 'stemWidth',0.3, 'facealpha',0.3);
        case 3 %connect two branching points
            segVec1 = diff([vertexCL(lineCL{splitLn(ctr-1)}(numPts(splitLn(ctr-1))),:) ; comPts(1,:)]);
            segVec2 = diff(comPts([1 end],:));
            mArrow3(vertexCL(lineCL{splitLn(ctr-1)}(numPts(splitLn(ctr-1))),:),comPts(1,:), 'stemWidth',0.3, 'facealpha',0.3);
    end
    angle(ctr) = acosd(dot(segVec1,segVec2)/vecnorm(segVec1')/vecnorm(segVec2'));

    text(comPts(1,1),comPts(1,2),comPts(1,3),['-\angle' num2str(180-angle(ctr),'%.0f')],'Interpreter','latex','FontSize',28);
    
 end
%% save data
save(['C:\Users\C20479\Documents\SNAKE\Lung_CT_cases\LIDC-Results\' LIDCcase '.mat'],'radius','normVec','carinaToPleura','angle')
    
