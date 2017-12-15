function [avenDWI,CVnDWI,aveDWI,CVDWI]=getSNR(DWI,nb0,PAR)

DIRfig='/home/schavez/MATLAB/prgms/multi-site-NIH-AV/DTI-QAfigs-June2017/'

close all

[Nx,Ny,numvols]=size(DWI)

nd=nb0*(nb0-1)/2
chknd=0;
for n=1:nb0-1

    for n2=n+1:nb0
        sprintf('%d %d',n,n2)
        chknd=chknd+1;
        Diffndwi(:,:,chknd)=DWI(:,:,n)-DWI(:,:,n2);
    end
end

% Mojdeh: added this for cases with one B0 only
if nb0==1
    nd=1;
end
%

for i=1:nd %compare the same number of DWI diff images
    imgnum=nb0+((i-1)*2)+1
    Diffdwi(:,:,i)=DWI(:,:,imgnum)-DWI(:,:,imgnum+1);
end


ccroi=makecirc(128,65,65,30);
e=edge(ccroi,'canny');

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    for i=1:nd
        tmp=Diffndwi(:,:,i);
            Nndwi{i}=double(tmp(find(ccroi)));
    end
end 

for i=1:nd %compare the same number of DWI diff images
    tmp=Diffdwi(:,:,i);
    Ndwi{i}=double(tmp(find(ccroi)));
end


% set window leveling values for images

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    tmp=Diffndwi(:,:,1);
else
    tmp=Diffdwi(:,:,1);
end
    

maxv=max(max(tmp));
minv=min(min(tmp));
maxabs=max(maxv,abs(minv));
cmin=-maxabs*1.05;
cmax=maxabs*1.05;

h1=figure(1)
%h1.Units='inches' % this isn't compatible with the set command below in older version of matlab (jdv)
set(h1, 'Visible', 'off');
set(h1, 'Position', [1000 691 1313 643]);
%set(h1, 'Position', [10 7 20 7]);
numr=2;
numc=nd;

if nd>10,
    numr=3;
    numc=ceil(nd/2);
end

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    for i=1:nd
        subplot(numr,numc,i)
        tmp=Diffndwi(:,:,i);
        tmp(find(e))=maxv*1.1;
        imagesc(tmp,[cmin,cmax])
        title(['Diffndwi #',num2str(i)])
        axis image
    end
end

for i=1:nd %compare the same number of DWI diff images
    subplot(numr,numc,(numc*(numr-1))+i)
    tmp=Diffdwi(:,:,i);
    tmp(find(e))=maxv*1.1;
    imagesc(tmp,[cmin,cmax])
    title(['Diffdwi #',num2str(i)])
    axis image
end



Ntotndwi=[];
Ntotdwi=[];
nd
for i=1:nd
    % Mojdeh: added the if statement to exclude cases with one B0 only 
    if nb0 ~= 1
        Ntotndwi=cat(1,Ntotndwi,Nndwi{i});
    end

    Ntotdwi=cat(1,Ntotdwi,Ndwi{i});
end


% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    [nndwi{1},x]=hist(Nndwi{1},50);

    for i=2:nd
        nndwi{i}=hist(Nndwi{i},x);
    end
end


for i=1:nd
    % Mojdeh: added the if statement for cases with one B0 only
    if nb0 == 1
        [ndwi{i},x]=hist(Ndwi{i},50);
    else
        ndwi{i}=hist(Ndwi{i},x);
    end 
end

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    ntot=hist(Ntotndwi,x);
end

ntot2=hist(Ntotdwi,x);

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    for i=1:nd
        fndwi{i}=fit(x',nndwi{i}','gauss1');
    end
end

for i=1:nd
    fdwi{i}=fit(x',ndwi{i}','gauss1');
end


% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    ftot=fit(x',ntot','gauss1');
end

ftot2=fit(x',ntot2','gauss1');

h2=figure(2)
set(h2,'Visible', 'off');
set(h2, 'Position', [1000 691 1313 643]);
subplot(1,3,1)

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    for i=1:nd
        stdndwi(i)=fndwi{i}.c1
        plot(i,stdndwi(i),'ro')
        hold on
    end
end

for i=1:nd
    stddwi(i)=fdwi{i}.c1
    plot(i,stddwi(i),'bo')
    hold on
end

title('Std results for various noise maps.. (red=ndwi)')
% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    stdallndwi=stdndwi(:)
end

stdalldwi=stddwi(:)

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    stdall=[stdallndwi', stdalldwi']
else
    stdall=[stdalldwi']
end

maxv=max(stdall)
minv=min(stdall)
axis([1 max(nd,4) minv*0.9 maxv*1.1])


subplot(1,3,2)
plot(x,ntot2,'b-') %see what you get if you use dwis
hold on
plot(ftot2)
legend off
stdn2=ftot2.c1
title(['DWI std(noise)=',num2str(stdn2)])
% for i=1:nd
%     plot(x,nndwi{i},'b-')
%     hold on
%     plot(fndwi{i})
%     legend off
% end
% for i=1:4
%     plot(x,ndwi{i},'k-')
%     hold on
%     plot(fdwi{i})
%     legend off
% end

subplot(1,3,3)

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    plot(x,ntot,'b-')
    hold on
    plot(ftot)
    legend off
    stdn=ftot.c1
    title(['nDWI std(noise)=',num2str(stdn)])
end

subplot(1,3,1)
plot([1 nd],[stdn2 stdn2],'b-')

% Mojdeh: added the if statement to exclude cases with one B0 only 
if nb0 ~= 1
    plot([1 nd],[stdn stdn],'r-')
end

% Mojdeh: added the if statement for cases with one B0 only
if nb0 == 1
    stdn = stdn2;
end

SNR=DWI/stdn;
for i=1:numvols
    dwi=DWI(:,:,i);
    SNRcc(i)=mean(dwi(find(ccroi)))/stdn;
end

h3=figure(3)
set(h3, 'Visible', 'off');
set(h3, 'Position', [1000 691 1313 643]);
for i=1:(2*nb0)
    subplot(2,nb0,i)
    tmp=DWI(:,:,i);% fig1name=strcat(DIRfig,'DiffImgs-',PAR)
% print('-f1',fig1name,'-djpeg')
%
% fig2name=strcat(DIRfig,'StdPlotsHist-',PAR)
% print('-f2',fig2name,'-djpeg')
%
% fig3name=strcat(DIRfig,'SNRImgs-',PAR);
% print('-f3',fig3name,'-djpeg')
%
% fig4name=strcat(DIRfig,'SNRplots-',PAR);
% print('-f4',fig4name,'-djpeg')

    tmp(find(e))=0;
    imagesc(tmp)
    axis image
    axis off
    title(['SNR=',num2str(SNRcc(i))])
end

h4=figure(4)
set(h4, 'Visible', 'off');
set(h4, 'Position', [1000 691 1313 643]);
plot(SNRcc,'k-')
hold on
plot(SNRcc(1:nb0),'ro')
plot([nb0+1:numvols],SNRcc(nb0+1:numvols),'bo')

avenDWI=mean(SNRcc(1:nb0))
stdnDWI=std(SNRcc(1:nb0))
CVnDWI=stdnDWI*100/avenDWI

aveDWI=mean(SNRcc(nb0+1:numvols))
stdDWI=std(SNRcc(nb0+1:numvols))
CVDWI=stdDWI*100/aveDWI

title(['ave(std)nDWI=',num2str(avenDWI,'%05.2f'),'(',num2str(stdnDWI,'%05.3f'),')',' ave(std)DWI=',num2str(aveDWI,'%05.2f'),'(',num2str(stdDWI,'%05.3f'),')'])
'prepare figures...'


% fig1name=strcat(DIRfig,'DiffImgs-',PAR)
% print('-f1',fig1name,'-djpeg')
%
% fig2name=strcat(DIRfig,'StdPlotsHist-',PAR)
% print('-f2',fig2name,'-djpeg')
%
% fig3name=strcat(DIRfig,'SNRImgs-',PAR);
% print('-f3',fig3name,'-djpeg')
%
% fig4name=strcat(DIRfig,'SNRplots-',PAR);
% print('-f4',fig4name,'-djpeg')

