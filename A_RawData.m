% A: Raw data handling and preparation


%% --- 1. BASIC PARAMETERS FOR RAWDATA PROCESSING ---

directory = input('Directory where raw data is stored: ','s');

%Define Scanner Serials: SCANNER 1 = UPRIVER, SCANNER 2 = DOWNRIVER
%RANDA: 11160041 Spreitgraben:13210487 (upriver as of 2014)
%RANDA: 11160042 Spreitgraben:12210590 (downriver as of 2014)

sc1 = input('Serial number of upriver scanner: '); 
sc2 = input('Serial number of downriver scanner: '); 

%% --- 2. READ RAWDATA ---

% Define the follwoing inputs: Filename and path (for later saving), angle
% width (centered around 90deg), Start and End file number, scan intervals
% and whether reflectance values should be read from the file or not.
name=input('Define filename:','s');
path=input('Define path to save file:','s');
awidth=input('Define angle width: ');
sf=input('Define start file number: ');
ef=input('Define end file number: ');
angleres=input('Define the scan steps (0.5 or 0.3333 degrees): ');
withrefl=input('Type 1 to include RSSI values, 0 otherwise: ');
mindist=input('Define threshold distance (mm): ');

%use function READLMSOLD [angle dist refl scan] = ReadLMSFiles(directory,
%angleresolution, with_reflection, ScannerSerials, AvgLines,MinDist,
%AngleWidth,FileNoFromTo)
[angle, dist, refl, scan] = ReadLMSFilesOld(directory, angleres, withrefl,...
    [sc1, sc2], 1, mindist, awidth,[sf ef]);

% Divide the measurements from the two scanners...
scrge{1} = find(scan(:,2) == sc1);
scrge{2} = find(scan(:,2) == sc2);

% ... and then mirror the data from the lower scanner. This is necessary
% because the lower scanner is mountend in the opposite way. 
dist(scrge{2},:) = dist(scrge{2},end:-1:1);
refl(scrge{2},:) = refl(scrge{2},end:-1:1);

% Save the extracted raw data into a variable
save(strcat(path,'/',name),'angle','dist','refl','scan','sc1','sc2','scrge')


%% --- 3. DEFINE ACROSS-BED EXTENT ---

% Visualizing the raw data is helpful in order to verify the data content
% and to decide on the extent of the x-axis that is needed to interpolate
% the raw data onto a regular grid for further processing. The x-axis needs
% to be symmetrical. Define the maximally necessary bed extent (e.g. 20m).

PlotEveryNthMeas = input('Plot only every nth-measurement (recommended):');

nlines1 = size(dist(scrge{1}(1:PlotEveryNthMeas:end),:),1);
nlines2 = size(dist(scrge{2}(1:PlotEveryNthMeas:end),:),1);


figure; 
% Plot first the "raw" data, x-axis = angle, y-axis = distance
subplot(2,1,1); hold on
plot(dist(scrge{1}(1:PlotEveryNthMeas:end),:)','b')
plot(dist(scrge{2}(1:PlotEveryNthMeas:end),:)','r')
xlabel('Winkel (# meas)'); ylabel('Distanz vom Laser (mm)'); 
title('Laser Rohdaten');


subplot(2,1,2); hold on
[xx1, yy1] = pol2cart(ones(nlines1,1)*(angle/180*pi), dist(scrge{1}...
    (1:PlotEveryNthMeas:end),:));
[xx2, yy2] = pol2cart(ones(nlines2,1)*(angle/180*pi), dist(scrge{2}...
    (1:PlotEveryNthMeas:end),:));
plot(xx1'/1000,-yy1'/1000,'*b');axis equal
plot(xx2'/1000,-yy2'/1000,'*r');axis equal
xlabel('Horizontaldistanz (m)'); ylabel('Vertikaldistanz (m)'); 
title('In kartesischen Koordinaten');

% Define x-axis for interpolation. The sampling-distance along the x-axis
% is defined as 50mm.
maxextent=input('Define maximum bed extent in meters:');
samplingx=50;
xax=-maxextent*1000:samplingx:maxextent*1000;

%Check data content: Uncomment and adjust color scaling if necessary
figure;
pcolor(yy2/1000);shading flat; colorbar; %caxis([10 12])
figure;
pcolor(yy1/1000);shading flat; colorbar; %caxis([5 6])


%% --- 4. CONSTRUCT INTERPOLATED DISTANCE MATRIX ---

name=input('Define filename:','s');
path=input('Define path to save file:','s');


% now convert into x/y coordinates
nlines1 = size(dist(scrge{1}(:),:),1);
nlines2 = size(dist(scrge{2}(:),:),1);

% prepare empty matrix
yxt{1} = zeros(size(dist(scrge{1},:),1), size(xax,2));
yxt{2} = zeros(size(dist(scrge{2},:),1), size(xax,2));

for ss=1:2
    [xx, yy] = pol2cart(ones(nlines1,1)*(angle/180*pi), dist(scrge{ss}(1:1:end),:));
	
    for ii=1:size(xx,1)
        thisx = xx(ii,:);
        thisy = yy(ii,:);
        % Remove NaNs
        rmvals = find(isnan(thisx) | isnan(thisy));
        thisx(rmvals) = []; thisy(rmvals) = [];
        
        yxt{ss}(ii,:) = interp1(thisx,thisy,xax);
    end
end

% Save yxt matrix to file
save(strcat(path,'/',name),'yxt','xax','samplingx');

%% --- 5. CROP MATRIX ---

% The yxt-variable constructed can contain large amounts of unused data in
% the spatial (across bed) as well as the temporal extent (before and after
% event). This section allows the user to crop the matrix to the spatial
% and temporal extens required. The generate plot can help to identify the
% necessary extent. 

name=input('Define filename:','s');
path=input('Define path to save file:','s');
raw_srate=input('Define scanner samplingrate: ');

%First: Define new spatial and temporal boundaries based on info from plot.
clrs = 'br';
figure;
subplot(2,1,1)
for ss=1:2
    plot(xax/1000,-yxt{ss}(floor(0.5*length(yxt{ss})),:)/1000,clrs(ss));...
    axis equal
    hold on
end
subplot(2,1,2)
for ss=1:2
    plot((1:length(yxt{ss}))/raw_srate,-yxt{ss}(:,floor(0.5*length(xax)))/1000,...
    clrs(ss))
    hold on
end

leftpoint=input('Enter left value (in m (get the sign right!)): ');
rightpoint=input('Enter right value (in m): ');
newaxis=((abs(min(xax/1000))+leftpoint)/(samplingx/1000))+1:1:...
    ((abs(min(xax/1000))+rightpoint)/(samplingx/1000)+1);
startpoint=input('Enter starting point in seconds:');
endpoint=input('Enter end point in seconds: ');

%Second:
%CROP YXT MATRIX
leftedge=min(newaxis);
rightedge=max(newaxis);
start=startpoint*raw_srate;
ending=endpoint*raw_srate;

for ss=1:2
    yxt{ss}=yxt{ss}(start:ending,leftedge:rightedge);
end

%compute new xaxis from the selected boundary positions
xax=leftpoint*1000:samplingx:rightpoint*1000;

save(strcat(path,'/',name),'yxt','xax','samplingx','raw_srate');

%% --- 6. CONSTRUCT REFLECTANCE VARIABLE ---

%If RSSI (received signal strength index) data is available (parameter
%'with reflectance' must be set to 1 when the raw data is read, then this
%can be extracted into a matrix here. Raw reflectance values are saved
%as 'reflectance', normalized intensities (scaled with distance) are saved
%as 'normint' (mormalized intensity).
name=input('Define filename:','s');
path=input('Define path to save file:','s');


reflectance(:,:,1)=refl(scrge{1}(1:1:end),:);
reflectance(:,:,2)=refl(scrge{2}(1:1:end),:);
[xx1,zz1]=pol2cart(ones(nlines1,1)*(angle/180*pi), dist(scrge{1}(1:1:end),:));
[xx2,zz2]=pol2cart(ones(nlines2,1)*(angle/180*pi), dist(scrge{2}(1:1:end),:));


%scaling

ndist=floor(min(zz1(:))); %standard distance from laser in mm;
normfact=(zz1/ndist).^2;

for i=1:2
normin(:,:,i)=reflectance(:,:,i).*normfact;
end

save(strcat(path,'/',name),'reflectance','normin','xax','samplingx','raw_srate');

%% --- 7. DATA REDUCTION: AVERAGE FILES ---

% Average a certain number of lines in order to reduce the amount of data.
% The new yxt matrix is then saved to file. 

name=input('Define filename:','s');
path=input('Define path to save file:','s');

Avlines=input('Number of lines to be averaged: ');

clear avg_yxt

for ss=1:2
    for m=1:size(yxt{ss},2);
        for n=1:Avlines:(floor(length(yxt{ss})/Avlines)*Avlines),
    avg_yxt{ss}(floor(n/Avlines)+1,m)=mean(yxt{ss}(n:n+(Avlines-1),m));
        end
    end
end

clear yxt
yxt=avg_yxt;

save(strcat(path,'/',name),'yxt','xax','samplingx','Avlines','raw_srate');


%% --- 8. NOISE REDUCTION: NOISE FILTERING ---

name=input('Define filename:','s');
path=input('Define path to save file:','s');
% wrw=input('Define filtering window length (rows, time dimension)');
% wcl=input('Define filtering window widths (columns, spatial dimension)');
wrw=3;
wcl=3;
window=[wrw wcl]; 

clear medfilt_yxt diff 

for ss=1:2
    medfilt_yxt{ss}=medfilt2(yxt{ss},window); % apply median filter
    diff{ss}=yxt{ss}-medfilt_yxt{ss}; %difference to original 
    difference{ss}=diff{ss};
    diff{ss}(diff{ss}<-3*nanstd(difference{ss}(:)))=0;%set filtering threshold

    yxt_filt{ss}=medfilt_yxt{ss}+diff{ss};%add to median filtered signal
    % count number of cells that were filtered out and compute percentage
    % of original file:
    diffilt{ss}=difference{ss}-diff{ss};
    count{ss}=find(diffilt{ss}~=0);
    count{ss}=length(count{ss})/(size(yxt{ss},1)*size(yxt{ss},2))*100;
    diffilt{ss}(diffilt{ss}==0)=NaN;
    meanout{ss}=nanmean(diffilt{ss}(:));
end

 clear yxt
 yxt=yxt_filt;


save(strcat(path,'/',name),'yxt','difference','xax','samplingx',...
    'count','diffilt','meanout','raw_srate');



