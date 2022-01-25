% B: Debris Flow Analysis


%% --- 1. LOAD YXT MATRIX AND BASIC PARAMETERS ---

%1. LOAD DATA WITH LOAD('PATH/FILENAME')

%2. SPECIFY BASIC PARAMETERS
dScanners = (1000-2*40-2*63)/1000;
%dScanners=0.84;
clrs = 'br';

AverageLines=input('How many lines are averaged in the file?'); % if no averaging has taken place, type 1
samplingrate = raw_srate/AverageLines; %HZ
scnnr=input('Which scanner is used for computation of the bed geometry?'); % 1 is upriver, 2 is downriver
date=input('Specify event date:','s'); % for plot titles
loc=input('Specifiy event location:','s'); % for plot titles

%% --- 2. DEFINE TIME SETTINGS TO ISOLATE EVENT ---

PlotEveryNthMeas = input('Plot only every nth-measurement (recommended):');

figure;
% Plot the whole
subplot(2,1,1); hold on; % Along x
for ss=1:2
    plot(xax/1000,-yxt{ss}(1:PlotEveryNthMeas:end,:)/1000,clrs(ss))
end
xlabel('Horizontal distance (m)'); ylabel('Vertical distance (m)');
subplot(2,1,2); hold on; % Along time
for ss=1:2
    plot((1:size(yxt{ss},1))/samplingrate,-yxt{ss}(:,floor(0.5*length(xax)))...
        /1000,clrs(ss)); % plot the center line
end
xlabel('Time (s)'); ylabel('Vertical distance (m)');
suptitle('Whole event');
legend('Upriver scanner', 'downriver scanner');


% ISOLATE EVENT VIA TIME INDICES:

lo_idx=input('Lower boundary [s]: '); % lower boundary in seconds
hi_idx=input('Upper boundary [s]: '); % upper boundary in seconds
tidx=lo_idx*samplingrate:hi_idx*samplingrate; % line index

lo_idx1=lo_idx-10; % Pre-event lower boundary
hi_idx1=lo_idx; % Pre-event upper boundary

%Incase a post-event bed geometry is used, uncomment the following block of
%code:
% % % 
% lo_idx1=hi_idx;
% hi_idx1=hi_idx+10;

tidx1 = lo_idx1*samplingrate:hi_idx1*samplingrate; % pre-event index

%% --- 3. COMPUTE BED GEOMETRY AND FLOW DEPTH ---

%The following block computes the flow depth. Within the specified time
%index from section 3, flow depth is computed as the difference between
%the defined bed-profile and the yxt-matrix values. The bedgeomety is an
%average of the 10 pre-event-seconds. If the bedgeometry is higher or equal
%to the flow stage, then flow depth equals zero. The FlowHeight output
%variable is in meters. 

bedgeom=nanmean(-yxt{scnnr}(tidx1,:)); %define average pre-event bed geometry

% new approach
bed=repmat(bedgeom,length(tidx),1);
FlowHeight=(-yxt{scnnr}(tidx,:)-bed)/1000;
FlowHeight(FlowHeight<0)=0;

% old aproach slow!!

% FlowHeight=zeros(length(tidx),length(xax)); 
% 
% for t=tidx;
%     for i=1:length(xax);
%           if bedgeom(i)>=-yxt{scnnr}(t,i)
%             FlowHeight(t-(min(tidx)-1),i)=0;
%                 
%           else FlowHeight(t-(min(tidx)-1),i)=(-yxt{scnnr}...
%                             (t,i)-bedgeom(i))/1000; 
%           end
%     end
% end

%% --- 4. BASIC PARAMETERS FOR FLOW VELOCITY COMPUTATION  ---

% Define size of correlation window through maximally and minimally
% expected fow velocities (in the temporal dimension). 
% vmax = input('Maximum flow velocity expected [m/s]: '); 
% vmin = input('Minimum flow velocity expected [m/s]: ');
vmax=20;
vmin=0.2;
% In order to be considered for the velocity computation, the correlation
% must exceed the following value. 
mincorrneeded = 0.1;    

% Define the lateral resolution of the corelation window. 
%reslateral = input('Define lateral resolution (mm)');
reslateral=1500;
% In cells this is xcorrwin.
xcorrwin = round(reslateral/samplingx);

% lowscidx defines the temporal range that needs to be searched for a
% possible correlation in the data from the lower scanner. 
lowscidx = 1+round(dScanners/vmax*samplingrate):...
    round(dScanners/vmin*samplingrate);

% And this results in a length of the correlation Window. 
corrwin  = 1:length(lowscidx);


% The stepsize defines the time-index steps between every correlation run.
% The idea is not to have to run a new correlation sequence for every line
% but more often than only every corrwin. 
stepsize = round(max(corrwin)/4);


%% --- 6. COMPUTE VELOCITIES ---
pad1=padarray(yxt{1},[0 xcorrwin],'both');
pad2=padarray(yxt{2},[0 xcorrwin],'both');

clear T X C sc1filt sc2filt nanflow FlowVelocity XFlowVelocity lowcorr N Quality i j
tnr = 0;

% Compute the correlation 
for i=1:stepsize:(length(tidx)-max(lowscidx))
    tnr = tnr+1;
    sc1 = pad1(tidx(i-1+corrwin),:);
    sc2 = pad2(tidx(i-1+lowscidx),:);
    

    
    taxis = 1:size(sc1,1);
  
    for k=1:size(sc1,2)
        pp = polyfit(taxis', sc1(:,k),2);
        sc1filt(:,k) = sc1(:,k) - polyval(pp,taxis');
        
        pp = polyfit(taxis', sc2(:,k),2);
        sc2filt(:,k) = sc2(:,k) - polyval(pp,taxis');
    end
    sc1filt(isnan(sc1filt)) = 0;
    sc2filt(isnan(sc2filt)) = 0;
    

     
    xnr=1;
 
    for m=xcorrwin+1:xcorrwin:(size(sc1,2)-2*xcorrwin)
        xnr = xnr+1;
        sc1fpart = sc1filt(:,m:(m+xcorrwin-1));
        sc2fpart = sc2filt(:,m-xcorrwin:(m+2*xcorrwin-1));
        

        xc2 = xcorr2(sc1fpart,sc2fpart);
        NN=sqrt(sum(dot(sc1fpart,sc1fpart))*sum(dot(sc2fpart,sc2fpart)));
        
     

        tcenterpoint = ceil(size(xc2,1)/2);
        xc2neg = xc2(1:tcenterpoint,:);
        maxcorr = max(max(xc2neg));
        C(tnr,xnr) = maxcorr;
        N(tnr,xnr)=max(NN);
        
          
        padxc2neg=padarray(xc2neg,[2 2],'both');
        maxi=max(padxc2neg(:));
        [a,b]=find(padxc2neg==maxi);
        
        if maxi==0
            qualcheck = zeros(5,5);
        else

        for i=-2:2
            for j=-2:2
                qual(i+3,j+3)=padxc2neg(a+i,b+j);
            end
        end

        for i=1:5
            for j=1:5
            if  qual(i,j)<0.7*maxi;
                qualcheck(i,j)=0;
            else qualcheck(i,j)=1;
            end
            end
        end
        end

        Quality(tnr,xnr)=sum(qualcheck(:));
        
        
        [tc, xc] = ind2sub(size(xc2neg), max(find(xc2neg==maxcorr)));
        % compute time (seconds):
        T(tnr,xnr) = (tc - tcenterpoint - (min(lowscidx)-1))/samplingrate; 
        % compute distance (mm):
        X(tnr,xnr) = (xc - ceil(size(xc2neg,2)/2))*samplingx; 

    end
end

Quality(Quality<18)=NaN; Quality(Quality>=18)=1;

Cnorm=(N/max(N(:)));

lowcorr = find(Cnorm<mincorrneeded);
Xfilt = X.*Quality; Tfilt = T.*Quality;
Xfilt(lowcorr)=NaN;
Tfilt(lowcorr)=NaN;


figure
subplot(1,2,1);
pcolor(xax([1 ((xcorrwin:xcorrwin:(size(sc1,2)-2*xcorrwin)))])/1000,...
(1:stepsize:(length(tidx)-max(lowscidx)))/samplingrate, -dScanners./Tfilt);...
shading flat;
xlabel('Horizontal distance (m)'); ylabel('Time (s)'); title('Speed (m/s)');...
colorbar

subplot(1,2,2);
pcolor(xax([1 ((xcorrwin:xcorrwin:(size(sc1,2)-2*xcorrwin)))])/1000, ...
(1:stepsize:(length(tidx)-max(lowscidx)))/samplingrate, Xfilt/1000); ...
shading flat;
xlabel('Horizontal distance (m)'); ylabel('Time (s)'); ...
title('Horizontal Shift (m)'); colorbar

FlowVelocity=(-dScanners./Tfilt); 
XFlowVelocity=((Xfilt/1000)./Tfilt);

%% --- 7. COMPUTE HYDROGRAPH ---

% The following block computes the discharge values for the selected
% event(s). The single steps are described in more detail below. 

clear NaNFlowHeight
clear NaNFlowHeightLR
clear ResizeFlowHeight

% 1. 
% Upsample flow velocities in order to match the resolution from the
% flow height matrix. The upsampling method that is used is a nearest
% neighbor algorithm. 

ResizeFlowVelocity=imresize(FlowVelocity,[size(FlowHeight,1) ...
    size(FlowHeight,2)],'nearest');

% 2. 
% Where the correlation values are below the defined threshold, no velocity
% information is available for the discharge computation. This circumstance
% is counteracted to a certain point by extrapolating flow velocities based
% on the relation with flow depth. 
% In this case, the flow velocity is computed, for every line, as the
% median flow velocity weighted with the relative flow depth (along the
% line). 

ExtrapVel=zeros(size(ResizeFlowVelocity,1),size(ResizeFlowVelocity,2));

for i=1:size(ResizeFlowVelocity,1)
    for j=1:size(ResizeFlowVelocity,2)
        if isnan(ResizeFlowVelocity(i,j))==0
            ExtrapVel(i,j)=ResizeFlowVelocity(i,j);
        else ExtrapVel(i,j)=nanmedian(ResizeFlowVelocity(i,:))*...
                (FlowHeight(i,j)/max(FlowHeight(i,:)));
        end
    end
end

% 3. For a first visual inspection and comparison between flow depth and
% flow velocity, substitute all flow depth values with NaN if the
% corresponding flow velocity cell is NaN

for i=1:size(FlowHeight,1)
    for j=1:size(FlowHeight,2)
        if isnan(ResizeFlowVelocity(i,j));
        NaNFlowHeight(i,j)=NaN;
        else NaNFlowHeight(i,j)=FlowHeight(i,j);
        end
    end
end

%4. Plot

figure
subplot(2,2,1);
pcolor(xax/1000,...
    (1:length(tidx))/samplingrate, ResizeFlowVelocity); 
shading flat;
xlabel('Horizontal distance (m)'); ylabel('Time (s)'); 
title('Speed (m/s)'); colorbar

subplot(2,2,2);
pcolor(xax/1000,...
    (1:length(tidx))/samplingrate, ...
    FlowHeight); 
shading flat;
xlabel('Horizontal distance (m)'); ylabel('Time (s)'); 
title('FlowHeight (m)'); colorbar

subplot(2,2,3);
pcolor(xax/1000,...
    (1:length(tidx))/samplingrate, NaNFlowHeight); 
shading flat;
xlabel('Horizontal distance (m)'); ylabel('Time (s)'); 
title('Flowdepth (m) filtered by Velocity'); colorbar

subplot(2,2,4);
pcolor(xax/1000,...
    (1:length(tidx))/samplingrate, ExtrapVel); 
shading flat;
xlabel('Horizontal distance (m)'); ylabel('Time (s)'); 
title('Extrapolated Flowvelocities'); colorbar

%5. Scatterplot to compare flowheight and flow velocities

ResizeFlowHeight=imresize(FlowHeight,[size(Tfilt,1) size(Tfilt,2)],...
    'bilinear');

for i=1:size(ResizeFlowHeight,1)
    for j=1:size(ResizeFlowHeight,2)
        if isnan(FlowVelocity(i,j));
        NaNFlowHeightLR(i,j)=NaN;
        else NaNFlowHeightLR(i,j)=ResizeFlowHeight(i,j);
        end
    end
end

figure;
vec_nanflow=reshape(NaNFlowHeightLR,1,numel(NaNFlowHeightLR));
vec_FlowVel=reshape(FlowVelocity,1,numel(FlowVelocity));
vec_Cnorm=reshape(Cnorm,1,numel(FlowVelocity));
scatter(vec_FlowVel,vec_nanflow,20,vec_Cnorm,'fill') 
caxis([0.1 0.4])
colorbar
title('Flow velocity vs. flow depth')
xlabel('flow velocity [ms^{-1}]')
ylabel('flow depth [m]')

% 6. Create event hydrograph for the selected event. 
%Replace -Inf Flow Velocities with NaN
ResizeFlowVelocity(ResizeFlowVelocity==-Inf)=NaN; 

Qdist=ExtrapVel.*(FlowHeight*(samplingx/1000)); %discharge for every line
Q=nansum(Qdist,2); %Sum discharges across every line
Q(Q==-Inf)=NaN; %Replace -Inf witht NaN
Q(Q==0)=NaN; % Replace 0s in Q with NaN
Qtot=nansum(Q)*(1/samplingrate); %Compute sum along whole tidx-timeline

%Use mean velocity to compute discharge curve:
Vmean=nanmean(ResizeFlowVelocity(:));
Qmean=nansum(Vmean.*(FlowHeight*(samplingx/1000)),2);
Qmeantot=nansum(Qmean)*(1/samplingrate);

% Fill Q vector with corresponding values from Qmed where Q=NaN
Qfill=Q;
Qfill(isnan(Qfill))=Qmean(isnan(Qfill)==1);
Qfilltot=sum(Qfill)*(1/samplingrate);


% 7. Plot the discharge hydrograph
%create a vector that only contains Qmed values, where no Q value is
%available for plotting hydrograph

Qmean_plot=Qmean;
Qmean_plot(isnan(Q)==0)=NaN;

figure
plot((1:length(tidx))/samplingrate,Q,'-b')
hold on
plot((1:length(tidx))/samplingrate,Qmean_plot,'-k')
title(strcat('Event hydrograph:',date,' at ',loc));
ylabel('Discharge [m^3s^{-1}]');
legend(['Total Discharge: ',num2str(round(Qfilltot)),'m^3'])

%% --- 8. WRITE OUTPUTS TO OUTPUT-STRUCTURE --

% This block allows the user to save a multitude of variables to a matlab
% structure. This makes it easy to analyse or plot the results later on,
% without having to rerun the calculations above. The whole structure is
% saved to file. 

name=input('Define filename:','s');
path=input('Define path to save file:','s');


fieldnumber=1;
%fieldnumber=length(experiment)+1;

% All the following variables are saved:

% Date= specified event date
% AveragedLines= Number of lines that were averaged
% Samplingrate= Samplingrate computed from original scan rate and number of
% averaged lines. 
% TIndx= Time index evaluated
% BedIndx= Time index used for bed geometry computation
% Bedgeometry= Vector containing the computed bed geometry
% Discharge= Discharge vector
% Qtot= Total discharge
% Qmax= Peak discharge
% Vmax= Maximum velocity
% Vmean= Average velocity
% Vmed= Median velocity
% FlowVel= Flow velocities matrix
% FlowDepth= Flow depth matrix
% Hmax= Maximum flow depth
% CorrQuality= Matrix of correlation quality values
% CorrCount= Number of cells that yielded a correlation above the specified
% correlation threshold.
% MinCorrNeeded= Specified correlation threshold.
% Comments= Field for various comments. 
% time= Time vector specified to isolate event, temporal axis 
% samplingx= spatial resolution in across bed direction (matrix
% interpolation steps)
% scnnr= Scanner number; Defines the scanner that is used to compute the
% bed geometry.
% xax= matrix x-axis (across bed axis)
% vmin,vmax, reslateral, lowscidx, corrwin, xcorrwin, stepsize= parameters
% as they were defined in section 4 for the velocity computation. 


% Define filed number as 1 (for first run) or make sure that stats file
% that is to be used is already loaded. Change file name if needed
% (shift+enter to change all the lines).



% most important values
experiment(fieldnumber).Date=date; 
experiment(fieldnumber).AveragedLines=AverageLines;
experiment(fieldnumber).Samplingrate=samplingrate; 
experiment(fieldnumber).TIndx=[lo_idx,hi_idx]; 
experiment(fieldnumber).BedIndx=[lo_idx1, hi_idx1]; 
experiment(fieldnumber).Bedgeometry=bedgeom; 
experiment(fieldnumber).Discharge=Q;
experiment(fieldnumber).Qtot=Qtot; 
experiment(fieldnumber).Qmax=max(Q);
experiment(fieldnumber).Qbest=Qfilltot;
experiment(fieldnumber).Vmax=max(ResizeFlowVelocity(:));
experiment(fieldnumber).Vmin=min(ResizeFlowVelocity(:));
experiment(fieldnumber).Vmean=nanmean(ResizeFlowVelocity(:)); 
experiment(fieldnumber).Vmed=nanmedian(ResizeFlowVelocity(:));
experiment(fieldnumber).FlowVel=ResizeFlowVelocity; 
experiment(fieldnumber).FlowDepth=FlowHeight; 
experiment(fieldnumber).Hmax=max(FlowHeight(:)); 
experiment(fieldnumber).CorrQuality=C; 
experiment(fieldnumber).CorrCount=sum(isfinite(FlowVelocity(:))); 
experiment(fieldnumber).MinCorrNeeded=mincorrneeded;
experiment(fieldnumber).comments='PRE event geometry';

%all general settings for analysis:
experiment(fieldnumber).time=tidx; 
experiment(fieldnumber).samplingx=samplingx;
experiment(fieldnumber).scnnr=scnnr; 
experiment(fieldnumber).xax=xax;

% all settings from correlation analysis
experiment(fieldnumber).vmin=vmin; 
experiment(fieldnumber).vmax=vmax; 
experiment(fieldnumber).reslateral=reslateral; 
experiment(fieldnumber).lowscidx=lowscidx; 
experiment(fieldnumber).corrwin=corrwin;
experiment(fieldnumber).xcorrwin=xcorrwin;
experiment(fieldnumber).stepsize=stepsize;
experiment(fieldnumber).dScanners=dScanners;

save(strcat(path,'/',name),'experiment');

%% --- 9. DATA ANIMATION ---


starttime=lo_idx; %Untere Zeitgrenze in Sekunden
endtime=hi_idx; %Obere Zeitgrenze in Sekunden
time=starttime*samplingrate:25:endtime*samplingrate;

writerObj = VideoWriter('final20140830sg','Uncompressed AVI');
open(writerObj);



for ii=time;
plot(xax/1000,-yxt{1}(ii,:)/1000,clrs(2),'Linewidth',1.1);
%Adjust axes:
axis([-15 15 -25 -12])
set(gca,'FontSize',18,'FontName','News Gothic MT');
%Adjust time location:
[hText]=text(-10,-6,(['time [s]=' num2str(ii/samplingrate)]));
set(hText,'FontSize',18,'FontName','News Gothic MT');
title('Guttannen Spreitgraben: 30.8.2014'); xlabel('Horizontal distance (m)'),...
ylabel('Distance to surface (m)');
frame=getframe(gcf);
writeVideo(writerObj,frame);
end
close(writerObj);

