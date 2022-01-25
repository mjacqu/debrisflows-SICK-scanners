%% --- 1. Binary identification ---
s=2; %define scanner number
% 
cax=(0.5*-size(yxt{1},2)*samplingx)+(0.5*samplingx):samplingx:...
     (0.5*size(yxt{1},2)*samplingx)-(0.5*samplingx);



clear flag

for t=tidx
   
   smu=-yxt{s}(t,:)/1000;
   
   poly=polyfit(cax/1000,smu,4);

   aprox=polyval(poly,cax/1000);
   test=gradient(aprox);
   
   count=0;
   
   for i=2:length(test)
       if((test(i) * test(i-1)) < 0)
           count = count+1;
           flag(t-(min(tidx)-1))=count;
       end
        
   end
   
  if poly(1)<0
      flag(t-(min(tidx)-1))=1;
  end

end



for i=1:length(flag)
    if flag(i)==1
        flag(i)=1;
    elseif flag(i)==3
            flag(i)=2;
        end
end



% figure;
% [hAxes,hBar,hLine]=plotyy(tidx/samplingrate,flag,...
%     tidx/samplingrate,-yxt{1}(tidx,...
%     floor(0.45*size(-yxt{1},2)))/1000,'bar','plot');
% set(hLine,'Color','b')
% bar_child=get(hBar,'Children');
% % no outlines around the bars:
% set(bar_child, 'EdgeColor', 'none'); 
% %assig value in flag to a 'direct' value:
% set(bar_child,'CDataMapping','direct'); 
% %direct value is in flag -> 1=red 2= green as in colormapvector:
% set(bar_child,'CData',flag) 
% mycolor=[1 0 0;0 1 0];
% colormap(mycolor)
% %set(hAxes(1),'xlim',[0 500]);
% %set(hAxes(2),'xlim',[0 500]);
% title('Scanner 1');
% ylabel(hAxes(1),'State: {\color{green}true} or {\color{red}false}',...
%     'Interpreter','Tex');
% ylabel(hAxes(2),'Distance from Scanner [m]');
% xlabel('Time [s]');

% --- 2. Curvature quantification ---

clear cW cF

for f=1:length(flag)
clear maximum minimum
    if flag(f)==2

raw=-yxt{s}(tidx(f),:);

polynomial=polyfit(cax/1000,raw,4); % 5 coeffiecients of 4th order polynomial
firstderiv=polyder(polynomial); % coefficients of 1st derivative
scndderiv=polyder(firstderiv); % coefficients of 2nd derivative
extrema=roots(firstderiv); %extrema are where the f'(x) is equal to zero (roots command)
min_max=polyval(scndderiv,extrema); %put x-coordinates of extrema into f''(x) to find out whether they are min or max


for i=1:length(min_max)
    if min_max(i)<0
        maximum=extrema(i);
    end
end

for i=1:length(min_max)
    if min_max(i)>0
        minimum(i)=extrema(i);
    end
end



h=abs(abs(polyval(polynomial,maximum))-abs((polyval(polynomial, minimum(1))...
    +polyval(polynomial,minimum(2)))/2)); 

%calculate width as difference between the x coordinates of the minima:
width=(abs(minimum(2)-minimum(1)))*1000; 
%compute curvature factor c_W as the ratio between height and width:
cW(f)=h/width;

%copute curvature factor c_H as ratio between height and max height:
Hmax=max(FlowHeight(f,:));
cH(f)=h/(Hmax*1000);
    else
        cW(f)=NaN;
        cH(f)=NaN;
    end

    
end
% 


figure;
subplot(2,1,1)
plot(tidx/samplingrate,cW,'k')
hold on
plot(tidx/samplingrate,cH,'m')
ylabel('{\color{black}h/w};{\color{magenta}h/H}')
xlim([480 540])

subplot(2,1,2)
[hAxes,hBar,hLine]=plotyy(tidx/samplingrate,flag,...
    tidx/samplingrate,-yxt{1}(tidx,...
    floor(0.45*size(-yxt{1},2)))/1000,'bar','plot');
set(hLine,'Color','b')
bar_child=get(hBar,'Children');
% no outlines around the bars:
set(bar_child, 'EdgeColor', 'none'); 
%assig value in flag to a 'direct' value:
set(bar_child,'CDataMapping','direct'); 
%direct value is in flag -> 1=red 2= green as in colormapvector:
set(bar_child,'CData',flag) 
mycolor=[1 0 0;0 1 0];
colormap(mycolor)
set(hAxes(1),'xlim',[480 540]);
set(hAxes(2),'xlim',[480 540]);
title('Scanner 1');
ylabel(hAxes(1),'State: {\color{green}true} or {\color{red}false}',...
    'Interpreter','Tex');
ylabel(hAxes(2),'Distance from Scanner [m]');
xlabel('Time [s]');


