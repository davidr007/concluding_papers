function plot_focmech(x,y,z,strike,dip,rake,fh,focsize,fcolor)
%
% function to add a focal mechanism(s) (lower hemisphere projection)
% to an existing plot
%
% Note that it is probably not very efficients because it starts with 
% an event that has strike = 0, dip = 90 and rake = 0 and then rotates everything.  
%
% INPUTS:
% x         x posn on plot
% y         y posn on plot
% z         z posn on plot
% strike    strike of event in degrees
% dip       dip of event in degrees
% rake      rake of event in degrees
% fh        figure handle
% focsize   size of circle
% fcolor     colour for focal mechanism - black if empty or not existing
%
% David Robinbson
% November 2006

tempfocsize = 1;    % The focal mechanism is created at this size first and then re-scaled
                    % In theory tempfocsize can be any number greater than
                    % 1 although this has not been tested
                 
% sort out default the colour
if nargin == 8 | isempty(fcolor)
    fcolor ='k';
end
    
% check input values
if focsize<=0
    error('focsize must be greater than 0')
end
if rake<0 | rake >90
    error('rake must be between 0 and 90 degrees')
end
if dip<0 | dip>90
    error('dip must be between 0 and 90 degrees')
end

%convert angles to radians
rstrike = strike*pi/180;
rdip = dip*pi/180;
rrake = rake*pi/180;

% turn warning off to avoid anoying divide by zero reminders 
warning off

if isempty(fh)
    figure
elseif length(fh) ==1
    figure(fh)
else
    figure(fh(1))
    axes(fh(2))
end

% create the focal circle with center of zero
circle_x = linspace(-tempfocsize,tempfocsize,1000);
circle_y_half = sqrt(tempfocsize^2-circle_x.^2);

%% rotate for dip first - i.e. we draw the great circle (see
%http://www.uwgb.edu/dutchs/STRUCTGE/SL112EquatNet.HTM)
% This corresponds to the great circle
% Note that the formulas do not work so well as dip->90
% we fix this at 90 with a separate if block and we fix those
% that -> 90 by a vertical stretch.
if dip == 90
    dgcxvec = zeros(1,1000);
    dgcyvec_upper = linspace(0,tempfocsize,1000);
    dgcyvec_lower = -linspace(0,tempfocsize,1000);
else
    drg = 2*tan(rdip)+tempfocsize*(90-dip)/90;  % dip - radius of great circle+
    dxg = -2*tan(rdip);    % dip -x position of great circle center (note that we add x and y later)
    dyg = 0;         % dip -y position of great circle center  (note that we add x and y later)
    dxmax = (90-dip)/90*tempfocsize;
    dgcxvec0 = linspace(0,dxmax,1000);
    dgcyvec0 = sqrt(drg^2-(dgcxvec0-dxg).^2);
    
    % Vertical Stretch - So long as tempfocsize>=1 the great
    % circle fails to reach the focmech circle
    % and stretching it out is simple. The vertical stretch is zero near
    % the center and max at the top
        vstretch = linspace(tempfocsize-max(dgcyvec0),0,1000);
        dgcxvec = dgcxvec0;
        dgcyvec_upper = dgcyvec0 +vstretch;
        dgcyvec_lower = -dgcyvec0-vstretch;
end
    
% Now we rotate for rake - i.e. we draw the small circle (see
%http://www.uwgb.edu/dutchs/STRUCTGE/SL112EquatNet.HTM)

%% Note that the formulas do not work so well as dip->90
% we fix this at 90 with a separate if block and we fix those
% that -> 90 by a vertical stretch.

if rake == 0
    rgcxvec_upper = linspace(0,tempfocsize,1000);
    rgcxvec_lower = -linspace(0,tempfocsize,1000);
    rgcyvec = zeros(1,1000);    
else
    rrg = 2*tan(pi/2-rrake)+tempfocsize*(rake)/90;  % rake radius of great circle+
    rxg = 0;    % rake x position of great circle center (note that we add x and y later)
    yg = -2*tan(pi/2-rrake);         % rake y position of great circle center  (note that we add x and y later)
    ymax = rake/90*tempfocsize;
    rgcyvec0 = linspace(0,ymax,1000);
    rgcxvec0 = sqrt(rrg^2-(rgcyvec0-yg).^2);
      
    % Horizontal version of vertical stretch discussed above. 
        hstretch = linspace(tempfocsize-max(rgcxvec0),0,1000);
        rgcyvec = rgcyvec0;
        rgcxvec_upper = rgcxvec0 +hstretch;
        rgcxvec_lower = -rgcxvec0 -hstretch;
end

%% Now we can re-set everything to the size (focsize) we want
% We can't do this re-size until all the stretching is done 
% because it doesn't work when focsize<1. 
[circle_x,circle_y_half] = scale_values(circle_x,circle_y_half,focsize);
[dgcxvec,dgcyvec_upper] = scale_values(dgcxvec,dgcyvec_upper,focsize);
[junk,dgcyvec_lower] = scale_values(dgcxvec,dgcyvec_lower,focsize);
[rgcxvec_upper,rgcyvec] = scale_values(rgcxvec_upper,rgcyvec,focsize);
[rgcxvec_lower,junk] = scale_values(rgcxvec_lower,rgcyvec,focsize);

%% Now we find the intersection of the two great circles
% Hint - we know the intersection occurs in the upper section of each gc
if dip ==90 &rake ==0
    intx = 0;
    inty = 0;
   % Lets pick the polygons
    P1_part1 = [0,focsize; 0,0];
    P1_part2 = [focsize,0];
    P1_part3 = [fliplr(circle_x(circle_x>=P1_part1(1,1)&circle_x<P1_part2(end,1)))', ...
             fliplr(circle_y_half(circle_x>P1_part1(1,1)&circle_x<P1_part2(end,1)))'];
    P1 = [P1_part1; P1_part2;P1_part3]; 
    
    P2_part1 = [0 0; 0 -focsize];
    P2_part2 = [];
    P2_part4 = [-focsize,0];
    P2_part3 = [fliplr(circle_x(circle_x<=P2_part1(end,1)&circle_x>P2_part4(1,1)))', ...
        fliplr(-circle_y_half(circle_x<=P2_part1(end,1)&circle_x>P2_part4(1,1)))'];
    P2_part5 = [];
    P2 = [P2_part1; P2_part2;P2_part3;P2_part4; P2_part5];
    
    
elseif rake ==0 % the intersection is at the end of gc1_upper
    intx = dgcxvec(end);
    inty = dgcyvec_upper(end);
    
    % Lets pick the polygons
    P1_part1 = [dgcxvec(dgcxvec<intx)',dgcyvec_upper(dgcxvec<intx)'];
    P1_part2 = [intx, inty; focsize,0];
    P1_part3 = [fliplr(circle_x(circle_x>=P1_part1(1,1)&circle_x<P1_part2(end,1)))', ...
             fliplr(circle_y_half(circle_x>P1_part1(1,1)&circle_x<P1_part2(end,1)))'];
    P1 = [P1_part1; P1_part2;P1_part3];
    
    P2_part1 = [dgcxvec(dgcxvec>intx)',dgcyvec_upper(dgcxvec>intx)'];
    P2_part2 = [fliplr(dgcxvec)', fliplr(dgcyvec_lower)'];
    P2_part4 = [-focsize,0;intx,0];
    P2_part3 = [fliplr(circle_x(circle_x<=P2_part2(end,1)&circle_x>P2_part4(1,1)))', ...
        fliplr(-circle_y_half(circle_x<=P2_part2(end,1)&circle_x>P2_part4(1,1)))'];
    P2_part5 = [];
    P2 = [P2_part1; P2_part2;P2_part3;P2_part4; P2_part5];
    
elseif dip ==90 % the intersection is at the end of gc2
    intx = rgcxvec_upper(end);
    inty = rgcyvec(end);
    
    % Lets pick the polygons
    P1_part1 = [0, focsize; intx, inty];
    P1_part2 = [fliplr(rgcxvec_upper(rgcxvec_upper>intx))',fliplr(rgcyvec(rgcxvec_upper>intx))'];
    P1_part3 = [fliplr(circle_x(circle_x>=P1_part1(1,1)&circle_x<P1_part2(end,1)))', ...
             fliplr(circle_y_half(circle_x>P1_part1(1,1)&circle_x<P1_part2(end,1)))'];
    P1 = [P1_part1; P1_part2;P1_part3];
    
    P2_part1 = [0,inty; 0 -focsize];
    P2_part2 = [];
    P2_part4 = [rgcxvec_lower',rgcyvec'];
    P2_part3 = [fliplr(circle_x(circle_x<=P2_part1(end,1)&circle_x>P2_part4(1,1)))', ...
        fliplr(-circle_y_half(circle_x<=P2_part1(end,1)&circle_x>P2_part4(1,1)))'];
    P2_part5 = [fliplr(rgcxvec_upper(rgcxvec_upper<intx))',fliplr(rgcyvec(rgcxvec_upper<intx))'];
    P2 = [P2_part1; P2_part2;P2_part3;P2_part4; P2_part5];
else
    tmp1y = dgcyvec_upper;
    tmp1x = dgcxvec;
    tmp2y = fliplr(rgcyvec);
    tmp2x = fliplr(rgcxvec_upper);
    newx = linspace(min([tmp1x,tmp2x]),max([tmp1x,tmp2x]),2000);
    newy1 = interp1(tmp1x,tmp1y,newx);
    newy2 = interp1(tmp2x,tmp2y,newx);
    tmpdiff = newy1-newy2;
    ind = find(tmpdiff<0);
    %plot(newx(ind(1)),newy1(ind(1)),'ro')
%     if ~isempty(ind) & ind(1)~=1
    if ~isempty(ind)
        intx = (newx(ind(1))+newx(ind(1)-1))/2;
        inty = (newy1(ind(1))+newy1(ind(1)-1))/2;
    elseif isempty(ind)         % if it is empty then new1y is close to vertical  
        ind1 = find(~isnan(newy1));   % and new1y->NaN before tempdiff<0. Hence take last non 
        intx = newx(ind1(end));       % NaN value of 
       
        inty = 0;
        %disp('gone in')
    end
       
%     elseif isempty(ind) | ind(1) ==1 
%         intx = newx(1);
%         inty = newy1(1);
%     end
    
    % Lets pick the polygons
    P1_part1 = [dgcxvec(dgcxvec<intx)',dgcyvec_upper(dgcxvec<intx)';intx, inty;];
    P1_part2 = [fliplr(rgcxvec_upper(rgcxvec_upper>intx))',fliplr(rgcyvec(rgcxvec_upper>intx))'];
    P1_part3 = [fliplr(circle_x(circle_x>=P1_part1(1,1)&circle_x<P1_part2(end,1)))', ...
             fliplr(circle_y_half(circle_x>P1_part1(1,1)&circle_x<P1_part2(end,1)))'];
    P1 = [P1_part1; P1_part2;P1_part3];
    
    P2_part1 = [dgcxvec(dgcxvec>intx)',dgcyvec_upper(dgcxvec>intx)'];
    P2_part2 = [fliplr(dgcxvec)', fliplr(dgcyvec_lower)'];
    P2_part4 = [rgcxvec_lower',rgcyvec'];
    P2_part3 = [fliplr(circle_x(circle_x<=P2_part2(end,1)&circle_x>P2_part4(1,1)))', ...
        fliplr(-circle_y_half(circle_x<=P2_part2(end,1)&circle_x>P2_part4(1,1)))'];
    P2_part5 = [fliplr(rgcxvec_upper(rgcxvec_upper<intx))',fliplr(rgcyvec(rgcxvec_upper<intx))'];
    P2 = [P2_part1; P2_part2;P2_part3;P2_part4; P2_part5; intx inty];
    
end


%% Now we can rotate for strike
% We must take every point and rotate it around the center (i.e. in xy-
% plane) by the angle = strike. 
[sgc1xvec_upper,sgc1yvec_upper] = srotate(dgcxvec,dgcyvec_upper,rstrike); % First we rotate gc1 upper
[sgc1xvec_lower,sgc1yvec_lower] = srotate(dgcxvec,dgcyvec_lower,rstrike); % Now we rotate gc1 lower
[sgc2xvec_upper,sgc2yvec_upper] = srotate(rgcxvec_upper,rgcyvec,rstrike); % Now we can repeat for gc2 - upper
[sgc2xvec_lower,sgc2yvec_lower] = srotate(rgcxvec_lower,rgcyvec,rstrike); % Finally we do it for gc2 - lower
[sintx,sinty] = srotate(intx,inty,rstrike);
[sP1x,sP1y] = srotate(P1(:,1),P1(:,2),rstrike);
[sP2x,sP2y] = srotate(P2(:,1),P2(:,2),rstrike);
    
   
% Now we do all of the plotting
% Note that we plot it in the x-y plane(z=0) first
% And then shift it up and down as required
warning on
fmh1 = plot(x+circle_x,y+circle_y_half,'color',fcolor);
fmtemp1 = length(get(fmh1,'Xdata'));
set(fmh1,'Zdata',z*ones(fmtemp1,1))
hold on
fmh2 = plot(x+circle_x,y-circle_y_half,'color',fcolor);
fmtemp2 = length(get(fmh2,'Xdata'));
set(fmh2,'Zdata',z*ones(fmtemp2,1))
fmh3 = fill(x+sP2x,y+sP2y,fcolor)
fmtemp3 = length(get(fmh3,'Xdata'));
set(fmh3,'EdgeColor',fcolor,'Zdata',z*ones(fmtemp3,1))
fmh4 = fill(x+sP1x,y+sP1y,fcolor)
fmtemp4 = length(get(fmh4,'Xdata'));
set(fmh4,'EdgeColor',fcolor,'Zdata',z*ones(fmtemp4,1))
if isempty(fh)
    axis equal
end
% plot(x+P1_part1(:,1), y+P1_part1(:,2),'b','linewidth',4)
% plot(x+P1_part2(:,1), y+P1_part2(:,2),'g','linewidth',4)
% plot(x+P1_part3(:,1), y+P1_part3(:,2),'y','linewidth',4)
% plot(x+P2_part5(:,1), y+P2_part5(:,2),'c')

function [xnew,ynew] = srotate(xorig,yorig,rotation);
    % used to rotate for strike
    radius = sqrt(xorig.^2+yorig.^2);  % strike - radius of gc1_upper
    phi2_tmp = atan(abs(yorig)./abs(xorig)); %angle of original points from the horizontal
    phi2 = zeros(size(phi2_tmp));
    phi2(xorig>0 & yorig>=0) = phi2_tmp(xorig>0 & yorig>=0);
    phi2(xorig<=0 & yorig>0) = pi-phi2_tmp(xorig<=0 & yorig>0);
    phi2(xorig<0 & yorig<=0) = pi+phi2_tmp(xorig<0 & yorig<=0);
    phi2(xorig>=0 & yorig<0) = 2*pi-phi2_tmp(xorig>=0 & yorig<0);
    phi1 = phi2-rotation;
    xnew = radius.*cos(phi1);
    ynew = radius.*sin(phi1);
    
    
 function [xnew,ynew] = scale_values(xorig,yorig,newsize);
    % used to re-scale values to focal mechanism size of interest
    radius = sqrt(xorig.^2+yorig.^2);  % strike - radius of gc1_upper
    phi2_tmp = atan(abs(yorig)./abs(xorig)); %angle of original points from the horizontal
    phi2 = zeros(size(phi2_tmp));
    phi2(xorig>0 & yorig>=0) = phi2_tmp(xorig>0 & yorig>=0);
    phi2(xorig<=0 & yorig>0) = pi-phi2_tmp(xorig<=0 & yorig>0);
    phi2(xorig<0 & yorig<=0) = pi+phi2_tmp(xorig<0 & yorig<=0);
    phi2(xorig>=0 & yorig<0) = 2*pi-phi2_tmp(xorig>=0 & yorig<0);
    xnew = (newsize.*radius).*cos(phi2);
    ynew = (newsize.*radius).*sin(phi2);
    