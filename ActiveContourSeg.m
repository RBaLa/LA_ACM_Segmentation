%Code to perform semi-automatic active contour segmentation.
%Author: Rahul Balamurugan
%Date last edited: 11 December 2020
%References:
%Ritwik Kumar (2020). Snakes: Active Contour Models 
%(https://www.mathworks.com/matlabcentral/fileexchange/28109-snakes-active-contour-models), 
%MATLAB Central File Exchange.
%
%Dirk-Jan Kroon (2020). Snake : Active Contour 
%(https://www.mathworks.com/matlabcentral/fileexchange/28149-snake-active-contour), 
%MATLAB Central File Exchange.
%
%
clear;
clc
%*************************************************************************
%Main Parameters:
alpha = 0.1; %internal weight, controls smoothness of contour
beta = 0.1; %internal weight, controls rigidity of contour
gamma = 1; %Time step
kappa = 2; %External weight, controls image force exerted on contour
delta = 0.6; %External weight, controls balloon force pushing the contour.. 
             %to expand
N = 1000; %Number of iterations
%*************************************************************************
%Read the CT image
I = imread('CT4.j2c');
im = I;
%*************************************************************************
%Preprocessing:
%Smooth the working copy of the image using Gaussian filter
siz = 5; %(Size of 1-D Gaussian Kernel)
sigma = 10; %(measure of std. deviation)
xk=-ceil(siz/2):ceil(siz/2); %(the kernel)
H = exp(-(xk.^2/(2*sigma^2)));
H = H/sum(H(:)); %(designed filter)
Hx=reshape(H,[length(H) 1]);
Hy=reshape(H,[1 length(H)]);
im=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
%Intensity threshold for better image contrast (strong edges)
im(im<100) = 0;
im(im>300) = 500;
im = im2double(im);
%*************************************************************************
%Get initial contour points:
figure
imshow(I); 
[yv,xv] = getpts; %Click to select a point, best 10-12 in a downward 
                  %triangle; double click to select last point.
P=[xv(:) yv(:)]';
close()
%Interpolate contour points to get smooth bounded initial contour:
n = length(P)+1;
P(:,n) = [P(1,1);P(2,1)];
t = 1:n;
ts = 1:0.1:n; %(finer spacing)
P = spline(t,P,ts);
%*************************************************************************
%Make contour clockwise to get outward normals:
% Area inside contour
P = P';
O=[P;P(1:2,:)];
area = 0.5*sum((O((1:size(P,1))+1,1) .* (O((1:size(P,1))+2,2) - O((1:size(P,1)),2))));
% If the area inside  the contour is positive, change from counter-clockwise to 
% clockwise
if(area>0), P=P(end:-1:1,:); end
P = P';
%*************************************************************************
%Internal contour force matrix:
b(1)=beta;
b(2)=-(alpha + 4*beta);
b(3)=(2*alpha + 6 *beta);
b(4)=b(2);
b(5)=b(1);
nPoints = length(P);
A=b(1)*circshift(eye(nPoints),2);
A=A+b(2)*circshift(eye(nPoints),1);
A=A+b(3)*circshift(eye(nPoints),0);
A=A+b(4)*circshift(eye(nPoints),-1);
A=A+b(5)*circshift(eye(nPoints),-2);
    %Find inv(A+tI):
B=inv(A + gamma*eye(nPoints));
%*************************************************************************
%External contour forces
[Gx,Gy] = gradient(im);
Eedge = -1*sqrt(Gx.^2+Gy.^2);
%*************************************************************************
%Fext1(image force attributed to gradient)
sigma2 = 20;
[xg,yg]=ndgrid(floor(-3*sigma2):ceil(3*sigma2),floor(-3*sigma2):ceil(3*sigma2));
DGaussx=-(xg./(2*pi*sigma2^4)).*exp(-(xg.^2+yg.^2)/(2*sigma2^2));
DGaussy=-(yg./(2*pi*sigma2^4)).*exp(-(xg.^2+yg.^2)/(2*sigma2^2));
Fx=imfilter(Eedge,DGaussx,'conv','symmetric');
Fy=imfilter(Eedge,DGaussy,'conv','symmetric');
Fext(:,:,1)=-Fx*2*sigma2^2;
Fext(:,:,2)=-Fy*2*sigma2^2;
%*************************************************************************
%Plot contour:
figure
title('Initial contour');
imshow(I,[]),hold on;
plot(P(2,:),P(1,:),'b.');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Move snake(iterative):
figure
title('moving snake contour');
imshow(I,[]),hold on;

for i = 1:N
    % Get image force on the contour points
    Fext1(:,1)=kappa*interp2(Fext(:,:,1),P(2,:),P(1,:));
    Fext1(:,2)=kappa*interp2(Fext(:,:,2),P(2,:),P(1,:));
    % Interp2, can give nan's if contour close to border
    Fext1(isnan(Fext1))=0;
%****************************************************
    %Fext2(balloon force)
    a=4;%(choosing points relatively close together for stable expansion)
    % From array to separate x,y
    xt=P(1,:); yt=P(2,:);
    % Derivatives of contour
    n=length(xt);
    f=(1:n)+a; f(f>n)=f(f>n)-n;
    l=(1:n)-a; l(l<1)=l(l<1)+n;
    dx=xt(f)-xt(l);
    dy=yt(f)-yt(l);
    % Normals of contourpoints
    L=sqrt(dx.^2+dy.^2);
    nx = -dy./L; 
    ny =  dx./L;
    Normal(:,1)=nx; 
    Normal(:,2)=ny;
    Fext2 = delta*Normal;%(balloon force)
%*****************************************************
    %New positions:
    vnewx = P(1,:)' + gamma*(Fext1(:,1) + Fext2(:,1));
    vnewy = P(2,:)' + gamma*(Fext1(:,2) + Fext2(:,2));
    temp1 = P(1,:);
    temp2 = P(2,:);
    P(1,:) = (B*vnewx)';
    P(2,:) = (B*vnewy)';
    %Stop expanding if contour point is at edge:
    Test = 0;
    for k = 1:length(temp1)
        if I(round(temp1(k)),round(temp2(k)))<100
            P(1,k) = temp1(k);
            P(2,k) = temp2(k);
            Test = 1;
        end
    end
%******************************************************
    plot(P(2,:),P(1,:),'b.');
    drawnow
end
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display image with segmented LA:
figure
title('CT image with segmented Left Atrium');
imshow(I,[]),hold on;
plot(P(2,:),P(1,:),'b.');
drawnow
%*************************************************************************
%Display segmented region as a binary image:
Isize = size(im);
P = P';
J=false(Isize+2);
% Loop through all line coordinates
x=round([P(:,1);P(1,1)]); x=min(max(x,1),Isize(1));
y=round([P(:,2);P(1,2)]); y=min(max(y,1),Isize(2));
for i=1:(length(x)-1)
   % Calculate the pixels needed to construct a line of 1 pixel thickness
   % between two coordinates.
   xp=[x(i) x(i+1)];  yp=[y(i) y(i+1)]; 
   dx=abs(xp(2)-xp(1)); dy=abs(yp(2)-yp(1));
   if(dx==dy)
     if(xp(2)>xp(1)), xline=xp(1):xp(2); else xline=xp(1):-1:xp(2); end
     if(yp(2)>yp(1)), yline=yp(1):yp(2); else yline=yp(1):-1:yp(2); end
   elseif(dx>dy)
     if(xp(2)>xp(1)), xline=xp(1):xp(2); else xline=xp(1):-1:xp(2); end
     yline=linspace(yp(1),yp(2),length(xline));
   else
     if(yp(2)>yp(1)), yline=yp(1):yp(2); else yline=yp(1):-1:yp(2); end
     xline=linspace(xp(1),xp(2),length(yline));   
   end
   % Insert all pixels in the fill image
   J(round(xline+1)+(round(yline+1)-1)*size(J,1))=1;
end
J=bwfill(J,1,1); J=~J(2:end-1,2:end-1);
figure
imshow(J,[]);
%Measure DSC value:
% J2 = logical(imread('CT4_gt.png')); %run after adding file to current folder
% DSC = dice(J,J2)
%END