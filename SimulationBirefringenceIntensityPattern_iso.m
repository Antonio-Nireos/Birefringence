%% Script simulating the birefringence intensity pattern for a alpha-BBO crystal with the optical axis oriented along the propagation direction
%Author: Antonio Perri CTO at NIREOS www.nireos.com 


clear all

%thickness in [m] of the crystal
t1=0.004;

xxx=linspace(-0.12,0.12,1000); %horizontal range of angles in radians
yyy=linspace(-0.12,0.12,1000); %Vertical range of angles radians


[XXX,YYY] = meshgrid(xxx,yyy); %Creating mesh
[w,i]=cart2pol(XXX,YYY); %converting to polar coordinates i=angle of incidence with respect to the crystal front face normal, w= polar angle with respect to the horizontal direction

%loading birefringence values of bbo from refractiveindex.info
load('birefringence.mat');

%calculating no and ne at the chosen wavelength
lambda=0.63; %Wavelength in um
[intens2d_r]=getBirefringencePattern(no,ne, w,i,lambda,t1);
lambda=0.550; %Wavelength in um
[intens2d_g]=getBirefringencePattern(no,ne, w,i,lambda,t1);
lambda=0.450; %Wavelength in um
[intens2d_b]=getBirefringencePattern(no,ne, w,i,lambda,t1);

%creating rgb image
rgb(:,:,1)=intens2d_r;
rgb(:,:,2)=intens2d_g;
rgb(:,:,3)=intens2d_b;

%Imaging result

image(xxx*180/pi,yyy*180/pi,rgb)

title(['Interference pattern']);
xlabel('Incidence Angle [deg]');
A=ylabel('Incidence Angle [deg]');
axis image
set(gca,'FontName','Courier New','FontWeight','bold');



function [intens2d]=getBirefringencePattern(no,ne, w,i,lambda,t1)
% Obtains the selected intensity pattern from the refractive indexes no and
% ne of the birefringent crystal for the specify angular ranges i and w,
% wavelength lambda and crystal thickness t1


%calculating no and ne at the chosen wavelength
[~,index]=min((ne(:,1)-lambda).^2);
a=1/ne(index,2);
b=1/no(index,2);

%defining arrays
deltap_1=zeros(size(w,1),size(w,2));
intens2d=deltap_1;

for jj=1:size(w,1)
    for kk=1:size(w,2)

%Computing phase difference between ordinary and extraordinary beams
deltap_1(jj,kk)=t1*(sqrt((1+(-a^2+b^2).*sin(i(jj,kk)).^2-b.^2.*sin(i(jj,kk)).^2)./b^2)-1/b*sqrt(1-sin(i(jj,kk)).^2*b^2));
%compute intensity patterns
intens2d(jj,kk)=0.5*((cos(w(jj,kk))).^4+(sin(w(jj,kk))).^4-2*(cos(w(jj,kk))).^2.*(sin(w(jj,kk))).^2+2*((cos(w(jj,kk))).^2.*(sin(w(jj,kk))).^2-0.5.*((cos(w(jj,kk))).^4+(sin(w(jj,kk))).^4)).*(cos((2*pi/lambda)*10^6*(deltap_1(jj,kk)))));


    end
end

end




