function [F2_pol_grd,I_as, f_p]=direct_fourier(phan,prj_angles,natural_0__crop_1__pad_2,outp);
%  Reconstruction using the Central Slice Theorem 
%  (Direct Fourier Tomographic Reconstruction)
%  Gianni Schena, Univ. Trieste, I, schenaATunits.it, may. 2018
%  input:
%  phan : a grayscale image e.g. phantom(128)
%  prj_angles : projection angles for the radon transform (degrees)
%  natural_0__crop_1__pad_2 for different ways to cropp or pad the sinogram
%  values 0 or 1 or 2 (see explanation in the code)
%  outp = 1 for generation of figures 0: no output
%  tested only with Matlab 2016a , requires image processing toolbox 
%%
% close all , clear all
%% some check on the input variables 
if nargin<3
% select important option 0 or 1 or 2 !!
display('missing selectable padding option ! ...then 2 ')
natural_0__crop_1__pad_2=3;
% 3 generous padding
% 2 padding sinogram before fft1
% 1 crops to keeps  the lengh of the sinogram to the same size of the image (side)
outp=1;
end

if nargin<2
display('projection angles missing ! ...then 0:179 ');
prj_angles=(0:0.5:179.5); % projection angles in degrees 

%%%%%
if 0 % random distributed angles (degrees) without repetition
ready = false;
nr_rand = 200; % number of random angles to generate
range = 180; % max angle of the range, [0 range]
while ~ready
   data  = rand(1, nr_rand).*range;
   ready = (numel(unique(data)) == nr_rand);
end
prj_angles=sort(data);
end
%%%%%%

outp=1;
natural_0__crop_1__pad_2=1;
end

if nargin <1
% alternative testing phantoms ,edit here ..   
%n=256; phan=phantom(n);
phan=imread('cameraman.tif');
%load('clown','X'); phan=X(:,:,1); 
%phan=imread('peppers.png'); phan=rgb2gray(phan);
%phan=imread('rice.png'); 
% just to make sure all is fine in the imput image
phan=double(phan); phan=phan./max(phan(:));
phan=imresize(phan,[128 128]); n=size(phan,1);
prj_angles=(0:0.5:179.5); % projection angles in degrees 
% select important option 0 or 1 or 2 !!
natural_0__crop_1__pad_2=2; % 2 for padding before fft1
% ( see bellow for values = 0 or 1 or  2) 
outp=1;
display('no input variables... then run demo example with default values! ')
end

if nargin==3 % check for the value in the range 0-2
    % and correct when out of range 
natural_0__crop_1__pad_2= int8(natural_0__crop_1__pad_2); % 
natural_0__crop_1__pad_2=max(natural_0__crop_1__pad_2,0); % 
natural_0__crop_1__pad_2=min(natural_0__crop_1__pad_2,3); % 
natural_0__crop_1__pad_2;% echo ?
outp=1;
end

if nargin==4 && ~(outp==1 || outp==0 )
    outp=1;
    display('warning: output parameter in input has been corrected')
end

if length(prj_angles)==1
    
display('warning: number of projection angles rather than their equispaced values' ),

if(prj_angles<1) ,
display('warning: projection angles wrong !, num angles less than 1') 
prj_angles=10;
display('warning: num of proj angles corrected to 10') 
end 

in_nprj=prj_angles;
prj_angles=linspace(0,180,in_nprj+1);  
prj_angles=prj_angles(1:end-1);
display('... then projection angle sequence has been computed')
end

if (length(prj_angles)>1 && any(diff(prj_angles)<=0)) || isempty(prj_angles)
 display('warning: projection angles in input not increasing or [], thus corrected with Nyquist number')
 n=size(phan,1);
 na_Nyquist=ceil((pi/2)*n); delta_theta=180/na_Nyquist;
 prj_angles=(0:delta_theta:180-delta_theta); 
 %prj_angles=(0:179);
 disp(['na:  ', num2str( length(prj_angles)), ' as dictated by Naquist principle']);
end

deltas_theta=diff(prj_angles(:));
%prj_angles
if all(deltas_theta~=deltas_theta(1))
display(['warning: not uniform delta angle ! :', num2str(deltas_theta(1))])   
display(['proceeding any way'])   
end

if prj_angles(end) ~=(180-deltas_theta(1))
display ('check last projection angle please: prj_angles(end) + delta =? 180'); 
display('e.g. [1:0.2:179.8] is fine')
display(['your last proj angle is : ' num2str(prj_angles(end))])
display(['your delta seems to be : ' num2str(deltas_theta(1))])
return
end

%% Interpolation method:the same for all the experiments

interp_m='cubic'; % for all but scatterinterpolant
interp_ms='linear'; % for scatteredInterpolant 

n=size(phan,1);
n_phan=n;
%% phantom image included within a disk , 0 outside 
mask=ones(n,n);
if (0) % (1) to generate the including disk 
    disk=logical(0*phan); mask=disk;% 
    szd=size(disk);
    disk(round(szd(1)./2),round(szd(2)./2))=1;
    disk=bwdist(logical(disk));
    mask(disk<=szd(1)./2)=1;
    phan=phan.*double(mask);
  % figure, imshow(phan,[]), title('original phantom in a disk');
else
 mask=ones(size(phan));   
end

%% radon transform (no noise) : create projections 
if 1 % use ifanbeam with infinite distance instead of radon
dtheta = prj_angles(2)-prj_angles(1);
[radon_prj] = fanbeam(phan,10000*size(phan,1),'FanSensorGeometry','line',...
'FanRotationIncrement',dtheta);
radon_prj=radon_prj(2:end-1,(1:round(size(radon_prj,2)./2)));

%size(radon_prj)
%round(size(radon_prj,2)./1),
else
radon_prj=radon(phan,prj_angles);
size(radon_prj)
end
proj_sino=radon_prj;

%%%%
%% create off-set artifacts ? (or use it to further centering the sinogram)
%proj_sino=padarray(proj_sino,1,0,'post'); % to create centre off-set artifacts
%proj_sino=padarray(proj_sino,1,0,'pre'); % to create centre off-set artifacts
%%%%%

% fan-beam acquisition geometry 
%%%aquis_geom='fan_beam';
aquis_geom='par'

if strcmp(aquis_geom,'fan_beam')
% fourier slice theorem does not hold with a fan beam geometry !
% however , if your source is VERY far from the virtual detector and you can assume quasi parallel rays ...  
DST=600;   % distance (virtual ) sensor to source focal spot 
DST_min=sqrt(size(phan,1)^2 + size(phan,2)^2)*2; % minimum reccomended by matlab  
DST=max(DST_min,DST); d_sensor1=1.0; aperture=15.; % fan aperture in degrees
%'Fan Rotation Increment' 
 [proj_sino,position,prj_angles]=...
     fanbeam(phan,DST,'FanSensorGeometry','line','FanSensorSpacing',d_sensor1); 
 [prj_angles(1,[1:2 3 end-1]), prj_angles(1,end),length(prj_angles)]
%  
if 1
  for i=1:180 % proj at alfa mirrows that at 180 + alfa 
  proj_sino(:,i)=(proj_sino(:,i)+proj_sino(end:-1:1,i+180))./2;
 % proj_sino(:,i)=min(proj_sino(:,i),proj_sino(end:-1:1,i+180));
  end
end

% only 0 -180 
prj_angles=prj_angles(1:round(length(prj_angles)/2));
proj_sino=proj_sino(:,1:length(prj_angles));
%length(prj_angles)
%more or less as in Genfire (Alan Pryor et al,2017, Equation [1] page 2)
len_pr=size(proj_sino,1);
Dy= abs((1:len_pr)-len_pr./2); % dist det.pixels dal centro detector
i_D_ipoten=1.0./ sqrt( Dy.^2 + DST.^2); % da focus a pixel sensore
alfa=atan(Dy/DST)/pi*180;
ind_in=alfa<=aperture;
sun_deno=sum(i_D_ipoten(ind_in));% scalare
p_wei=i_D_ipoten./sun_deno; % same length as the projection
% comment oout the next 2 if you do not like the weight
%proj_sino=repmat(p_wei(:),[1,size(proj_sino,2)]).*proj_sino;
%proj_sino(alfa>aperture,:)=0.0;  
proj_sino(~ind_in,:)=0.0;  
end % here ends the fan beam section ( again .. the the DFT does not hold with divergent beams)
%

%%%% DIGITAL FILTER the projections ( rather than analog filter)
% digital filter on the radon transform (rather than analog filter window)
% requires curvefitting toolbox
if 0 % (0) for no digital filter  
[ ' digial filtering of the projection ; rather than the usual filter in the frequency domain!']
order=5;
framelen=7; % odd or it is reduced by 1 to make it odd
%proj_sino = sgolayfilt(proj_sino,order,framelen) % requires signal proc toolbox
% see also on matlab FEX ''Fast Savitzky Golay filter as multi-threaded C-Mex''
for ipr=1:size(proj_sino,2)% filter each columns of the sinog separately.
proj_sino(:,ipr) = smooth(proj_sino(:,ipr),framelen,'sgolay',order);
end

end
%%%%
%% AVAILABLE ANALOG low-pass FILTERS for filtering the projections
%    filter='none',% if you use a digital filter ! probably 'none' is your filter
   filter= 'ram-lak'
%    filter= 'shepp-logan'
%    filter= 'cosine'
%    filter='connes'  % connes windows 
%    filter= 'hamming'
%    filter= 'hann'
%    filter= 'blackman',
%    filter= 'welch'
%    filter='parzen' 
%    filter= 'blackmanharris' 
%    filter='expwin'
%    filter='kaiser'
%    filter='Lanczos'
%    filter= 'udf' % user defined filter 



%% crop projection (see also option [n] in fft)
% dxp>0 is to crop sino to make its length equal to the side
% of the image , dx=0 id you do not like this option
% dxp<0 to pad the sinogram colums with 0 

if mod(length(proj_sino(:,1)),2)==1 % odd to even
   proj_sino= [proj_sino;zeros(1,size(proj_sino,2))]; % 
   % adds one rows of 0 to make it even 
   
   if exist('p_wei','var') 
   p_wei=   [p_wei,0.0]; %
   % proj weinghts same saze as projection 
   % often not used 
   end 
   
   
end

dxp=length(proj_sino(:,1))-size(phan,1) ; 

if natural_0__crop_1__pad_2==1 % crop sinogram length to image side 
dx_up=ceil(dxp./2);
dx_dw=dxp-dx_up; %nc=dx-dx_sx-dx_dx;
nx=length(proj_sino(dx_dw+1:end-dx_up,1));
proj_sino=proj_sino(dx_up+1:end-dx_dw,:); % !! crop
L_pad=0;
display('sinogram length reduced to image size: pad option:1 ')
end

if natural_0__crop_1__pad_2==2  % pad with 0-s the sinogram before fft1
L_pad=max(16,round(n./2)); %  padding the sinogram to double length
proj_sino=padarray( proj_sino,[L_pad 0],0,'both');  
disp('sinogram padded before fft1, pad option:2')
L_pad=L_pad+ceil(((n.*sqrt(2)+2)-n)/2)+1, % to account for the pad added by radon
end

if natural_0__crop_1__pad_2==3  % pad with 0-s the sinogram before fft1
    % i.e very generous paddigg 
L_pad=round(max(64,1.0*n)); %  padding the sinogram to more than double length
proj_sino=padarray( proj_sino,[L_pad 0],0,'both');  
disp('sinogram padded before fft1, pad option:3 (generous padding)')
L_pad=L_pad+ceil(((n.*sqrt(2)+2)-n)/2)+1, % to account for the pad added by radon
end

if natural_0__crop_1__pad_2==0 
    % do nothing with the sinog length
    % keep the sinogram with the length that results from radon:
    % i.e. with the 0 padding given by radon function
    disp('sinogram padded by radon tranform ')
    L_pad=ceil(((n.*sqrt(2)+2)-n)/2)+1,  %pad added by radon
end
% 
%% fft1 sinogram 
f_p=ifftshift(proj_sino,1);
n_p=size(proj_sino,1);
ncr=n;
f_p=fft(f_p,[ ],1); % fourier transf for radon ,[n]: crop here instead 
%%f_p=fftshift(f_p,1); % put DC central 
nfp=length(f_p(:,1));


% filtering 
if exist('filter','var') && ~strcmp(filter,'none') 
    [ 'Fourier filtering with : ',filter]
    H=filtering(filter,nfp);
    H=ifftshift(H); % since I ifftshift the original projction prior fft
    
    % all at onece !
    f_p=f_p.*repmat(H,[1,size(f_p,2)]); % do filter !
    
   % one by one may be iwth different filters 
%     for idf=1:1:size(f_p,2) % each angle
%     f_p(:,idf)=f_p(:,idf).*H(:);
%     end
    
else
     filter
     [ 'no analog filter for projections ']
end

f_p=fftshift(f_p,1); % put DC central after filtering 
%% return to original sinogram length after fft 
if natural_0__crop_1__pad_2==2 % padded
if(0) % almost useless operation ! avoided 
dxp=length(f_p(:,1))-size(phan,1) ; 
dx_up=ceil(dxp./2);
dx_dw=dxp-dx_up; %nc=dx-dx_sx-dx_dx;
nx=length(f_p(dx_dw+1:end-dx_up,1));
f_p=f_p(dx_up+1:end-dx_dw,:); % !! crop
end
end
% final fft1 length 
nfp=length(f_p(:,1)); % final length of the fft1 
%% un-necessary check on the ifft2(fft2(Image))
if(0) % fft2 image :central slice
ft2_phan=ifftshift(phan);
ft2_phan=fft2(ft2_phan,n,n); %fft2 phanntom
ft2_phan=fftshift(ft2_phan);   
phan_reco=fftshift(ifft2(ifftshift(ft2_phan)));
figure, imshow(phan_reco,[]); title('phantom ifft2(fft2)')
end
%% interpolation method : possible selections 
%intertp_m='cubic'; % overwriting the earlier definition

%% gridding for polar to cartesian 
omega_sino=(-(nfp-1)/2:(nfp-1)/2).*(2*pi/size(f_p,1));
theta=prj_angles*pi/180;
omega_image=omega_sino; % let assume same resolution sino % image
% grid that copies(replicates) fft1( sino(:,theta))  coordinates
[theta_grid, omega_grid] = meshgrid(theta,omega_sino); 
% grid image target 
[omega_grid_x, omega_grid_y] = meshgrid(omega_image, omega_image);
[coo_th_fft2, coo_r_fft2] = cart2pol(omega_grid_x,omega_grid_y);
%%%
 coo_r_fft2=coo_r_fft2.*sign(coo_th_fft2); % if theta >pi
 coo_th_fft2(coo_th_fft2<0)=coo_th_fft2(coo_th_fft2<0)+pi;
% %%%%% 

%% put the fft1(sino) radially in the prepared fft2 grid 
Fourier2_radial = interp2(theta_grid,omega_grid,f_p,coo_th_fft2,coo_r_fft2,interp_m,(0+1i.*0));
%'makima'
%Fourier2_radial = interp2(theta_grid,omega_grid,f_p,coo_th_fft2,coo_r_fft2,'makima',(0+1i.*0));


[I_r,I_a,I_i,LIa_1,err11] = Ifft2_2_Img(Fourier2_radial,mask,L_pad,[],phan); 

%% alternative 1 , reconstruction
OMEGA=repmat(omega_sino',[1, length(theta)]);
thetan=theta;
%thetan=theta+0.01.*(rand(1,length(theta)) -0.5);
x=bsxfun(@times, OMEGA, cos(thetan));
y=bsxfun(@times, OMEGA, sin(thetan));

%interp_m='cubic';
Extrapol_Method='nearest';
Fscat=scatteredInterpolant(x(:),y(:),f_p(:),interp_ms,Extrapol_Method); 
%Fcart=TriScatteredInterp(x(:),y(:),f_p(:),interp); 
F2_rad = Fscat(omega_grid_x,omega_grid_y);
[I_r1,I_a1,I_i1,LI_a1,err_sca_13] = Ifft2_2_Img(F2_rad,mask,L_pad,[],phan);

%% alternative 2 ,

Extrapol_Method='nearest';
%Extrapol_Method='';
%
omega_sino=(-(nfp-1)/2:1:(nfp-1)/2).*(2*pi/size(f_p,1));
[theta_grid, omega_grid] = ndgrid(theta,omega_sino); 
f_p1=f_p'; %f_p1=f_p1(:,:);
Fscat = griddedInterpolant(theta_grid,omega_grid,f_p1,interp_m,Extrapol_Method); 
%Fscat = griddedInterpolant(theta_grid,omega_grid,f_p1,'makima',Extrapol_Method); 
[omega_grid_x, omega_grid_y] = ndgrid(omega_sino, omega_sino);
[coo_th_fft2, coo_ro_fft2] = cart2pol(omega_grid_y,omega_grid_x);
%%%%
coo_ro_fft2=coo_ro_fft2.*sign(coo_th_fft2); % cambia segno se theta >pi
% in_0=find(coo_th_fft2(:)==0)
coo_th_fft2(coo_th_fft2<0)=coo_th_fft2(coo_th_fft2<0)+pi;
%%%%%
F2_rad_2=Fscat(coo_th_fft2, coo_ro_fft2);
%F2_rad=flipud(F2_rad);
%F2_rad_2=F2_rad_2.*mask;
[I_r_a2,I_a_a2,I_i_a2,LI_a_a2,err_gi_23] = Ifft2_2_Img(F2_rad_2,mask,L_pad,[],phan); 
%%%

%% alternatice 3, reconstruction based on Polar Gridding codes by Miki Elad
% n_circles x 2 = data size(f_pt,2)
n_cir=ceil(size(f_p,1)./2)+1; % concentric circles, include the 0 radius circle
n_ry=length(theta).*2; % rays : twice the number of diameters
radius=nfp./2-0.5; % max radius
Param=[n_cir,n_ry,radius]; % number of circles , rays and the max radius
%figure, [XP, YP]=Create_Grid_Polar('P',Param,'.b');
[XP, YP]=Create_Grid_Polar('P',Param,[]);
%
XP=XP(:,2:end); YP=YP(:,2:end); % remove 0,0 circle
% allocable fft data twice the n_cir 
XPT=[ XP(n_ry/2+1:n_ry, end:-1:1),  XP(1:n_ry/2,1:end)];
YPT=[ YP(n_ry/2+1:n_ry, end:-1:1),  YP(1:n_ry/2,1:end)];

interp_ms='linear';   
Extrapol_Method='nearest';
f_p_tr=f_p'; % transposed 
F2_rad = scatteredInterpolant(XPT(:),YPT(:),f_p_tr(:),interp_ms,Extrapol_Method); 

[F2_x , F2_y]=meshgrid(-(nfp-1)/2:1:(nfp-1)/2);
F2_pol_grd=F2_rad(F2_x,F2_y);
[I_r,I_a,I_i,LIa,err_pg] = Ifft2_2_Img(F2_pol_grd,mask,L_pad,[],phan); 

%% alternative 4, Exploting bwdist to allocate radially the fft1 of the sinogram in the fft2
ovs=1; % oversampling coef. ovs>=1 e.g 1. , 2 ;
nfp2=round(nfp*ovs);
if mod(nfp2,2)==1, nfp2=nfp2+1; end;  n2=round(nfp2/2); %
Radl=zeros(nfp2,nfp2); Radl(n2,n2)=1; Radl=double(bwdist(Radl,'euclidean')); Radl=Radl./ovs;
%id_core=find(Radl<nfp2./4); % core of the fft2 coef matrix
Radl(n2,n2)=0.5; Radl(n2+1:end,:)=-Radl(n2+1:end,:);
tmp_vect=-(nfp2-1)/2:+1:(nfp2-1)/2;
[Mx,My]=meshgrid(tmp_vect,tmp_vect);
Tht=atan2(-My,Mx); Tht(Tht<0)=Tht(Tht<0)+pi;
[th,ro]=meshgrid(theta,-(nfp-1)/2:1:(nfp-1)/2);
%interp_ms='nearest'; %  neighbor interpolation

if (1)
ro(abs(ro)> n2)=NaN; th(isnan(ro))=NaN; f_p1=f_p; f_p1(isnan(ro))=NaN;
ro = ro(~isnan(ro));th = th(~isnan(th)); f_p1 = f_p1(~isnan(f_p1));
F2=scatteredInterpolant(ro(:),th(:),f_p1(:),interp_ms,Extrapol_Method); %
%F2=scatteredInterpolant(ro(:),th(:),f_p1(:),interp_ms); %
%Radl(Radl>n2)=NaN;  Tht(isnan(Radl))=NaN;
%Radl=Radl(~isnan(Radl));  Tht=Tht(~isnan(Tht));
F2_bwd=F2(Radl,Tht);
%whos F2_bwd
else % TriScatteredInterp will be removed from matlab
F2=TriScatteredInterp(ro(:),th(:),real(f_p(:)),'nearest'); %
r_F2_bwd=F2(double(Radl),double(Tht));
F2=TriScatteredInterp(ro(:),th(:),imag(f_p(:)),'nearest'); %
i_F2_bwd=F2(double(Radl),double(Tht));
F2_bwd=r_F2_bwd+1i.*i_F2_bwd  ; F2_bwd(isnan(F2_bwd)) = 0 ;
end
%F2_bwd=imresize(F2_bwd,1./ovs,'bilinear', 'Antialiasing',true);
if (ovs ~=1) , F2_bwd=imresize(F2_bwd,1./ovs,'bicubic', 'Antialiasing',true); end

[I_r,I_a,I_i,LIabw,err_bwd] = Ifft2_2_Img(F2_bwd,mask,L_pad,[],phan); % ifft2

%% alternative ''very simple code''
s=-(nfp-1)/2:+1:(nfp-1)/2; xi=zeros(length(s),length(theta)); yi=xi;

for it=1:length(theta)
xi(:,it)=s(:)*sin(theta(it));
yi(:,it)=s(:)*cos(theta(it));
end
N=nfp/2;

if (1)
% shows raw reconstruction without any interpolation
occurences=population_ghost(xi,yi,f_p,theta, size(xi,1),L_pad,n_phan,phan);
end

%%return

%Matlab interpolation with griddata
[Mx,My]=meshgrid(s,s);

% index pixels low frequency undergoing NN  interpolation process 
% index pixels low frequency undergoing NN  interpolation process
sl2h=0.3; % [0 , 0.5]
id_core0= find(sqrt(xi(:).^2 + yi(:).^2)<=sl2h*(nfp+2)); % L
ind_ex  = find(sqrt(xi(:).^2 + yi(:).^2)>=sl2h*(nfp-2)); % H
lc0=length(id_core0);
lx=length(x(:));

if lx>lc0
% possibility to use different interpolation methods for low & high feq.
% high freq
'different interpolation method for high & low freq., subplot 21 '

F2_simple = ...
    griddata(xi(ind_ex),yi(ind_ex),f_p(ind_ex),My,Mx,'linear'); % H
%whos F2_simple
%%F2_simple(isnan(F2_simple)) = 0; % extrapolation

% low freq ... core 
if lc0>0
    F2_simple0 = ...
        griddata(xi(id_core0),yi(id_core0),f_p(id_core0),My,Mx,'cubic');%L
   % whos F2_simple0
    F2_simple(~isnan(F2_simple0))=F2_simple0(~isnan(F2_simple0));
    F2_simple(isnan(F2_simple)) = 0; % extrapolation
end


else
 F2_simple = griddata(xi(:),yi(:),f_p(:),My,Mx,'cubic');  
 F2_simple(isnan(F2_simple)) = 0; % extrapolation
end 

%F2_simple(isnan(F2_simple)) = 0; % extrapolation
[I_r,I_as,I_i,LIab_simple,Err_31] = Ifft2_2_Img(F2_simple,mask,L_pad,[],phan); % ifft2

%  tinterp - an alternative to griddata
% version 1.0 (3.58 KB) by Darren Engwirda
% Scattered data
p = [xi(:),yi(:)];     % xy vector
clear xi; clear xi;
t = delaunay(p);   %  delaunayn ! % Triangulation OLD !
%t= delaunayTriangulation(p);
%'delaunaytriangulation'


% Interpolate onto Cartesian grid
F2_Engwirda  = tinterp(p,t,f_p(:),My,Mx,'quadratic');         % Quadratic interpolation via tinterp
F2_Engwirda(isnan(F2_Engwirda)) = 0; % extrapolation 
[I_r,I_a,I_i,LIab_Engwirda,err_DE] = Ifft2_2_Img(F2_Engwirda,mask,L_pad,[],phan); % ifft2    

LIab_rad_base=F2_simple*0.;
rb_methd='regularized';
% if(0) % 
% tic 
% %coord=[xi(:), yi(:)]'; %size(coord)
% values=f_p(:)'; %size(values)
% %zGrid = regularizeNd(x, y, xGrid, smoothness, interpMethod, solver)
% smoothness=[ 0.40E-10 , 0.40E-10 ]; interpMethod='cubic'; solver='\';
% xGrid = cell(1,2); xGrid={s',s'};
% %whos xGrid
% zGrid = regularizeNd([(xi(:)), (yi(:))],values(:),xGrid,smoothness,interpMethod,solver);
% % create girrdedInterpolant function
% F6 = griddedInterpolant(xGrid, zGrid, 'cubic','none');
% %F6 = griddedInterpolant(xGrid, zGrid, 'makima','none');
% F6_fft2=F6(xGrid);
% F6_fft2(isnan(F6_fft2)) = 0.; % extrapolation
% %whos zGrid
% toc 
% [I_r,I_a,I_i,LIab_rad_base] = Ifft2_2_Img(F6_fft2,mask,L_pad,occurences); % ifft2
% end % if radial basis 
%% end simple 

end %end of the main function 
%%
function [GridX,GridY]=Create_Grid_Polar(Kind,Param,Show)
% Written by Miki Elad on March 20th 2005, based on programs written by Dave Donoho.
%=====================================================================
% This function creates a set of points in the plane froming either Cartesian, Polar, or 
% Recto-Polar grids (several variations). A set of parameters dictate the span of these 
% points, their number, and more. 
%
% Synopsis: [GridX,GridY]=Create_Grid(Kind,Param,Show)
%
% Input: 
%    Kind - 'C' for Cartesian, 'P' for Exact Polar, 'R' for exact Recto-Polar grid,  
%               'D' for Donoho and Averbuch Recto-Polar grid, 'S'  for variation on 'D'
%               Recto-Polar that has uniform angle sampling, and 'X' for the Polar grid that 
%               is obtained from the 'S' Recto-Polar by a uniform steps along all rays.
%    Param - A set of parameters dictating the grid size and coverage

%                  'P' - Param=[nr,nt,range] - number of circles/rays and the maximal 
%                          radius to use

%    Show - if non-empty, show a graph of the grid points. Show is used as the 
%                 plot string
%
% Output:
%    GridX,GridY - two 2D arrays containing the grid points
%
% Remark: Exact grids are grid that are pure rotational. Since Donoho-Averbuch grid 
%                has an even number of points along each ray, it is not exact. Thus, 'C', 'P', 
%               and 'R' are exact, while 'D', 'S', and 'X' are not so.
% The sections on 'S' and 'X" are taken from 'Grid_Evolution.m' file.
%
% Examples:
%       figure(1); clf; 

%  

%[XP,YP]=Create_Grid('P',[8,8,1],'.b');

% Written by Miki Elad on March 20th 2005, based on programs written by Dave Donoho.
%=====================================================================

switch Kind

case 'P' % Polar
    
    if nargin==1
        Param=[8,8,1];
    end;
    Nr=Param(1);
    Nt=Param(2);
    MaxRadius=Param(3);    
    Rpoints=0:MaxRadius/Nr:MaxRadius;
    Rpoints=Rpoints(1:end-1);
    Tpoints=0:2*pi/Nt:2*pi;
    Tpoints=Tpoints(1:end-1);
    [GridR,GridT]=meshgrid(Rpoints,Tpoints);
    GridX=GridR.*cos(GridT);
    GridY=GridR.*sin(GridT);
    
    if Show % display the grid
        plot(GridX,GridY,Show);
        axis equal;
        axis([-MaxRadius MaxRadius -MaxRadius MaxRadius]);
    end;
    
end
end % end function 
%%
%function [I_r,I_a,I_i,Lg_I_a] = Ifft2_2_Img(I_fft2,dammy,crop,bw_occurences),
function [Lg_I_a,I_r,I_i,I_a,Error] = Ifft2_2_Img(I_fft2,dammy,crop,bw_occurences,phant)
%I_r=real(target); % real , I_a= abs(target); % module 
%I_i=imag(target); imag , Lg_I_a=log(1+I_a); % for better image visibility 
% input : fft2 , output : images ; 
% crop > 1 crop reconstructed image to original phantom size 

%I_fft2(isnan(I_fft2))=0; % just in case: NaNs left by the Extrapol.

 if nargin >3
 %I_fft2=I_fft2.*double(bw_occurences);% in general better results without it
 end


if ( 0 ) % (1) to experiment with half Fourier space only 
    
fs1 = size(I_fft2,1); fs2 = size(I_fft2,2); 
ms1 = round(size(I_fft2,1)./2)+0; 
ms2 = round(size(I_fft2,2)./2)+0;
if(1) % use upper half of the fft2 
I_fft2=[I_fft2(1:ms1,:); 0*I_fft2(ms1+1:fs1,:)];
else % use left half of the fft2
I_fft2=[I_fft2(:,1:ms2), 0.*I_fft2(:,ms2+1:fs2)];
end

end % experiments with half F2


if (1) %  (1) cuts high frequencies & mitigates noise in the reconstruction
    sz=size(I_fft2,1); mask=ones(sz,sz);
    % one of the following filterng strategies 
    %%
    if 0 % rough low pass mask rather than a circles or ellipses ()
        nrzo=round(0.3*sz); % parameter
        mask(( [1:nrzo, (sz-nrzo+1:sz)]),( [1:nrzo, (sz-nrzo+1:sz)]))=0;
     %   figure, imshow(mask,[]),
     I_fft2=I_fft2.*complex(real(mask),ones(size(mask)));
    end
    %%
    if 0, % ideal low pass filter 
        mask=bwcircle(sz,sz.*0.95,sz.*0.95); % disk 
        % (size, diameter vertical ,diam.  horizontal )
     %   figure, imshow(mask,[])
        % used to zeroing the high frequencies in the fft2
     I_fft2=I_fft2.*complex(real(mask),ones(size(mask)));   
    end
    %%
    if 1,
   %%% I_fft2=butterworth_2D(I_fft2,0.1*sz,0.9.*sz,2); % this code has not been scrutinazed 
    Dl=0.6*sz;
    kl = butter_lp_kernel(sz, Dl, 2);
    I_fft2=I_fft2.*kl;
    end
    %%
    
end

% if exist('mask','var')
% %     [' mask exists ']
% % I_fft2=ifftshift(I_fft2);    
% % I_fft2=I_fft2.*complex(real(mask),ones(size(mask)));
% % target=fftshift(ifft2(I_fft2));    
% % else
% %     [' mask does not exists ']
target=fftshift(ifft2(ifftshift(I_fft2))); 
%end

target=(target-min(target(:)))./(max(target(:))-min(target(:))); % nolrmalize [0 1]
%target = imadjust(target,[min(target(:)); max(target(:))],[0; 1])

if crop>0
  
    target=target(crop+1:end-crop,crop+1:end-crop);
    
end

%sz_after=size(target)


I_r=real(target); % real 
I_a= abs(target); % module 
I_i=imag(target);
Lg_I_a=log(1+I_a); % for better image visibility 
%Lg_I_a=I_r;
%Lg_I_a=I_i;

Error=NaN;
if 1
    % computer reconstruction error
    flag=exist('phant','var') && ~isempty(phant);
    if flag
  
 %  Error=  immse(I_r,single(phantom));
 %  fprintf('\n The Mean-Squared Error is %0.4f', Error);
   Error=  ssim(single(I_r),single(phant));
   ['The Similarity Index is : ', num2str(Error)],
   
    end

end



end
%


function r=dist2YA(A)
A=A';
c1 = bsxfun(@minus, A(:,1), A(:,1)');
c2 = bsxfun(@minus, A(:,2), A(:,2)');
r = sqrt(c1.^2 + c2.^2); %sym by def
end

function fi = tinterp(p,t,f,xi,yi,method)

% Interpolation of scattered data.
%
%   fi = tinterp(p,t,f,xi,yi,method);
%
% INPUTS:
%
%   p = [x1,y1; x2,y2; etc]             - an Nx2 matrix defining the [x,y] 
%                                         coordinates of the scattered vertices 
%
%   t = [n11,n12,n13; n21,n22,n23; etc] - an Mx3 matrix defining the
%                                         triangulation of the points in p, as 
%                                         returned by delaunayn(p)
%
%   f = [f1; f2; etc]                   - an Nx1 vector defining the value
%                                         of the function at the scattered 
%                                         vertices
%
%   xi,yi                               - Vectors or matrices of identical
%                                         size, defining the [x,y] coordinates
%                                         where the function is to be
%                                         interpolated
%
%   method = 'linear','quadratic'       - Type of interpolant, method is
%                                         optional and if not passed the 
%                                         'linear' method is used
%
% OUTPUT:
%
%   fi                                  - The interpolation evaluated at
%                                         all points in xi,yi. 
%
% The 'linear' method should return identical results to griddata. The
% 'quadratic' method is NEW! (I think) and should return more accurate results 
% when compared with the 'cubic' method in gridata (See tester.m).  
%
% Example:
%
%   p = rand(100,2);
%   t = delaunayn(p);
%   f = exp(-20*sum((p-0.5).^2,2));
%
%   figure, trimesh(t,p(:,1),p(:,2),f), title('Exact solution on scattered data')
%
%   x     = linspace(0,1,50);
%   [x,y] = meshgrid(x);
%   fL    = tinterp(p,t,f,x,y);
%   fQ    = tinterp(p,t,f,x,y,'quadratic');
%
%   figure, mesh(x,y,fL), title('Linear interpolation')
%   figure, mesh(x,y,fQ), title('Quadratic interpolation')
%
% See also TESTER, DELAUNAYN, GRIDDATA

% Any comments? Let me know:
%
%   d_engwirda@hotmail.com
%
% Darren Engwirda - 2006


% Check I/O
if nargin<6
    method = 'linear';
elseif nargin>6
    error('Too many inputs');
end
if nargout~=1
    error('Wrong number of outputs');
end

% Force lowercase
method = lower(method);

% Error checking
if size(p,2)~=2
    error('p must be an Nx2 vector as passed to delaunayn');
end
if size(t,2)~=3
    error('t must be an Mx3 matrix containing the triangulation of p as per delaunayn');
end
if (max(t(:))>size(p,1))||(min(t(:))<=0)
    error('t is not a valid triangulation for p');
end
if (size(f,2)~=1)||(size(f,1)~=size(p,1))
    error('f must be an Nx1 vector defining the function f(p) at the vertices')
end
if ~strcmp(method,'linear')&&~strcmp(method,'quadratic')
    error('method must be either ''linear'' or ''quadratic'' ');
end
sizx = size(xi); sizy = size(yi);
if any(sizx-sizy)
    error('xi & yi must be the same size');
end
pi = [xi(:),yi(:)];

% Allocate output
fi = zeros(size(pi,1),1);

% Find enclosing triangle of points in pi

%i = tsearch(p(:,1),p(:,2),t,pi(:,1),pi(:,2)); % old 
% ( vertices, triangles, interrogation/spatial query) 

 i = tsearchn(p,t,pi); % OLD function in connection with delaunay
% [i, ~] = pointLocation(t, pi); NEW in connection with
% trialgolarizatiodelaunay

% tic
% [~,i,~]=findtria(p,t,pi);
% toc
% try Findtria on FEX as a faster alternative to tsearchn

% Deal with points oustide convex hull
in = ~isnan(i); 

% Keep internal points
pin = pi(in,:); 
tin = t(i(in),:);

% Corner nodes
t1 = tin(:,1); t2 = tin(:,2); t3 = tin(:,3);

% Linear shape functions
dp1 = pin-p(t1,:);
dp2 = pin-p(t2,:);
dp3 = pin-p(t3,:);
A3  = abs(dp1(:,1).*dp2(:,2)-dp1(:,2).*dp2(:,1));
A2  = abs(dp1(:,1).*dp3(:,2)-dp1(:,2).*dp3(:,1));
A1  = abs(dp3(:,1).*dp2(:,2)-dp3(:,2).*dp2(:,1));

% Interpolant
if strcmp(method,'linear')      % Linear interpolation
    
    % Take weighted average of constant extrapolations
    fi(in) = (A1.*f(t1)+A2.*f(t2)+A3.*f(t3)) ./ (A1+A2+A3);
    
else                            % Quadratic interpolation
    
    % Form gradients
    [fx,fy] = tgrad(p,t,f);
    
    % Take weighted average of linear extrapolations
    fi(in) = ( A1.*(f(t1) + dp1(:,1).*fx(t1) + dp1(:,2).*fy(t1)) + ...
               A2.*(f(t2) + dp2(:,1).*fx(t2) + dp2(:,2).*fy(t2)) + ...
               A3.*(f(t3) + dp3(:,1).*fx(t3) + dp3(:,2).*fy(t3)) ) ./ (A1+A2+A3);
           
end

% Deal with outside points
if any(~in)
    fi(~in) = NaN;
end

% Reshape output
fi = reshape(fi,sizx);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fnx,fny] = tgrad(p,t,f)

% Approximate [df/dx,df/dy] at the vertices.

% Nodes
t1 = t(:,1); t2 = t(:,2); t3 = t(:,3);

% Evaluate centroidal gradients (piecewise-linear interpolants)
x23 = p(t2,1)-p(t3,1);  y23 = p(t2,2)-p(t3,2);
x21 = p(t2,1)-p(t1,1);  y21 = p(t2,2)-p(t1,2);

% Centroidal values - [df/dx,df/dy]
den =  (x23.*y21-x21.*y23);
fcx =  (y23.*f(t1)+(y21-y23).*f(t2)-y21.*f(t3)) ./ den;
fcy = -(x23.*f(t1)+(x21-x23).*f(t2)-x21.*f(t3)) ./ den;

% Calculate simplex quality
q = quality(p,t);

% Form nodal gradients.
% Take the quality weighted average of the neighbouring centroidal values.
% Low quality simplices ("thin" triangles) can contribute inaccurate
% gradient information - quality weighting seems to limit this effect.
fnx = 0*f; 
fny = fnx; 
den = fnx;
for k = 1:size(t,1)
    % Nodes
    n1 = t1(k); n2 = t2(k); n3 = t3(k);
    % Current values
    qk = q(k); qfx = qk*fcx(k); qfy = qk*fcy(k); 
    % Average to n1
    fnx(n1) = fnx(n1)+qfx;
    fny(n1) = fny(n1)+qfy;
    den(n1) = den(n1)+qk;
    % Average to n2
    fnx(n2) = fnx(n2)+qfx;
    fny(n2) = fny(n2)+qfy;
    den(n2) = den(n2)+qk;
    % Average to n3
    fnx(n3) = fnx(n3)+qfx;
    fny(n3) = fny(n3)+qfy;
    den(n3) = den(n3)+qk;
end
fnx = fnx./den;
fny = fny./den;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = quality(p,t)

% Approximate simplex quality

% Nodes
p1 = p(t(:,1),:); p2 = p(t(:,2),:); p3 = p(t(:,3),:);

% Approximate quality
d12 = p2-p1;
d13 = p3-p1;
q   = 3.4641*abs(d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))...
                            ./(sum(d12.^2,2)+sum(d13.^2,2)+sum((p3-p2).^2,2));
end


function pop = population_ghost(x,y,ft,theta,N,crop,n_phan,phan);
% shows ffft2 positions populated without/before any kind of interpolation
% shows also rough reconstruction based on non interpolated fft2
pop=zeros(N,N); val_ft=pop; I_fft2=pop;
size(ft)

%%
for ic=1:size(x,2);
xr=round(x(:,ic)+N/2); xr=max(1,xr);
yr=round(y(:,ic)+N/2); yr=max(1,yr);

linearInd = sub2ind([N,N],xr,yr);

pop(linearInd)=1;
val_ft(linearInd)=ft(:,ic);
end


%%%

if 0
filtertype='lp'; % lp : low pass
lambda=4;    
%Zf = filt2(val_ft,1,lambda,filtertype) ; % from matlab FEX
%Zfr=circ_filt(real(val_ft));  Zfi=circ_filt(imag(val_ft)); FEX
%%Zfi=imag(val_ft);
Zf=Zfr+1.i*Zfi;
else
Zfr=medfilt2(real(val_ft)); Zfi=medfilt2(imag(val_ft)); Zf=Zfr+1.i*Zfi;
Zf(val_ft~=0)=val_ft(val_ft~=0); % keep original values
end

figure
subplot(2,4,1)
imshow(pop,[]);

title(' fft2 population w/o interpolation')
subplot(2,4,2)
imshow(log(1+abs(real(val_ft))),[]);%log(1+abs(fft2 FBP))
sparcity_fft2=nnz(pop)/prod(size(pop))
title([' real fft2 coef.s sparcity :', num2str( sparcity_fft2)])

subplot(2,4,3)
I_fft2=val_ft;

if ( 0 ) % (1) to experiment with half Fourier space only

fs1 = size(I_fft2,1); fs2 = size(I_fft2,2);
ms1 = round(size(I_fft2,1)./2)+0;
ms2 = round(size(I_fft2,2)./2)+0;

if(0) % use upper half of the fft2
I_fft2=[I_fft2(1:ms1,:); 0*I_fft2(ms1+1:fs1,:)];
else % use left half of the fft2
I_fft2=[I_fft2(:,1:ms2), 0.*I_fft2(:,ms2+1:fs2)];
end

end % experiments with half Fft2

target=fftshift(ifft2(ifftshift(I_fft2)));

size_target=size(target)
if crop>0,  target=target(crop+1:end-crop,crop+1:end-crop); end
size_target=size(target)

target=flipud(target);
imshow(abs(target),[ ]); title('ifft2 before any interpol. (min energy)')
error=  ssim(single(real(target)),single(phan));
xlabel(['Struc Simil index:', num2str(error)]);

%%

%%
if 1 % duplication of rays
f_p1=ft.*0; f_p1(:,end)=[];
nfp=size(ft,1);
s=-(nfp-1)/2:+1:(nfp-1)/2;

for it=1:length(theta)-1
theta1=(theta(it)+theta(it+1)).*0.5;
xi1(:,it)=s(:)*sin(theta1);
yi1(:,it)=s(:)*cos(theta1);
f_p1(:,it)=min(ft(:,it),ft(:,it+1));

xr1=round(xi1(:,it)+N/2); xr1=max(1,xr1);
yr1=round(yi1(:,it)+N/2); yr=max(1,yr1);

linearInd = sub2ind([N,N],xr1,yr1);
% put fft1 into the polar cartesian matrix
% ... keep putting in for it ha sbeen used already
val_ft(linearInd)=f_p1(:,it);
end

pause(1)

I_fft2=val_ft;
target=fftshift(ifft2(ifftshift(I_fft2)));
if crop>0,  target=target(crop+1:end-crop,crop+1:end-crop); end

size_target=size(target)

subplot(2,4,4)
target=flipud(target);
imshow(abs(target),[ ]);
title('raw ifft2 with rays duplicated')
sparcity_fft2=nnz(val_ft)/prod(size(pop))
end

subplot(2,4,7)
imshow(phan,[]); title(['phantom: ', num2str(n_phan)])

subplot(2,4,8)
%Zf!!fft2 filtrata 
target=fftshift(ifft2(ifftshift(Zf))); 
if crop>0,  target=target(crop+1:end-crop,crop+1:end-crop); end
target=flipud(target);
size_target=size(target)
imshow(abs(target),[ ]); title('2d low pass filtering')

subplot(2,4,5)

magni=1; % [suggested : magni=2]
% less than 1 produces smooth 'foggy' results 
% larger than 1 'sandy' results 
% for very large values ... artifacts  
%xedges=(-N/2:magni:+N/2);
%lenght_xedges=length(xedges)
xedges=(-(nfp-1)/2:(nfp-1)/2);
%xedges = powmspace_m(-round(N./2), round(N./2), 1, length(xedges));
lenght_xedges=length(xedges)
yedges=xedges; 
end % end population


function w = parzenwin_1(n)
%PARZENWIN Parzen window.
%   PARZENWIN(N) returns the N-point Parzen (de la Valle-Poussin) window in a column vector.
%
%   See also BARTHANNWIN, BARTLETT, BLACKMANHARRIS, BOHMANWIN,
%            FLATTOPWIN, NUTTALLWIN, RECTWIN, TRIANG, WINDOW.

%   Reference:
%     [1] fredric j. harris [sic], On the Use of Windows for Harmonic
%         Analysis with the Discrete Fourier Transform, Proceedings of
%         the IEEE, Vol. 66, No. 1, January 1978


%   Author(s): V.Pellissier
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2002/11/21 15:46:25 $

% Check for valid window length (i.e., n < 0)


% Index vectors
k = -(n-1)/2:(n-1)/2;
k1 = k(k<-(n-1)/4);
k2 = k(abs(k)<=(n-1)/4);

% Equation 37 of [1]: window defined in three sections
w1 = 2 * (1-abs(k1)/(n/2)).^3;
w2 = 1 - 6*(abs(k2)/(n/2)).^2 + 6*(abs(k2)/(n/2)).^3;
w = [w1 w2 w1(end:-1:1)]';

end

function [H]=filtering(filter,nfp);

 d=0.95;
 filt=[2*(0:nfp/2)./nfp]';
 w = [2*pi*(0:length(filt)-1)/nfp]';   % frequency axis up to Nyquist
% 
%filter
    switch filter
        case 'none'
            filt=[ones(1,nfp/2+1)]'; % thst is : Do nothing
             'projection filter: none '
        case 'ram-lak'
            filt=[2*(0:nfp/2)./nfp]';
            'projection filter: ram-lak '
        case 'shepp-logan'
            % be careful not to divide by 0:
            filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
             'projection filter: shepp-logan '
        case 'cosine'
           
            filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
        case 'hamming'
            filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
            'projection filter: hamming '
        case 'hann'
            filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
             'projection filter: hann '
        case 'blackman'
            filt(2:end) = filt(2:end) .*(0.42 - 0.5*cos(w(2:end)./d) + ...
                0.08*cos(2.*w(2:end)./d));
             'projection filter: blackman '
        case 'welch'
            filt(2:end) = filt(2:end) .*(1-(w(2:end)/d).^2);
            'projection filter : welch'
        case 'connes'    
            filt(2:end) = filt(2:end) .*(1-(w(2:end)/d).^2).^2;
            
            'projection filter : connes '
        case 'blackmanharris' % Blackman-Harris window
             % check !!
            t=w(2:end)./d;
            filt(2:end)= filt(2:end).*...
                (  0.35875-.48829*cos(t)+.14128*cos(2*t)-.01168*cos(3*t));  % half window
        case 'udf' % user defined 
            %filt=[ones(1,nfp/2+1)]'; % that is : Do nothing
              filt=(2*(0:nfp/2)./nfp)'; % r-l
              filt( filt < 0.2 ) = 0.2;
             'projection filter : user defined filter'
        case 'expwin'
        beta=1; N=nfp; M = (N-1)/2; 
        for k = 0:M
        n = k-M;
        filt(k+1) = exp(beta*sqrt(1-4*n*n/(N-1)/(N-1)))/exp(beta);
        end
            
        case 'parzen'
            filt=parzenwin_1(nfp);
            'projection filter: parzen '
        case 'kaiser'
            %  Joe Henning - Dec 2013
            beta=1 ;% 1 kaiser,0 Rectangular,5 Hann,6 Hamming,12.2 Blackman-Harris
            p=1;  N=nfp; M = (N-1)/2;
            w = [];
            for k = 0:M
                n = k-M;
                w(k+1) = besseli(0,beta*sqrt(1-4*n*n/(N-1)/(N-1)))/besseli(0,beta);
                % w(N-k) = w(k+1);
            end
            w = w.^p;  w = w(:)/max(w);
            'projection filter: kaiser '
        case 'Lanczos'
            n=length(filter) ;
            if ~rem(n,2)
                % Even length window
                m = n/2;
                x = (0:m-1)'/(n-1);
                w = sinc(2*x - 1);
                %  w = [w; w(end:-1:1)];
            else
                % Odd length window
                m = (n+1)/2;
                x = (0:m-1)'/(n-1);
                w = sinc(2*x - 1);
                %  w = [w; w(end-1:-1:1)];
            end
            filt=w;
            
            
        otherwise
            error('Invalid selection: the filter selected in not listed.');
    end
    
    %%filt
    
    H=[filt;filt(end-1:-1:1)];
    H(nfp+1:end)=[];
    % show frquency filter 
    % figure, plot(H,'r+'), title(['frequency filter: ', filter])

end
% SINC   Sinc function.
%    sinc(X) returns sin(pi * x) / (pi * x) for elements of X. X may be 
%    vector or matrix. sinc(0) = 1.

% Author: Gene Dial, GeoEye, 2011-07-28.

function y = sinc(x)
y = ones(size(x));
isNonZero = x~=0;
y(isNonZero) = sin(pi*x(isNonZero))./(pi*x(isNonZero));
end

function mask= bwcircle(n,dh,dw)
if nargin==1 , dh=n; end
% mask is 1 within the circle of diameter d,
% 0 outside

if nargin==2 % circle 
disk=false(n,n); mask=disk;%
r=dh/2;
disk(round(n/2),round(n/2))=1;
disk=bwdist(logical(disk),'euclidean');
mask(disk<=r)=true;
end

if nargin==3
% mak with an ellipse 
% dw &  dh are diameters vertical and horizontal
h=dh/2; % vertical 
w=dw/2; % hor
[X, Y ]=meshgrid(1:n,1:n);
mask=((X-n./2)/w).^2 + ((Y-n./2)/h).^2 <=true;
end

%figure, imshow(mask,[]);
wedge=0; %%%bow-tie type wedge 
% (1)
if wedge
u=n/10; 
yc=[n/2,n/2+u,n/2-u,n/2];
xc=[n/2,n,n,n/2];
BW=~poly2mask(xc,yc,n,n);
BW= bitand( BW , BW(:,end:-1:1));
%figure, imshow(BW,[]);
% intersect circle with Horiz. bow tie
mask=bitand((mask) ,BW);
% also Vertical bow tie 
% BW=rot90(BW); 
% mask=bitand(mask,BW);
%mask=~mask;
end
%
%figure, imshow(mask,[]);title('mask for freq. cut-off')
%pause
end

function kh = butter_hp_kernel(N, Dh, no),
% butterworth high pass filter kernel 
D=false(N,N);
r=N/2;
D(round(r),round(r))=1;
D=bwdist(logical(D),'euclidean');
denom = 1+(sqrt(2)-1).*(D./Dh).^(2*no);
kh = 1-1./(denom);
%figure, imshow(kh,[]); title('high pass butterw')
end


function kl = butter_lp_kernel(N, Dl, no),
% butterworth low pass filter kernel 
% no: order , N : size of the matrix filter
D=false(N,N);
r=N/2;
D(round(r),round(r))=1;
D=bwdist(logical(D),'euclidean');
denom = 1+(sqrt(2)-1).*(Dl./D).^(2*no);
kl = 1-1./(denom);
%figure, imshow(kl,[]); title ('low pass butterw')
['compute butterworth low pass filter kernel'],
end







