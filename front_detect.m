
% Example code for frontal detection following algorithm
% of Simmonds et al. (2012) and Schemm et al. (2015)

% Author: Erik Swenson

% Input ERA-Interim (from file 'erai_front.nc'):
%   longitude ('lon'), latitude ('lat'),
%   10 m winds on 18Z10May2017 ('u1' and 'v1')
%          and on 00Z11May2017 ('u2' and 'v2')

% Output (to same file 'erai_front.nc'):
%   frontal detection on 00Z11May2017 ('front')

%   front==1 at grid points where front is detected
%   front==0 everywhere else


clearvars;

% Reads in necessary fields
fname = 'erai_front.nc';
lon = ncread(fname,'lon'); nlon = length(lon);
lat = ncread(fname,'lat'); nlat = length(lat);
% 10 m winds 6 hours prior (1) and current (2)
U1 = ncread(fname,'u1');
V1 = ncread(fname,'v1');
U2 = ncread(fname,'u2');
V2 = ncread(fname,'v2');

% Grid spacing (km)
dx = cos(lat*pi/180.0)*2*pi*6370/360*(lon(nlon)-lon(1))/(nlon-1);
dy = ((lat(2)-lat(1))/180.0)*6370*pi;

front = zeros(nlon*nlat,1);
shift = zeros(nlon,nlat);

% Must be westerly with meridional wind shift from
% negative to positive with at least 2 m/s difference:
% shift==1 where detection criteria satisfied
shift(find(U1>0.0&V1<0.0&U2>0.0&V2>0.0&(V2-V1)>2.0)) = 1;

% label groups of contigious grid points (objects)
% with connectivity = 8 as potential fronts
obj = watershed(-shift,8);
obj(find(shift==0)) = 0;

% Loop over objects
   for z = 1:max(max(obj))
% Location indices
      iobj = find(obj==z); nobj = size(iobj,1);
% Longitude and latitude indices
      [ilon,ilat] = ind2sub(size(obj),iobj);
% Mean latitude
      mlat = mean(lat(ilat));
      [dum,imlat] = min((lat-mlat).^2);
% Area
      area = dy.*sum(dx(ilat));
% Extent (+/- 2 sigma of ellipse with same area)
      xlen = dx(imlat)*4*std(ilon); ylen = dy*4*std(ilat);
% Length
      length = sqrt(ylen^2+xlen^2);
% Length must be at least 500 km to be front
      if (length>500.0)
         front(find(obj==z)) = 1.0;
      end
   end

% Writes out frontal detection 
front = reshape(front,nlon,nlat);
ncwrite(fname,'front',front);


