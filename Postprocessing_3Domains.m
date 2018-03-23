% ------------------------------------------------------------------------------------------
% Extracting netcdf variables from WRF-runs and saves variables of interest
% to a .mat-file.
% Explanation of the individual WRF variables can be found by entering
% "ncdisp(path)" in MATLAB. Or by "ncdump -h <filename>" in a terminal.
%
% NetCDF variables (subscripts denote domain numbers):
%   (u,v,w) : 3-D wind velocity
%   PHB : base-state geopotential
%   PH  : perturbation geopotential
%   T2  : Temperature at 2 meters a.g.l.
%   lon : degrees east
%   lat : degrees north
%   PSFC: Surface pressure [Pa]
%   HGT : Terrain height [m]
%   znw : height levels (eta-levels)
%   znu : half height levels
%
% Last edited: 21.March.2018, Torgeir
% ------------------------------------------------------------------------------------------

catalog = 'WRF';                                % catalog containing wrfout-file
dom1    = 'wrfout_d01_2014-12-31_00:00:00';     % name of wrfout-file
dom2    = 'wrfout_d02_2014-12-31_00:00:00';
dom3    = 'wrfout_d03_2014-12-31_00:00:00';
path1   = [catalog '/' dom1];                   % path to outfile
path2   = [catalog '/' dom2];
path3   = [catalog '/' dom3];

% ----------- WRF variables ----------------------------------------------------------------
% 4-D variables: (west-east, north-south, height-level, time)
u1  = ncread(path1, 'U');       u2  = ncread(path2, 'U');       u3  = ncread(path3, 'U');
v1  = ncread(path1, 'V');       v2  = ncread(path2, 'V');       v3  = ncread(path3, 'V');
w1  = ncread(path1, 'W');       w2  = ncread(path2, 'W');       w3  = ncread(path3, 'W');
PHB1= ncread(path1, 'PHB');     PHB2= ncread(path2, 'PHB');     PHB3= ncread(path3, 'PHB');
PH1 = ncread(path1, 'PH');      PH2 = ncread(path2, 'PH');      PH3 = ncread(path3, 'PH');

% 3-D variables: (west-east, south-north, time)
T2_1 = ncread(path1, 'T2');     T2_2 = ncread(path2, 'T2');     T2_3 = ncread(path3, 'T2');
lon1  = ncread(path1, 'XLONG'); lon2  = ncread(path2, 'XLONG'); lon3  = ncread(path3, 'XLONG');
lat1  = ncread(path1, 'XLAT');  lat2  = ncread(path2, 'XLAT');  lat3  = ncread(path3, 'XLAT');
PSFC1 = ncread(path1, 'PSFC');  PSFC2 = ncread(path2, 'PSFC');  PSFC3 = ncread(path3, 'PSFC');
HGT1  = ncread(path1, 'HGT');   HGT2  = ncread(path2, 'HGT');   HGT3  = ncread(path3, 'HGT');

% 2-D variables: (height-level, time)
znw1  = ncread(path1, 'ZNW');   znw2  = ncread(path2, 'ZNW');   znw3  = ncread(path3, 'ZNW');
znu1  = ncread(path1, 'ZNU');   znu2  = ncread(path2, 'ZNU');   znu3  = ncread(path3, 'ZNU');
% ------------------------------------------------------------------------------------------

% Simulation time;
time1          = ncread(path1, 'Times')';
Timestring1    = datenum(time1)';
time2          = ncread(path2, 'Times')';
Timestring2    = datenum(time2)';
time3          = ncread(path3, 'Times')';
Timestring3    = datenum(time3)';

% Actual longitude and latitude of site
R_lon = 20.6804;
R_lat = 69.1867;

% Obtaining the latitude and longitude indexes of the Rieppi site
lonxy1 = lon1(:,:,1);   lonxy2 = lon2(:,:,1);   lonxy3 = lon3(:,:,1);
latxy1 = lat1(:,:,1);   latxy2 = lat2(:,:,1);   latxy3 = lat3(:,:,1);

[nx1,ny1] = size(squeeze(lonxy1));
[nx2,ny2] = size(squeeze(lonxy2));
[nx3,ny3] = size(squeeze(lonxy3));

dist1 = zeros(nx1,ny1); dist2 = zeros(nx2,ny2); dis3t = zeros(nx3,ny3);
r = 6.370e6;        % Radius of the Earth
for i = 1:nx1
   for j = 1:ny1
        dist1(i,j) = r*sqrt(cos((R_lat+latxy1(i,j))/2)^2*(R_lon-lonxy1(i,j))^2 ...
            + (R_lat-latxy1(i,j))^2);
   end
end

for i = 1:nx2
   for j = 1:ny2
        dist2(i,j) = r*sqrt(cos((R_lat+latxy2(i,j))/2)^2*(R_lon-lonxy2(i,j))^2 ...
            + (R_lat-latxy2(i,j))^2);
   end
end

for i = 1:nx3
   for j = 1:ny3
        dist3(i,j) = r*sqrt(cos((R_lat+latxy3(i,j))/2)^2*(R_lon-lonxy3(i,j))^2 ...
            + (R_lat-latxy3(i,j))^2);
   end
end

% Finding the indexes of the closest grid points
[~,index_min1] = min(dist1(:));
[lat_index1,lon_index1] = ind2sub(size(dist1),index_min1);

[~,index_min2] = min(dist2(:));
[lat_index2,lon_index2] = ind2sub(size(dist2),index_min2);

[~,index_min3] = min(dist3(:));
[lat_index3,lon_index3] = ind2sub(size(dist3),index_min3);

% Saving to .mat file for all three domains
outDataPath = '/global/work/blasterdalen/WRF_narvik/WRFnygfj/DataExtracts/'

uu1  = squeeze(u1(lon_index1,lat_index1,:,:));
vv1  = squeeze(v1(lon_index1,lat_index1,:,:));
ww1  = squeeze(w1(lon_index1,lat_index1,:,:));
ph1  = squeeze(PH1(lon_index1,lat_index1,:,:));
phb1 = squeeze(PHB1(lon_index1,lat_index1,:,:));
psfc1= squeeze(PSFC1(lon_index1,lat_index1,:));
t2_1 = squeeze(T2_1(lon_index1,lat_index1,:))';
hgt1 = squeeze(HGT1(lon_index1,lat_index1,:));
save('WRF_outputs/Simulation_December/D01_31_dec.mat','Timestring1','uu1','vv1','ww1','ph1','phb1','psfc1','t2_1','znw1','znu1','hgt1','lonxy1','latxy1')

uu2  = squeeze(u2(lon_index2,lat_index2,:,:));
vv2  = squeeze(v2(lon_index2,lat_index2,:,:));
ww2  = squeeze(w2(lon_index2,lat_index2,:,:));
ph2  = squeeze(PH2(lon_index2,lat_index2,:,:));
phb2 = squeeze(PHB2(lon_index2,lat_index2,:,:));
psfc2= squeeze(PSFC2(lon_index2,lat_index2,:));
t2_2 = squeeze(T2_2(lon_index2,lat_index2,:))';
hgt2 = squeeze(HGT2(lon_index2,lat_index2,:));
save('WRF_outputs/Simulation_December/D02_31_dec.mat','Timestring2','uu2','vv2','ww2','ph2','phb2','psfc2','t2_2','znw2','znu2','hgt2','lonxy2','latxy2')

uu3  = squeeze(u3(lon_index3,lat_index3,:,:));
vv3  = squeeze(v3(lon_index3,lat_index3,:,:));
ww3  = squeeze(w3(lon_index3,lat_index3,:,:));
ph3  = squeeze(PH3(lon_index3,lat_index3,:,:));
phb3 = squeeze(PHB3(lon_index3,lat_index3,:,:));
psfc3= squeeze(PSFC3(lon_index3,lat_index3,:));
t2_3 = squeeze(T2_3(lon_index3,lat_index3,:))';
hgt3 = squeeze(HGT3(lon_index3,lat_index3,:));
save('WRF_outputs/Simulation_December/D03_31_dec.mat','Timestring3','uu3','vv3','ww3','ph3','phb3','psfc3','t2_3','znw3','znu3','hgt3','lonxy3','latxy3')

quit
