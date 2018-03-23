% ------------------------------------------------------------------------------------------
% Script for plotting WRF-output variables
%
% Last edited: 22.March.2018, Torgeir
% ------------------------------------------------------------------------------------------

addpath m_map/

% Identify WRF-outfile
wrfpath = '/global/work/blasterdalen/WRFsimMasterThesis/WRF_simulations/WRF_Dec/';
yy = 2014;
mm = 12
dd = 21;
HH = 00;
MM = 00;
SS = 00;
domain = 3;
fileID   = strcat(wrfpath, 'wrfout_d0', num2str(domain), '_', ...
                  num2str(yy), '-', num2str(mm), '-', num2str(dd), '_', ...
                  num2str(HH), ':', num2str(MM), ':', num2str(SS));

% Identify geographical file
geopath = '/home/blasterdalen/Ymse/geo_em.d03.nc';

% time of interest
time = 72;          % 12:00 because outfile has 144 frames

% Extracting 3-D variables (time, south-north, west-east) from WRF
lon = ncread(fileID, 'XLONG');
lat = ncread(fileID, 'XLAT');

% 4-D variables (time, eta, lon, lat)
eta = 4 	% index of 4 gives (normally) height 60 meter a.g.l.
U = ncread(filID, 'U');
V = ncread(filID, 'V');

u = squeeze(U(time, eta, 0:300,0:300));
v = squeeze(V(time, eta, 0:300,0:300));

% Extracting top. height from geo-file
lon_geo = ncread(gepath, 'XLONG_M');
lat_geo = ncread(gepath, 'XLAT_M');
HGT = ncread(geopath, 'HGT_M');

% The terrain in the geo-file is static, therefore we only need at time 1
geoLon = lon_geo(1, :, :);
geoLat = lat_geo(1, :, :);
hgt = HGT(1, :, :);



fig1 = figure(1);
hold all
m_proj('lambert','long',[20 21.5],'lat',[68.9 69.3]);
m_contour(geoLon, geoLat, hgt, 0:200:1400, ...
    'color',[.65 .65 .65],'ShowText','on')

WindField = sqrt(u.^2 + v.^2);
m_contourf(lon, lat, WindField)
colormap jet
colorbar

m_grid('xtick', 20:.5:21.5, 'tickdir', 'out', 'ytick', ...
    68:.5:69.5, 'linest', ':')
m_coast('patch',[.99 .99 .99],'edgecolor','none');


[xR, yR] = m_ll2xy(20.6804, 69.1867);
plot(xR, yR, 'rx')

% Save to file
print('-dpng', 'WindField_d03_60m', '-painters')
