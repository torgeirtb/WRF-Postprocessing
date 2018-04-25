%:------------------------------------------------------------------------------------------
% Extracting netcdf variables from WRF-runs and saves variables of interest
% to a .mat-file.
% 
% - NOTE: -
% This script is written to be memory efficient than the script 
% "PostProcessing_1domain.m", and requires that the
% point of interest is the centerpoint of the domain!
% ---------
% 
% Explanation of the individual WRF variables can be found by entering
% "ncdisp(path)" in MATLAB. Or by "ncdump -h <filename>" in a terminal.
%
% Some NetCDF variables:
%   XLAT  : degrees north, South is negative. [degrees]
%           (west-east, south-north, time). 
%   XLONG : degrees east, West is negative.             --
%   U     : x-wind component [m/s]
%           (west-east - staggered, south-north, eta, time). 
%   V     : y-wind component [m/s]
%           (west-east, south-north - staggered, eta, time). 
%   W     : z-wind component      --
%   T2    : Temperature at 2 meters a.g.l. [K]
%   T     : Perturbation potential temperature [K] 
%   PSFC: Surface pressure [Pa]
%   HGT : Terrain height [m]
%    
% Order of arguments (in MATLAB):
% 3-D: (west-east, south-north, time) 
% 4-D: (west-east, south-north, bottom-top, time) 
% 
% 
% Last edited: 24.April.2018, Torgeir
%:------------------------------------------------------------------------------------------


% Identify WRF-outfile and where to store extracted data
wrfpath = '/global/work/blasterdalen/WindCoE/WRF_rawoutput/';
yyyy = 2014;
mm = 03;
dd = 20;
HH = 00;
MM = 00;
SS = 00;
domain = 2;
fileID = strcat(wrfpath, 'wrfout_d0', num2str(domain), '_', ...
                num2str(yyyy), '-', num2str(mm, '%02d'), '-', ...
                num2str(dd, '%02d'), '_', num2str(HH, '%02d'), ':', ...
                num2str(MM, '%02d'), ':', num2str(SS, '%02d'));

outDataPath = '../WRF_dataextracts/';

% Eta-levels of interest. 1-6 bives approx. 0 - 150 m 
eta = 1:6;     


% ---- Destaggering wind components ---------------------------------------
% x-component, West - East destaggering 
U  = ncread(fileID, 'U');

% Obtain two center points (icp = West-East, jcp = South-North)
icp = [round(size(U, 1)*.5) round(size(U, 1)*.5) + 1];
jcp = round(size(U, 2)*.5);
% Iterate over all timesteps at different eta-levels to obtain destaggered
% x-wind component
u = NaN(length(eta), size(U, 4));
for t = 1:length(u)
    for z = 1:length(eta)
        u(z, t) = 0.5*(U(icp(1), jcp, z, t) + U(icp(2), jcp, z, t));
    end
end
clear *cp t U z


% y-component, South - North destaggering 
V  = ncread(fileID, 'V');
% Obtain two center points
icp = round(size(V, 1)*.5);
jcp = [round(size(V, 2)*.5) round(size(V, 2)*.5)+1];
% Destaggering the y-component
v = NaN(length(eta), size(V, 4));
for t = 1:length(v)
    for z = 1:length(eta)
        v(z, t) = 0.5*(V(icp, jcp(1), z, t) + V(icp, jcp(2), z, t));
    end
end
clear *cp t V z

% z-component, bottom - top destaggering 
W  = ncread(fileID, 'W');
% Obtain two center points
icp = round(size(W, 1)*.5);          % west-east centerpoints
jcp = round(size(W, 2)*.5);          % south-north centerpoint
% Destaggering the z-components
w = NaN(length(eta), size(W, 4));
for t = 1:length(w)
    for z = 1:length(eta)
        w(z, t) = 0.5*(W(icp, jcp, z, t) + W(icp, jcp, z + 1, t));
    end
end
clear *cp t W z


% ---- Untaggered points (allready in theta-coordinates) ------------------
% 3-D variables 
XLAT = ncread(fileID, 'XLAT');
XLON = ncread(fileID, 'XLONG');
T2 = ncread(fileID, 'T2');
PSFC = ncread(fileID, 'PSFC');
HGT  = ncread(fileID, 'HGT');

% 4-D variables 
QVAPOR = ncread(fileID, 'QVAPOR');
T = ncread(fileID, 'T');
P = ncread(fileID, 'P');
PB = ncread(fileID, 'PB');


% Obtain two center points
icp = round(size(T2, 1)*.5);          % west-east centerpoints
jcp = round(size(T2, 2)*.5);          % south-north centerpoint

lat = squeeze(XLAT(icp, jcp, :));
lon = squeeze(XLON(icp, jcp, :));
t2 = squeeze(T2(icp, jcp, :));
psfc = squeeze(PSFC(icp, jcp, :));
hgt = squeeze(HGT(icp, jcp, :));
qvapor = squeeze(QVAPOR(icp, jcp, eta, :));
thetap = squeeze(T(icp, jcp, eta, :));
pp = squeeze(P(icp, jcp, eta, :)); 
pb = squeeze(PB(icp, jcp, eta, :));

% Tidy up
clear *cp XLAT XLON T2 PSFC QVAPOR T P PB


% ---- Simulation time; ---------------------------------------------------
timestring = ncread(fileID, 'Times')';
timenum = datenum(timestring)';

% Save extracted data 
save(strcat(outDataPath, 'Simulation_', num2str(dd, '%02d'), ...
            num2str(mm, '%02d'), num2str(yyyy), ...
            'Domain', num2str(domain), '.mat'), ...
            'domain', 'timestring', 'timenum', 'u', 'v', 'w', ...
            'lat', 'lon', 't2', 'psfc', 'hgt', 'qvapor', 'thetap', 'pp', 'pb')

quit
