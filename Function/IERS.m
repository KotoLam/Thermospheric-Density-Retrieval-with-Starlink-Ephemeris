%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%  
% Last modified:   2018/02/01   Meysam Mahooti
% 
% Modified: linear interpolation only
%--------------------------------------------------------------------------

function [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] ...
    = IERS(eopdata,mjd_UTC)

Arcs = 3600*180/pi;         % Arcseconds per radian

% linear interpolation
mjd = floor(mjd_UTC);
idx = find(mjd==eopdata(4,:),1,'first');
preeop = eopdata(:,idx);
nexteop = eopdata(:,idx+1);
fixf = mjd_UTC-floor(mjd_UTC);
% Setting of IERS Earth rotation parameters
% (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
TAI_UTC = preeop(13);         % TAI-UTC time difference [s]

x_pole  = x_pole/Arcs;  % Pole coordinate [rad]
y_pole  = y_pole/Arcs;  % Pole coordinate [rad]
dpsi    = dpsi/Arcs;    % [rad]
deps    = deps/Arcs;    % [rad]
dx_pole = dx_pole/Arcs; % Pole coordinate [rad]
dy_pole = dy_pole/Arcs; % Pole coordinate [rad]


