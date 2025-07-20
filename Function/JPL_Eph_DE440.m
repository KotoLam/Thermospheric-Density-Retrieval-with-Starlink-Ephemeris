%--------------------------------------------------------------------------
%
% JPL_Eph_DE440: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Input:
%   JD_TDB         Julian date of TDB
%
% Output:
%   r_Earth and r_SunSSB (solar system barycenter (SSB)), r_Mercury, r_Venus,
%   r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon,
%   r_Sun (geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Last modified:   2023/12/28   Meysam Mahooti
% 
% Modified: Only need Moon and Sun position
%--------------------------------------------------------------------------
function [r_Moon,r_Sun] = JPL_Eph_DE440(JD_TDB,DE440Coeff)


i = find(DE440Coeff(:,1)<=JD_TDB & JD_TDB<=DE440Coeff(:,2),1,'first'); 
PCtemp = DE440Coeff(i,:);

t1 = PCtemp(1); % JD_TDB at start of interval

dt = JD_TDB - t1;

temp = (231:13:270);
Cx_EarthMoon = PCtemp(temp(1):temp(2)-1);
Cy_EarthMoon = PCtemp(temp(2):temp(3)-1);
Cz_EarthMoon = PCtemp(temp(3):temp(4)-1);
temp = temp+39;
Cx = PCtemp(temp(1):temp(2)-1);
Cy = PCtemp(temp(2):temp(3)-1);
Cz = PCtemp(temp(3):temp(4)-1);
Cx_EarthMoon = [Cx_EarthMoon,Cx];
Cy_EarthMoon = [Cy_EarthMoon,Cy];
Cz_EarthMoon = [Cz_EarthMoon,Cz];    
if (0<=dt && dt<=16)
    j=0;
    JD0 = t1;
elseif(16<dt && dt<=32)
    j=1;
    JD0 = t1+16*j;
end
r_EarthMoon = 1e3*Cheb3D(JD_TDB, 13, JD0, JD0+16, Cx_EarthMoon(13*j+1:13*j+13),...
                     Cy_EarthMoon(13*j+1:13*j+13), Cz_EarthMoon(13*j+1:13*j+13))';

temp = (441:13:480);
Cx_Moon = PCtemp(temp(1):temp(2)-1);
Cy_Moon = PCtemp(temp(2):temp(3)-1);
Cz_Moon = PCtemp(temp(3):temp(4)-1);
for i=1:7
    temp = temp+39;
    Cx = PCtemp(temp(1):temp(2)-1);
    Cy = PCtemp(temp(2):temp(3)-1);
    Cz = PCtemp(temp(3):temp(4)-1);   
    Cx_Moon = [Cx_Moon,Cx];
    Cy_Moon = [Cy_Moon,Cy];
    Cz_Moon = [Cz_Moon,Cz];    
end
if (0<=dt && dt<=4)
    j=0;
    JD0 = t1;
elseif(4<dt && dt<=8)
    j=1;
    JD0 = t1+4*j;
elseif(8<dt && dt<=12)
    j=2;
    JD0 = t1+4*j;
elseif(12<dt && dt<=16)
    j=3;
    JD0 = t1+4*j;
elseif(16<dt && dt<=20)
    j=4;
    JD0 = t1+4*j;
elseif(20<dt && dt<=24)
    j=5;
    JD0 = t1+4*j;
elseif(24<dt && dt<=28)
    j=6;
    JD0 = t1+4*j;
elseif(28<dt && dt<=32)
    j=7;
    JD0 = t1+4*j;
end
r_Moon = 1e3*Cheb3D(JD_TDB, 13, JD0, JD0+4, Cx_Moon(13*j+1:13*j+13),...
                    Cy_Moon(13*j+1:13*j+13), Cz_Moon(13*j+1:13*j+13))';

temp = (753:11:786);
Cx_Sun = PCtemp(temp(1):temp(2)-1);
Cy_Sun = PCtemp(temp(2):temp(3)-1);
Cz_Sun = PCtemp(temp(3):temp(4)-1);
temp = temp+33;
Cx = PCtemp(temp(1):temp(2)-1);
Cy = PCtemp(temp(2):temp(3)-1);
Cz = PCtemp(temp(3):temp(4)-1);   
Cx_Sun = [Cx_Sun,Cx];
Cy_Sun = [Cy_Sun,Cy];
Cz_Sun = [Cz_Sun,Cz];
if (0<=dt && dt<=16)
    j=0;
    JD0 = t1;
elseif(16<dt && dt<=32)
    j=1;
    JD0 = t1+16*j;
end
r_Sun = 1e3*Cheb3D(JD_TDB, 11, JD0, JD0+16, Cx_Sun(11*j+1:11*j+11),...
                   Cy_Sun(11*j+1:11*j+11), Cz_Sun(11*j+1:11*j+11))';

EMRAT = 81.3005682214972154; % DE440
EMRAT1 = 1/(1+EMRAT);
r_Earth = r_EarthMoon-EMRAT1*r_Moon;

r_Sun = -r_Earth+r_Sun;

end

