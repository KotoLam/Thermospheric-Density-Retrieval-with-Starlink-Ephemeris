function V_Moon = Potential_Moon(r_Sat,r_Moon)
r_Sat = r_Sat(:);
r_Moon = r_Moon(:);
GM_Earth   = 0.3986004415E+15;         	 % [m^3/s^2]; EGM2008
GM_Moon    = GM_Earth/81.3005682214972154; % [m^3/s^2]; DE440
s = norm(r_Sat - r_Moon);   % 卫星到月亮距离
V_Moon = GM_Moon*(1/s-1/norm(r_Moon)-dot(r_Sat,r_Moon)/norm(r_Moon)^3);

end