function V_Sun = Potential_Sun(r_Sat,r_Sun)
r_Sat = r_Sat(:);
r_Sun = r_Sun(:);
GM_Sun     = 132712440041.279419e9; 			    % [m^3/s^2]; DE440
s = norm(r_Sat - r_Sun);   % 卫星到太阳距离
V_Sun = GM_Sun*(1/s-1/norm(r_Sun)-dot(r_Sat,r_Sun)/norm(r_Sun)^3);

end