cd("Function")
try
mex iauPnm06a_mex.c pnm06a.c fw2m.c nut06a.c pfw06.c rx.c ry.c rz.c ir.c nut00a.c fae03.c faf03.c faju03.c fal03.c fama03.c fame03.c faom03.c fapa03.c fasa03.c faur03.c fave03.c obl06.c
disp('iauPnm06a_mex build sucess');
mex iauGst06_mex.c gst06.c bpn2xy.c eors.c s06.c era00.c anp.c fad03.c fae03.c faf03.c fal03.c falp03.c faom03.c fapa03.c fave03.c
disp('iauGst06_mex build sucess');
mex iauSp00_mex.c sp00.c
disp('iauSp00_mex build sucess');
mex iauPom00_mex.c pom00.c rx.c ry.c rz.c ir.c
disp('iauPom00_mex build sucess');
mex iauGmst06_mex.c gmst06.c era00.c anp.c
disp('iauGmst06_mex build sucess');
mex Legendre_mex.c Legendre.c
disp('Legendre_mex build sucess');
catch ME
    disp(['错误消息：', ME.message]);
end
cd("..")