function OFR2019_1035Figure4
for lw = [1 3]

dvec = 0:100:30000;
Pressure = 1.01e5 + dvec*9.8*2700;

fs = 24;
if lw == 1
	qs = 90e-3;
	figure(1)
	clf; hold on
	subplot(2,2,1)
	hold on
	text(5350,0,'Lake and Marine Sediment','color','g','rotation',-90,'fontsize',8)
	text(5700,0,'Limestone','color','b','rotation',-90,'fontsize',8)
	text(5910,0,'Sandstone','color','r','rotation',-90,'fontsize',8)
	text(6100,0,'Granite','color','m','rotation',-90,'fontsize',8)
	text(6800,0,'Basalt','color','k','rotation',-90,'fontsize',8)
	text(5000,36,'Thin lines indicate moderate geotherm','color','k','fontsize',12)
	text(5000,38,'Thick lines indicate hot geotherm','color','k','fontsize',12)
	text(5000,40,'Dashed lines indicate prediction with Brochers (2005) relations','color','k','fontsize',12)
else
	qs = 120e-3;
end
qm = 20e-3;
k = 2.5;
Ts = 0;
hr = 20e3;
Temp = 273.15 + Ts + qm*dvec/k + (qs - qm)*hr*(1-exp(-dvec/hr))/k;

dvec = dvec/1e3;

Rock(11).MinIndex = [1     5     7    32    55    60    62    66];
Rock(11).VolFrac = [0.3000    0.0500    0.0500    0.1200    0.1200    0.1200    0.1200    0.1200];
Rock(53).MinIndex = [1     6     7    55    60    62];
Rock(53).VolFrac = [0.5500    0.0750    0.0750    0.1000    0.1000    0.1000];
Rock(67).MinIndex = [1     5    51    66];
Rock(67).VolFrac = [0.0500    0.0500    0.8000    0.1000];
Rock(107).MinIndex = [1     5     7    14    15    25];
Rock(107).VolFrac = [0.0175    0.5425    0.1400    0.1000    0.1000    0.1000];
Rock(129).MinIndex = [ 1     5     7    27    28    29];
Rock(129).VolFrac = [0.3500    0.1750    0.3500    0.03125    0.03125    0.0625];
% Following units are in the top 10 of those represented in the NCM geologic framework
% Lake and Marine sediment, 11
Comp.Min = Rock(11).MinIndex;
Comp.Fr = Rock(11).VolFrac;
MV1 = MinVel(Comp,Pressure,Temp);
P1 = Brocher2005(MV1.Vp/1e3);
% Sandstone, 53
Comp.Min = Rock(53).MinIndex;
Comp.Fr = Rock(53).VolFrac;
MV2 = MinVel(Comp,Pressure,Temp);
P2 = Brocher2005(MV2.Vp/1e3);
% Limestone, 67
Comp.Min = Rock(67).MinIndex;
Comp.Fr = Rock(67).VolFrac;
MV3 = MinVel(Comp,Pressure,Temp);
P3 = Brocher2005(MV3.Vp/1e3);
% Basalt, 107
Comp.Min = Rock(107).MinIndex;
Comp.Fr = Rock(107).VolFrac;
MV4 = MinVel(Comp,Pressure,Temp);
P4 = Brocher2005(MV4.Vp/1e3);
% Granite, 129
Comp.Min = Rock(129).MinIndex;
Comp.Fr = Rock(129).VolFrac;
MV5 = MinVel(Comp,Pressure,Temp);
P5 = Brocher2005(MV5.Vp/1e3);

figure(1)
subplot(2,2,1)
hold on
plot(MV1.Vp,dvec,'g','linewidth',lw);
plot(MV2.Vp,dvec,'r','linewidth',lw);
plot(MV3.Vp,dvec,'b','linewidth',lw);
plot(MV5.Vp,dvec,'m','linewidth',lw);
plot(MV4.Vp,dvec,'k','linewidth',lw);
set(gca,'ydir','reverse');
box on
grid on
xlabel('\itV\rm_P (m/s)')
ylabel('Depth (km)')
set(gcf,'renderer','painters')

subplot(2,2,2)
hold on
plot(MV1.Vs,dvec,'g','linewidth',lw);
plot(MV2.Vs,dvec,'r','linewidth',lw);
plot(MV3.Vs,dvec,'b','linewidth',lw);
plot(MV5.Vs,dvec,'m','linewidth',lw);
plot(MV4.Vs,dvec,'k','linewidth',lw);
plot(P1.Vs*1e3,dvec,'g-.','linewidth',lw);
plot(P2.Vs*1e3,dvec,'r-.','linewidth',lw);
plot(P3.Vs*1e3,dvec,'b-.','linewidth',lw);
plot(P5.Vs*1e3,dvec,'m-.','linewidth',lw);
plot(P4.Vs*1e3,dvec,'k-.','linewidth',lw);
set(gca,'ydir','reverse');
box on
grid on
xlabel('\itV\rm_S (m/s)')
ylabel('Depth (km)')
set(gcf,'renderer','painters')

subplot(2,2,3)
hold on
plot(MV1.Vp./MV1.Vs,dvec,'g','linewidth',lw);
plot(P1.VpVs1,dvec,'g-.','linewidth',lw);
plot(MV2.Vp./MV2.Vs,dvec,'r','linewidth',lw);
plot(P2.VpVs1,dvec,'r-.','linewidth',lw);
plot(MV3.Vp./MV3.Vs,dvec,'b','linewidth',lw);
plot(P3.VpVs1,dvec,'b-.','linewidth',lw);
plot(MV5.Vp./MV5.Vs,dvec,'m','linewidth',lw);
plot(P5.VpVs1,dvec,'m-.','linewidth',lw);
plot(MV4.Vp./MV4.Vs,dvec,'k','linewidth',lw);
plot(MV4.Vp./MV4.Vs,dvec,'k','linewidth',lw);
plot(P4.VpVs1,dvec,'k-.','linewidth',lw);
set(gca,'ydir','reverse');
box on
grid on
xlabel('\itV\rm_P/\itV\rm_S ratio')
ylabel('Depth (km)')
set(gcf,'renderer','painters')

subplot(2,2,4)
hold on
plot(MV1.p,dvec,'g','linewidth',lw);
plot(P1.Den*1e3,dvec,'g-.','linewidth',lw);
plot(MV2.p,dvec,'r','linewidth',lw);
plot(P2.Den*1e3,dvec,'r-.','linewidth',lw);
plot(MV3.p,dvec,'b','linewidth',lw);
plot(P3.Den*1e3,dvec,'b-.','linewidth',lw);
plot(MV5.p,dvec,'m','linewidth',lw);
plot(P5.Den*1e3,dvec,'m-.','linewidth',lw);
plot(MV4.p,dvec,'k','linewidth',lw);
plot(P4.Den*1e3,dvec,'k-.','linewidth',lw);
set(gca,'ydir','reverse');
box on
grid on
xlabel('Density (kg/m^3)')
ylabel('Depth (km)')
set(gcf,'renderer','painters')

end

function P = Brocher2005(Vp)
%function P = Brocher2005(Vp)
%
%       Vp in km/s
%
%       Outut in km/s and gm/cc
%
% Brocher's 2005 relations
% Brocher, T.M., 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust: Bulletin of the Seismological Society of America, v. 95, no. 6, p. 2081–2092.
%

%Vp in km/2; Den in g/cc
% For Vp between 1.5 and 8 km/s
P.Den = 1.6612*Vp - 0.4721*Vp.^2 + 0.0671*Vp.^3 - 0.0043*Vp.^4 + 0.000106*Vp.^5;
P.Vs = 0.7858 - 1.2344.*Vp + 0.7949*Vp.^2 - 0.1238*Vp.^3 + 0.0064*Vp.^4;
P.v = 0.8835 - 0.315*Vp + 0.0491*Vp.^2 - 0.0024*Vp.^3;
P.VpVs1 = Vp./P.Vs;
P.VpVs2 = sqrt((2*P.v - 2)./(2*P.v - 1));
