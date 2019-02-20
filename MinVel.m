function MV = MinVel(Comp, P, T)
%function MV = MinVel(Comp, Pressure, Temperature)
%
%	Written by Oliver Boyd
%	V0.1 released Nov. 16th 2007, last updated 12/11/2018
%	Mineral database updated 12/10/2018
%
%	This program calculates the velocity of the mineral assemblage 
%	given by Comp at a given pressure and temperature (which may be
%	vectors). It returns the geophysical parameters in the
%	structure MV. The velocities are expressed as Voigt, Reuss,
%       and Voigt-Reuss-Hill averages.
%
%	The data required for this analysis is taken from Hacker and Abers (2003),
%	updated by Abers and Hacker in 2016, and expanded by Boyd in 2018.
%	The moduli at pressure and temperature are calculated based on the
%	procedures of Hacker and Abers (2004), Bina and Helffrich (1992) and
%	Holland and Powell (1998) as outlined in the supplementary section of 
%	Boyd et al. (2004) with updates by Abers and Hacker (2016) for quartz.
%
%	OUTPUT (SI Units)
%		MV - structure containing
%			Vpv - P-wave velocity, Voigt average
%			Vpr - P-wave velocity, Reuss average
%			Vsv - S-wave velocity, Voigt average
%			Vsr - S-wave velocity, Reuss average
%			p   - Density
%			a   - Thermal Expansivity
%			Kv  - Bulk modulus, Voigt average
%			Gv  - Shear modulus, Voigt average
%			Kr  - Bulk modulus, Reuss average
%			Gr  - Shear modulus, Reuss average
%			Vp  - P-wave velocity, Voigt-Reuss-Hill average
%			Vs  - S-wave velocity, Voigt-Reuss-Hill average
%			K   - Bulk modulus, Voigt-Reuss-Hill average
%			G   - Shear modulus, Voigt-Reuss-Hill average
%
%	Additional outputs based on Voigt-Reuss-Hill average
%			l   - Lambda
%			v   - Poisson's ratio
%			E   - Young's modulus
%
%	INPUTS
%	Pressure - desired pressure or vector of pressures (Pa)
%	Temperature - desired temperature or vector of temperatures (K)
%
%	Comp - a composition structure with the following fields: 
%		Min - The mineral index vector.
%		Fr - Volume fraction for each mineral in Min (0 to 1).
%
%	    Quartz
%		 1. Alpha Quartz                  
%		 2. Beta Quartz                   
%		 3. Coesite                       
%   	    Feldspar group
%             Plagioclase
%		 4. High Albite                   
%		 5. Low Albite                    
%		 6. Anorthite                     
%
%		 7. Orthoclase                    
%		 8. Sanidine                      
%	    Garnet structural group
%		 9. Almandine                     
%		10. Grossular                     
%		11. Pyrope                        
%	    Olivine group
%		12. Forsterite                    
%		13. Fayalite                      
%	    Pyroxene group
%		14. Diopside                      
%		15. Enstatite                     
%		16. Ferrosilite                   
%		79. Mg-Tschermak                  
%		17. Jadeite                       
%		18. Hedenbergite                  
%		80. Acmite                        
%		81. Ca-Tschermak                  
%	    Amphibole supergroup
%		19. Glaucophane                   
%		20. Ferroglaucophane              
%		21. Tremolite                     
%		22. Ferroactinolite               
%		23. Tshermakite                   
%		24. Pargasite                     
%		25. Hornblende                    
%		26. Anthophyllite                 
%	    Mica group
%		27. Phlogopite                    
%		28. Annite                        
%		29. Muscovite                     
%		30. Celadonite                    
%	    Other
%		31. Talc                          
%		32. Clinochlore                   
%		33. Daphnite                      
%		34. Antigorite                    
%		35. Zoisite                       
%		36. Clinozoisite                  
%		37. Epidote                       
%		38. Lawsonite                     
%		39. Prehnite                      
%		40. Pumpellyite                   
%		41. Laumontite                    
%		42. Wairakite                     
%		43. Brucite                       
%		44. Clinohumite                   
%		45. Phase A                       
%		46. Sillimanite                   
%		47. Kyanite                       
%		48. Spinel                        
%		49. Hercynite                     
%		50. Magnetite                     
%		51. Calcite                       
%		52. Aragonite                     
%		82. Magnesite                     
%		83. En79Fs09Ts12                  
%		84. Di75He9Jd3Ts12                
%		85. ilmenite                      
%		86. cordierite                    
%		87. scapolite (meionite)          
%		88. rutile                        
%		89. sphene                        
%		53. Corundum                      
%		54. Dolomite                      
%		74. Halite                        
%		77. Pyrite                        
%		78. Gypsum   
%		90. Anhydrite   
%		 0. Water   
%		-1. Ice   
%	    Clays
%		55. Montmorillonite (Saz-1)
%		56. Montmorillonite (S Wy-2)
%		57. Montmorillonite (STX-1)
%		58. Montmorillonite (S Wy-1)
%		59. Montmorillonite (Shca-1)
%		60. Kaolinite (Kga-2)
%		61. Kaolinite (Kga-1b)
%		62. Illite (IMT-2)
%		63. Illite (ISMT-2)
%		66. Smectite (S Wa-1)
%		70. Montmorillonite (S YN-1)
%		71. Chrysotile                    
%		72. Lizardite                     
%		76. Dickite                       
%

load('MineralPhysicsDatabase_181107');

if abs(sum(Comp.Fr) - 1) > 1e-6
	disp('Composition does not sum to one. - Exiting')
	Comp
	fprintf('Composition sums to: %f\n',sum(Comp.Fr));
	MV = [];
	return
end
nPar = length(Par(1,:));
nMin = length(Par(:,1));
nPT = length(P);
p = zeros(nMin,nPT);
Vp = zeros(nMin,nPT);
Vs = zeros(nMin,nPT);
a = zeros(nMin,nPT);

MinIndex = Par(:,1);
gind = Comp.Min;
nMin = length(gind);
for j = 1:nMin
	k = find(MinIndex == Comp.Min(j));
	cmp(j) = Comp.Fr(j);
	if Comp.Min(j) == 1 % Quartz
		k = [1 2];
		[p(j,:), Vp(j,:), Vs(j,:), a(j,:), K(j,:), G(j,:), V(j,:)] = GetV(P,T, ...
			Par(k,23), Par(k,25), Par(k,19), ...
			Par(k,6), Par(k,9), ...
			Par(k,16), Par(k,21), ...
			Par(k,11), Par(k,14), Par(k,4));
	else
		[p(j,:), Vp(j,:), Vs(j,:), a(j,:), K(j,:), G(j,:), V(j,:)] = GetV(P,T, ...
			Par(k,23), Par(k,25), Par(k,19), ...
			Par(k,6), Par(k,9), ...
			Par(k,16), Par(k,21), ...
			Par(k,11), Par(k,14), Par(k,4));
	end
end

MV = GetMV(p, K, G, Vp, Vs, a, cmp);

return

function [pTP, Vp, Vs, aTP, Ks, GTP, VTP] = GetV(P, T, th, dt, gamma, p, a, G, dGdPr, K, dKdPr, Vo)

% GetV returns the density, p, P and S velocity, Vp and Vs, and thermal
% expansivity, a, at a given and temperature and pressure, T and P. The postscript
% TP stands for temperature and pressure.
% It requires the following parameters:
% 
%	P - Pressure
%	T - Temperature
%	a0 - Thermal expansivity
%	p - Density
%	G - Shear Modulus
%	K - Bulk Modulus
%	d.dP - Derivative with respect to pressure
%	dt = isothermal second gruneisen parameter
%	th - first gruneisen parameter
%	gamma - dlnGdp at constant pressure
%	Vo - Volume at STP
%
% Equations based on Hacker and Abers (2004) with modifications to equations in GetF by Boyd et al (2004).
% Further modifcation for quartz adapted from Abers and Hacker (2016)
%

if length(K) == 2 % Dealing with quartz
	Tab = 273 + 574.3 + 0.2559*P/1e6 - 6.406e-6*(P/1e6).^2;

	dt = dt(1,:)*ones(size(T));
	gamma = gamma(1,:)*ones(size(T));

	j = find(T <= Tab); % alpha
	if length(j) > 0
		%%%% Expansion fit: form of Carpenter et al. 1998
		Ttr = Tab(j);
		Tc = Ttr - 3.1;    % values from Carpenter et al.
		Vso = -5.4E-3;   % Slightly adjusted from -5.1E-3, to better fit Ohno density data
		Vs = (2/3)*Vso*(1 + sqrt(1 - 0.75*(T(j) - Tc)./(Ttr - Tc)));  % Carpenter et al eq 25
		Vsstp = (2/3)*Vso*(1 + sqrt(1 - 0.75*(298 - Tc)./(Ttr - Tc)));
		phi(j) = log(1 + Vs) - log(1 + Vsstp);  % integrated therm. expans.;
		% Derivative gives alpha(T) thermal expansion
		aT(j) = -0.25*Vso./((1 + Vs).*(Ttr - Tc).*sqrt(1 - 0.75*(T(j) - Tc)./(Ttr - Tc)));
	
		pT(j) = p(1)*exp(-phi(j));
		VT(j) = Vo(1)*exp(phi(j));

 		% regressions by Abersw and Hacker (2016) on Ohno et al. 2006 data
		aktp0 = (11.384 + 0.6478*log(Tab(j) - 298) + 0.5878*(log(Tab(j) - 298).^2)) ; %Ks(P,25C)
		amup0 = (41.721 - 0.034879*log(Tab(j) - 298) + 0.08424*(log(Tab(j) - 298).^2)); %G(P,25C)
		KT(j) = (K(1)./aktp0).*(11.384 + 0.6478*log(Tab(j) - T(j)) + 0.5878*(log(Tab(j) - T(j)).^2)) ; %Ks(0,T)
		GT(j) = (G(1)./amup0).*(41.721 - 0.034879*log(Tab(j) - T(j)) + 0.08424*(log(Tab(j) - T(j)).^2)); % G(0,T)

		dKdP(j) = dKdPr(1);
		dGdP(j) = dGdPr(1);
		f(j) = GetF(P(j),KT(j),dKdPr(1));
		thg(j) = th(1);
		p0(j) = p(1);
		pP(j) = p(1)*(1 + 2*f(j)).^(3/2);
	end
	j = find(T > Tab); % beta
	if length(j) > 0
		% I was uncomfortable assuming phi = 0 and aT = 0. I therefor maintained the original formulation
		% but normalized dt and gamma by 1/phi to get the same fit.

		% This is the orginal method in Abers and Hacker (2016)
		%dt(j) = -0.0963*log(T(j) - Tab(j)) + 0.0072; %fit to Ohno et al. (2006) & Ackerman et al. (1974)
		%gamma(j) = 0.0002 - 0.0002*(T(j) - Tab(j)); %fit to Ohno et al. (2006)   dGdT in Excel macro
		%phi(j) = ones(1,length(j));
		%aT(j) = zeros(1,length(j));
		%pT(j) = p(2);
		%VT(j) = Vo(2);

		phi(j) = a(2)*((T(j)-298) - 20*(sqrt(T(j)) - sqrt(298)));
		aT(j) = a(2)*(1 - 10./sqrt(T(j)));
		dt(j) = (1./phi(j)).*(-0.0963*log(T(j) - Tab(j)) + 0.0072); %fit to Ohno et al. (2006) & Ackerman et al. (1974)  
		gamma(j) = (1./phi(j)).*(0.0002 - 0.0002*(T(j) - Tab(j))); %fit to Ohno et al. (2006) 
		pT(j) = p(2)*exp(-phi(j));
		VT(j) = Vo(2)*exp(phi(j));
		KT(j) = K(2)*exp(-dt(j).*phi(j));
		GT(j) = G(2)*exp(-gamma(j).*phi(j));
		dKdP(j) = dKdPr(2);
		dGdP(j) = dGdPr(2);
		f(j) = GetF(P(j),KT(j),dKdPr(2));
		thg(j) = th(2);
		p0(j) = p(2);
		pP(j) = p(2)*(1 + 2*f(j)).^(3/2);
	end
else
	phi = a*((T-298) - 20*(sqrt(T) - sqrt(298)));
	aT = a*(1 - 10./sqrt(T));
	pT = p*exp(-phi);
	VT = Vo*exp(phi);
	KT = K*exp(-dt*phi);
	GT = G*exp(-gamma*phi);
	dKdP = dKdPr;
	dGdP = dGdPr;
	f = GetF(P,KT,dKdP);
	thg = th;
	p0 = p;
	pP = p*(1 + 2*f).^(3/2);
end

VTP = VT.*(1 + 2*f).^(-3/2);

%KTP = KT.*((1+2*f).^(5/2)).*(1 - f.*(5 - 3*dKdP));
%GTP = GT.*((1+2*f).^(5/2)).*(1 - f.*(5 - 3*dGdP*KT./GT));

%Assume ddKddP is 1/KT and ddGddP is 1/GT as suggested by Bina and Helffrich, 1992
KTP = KT.*((1+2*f).^(5/2)).*(1 - f.*(5 - 3*dKdP) + (f.^2).*(9 + (3*dKdP - 7).*(3*dKdP - 5))/2);
if GT > 0
	GTP = GT.*((1+2*f).^(5/2)).*(1 - f.*(5 - 3*dGdP.*KT./GT) + (f.^2).*(9*KT.^2./GT.^2 + 9*(dKdP - 4).*dGdP.*KT./GT + 35)/2);
else
	GTP = 0;
end

aTP = aT.*((pP./p0).^-dt);
pTP = pP.*pT./p0;
Ks = KTP.*(1 + T.*thg.*aTP);

Vp = sqrt((Ks + 4*GTP/3)./pTP);
Vs = sqrt(GTP./pTP);

return

function f = GetF(P,KT,dKdP)

% From Bina and Helffrich, 1992. Assume second derivative of K with respect to P,
% ddKddP, is 1/dKdP as suggested by Bina and Helffrich, else
%
% y = (9*dKdP*ddKddP + 4*e*(4 - 3*dKdP) + 5*(3*dKdP - 5))/6;
%
% Took Eq from Bina and Helfrich, Turned into polynomial, and used matlab roots
% command to get the roots of the polynomial. I use solution closest to First
% order solution.
%

e = 0.75*(4-dKdP);
y = (9 + 4*e*(4 - 3*dKdP) + 5*(3*dKdP - 5))/6;

y2 = y^2;
e2 = e^2;
ey = e*y;

C(1) = 32*y2;
C(2) = 80*y2 - 128*ey;
C(3) = 80*y2 - 320*ey + 64*y + 128*e2;
C(4) = 40*y2 - 320*ey + 160*y + 320*e2 - 128*e;
C(5) = 10*y2 - 160*ey + 160*y + 320*e2 - 320*e + 32;
C(6) = y2 - 40*ey + 80*y + 160*e2 - 320*e + 80;
C(7) = -4*ey + 20*y + 40*e2 - 160*e + 80;
C(8) = 2*y + 4*e2 - 40*e + 40;
C(9) = 10 - 4*e;
C(10) = 1;
C(11) = 0;
for j = 1:length(P)
	C(12) = -((P(j)/KT(j))^2)/9;

	F1 = sqrt(-C(12));
	F = roots(C);
	k = find(imag(F) == 0 & abs(F-F1) == min(abs(F-F1)));
	f(j) = F(k(1));
end
return

function MV = GetMV(p, K, G, Vp, Vs, a, cmp)

% GetMV calculates the Voigt-Reuss averages (of velocities but should be moduli
% though VRH will be about the same) of the s and p-wave velocity (Anderson 1997).

MV.Vpv = zeros(1,length(Vp(1,:)));
MV.Vpr = zeros(1,length(Vp(1,:)));
MV.Vsv = zeros(1,length(Vp(1,:)));
MV.Vsr = zeros(1,length(Vp(1,:)));
MV.p = zeros(1,length(Vp(1,:)));
MV.a = zeros(1,length(Vp(1,:)));

MV.p = cmp*p(1:length(cmp),:);
MV.a = cmp*a(1:length(cmp),:);

% Voigt Average
%MV.Vpv = (Wt.*cmp)*Vp(1:length(cmp),:)/WT;
%MV.Vsv = (Wt.*cmp)*Vs(1:length(cmp),:)/WT;
MV.Kv = cmp*K(1:length(cmp),:);
MV.Gv = cmp*G(1:length(cmp),:);
MV.Vpv = sqrt((MV.Kv + 4*MV.Gv/3)./MV.p);
MV.Vsv = sqrt(MV.Gv./MV.p);

% Reuss Average
%MV.Vpr = WT./((Wt.*cmp)*(1./Vp(1:length(cmp),:)));
%MV.Vsr = WT./((Wt.*cmp)*(1./Vs(1:length(cmp),:)));
MV.Kr = 1./(cmp*(1./K(1:length(cmp),:)));
MV.Gr = 1./(cmp*(1./G(1:length(cmp),:)));
MV.Vpr = sqrt((MV.Kr + 4*MV.Gr/3)./MV.p);
MV.Vsr = sqrt(MV.Gr./MV.p);

% Voigt Reuss Hill Average
MV.Vp = (MV.Vpv + MV.Vpr)/2;
MV.Vs = (MV.Vsv + MV.Vsr)/2;
MV.K = (MV.Kv + MV.Kr)/2;
MV.G = (MV.Gv + MV.Gr)/2;

MV.l = MV.K-2*MV.G/3;
MV.v = MV.l./(2*(MV.l+MV.G));
MV.E = 2*MV.G.*(1+MV.v);
return
