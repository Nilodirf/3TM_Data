% Script to extract the phonon DOS from literature (based on a script by
% Daniela Zahn).

% Controls: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x_scale = 18;             % THz
THz_to_meV_factor = 4.13567; % multiply THz to get meV
x_range = [2.65 3.8];     % THz
num_of_low_energy_pts = 100;
temperatures = [0:10:290 linspace(296,315,10000) 310:10:1000]; % kelvin
cutoff_opt_acc = 44;      % meV
% material properties:
atoms_per_unit_cell = 8;  % 4 x Ni and 4 x O
density = 6.67;           % g/cm^3      <-fits to the uc with 4+4 ions
   density = density*1000;% kg/m^3 
   Ni_mass = 58.6934;     % unit: u
   O_mass  = 15.999 ;     % unit: u
   u = 1.66e-27;          % atomic mass unit in kg
M = u*(Ni_mass+O_mass)*4; % mass of the UC in kg
N_a=6.0221409e+23;        % Avogadros number
UC_to_volume =  density / M;  % conversion factor
UC_to_mol = 2*N_a / atoms_per_unit_cell; % makes it fit to the Coy 1976 paper
% load the clean PNG file: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fullFilePath = 'pDOS_clean.png'; % the original image has x axis in units of THz
[x,y] = png_to_data(fullFilePath,x_scale,0.01824,1,0);
y(isnan(y)) = 0;

% fit the lower end to zero with power law: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inds = yww_indsBetween(x,x_range(1),x_range(2)); % indexes to fit
ft = fit(x(inds)',y(inds)','a*x^b');
x_low = 0:x_range(1)/num_of_low_energy_pts:x_range(1) ; % x points to simulate directly below the fit range.
y_low = ft.a*x_low.^ft.b;

% plot this
yww_freshFigure(1,[]); 
plot(x,y,'displayname','data'); 
plot(x_low  ,y_low              ,'-r','displayname','simulated region'); 
plot(x(inds),ft.a.*x(inds).^ft.b,'-r','linewidth',3,'displayname',['fit, power: ' num2str(ft.b,'%1.2f')  ]);
xlabel('THz'); xlim([0 18]); legend('show');

%% Final DOS arrangement: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% make new variables instead of x and y:
I1 = min(find(inds == 1)); % The first point in the "good" region.
X = [x_low x(I1:end)];
Y = [y_low y(I1:end)];

% fix axes to proper units: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X = X .*THz_to_meV_factor; % convert the x axis from THz to meV
normalizationFactor = trapz(X,Y);
Y = Y ./ normalizationFactor * 3 * atoms_per_unit_cell;

% plot this too
yww_freshFigure(2,[]); 
plot(X,Y,'displayname','data'); 
xlabel('meV'); xlim([0 18*THz_to_meV_factor]); 
ylabel('DOS (states / meV / unit cell)')
title('NiO Phonon density of states');% legend('show');

% save this daddy
XY = [X' Y'];
save('NiO_phonon_DOS_meV.mat','XY');

%% Lattice heat capacity: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[c_p,T] = c_p_from_vDOS_yww(XY,temperatures); % J/UC/K
% c_p = c_p * UC_to_mol ; % convert from per mol to per m^3
c_p = c_p * UC_to_volume ; % convert from per mol to per m^3

yww_freshFigure(3,[]); 
plot(T,c_p,'.','displayname','Lattice heat capacity'); 
xlabel('T (K)'); %xlim([0 18*THz_to_meV_factor]); 
ylabel('Lattice heat capacity (J/m^3K)')
title('Lattice heat capacity');% legend('show');

% next step: separate heat capacities for optical and acoutsic phonons
% acoustic phonons
X1 = X(X < cutoff_opt_acc);  Y1 = Y(X < cutoff_opt_acc); 
[c_p1,T1] = c_p_from_vDOS_yww([X1' Y1'],temperatures); % J/mol
c_p1 = c_p1 * UC_to_volume ; % convert from per mol to per m^3
plot(T1,c_p1,'.','displayname','Acoustic phonons heat capacity'); 

% optical phonons
X2 = X(X > cutoff_opt_acc);  Y2 = Y(X > cutoff_opt_acc); 
[c_p2,T2] = c_p_from_vDOS_yww([X2' Y2'],temperatures); % J/mol
c_p2 = c_p2 * UC_to_volume ; % convert from per mol to per K*m^3
plot(T2,c_p2,'.','displayname','Optical phonons heat capacity'); 

legend('show','location','southeast')

total_lattice = [T' c_p']; acoustic = [T1' c_p1']; optical = [T2' c_p2'];
save('NiO_heat_capacities.mat','total_lattice','acoustic','optical');
