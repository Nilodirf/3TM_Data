function [c_p,T] = c_p_from_vDOS_yww(vDOS,temperature)
% function calculates molar phonon heat capacity from vDOS
%inputs: vDOS in meV (2D array), temperatures in K (1D array)
%c_p is in J/(UC*K)

%calculate energy as function of temperature (in meV)
for t=1:length(temperature)
    %in total
    thermal_occupation = (vDOS(:,2))./(exp((vDOS(:,1)./(8.617e-2*temperature(t))))-1);
    thermal_occupation(isnan(thermal_occupation)) = 0;
    energy(t)=trapz(vDOS(:,1),thermal_occupation.*vDOS(:,1));
end

%convert energy from meV to J
e=1.60217e-19;  % elementary charge
energy=energy*1e-3*e;

%step size (for differentiation)
c_p = diff(energy)./diff(temperature);
T = temperature(1:end-1) + diff(temperature)/2;


end

