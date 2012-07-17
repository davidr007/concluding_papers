function [energy] = compute_energy(event,time)

%Inputs: 
% event     [vector 1Xn] seismogram
% time      [vector 1Xn] time
%
% OUTPUTS:
% energy    [scalar] = int(dt x energy.^2, time(1),time(end))

%error('BE CAREFUL - I Think you should normalise by time(end)-time(1)')
energy = trapz(time,event.^2);
