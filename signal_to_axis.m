function [ axis ] = signal_to_axis( signal )
%SIGNAL_TO_AXIS Returns an axis of values found in the input signal
%   Warning: the input sweep should come from an ADC, otherwise
%   the function will not find the discrete steps

%Find the minimum step
[ sweep_step ] = signal_to_step( signal );
%Find min and max
sweep_max = max(signal);
sweep_min = min(signal);
%Find the proper amount of elements
Nsweep = (sweep_max-sweep_min)/sweep_step + 1;
%Re-create the axis
axis = linspace(sweep_min,sweep_max,Nsweep);

end

