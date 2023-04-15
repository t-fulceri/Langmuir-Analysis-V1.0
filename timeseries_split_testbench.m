f = figure;
movegui(f,'center');
set(gcf,'Renderer','painters');
plot(voltage)

%Calculate resolution of voltage
dV = diff(voltage);
delta_voltage = min(abs(dV(dV~=0)));

ind_zerocross = find(abs(voltage(500:end)) < delta_voltage);

T_vec = diff(ind_zerocross);

f = figure;
movegui(f,'center');
set(gcf,'Renderer','painters');
hist(T_vec)

test = 1;