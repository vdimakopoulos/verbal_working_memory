function plot_leadfield(elec_realigned,sourcemodel_and_leadfield,headmodel,mesh_brain)
figure('units', 'normalized', 'outerposition', [0 0 0.5 0.5])
source_index = 1; %% a superficial sources
sensory_dipole_current = 100e-9; % Am (realistic)

n_sensors = length(elec_realigned.label);

inside_sources = find(sourcemodel_and_leadfield.inside);
inside_index = inside_sources(source_index);
lead = sourcemodel_and_leadfield.leadfield{inside_index};
xs = zeros(1, n_sensors);
ys = zeros(1, n_sensors);
zs = zeros(1, n_sensors);
voltages = zeros(1, n_sensors);
titles = {'Lead field (x)' 'Lead field (y)' 'Lead field (z)'};

% get the xyz and norm

for sensor_index = 1:n_sensors
    this_x = lead(sensor_index, 1);
    this_y = lead(sensor_index, 2);
    this_z = lead(sensor_index, 3);
    this_norm = norm(lead(sensor_index, :));
    xs(sensor_index) = this_x * sensory_dipole_current;
    ys(sensor_index) = this_y * sensory_dipole_current;
    zs(sensor_index) = this_z * sensory_dipole_current;
    voltages(sensor_index) = this_norm * sensory_dipole_current;
end

% plot xyz
axes = {xs ys zs};

for axis_index = 1:3
    this_axis = axes{axis_index};
    subplot(1, 3, axis_index)
    hold on
    ft_plot_topo3d(elec_realigned.chanpos, this_axis, 'facealpha', 0.8)
    if strcmp(headmodel.type, 'dipoli')
        caxis([-10e-6, 10e-6])
    end
    c = colorbar('location', 'southoutside');
    c.Label.String = 'Lead field (V)';
    axis tight
    ft_plot_mesh(mesh_brain, 'facealpha', 0.10);
    ft_plot_sens(elec_realigned, 'elecsize', 20);
    title(titles{axis_index})
    plot3(sourcemodel_and_leadfield.pos(inside_index, 1), ...
        sourcemodel_and_leadfield.pos(inside_index, 2), ...
        sourcemodel_and_leadfield.pos(inside_index, 3), 'bo', ...
        'markersize', 20, 'markerfacecolor', 'r')
end

% plot norm

figure('units', 'normalized', 'outerposition', [0 0 0.5 0.85])
hold on
ft_plot_topo3d(elec_realigned.chanpos, voltages, 'facealpha', 0.8)
if strcmp(headmodel.type, 'dipoli')
    caxis([0, 10e-6])
end
c = colorbar('location', 'eastoutside');
c.Label.String = 'Lead field (V)';
axis tight
ft_plot_mesh(mesh_brain, 'facealpha', 0.10);
ft_plot_sens(elec_realigned, 'elecsize', 20);
title('Leadfield magnitude')
plot3(sourcemodel_and_leadfield.pos(inside_index, 1), ...
    sourcemodel_and_leadfield.pos(inside_index, 2), ...
    sourcemodel_and_leadfield.pos(inside_index, 3), 'bo', ...
    'markersize', 20, 'markerfacecolor', 'r')

view(-90, 0)



figure
hold on
ft_plot_sens(elec_realigned, 'elecsize', 40);
ft_plot_headmodel(headmodel, 'facealpha', 0.5);
view(90, 0)

source_index = 2;
sensory_dipole_current = 1e-1; % Am (unrealistic)
n_sensors = length(elec_realigned.label);

inside_sources = find(sourcemodel_and_leadfield.inside);
inside_index = inside_sources(source_index);
lead = sourcemodel_and_leadfield.leadfield{inside_index};

for sensor_index = 1:n_sensors
    this_pos = elec_realigned.chanpos(sensor_index, :);
    this_lead = lead(sensor_index, :);
    this_voltage = this_lead * sensory_dipole_current;
    quiver3(this_pos(1), this_pos(2), this_pos(3), ...
            this_voltage(1), this_voltage(2), this_voltage(3), 'r', ...
            'linewidth', 3, 'markersize', 100);
end

plot3(sourcemodel_and_leadfield.pos(inside_index, 1), ...
      sourcemodel_and_leadfield.pos(inside_index, 2), ...
      sourcemodel_and_leadfield.pos(inside_index, 3), 'o', ...
      'markersize', 40, 'markerfacecolor', 'r')
end