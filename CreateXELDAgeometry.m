%% This creates the geometry for ray tracing in XELDA
%
%
function [surface_list, photonspecs] = CreateXELDAgeometry()

grids_names = {'Cathode','Gate','Anode','Top'};
grids_origin = zeros(4,3);
grids_pitch = zeros(4,1) + .2; % distance between opposite sides of hexagon
grids_wirerad = [.0025; .0025; .006; .006];
grids_orientation = zeros(4,1);
grid_phase = 'llgg';
grids_origin(3:4,1) = grids_pitch(3:4) / sqrt(3);
grids_origin(:,3) = [-1.9185; -0.6435; 0; 0.6470];
z_zero_gridix = 1;
z_offset = -grids_origin(z_zero_gridix, 3);
grids_origin(:,3) = grids_origin(:,3) + z_offset;

n_xenon = 1.69; % E. Grace, triple point
n_gxenon = 1.0011262;
n_sio2 = 1.597;
n_ptfe = 1.79;

abslength_xenon = 460;
abslength_gxenon = 46000;
abslength_sio2 = 100;

scatlength_xenon = 29;% Coimbra, Seidel 2002  % 35; % E/ Grace, triple point
scatlength_gxenon = 1000;
scatlength_sio2 = 1000;


abs_ptfe = 0;%.01;
ptfe_mode = 'unified';
ptfe_uparams = [0, .978, 1, 0, 0];

abs_gasptfe = 0;%.01;
gasptfe_mode = 'unified';
gasptfe_uparams = [0, .7, 1, 0, 0];

liquidlevel = 0.5;

z_toppmt = 1.6430 + z_offset;
z_toppmtwin = 1.5230 + z_offset;
z_topptfe = 1.2880 + z_offset;
z_topvent = 0.9880 + z_offset;
z_tpctop = grids_origin(3,3) - grids_wirerad(3);
z_midvent = -0.4832 + z_offset;
z_tpcbot = grids_origin(1,3) + grids_wirerad(1);
z_notch = -2.0004 + z_offset;
z_botventtop = -3.1364 + z_offset;
z_botventbot = -3.5364 + z_offset;
z_botpmtwin = -3.7727 + z_offset;
z_botpmt = -3.8727 + z_offset;

r_gasgap = 3.3782;
r_tpc = 3.175;
r_notch = 3.4179;
r_botpmt = 3.2;

w_bigvent = .635;
w_midvent = .1;

phi_bigvent = 15*pi/180;
phi_midvent = 150*pi/180;

pmtwin_cornercylrad = 0.2540;
pmt_center = 1.35505;
pmt_width = 2.05;
pmtwin_tilt = pi/4;
pmtwin_botwidth = 2.5263;


%% derived quantities
grids_hexside = grids_pitch / sqrt(3);
z_liquidlevel = grids_origin(2,3) + liquidlevel*diff(grids_origin(2:3,3));
pmtwin_cornertilt = atan(sqrt(2)*tan(pmtwin_tilt));
pmtwin_cornereffrad = pmtwin_cornercylrad * cos(pmtwin_cornertilt) * sqrt(2/(1 + (cos(pmtwin_cornertilt))^2));
pmtwin_centeroffset = pmtwin_cornercylrad^2 / pmtwin_cornereffrad;
pmtwin_bevelextent = pmtwin_centeroffset - pmtwin_cornereffrad;
pmtwin_topwidth = pmtwin_botwidth - 2*tan(pmtwin_tilt)*(z_toppmtwin-z_topptfe);
pmtwin_botinner = pmtwin_botwidth - 2*pmtwin_cornereffrad;
pmtwin_topinner = pmtwin_topwidth - 2*pmtwin_cornereffrad;
pmt_halfwidth = .5*pmt_width;
pmtwin_tantilt = tan(pmtwin_tilt);

pmt_centers = zeros(4,2);
cyl_corner_centers = zeros(4,4,3);
cyl_corner_centers(:,:,3) = z_topptfe;
for i_pmt=1:4
    pmt_centers(i_pmt, 1) = pmt_center * (1 - 2*(i_pmt==2 | i_pmt==3));
    pmt_centers(i_pmt, 2) = pmt_center * (1 - 2*(i_pmt==1 | i_pmt==3));
    for i_corner=1:4
        cyl_corner_centers(i_pmt,i_corner, 1) = ...
            pmt_center * (1 - 2*(i_pmt==2 | i_pmt==3)) + ...
            (.5*pmtwin_botwidth - pmtwin_centeroffset) * (1 - 2*(i_corner==2 | i_corner==3));
        cyl_corner_centers(i_pmt,i_corner, 2) = ...
            pmt_center * (1 - 2*(i_pmt==1 | i_pmt==3)) + ...
            (.5*pmtwin_botwidth - pmtwin_centeroffset) * (1 - 2*(i_corner==1 | i_corner==3));
    end
end

cyl_axes = zeros(4,3);
cyl_axes(:,3) = cos(pmtwin_cornertilt);
for i_corner=1:4
    cyl_axes(i_corner,1) = ...
        sin(pmtwin_cornertilt) * sqrt(.5) * (1 - 2*(i_corner==1 | i_corner==4));
    cyl_axes(i_corner,2) = ...
        sin(pmtwin_cornertilt) * sqrt(.5) * (1 - 2*(i_corner==2 | i_corner==4));
end
cyl_levelaxes = cyl_axes ./ repmat(cyl_axes(:,3),1,3);
cyl_signaxes = sign(cyl_levelaxes);

bevel_centers = zeros(4,4,3);
bevel_centers(:,:,3) = z_topptfe;
for i_pmt=1:4
    for i_side=1:4
        bevel_centers(i_pmt,i_side, 1) = ...
            pmt_center * (1 - 2*(i_pmt==2 | i_pmt==3)) + ...
            .5*pmtwin_botwidth*((-1)^(1+i_side))*(i_side<3);
        bevel_centers(i_pmt,i_side, 2) = ...
            pmt_center * (1 - 2*(i_pmt==1 | i_pmt==3)) + ...
            .5*pmtwin_botwidth*((-1)^i_side)*(i_side>2);
    end
end

bevel_normals = zeros(4,3);
bevel_normals(:,3) = -cos(pmtwin_tilt);
bevel_normals(:,1) = [-1; 1; 0; 0] * sin(pmtwin_tilt);
bevel_normals(:,2) = [0; 0; 1; -1] * sin(pmtwin_tilt);

bevel_leveltan = zeros(4,3);
bevel_leveltan(1:2,2) = 1;
bevel_leveltan(3:4,1) = 1;

cphi_bigvent = cos(phi_bigvent);
sphi_bigvent = sin(phi_bigvent);
cphi_midvent = cos(phi_midvent);
sphi_midvent = sin(phi_midvent);

%% make empty surface struct array
surface_list = struct( ...
    'description', {}, ...
    'intersect_function', {}, ...
    'inbounds_function', {}, ...
    'n_outside', {}, ...
    'n_inside', {}, ...
    'surface_type', {}, ...
    'absorption', {}, ...
    'abslength_outside', {}, ...
    'abslength_inside', {}, ...
    'rayleigh_outside', {}, ...
    'rayleigh_inside', {}, ...
    'unifiedparams', {});

%%
% Make some hexagonal grids by magic (or trig)
surface_list(end+1).description = 'PMT: N';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_toppmt], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (abs(p(:,1,:) - pmt_center) <= pmt_halfwidth) & (abs(p(:,2,:) + pmt_center) <= pmt_halfwidth), ...
    size(p, 1), []));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'PMT: S';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_toppmt], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (abs(p(:,1,:) + pmt_center) <= pmt_halfwidth) & (abs(p(:,2,:) - pmt_center) <= pmt_halfwidth), ...
    size(p, 1), []));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'PMT: E';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_toppmt], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (abs(p(:,1,:) + pmt_center) <= pmt_halfwidth) & (abs(p(:,2,:) + pmt_center) <= pmt_halfwidth), ...
    size(p, 1), []));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'PMT: W';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_toppmt], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (abs(p(:,1,:) - pmt_center) <= pmt_halfwidth) & (abs(p(:,2,:) - pmt_center) <= pmt_halfwidth), ...
    size(p, 1), []));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'PMT: B';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_botpmt], [0, 0, -1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,1,:).^2 + p(:,2,:).^2) <= (r_botpmt^2), ...
    size(p, 1), []));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'TopPMTwindow: Missed';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_toppmt], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (abs(p(:,1,:)) <= (pmt_center - pmt_halfwidth)) | ...
    (abs(p(:,2,:)) <= (pmt_center - pmt_halfwidth)) | ...
    (abs(p(:,1,:)) >= (pmt_center + pmt_halfwidth)) | ...
    (abs(p(:,2,:)) >= (pmt_center + pmt_halfwidth)), ...
    size(p, 1), []));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'TopPMTwindow: yz shade';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_toppmt], [1, 0, 0]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    p(:, 3, :) >= z_toppmtwin, ...
    size(p, 1), []));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'TopPMTwindow: xz shade';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_toppmt], [0, 1, 0]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    p(:, 3, :) >= z_toppmtwin, ...
    size(p, 1), []));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'TopPMTwindow: xe-to-silica';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_toppmtwin], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(true(size(p,1),size(p,3)));
surface_list(end).n_outside = n_sio2;
surface_list(end).n_inside = n_gxenon;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_gxenon;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_gxenon;
surface_list(end).unifiedparams = zeros(1,5);

for i_pmt=1:4
    for i_side=1:4
        surface_list(end+1).description = sprintf('TopPMT flat bevel, PMT %d, side %d',i_pmt,i_side);
        surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
            reshape(bevel_centers(i_pmt, i_side, :),[1,3]), bevel_normals(i_side,:));
        surface_list(end).inbounds_function = @(p)(reshape( ...
            (p(:, 3, :) <= z_toppmtwin) & ...
            (p(:, 3, :) >= z_topptfe) & ...
            (abs((p(:,1,:)-pmt_centers(i_pmt,1)).*bevel_leveltan(i_side,1) + (p(:,2,:)-pmt_centers(i_pmt,2)).*bevel_leveltan(i_side,2)) <= ...
            (0.5 * pmtwin_botinner - (p(:, 3, :) - z_topptfe)*pmtwin_tantilt)), ...
            size(p, 1), []));
        surface_list(end).n_outside = n_gxenon;
        surface_list(end).n_inside = n_ptfe;
        surface_list(end).surface_type = gasptfe_mode;
        surface_list(end).absorption = abs_gasptfe;
        surface_list(end).abslength_outside = abslength_gxenon;
        surface_list(end).abslength_inside = abslength_gxenon;
        surface_list(end).rayleigh_outside = scatlength_gxenon;
        surface_list(end).rayleigh_inside = scatlength_gxenon;
        surface_list(end).unifiedparams = gasptfe_uparams;
    end
    for i_corner=1:4
        surface_list(end+1).description = sprintf('TopPMT cyl bevel, PMT %d, corner %d',i_pmt,i_side);
        surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
            reshape(cyl_corner_centers(i_pmt, i_corner, :),[1,3]), cyl_axes(i_corner,:), pmtwin_cornercylrad);
        surface_list(end).inbounds_function = @(p)(reshape( ...
            (p(:, 3, :) <= z_toppmtwin) & ...
            (p(:, 3, :) >= z_topptfe) & ...
            (((p(:,1,:) - (cyl_corner_centers(i_pmt,i_corner,1) + cyl_levelaxes(i_corner,1) .* (p(:,3,:) - z_topptfe))) * (-cyl_signaxes(i_corner, 1))) >= pmtwin_bevelextent) & ...
            (((p(:,2,:) - (cyl_corner_centers(i_pmt,i_corner,2) + cyl_levelaxes(i_corner,2) .* (p(:,3,:) - z_topptfe))) * (-cyl_signaxes(i_corner, 2))) >= pmtwin_bevelextent), ...
            size(p, 1), []));
        surface_list(end).n_outside = n_ptfe;
        surface_list(end).n_inside = n_gxenon;
        surface_list(end).surface_type = gasptfe_mode;
        surface_list(end).absorption = abs_gasptfe;
        surface_list(end).abslength_outside = abslength_gxenon;
        surface_list(end).abslength_inside = abslength_gxenon;
        surface_list(end).rayleigh_outside = scatlength_gxenon;
        surface_list(end).rayleigh_inside = scatlength_gxenon;
        surface_list(end).unifiedparams = gasptfe_uparams;
    end
end

surface_list(end+1).description = 'Top PTFE surface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_topptfe], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (abs(p(:,1,:)) <= (pmt_center - .5*pmtwin_botwidth)) | ...
    (abs(p(:,2,:)) <= (pmt_center - .5*pmtwin_botwidth)) | ...
    (abs(p(:,1,:)) >= (pmt_center + .5*pmtwin_botwidth)) | ...
    (abs(p(:,2,:)) >= (pmt_center + .5*pmtwin_botwidth)) | ...
    ( ...
    (((abs(p(:,1,:)) <= (pmt_center - .5*pmtwin_botinner)) & (abs(p(:,2,:)) <= (pmt_center - .5*pmtwin_botinner))) | ...
    ((abs(p(:,1,:)) >= (pmt_center + .5*pmtwin_botinner)) & (abs(p(:,2,:)) <= (pmt_center - .5*pmtwin_botinner))) | ...
    ((abs(p(:,1,:)) <= (pmt_center - .5*pmtwin_botinner)) & (abs(p(:,2,:)) >= (pmt_center + .5*pmtwin_botinner))) | ...
    ((abs(p(:,1,:)) >= (pmt_center + .5*pmtwin_botinner)) & (abs(p(:,2,:)) >= (pmt_center + .5*pmtwin_botinner)))) & ...
    (((p(:,1,:)-cyl_corner_centers(1,1,1)).^2 + (p(:,2,:)-cyl_corner_centers(1,1,2)).^2 - ((p(:,1,:)-cyl_corner_centers(1,1,1)).*cyl_axes(1,1)+(p(:,2,:)-cyl_corner_centers(1,1,2)).*cyl_axes(1,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(2,1,1)).^2 + (p(:,2,:)-cyl_corner_centers(2,1,2)).^2 - ((p(:,1,:)-cyl_corner_centers(2,1,1)).*cyl_axes(1,1)+(p(:,2,:)-cyl_corner_centers(2,1,2)).*cyl_axes(1,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(3,1,1)).^2 + (p(:,2,:)-cyl_corner_centers(3,1,2)).^2 - ((p(:,1,:)-cyl_corner_centers(3,1,1)).*cyl_axes(1,1)+(p(:,2,:)-cyl_corner_centers(3,1,2)).*cyl_axes(1,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(4,1,1)).^2 + (p(:,2,:)-cyl_corner_centers(4,1,2)).^2 - ((p(:,1,:)-cyl_corner_centers(4,1,1)).*cyl_axes(1,1)+(p(:,2,:)-cyl_corner_centers(4,1,2)).*cyl_axes(1,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(1,2,1)).^2 + (p(:,2,:)-cyl_corner_centers(1,2,2)).^2 - ((p(:,1,:)-cyl_corner_centers(1,2,1)).*cyl_axes(2,1)+(p(:,2,:)-cyl_corner_centers(1,2,2)).*cyl_axes(2,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(2,2,1)).^2 + (p(:,2,:)-cyl_corner_centers(2,2,2)).^2 - ((p(:,1,:)-cyl_corner_centers(2,2,1)).*cyl_axes(2,1)+(p(:,2,:)-cyl_corner_centers(2,2,2)).*cyl_axes(2,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(3,2,1)).^2 + (p(:,2,:)-cyl_corner_centers(3,2,2)).^2 - ((p(:,1,:)-cyl_corner_centers(3,2,1)).*cyl_axes(2,1)+(p(:,2,:)-cyl_corner_centers(3,2,2)).*cyl_axes(2,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(4,2,1)).^2 + (p(:,2,:)-cyl_corner_centers(4,2,2)).^2 - ((p(:,1,:)-cyl_corner_centers(4,2,1)).*cyl_axes(2,1)+(p(:,2,:)-cyl_corner_centers(4,2,2)).*cyl_axes(2,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(1,3,1)).^2 + (p(:,2,:)-cyl_corner_centers(1,3,2)).^2 - ((p(:,1,:)-cyl_corner_centers(1,3,1)).*cyl_axes(3,1)+(p(:,2,:)-cyl_corner_centers(1,3,2)).*cyl_axes(3,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(2,3,1)).^2 + (p(:,2,:)-cyl_corner_centers(2,3,2)).^2 - ((p(:,1,:)-cyl_corner_centers(2,3,1)).*cyl_axes(3,1)+(p(:,2,:)-cyl_corner_centers(2,3,2)).*cyl_axes(3,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(3,3,1)).^2 + (p(:,2,:)-cyl_corner_centers(3,3,2)).^2 - ((p(:,1,:)-cyl_corner_centers(3,3,1)).*cyl_axes(3,1)+(p(:,2,:)-cyl_corner_centers(3,3,2)).*cyl_axes(3,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(4,3,1)).^2 + (p(:,2,:)-cyl_corner_centers(4,3,2)).^2 - ((p(:,1,:)-cyl_corner_centers(4,3,1)).*cyl_axes(3,1)+(p(:,2,:)-cyl_corner_centers(4,3,2)).*cyl_axes(3,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(1,4,1)).^2 + (p(:,2,:)-cyl_corner_centers(1,4,2)).^2 - ((p(:,1,:)-cyl_corner_centers(1,4,1)).*cyl_axes(4,1)+(p(:,2,:)-cyl_corner_centers(1,4,2)).*cyl_axes(4,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(2,4,1)).^2 + (p(:,2,:)-cyl_corner_centers(2,4,2)).^2 - ((p(:,1,:)-cyl_corner_centers(2,4,1)).*cyl_axes(4,1)+(p(:,2,:)-cyl_corner_centers(2,4,2)).*cyl_axes(4,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(3,4,1)).^2 + (p(:,2,:)-cyl_corner_centers(3,4,2)).^2 - ((p(:,1,:)-cyl_corner_centers(3,4,1)).*cyl_axes(4,1)+(p(:,2,:)-cyl_corner_centers(3,4,2)).*cyl_axes(4,2)).^2) >= (pmtwin_cornercylrad^2)) & ...
    (((p(:,1,:)-cyl_corner_centers(4,4,1)).^2 + (p(:,2,:)-cyl_corner_centers(4,4,2)).^2 - ((p(:,1,:)-cyl_corner_centers(4,4,1)).*cyl_axes(4,1)+(p(:,2,:)-cyl_corner_centers(4,4,2)).*cyl_axes(4,2)).^2) >= (pmtwin_cornercylrad^2))), ...
    size(p, 1), []));
surface_list(end).n_outside = n_ptfe;
surface_list(end).n_inside = n_gxenon;
surface_list(end).surface_type = gasptfe_mode;
surface_list(end).absorption = abs_gasptfe;
surface_list(end).abslength_outside = abslength_gxenon;
surface_list(end).abslength_inside = abslength_gxenon;
surface_list(end).rayleigh_outside = abslength_gxenon;
surface_list(end).rayleigh_inside = scatlength_gxenon;
surface_list(end).unifiedparams = gasptfe_uparams;

for i_g=1:length(grids_names)
    surface_list(end+1).description = sprintf('Grid:  %s',grids_names{i_g});
    surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
        grids_origin(i_g,:), [0, 0, 1]);
    surface_list(end).inbounds_function = @(p)(reshape( ...
        ((mod((p(:,1,:)-grids_origin(i_g,1))*cos(grids_orientation(i_g)) + (p(:,2,:)-grids_origin(i_g,2))*sin(grids_orientation(i_g)), 3*grids_hexside(i_g)) <= grids_hexside(i_g)) & ...
        (abs(mod((p(:,1,:)-grids_origin(i_g,1))*sin(grids_orientation(i_g)) - (p(:,2,:)-grids_origin(i_g,2))*cos(grids_orientation(i_g)) + .5*grids_pitch(i_g), grids_pitch(i_g)) - .5*grids_pitch(i_g)) <= grids_wirerad(i_g))) | ...
        ((mod((p(:,1,:)-grids_origin(i_g,1))*cos(grids_orientation(i_g) + (2*pi/3)) + (p(:,2,:)-grids_origin(i_g,2))*sin(grids_orientation(i_g) + (2*pi/3)), 3*grids_hexside(i_g)) <= grids_hexside(i_g)) & ...
        (abs(mod((p(:,1,:)-grids_origin(i_g,1))*sin(grids_orientation(i_g) + (2*pi/3)) - (p(:,2,:)-grids_origin(i_g,2))*cos(grids_orientation(i_g) + (2*pi/3)) + .5*grids_pitch(i_g), grids_pitch(i_g)) - .5*grids_pitch(i_g)) <= grids_wirerad(i_g))) | ...
        ((mod((p(:,1,:)-grids_origin(i_g,1))*cos(grids_orientation(i_g) - (2*pi/3)) + (p(:,2,:)-grids_origin(i_g,2))*sin(grids_orientation(i_g) - (2*pi/3)), 3*grids_hexside(i_g)) <= grids_hexside(i_g)) & ...
        (abs(mod((p(:,1,:)-grids_origin(i_g,1))*sin(grids_orientation(i_g) - (2*pi/3)) - (p(:,2,:)-grids_origin(i_g,2))*cos(grids_orientation(i_g) - (2*pi/3)) + .5*grids_pitch(i_g), grids_pitch(i_g)) - .5*grids_pitch(i_g)) <= grids_wirerad(i_g))) | ...
        ((mod((p(:,1,:)-grids_origin(i_g,1))*cos(grids_orientation(i_g)) + (p(:,2,:)-grids_origin(i_g,2))*sin(grids_orientation(i_g)) + 1.5*grids_hexside(i_g), 3*grids_hexside(i_g)) <= grids_hexside(i_g)) & ...
        (abs(mod((p(:,1,:)-grids_origin(i_g,1))*sin(grids_orientation(i_g)) - (p(:,2,:)-grids_origin(i_g,2))*cos(grids_orientation(i_g)), grids_pitch(i_g)) - .5*grids_pitch(i_g)) <= grids_wirerad(i_g))) | ...
        ((mod((p(:,1,:)-grids_origin(i_g,1))*cos(grids_orientation(i_g) + (2*pi/3)) + (p(:,2,:)-grids_origin(i_g,2))*sin(grids_orientation(i_g) + (2*pi/3)) + 1.5*grids_hexside(i_g), 3*grids_hexside(i_g)) <= grids_hexside(i_g)) & ...
        (abs(mod((p(:,1,:)-grids_origin(i_g,1))*sin(grids_orientation(i_g) + (2*pi/3)) - (p(:,2,:)-grids_origin(i_g,2))*cos(grids_orientation(i_g) + (2*pi/3)), grids_pitch(i_g)) - .5*grids_pitch(i_g)) <= grids_wirerad(i_g))) | ...
        ((mod((p(:,1,:)-grids_origin(i_g,1))*cos(grids_orientation(i_g) - (2*pi/3)) + (p(:,2,:)-grids_origin(i_g,2))*sin(grids_orientation(i_g) - (2*pi/3)) + 1.5*grids_hexside(i_g), 3*grids_hexside(i_g)) <= grids_hexside(i_g)) & ...
        (abs(mod((p(:,1,:)-grids_origin(i_g,1))*sin(grids_orientation(i_g) - (2*pi/3)) - (p(:,2,:)-grids_origin(i_g,2))*cos(grids_orientation(i_g) - (2*pi/3)), grids_pitch(i_g)) - .5*grids_pitch(i_g)) <= grids_wirerad(i_g))) ...
        , size(p,1), []));
    surface_list(end).n_outside = (n_xenon*(grid_phase(i_g)=='l')) + (n_gxenon*(grid_phase(i_g)=='g'));
    surface_list(end).n_inside = (n_xenon*(grid_phase(i_g)=='l')) + (n_gxenon*(grid_phase(i_g)=='g'));
    surface_list(end).surface_type = 'normal';
    surface_list(end).absorption = 1; % call grids black
    surface_list(end).abslength_outside = (abslength_xenon*(grid_phase(i_g)=='l')) + (abslength_gxenon*(grid_phase(i_g)=='g')); % photons shouldn't get there
    surface_list(end).abslength_inside = (abslength_xenon*(grid_phase(i_g)=='l')) + (abslength_gxenon*(grid_phase(i_g)=='g')); % xenon
    surface_list(end).rayleigh_outside = (scatlength_xenon*(grid_phase(i_g)=='l')) + (scatlength_gxenon*(grid_phase(i_g)=='g')); % photons shouldn't get there
    surface_list(end).rayleigh_inside = (scatlength_xenon*(grid_phase(i_g)=='l')) + (scatlength_gxenon*(grid_phase(i_g)=='g'));
    surface_list(end).unifiedparams = zeros(1,5);
end

surface_list(end+1).description = 'gas region wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], r_gasgap);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:, 3, :) <= z_topptfe) & ...
    (p(:, 3, :) >= z_tpctop) & ...
    ((p(:,3,:) <= z_topvent) | ...
    ((abs(p(:, 1, :)*cphi_bigvent + p(:, 2, :)*sphi_bigvent)>=(.5*w_bigvent)) & ...
    (abs(p(:, 1, :)*sphi_bigvent - p(:, 2, :)*cphi_bigvent)>=(.5*w_bigvent)))), ...
    size(p, 1), []));
surface_list(end).n_outside = n_ptfe;
surface_list(end).n_inside = n_gxenon;
surface_list(end).surface_type = gasptfe_mode;
surface_list(end).absorption = abs_gasptfe;
surface_list(end).abslength_outside = abslength_gxenon;
surface_list(end).abslength_inside = abslength_gxenon;
surface_list(end).rayleigh_outside = scatlength_gxenon;
surface_list(end).rayleigh_inside = scatlength_gxenon;
surface_list(end).unifiedparams = gasptfe_uparams;

surface_list(end+1).description = 'gas region wall vents';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], r_gasgap);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) >= z_topvent) & ...
    (p(:,3,:) <= z_topptfe) & ...
    ((abs(p(:, 1, :)*cphi_bigvent + p(:, 2, :)*sphi_bigvent)<=(.5*w_bigvent)) | ...
    (abs(p(:, 1, :)*sphi_bigvent - p(:, 2, :)*cphi_bigvent)<=(.5*w_bigvent))), ...
    size(p, 1), []));
surface_list(end).n_outside = n_gxenon;
surface_list(end).n_inside = n_gxenon;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_gxenon;
surface_list(end).abslength_inside = abslength_gxenon;
surface_list(end).rayleigh_outside = scatlength_gxenon;
surface_list(end).rayleigh_inside = scatlength_gxenon;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'tpc region wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], r_tpc);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) >= z_liquidlevel) & ...
    (((p(:, 3, :) <= z_notch) & ...
    (p(:, 3, :) >= z_botpmtwin) & ...
    ((p(:,3,:) <= z_botventbot) | (p(:,3,:) >= z_botventtop) | ...
    ((abs(p(:, 1, :)*cphi_bigvent + p(:, 2, :)*sphi_bigvent)>=(.5*w_bigvent)) & ...
    (abs(p(:, 1, :)*sphi_bigvent - p(:, 2, :)*cphi_bigvent)>=(.5*w_bigvent))))) | ...
    ((p(:,3,:) <= z_tpctop) & ...
    (p(:,3,:) >= z_tpcbot) & ...
    ((p(:,3,:) <= z_midvent) | ...
    (abs(p(:,1,:)*sphi_midvent - p(:,2,:)*cphi_midvent) >= (.5*w_midvent)) | ...
    ((p(:,1,:)*cphi_midvent + p(:,2,:)*sphi_midvent) <= 0)))), ...
    size(p, 1), []));
surface_list(end).n_outside = n_ptfe;
surface_list(end).n_inside = n_gxenon;
surface_list(end).surface_type = gasptfe_mode;
surface_list(end).absorption = abs_gasptfe;
surface_list(end).abslength_outside = abslength_gxenon;
surface_list(end).abslength_inside = abslength_gxenon;
surface_list(end).rayleigh_outside = scatlength_gxenon;
surface_list(end).rayleigh_inside = scatlength_gxenon;
surface_list(end).unifiedparams = gasptfe_uparams;

surface_list(end+1).description = 'tpc region wall vents';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], r_tpc);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) >= z_liquidlevel) & ...
    (((p(:,3,:) >= z_botventbot) & ...
    (p(:,3,:) <= z_botventtop) & ...
    ((abs(p(:, 1, :)*cphi_bigvent + p(:, 2, :)*sphi_bigvent)<=(.5*w_bigvent)) | ...
    (abs(p(:, 1, :)*sphi_bigvent - p(:, 2, :)*cphi_bigvent)<=(.5*w_bigvent)))) | ...
    ((p(:,3,:) >= z_midvent) & ...
    (p(:,3,:) <= z_tpctop) & ...
    (abs(p(:,1,:)*sphi_midvent - p(:,2,:)*cphi_midvent) <= (.5*w_midvent)) & ...
    ((p(:,1,:)*cphi_midvent + p(:,2,:)*sphi_midvent) >= 0))), ...
    size(p, 1), []));
surface_list(end).n_outside = n_gxenon;
surface_list(end).n_inside = n_gxenon;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_gxenon;
surface_list(end).abslength_inside = abslength_gxenon;
surface_list(end).rayleigh_outside = scatlength_gxenon;
surface_list(end).rayleigh_inside = scatlength_gxenon;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'tpc region wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], r_tpc);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) <= z_liquidlevel) & ...
    (((p(:, 3, :) <= z_notch) & ...
    (p(:, 3, :) >= z_botpmtwin) & ...
    ((p(:,3,:) <= z_botventbot) | (p(:,3,:) >= z_botventtop) | ...
    ((abs(p(:, 1, :)*cphi_bigvent + p(:, 2, :)*sphi_bigvent)>=(.5*w_bigvent)) & ...
    (abs(p(:, 1, :)*sphi_bigvent - p(:, 2, :)*cphi_bigvent)>=(.5*w_bigvent))))) | ...
    ((p(:,3,:) <= z_tpctop) & ...
    (p(:,3,:) >= z_tpcbot) & ...
    ((p(:,3,:) <= z_midvent) | ...
    (abs(p(:,1,:)*sphi_midvent - p(:,2,:)*cphi_midvent) >= (.5*w_midvent)) | ...
    ((p(:,1,:)*cphi_midvent + p(:,2,:)*sphi_midvent) <= 0)))), ...
    size(p, 1), []));
surface_list(end).n_outside = n_ptfe;
surface_list(end).n_inside = n_xenon;
surface_list(end).surface_type = ptfe_mode;
surface_list(end).absorption = abs_ptfe;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = scatlength_xenon;
surface_list(end).rayleigh_inside = scatlength_xenon;
surface_list(end).unifiedparams = ptfe_uparams;

surface_list(end+1).description = 'tpc region wall vents';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], r_tpc);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) <= z_liquidlevel) & ...
    (((p(:,3,:) >= z_botventbot) & ...
    (p(:,3,:) <= z_botventtop) & ...
    ((abs(p(:, 1, :)*cphi_bigvent + p(:, 2, :)*sphi_bigvent)<=(.5*w_bigvent)) | ...
    (abs(p(:, 1, :)*sphi_bigvent - p(:, 2, :)*cphi_bigvent)<=(.5*w_bigvent)))) | ...
    ((p(:,3,:) >= z_midvent) & ...
    (p(:,3,:) <= z_tpctop) & ...
    (abs(p(:,1,:)*sphi_midvent - p(:,2,:)*cphi_midvent) <= (.5*w_midvent)) & ...
    ((p(:,1,:)*cphi_midvent + p(:,2,:)*sphi_midvent) >= 0))), ...
    size(p, 1), []));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = n_xenon;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = scatlength_xenon;
surface_list(end).rayleigh_inside = scatlength_xenon;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'tpc region notch wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], r_notch);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:, 3, :) >= z_notch) & ...
    (p(:, 3, :) <= z_tpcbot), ...
    size(p, 1), []));
surface_list(end).n_outside = n_ptfe;
surface_list(end).n_inside = n_xenon;
surface_list(end).surface_type = ptfe_mode;
surface_list(end).absorption = abs_ptfe;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = scatlength_xenon;
surface_list(end).rayleigh_inside = scatlength_xenon;
surface_list(end).unifiedparams = ptfe_uparams;

surface_list(end+1).description = 'Bot PMT window side wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0, 0, 0], [0, 0, 1], r_botpmt);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:,3,:) >= z_botpmt) & ...
    (p(:,3,:) <= z_botpmtwin), ...
    size(p, 1), []));
surface_list(end).n_outside = n_gxenon;
surface_list(end).n_inside = n_gxenon;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
surface_list(end).abslength_outside = abslength_sio2;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_sio2;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'TPC shelf';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_tpctop], [0, 0,1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:, 1, :).^2 + p(:, 2, :).^2) >= (r_tpc^2), ...
    size(p, 1), []));
surface_list(end).n_outside = n_gxenon;
surface_list(end).n_inside = n_ptfe;
surface_list(end).surface_type = gasptfe_mode;
surface_list(end).absorption = abs_gasptfe;
surface_list(end).abslength_outside = abslength_gxenon;
surface_list(end).abslength_inside = abslength_gxenon;
surface_list(end).rayleigh_outside = scatlength_gxenon;
surface_list(end).rayleigh_inside = scatlength_gxenon;
surface_list(end).unifiedparams = gasptfe_uparams;

surface_list(end+1).description = 'Liquid surface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_liquidlevel], [0, 0,1]);
surface_list(end).inbounds_function = @(p)(true(size(p,1),size(p,3)));
surface_list(end).n_outside = n_gxenon;
surface_list(end).n_inside = n_xenon;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = abslength_gxenon;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = scatlength_gxenon;
surface_list(end).rayleigh_inside = scatlength_xenon;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'TPC notch top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_tpcbot], [0, 0, -1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:, 1, :).^2 + p(:, 2, :).^2) >= (r_tpc^2), ...
    size(p, 1), []));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = n_ptfe;
surface_list(end).surface_type = ptfe_mode;
surface_list(end).absorption = abs_ptfe;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = scatlength_xenon;
surface_list(end).rayleigh_inside = scatlength_xenon;
surface_list(end).unifiedparams = ptfe_uparams;

surface_list(end+1).description = 'TPC notch bot';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_notch], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:, 1, :).^2 + p(:, 2, :).^2) >= (r_tpc^2), ...
    size(p, 1), []));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = n_ptfe;
surface_list(end).surface_type = ptfe_mode;
surface_list(end).absorption = abs_ptfe;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = scatlength_xenon;
surface_list(end).rayleigh_inside = scatlength_xenon;
surface_list(end).unifiedparams = ptfe_uparams;

surface_list(end+1).description = 'bot PMT silica to xenon';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_botpmtwin], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(true(size(p,1),size(p,3)));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = n_sio2;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_sio2;
surface_list(end).rayleigh_outside = scatlength_xenon;
surface_list(end).rayleigh_inside = scatlength_sio2;
surface_list(end).unifiedparams = zeros(1,5);

surface_list(end+1).description = 'TPC gass wall top';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0, 0, z_topptfe], [0, 0, 1]);
surface_list(end).inbounds_function = @(p)(reshape( ...
    (p(:, 1, :).^2 + p(:, 2, :).^2) >= (r_gasgap^2), ...
    size(p, 1), []));
surface_list(end).n_outside = n_xenon;
surface_list(end).n_inside = n_ptfe;
surface_list(end).surface_type = ptfe_mode;
surface_list(end).absorption = abs_ptfe;
surface_list(end).abslength_outside = abslength_xenon;
surface_list(end).abslength_inside = abslength_xenon;
surface_list(end).rayleigh_outside = scatlength_xenon;
surface_list(end).rayleigh_inside = scatlength_xenon;
surface_list(end).unifiedparams = ptfe_uparams;

%% Use the photonspecs structure to determine where to put photons
photonspecs = struct();
photonspecs.tpc_bot = grids_origin(1,3);
photonspecs.tpc_top = grids_origin(2,3);
photonspecs.liquidlevel = z_liquidlevel;
photonspecs.s2_top = grids_origin(3,3);
photonspecs.tpc_rad = r_tpc;

l1 = grids_pitch(2) * [0; 1];
l2 = grids_pitch(2) * [-sqrt(3)*.5; .5];
grid_rot = [cos(grids_orientation(2)), -sin(grids_orientation(2)) ; sin(grids_orientation(2)), cos(grids_orientation(2))];
l1 = grid_rot * l1;
l2 = grid_rot * l2;

gridsize = ceil(10*r_tpc/grids_pitch(2));
steps1 = repmat(-gridsize:gridsize, 2*gridsize+1, 1);
steps2 = repmat((-gridsize:gridsize)', 1, 2*gridsize+1);

xyarray = repmat(steps1(:), 1, 2) .* repmat(l1', numel(steps1), 1) + ...
    repmat(steps2(:), 1, 2) .* repmat(l2', numel(steps2), 1) + ...
    repmat(grids_origin(2,1:2) - grids_hexside(2)*[cos(grids_orientation(2)),sin(grids_orientation(2))], numel(steps1), 1);
xyarraycut = sum(xyarray.^2, 2) < r_tpc^2;

photonspecs.s2_xyarray = xyarray(xyarraycut, :);



