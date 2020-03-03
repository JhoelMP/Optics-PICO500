  %% CommentedOpticsCode.m

% This code crend acts as the hub to call all the other required functions. It
% is adaptedates all the graphs used in simulating the optics of
% PICO-40L a from the code called pressure_vessel_optics_simulator.m in
% Eric's user code in the file couppjunk.

%%

% Starting with the camera specs
calidad=1/8
% The resolution is in pixels, [horizontal vertical]
geospecs.cam_resolution = round([1200 1920]*calidad-1) %[1088 2048]; 
% cam_pixelpitch is the pixel pitch of the sensor, which is the distance
% between the center of two pixels
geospecs.cam_pixelpitch = 0.00059/calidad
% cam_f is the focal length of the lens

% cam_barreld this is the radial distortion
geospecs.cam_barreld = 0;
% cam_lenstype, not exactly sure what this represents, see
% GenerateRaysFromCameraCommented for how it effects the code
geospecs.cam_lenstype = 'theta';
% cam_sensorsize is the physical size of the sensor
geospecs.cam_sensorsize = geospecs.cam_resolution .* geospecs.cam_pixelpitch;

%% Create ray lists coming out of the camera
% This function creates a map of the pixels in the camera sensor and the
% direction of the ray from each pixel as it leaves the camera
[raydirectionsupper, pixelmap] = GenerateRaysFromCameraCommented(geospecs.cam_resolution, geospecs.cam_pixelpitch, .5*(1+geospecs.cam_resolution), geospecs.cam_f, ...
    geospecs.camupper_pitch, geospecs.camupper_yaw, geospecs.camupper_roll, geospecs.cam_barreld, geospecs.cam_lenstype);
    
% rays amends the ray directions from above to have 0,0,1,1,0,0,0 after
% each entry, note that each line of rays{1} corresponds to the pixel held
% in the same line of pixelmap/pixels{1}
rays{1} = [raydirectionsupper repmat([0 0 1 1 0 0 0],size(raydirectionsupper,1),1)];
pixels{1} = pixelmap;
% ray_startingpoints gives each of the rays their geometric starting
% position according to the geospecs outlined above
ray_startingpoints{1} = repmat([geospecs.camupper_x (geospecs.camupper_y - cos(geospecs.vp_upperangle)) (geospecs.camupper_z - sin(geospecs.vp_upperangle) + geospecs.camupper_zoffset)],size(raydirectionsupper,1),1);

% We repeat everything we just did for the lower camera
[raydirectionslower, pixelmap] = GenerateRaysFromCameraCommented(geospecs.cam_resolution, geospecs.cam_pixelpitch, .5*(1+geospecs.cam_resolution), geospecs.cam_f, ...
    geospecs.camlower_pitch, geospecs.camlower_yaw, geospecs.camlower_roll, geospecs.cam_barreld, geospecs.cam_lenstype);
    
rays{2} = [raydirectionslower repmat([0 0 1 1 0 0 0],size(raydirectionslower,1),1)];
pixels{2} = pixelmap;
ray_startingpoints{2} = repmat([geospecs.camlower_x (geospecs.camlower_y - cos(geospecs.vp_lowerangle)) (geospecs.camlower_z - sin(geospecs.vp_lowerangle) + geospecs.camlower_zoffset)],size(raydirectionslower,1),1);

% Again we repeat but this time it is for the upper camera that is not on
% the x-axis
[raydirections2upper, pixelmap2] = GenerateRaysFromCameraCommented(geospecs.cam_resolution, geospecs.cam_pixelpitch, .5*(1+geospecs.cam_resolution), geospecs.cam_f, ...
    geospecs.camupper_pitch, (geospecs.camupper_yaw - geospecs.vp_phi), geospecs.camupper_roll, geospecs.cam_barreld, geospecs.cam_lenstype);
    
rays{3} = [raydirections2upper repmat([0 0 1 1 0 0 0],size(raydirections2upper,1),1)];
pixels{3} = pixelmap2;
ray_startingpoints{3} = repmat([((geospecs.camupper_y - cos(geospecs.vp_upperangle)) * sin(-geospecs.vp_phi)) ((geospecs.camupper_y - cos(geospecs.vp_upperangle)) * cos(geospecs.vp_phi)) (geospecs.camupper_z - sin(geospecs.vp_upperangle))],size(raydirections2upper,1),1);

% Again we repeat but this time it is for the lower camera that is not on
% the x-axis
[raydirections2lower, pixelmap2] = GenerateRaysFromCameraCommented(geospecs.cam_resolution, geospecs.cam_pixelpitch, .5*(1+geospecs.cam_resolution), geospecs.cam_f, ...
    geospecs.camlower_pitch, (geospecs.camlower_yaw - geospecs.vp_phi), geospecs.camlower_roll, geospecs.cam_barreld, geospecs.cam_lenstype);

rays{4} = [raydirections2lower repmat([0 0 1 1 0 0 0],size(raydirections2lower,1),1)];
pixels{4} = pixelmap2;
ray_startingpoints{4} = repmat([((geospecs.camlower_y - cos(geospecs.vp_lowerangle)) * sin(-geospecs.vp_phi)) ((geospecs.camlower_y - cos(geospecs.vp_lowerangle)) * cos(geospecs.vp_phi)) (geospecs.camlower_z - sin(geospecs.vp_lowerangle))],size(raydirections2lower,1),1);

fprintf('Rays and pixels created\n')
%%  Y-Z cross section of chamber
% We first create the number of rays that will be sent out from the camera
numrays = 256; %default was 64
% Depending on the lens type there is a different max angle at which rays 
% can be released which is calculated below
switch geospecs.cam_lenstype
    case 'tan'
        max_angle = atan(.5*geospecs.cam_sensorsize(2)*(1+sum(geospecs.cam_barreld.*(.5*geospecs.cam_sensorsize(2)/geospecs.cam_f).^(2*(1:length(geospecs.cam_barreld)))))/geospecs.cam_f);
    case 'theta'
        max_angle = .5*geospecs.cam_sensorsize(2)/geospecs.cam_f;
end

% This variable gives the starting coordinates for the camera center to
% each of the rays being created for the upper camera
yz_startingpointsupper = repmat(ray_startingpoints{1}(1,:), numrays, 1);
yz_startingpointslower = repmat(ray_startingpoints{2}(1,:), numrays, 1);
yz_startingpoints = [yz_startingpointsupper; yz_startingpointslower];

% td creates an array of angles from the - max angle to the + max angle and
% adjusts by the camera's pitch
tdupper = linspace(-max_angle, max_angle, numrays)+geospecs.camupper_pitch;
tdlower = linspace(-max_angle, max_angle, numrays)+geospecs.camlower_pitch;
% yz_rays gives a 10 entry row to each ray to be used in RayTracerRetro
yz_raysupper = repmat([0 0 0 1 0 0 1 0 0 0],numrays,1);
yz_rayslower = repmat([0 0 0 1 0 0 1 0 0 0],numrays,1);
% This replaces the 2nd column of each ray with the cos value for the
% corresponding angle in td
yz_raysupper(:,2) = cos(tdupper);
yz_rayslower(:,2) = cos(tdlower);
% This replaces the 3rd column of each ray with the sin value of the
% corresponding angle in td
yz_raysupper(:,3) = sin(tdupper);
yz_rayslower(:,3) = sin(tdlower);
% Note that now the 2nd and 3rd columns give the direction of the ray in the yz plane

yz_rays = [yz_raysupper; yz_rayslower];

% max_scatters gives the number of scatters (interactions with a surface)
% a ray can go through before the trace stops
max_scatters = 16;
% yz_output is the end result of tracing the rays, for more information see
% RayTracerRetro
yz_output = RayTracerRetro(yz_startingpoints, yz_rays, surface_list, max_scatters, 1e-5, [0 100]);

% Scatter points as a variable gives the location of the ray at each point
% it interacts with a surface
scatter_points = zeros([max_scatters+2, 2 * numrays, 3]);
scatter_points(:) = NaN;
scatter_points(1,:,:) = yz_startingpoints;

for n=1:length(yz_output)
    scatter_points(n+1,abs(yz_output(n).ray_index),:) = yz_output(n).intersection_point;
end

% The following plots the rays in the yz-plane
plotrays = reshape(scatter_points,[],3);
close all
figure(101);
clf;
set(gca,'fontsize',16);
hold on
axis equal
yz_lines = SectionPlotter(surface_list,yz_startingpoints(1,:),[1 0 0],10000,101);
plot(plotrays(1:(numrays*(max_scatters+2)),2),plotrays(1:(numrays*(max_scatters+2)),3),'-or');
plot(plotrays(1+(numrays*(max_scatters+2)):size(plotrays,1),2),plotrays(1+(numrays*(max_scatters+2)):size(plotrays,1),3),'-og');
title('YZ Cross-section')
xlabel('Y Distance (cm)')
ylabel('Z Distance (cm)')
axis([-75 75 -1 210]) 
fig1 = gcf;
savefig(fig1,strcat('YZ_',char(datetime('now','Format','MM_dd_hh_mm'))));

fprintf('YZ cross-section figure displayed\n')
%%  X-Y cross section of chamber
% This acts essentially the same as the yz cross section code above but in
% the xy plane instead
numrays = 256;
switch geospecs.cam_lenstype
    case 'tan'
        max_angle = atan(.5*geospecs.cam_sensorsize(1)*(1+sum(geospecs.cam_barreld.*(.5*geospecs.cam_sensorsize(1)/geospecs.cam_f).^(2*(1:length(geospecs.cam_barreld)))))/geospecs.cam_f);
    case 'theta'
        max_angle = .5*geospecs.cam_sensorsize(1)/geospecs.cam_f;
end

xy_startingpointsleft = repmat(ray_startingpoints{1}(1,:), numrays, 1);
xy_startingpointsright = repmat(ray_startingpoints{3}(1,:), numrays, 1);
xy_startingpoints = [xy_startingpointsleft; xy_startingpointsright];
%xy_startingpoints_ = [xy_startingpointsleft; xy_startingpointsright]+[zeros(1,128); zeros(1,128); 15*ones(1,128)]';

tdleft = linspace(-max_angle, max_angle, numrays) - geospecs.camupper_yaw;
tdright = linspace(-max_angle, max_angle, numrays) - geospecs.camupper_yaw - geospecs.vp_phi;

xy_raysleft = repmat([0 0 0 1 0 0 1 0 0 0],numrays,1);
xy_raysright = repmat([0 0 0 1 0 0 1 0 0 0],numrays,1);
xy_raysleft(:,2) = cos(tdleft);
xy_raysright(:,2) = cos(tdright);
xy_raysleft(:,1) = sin(tdleft);
xy_raysright(:,1) = sin(tdright);
xy_rays = [xy_raysleft;xy_raysright];

max_scatters = 16;
xy_output = RayTracerRetro(xy_startingpoints, xy_rays, surface_list, max_scatters, 1e-5, [0 100]);

scatter_points = zeros([max_scatters+2, 2*numrays, 3]);
scatter_points(:) = NaN;
scatter_points(1,:,:) = xy_startingpoints;

for n=1:length(xy_output)
    scatter_points(n+1,abs(xy_output(n).ray_index),:) = xy_output(n).intersection_point;
end

plotrays = reshape(scatter_points,[],3);
figure(102);
clf;
set(gca,'fontsize',16);
hold on
axis equal
SectionPlotter(surface_list,xy_startingpoints(1,:),[0 0 1],10000,102);
plot(plotrays(1:(numrays*(max_scatters+2)),1),plotrays(1:(numrays*(max_scatters+2)),2),'-or');
plot(plotrays(1+(numrays*(max_scatters+2)):size(plotrays,1),1),plotrays(1+(numrays*(max_scatters+2)):size(plotrays,1),2),'-og');
title('XY Cross-section')
xlabel('X Distance (cm)')
ylabel('Y Distance (cm)')
axis([-50 60 -60 50])
fig2 = gcf;
savefig(fig2,strcat('XY_',char(datetime('now','Format','MM_dd_hh_mm'))));

fprintf('XY cross-section figure displayed\n')

%% run RayTracerRetro test ray
% Test pixel uses the pixel at the centre of the camera sensor
test_pixel1 = fix(.5*max(pixels{1}, [], 1));
% This gives the index corresponding to the test_pixel
test_ix1 = find(all(pixels{1}==repmat(test_pixel1,size(pixels{1},1),1),2));

% The following runs the test ray
raytracer_test1 = RayTracerRetro(ray_startingpoints{1}(test_ix1,:), ...
    rays{1}(test_ix1,:), surface_list, 100, 1e-5, [1e-5 100]);

% This sets out_of_cf3i_index to be the scatter number when the ray leaves
% the target fluid
out_of_cf3i_index1 = 0;

out_of_target_surfaces1 = [-1, -3, -5, 8, 10, 12];

for sn=1:length(raytracer_test1)
    if any(raytracer_test1(sn).surface_index(1)==out_of_target_surfaces1)
        out_of_cf3i_index1 = sn;
        break
    end
end

% This sets first_retro_surface to be the scatter number when the ray hits
% the retroreflector
first_retro_surface1 = 0;

for sn=1:length(surface_list)
    if strcmp(surface_list(sn).surface_type,'retro')
        first_retro_surface1 = sn;
        break
    end
end

% We now run test pixel 2 for the lower camera
test_pixel2 = fix(.5*max(pixels{2}, [], 1));
test_ix2 = find(all(pixels{2}==repmat(test_pixel2,size(pixels{2},1),1),2));

raytracer_test2 = RayTracerRetro(ray_startingpoints{2}(test_ix2,:), ...
    rays{2}(test_ix2,:), surface_list, 100, 1e-5, [1e-5 100]);

out_of_cf3i_index2 = 0;

out_of_target_surfaces2 = [-1, -3, -5, 8, 10, 12];

for sn=1:length(raytracer_test2)
    if any(raytracer_test2(sn).surface_index(1)==out_of_target_surfaces2)
        out_of_cf3i_index2 = sn;
        break
    end
end

first_retro_surface2 = 0;

for sn=1:length(surface_list)
    if strcmp(surface_list(sn).surface_type,'retro')
        first_retro_surface2 = sn;
        break
    end
end

fprintf('Test ray run\n')
%% run RayTracerRetro
raytracer_output1 = RayTracerRetro(ray_startingpoints{1}, ...
    rays{1}, surface_list, length(raytracer_test1) + 5, 1e-5, [1e-5 100]);

% Not 100% sure why this warning is here, it seems to be inevitable by
% input of RayTracerRetro
if (length(raytracer_output1)-length(raytracer_test1)) > 1
    disp('Might need bigger scatter history1...');
end

% This is for the lower camera
raytracer_output2 = RayTracerRetro(ray_startingpoints{2}, ...
    rays{2}, surface_list, length(raytracer_test2) + 5, 1e-5, [1e-5 100]);
if (length(raytracer_output2)-length(raytracer_test2)) > 1
    disp('Might need bigger scatter history2...');
end


%% Generic discrete image w/o reflections
surface_history1 = zeros(size(rays{1},1), max_scatters);
for scatternum = 1:length(raytracer_output1)
    % This cuts all rays that have reflected at some point
    noreflectionscut = raytracer_output1(scatternum).ray_index > 1;
    % This fills surface history with just the surfaces from scatters
    % before reflection is taken into account
    surface_history1(raytracer_output1(scatternum).ray_index(noreflectionscut), scatternum) = ...
        raytracer_output1(scatternum).surface_index(noreflectionscut);
end

% Unsure of what this does
[scatter_history_legend1, ~, scatter_history_ix1] = unique(surface_history1, 'rows');

% Takes the maximum of the pixels column giving the resolution of the
% sensor
resolution1 = max(pixels{1}, [], 1);

% Random numbers between 0 and 1 for the column size of
% scatter_history_legend in rows by 3 columns
cmap = rand(size(scatter_history_legend1,1),3);

% Then take the random number from the array to put into the image as a
% colour for a certain surface defined by scatter_history_ix
scatter_history_image1 = cmap(scatter_history_ix1,:);

% Reshape it into RGB with the proper number of pixels
scatter_history_image1 = reshape(scatter_history_image1,[resolution1 3]);

% We repeat for cam2
surface_history2 = zeros(size(rays{2},1), max_scatters);
for scatternum = 1:length(raytracer_output2)
    noreflectionscut = raytracer_output2(scatternum).ray_index > 0;
    surface_history2(raytracer_output2(scatternum).ray_index(noreflectionscut), scatternum) = ...
        raytracer_output2(scatternum).surface_index(noreflectionscut);
end
[scatter_history_legend2, ~, scatter_history_ix2] = unique(surface_history2, 'rows');
resolution2 = max(pixels{2}, [], 1);
scatter_history_image2 = cmap(scatter_history_ix2,:);
scatter_history_image2 = reshape(scatter_history_image2,[resolution2 3]);

% Creates the figure of the detector without reflections
figure(111);
clf;

subplot(1,2,1)
set(gca,'fontsize',16);
image(permute(scatter_history_image1,[2 1 3]));
axis image
title('Upper Camera')
xlabel('X Pixels')
ylabel('Y Pixels')

subplot(1,2,2)
set(gca,'fontsize',16);
image(permute(scatter_history_image2,[2 1 3]));
axis image
title('Lower Camera')
xlabel('X Pixels')
ylabel('Y Pixels')

fig3 = gcf;
savefig(fig3,strcat('detector_',char(datetime('now','Format','MM_dd_hh_mm'))));

fprintf('Figure of detector without reflections created.\n')
%% plot some results
% This creates a cut for all rays which do not leave the target liquid
out_of_cf3i_cut1 = (raytracer_output1(out_of_cf3i_index1).ray_index>0) & ...
    ( (raytracer_output1(out_of_cf3i_index1).surface_index == -1) | ...
    (raytracer_output1(out_of_cf3i_index1).surface_index == -3) | ...
    (raytracer_output1(out_of_cf3i_index1).surface_index == -5) | ...
    (raytracer_output1(out_of_cf3i_index1).surface_index == 8) | ...
    (raytracer_output1(out_of_cf3i_index1).surface_index == 10) | ...
    (raytracer_output1(out_of_cf3i_index1).surface_index == 12) );

out_of_cf3i_cut2 = (raytracer_output2(out_of_cf3i_index2).ray_index>0) & ...
        ( (raytracer_output2(out_of_cf3i_index2).surface_index == -1) | ...
        (raytracer_output2(out_of_cf3i_index2).surface_index == -3) | ...
        (raytracer_output2(out_of_cf3i_index2).surface_index == -5) | ...
        (raytracer_output2(out_of_cf3i_index2).surface_index == 8) | ...
        (raytracer_output2(out_of_cf3i_index2).surface_index == 10) | ...
        (raytracer_output2(out_of_cf3i_index2).surface_index == 12) );
    
% Chooses all rays up until the point where they leave the jar of c3f8/cf3i
cf3i_raylist1 = raytracer_output1(out_of_cf3i_index1).ray_index(out_of_cf3i_cut1);
cf3i_raylist2 = raytracer_output2(out_of_cf3i_index2).ray_index(out_of_cf3i_cut2);

colorlist = 'kgbcyrm';
% First figure gives all surfaces where the c3f8 is
figure(7);
clf;
subplot(1,2,1)
set(gca,'fontsize',16);
for n=1:out_of_cf3i_index1
    % raycut takes only the rays in cf3i_raylist
    raycut1 = ismember(raytracer_output1(n).ray_index,cf3i_raylist1);
    plot3(raytracer_output1(n).intersection_point(raycut1,1), ...
        raytracer_output1(n).intersection_point(raycut1,2), ...
        raytracer_output1(n).intersection_point(raycut1,3), ...
        'o', 'color', colorlist(n), 'markerfacecolor', colorlist(n), ...
        'markersize', 4);
    hold on;
end

subplot(1,2,2)
set(gca,'fontsize',16);
for n=1:out_of_cf3i_index2
    % raycut takes only the rays in cf3i_raylist
    raycut2 = ismember(raytracer_output2(n).ray_index,cf3i_raylist2);
    plot3(raytracer_output2(n).intersection_point(raycut2,1), ...
        raytracer_output2(n).intersection_point(raycut2,2), ...
        raytracer_output2(n).intersection_point(raycut2,3), ...
        'o', 'color', colorlist(n), 'markerfacecolor', colorlist(n), ...
        'markersize', 4);
    hold on;
end

% This takes rays that hit all the other aspects of the inside of the
% pressure vessel such as the inner vessel, the bellows, and the
% retroreflector
quartz_raylist1 = raytracer_output1(out_of_cf3i_index1-2).ray_index( ...
    ((raytracer_output1(out_of_cf3i_index1-2).ray_index>0) & ...
    ( (raytracer_output1(out_of_cf3i_index1-2).surface_index==2) | ...
    (raytracer_output1(out_of_cf3i_index1-2).surface_index==4) | ...
    (raytracer_output1(out_of_cf3i_index1-2).surface_index==6) ) ) | ...
    ((raytracer_output1(out_of_cf3i_index1-2).ray_index>0) & ...
    ((raytracer_output1(out_of_cf3i_index1-2).surface_index==8) | ...
    (raytracer_output1(out_of_cf3i_index1-2).surface_index==10) | ...
    (raytracer_output1(out_of_cf3i_index1-2).surface_index==12) )) );
quartz_raylist1 = quartz_raylist1(~ismember(quartz_raylist1,cf3i_raylist1));
bellows_raylist1 = raytracer_output1(out_of_cf3i_index1).ray_index( ...
    (raytracer_output1(out_of_cf3i_index1).ray_index>0) & ...
    (raytracer_output1(out_of_cf3i_index1).surface_index==13) );
reflector_raylist1 = raytracer_output1(out_of_cf3i_index1).ray_index( ...
    (raytracer_output1(out_of_cf3i_index1).ray_index>0) & ...
    (raytracer_output1(out_of_cf3i_index1).surface_index<=-first_retro_surface1) );

quartz_raylist2 = raytracer_output2(out_of_cf3i_index2-2).ray_index( ...
    ((raytracer_output2(out_of_cf3i_index2-2).ray_index>0) & ...
    ( (raytracer_output2(out_of_cf3i_index2-2).surface_index==2) | ...
    (raytracer_output2(out_of_cf3i_index2-2).surface_index==4) | ...
    (raytracer_output2(out_of_cf3i_index2-2).surface_index==6) ) ) | ...
    ((raytracer_output2(out_of_cf3i_index2-2).ray_index>0) & ...
    ((raytracer_output2(out_of_cf3i_index2-2).surface_index==8) | ...
    (raytracer_output2(out_of_cf3i_index2-2).surface_index==10) | ...
    (raytracer_output2(out_of_cf3i_index2-2).surface_index==12) )) );
quartz_raylist2 = quartz_raylist2(~ismember(quartz_raylist2,cf3i_raylist2));
bellows_raylist2 = raytracer_output2(out_of_cf3i_index2).ray_index( ...
    (raytracer_output2(out_of_cf3i_index2).ray_index>0) & ...
    (raytracer_output2(out_of_cf3i_index2).surface_index==13) );  
reflector_raylist2 = raytracer_output2(out_of_cf3i_index2).ray_index( ...
    (raytracer_output2(out_of_cf3i_index2).ray_index>0) & ...
    (raytracer_output2(out_of_cf3i_index2).surface_index<=-first_retro_surface2) );

reflector_raylist1 = reflector_raylist1(~ismember(reflector_raylist1,cf3i_raylist1));
reflector_raylist1 = reflector_raylist1(~ismember(reflector_raylist1,quartz_raylist1));

reflector_raylist2 = reflector_raylist2(~ismember(reflector_raylist2,cf3i_raylist2));
reflector_raylist2 = reflector_raylist2(~ismember(reflector_raylist2,quartz_raylist2));

if prod(resolution1) ~= size(pixels{1},1) || ...
        (resolution1(1) > 1 && diff(pixels{1}(1:2,1)) ~= 1) || ...
        (resolution1(1) > 1 && diff(pixels{1}(1:2,2)) ~= 0)
    disp('We have a pixel numbering problem');
end

if prod(resolution2) ~= size(pixels{2},1) || ...
        (resolution2(1) > 1 && diff(pixels{2}(1:2,1)) ~= 1) || ...
        (resolution2(1) > 1 && diff(pixels{2}(1:2,2)) ~= 0)
    disp('We have a pixel numbering problem');
end

% The following constructs the image
composite_image1 = ones([resolution1 3]);
r_image1 = zeros(resolution1);
g_image1 = zeros(resolution1);
b_image1 = zeros(resolution1);

r_image1(quartz_raylist1) = 1;
r_image1(bellows_raylist1) = 1;
b_image1(bellows_raylist1) = 1;
b_image1(cf3i_raylist1) = 1;
r_image1(reflector_raylist1) = 1;
g_image1(reflector_raylist1) = 1;
b_image1(reflector_raylist1) = 1;

composite_image1(:,:,1) = (r_image1);
composite_image1(:,:,2) = (g_image1);
composite_image1(:,:,3) = (b_image1);

composite_image2 = ones([resolution2 3]);
r_image2 = zeros(resolution2);
g_image2 = zeros(resolution2);
b_image2 = zeros(resolution2);

r_image2(quartz_raylist2) = 1;
r_image2(bellows_raylist2) = 1;
b_image2(bellows_raylist2) = 1;
b_image2(cf3i_raylist2) = 1;
r_image2(reflector_raylist2) = 1;
g_image2(reflector_raylist2) = 1;
b_image2(reflector_raylist2) = 1;

composite_image2(:,:,1) = (r_image2);
composite_image2(:,:,2) = (g_image2);
composite_image2(:,:,3) = (g_image2);

figure(9);
clf;
subplot(1,2,1)
set(gca,'fontsize',16);
image(permute(composite_image1,[2 1 3]));
axis image
title('Upper Camera')
xlabel('X Pixels')
ylabel('Y Pixels')
subplot(1,2,2)
set(gca,'fontsize',16);
image(permute(composite_image2,[2 1 3]));
axis image
title('Lower Camera')
xlabel('X Pixels')
ylabel('Y Pixels')
fig4 = gcf;
savefig(fig4,strcat('c3f8_cut_',char(datetime('now','Format','MM_dd_hh_mm'))));

%% pixel density figure
% The following section takes the identity of the rays as discussed above
% and calculates the final density of pixels in the final image
cf3i_pixelcut1 = false(resolution1);
cf3i_pixelcut1(cf3i_raylist1) = true;
pixelid1 = 1:numel(cf3i_pixelcut1);

pixelid_down1 = pixelid1 + resolution1(1);
pixelid_down1(pixelid1 > (resolution1(2)-1)*resolution1(1)) = length(pixelid1)+1;

pixelid_up1 = pixelid1 - resolution1(1);
pixelid_up1(pixelid1 <= resolution1(1)) = length(pixelid1)+1;

pixelid_left1 = pixelid1 - 1;
pixelid_left1(mod(pixelid1,resolution1(1))==1) = length(pixelid1)+1;

pixelid_right1 = pixelid1 + 1;
pixelid_right1(mod(pixelid1,resolution1(1))==0) = length(pixelid1)+1;

out_coordinates1 = zeros(length(pixelid1)+1,3);
out_direction1 = zeros(length(pixelid1)+1,3);

out_coordinates1(cf3i_raylist1,:) = raytracer_output1(out_of_cf3i_index1).intersection_point(out_of_cf3i_cut1,:);
left_coordinates1 = out_coordinates1(pixelid_left1,:);
right_coordinates1 = out_coordinates1(pixelid_right1,:);
up_coordinates1 = out_coordinates1(pixelid_up1,:);
down_coordinates1 = out_coordinates1(pixelid_down1,:);
out_coordinates1 = out_coordinates1(1:(end-1),:);

out_direction1(cf3i_raylist1,:) = raytracer_output1(out_of_cf3i_index1).incoming_ray(out_of_cf3i_cut1,1:3);
left_direction1 = out_direction1(pixelid_left1,:);
right_direction1 = out_direction1(pixelid_right1,:);
up_direction1 = out_direction1(pixelid_up1,:);
down_direction1 = out_direction1(pixelid_down1,:);
out_direction1 = out_direction1(1:(end-1),:);

out_dot_left1 = sum(out_direction1 .* left_direction1, 2);
out_dot_right1 = sum(out_direction1 .* right_direction1, 2);
out_dot_up1 = sum(out_direction1 .* up_direction1, 2);
out_dot_down1 = sum(out_direction1 .* down_direction1, 2);

out_dot_left1(out_dot_left1==0) = inf;
out_dot_right1(out_dot_right1==0) = inf;
out_dot_up1(out_dot_up1==0) = inf;
out_dot_down1(out_dot_down1==0) = inf;

new_left1 = left_coordinates1 + left_direction1 .* ...
    repmat(sum(out_direction1 .* (out_coordinates1 - left_coordinates1), 2)./out_dot_left1,1,3);
new_right1 = right_coordinates1 + right_direction1 .* ...
    repmat(sum(out_direction1 .* (out_coordinates1 - right_coordinates1), 2)./out_dot_right1,1,3);
new_up1 = up_coordinates1 + up_direction1 .* ...
    repmat(sum(out_direction1 .* (out_coordinates1 - up_coordinates1), 2)./out_dot_up1,1,3);
new_down1 = down_coordinates1 + down_direction1 .* ...
    repmat(sum(out_direction1 .* (out_coordinates1 - down_coordinates1), 2)./out_dot_down1,1,3);

area1 = zeros(size(out_coordinates1,1),4);
area1(:,1) = sqrt(sum(cross(out_coordinates1-new_up1,out_coordinates1-new_left1).^2,2));
area1(:,2) = sqrt(sum(cross(out_coordinates1-new_up1,out_coordinates1-new_right1).^2,2));
area1(:,3) = sqrt(sum(cross(out_coordinates1-new_down1,out_coordinates1-new_left1).^2,2));
area1(:,4) = sqrt(sum(cross(out_coordinates1-new_down1,out_coordinates1-new_right1).^2,2));

area1(sum(out_coordinates1.^2,2)==0, :) = 0;
area1(sum(up_coordinates1.^2,2)==0, [1 2]) = 0;
area1(sum(down_coordinates1.^2,2)==0, [3 4]) = 0;
area1(sum(left_coordinates1.^2,2)==0, [1 3]) = 0;
area1(sum(right_coordinates1.^2,2)==0, [2 4]) = 0;

mean_area1 = sum(area1,2);
count_area1 = sum(area1>0,2);
mean_area1(mean_area1>0) = mean_area1(mean_area1>0) ./ count_area1(mean_area1>0);

pixels_per_mm2_1 = zeros(size(mean_area1));
pixels_per_mm2_1(mean_area1>0) = .01 ./ mean_area1(mean_area1>0);

% Things above 15 are breaking the code so we put them to 0
set_to_0_ix = pixels_per_mm2_1>10; 
pixels_per_mm2_1(set_to_0_ix) = 0;

ncolors = 200;
binedges1 = linspace(min(pixels_per_mm2_1(pixels_per_mm2_1>0)), max(pixels_per_mm2_1(pixels_per_mm2_1>0)), ncolors);
cmap = hsv(ncolors);
[~, bins] = histc(pixels_per_mm2_1, binedges1);

new_b_image1 = b_image1;
new_g_image1 = g_image1;
new_r_image1 = r_image1;

new_b_image1(pixels_per_mm2_1>0) = cmap(bins(pixels_per_mm2_1>0),3);
new_g_image1(pixels_per_mm2_1>0) = cmap(bins(pixels_per_mm2_1>0),2);
new_r_image1(pixels_per_mm2_1>0) = cmap(bins(pixels_per_mm2_1>0),1);
composite_image1(:,:,1) = (new_r_image1);
composite_image1(:,:,2) = (new_g_image1);
composite_image1(:,:,3) = (new_b_image1);

% Lower Camera Pixel Density
cf3i_pixelcut2 = false(resolution2);
cf3i_pixelcut2(cf3i_raylist2) = true;
pixelid2 = 1:numel(cf3i_pixelcut2);
pixelid_down2 = pixelid2 + resolution2(1);
pixelid_down2(pixelid2 > (resolution2(2)-1)*resolution2(1)) = length(pixelid2)+1;
pixelid_up2 = pixelid2 - resolution2(1);
pixelid_up2(pixelid2 <= resolution2(1)) = length(pixelid2)+1;
pixelid_left2 = pixelid2 - 1;
pixelid_left2(mod(pixelid2,resolution2(1))==1) = length(pixelid2)+1;
pixelid_right2 = pixelid2 + 1;
pixelid_right2(mod(pixelid2,resolution2(1))==0) = length(pixelid2)+1;
out_coordinates2 = zeros(length(pixelid2)+1,3);
out_direction2 = zeros(length(pixelid2)+1,3);
out_coordinates2(cf3i_raylist2,:) = raytracer_output2(out_of_cf3i_index2).intersection_point(out_of_cf3i_cut2,:);
left_coordinates2 = out_coordinates2(pixelid_left2,:);
right_coordinates2 = out_coordinates2(pixelid_right2,:);
up_coordinates2 = out_coordinates2(pixelid_up2,:);
down_coordinates2 = out_coordinates2(pixelid_down2,:);
out_coordinates2 = out_coordinates2(1:(end-1),:);
out_direction2(cf3i_raylist2,:) = raytracer_output2(out_of_cf3i_index2).incoming_ray(out_of_cf3i_cut2,1:3);
left_direction2 = out_direction2(pixelid_left2,:);
right_direction2 = out_direction2(pixelid_right2,:);
up_direction2 = out_direction2(pixelid_up2,:);
down_direction2 = out_direction2(pixelid_down2,:);
out_direction2 = out_direction2(1:(end-1),:);
out_dot_left2 = sum(out_direction2 .* left_direction2, 2);
out_dot_right2 = sum(out_direction2 .* right_direction2, 2);
out_dot_up2 = sum(out_direction2 .* up_direction2, 2);
out_dot_down2 = sum(out_direction2 .* down_direction2, 2);
out_dot_left2(out_dot_left2==0) = inf;
out_dot_right2(out_dot_right2==0) = inf;
out_dot_up2(out_dot_up2==0) = inf;
out_dot_down2(out_dot_down2==0) = inf;
new_left2 = left_coordinates2 + left_direction2 .* ...
    repmat(sum(out_direction2 .* (out_coordinates2 - left_coordinates2), 2)./out_dot_left2,1,3);
new_right2 = right_coordinates2 + right_direction2 .* ...
    repmat(sum(out_direction2 .* (out_coordinates2 - right_coordinates2), 2)./out_dot_right2,1,3);
new_up2 = up_coordinates2 + up_direction2 .* ...
    repmat(sum(out_direction2 .* (out_coordinates2 - up_coordinates2), 2)./out_dot_up2,1,3);
new_down2 = down_coordinates2 + down_direction2 .* ...
    repmat(sum(out_direction2 .* (out_coordinates2 - down_coordinates2), 2)./out_dot_down2,1,3);
area2 = zeros(size(out_coordinates2,1),4);
area2(:,1) = sqrt(sum(cross(out_coordinates2-new_up2,out_coordinates2-new_left2).^2,2));
area2(:,2) = sqrt(sum(cross(out_coordinates2-new_up2,out_coordinates2-new_right2).^2,2));
area2(:,3) = sqrt(sum(cross(out_coordinates2-new_down2,out_coordinates2-new_left2).^2,2));
area2(:,4) = sqrt(sum(cross(out_coordinates2-new_down2,out_coordinates2-new_right2).^2,2));
area2(sum(out_coordinates2.^2,2)==0, :) = 0;
area2(sum(up_coordinates2.^2,2)==0, [1 2]) = 0;
area2(sum(down_coordinates2.^2,2)==0, [3 4]) = 0;
area2(sum(left_coordinates2.^2,2)==0, [1 3]) = 0;
area2(sum(right_coordinates2.^2,2)==0, [2 4]) = 0;
mean_area2 = sum(area2,2);
count_area2 = sum(area2>0,2);
mean_area2(mean_area2>0) = mean_area2(mean_area2>0) ./ count_area2(mean_area2>0);
pixels_per_mm2_2 = zeros(size(mean_area2));
pixels_per_mm2_2(mean_area2>0) = .01 ./ mean_area2(mean_area2>0);
set_to_0_ix = find(pixels_per_mm2_2>15); 
pixels_per_mm2_2(set_to_0_ix) = 0;
ncolors = 200;
binedges2 = linspace(min(pixels_per_mm2_2(pixels_per_mm2_2>0)), max(pixels_per_mm2_2(pixels_per_mm2_2>0)), ncolors);
[~, bins] = histc(pixels_per_mm2_2, binedges2);
new_b_image2 = b_image2;
new_g_image2 = g_image2;
new_r_image2 = r_image2;
new_b_image2(pixels_per_mm2_2>0) = cmap(bins(pixels_per_mm2_2>0),3);
new_g_image2(pixels_per_mm2_2>0) = cmap(bins(pixels_per_mm2_2>0),2);
new_r_image2(pixels_per_mm2_2>0) = cmap(bins(pixels_per_mm2_2>0),1);
composite_image2(:,:,1) = (new_r_image2);
composite_image2(:,:,2) = (new_g_image2);
composite_image2(:,:,3) = (new_b_image2);

% Creating the pixel density picture
figure(10);
clf;
subplot(1,2,1)
set(gca,'fontsize',16);
image(permute(composite_image1,[2 1 3]));
title('Upper Camera Pixel Density');
colormap(cmap);
caxis(binedges1([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges1([1 end]));
set(ch,'Ydata',binedges1([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Pixels per mm^2');
axis image
subplot(1,2,2)
set(gca,'fontsize',16);
image(permute(composite_image2,[2 1 3]));
title('Lower Camera Pixel Density');
colormap(cmap);
caxis(binedges2([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges2([1 end]));
set(ch,'Ydata',binedges2([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Pixels per mm^2');
axis image
fig5 = gcf;
savefig(fig5,strcat('pixels_',char(datetime('now','Format','MM_dd_hh_mm'))));

%% get pixel ray endpoint data

% found_endpoints and final_rays work to take the last scatter of the ray
% i.e. where it hits an absorbing material, and identify the material,
% distance travelled, and other charecteristics of the end point
found_endpoints1 = false(size(pixelid1));

final_rays1 = struct( ...
    'incoming_ray', zeros(length(pixelid1), 10), ...
    'intersection_point', zeros(length(pixelid1), 3), ...
    'surface_normal', zeros(length(pixelid1), 3), ...
    'ray_index', zeros(length(pixelid1), 1), ...
    'surface_index', zeros(length(pixelid1), 1) );

for n=length(raytracer_output1):-1:1
    if all(found_endpoints1)
        break
    end
    
    [pixel_present1, pixel_ix1] = ismember(pixelid1, raytracer_output1(n).ray_index);
    these_endpoints1 = ~found_endpoints1 & pixel_present1;

    found_endpoints1 = found_endpoints1 | these_endpoints1;
    final_rays1.incoming_ray(pixelid1(these_endpoints1),:) = raytracer_output1(n).incoming_ray(pixel_ix1(these_endpoints1),:);
    final_rays1.intersection_point(pixelid1(these_endpoints1),:) = raytracer_output1(n).intersection_point(pixel_ix1(these_endpoints1),:);
    final_rays1.surface_normal(pixelid1(these_endpoints1),:) = raytracer_output1(n).surface_normal(pixel_ix1(these_endpoints1),:);
    final_rays1.ray_index(pixelid1(these_endpoints1)) = raytracer_output1(n).ray_index(pixel_ix1(these_endpoints1));
    final_rays1.surface_index(pixelid1(these_endpoints1)) = raytracer_output1(n).surface_index(pixel_ix1(these_endpoints1));
    
end

% Calculating the final angle of the ray to the surface normal
final_rays1.cos_incident = -sum(final_rays1.incoming_ray(:,1:3) .* final_rays1.surface_normal, 2);
final_rays1.sin_incident = sqrt(1 - final_rays1.cos_incident.^2);
final_rays1.theta_incident = acos(final_rays1.cos_incident);

% CAM LOWER
found_endpoints2 = false(size(pixelid2));
final_rays2 = struct( ...
    'incoming_ray', zeros(length(pixelid2), 10), ...
    'intersection_point', zeros(length(pixelid2), 3), ...
    'surface_normal', zeros(length(pixelid2), 3), ...
    'ray_index', zeros(length(pixelid2), 1), ...
    'surface_index', zeros(length(pixelid2), 1) );
for n=length(raytracer_output2):-1:1
    if all(found_endpoints2)
        break
    end
    
    hits_ring=zeros(size(raytracer_output2(n).surface_index));
    hits_ring(abs(raytracer_output2(n).surface_index)<3 & abs(raytracer_output2(n).surface_index)>0)=1;
    sum(hits_ring)
    
    [pixel_present2, pixel_ix2] = ismember(pixelid2, raytracer_output2(n).ray_index);
    these_endpoints2 = ~found_endpoints2 & pixel_present2;
    found_endpoints2 = found_endpoints2 | these_endpoints2;
    final_rays2.incoming_ray(pixelid2(these_endpoints2),:) = raytracer_output2(n).incoming_ray(pixel_ix2(these_endpoints2),:);
    final_rays2.intersection_point(pixelid2(these_endpoints2),:) = raytracer_output2(n).intersection_point(pixel_ix2(these_endpoints2),:);
    final_rays2.surface_normal(pixelid2(these_endpoints2),:) = raytracer_output2(n).surface_normal(pixel_ix2(these_endpoints2),:);
    final_rays2.ray_index(pixelid2(these_endpoints2)) = raytracer_output2(n).ray_index(pixel_ix2(these_endpoints2));
    final_rays2.surface_index(pixelid2(these_endpoints2)) = raytracer_output2(n).surface_index(pixel_ix2(these_endpoints2));
end
final_rays2.cos_incident = -sum(final_rays2.incoming_ray(:,1:3) .* final_rays2.surface_normal, 2);
final_rays2.sin_incident = sqrt(1 - final_rays2.cos_incident.^2);
final_rays2.theta_incident = acos(final_rays2.cos_incident);

% Building the ray to retroreflector angle figure
ncolors = 200;
cmap = hsv(ncolors);
binedges1 = linspace(0, 90, ncolors);
[~, bins] = histc(final_rays1.theta_incident*180/pi, binedges1);
bins(bins==0) = length(binedges1);

new_b_image1 = b_image1;
new_g_image1 = g_image1;
new_r_image1 = r_image1;
new_r_image1(new_r_image1==1) = .5;

new_b_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),3);
new_g_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),2);
new_r_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),1);

miss_retroreflector1 = abs(final_rays1.surface_index)<first_retro_surface1;
    
new_b_image1(cf3i_pixelcut1(:) & miss_retroreflector1(:)) = 0;
new_g_image1(cf3i_pixelcut1(:) & miss_retroreflector1(:)) = 0;
new_r_image1(cf3i_pixelcut1(:) & miss_retroreflector1(:)) = 0;

composite_image1(:,:,1) = (new_r_image1);
composite_image1(:,:,2) = (new_g_image1);
composite_image1(:,:,3) = (new_b_image1);

%Cam Lower
binedges2 = linspace(0, 90, ncolors);
[~, bins] = histc(final_rays2.theta_incident*180/pi, binedges2);
bins(bins==0) = length(binedges2);
new_b_image2 = b_image2;
new_g_image2 = g_image2;
new_r_image2 = r_image2;
new_r_image2(new_r_image2==1) = .5;
new_b_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),3);
new_g_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),2);
new_r_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),1);
miss_retroreflector2 = abs(final_rays2.surface_index)<first_retro_surface2;
new_b_image2(cf3i_pixelcut2(:) & miss_retroreflector2(:)) = 0;
new_g_image2(cf3i_pixelcut2(:) & miss_retroreflector2(:)) = 0;
new_r_image2(cf3i_pixelcut2(:) & miss_retroreflector2(:)) = 0;
composite_image2(:,:,1) = (new_r_image2);
composite_image2(:,:,2) = (new_g_image2);
composite_image2(:,:,3) = (new_b_image2);

figure(11);
clf;
subplot(1,2,1)
image(permute(composite_image1,[2 1 3]));

set(gca,'fontsize',16);
title('Ray to Retroreflector Angle in Upper Camera');

colormap(cmap);
axis image
caxis(binedges1([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges1([1 end]));
set(ch,'Ydata',binedges1([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Incident Angle (degrees)');

subplot(1,2,2)
image(permute(composite_image2,[2 1 3]));

set(gca,'fontsize',16);
title('Ray to Retroreflector Angle in Lower Camera');

colormap(cmap);
axis image
caxis(binedges2([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges2([1 end]));
set(ch,'Ydata',binedges2([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Incident Angle (degrees)');
fig6 = gcf;
savefig(fig6,strcat('angle_',char(datetime('now','Format','MM_dd_hh_mm'))));

%% Build corrected amplitude figure
% Uses retroreflector data to correct intensity for different incident
% angles

ncolors = 200;
binedges1 = linspace(0, 1, ncolors);
cmap = cool(ncolors);
retro_mu = 0;
retro_sigma = 61.44529266;
[~, bins] = histc(final_rays1.incoming_ray(:,7).*exp(-(retro_mu-final_rays1.theta_incident*180/pi).^2/(2*retro_sigma^2)), binedges1);
bins(bins==0) = length(binedges1);

new_b_image1 = b_image1;
new_g_image1 = g_image1;
new_r_image1 = r_image1;

new_b_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),3);
new_g_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),2);
new_r_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),1);
composite_image1(:,:,1) = (new_r_image1);
composite_image1(:,:,2) = (new_g_image1);
composite_image1(:,:,3) = (new_b_image1);

%Lower cam
binedges2 = linspace(0, 1, ncolors);
[~, bins] = histc(final_rays2.incoming_ray(:,7).*exp(-(retro_mu-final_rays2.theta_incident*180/pi).^2/(2*retro_sigma^2)), binedges2);
bins(bins==0) = length(binedges2);
new_b_image2 = b_image2;
new_g_image2 = g_image2;
new_r_image2 = r_image2;
new_b_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),3);
new_g_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),2);
new_r_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),1);
composite_image2(:,:,1) = (new_r_image2);
composite_image2(:,:,2) = (new_g_image2);
composite_image2(:,:,3) = (new_b_image2);

figure(33);
clf;
subplot(1,2,1)
image(permute(composite_image1,[2 1 3]));

set(gca,'fontsize',16);
title('Corrected Intensity Hitting Retroreflector');

colormap(cmap);
axis image
caxis(binedges1([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges1([1 end]));
set(ch,'Ydata',binedges1([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Reflected Amplitude (0..1)');

subplot(1,2,2)
image(permute(composite_image2,[2 1 3]));
set(gca,'fontsize',16);
title('Corrected Intensity Hitting Retroreflector');
colormap(cmap);
axis image
caxis(binedges2([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges2([1 end]));
set(ch,'Ydata',binedges2([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Reflected Amplitude (0..1)');
fig7 = gcf;
savefig(fig7,strcat('intensity_',char(datetime('now','Format','MM_dd_hh_mm'))));

% Next we build the ray amplitude figure, this gives the amplitude of the
% ray when it hits the retroreflector, as a ratio of starting amplitude and
% final amplitude
ncolors = 200;
binedges1 = linspace(0, 1, ncolors);
cmap = cool(ncolors);
[~, bins] = histc(final_rays1.incoming_ray(:,7), binedges1);
bins(bins==0) = length(binedges1);

new_b_image1 = b_image1;
new_g_image1 = g_image1;
new_r_image1 = r_image1;

new_b_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),3);
new_g_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),2);
new_r_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),1);
composite_image1(:,:,1) = (new_r_image1);
composite_image1(:,:,2) = (new_g_image1);
composite_image1(:,:,3) = (new_b_image1);

%Lower cam
binedges2 = linspace(0, 1, ncolors);
[~, bins] = histc(final_rays2.incoming_ray(:,7), binedges2);
bins(bins==0) = length(binedges2);
new_b_image2 = b_image2;
new_g_image2 = g_image2;
new_r_image2 = r_image2;
new_b_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),3);
new_g_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),2);
new_r_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),1);
composite_image2(:,:,1) = (new_r_image2);
composite_image2(:,:,2) = (new_g_image2);
composite_image2(:,:,3) = (new_b_image2);

figure(12);
clf;
subplot(1,2,1)
image(permute(composite_image1,[2 1 3]));

set(gca,'fontsize',16);
title('Ray Amplitude Hitting Retroreflector');

colormap(cmap);
axis image
caxis(binedges1([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges1([1 end]));
set(ch,'Ydata',binedges1([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Incident Amplitude (0..1)');

subplot(1,2,2)
image(permute(composite_image2,[2 1 3]));
set(gca,'fontsize',16);
title('Ray Amplitude Hitting Retroreflector');
colormap(cmap);
axis image
caxis(binedges2([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges2([1 end]));
set(ch,'Ydata',binedges2([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Incident Amplitude (0..1)');
%%

ncolors = 2*length(surface_list) + 1;
binedges1 = (-.5-length(surface_list)):(.5+length(surface_list));
cmap = rand(ncolors,3);
[~, bins] = histc(final_rays1.surface_index, binedges1);
bins(bins==0) = length(binedges1);

new_b_image1 = b_image1;
new_g_image1 = g_image1;
new_r_image1 = r_image1;

new_b_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),3);
new_g_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),2);
new_r_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),1);
composite_image1(:,:,1) = (new_r_image1);
composite_image1(:,:,2) = (new_g_image1);
composite_image1(:,:,3) = (new_b_image1);
figure(21);
clf;
image(permute(composite_image1,[2 1 3]));

set(gca,'fontsize',16);
title('Ray Amplitude Hitting Retroreflector');

colormap(cmap);
axis image
caxis(binedges1([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges1([1 end]));
set(ch,'Ydata',binedges1([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Final Surface (index)');


%% get ray lengths

% These last calculations give the length of the ray from leaving the
% camera to when it hits the retroreflector and then builds the figure for
% it
ray_length1 = zeros(size(pixelid1));
ray_optical_length1 = zeros(size(pixelid1));
ray_angular_length1 = zeros(size(pixelid1));
for n=1:length(raytracer_output1)
    [pixel_present1, pixel_ix1] = ismember(pixelid1, raytracer_output1(n).ray_index);
    ray_length1(pixel_present1) = ray_length1(pixel_present1) + raytracer_output1(n).distance_traveled(pixel_ix1(pixel_present1))';
    ray_optical_length1(pixel_present1) = ray_optical_length1(pixel_present1) + ...
        (raytracer_output1(n).distance_traveled(pixel_ix1(pixel_present1))' .* raytracer_output1(n).n_incident(pixel_ix1(pixel_present1))');
    ray_angular_length1(pixel_present1) = ray_angular_length1(pixel_present1) + ...
        (raytracer_output1(n).distance_traveled(pixel_ix1(pixel_present1))' ./ raytracer_output1(n).n_incident(pixel_ix1(pixel_present1))');
end

ncolors = 200;
binedges1 = linspace(min(ray_angular_length1(cf3i_pixelcut1(:))), max(ray_angular_length1(cf3i_pixelcut1(:))), ncolors);
cmap = cool(ncolors);
[~, bins] = histc(ray_angular_length1, binedges1);
bins(bins==0) = length(binedges1);

new_b_image1 = b_image1;
new_g_image1 = g_image1;
new_r_image1 = r_image1;

new_b_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),3);
new_g_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),2);
new_r_image1(cf3i_pixelcut1(:)) = cmap(bins(cf3i_pixelcut1(:)),1);
composite_image1(:,:,1) = (new_r_image1);
composite_image1(:,:,2) = (new_g_image1);
composite_image1(:,:,3) = (new_b_image1);

%Lower cam
ray_length2 = zeros(size(pixelid2));
ray_optical_length2 = zeros(size(pixelid2));
ray_angular_length2 = zeros(size(pixelid2));
for n=1:length(raytracer_output2)
    [pixel_present2, pixel_ix2] = ismember(pixelid2, raytracer_output2(n).ray_index);
    ray_length2(pixel_present2) = ray_length2(pixel_present2) + raytracer_output2(n).distance_traveled(pixel_ix2(pixel_present2))';
    ray_optical_length2(pixel_present2) = ray_optical_length2(pixel_present2) + ...
        (raytracer_output2(n).distance_traveled(pixel_ix2(pixel_present2))' .* raytracer_output2(n).n_incident(pixel_ix2(pixel_present2))');
    ray_angular_length2(pixel_present2) = ray_angular_length2(pixel_present2) + ...
        (raytracer_output2(n).distance_traveled(pixel_ix2(pixel_present2))' ./ raytracer_output2(n).n_incident(pixel_ix2(pixel_present2))');
end
binedges2 = linspace(min(ray_angular_length2(cf3i_pixelcut2(:))), max(ray_angular_length2(cf3i_pixelcut2(:))), ncolors);
[~, bins] = histc(ray_angular_length2, binedges2);
bins(bins==0) = length(binedges2);
new_b_image2 = b_image2;
new_g_image2 = g_image2;
new_r_image2 = r_image2;
new_b_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),3);
new_g_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),2);
new_r_image2(cf3i_pixelcut2(:)) = cmap(bins(cf3i_pixelcut2(:)),1);
composite_image2(:,:,1) = (new_r_image2);
composite_image2(:,:,2) = (new_g_image2);
composite_image2(:,:,3) = (new_b_image2);

figure(13);
clf;
subplot(1,2,1)
image(permute(composite_image1,[2 1 3]));

set(gca,'fontsize',16);
title('Angular Distance to Retroreflector');

colormap(cmap);
axis image
caxis(binedges1([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges1([1 end]));
set(ch,'Ydata',binedges1([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Optical Distance to Retroreflector (cm)');

subplot(1,2,2)
image(permute(composite_image2,[2 1 3]));

set(gca,'fontsize',16);
title('Angular Distance to Retroreflector');

colormap(cmap);
axis image
caxis(binedges2([1 end]));
cb = colorbar;
ch = get(cb,'children');
set(cb,'ylim',binedges2([1 end]));
set(ch,'Ydata',binedges2([1 end]));
set(cb,'fontsize',16);
ylabel(cb,'Optical Distance to Retroreflector (cm)');
fig8 = gcf;
savefig(fig8,strcat('distance_',char(datetime('now','Format','MM_dd_hh_mm'))));

min_ray_optical_length1 = min(ray_optical_length1(cf3i_pixelcut1(:)));
min_ray_optical_length2 = min(ray_optical_length2(cf3i_pixelcut2(:)));
