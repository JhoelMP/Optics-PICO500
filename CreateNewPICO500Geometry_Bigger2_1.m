% function [surface_list rays ray_startingpoints pixels] = CreatePICO40Geometry()
%
% This function creates a structure array of surfaces to be used by
% RayTracer.  Follow this architecture to create any geometry you like.
%
% Each surface has six fields.  The first (intersect_function) defines the
% geometry, and is a function handle to an anonymous function, calling a
% RayToXXXXX function with the appropriate geometry inputs.  For example, 
%
%  @(sp,indir)RayToCylinder(sp,indir, [0 0 0], [0 0 1], 10) 
%
% defines a cylinder on the z-axis with radius 10.  See all RayToXXXXX
% functions in the RayTracing directory for other possible shapes (and
% create your own if you desire).
%
% The second field (inbounds_function) defines the bounds of the surface,
% and is a function handle to an anonymous function that inputs an N-by-3-by-M 
% array and outputs an N-by-M logical.  It can be assumed that all input
% points are on the surface defined by intersect_function, giving true if
% the input point is contained within the bounds, false otherwise.  For
% example,
%
%  @(p)(reshape( ...
%      (p(:,3,:)>20) & (p(:,3,:)<80) & (atan2(p(:,2,:),p(:,1,:))>0), ...
%      size(p,1), [] ));
%
% would cut the above cylinder in half along the xz plane and truncate it 
% at 20 and 80 in the z-coordinate.  (Generically, p is an N-by-3-by-M
% matrix, and the output of the function should be an N-by-M boolean.)
%
% The third and fourth fields are n_outside and n_inside, and give the
% indices of refraction on the two sides of the surface.  What is 'inside'
% and 'outside' is defined in the RayToXXXXX function -- for spheres and
% cylinders, this is obvious, for planes less so.  See the documentation
% for each RayToXXXXX function for details.  Also, setting n to inf makes
% that side a perfect conductor (in terms of calculating reflection and
% polarization, at least).
%
% The fifth field is surface type, and may be 'normal', 'diffuse', or
% 'retro'.  For 'diffuse' and 'retro' surfaces, the normal direction
% returned by the intersect_function is replaced by a random direction
% within pi/4 of normal or the reverse of the incoming ray, respectively.
%
% The sixth field is an absorption coefficient -- RayTracer will multiply
% the intensity of both the reflected and refracted rays coming from this
% surface by 1-absorption.
%
% 12/16/09, CED

%% Geospecs
geospecs = struct();
% PICO 40 geospecs, all specs are in cm unless otherwise specified, these
% dimensions are taken from or estimated from docdb entries #1867, #1537,
% and #1903
% All commented values are the previous values

% New dimensions have been measured from the Solidworks assembly

% Indices of refraction
geospecs.n_jar = 1.456; % checked quartz on online database
geospecs.n_target = 1.16;
%geospecs.n_hydraulic = 1.434;
geospecs.n_hydraulic = 1.4418; % mineral oil, not hydraulic fluid
geospecs.n_water = 1.33;
geospecs.n_air = 1;
geospecs.n_pressurewindow = 1.456;

% Note that all positions are relative to the xy-middle of the pressure
% vessel on the top of the base flange as (0,0,0) (x,y,z)
% The y-axis passes through the center of one set of viewports
% Positive y is away from the viewports
% The x-axis is perpendicular to the y-axis such that the xy-plane is the
% plane of the top of the base flange
% Positive x is to the right after aligning according to positive y
% The z-axis is the cylindrical center of all vessels
% Positve z is upwards from the base flange

% Outer inner vessel specs
geospecs.outervessel_cylrad  = 23*2; %460mm/2=23 cm
geospecs.outervessel_cylthick = 0.7*2;
geospecs.outervessel_knucklerad = 11*2;
geospecs.outervessel_sphererad = 46*2;
%geospecs.outervessel_cyllength =%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geospecs.outervessel_cyllength = 78*2.54;
% This is the z-coordinate to start the cylinder it takes into account the
% size of the bellows and such
%geospecs.outervessel_startcyl = 60.08; %(22.63*2.54) + 2.553;
geospecs.outervessel_startcyl = (9.98+5.61)*2.54;  %?????????????????????????

% Inner inner vessel specs
geospecs.innervessel_cylrad  = 21.4*2;
geospecs.innervessel_cylthick = 0.7*2;
geospecs.innervessel_knucklerad = 10*2;
geospecs.innervessel_sphererad = 42.8*2;
%geospecs.innervessel_cyllength =%%%%%%%%%%%%%%%%%%%%%%%%
geospecs.innervessel_cyllength = 57.9;
% This is the z-coordinate to start the cylinder it takes into account the
% size of the bellows and such
%geospecs.innervessel_startcyl = 32.809; %(36.7*2.54) - 67.3 + 3.823;
geospecs.innervessel_startcyl =0;  %?????????????????????????


% Pressure vessel specs
geospecs.pv_outerrad = 94.66/2*2.54;
%geospecs.pv_cyllength = 69.8*2.54+6.21*2.54;
geospecs.pv_cyllength = 115.96*2.54;
%geospecs.pv_topdomerad = ((geospecs.pv_outerrad .^2) + (15^2))/(2*15);
geospecs.pv_topdomerad = (96+9.98+9.98)*2.54;
geospecs.pv_domeheight =(9.98)*2.54;
geospecs.reflector_topdomerad = (((geospecs.pv_outerrad - 3) .^2) + (12^2))/(2*12);

inclinacion=0

geospecs.cam_f = 0.5;

% Lastly we put in the cameras pitch, yaw, and roll as found from the
% technical drawings (docdb 1867) for the top and bottom cameras in radians
geospecs.camupper_pitch = deg2rad(-inclinacion); %geospecs.vp_upperangle;
geospecs.camupper_yaw = deg2rad(0);
geospecs.camupper_roll = deg2rad(0);
geospecs.camlower_pitch = deg2rad(inclinacion); %geospecs.vp_lowerangle;
geospecs.camlower_yaw = deg2rad(0);
geospecs.camlower_roll = deg2rad(0);





% Viewport specsdeg
geospecs.vp_outerrad = (6 * 2.54);
geospecs.vp_upperangle = deg2rad(-inclinacion);
geospecs.vp_lowerangle = deg2rad(inclinacion);
%geospecs.vp_upperheight = (47.25*2.54) + (((9.02+18)*2.54)*tan(deg2rad(14)));
%geospecs.vp_lowerheight = (47.25*2.54) - (((9.02+18)*2.54)*tan(deg2rad(14)));
%geospecs.vp_upperheight = 200 + geospecs.pv_outerrad*tan(deg2rad(14));
%geospecs.vp_lowerheight =  200- geospecs.pv_outerrad*tan(deg2rad(14));



geospecs.vp_upperheight=  74*2.54;
geospecs.vp_lowerheight=  42*2.54;


geospecs.vp_phi = deg2rad(90);
geospecs.vp_winthick = 4.921; %2*2.54;

% General specs
% This is the height from the top of the pressure vessel base flange to the
% top of the outer inner vessel
geospecs.totalinnerheight = 61.71*2.54; % this doesn't appear anywhere in the PICO40L code
% bellowsrad is used in creating the geometry such that the bellows absorb
% everything going into it, thus we use the outer inner vessel radius
geospecs.bellowsrad = geospecs.outervessel_cylrad;

%Jugar con barrel,lens type, geospec.cam_f, etc


% Next we put down the physical position of the camera, note that this is
% where the lens meets the viewport
geospecs.camupper_x = 0;
geospecs.camupper_y = (-geospecs.pv_outerrad) + (-3.32 * 2.54);
%geospecs.camupper_z = (47.25*2.54) + 7.56 * 2.54;
geospecs.camupper_z = geospecs.vp_upperheight + ...
    (geospecs.vp_outerrad*tan(-geospecs.vp_upperangle)+geospecs.vp_winthick)*sin(-geospecs.vp_upperangle);
geospecs.camlower_x = 0;
geospecs.camlower_y = (-geospecs.pv_outerrad) + (-3.32 * 2.54);
%geospecs.camlower_z = (47.25*2.54) - 7.56 * 2.54;
geospecs.camlower_z = geospecs.vp_lowerheight - ...
    (geospecs.vp_outerrad*tan(geospecs.vp_lowerangle)+geospecs.vp_winthick)*sin(geospecs.vp_lowerangle);
geospecs.camupper_zoffset = 0;
geospecs.camlower_zoffset = 0;

%% Moving onwards and making the list of all surfaces1
[surface_list] = CreateNewPICO40LGeometry1(geospecs);

fprintf('Geospecs structure created\n')
 
function [surface_list] = CreateNewPICO40LGeometry1(geospecs)   %Name Fixed according to the file name

%% set defaults
if nargin<1 || isempty(geospecs) || ~isstruct(geospecs)
    geospecs = struct();
end

%% define surface_list structure
surface_list = struct( ...
    'description', {}, ...
    'intersect_function', {}, ...
    'inbounds_function', {}, ...
    'n_outside', {}, ...
    'n_inside', {}, ...
    'surface_type', {}, ...
    'absorption', {});

%% indices of refraction and dimensions used below


%% extract geospecs
fn = fieldnames(geospecs);
for n=1:length(fn)
    if ~isempty(geospecs.(fn{n}))
        eval([fn{n} '=geospecs.(fn{n});']);
    end
end

%% derived dimensions

% The following variables give some distances for the outer and then the inner vessel as a 1x2
% array where the the 1st element is the outer distance and the 2nd element
% is the inside distance

% t gives the change to radii, 0 for outside, jar thickness for inside
touter = [0 outervessel_cylthick];

% Adjusted cylinder radii
router1 = outervessel_cylrad - touter;
% Adjusted knuckle radii
router2 = outervessel_knucklerad - touter;
% Adjusted sphere radii
router3 = outervessel_sphererad - touter;

% s gives the distance from the cylinder axis to the point where the torus
% and sphere meet
souter = router3.*(router1-router2)./(router3-router2);

% z gives the z-distance from the end of the cylinder to the point where the
% sphere and the torus meet
zouter = router2 .* sqrt(1 - (souter./router3).^2);

% d gives the z-distance from the end of the cylinder to the point where
% the center of the sphere is
douter = - router3 .* zouter .* ((1./router3)-(1./router2));

% t gives the change to radii, 0 for outside, jar thickness for inside
tinner = [0 innervessel_cylthick];

% Adjusted cylinder radii
rinner1 = innervessel_cylrad - tinner;
% Adjusted knuckle radii
rinner2 = innervessel_knucklerad - tinner;
% Adjusted sphere radii
rinner3 = innervessel_sphererad - tinner;

% s gives the distance from the cylinder axis to the point where the torus
% and sphere meet
sinner = rinner3.*(rinner1-rinner2)./(rinner3-rinner2);

% z gives the z-distance from the end of the cylinder to the point where the
% sphere and the torus meet
zinner = rinner2 .* sqrt(1 - (sinner./rinner3).^2);

% d gives the z-distance from the end of the cylinder to the point where
% the center of the sphere is
dinner = - rinner3 .* zinner .* ((1./rinner3)-(1./rinner2));

%% surface list

% Each section through here is its own individual surface, see the
% explanation at the beginning of the function for an explanation of what
% each surface_list.XXXXX does

%1
surface_list(end+1).description = 'inside surface of quartz outer vessel cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, [0 0 0], [0 0 1], outervessel_cylrad - outervessel_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=outervessel_startcyl) & (p(:,3,:)<=(outervessel_startcyl+outervessel_cyllength)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%2
surface_list(end+1).description = 'outside surface of quartz outer vessel cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, [0 0 0], [0 0 1], outervessel_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=outervessel_startcyl) & (p(:,3,:)<=(outervessel_startcyl+outervessel_cyllength)), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%3
surface_list(end+1).description = 'inside surface of quartz outer vessel hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 (outervessel_startcyl + outervessel_cyllength - douter(2))], router3(2)); %rinner changed for router
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(outervessel_startcyl+outervessel_cyllength+zouter(2))), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%4
surface_list(end+1).description = 'outside surface of quartz outer vessel hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 (outervessel_startcyl + outervessel_cyllength - douter(1))], router3(1));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(outervessel_startcyl+outervessel_cyllength+zouter(1))), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%5
surface_list(end+1).description = 'inside surface of quartz outer vessel knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 (outervessel_startcyl + outervessel_cyllength)], [0 0 1], outervessel_cylrad-outervessel_knucklerad, outervessel_knucklerad-outervessel_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(outervessel_startcyl+outervessel_cyllength) ...
	& p(:,3,:)<=(outervessel_startcyl+outervessel_cyllength+zouter(2)) & (p(:,1,:).^2+p(:,2,:).^2)>((outervessel_cylrad-outervessel_knucklerad)^2) ), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;


%6
surface_list(end+1).description = 'outside surface of quartz outer vessel knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 (outervessel_startcyl + outervessel_cyllength)], [0 0 1], outervessel_cylrad-outervessel_knucklerad, outervessel_knucklerad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(outervessel_startcyl+outervessel_cyllength) ...
	& p(:,3,:)<=(outervessel_startcyl+outervessel_cyllength+zouter(1)) & (p(:,1,:).^2+p(:,2,:).^2)>((outervessel_cylrad-outervessel_knucklerad)^2) ), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%7
surface_list(end+1).description = 'inside surface of quartz inner vessel cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, [0 0 0], [0 0 1], innervessel_cylrad - innervessel_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=innervessel_startcyl) & (p(:,3,:)<=(innervessel_startcyl+innervessel_cyllength)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%8
surface_list(end+1).description = 'outside surface of quartz inner vessel cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, [0 0 0], [0 0 1], innervessel_cylrad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=innervessel_startcyl) & (p(:,3,:)<=(innervessel_startcyl+innervessel_cyllength)), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%9
surface_list(end+1).description = 'inside surface of quartz inner vessel hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 (innervessel_startcyl + innervessel_cyllength - dinner(2))], rinner3(2));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(innervessel_startcyl+innervessel_cyllength+zinner(2))), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%10
surface_list(end+1).description = 'outside surface of quartz inner vessel hemisphere';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0 0 (innervessel_startcyl + innervessel_cyllength - dinner(1))], rinner3(1));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(innervessel_startcyl+innervessel_cyllength+zinner(1))), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%11
surface_list(end+1).description = 'inside surface of quartz inner vessel knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 (innervessel_startcyl + innervessel_cyllength)], [0 0 1], innervessel_cylrad-innervessel_knucklerad, innervessel_knucklerad-innervessel_cylthick);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(innervessel_startcyl+innervessel_cyllength) ...
	& p(:,3,:)<=(innervessel_startcyl+innervessel_cyllength+zinner(2)) & (p(:,1,:).^2+p(:,2,:).^2)>((innervessel_cylrad-innervessel_knucklerad)^2) ), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%12
surface_list(end+1).description = 'outside surface of quartz inner vessel knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 (innervessel_startcyl + innervessel_cyllength)], [0 0 1], innervessel_cylrad-innervessel_knucklerad, innervessel_knucklerad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(innervessel_startcyl+innervessel_cyllength) ...
	& p(:,3,:)<=(innervessel_startcyl+innervessel_cyllength+zinner(1)) & (p(:,1,:).^2+p(:,2,:).^2)>((innervessel_cylrad-innervessel_knucklerad)^2) ), size(p,1), [] ));
surface_list(end).n_outside = n_target;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%13
surface_list(end+1).description = 'Bellows cylinder (approx)';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], bellowsrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (p(:,3,:) < outervessel_startcyl) & ...
    (p(:,3,:) > 0) ), ...
    size(p,1), []));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = 1;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;


%14
% Note the excessive coding here takes into account the angled viewport
% entry into the pressure vessel and cuts out holes for it
surface_list(end+1).description = 'PV - cylinder wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], pv_outerrad);
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,3,:)>0) & (p(:,3,:)<(pv_cyllength + pv_domeheight - (pv_topdomerad - sqrt(pv_topdomerad^2-pv_outerrad^2)))) ) &  ...
    ((((p(:,3,:)-vp_upperheight + ...
	((pv_outerrad - (((p(:,1,:) .^ 2) + (pv_outerrad .^2)) .^ (1/2))) .* (tan(vp_upperangle)))).^2 + ...
	(p(:,1,:).^2)) >= ((vp_outerrad/(cos(vp_upperangle)))^2)) | (p(:,2,:)>0)) & ...
	((((p(:,3,:)-vp_lowerheight + ...
	((pv_outerrad - (((p(:,1,:) .^ 2) + (pv_outerrad .^2)) .^ (1/2))) .* (tan(vp_lowerangle)))).^2 + ...
	(p(:,1,:).^2)) >= ((vp_outerrad/(cos(vp_lowerangle)))^2)) | (p(:,2,:)>0)) & ...
	((((p(:,3,:)-vp_upperheight + ...
	((pv_outerrad - (((p(:,1,:) .^ 2) + (pv_outerrad .^2)) .^ (1/2))) .* (tan(vp_upperangle)))).^2 + ...
	(p(:,1,:).*cos(vp_phi) + p(:,2,:).*sin(vp_phi)).^2) >= ((vp_outerrad/(cos(vp_upperangle)))^2)) | ...
	((p(:,2,:).*cos(vp_phi) - p(:,1,:).*sin(vp_phi))>0))& ...
	((((p(:,3,:)-vp_lowerheight + ...
	((pv_outerrad - (((p(:,1,:) .^ 2) + (pv_outerrad .^2)) .^ (1/2))) .* (tan(vp_lowerangle)))).^2 + ...
	(p(:,1,:).*cos(vp_phi) + p(:,2,:).*sin(vp_phi)).^2) >= ((vp_outerrad/(cos(vp_lowerangle)))^2)) | ...
	((p(:,2,:).*cos(vp_phi) - p(:,1,:).*sin(vp_phi))>0)), ...
    size(p,1), []));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%15
% Note the top dome is simplified to the a hemisphere on top of the vessel
surface_list(end+1).description = 'PV - top dome';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    [0,0,(pv_cyllength + pv_domeheight - (geospecs.pv_topdomerad))],geospecs.pv_topdomerad);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>(pv_cyllength + pv_domeheight - (pv_topdomerad - sqrt(pv_topdomerad^2-pv_outerrad^2)))),  ...
    size(p,1), []));
surface_list(end).n_outside = n_water;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%16
% Note the bottom dome is simplified to the a flat plane which totally absorbs
% all light, it is placed as the top of the base flange, i.e. the plane
% containing (0,0,0)

surface_list(end+1).description = 'PV - bottom dome';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0,0,0],[0,0,1]);
surface_list(end).inbounds_function = @(p)(reshape( (((p(:,1,:) .^2) + (p(:,2,:).^2)) <= (pv_outerrad.^2)),  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_water;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%17

surface_list(end+1).description = 'PV - Copper ring';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0,0,(innervessel_startcyl+innervessel_cyllength)],[0,0,1]);
surface_list(end).inbounds_function = @(p)(reshape( (((p(:,1,:) .^2) + (p(:,2,:).^2)) <= ((pv_outerrad).^2) & ((p(:,1,:) .^2) + (p(:,2,:).^2)) >= ((outervessel_cylrad + 1).^2)),  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;


%18
%{
surface_list(end+1).description = 'PV - Copper ring cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], (outervessel_cylrad + 1));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=(outervessel_startcyl-2)) & (p(:,3,:)<=(innervessel_startcyl+innervessel_cyllength)),  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = inf;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
%}
%19
surface_list(end+1).description = 'VP1upper - window casing';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 camupper_y camupper_z], [0 1 tan(vp_upperangle)], vp_outerrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,3,:) - (camupper_z - (vp_outerrad .* cos(vp_upperangle)))) <= (p(:,2,:) - (camupper_y + (vp_outerrad .* sin(vp_upperangle))))/(tan(-vp_upperangle)))   & ...
	(p(:,2,:) <= 0) & p(:,2,:) >= (camupper_y + (vp_outerrad .* sin(vp_upperangle))) & ...
    (((p(:,1,:) .^2) + (p(:,2,:) .^2)) >= (pv_outerrad .^ 2)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%20
surface_list(end+1).description = 'VP1upper - glass-air interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 camupper_y camupper_z], [0 1 tan(vp_upperangle)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)) .^ 2) + ((p(:,2,:) - camupper_y) .^ 2) + ((p(:,3,:) - camupper_z) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%21
surface_list(end+1).description = 'VP1upper - glycol-glass interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 (camupper_y + (vp_winthick*cos(geospecs.vp_upperangle))) (camupper_z + (vp_winthick*sin(geospecs.vp_upperangle)))], [0 1 tan(vp_upperangle)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)) .^ 2) + ((p(:,2,:) - (camupper_y + (vp_winthick*cos(geospecs.vp_upperangle)))) .^ 2) + ((p(:,3,:) - (camupper_z + (vp_winthick*sin(geospecs.vp_upperangle)))) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%22
surface_list(end+1).description = 'VP1lower - window casing';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 camlower_y camlower_z], [0 1 tan(vp_lowerangle)], vp_outerrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((p(:,3,:) - (camlower_z + (vp_outerrad .* cos(vp_lowerangle)))) >= (p(:,2,:) - (camlower_y - (vp_outerrad .* sin(vp_lowerangle))))/(tan(-vp_lowerangle)))   & ...
	(p(:,2,:) <= 0) & ...
    (((p(:,1,:) .^2) + (p(:,2,:) .^2)) >= (pv_outerrad .^ 2)) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%23
surface_list(end+1).description = 'VP1lower - glass-air interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 camlower_y camlower_z], [0 1 tan(vp_lowerangle)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)) .^ 2) + ((p(:,2,:) - camlower_y) .^ 2) + ((p(:,3,:) - camlower_z) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%24
surface_list(end+1).description = 'VP1lower - glycol-glass interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 (camlower_y + (vp_winthick*cos(vp_lowerangle))) (camlower_z + (vp_winthick*sin(vp_lowerangle)))], [0 1 tan(vp_lowerangle)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)) .^ 2) + ((p(:,2,:) - (camlower_y + (vp_winthick*cos(vp_lowerangle)))) .^ 2) + ((p(:,3,:) - (camlower_z + (vp_winthick*sin(vp_lowerangle)))) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%25
surface_list(end+1).description = 'VP2upper - window casing';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, [(camupper_y * sin(-vp_phi)) (camupper_y * cos(vp_phi)) camupper_z], ...
	[(cos(vp_upperangle)*sin(-vp_phi)) (cos(vp_upperangle)*cos(vp_phi)) (sin(vp_upperangle))], vp_outerrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
	(((cos(vp_upperangle)*sin(-vp_phi)) .* (p(:,1,:))) + ((cos(vp_upperangle)*cos(vp_phi)) .* (p(:,2,:))) + ((sin(vp_upperangle)) .* (p(:,3,:))) >= ...
	((cos(vp_upperangle)*sin(-vp_phi)) .* (camupper_y * sin(-vp_phi))) + ((cos(vp_upperangle)*cos(vp_phi)) .* (camupper_y * cos(vp_phi))) + ((sin(vp_upperangle)) .* (camupper_z))) & ... 	
	(p(:,2,:) <= 0) & p(:,2,:) >= (camupper_y + (vp_outerrad .* sin(vp_upperangle))) & ...
    (((p(:,1,:) .^2) + (p(:,2,:) .^2)) >= (pv_outerrad .^ 2)) ), ...
     size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%26
surface_list(end+1).description = 'VP2upper - glass-air interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [(camupper_y * sin(-vp_phi)) (camupper_y * cos(vp_phi)) camupper_z], [(cos(vp_upperangle)*sin(-vp_phi)) (cos(vp_upperangle)*cos(vp_phi)) (sin(vp_upperangle))]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((((p(:,1,:)) - (camupper_y * sin(-vp_phi))) .^ 2) + ((p(:,2,:) - (camupper_y * cos(vp_phi))) .^ 2) + ((p(:,3,:) - camupper_z) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%27
surface_list(end+1).description = 'VP2upper - glycol-glass interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [((camupper_y * sin(-vp_phi)) + (vp_winthick*cos(vp_upperangle) *sin(-vp_phi))) ...
	((camupper_y * cos(vp_phi)) + (vp_winthick*cos(vp_upperangle) * cos(vp_phi))) ...
	(camupper_z + (vp_winthick*sin(vp_upperangle)))], ...
	[(cos(vp_upperangle)*sin(-vp_phi)) (cos(vp_upperangle)*cos(vp_phi)) (sin(vp_upperangle))]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((((p(:,1,:)) - ((camupper_y * sin(-vp_phi)) + (vp_winthick*cos(vp_upperangle) *sin(-vp_phi)))) .^ 2) + ...
	((p(:,2,:) - ((camupper_y * cos(vp_phi)) + (vp_winthick*cos(vp_upperangle) * cos(vp_phi)))) .^ 2) + ...
	((p(:,3,:) - (camupper_z + (vp_winthick*sin(vp_upperangle)))) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%28
surface_list(end+1).description = 'VP2lower - window casing';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [(camlower_y * sin(-vp_phi)) (camlower_y * cos(vp_phi)) camlower_z], ...
	[(cos(vp_lowerangle)*sin(-vp_phi)) (cos(vp_lowerangle)*cos(vp_phi)) (sin(vp_lowerangle))], vp_outerrad);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((cos(vp_lowerangle)*sin(-vp_phi)) .* (p(:,1,:))) + ((cos(vp_lowerangle)*cos(vp_phi)) .* (p(:,2,:))) + ((sin(vp_lowerangle)) .* (p(:,3,:))) >= ...
	((cos(vp_upperangle)*sin(-vp_phi)) .* (camlower_y * sin(-vp_phi))) + ((cos(vp_upperangle)*cos(vp_phi)) .* (camlower_y * cos(vp_phi))) + ((sin(vp_upperangle)) .* (camlower_z))) & ... 	
	(p(:,2,:) <= 0) & p(:,2,:) >= (camupper_y + (vp_outerrad .* sin(vp_upperangle))) & ...
    (((p(:,1,:) .^2) + (p(:,2,:) .^2)) >= (pv_outerrad .^ 2)) ), ...
     size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;

%29
surface_list(end+1).description = 'VP2lower - glass-air interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [(camlower_y * sin(-vp_phi)) (camlower_y * cos(vp_phi)) camlower_z], [(cos(vp_lowerangle)*sin(-vp_phi)) (cos(vp_lowerangle)*cos(vp_phi)) (sin(vp_lowerangle))]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((((p(:,1,:)) - (camlower_y * sin(-vp_phi))) .^ 2) + ((p(:,2,:) - (camlower_y * cos(vp_phi))) .^ 2) + ((p(:,3,:) - camlower_z) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%30
surface_list(end+1).description = 'VP2lower - glycol-glass interface';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [((camlower_y * sin(-vp_phi)) + (vp_winthick*cos(vp_lowerangle) *sin(-vp_phi))) ...
	((camlower_y * cos(vp_phi)) + (vp_winthick*cos(vp_lowerangle) * cos(vp_phi))) ...
	(camlower_z + (vp_winthick*sin(vp_lowerangle)))], ...
	[(cos(vp_lowerangle)*sin(-vp_phi)) (cos(vp_lowerangle)*cos(vp_phi)) (sin(vp_lowerangle))]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    ((((p(:,1,:)) - ((camlower_y * sin(-vp_phi)) + (vp_winthick*cos(vp_lowerangle) *sin(-vp_phi)))) .^ 2) + ...
	((p(:,2,:) - ((camlower_y * cos(vp_phi)) + (vp_winthick*cos(vp_lowerangle) * cos(vp_phi)))) .^ 2) + ...
	((p(:,3,:) - (camlower_z + (vp_winthick*sin(vp_lowerangle)))) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

%% Geometry that defines the camera surfaces    ((((((((((THIS CODE IS NOT RUNNING)))))))))))
%{
surface_list(end+1).description = 'Lower LED ring 1';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 (camlower_y + (vp_winthick*cos(vp_lowerangle))) (camlower_z + (vp_winthick*sin(vp_lowerangle)))], ...
    [0 1 tan(vp_lowerangle)]);
surface_list(end).inbounds_function = @(p)(reshape( ... 
    (p(:,1,:)) .^2 + (p(:,3,:)-camlower_z-(vp_winthick*sin(vp_lowerangle))) .^2 < (vp_outerrad)^2 & ...
    (p(:,1,:)) .^2 + (p(:,3,:)-camlower_z-(vp_winthick*sin(vp_lowerangle))) .^2 > 10^2, ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
%}
%{
surface_list(end+1).description = 'Camera 1 lower';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 camlower_y-1 camlower_z-1], [0 1 tan(vp_lowerangle)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)) .^ 2) + ((p(:,2,:) - (camlower_y-1)) .^ 2) + ((p(:,3,:) - (camlower_z-1)) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Camera 1 upper';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 camupper_y-1 camupper_z-1], [0 1 tan(vp_upperangle)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)) .^ 2) + ((p(:,2,:) - (camupper_y-1)) .^ 2) + ((p(:,3,:) - (camupper_z-1)) .^ 2) <= (vp_outerrad .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'Test rectangle';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 camupper_y-1+5 camupper_z-1+10], [0 1 tan(vp_upperangle)]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    p(:,1,:) <= 3 & p(:,1,:)>=-3), ...
    size(p,1), [] ));
surface_list(end).n_outside = n_pressurewindow;
surface_list(end).n_inside = n_air;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
%}
%}
%% Geometry that defines the retroreflector   ((((((SOME PART OF THIS CODE IS NOT RUNNING))))))
% Next we create the reflector which is just a lining of the pressure
% vessel with a cut out for where the cameras are, it is a distance of 3
% cm from the walls of the pressure vessel


retro_rad = pv_outerrad-0.01;
%{
surface_list(end+1).description = 'Reflector - just above copper ring';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0,0,(innervessel_startcyl+innervessel_cyllength+1)],[0,0,1]);
surface_list(end).inbounds_function = @(p)(reshape( (((p(:,1,:) .^2) + (p(:,2,:).^2)) ...
    <= ((retro_rad).^2) & ((p(:,1,:) .^2) + (p(:,2,:).^2)) >= ((outervessel_cylrad + 1).^2)),  ...
    size(p,1), []));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 0.5;
height_adj = -10;

%}
refl='hemi';
cone_slope=1.5;
dome_angle=0;
switch refl
    case 'cone'
        height_adj=0;
    case 'flat'
        height_adj=10;
    case 'hemi'
        height_adj=-10;
end






ref_angle=deg2rad(8);
%cyl_height=(pv_cyllength) - ((retro_rad) .* (tan(deg2rad(30))))+height_adj; %subtract 5 for dome roof, add 10 for flat roof
cyl_height =pv_cyllength+innervessel_startcyl;
surface_list(end+1).description = 'Reflector 3/4 cylinder wall';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], [0 0 1], (pv_outerrad-0.01));
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,3,:)>=innervessel_startcyl) & ...
    (p(:,3,:)<(cyl_height)) & ...
    (((p(:,1,:) < -pv_outerrad*sin(ref_angle) ) | (p(:,2,:) > -pv_outerrad*cos(vp_phi+ref_angle))) | p(:,3,:) < innervessel_cyllength+innervessel_startcyl+2.54) ), ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;



%%%%%%%%%%%%%%%%%%

ver1 = 50;
hor1 = 20;
bottom_slope = ver1/hor1;


lower_bottom = 71.5*2.54+innervessel_cyllength+innervessel_startcyl; %parte baja
bottom_height = lower_bottom + retro_rad*bottom_slope; %vertice , retro_rad=radio base del cono
%bottom_height = cyl_height + (retro_rad)*bottom_slope;
surface_list(end+1).description = 'Reflector - bottom cone section';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    [1,0,0;0,1,0;0,0,-bottom_slope^-2],[0,0,(2 *(bottom_slope^-2).* bottom_height)],-(bottom_height .^ 2) .* bottom_slope^-2);
%surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2 + p(:,2,:).^2 <= (retro_rad .^ 2)) & (p(:,3,:) < cyl_height + ver1) & ...
    %(((p(:,1,:) < -pv_outerrad*sin(ref_angle) ) | (p(:,2,:) > -pv_outerrad*cos(vp_phi+ref_angle))))),  ...
    %size(p,1), []));
surface_list(end).inbounds_function = @(p)(reshape( ((p(:,1,:).^2 + p(:,2,:).^2 <= (retro_rad .^ 2)) & (p(:,3,:) < lower_bottom + ver1) ...
    ), size(p,1), [])); %((p(:,1,:) < -pv_outerrad*sin(ref_angle) ) | (p(:,2,:) > -pv_outerrad*cos(vp_phi+ref_angle)))
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;



ver2 = 20;
hor2 = 20;
top_slope = ver2/hor2;
retro_rad1=retro_rad - 50
spacing = 0;
top_height = lower_bottom + ver1 + (retro_rad1 - hor1)*top_slope+spacing ;
plane_height = ver1 + ver2 + lower_bottom + spacing;
surface_list(end+1).description = 'Reflector - top cone section';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    [1,0,0;0,1,0;0,0,-top_slope^-2],[0,0,(2*top_slope^-2 .* top_height)],-(top_height .^ 2) .* top_slope^-2);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2+p(:,2,:).^2<= ((retro_rad1) .^ 2) & (p(:,3,:)<plane_height-spacing)),  ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;

%{

ver3 = 20;
hor3 = 35;
top_slope = ver3/hor3;
retro_rad1=retro_rad - 65
spacing = 0;
top_height = lower_bottom + ver1 + (retro_rad1 - hor1)*top_slope+spacing ;
plane_height = ver1 + ver3 + lower_bottom + spacing;
surface_list(end+1).description = 'Reflector - top cone section';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    [1,0,0;0,1,0;0,0,-top_slope^-2],[0,0,(2*top_slope^-2 .* top_height)],-(top_height .^ 2) .* top_slope^-2);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2+p(:,2,:).^2<= ((retro_rad1) .^ 2) & (p(:,3,:)<plane_height-spacing)),  ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;

%}


%{

ver2 = 15;
hor2 = 20;
top_slope = ver2/hor2;
spacing = 5;
top_height = 4*2.54;
top_height = lower_bottom + ver1 + (retro_rad - hor1)*top_slope+spacing;
plane_height = ver1 + ver2 + lower_bottom + spacing;
surface_list(end+1).description = 'Reflector - Connecting Part';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    [1,0,0;0,1,0;0,0,-top_slope^-2],[0,0,(2*top_slope^-2 .* top_height)],-(top_height .^ 2) .* top_slope^-2);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:).^2+p(:,2,:).^2<= ((retro_rad-0.1) .^ 2) & (p(:,3,:)<plane_height-spacing)),  ...
    size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
%}
%{
plane_rad = 20;
surface_list(end+1).description = 'Reflector - top plane';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0,0,plane_height],[0,0,1]);
surface_list(end).inbounds_function = @(p)(reshape( p(:,1,:) .^2 + ...
    p(:,2,:) .^2 <= (plane_rad) .^2, size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;
%}

ver3 = 25;
hor3 =20;
lower_slope = ver3/hor3;
lower_height = 0;
surface_list(end+1).description = 'Reflector - bottom cone section';
surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
    [1,0,0;0,1,0;0,0,-lower_slope^-2],[0,0,(2 *(lower_slope^-2).* lower_height)],(lower_height .^ 2) .* lower_slope^-2);
surface_list(end).inbounds_function = @(p)(reshape( p(:,3,:) > 35 & p(:,3,:)<45,size(p,1), []));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;

%specify where the retroreflector surface should end, given in radians from
%the centre of both viewports. The offset is chosen to be equal on both
%sides

%{
switch refl
    case 'cone'
        cone_height=(retro_rad)/sqrt(cone_slope)+cyl_height;
        surface_list(end+1).description = 'Reflector - top cone';
        surface_list(end).intersect_function = @(sp,indir)RayToQuadsurface(sp,indir, ...
            [1,0,0;0,1,0;0,0,-cone_slope],[0,0,(2*cone_slope .* cone_height)],-(cone_height .^ 2) .* cone_slope);
        surface_list(end).inbounds_function = @(p)(reshape( (((((p(:,1,:)).^2) + ((p(:,2,:)).^2)) <= ((retro_rad) .^ 2)) & (p(:,3,:)<cone_height)),  ...
            size(p,1), []));
        surface_list(end).n_outside = inf;
        surface_list(end).n_inside = n_hydraulic;
        surface_list(end).surface_type = 'retro';
        surface_list(end).absorption = 0;
    case 'flat'
        plane_height=cyl_height;
        surface_list(end+1).description = 'Reflector - top plane';
        surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
            [0,0,plane_height],[0,0,1]);
        surface_list(end).inbounds_function = @(p)(reshape( p(:,1,:) .^2 + ...
            p(:,2,:) .^2 <= (retro_rad) .^2, size(p,1), []));
        surface_list(end).n_outside = inf;
        surface_list(end).n_inside = n_hydraulic;
        surface_list(end).surface_type = 'retro';
        surface_list(end).absorption = 0;
    case 'hemi'
        %dome_height=cyl_height-(retro_rad)*tan(dome_angle);
        surface_list(end+1).description = 'Reflector - top hemisphere';
        surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
            [0,0,146.5],pv_outerrad-1);
        surface_list(end).inbounds_function = @(p)(reshape( ( ...
    p(:,3,:) > 94 & p(:,3,:) < 200 & ((p(:,1,:) < -pv_outerrad*sin(ref_angle) ) ...
    | (p(:,2,:) > -pv_outerrad*cos(vp_phi+ref_angle)))),size(p,1), []));
        %p(:,2,:) > -(retro_rad)*sqrt(5/9)
        %(acos((p(:,3,:)-dome_height) ./sqrt((p(:,1,:)).^2+(p(:,2,:)).^2+(p(:,3,:)-dome_height).^2))).^2 + (atan(p(:,2,:) ./ p(:,1,:))).^2 > (60*pi/180)^2
        surface_list(end).n_outside = inf;
        surface_list(end).n_inside = n_hydraulic;
        surface_list(end).surface_type = 'retro';
        surface_list(end).absorption = 0;
%}
%{
surface_list(end+1).description = 'Temp sensor plate 1';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 pv_outerrad - 5 camupper_z-1], [0 1 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)) .^ 2) + ((p(:,3,:) - (camupper_z-1)) .^ 2) <= (5 .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;
%}
%{
surface_list(end+1).description = 'Temp sensor plate 2';
surface_list(end).intersect_function = @(sp,indir)RayToPlane(sp,indir, ...
    [0 pv_outerrad - 5 camlower_z-1], [0 1 0]);
surface_list(end).inbounds_function = @(p)(reshape( ( ...
    (((p(:,1,:)) .^ 2) + ((p(:,3,:) - (camlower_z-1)) .^ 2) <= (5 .^ 2))), ...
    size(p,1), [] ));
surface_list(end).n_outside = inf;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'retro';
surface_list(end).absorption = 1;
%}

%{
surface_list(end+1).description = 'Source Tube';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [12*2.54*sind(30) -12*2.54*cosd(30) 0], [0 0 1], 2.54);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:) > 120), ...
    size(p,1), []));
surface_list(end).n_outside = 1;
surface_list(end).n_inside = n_hydraulic;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
%}

%{
t_fm = 0.3;
l1_fm = 1;
l2_fm = 0.7;
r_fm = outervessel_cylrad+0.1;
surface_list(end+1).description = 'Fiducial Marks';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, [0 0 0], [0 0 1], r_fm);
surface_list(end).inbounds_function = @(p)(reshape( (p(:,3,:)>=outervessel_startcyl)...
    &( (p(:,3,:)<=(outervessel_startcyl+outervessel_cyllength))...
    & (p(:,3,:)>=outervessel_startcyl+50-t_fm/2 & p(:,3,:)<=outervessel_startcyl+50+t_fm/2)...
    &( ((p(:,1,:)>=-r_fm*sin(pi/6+l1_fm/(2*r_fm)) & p(:,1,:)<=-r_fm*sin(pi/6-l1_fm/(2*r_fm)))...
    & (p(:,2,:)>=0))...
    | ((p(:,1,:)<=-r_fm*sin(5*pi/6+l1_fm/(2*r_fm)) & p(:,1,:)>=-r_fm*sin(5*pi/6-l1_fm/(2*r_fm)))...
    & (p(:,2,:)<=0))...
    | ((p(:,2,:)<=r_fm*sin(l1_fm/(2*r_fm)) & p(:,2,:)>=-r_fm*sin(l1_fm/(2*r_fm)))...
    & (p(:,1,:)>=0)) )...
    )...
    |( (p(:,3,:)<=(outervessel_startcyl+outervessel_cyllength))...
    & (p(:,3,:)>=outervessel_startcyl+50-t_fm/2-l2_fm & p(:,3,:)<=outervessel_startcyl+50-t_fm/2)...
    &( ((p(:,1,:)>=-r_fm*sin(pi/6+t_fm/(2*r_fm)) & p(:,1,:)<=-r_fm*sin(pi/6-t_fm/(2*r_fm)))...
    & (p(:,2,:)>=0))...
    | ((p(:,1,:)<=-r_fm*sin(5*pi/6+t_fm/(2*r_fm)) & p(:,1,:)>=-r_fm*sin(5*pi/6-t_fm/(2*r_fm)))...
    & (p(:,2,:)<=0))...
    | ((p(:,2,:)<=r_fm*sin(t_fm/(2*r_fm)) & p(:,2,:)>=-r_fm*sin(t_fm/(2*r_fm)))...
    & (p(:,1,:)>=0)) )...
    )...
    , size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 1;
%}
end


