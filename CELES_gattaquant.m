% Generate E field data inside a scattering medium for multiple randomly
% oriented plane-wave illumination conditions
%
% Author: Lorenzo Pattelli

clear, close all

addpath(genpath(' path / to / CELES '))

%% load target structure (scatterers + gattaquant nanorulers)
packingfile = 'spherical_target_120k.mat';

load(packingfile)

%% set up and run CELES simulations

% number of independent, randomly oriented plane-waves to simulate
repetitions = 250;

for iii=1:repetitions
    % generate random rotation matrix
    theta = unifrnd(0,2*pi);
    phi = unifrnd(0,2*pi);
    u = unifrnd(0,1);
    V = [sqrt(u)*cos(phi),sqrt(u)*sin(phi),sqrt(1-u)].';
    RM = (2*(V * transpose(V)) - eye(3))*[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];

    % rotate particles and rods
    rotpos = pos*RM;
    rotrod(:,1:3) = rodpos(:,1:3)*RM;
    rotrod(:,4:6) = rodpos(:,4:6)*RM;

    % define plotting points along the nanorulers
    resolution = 5; % [nm]
    npoints = Lrods/resolution +1;
    x = zeros(Nrods, npoints);
    y = zeros(Nrods, npoints);
    z = zeros(Nrods, npoints);
    for i=1:Nrods
        x(i,:) = linspace(rotrod(i,1),rotrod(i,4),npoints).';
        y(i,:) = linspace(rotrod(i,2),rotrod(i,5),npoints).';
        z(i,:) = linspace(rotrod(i,3),rotrod(i,6),npoints).';
    end

    % simulation parameters
    wl = 532;       % wavelength
    ns = 1.963;     % ZnO
    nh = 1.334;     % water-agarose gel
    lmax = 3;
    refridx = ones(length(pos),1)*ns;

    % configure simulation
    prtcls = celes_particles('positionArray', rotpos, 'refractiveIndexArray', refridx, 'radiusArray', prad);
    ifld = celes_initialField('polarAngle', 0, 'azimuthalAngle', 0, 'polarization', 'TE', 'beamWidth', Inf, 'focalPoint', [0,0,0]);
    input = celes_input('wavelength', wl, 'mediumRefractiveIndex', nh, 'particles', prtcls, 'initialField', ifld);
    precnd = celes_preconditioner('type', 'blockdiagonal', 'partitionEdgeSizes', [600,600,600]);
    solver = celes_solver('type', 'GMRES', 'tolerance', 1e-3, 'maxIter', 2000, 'restart', 200, 'preconditioner', precnd);
    numrcs = celes_numerics('lmax', lmax, 'polarAnglesArray', 0:pi/5e3:pi, 'azimuthalAnglesArray', 0:pi/1e2:2*pi, 'gpuFlag', true, 'particleDistanceResolution', 0.02, 'solver', solver);
    output = celes_output('fieldPoints', [x(:),y(:),z(:)], 'fieldPointsArrayDims', size(x));

    simul = celes_simulation('input', input, 'numerics', numrcs, 'output', output);

    % run simulation
    simul.run;

    % evaluate field at output.fieldPoints
    simul.evaluateFields;

    Erod = gather(simul.output.totalField)/RM; % counter-rotate field
    Erod = reshape(Erod, [size(x),3]);
    % Irod = reshape(sum(abs(Erod).^2,2), size(x));

    save([packingfile(1:end-4), '_Erods', int2str(Lrods), '_', int2str(iii), '.mat'], 'Erod', 'RM', 'theta', 'phi', 'u', 'rodpos', 'resolution')

    % let's plot & save a slice, too
    xarray = -2000:4:2000; zarray = xarray;
    [x,z] = meshgrid(xarray, zarray); y = zeros(size(x));
    fprintf('evaluating slice (corresponding to the pre-rotated y=0 plane)\n')
    simul.output.fieldPoints = [x(:),y(:),z(:)]*RM;
    simul.output.fieldPointsArrayDims = size(x);
    simul.evaluateFields;

    Eslice = gather(simul.output.totalField)/RM; % counter-rotate field
    Eslice = reshape(Eslice,[size(x),3]);
    Ipts = [x(:),y(:),z(:)]; % save directly the original points rather than counter-rotating them
    Iptsr = reshape(Ipts,[size(x),3]);
    intpts = simul.output.internalIndices;
    save([packingfile(1:end-4), '_Eslice_', int2str(iii), '.mat'], 'Eslice', 'RM', 'theta', 'phi', 'u', 'Iptsr', 'intpts')

    clear xarray zarray x y z Eslice Ipts Iptsr intpts
end
