% Generate E field data in the gap between two particles for multiple
% randomly oriented plane-wave illumination conditions and check the
% resulting k-vector orientation distribution
%
% Author: Lorenzo Pattelli

clear, close all

addpath(genpath(' path / to / CELES '))

%% define two-particle system and simulate several illumination conditions

rho = 1265; % [um^-3] number density of particles
a = 1000*rho^(-1/3); % [nm] typical average separation
prad = [40; 44]; % pick two plausible radii
gapwidth = a - sum(prad);
pos = [0, 0, -(gapwidth/2 + prad(1)); 0, 0, gapwidth/2 + prad(2)];

Nrealiz = 150000;

% initialize arrays to store computed field values
Esca_pool = zeros(Nrealiz,3);
Etot_pool = zeros(Nrealiz,3);

% define plotting point
x = 0; y = 0; z = 0;

% configure parameters
wl = 532;
ns = 1.963;
nh = 1.334;
lmax = 3;
refridx = ones(size(pos,1),1)*ns;

for iii=1:Nrealiz
    % generate rotation matrix
    theta = unifrnd(0,2*pi);
    phi = unifrnd(0,2*pi);
    u = unifrnd(0,1);
    V = [sqrt(u)*cos(phi),sqrt(u)*sin(phi),sqrt(1-u)].';
    RM = (2*(V * transpose(V)) - eye(3))*[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];

    % rotate particles
    rotpos = pos*RM;

    % configure simulation
    prtcls = celes_particles('positionArray', rotpos, 'refractiveIndexArray', refridx, 'radiusArray', prad);
    ifld = celes_initialField('polarAngle', 0, 'azimuthalAngle', 0, 'polarization', 'TE', 'beamWidth', Inf, 'focalPoint', [0,0,0]); % plane wave
    input = celes_input('wavelength', wl, 'mediumRefractiveIndex', nh, 'particles', prtcls, 'initialField', ifld);
    precnd = celes_preconditioner('type', 'none'); % no preconditioning needed in this case (only 2 particles!)
    solver = celes_solver('type', 'GMRES', 'tolerance', 1e-4, 'maxIter', 100, 'restart', 10, 'preconditioner', precnd);
    numrcs = celes_numerics('lmax', lmax, 'polarAnglesArray', 0:pi/5e4:pi, 'azimuthalAnglesArray', 0:pi/1e3:2*pi, 'gpuFlag', true, 'particleDistanceResolution', 0.01, 'solver', solver);
    output = celes_output('fieldPoints', [x(:),y(:),z(:)], 'fieldPointsArrayDims', size(x));

    simul = celes_simulation('input', input, 'numerics', numrcs, 'output', output);

    % run simulation
    simul.run;

    % evaluate field at output.fieldPoints
    simul.evaluateFields;

    Esca = gather(simul.output.scatteredField)/RM;  Esca_pool(iii,:) = reshape(Esca, [numel(x),3]); % counter-rotate field
    Etot = gather(simul.output.totalField)/RM;      Etot_pool(iii,:) = reshape(Etot, [numel(x),3]); % counter-rotate field
end

%% recombine all simulations applying random phases

Nrealiz_deph = 50000;
Epool = zeros(Nrealiz_deph,3);

for rr = 1:Nrealiz_deph
    sf = rand([Nrealiz,1]);
    dephasing = unifrnd(0,2*pi,[Nrealiz,1]);
    Epool(rr,:) = sum(sf.*Etot_pool.*exp(1i*dephasing));
end

%% determine the propagation direction of the resulting field

Qpr_deph = zeros(Nrealiz_deph, 3);

for ii=1:length(Epool)
    E = Epool(ii,:);
    [~, ~, ~, Qs] = field_pol(E, false); % set plot flag = false
    Qpr_deph(ii,:) = Qs(:,3)/norm(E);
end

%% plot final distribution

% in-plane directions are symmetric, z is not because
% the two particles have different sizes
Qpr_deph_flip = [abs(Qpr_deph(:,1:2)), Qpr_deph(:,3)];

figure('Renderer', 'painters', 'Position', [100 100 900 150])
histogram(Qpr_deph_flip(:,3),-1:0.02:1,'normalization','pdf'),
axis([-1 1 0 2])
figure('Position', [100 100 400 900])
scatter3(Qpr_deph_flip(:,1), Qpr_deph_flip(:,2), Qpr_deph_flip(:,3), 1, '.'), axis equal
axis([0 1 0 1 -1 1]), view([70 10])
