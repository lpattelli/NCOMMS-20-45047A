function plot_intensity(ax, pArr, rArr, I, fldPoints, edgecolor)
%PLOT_INTENSITY Utility function to plot an intensity pattern generated
%with CELES.
%
% This function can be used for convenience instead of the plot_intensity
% function bundled with CELES, which requires a Simulation object as input
%
% Author: Lorenzo Pattelli

hold(ax,'on')

if all(fldPoints(:,:,1) == fldPoints(1,1,1))    % fldPoints are on the yz plane
    perpdim = 1;                                % 1->x is the perp. direction
    draw_image(ax, fldPoints, I, perpdim, pArr, rArr, edgecolor)
    xlabel('y')
    ylabel('z')
    titlestring = sprintf('x = %g nm',fldPoints(1,1,1));
elseif all(fldPoints(:,:,2) == fldPoints(1,1,2))% fldPoints are on the xz plane
    perpdim = 2;                                % 2->y is the perp. direction
    draw_image(ax, fldPoints, I, perpdim, pArr, rArr, edgecolor)
    xlabel('x')
    ylabel('z')
    titlestring = sprintf('y = %g nm',fldPoints(1,1,2));
elseif all(fldPoints(:,:,3) == fldPoints(1,1,3))% fldPoints are on the xy plane
    perpdim = 3;                                % 3->z is the perp. direction
    draw_image(ax, fldPoints, I, perpdim, pArr, rArr, edgecolor)
    xlabel('x')
    ylabel('y')
    titlestring = sprintf('z = %g nm',fldPoints(1,1,3));
else
    error('fieldPoint must define an xy, xz or yz-like plane')
end

ax.DataAspectRatio = [1,1,1];
title(['Intensity map at ', titlestring])
hold(ax,'off')
end

function draw_image(ax, fldP, fld, perpdim, pArr, rArr, edgecolor)
xy = setdiff([1,2,3], perpdim);                 % here xy are the in-plane dimensions
x = fldP(:,:,xy(1));
y = fldP(:,:,xy(2));
imagesc(x(1,:), y(:,1), fld)                    % plot field on a xy plane
dist = abs(pArr(:,perpdim) - fldP(1,1,perpdim));% particle distances from xy plane
idx = find(dist<rArr);                          % find particles intersecting the plane
rArr(idx) = sqrt(rArr(idx).^2 - dist(idx).^2);  % overwrite radius of the intersection circle
for i=1:length(idx)
    rectangle(ax, ...
             'Position', [pArr(idx(i),xy)-rArr(idx(i)), [2,2]*rArr(idx(i))], ...
             'Curvature', [1 1], ...
             'FaceColor', 'none', ...
             'EdgeColor', edgecolor*[1,1,1], ...
             'LineWidth', 0.2)
end
axis([min(x(:)),max(x(:)),min(y(:)),max(y(:))]) % set axis tight to fldPoints
end
