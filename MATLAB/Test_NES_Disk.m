clear; close all; clc;

diskTri = diskTriangulation(30);
F = diskTri.ConnectivityList;
x = diskTri.Points;
zeroID = knnsearch(x, [0 0]);

initV = [x, -1e-3 .* (exp(-sum(x.^2, 2))-1)];
tarV = [x, -exp(-sum(x.^2, 2))+1];

[tarL, tarAng] = calculateEdgeLengthsAndAngles(F, tarV);

h = 0.1;
nu = 1/3;

movV = minimizeElasticEnergy( F, initV, tarL, ...
    'TargetAngles', tarAng, ...
    'Thickness', h, 'Poisson', nu, ...
    'MaxIterations', 10000, 'iterDisplay', true, ...
    'LineSearchMethod', 'Backtracking', ...
    'TargetVertices', zeroID, 'TargetLocations', initV(zeroID, :), ...
    'Alpha', 10 );

[movL, movAng] = calculateEdgeLengthsAndAngles(F, movV);

lengthErr = abs(movL-tarL) ./ abs(tarL);
angleErr = abs(movAng - tarAng) ./ abs(tarAng);

fprintf('Max edge length error = %0.5e\n', max(lengthErr));
fprintf('Median edge length error = %0.5e\n', median(lengthErr));

fprintf('Max bending angle error = %0.5e\n', max(angleErr));
fprintf('Median bending angle error = %0.5e\n', median(angleErr));


subplot(1,3,1)
trisurf(triangulation(F, initV));
hold on
scatter3(initV(zeroID, 1), initV(zeroID, 2), initV(zeroID, 3), 'filled', 'r');
hold off
axis equal tight

subplot(1,3,2)
trisurf(triangulation(F, movV));
hold on
scatter3(movV(zeroID, 1), movV(zeroID, 2), movV(zeroID, 3), 'filled', 'r');
hold off
axis equal tight

subplot(1,3,3)
trisurf(triangulation(F,tarV));
hold on
scatter3(tarV(zeroID, 1), tarV(zeroID, 2), tarV(zeroID, 3), 'filled', 'r');
hold off
axis equal tight

