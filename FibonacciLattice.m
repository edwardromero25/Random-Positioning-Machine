classdef FibonacciLattice
    properties
        ID
        x
        y
        z
        pathCoords
        num_points = 1000
    end

    methods
        function obj = FibonacciLattice(ID, x, y, z)
            obj.ID = ID;
            obj.x = x;
            obj.y = y;
            obj.z = z;
            obj.pathCoords = [x(:), y(:), z(:)];
        end

        function score = getDistribution(obj)
            [Xs, Ys, Zs] = obj.createSphere();
            sphereCoords = [Xs(:), Ys(:), Zs(:)];
            score = obj.getDistributionNum(sphereCoords);
        end
    end

    methods (Access = private)
        function [Xs, Ys, Zs] = createSphere(obj)
            golden_r = (sqrt(5.0) + 1.0) / 2.0;
            golden_a = (2.0 - golden_r) * (2.0 * pi);

            Xs = zeros(1, obj.num_points);
            Ys = zeros(1, obj.num_points);
            Zs = zeros(1, obj.num_points);

            for i = 0:obj.num_points-1
                ys = 1 - (i / (obj.num_points - 1)) * 2;
                radius = sqrt(1 - ys^2);
                theta = golden_a * i;

                xs = cos(theta) * radius;
                zs = sin(theta) * radius;

                Xs(i+1) = xs;
                Ys(i+1) = ys;
                Zs(i+1) = zs;
            end
        end

        function octants = splitSphere(~, sphereCoords)
            octants = struct('posI', [], 'posII', [], 'posIII', [], 'posIV', [], ...
                             'negI', [], 'negII', [], 'negIII', [], 'negIV', []);

            for i = 1:size(sphereCoords, 1)
                row = sphereCoords(i, :);
                if row(3) > 0
                    if row(2) > 0
                        if row(1) > 0
                            octants.posI = [octants.posI; row];
                        else
                            octants.posII = [octants.posII; row];
                        end
                    elseif row(1) > 0
                        octants.posIV = [octants.posIV; row];
                    else
                        octants.posIII = [octants.posIII; row];
                    end
                else
                    if row(2) > 0
                        if row(1) > 0
                            octants.negI = [octants.negI; row];
                        else
                            octants.negII = [octants.negII; row];
                        end
                    elseif row(1) > 0
                        octants.negIV = [octants.negIV; row];
                    else
                        octants.negIII = [octants.negIII; row];
                    end
                end
            end
        end

        function oct = getPathOctant(~, pathRow)
            if pathRow(3) > 0
                if pathRow(2) > 0
                    if pathRow(1) > 0
                        oct = 'posI';
                    else
                        oct = 'posII';
                    end
                elseif pathRow(1) > 0
                    oct = 'posIV';
                else
                    oct = 'posIII';
                end
            else
                if pathRow(2) > 0
                    if pathRow(1) > 0
                        oct = 'negI';
                    else
                        oct = 'negII';
                    end
                elseif pathRow(1) > 0
                    oct = 'negIV';
                else
                    oct = 'negIII';
                end
            end
        end

        function dist = getDistanceBetween(~, pathTupleCoords, sphereTupleCoords)
            diffX = pathTupleCoords(1) - sphereTupleCoords(1);
            diffY = pathTupleCoords(2) - sphereTupleCoords(2);
            diffZ = pathTupleCoords(3) - sphereTupleCoords(3);
            dist = sqrt(diffX^2 + diffY^2 + diffZ^2);
        end

        function score = getDistributionNum(obj, sphereCoords)
            octants = obj.splitSphere(sphereCoords);
            pathMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
            repeatTime = -1;

            for i = 1:size(obj.pathCoords, 1)
                pathRow = obj.pathCoords(i, :);
                pathOctant = obj.getPathOctant(pathRow);
                sphereCoordsSplit = octants.(pathOctant);
                distDict = zeros(size(sphereCoordsSplit, 1), 1);
                repeatTime = repeatTime + 1;

                for j = 1:size(sphereCoordsSplit, 1)
                    distDict(j) = obj.getDistanceBetween(pathRow, sphereCoordsSplit(j, :));
                end

                [~, idx] = sort(distDict);
                top3 = sphereCoordsSplit(idx(1:3), :);
                
                segmentKey = sprintf('%.8f,', top3');
                segmentKey(end) = []; 

                if isKey(pathMap, segmentKey)
                    pathMap(segmentKey) = [pathMap(segmentKey), repeatTime];
                else
                    pathMap(segmentKey) = [repeatTime];  
                end
            end

            score = length(pathMap.keys);
        end
    end
end
