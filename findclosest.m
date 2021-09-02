% findclosest.m
% This file takes a single value ("ritzval") and finds the index
% ("closest") and distance ("dist") to the closest item in list "ee". 
%
% Last edited by: Erin Carson, 2021

function [closest, dist] = findclosest(ritzval, ee)

closest = 1;
dist = norm(ritzval - ee(1));
for i = 2:numel(ee)
    disttmp = norm(ritzval - ee(i));
    if disttmp < dist
        closest = i;
        dist = disttmp;
    end
end


end