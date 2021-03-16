% function to calcluate the midpoint between two given lat/lon
% assumes input is in degrees

function [Lat,Lon] = LatLonMidpoint(lat1,lon1,lat2,lon2)

X = cosd(lat2) * cosd(lon2 - lon1);
Y = cosd(lat2) * sind(lon2 - lon1);

Lat = atan2d(sind(lat1) + sind(lat2), sqrt( (cosd(lat1) + X)^2 + Y^2 ));
Lon = lon1 + atan2d(Y, cosd(lat1) + X);

end