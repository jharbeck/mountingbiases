function ecef2geodetic,x,y,z
;function to convert ECEF to lat, lon, and altitude


;[phi, lambda, h] = ecef2geodetic(p(:,1), p(:,2), p(:,3), [R sqrt( 1 - ( 1 - f )^2 )]);



;************Constants*************
f = 1/298.257223563D  ;flattening term from WGS84 definition
er = 6378137.0D ;Earth equatorial radius in meters
el = 0.0818191908426215D ;8.1819190842622e-2  ;WGS84 constant
;**********************************

;Get pixel position in lat, lon, altitude, old way
;x=pixel_ecef[0]
;y=pixel_ecef[1]
;z=pixel_ecef[2]
;bcam = sqrt(er^2*(1-el^2))
;epcam = sqrt((er^2-bcam^2)/bcam^2)
;pcam = sqrt(x^2+y^2)
;thcam = atan(er*z,bcam*pcam)
;
;img_lat[ii,jj] = atan( (z+epcam^2*bcam*sin(thcam)^3), (pcam-el^2*er*cos(thcam)^3) )
;img_lon[ii,jj] = atan(y,x)










; Ellipsoid constants
;a  = er ;ellipsoid(1);       ; Semimajor axis
;e2 = ellipsoid(2) ^ 2;   ; Square of first eccentricity
a = er
e2 = el ^2.0D
ep2 = e2 / (1 - e2);     ; Square of second eccentricity
f = 1 - sqrt(1 - e2);    ; Flattening
b = a * (1 - f);         ; Semiminor axis

; Longitude
lambda = atan(y,x);

; Distance from Z-axis
rho = sqrt( ABS(x)^2.0D + ABS(y)^2.0D) ;hypot(x,y);

; Bowring's formula for initial parametric (beta) and geodetic (phi) latitudes
beta = atan(z, (1 - f) * rho);
phi = atan(z   + b * ep2 * sin(beta)^3.0, rho - a * e2  * cos(beta)^3.0);

; Fixed-point iteration with Bowring's formula
; (typically converges within two or three iterations)
betaNew = atan((1 - f)*sin(phi), cos(phi));
count = 0;
while beta NE betaNew AND count < 5 DO BEGIN
  beta = betaNew;
  phi = atan(z   + b * ep2 * sin(beta)^3,  rho - a * e2  * cos(beta)^3.0);
  betaNew = atan((1 - f)*sin(phi), cos(phi));
  count = count + 1;
endwhile

; Calculate ellipsoidal height from the final value for latitude
sinphi = sin(phi);
N = a / sqrt(1 - e2 * sinphi^2.0);
h = rho * cos(phi) + (z + e2 * N * sinphi) * sinphi - N;


llh = [phi, lambda, h]  ;lat and lon in radians, height in meters
return,llh
END