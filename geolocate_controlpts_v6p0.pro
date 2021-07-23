function geolocate_ControlPts_v6p0,img_cols,img_rows,ctrl_pixels,ctrl_pixels_elev,lat_in,lon_in,alt_in,pitch_in,roll_in,heading_in,focal_length,sensor_x,sensor_y,camera_offset,pitch_bias=pitch_bias,roll_bias=roll_bias,heading_bias=heading_bias,xcorr=xcorr,ycorr=ycorr


  ;This program gets the latitude and longitude for a set of control points (instead of a whole image array)

  
  ;ctrl_pixels is the pixel (i,j) location of the control points in the image
  ;img_cols and img_rows are the total number of columns and rows in the image

  ;Version 1p0, 02/04/2016 first version of the code
  ;Version 2p0, 02/23/2016 added in lens distortion correction capability
  ;Version 3p0, 05/22/2016 fixed error with camera pointing vector not extending to the surface for non-zero pitch and roll, mounting biases now
  ;                         taken into account properly through a coordinate transformation
  ;Version 3p2, 02/17/2017 fixed error in total path length and angle calculation
  ;Version 4p0, 10/17/2017 fixed error in aircraft pointing and mounting bias matrix multiplication order, added in function to calculate
  ;                         ECEF to lat/lon/altitude conversion which fixes a previous error in altitude determination
  ;Version 4p1, 10/25/2017 fixed error in determination of the angular component of each image pixel, previous version assumed linear relationship
  ;                         between nadir and FOV, new version assumes image pixels have a uniform size for a flat undistorted image
  ;Version 5p0, 07/17/2018 Switched from flat plane trigonometry to a projection camera model, needed to make the code applicable for larger pitch and roll values
  ;Version 5p3 07/20/2018 Switched from flat plane trigonometry to a projection camera model, needed to make the code applicable for larger pitch and roll values, fixed many errors from 5p0
  ;version 6p0, 2/25/2020 Updated to include elevation input for input ctrl_pixels

  ;****************Input*************
  gps_lat = lat_in  ;geodetic latitude of gps antenna in degrees
  gps_lon = lon_in  ;geodetic longitude of gps antenna in degrees
  gps_alt = alt_in  ;altitude (above ground) of gps antenna in meters
  p = pitch_in  ;pitch of airplane in degrees
  r = roll_in ;roll of airplane in degrees
  h = heading_in  ;heading of airplane in degrees
  ;camera_offset = [-3.89, 0.0, 4.05]   ;offset for C-130 from Jim and Craig, 3-d distance vector between gps antenna and camera boresight (in NED coordinate system), positive x-axis in the direction of North. column 1 is north (positive value is north), co
  ;camera_offset = [15, 0.0, 3.0] ;offset for NOAA43 P-3 from Matt and Dana
  ;camera_offset = [-4.148D,0.004D,3.104D] ;offset for P-3B for 2017, from 20170311 ancillary file

  ;fov_y = 22.0D*2 ; ;45.627 originally   ;Full angle field of view (degrees)
  ;fov_x = 33.0D*2 ; 66.0 original    ;field of view (degrees) ;old fov: 23.962*2
  ;focal_length = 28e-3  ;lens focal length
  ;**********************************

  ;sensor_x = 2.0*focal_length*tan(fov_x/2.0*!pi/180)  ;physical size of the camera sensor derived from the effective field of view and focal length, can be used if field of view is known but sensor size is not
  ;sensor_y = 2.0*focal_length*tan(fov_y/2.0*!pi/180)

  ;************Keywords**************
  if n_elements(pitch_bias) eq 0 then pitch_bias = 0.0
  if n_elements(roll_bias) eq 0 then roll_bias = 0.0
  if n_elements(heading_bias) eq 0 then heading_bias = 0.0
  ;**********************************

  ;************Constants*************
  f = 1/298.257223563D  ;flattening term from WGS84 definition
  er = 6378137.0D ;Earth equatorial radius in meters
  el = 0.0818191908426215D  ;WGS84 constant
  ;**********************************


  ;***Convert all variables to radians***
  gps_lat = gps_lat*!pi/180.0D
  gps_lon = gps_lon*!pi/180.0D
  p = 1.0*p*!pi/180.0D
  r = 1.0*r*!pi/180.0D
  h = h*!pi/180.0D
  pitch_bias = 1.0*pitch_bias*!pi/180.0D
  roll_bias = 1.0*roll_bias*!pi/180.0D
  heading_bias = heading_bias*!pi/180.0D
  ;**************************************

  ;***********Convert gps lat/lon to Earth-Centered Earth-Fixed (ECEF) coordinates*************
  glat = atan( (1.0-f)^2.0 * tan(gps_lat)) ;geocentric latitude at mean sea-level
  rs = sqrt( er^2 / (1+(1.0/(1-f)^2 - 1.0)*(sin(glat))^2)  )  ;radius at surface point
  p_ecef = dblarr(3)
  p_ecef[0] = rs*cos(glat)*cos(gps_lon) + gps_alt*cos(gps_lat)*cos(gps_lon)  ;x coordinate
  p_ecef[1] = rs*cos(glat)*sin(gps_lon) + gps_alt*cos(gps_lat)*sin(gps_lon)  ;y coordinate
  p_ecef[2] = rs*sin(glat) + gps_alt*sin(gps_lat) ;z coordinate
  ;********************************************************************************************


  ;**********Get coordinate rotation matrix to align camera coordinate system with aircraft coordinate system with the x-axis in the direction of North*************
  cp = cos(p)
  sp = sin(p)
  cr = cos(r)
  sr = sin(r)
  ch = cos(h)
  sh = sin(h)

  T = [ [ch*cp, ch*sp*sr-sh*cr, ch*sp*cr+sh*sr],$ ;The T(heading,pitch,roll) rotation aligns the camera’s coordinate system with the aircraft coordinate system, with the x-axis in the direction of the north
    [sh*cp, sh*sp*sr+ch*cr, sh*sp*cr-ch*sr],$
    [-1.0*sp, cp*sr, cp*cr]]

  st = sin(gps_lat)
  ct = cos(gps_lat)
  sl = sin(gps_lon)
  cl = cos(gps_lon)


  NED_R = [ [-1.0*st*cl, -1.0*sl, -1.0*ct*cl],$ ;R(theta,lambda) rotation transforms the camera location direction vector to ECEF geocentric coordinates
    [-1.0*st*sl, cl, -1.0*ct*sl],$
    [ct, 0, -1.0*st] ]

  camera_ecef = NED_R ## T ## camera_offset + p_ecef  ;camera position in ECEF

  ;******************************************************************************************************************************************************************


  ;*****************Get camera position in lat, lon, altitude*********************
  llh = ecef2geodetic(camera_ecef[0],camera_ecef[1],camera_ecef[2])
  camera_lat = llh[0]
  camera_lon = llh[1]
  camera_alt = llh[2]


  ;print,'Camera lat, lon, alt: ', camera_lat*180.0/!pi,camera_lon*180.0/!pi,camera_alt
  ;print,'Distance between gps and camera points (meters): ',distance(gps_lat*180.0/!pi,gps_lon*180.0/!pi,camera_lat*180.0/!pi,camera_lon*180.0/!pi)*1000.0

  ;**********************************************************************************



  center_x = float(img_cols/2)
  center_y = float(img_rows/2)




  camera_r_vec = dblarr(3)  ;camera pointing vector in NED coordinates, column 1 is north, column 2 is east, column 3 is down



  st = sin(camera_lat)
  ct = cos(camera_lat)
  sl = sin(camera_lon)
  cl = cos(camera_lon)


  NED_R_camera =  [ [-1.0*st*cl, -1.0*sl, -1.0*ct*cl,0],$ ;rotation matrix to transform from North-East-Down coordinate system (NED) to ECEF
    [-1.0*st*sl, cl, -1.0*ct*sl,0],$
    [ct, 0, -1.0*st,0],$
    [0,0,0,1] ]



  ;**********Get coordinate rotation matrix to align camera coordinate system with aircraft coordinate system with the x-axis in the direction of North for camera with mounting bias*************
 
  cp = cos(-1.0*p)	;negative signs needed to compensate for coordinate system change in image plane
  sp = sin(-1.0*p)
  cr = cos(-1.0*r)
  sr = sin(-1.0*r)
  ch = cos(h)
  sh = sin(h)   


  T = [ [ch*cp, cp*sh, -1.0*sp,0],$ ;The T(heading,pitch,roll) rotation aligns the camera coordinate system with the aircraft coordinate system, with the x-axis in the direction of the north, from Barber and Redding, 2006 eqn 3 (transpose of Euler rotation matrix)
      [sr*sp*ch - cr*sh, sh*sp*sr+ch*cr, sr*cp,0],$
      [cr*sp*ch + sr*sh, cr*sp*sh-sr*ch, cr*cp,0],$
      [0,0,0,1] ]



     
     
  cp = cos(-1.0*pitch_bias)	;negative signs needed to compensate for coordinate system change in image plane
  sp = sin(-1.0*pitch_bias)
  cr = cos(-1.0*roll_bias)
  sr = sin(-1.0*roll_bias)
  ch = cos(heading_bias)
  sh = sin(heading_bias)
     



  T_mounting = [ [ch*cp, cp*sh, -1.0*sp,0],$ ;Aligns the camera coordinate system with the aircraft coordinate system, accounts for mounting bias of the instrument
    [sr*sp*ch - cr*sh, sh*sp*sr+ch*cr, sr*cp,0],$
    [cr*sp*ch + sr*sh, cr*sp*sh-sr*ch, cr*cp,0],$
    [0,0,0,1] ]




  ;********************Geolocate individual camera pixels**************************

;  pixel_ecef = image_data*0.0D

  sz_ctrl = size(ctrl_pixels,/dimensions)
  tot_points = sz_ctrl(0)
  

  
  img_lat = dblarr(tot_points)
  img_lon = dblarr(tot_points)
  img_elev = dblarr(tot_points)



  FOR iip=0L,tot_points-1L DO BEGIN

      ii = ctrl_pixels(iip,0)
      jj = ctrl_pixels(iip,1)


      ;camera pointing vector is in NED coordinates, need to transform to ECEF then add to camera position to get pixel location

      ;**********This part depends on camera orientation, assumes camera is oriented such that if the airplane is traveling due north the image (0,0) index is the northwest corner

      IF n_elements(xcorr) EQ 0 THEN BEGIN  ;no lens distortion
        ;get pixel distance from the center point
        
        yp_dist = jj - center_y*2.0 ;center_y - jj ;y distance of pixel
        xp_dist = ii -  center_x*2.0 ;center_x - ii ;x distance of pixel

      ENDIF ELSE BEGIN  ;else correct for lens distortion
        ;get pixel distance from the center point with lens distortion, xcorr and ycorr are the apparent pixel locations due to distortion

        yp_dist = ycorr(ii,jj) - center_y*2.0 ;y distance of pixel
        xp_dist = xcorr(ii,jj) - center_x*2.0 ;x distance of pixel

      ENDELSE



      ;****************************************************
      
      sx = 1.0*sensor_x/img_cols/2.0 ;conversion factor from pixels to meters with units of meters/pixel, see eqn. 6 of Barber and Redding
      sy = 1.0*sensor_y/img_rows/2.0 ;conversion factor from pixels to meters, I think factor of 2 needs to be in since it is from the middle of the sensor
      
      cx = xp_dist
      cy = yp_dist
      
      

      Cam = [ [0, -1.0*focal_length/sx, -1.0*cx, 0],$    ;-1 sign change only from Barber and Redding et al., 2006 eqn 7
        [1.0*focal_length/sy, 0, -1.0*cy, 0],$
        [0, 0, 1, 0],$
        [0, 0, 0, 1]]
        
        

      q=dblarr(4)	;pixel vector
      q[0] = 1.0*ii
      q[1] = 1.0*jj
      q[2] = 1.0
      q[3] = 1.0
      
      
      
                    
                    
      sensor_m = [ [1,0,0,0],$    ;doing this to put camera North and east coords as 0 and only keep altitude, the right hand column vector is the x,y,z position of the camera, equivalent to T_i_v transformation in Barber and Redding, 2006
      [0,1,0,0],$
      [0,0,1,-1.0*camera_alt],$
      [0,0,0,1]]
      
      
     
      q_I_obj = invert(Cam ## T_mounting ## T ## sensor_m  ) ## q   ;find the vector to the surface
      

   
     ;llh_q = ecef2geodetic(q_I_obj[0],q_I_obj[1],q_I_obj[2]) ;if in ECEF coordinates
     ;z_I_obj = llh_q[2]  ;altitude
     z_I_obj = q_I_obj[2]  ;Used when using only NED coordinates
      
      
      p_I_cc = [0.0, 0.0, camera_alt, 1] ;temp for now to be north, east, and down
     ;p_I_cc = invert(T ## T_mounting ## sensor_m  ) ## [0, 0 , 0, 1]  ;from Barber et al., but not actually needed since already calculated
     ;p_I_cc = [ camera_ecef[0],camera_ecef[1],camera_ecef[2],1 ] ;In ECEF coordinates
      
      ;elev = 0.0  ;should be equal to the point where the optical axis intersects the terrain
      elev = ctrl_pixels_elev[p] 
      Z_I_cc = p_I_cc[2]  
        
      
      scale_factor =  ( elev - Z_I_cc) / (Z_I_obj - z_I_cc)  ; = lambda, or distance along camera's optical axis to object in image used to scale the vector to get to the surface
      

      
      
      p_I_obj = p_I_cc[0:3] + scale_factor*( q_I_obj[0:3] - p_I_cc[0:3])
      p_I_obj[2] = camera_alt - elev	;need to do this since p_I_obj z component is initially zero


      
      p_I_obj_ecef = NED_R_camera ## p_I_obj + camera_ecef	;put the vector to the surface in ECEF coordinates and add to the camera position to get the final pixel position in ECEF coordinates
      

      llh = ecef2geodetic(p_I_obj_ecef[0],p_I_obj_ecef[1],p_I_obj_ecef[2])	;convert from ECEF to lat/lon/altitude

            
           tot_path_length = sqrt(p_I_obj[0]^2+p_I_obj[1]^2+p_I_obj[2]^2)	;total path length
           tot_path_angle = acos( (camera_alt-elev)/tot_path_length)		;angle between the camera and the surface
     
     

           img_lat[iip] = llh[0]*180.0/!pi  ;get in degrees
           img_lon[iip] = llh[1]*180.0/!pi  ;get in degrees
           img_elev[iip] = llh[2] ;already in meters



      

      ENDFOR



  final_data = create_struct('lat',img_lat,'lon',img_lon,'elevation',img_elev)
  return,final_data

END
