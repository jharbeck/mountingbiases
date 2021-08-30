PRO CAMBOT_georeferencing_v2p0_mounting_biases

;Version 0p0 written by Al
;Version 0p1 updated by Jeremy to work on his system
;Version 0p2 adding in mounting bias section
;Version 0p3 adding in time adjustment section
;Version 0p4 adding in FoV adjustment section
;Version 0p5 updated image scenes
;Version 0p6 updated georeferencing function call
;Version 0p7 updated ac subsetting section
;Version 0p7p1 updated with additional image gridding/display sections as well as pseudo GCPs
;Version 0p7p2 adding in elevation changes to for loops section
;Version 0p7p3 switched elevation changes for focal length changes
;Version 0p7p4 updating FOR loops for tighter offsets/quicker runs, added in lens corrections
;Version 0p7p5 moving to parallel code
;Version 0p7p6 moving to multi-image support
;Version 0p8 updated to be multi campaign
;Version 0p9 updated to ingest GCP data from external files
;Version 0p10 updated to include altitude->range fix and additional ANT2018 specifics
;Version 0p10p1 updated to fix init_time issue
;Version 0p11 updated focal length to correct value given shortened CAMBOT adaptor tube, removed pulling multiple
;              lines from ancillary file and use of interp_attitude function, as we now know there is no image 
;              offset. Remove use of POS/AV data; cleaned up code.
;Version 0p12 added in lats/lons output to save files, for inter-image comparison
;Version 1p0 adds in z values for corner points and fixes the issue with using range instead of antenna altitude for the elevation
;Version 1p1 a temp version
;Version 2p0 combines updates for mounting biases for multiple campaigns

allgeo = tic('all geolocation')
firstpart = tic('first part')


;set first_run to 1 to turn off save file generation and run entire code only on the first image. Run this initially to ensure 
;CAMBOT camera image rotation/orientation makes sense and GCPs are aligned with correct X/Y locations in image, as well as get
;gross adjustments in roll, pitch and yaw. Set initial roll, pitch and yaw FOR loop values large and gross.
first_run = 0  

;~~~Select campaign to run this code for~~~
;Arctic2018
;Antarctic2018
;Arctic2019
;ArcticSummer2019
;Antarctic2019
campaign = ''

sl = path_sep()
if strcmp(sl,'\') eq 1 then begin
  base_dir1 = 'Y:' ;code locations
  base_dir2 = 'X:' ;data storage
  base_dir3 = 'W:'
endif else begin
  base_dir1 = '/data/users'
  base_dir2 = '/data/derived'
  base_dir3 = '/scratch/icebridgedata'
endelse

;Setup directories
work_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'georeferencing'+sl
lens_corr_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl

;Campaign specific variables
CASE 1 OF
  strcmp(campaign,'Arctic2018'): begin
    date = '20180314'
    cambot_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20180314'+sl+'ramp_passes'+sl+'all'+sl ;location of all images during test flights (or just a test image. Used for getting image size)
    GCP_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20180314'+sl+'ramp_passes'+sl ;location of GCP files
    GCP_LatLon_file = GCP_dir+'WFF_GCP_LatLon.csv' ;file containing lat,lon,elev locations for all GCPs
    GCP_XY_file = GCP_dir+'' ;file with row/col locations of GCPs for select images
    lens_corr_file = lens_corr_dir+'Lens_corrections_CAMBOT_15814818.h5' ;lens correction file
    cambot_anc_file = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20180314'+sl+date+'_CAMBOT_ancillary.csv'
  end
  
  strcmp(campaign,'Antarctic2018'): begin
    date = '20181002'
    cambot_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF2-20181002'+sl+'ramp_passes'+sl+'all'+sl
    GCP_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF2-20181002'+sl+'ramp_passes'+sl
    GCP_LatLon_file = GCP_dir+'Palmdale_GCP_LatLon.csv'
    GCP_XY_file = GCP_dir+'CAMBOT_GCP_XY_Palmdale_orthogonal.csv'
    lens_corr_file = lens_corr_dir+'Lens_corrections_CAMBOT_51500462.h5'
    cambot_anc_file = base_dir2+sl+strmid(date,0,4)+'_fall_quicklook'+sl+'CAMBOT'+sl+date+'_CAMBOT_ancillary_precise.csv'
  end
    
  strcmp(campaign,'Arctic2019'): begin
    date = '20190327'
    cambot_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20190327'+sl+'ramp_passes'+sl+'all'+sl
    GCP_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20190327'+sl+'ramp_passes'+sl
    GCP_LatLon_file = GCP_dir+'WFF_GCP_LatLon.csv'
    GCP_XY_file = GCP_dir+'CAMBOT_GCP_XY_WFF_20190327_Try11.csv'
    lens_corr_file = lens_corr_dir+'Lens_corrections_CAMBOT_51500462.h5'
    cambot_anc_file = base_dir3+sl+'CAMBOT_JPEG'+sl+'ancillary_files'+sl+date+'_CAMBOT_ancillary_precise.csv'
  end
  
  strcmp(campaign,'ArcticSummer2019'): begin
    date = '20190819'
    cambot_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20190819'+sl+'ramp_passes'+sl+'all'+sl
    GCP_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20190819'+sl+'ramp_passes'+sl
    GCP_LatLon_file = GCP_dir+'KLBX_GCP_LatLon.csv'
    GCP_XY_file = GCP_dir+'CAMBOT_GCP_XY_KLBX_20190819.csv'
    lens_corr_file = lens_corr_dir+'Lens_corrections_CAMBOT_51500462.h5'
    cambot_anc_file = base_dir2+sl+strmid(date,0,4)+'_summer_quicklook'+sl+'CAMBOT'+sl+date+'_CAMBOT_ancillary_precise.csv'
  end
    
  strcmp(campaign,'Antarctic2019'): begin
    date = '20191017'
    cambot_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20191017'+sl+'ramp_passes'+sl+'all'+sl
    GCP_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20191017'+sl+'ramp_passes'+sl
    GCP_LatLon_file = GCP_dir+'KLBX_GCP_LatLon.csv'
    GCP_XY_file = GCP_dir+'CAMBOT_GCP_XY_KLBX_20191017.csv'
    lens_corr_file = lens_corr_dir+'Lens_corrections_CAMBOT_51500462.h5'
    cambot_anc_file = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20191017'+sl+date+'_CAMBOT_ancillary_precise.csv'
  end
  else: STOP,'No campaign correctly identified'
endcase

;output directories
worker_output_dir = cambot_dir+'worker_output'+sl
out_dir = GCP_dir+'mounting_bias_runs'+sl


;check directories/files
if file_test(GCP_dir,/dir) eq 0 then STOP,'cannot find GCP directory'
if file_test(GCP_LatLon_file) eq 0 then STOP,'cannot find GCP LatLon file'
if file_test(work_dir,/dir) eq 0 then STOP,'cannot find work directory'
if file_test(lens_corr_file) eq 0 then STOP,'cannot find lens correction file'
if file_test(cambot_dir,/dir) eq 0 then STOP,'cannot find CAMBOT directory'
if file_test(cambot_anc_file) eq 0 then STOP,'cannot find CAMBOT ancillary file'
if file_test(out_dir) eq 0 then file_mkdir,out_dir

verbose = 0
t1 = systime(1)

;----------------------------------------------------------------------------------
;                                import GCP files
;----------------------------------------------------------------------------------
result = query_ascii(GCP_LatLon_file,info)
GCP_coords = replicate(-99999.0,info.lines-1,3) ;latlon coordinates (& elevation)
GCP_pixels_temp = replicate(-99999.0,info.lines-1) ;xy coordinates
GCP_names = replicate('',info.lines-1)
tmp = ''
openr,lun,GCP_LatLon_file,/get_lun
readf,lun,tmp ;header
for n=0,info.lines-2 do begin
  readf,lun,tmp 
  tmp2 = strsplit(tmp,',',/extract)
  GCP_names[n] = strcompress(tmp2[0],/remove_all)
  GCP_coords[n,*] = float(tmp2[1:3])
endfor
free_lun,lun


if file_test(GCP_XY_file) eq 0 then STOP,'cannot find GCP XY file'
result = query_ascii(GCP_XY_file,info2)

;read through file one time first and find how many images we are working with
GCP_XY_names = replicate('',info2.lines-1)
tmp = ''
openr,lun,GCP_XY_file,/get_lun
readf,lun,tmp ;header
for n=0,info2.lines-2 do begin
  readf,lun,tmp
  tmp2 = strsplit(tmp,',',/extract)
  if strcmp(strmid(tmp2[0],0,1),'-') eq 0 then GCP_XY_names[n] = tmp2[0] 
endfor
free_lun,lun

tempsel = where(GCP_XY_names ne '')
GCP_XY_names = GCP_XY_names[tempsel]
sel_images = uniq(GCP_XY_names,sort(GCP_XY_names))
if strcmp(sel_images[0],'') eq 1 then STOP,'cannot find image names'
if n_elements(sel_images) lt 1 then STOP,'cannot find image names' 
image_names = GCP_XY_names[sel_images]
nimages = n_elements(image_names)

;define and build a structure to for each individual image to use
GCP_struct_def2 = {GCP_struct2, $
                  X:GCP_pixels_temp, $
                  Y:GCP_pixels_temp, $
                  Z:GCP_pixels_temp}
GCP_pixels_struct = replicate(GCP_struct_def2,nimages)

openr,lun,GCP_XY_file,/get_lun
readf,lun,tmp ;header
for n=0,info2.lines-2 do begin
  readf,lun,tmp
  tmp2 = strsplit(tmp,',',/extract)
  if strcmp(strmid(tmp2[0],0,1),'-') eq 0 then begin
    sel_image = where(tmp2[0] eq image_names)
    if sel_image eq -1 then STOP,'image name matching issue '
    sel_point = where(tmp2[1] eq GCP_names)
    if sel_point eq -1 then STOP,'point name matching issue'
    
    GCP_pixels_struct[sel_image].X[sel_point] = tmp2[2]
    GCP_pixels_struct[sel_image].Y[sel_point] = tmp2[3] 
  endif
  
endfor
free_lun,lun
;==================================================================================


;----------------------------------------------------------------------------------
;                  start FOR loop going through all images here
;----------------------------------------------------------------------------------

if first_run eq 1 then nimages=1
FOR ni = 0,nimages-1 DO BEGIN
  
  ;assign variables to daily ones
  file_in = image_names[ni]+'.jpg'
  GCP_pixels = replicate(-99999.0,info.lines-1,2)
  GCP_pixels[*,0] = GCP_pixels_struct[ni].X
  GCP_pixels[*,1] = GCP_pixels_struct[ni].Y
  
  ;find active GCPs from entered pixel locations above and assign them to variables used below
  cp_sel = where(GCP_pixels[*,0] gt -99999,npt)
  ctrl_pixels = GCP_pixels[cp_sel,*]
  ctrl_coords = GCP_coords[cp_sel,*]
  this_GCP_names = GCP_names[cp_sel]
  ctrl_pixels_elev = ctrl_coords[*,2] ;assign GCP elevtions to each pixel
  
  ;find and read in CAMBOT image
  label = file_basename(file_in,'.jpg')
  
  file_in_full = cambot_dir + file_in
  if file_test(file_in_full) eq 0 then STOP,'cannot find input image'
  img_in = read_image(file_in_full)
  ;help, img
  
  ;; just use a single band
  img = reform(img_in(0, *, *))
  
  ;do this to determine x&y pixel locations to match with lat/lon values. This is with the lower, long edge of the CAMBOT image in the direction of flight.
  ;test = image(img_in) 
  
  ;----------------------------------------------------------------------------
  ;rotate input image and GCP coordinates so that flight direction is upwards
  ;----------------------------------------------------------------------------
  ;;swap fov_x/y accordingly
  ;;also orient the image so that up is azimuth 0 (north)
  ;img = rotate(img, 3) ;DMS
  ;img = rotate(img, 2) ;CAMBOT Arctic 2017
  case 1 of 
    strcmp(campaign,'Arctic2018'): begin
                                    ;***no rotation
                                   end
    
    strcmp(campaign,'Antarctic2018'): begin  
                                        img_size = size(img)
                                        xsize = img_size[1]
                                        ysize = img_size[2]
        
                                        img = rotate(img, 2) ;image
                                        ctrl_pixels[*,0] = xsize-ctrl_pixels[*,0] ;GCP pixel locations to match image rotation
                                        ctrl_pixels[*,1] = ysize-ctrl_pixels[*,1]
                                      end
    strcmp(campaign,'Arctic2019'): begin 
                                    ;CAMBOT Arctic 2019
                                    ;***no rotation
                                   end
    
    strcmp(campaign,'ArcticSummer2019'): begin 
                                          ;CAMBOT Arctic Summer 2019
                                          img_size = size(img)
                                          xsize = img_size[1]
                                          ysize = img_size[2]
                                      
                                          img = rotate(img, 2) ;image
                                          ctrl_pixels[*,0] = xsize-ctrl_pixels[*,0] ;GCP pixel locations to match image rotation
                                          ctrl_pixels[*,1] = ysize-ctrl_pixels[*,1]
                                        end
    
    strcmp(campaign,'Antarctic2019'): begin 
                                        ;CAMBOT Antarctic 2019
                                        img_size = size(img)
                                        xsize = img_size[1]
                                        ysize = img_size[2]
                                    
                                        img = rotate(img, 2) ;image
                                        ctrl_pixels[*,0] = xsize-ctrl_pixels[*,0] ;GCP pixel locations to match image rotation
                                        ctrl_pixels[*,1] = ysize-ctrl_pixels[*,1]
                                      end
    
  ;============================================================================
  
  ;----------------------------------------------------------------------------
  ;           correct the longitude mirror flip issue (all campaigns)
  ;----------------------------------------------------------------------------
  ;applied for all campaigns so far...
  img = reverse(img, 1)
  
  img_size = size(img)
  xsize = img_size[1]
  ysize = img_size[2]
  
  ;flip x-coordinates to match reverse above
  ctrl_pixels[*,0] = xsize-ctrl_pixels[*,0]
  ;============================================================================
      
  ;----------------------------------------------------------------------------
  ;
  ;----------------------------------------------------------------------------
  img_ncol = xsize
  img_nrow = ysize
  
  ;pull out timestamp from CAMBOT filename
  file_time = strmid(file_in,9,11)
  
  ;pull lens correction data out from the file provided
  file_id = H5F_OPEN(lens_corr_file)
  ID = H5D_OPEN(file_id,'xcorr')
  xcorr = H5D_READ(ID)
  ID = H5D_OPEN(file_id,'ycorr')
  ycorr = H5D_READ(ID)
  
  ;transpose corrections from row (MATLAB) to column major (IDL)
  xcorr = transpose(xcorr)
  ycorr = transpose(ycorr)
  
  ;read-in CAMBOT ancillary data file
  anc_data =  read_cambot_ancillary_file_v1p0(cambot_anc_file)
  asel = (where(anc_data.filename eq file_in,nsel))[0]
  if nsel ne 1 then STOP,'cannot find a filename match in the ancillary file'

  unix_offset = (julday(strmid(date,4,2),strmid(date,6,2),strmid(date,0,4))-julday(1,1,1970))*24.0*3600.0
  
  ac_time = anc_data.posix_time[asel]-unix_offset
  ac_lat = anc_data.lat[asel]
  ac_lon = anc_data.lon[asel]
  ac_elev = anc_data.range[asel] ;selects range from ancillary file
  ;ac_elev = anc_data.alt[asel] ;selects altitude from ancillary file
  ac_roll = anc_data.roll[asel]
  ac_pitch = anc_data.pitch[asel]
  ac_azim = anc_data.heading[asel]
  
  ;; - the aircraft heading is positive clockwise (as one would expect)
  ;; - the function geolocate_image_v4p1 seems to implement the heading as positive counter-clockwise
  ;; - multiply the heading by -1 for use with the function
  ;; - a headin of zero points the top of the image to the south, so add 180 to make it point north
  ;; - the top of the image is not the direction of travel.  Add 90 (for dms) so that the image is rotated to the direction of travel
  
  label_ver = '_00'
  
  ;----------------------------------------------------------------------------------
  ;                      camera settings and mounting bias values 
  ;----------------------------------------------------------------------------------
  
  ;lens values from pointsinfocus.com lens calculator with manual sensor size value of: 26.93mm x 17.95mm
  fov_x = 51.4
  fov_y = 35.5
  sensor_x = 0.02693 ;sensor x dimension in meters
  sensor_y = 0.01795
  
  strcmp(campaign,'Arctic2018'):       camera_offset = [-4.148D,0.25D,3.2D] ;offset for P-3B for 2018, from 20180314 ancillary file
  strcmp(campaign,'Antarctic2018'):    camera_offset = [-2.7D,-0.4D,3.95D] ;offset for DC-8 for ANT 2018, from 20181002 ancillary file
  strcmp(campaign,'ArcticSummer2019'): camera_offset = [-4.148D,0.25D,3.2D] ;offset for P-3B for 2019, from 20190327 ancillary file
  strcmp(campaign,'Antarctic2019'):    camera_offset = [-4.463D,0.092D,2.042D] ;offset for G-V for summer 2019, from 20190819 ancillary file 
   
  focal_length =  0.02859 ;lens focal length (Arctic 2019 : calc from min mean dist, try 13i
  
  
  
  ;~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~   
  ;These values are all manually entered and updated as each iteration gets us closer to the
  ;   actual mounting bias value for a specific campaign
  ;~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

  ;mounting bias values (start with zeros and update as necessary)
  pitch_bias = 0.0 ;-4.09 ;-4.08
  roll_bias = 0.0 ;-1.10 ;-1.11
  heading_bias = 0.00 ;0.50 ;0.46
  
  pitch_bias_init = pitch_bias
  roll_bias_init = roll_bias
  heading_bias_init = heading_bias
  
  ;setup +/- values we want to iterate through for this set of mounting biases 
  count_iter = 0L
  
  pm_check = 5.0 ;how far +/- in degrees we want to search out in our loops for roll & pitch mounting bias
  dpm = 0.1
  
  hm_check = 5.0 ; ;how far +/- in degrees we want to search in our loop for heading bias offsets
  dhm = 0.1
  
  bias_len = ulong((pm_check/dpm*2+1)^2 * (hm_check/dhm*2+1) * 1.05) ;multiply by 1.2 just to give us some buffer. change ^2 to ^3 if you add in heading. 
  
  mcp_dist = fltarr(bias_len)
  mcp_dist_all = fltarr(bias_len,npt)
  pitch_vals = fltarr(bias_len)
  roll_vals = fltarr(bias_len)
  heading_vals = fltarr(bias_len)
  
  pcnt_update = 16
  update_point = long(bias_len/pcnt_update)
  
  print,'Beginning mounting bias loop values generation'
      
  
  ;now iterate through this timing offset with all possibilities of mounting biases
  FOR hh=-1.0*hm_check,hm_check,dhm DO BEGIN ;heading offset
    FOR rr=-1.0*pm_check,pm_check,dpm DO BEGIN
      FOR pp=-1.0*pm_check,pm_check,dpm DO BEGIN
  
        pitch_vals[count_iter] = pp ;pitch_bias_in
        roll_vals[count_iter] = rr ;roll_bias_in
        heading_vals[count_iter] = hh ;heading_bias_in
        count_iter++
  
      ENDFOR  ;end loop over pitch
    ENDFOR  ;end loop over roll
  ENDFOR ;end loop over heading
  
  ;truncate arrays to valid values
  mcp_dist = mcp_dist[0L:count_iter-1L]
  mcp_dist_all = mcp_dist_all[0L:count_iter-1,*]
  mcp_dist_3D = mcp_dist
  mcp_dist_all_3D = mcp_dist_all
  mcp_lats = mcp_dist_all
  mcp_lons = mcp_dist_all
  pitch_vals = pitch_vals[0L:count_iter-1L]
  roll_vals = roll_vals[0L:count_iter-1L]
  heading_vals = heading_vals[0L:count_iter-1L]
  
  ;add-in ititial values here so we're only working with the full values from now on
  pitch_vals = pitch_vals+pitch_bias_init
  roll_vals = roll_vals+roll_bias_init
  heading_vals = heading_vals+heading_bias_init
  toc,firstpart
  
  ;---------------------------------------------------------------------------------------------------
  ;                                 parallel part of the code
  ;---------------------------------------------------------------------------------------------------
  proj_info = map_proj_init('Polar Stereographic',center_longitude=ac_lon,true_scale_latitude=ac_lat)
  
  ;save structures to an IDL save file, for use later in the child processes (as I can't get the structure passing to work)
  proj_file = work_dir+'temp_proj_file_'+strcompress(floor(systime(/seconds)),/rem)+'.sav'
  save,filename=proj_file,proj_info
  
  data_index = bytarr(count_iter)
  ntodo = count_iter-1L
  go_parallel = 1
  
  IF go_parallel eq 1 THEN BEGIN
  
    parallelsection = tic('Parallel section')
    nthreads = 16
    
    SPLIT_FOR_v2,0,ntodo,ctvariable_name='ci',commands=[$
      'pitch_bias_in = pitch_vals[ci]',$
      'roll_bias_in = roll_vals[ci]',$
      'heading_bias_in = heading_vals[ci]',$
      'floc = geolocate_ControlPts_v6p0(img_ncol,img_nrow,ctrl_pixels,ctrl_pixels_elev,ac_lat,ac_lon,ac_elev,ac_pitch,ac_roll,ac_azim,focal_length,sensor_x,sensor_y,camera_offset,pitch_bias=pitch_bias_in,roll_bias=roll_bias_in,heading_bias=heading_bias_in,xcorr=xcorr,ycorr=ycorr)',$
      'szl = size(floc.lat,/dim)',$
      'mcp_lats[ci,*] = floc.lat',$
      'mcp_lons[ci,*] = floc.lon',$
      'c_dist = fltarr(npt)',$
      'FOR pt=0,npt-1 DO c_dist[pt] = map_2points(floc.lon[pt],floc.lat[pt],ctrl_coords[pt,1],ctrl_coords[pt,0],/meters)',$
      'mcp_dist_all[ci,*] = c_dist',$
      'mcp_dist[ci] = mean(c_dist)',$
      'restore,proj_file',$
      'floc_xy = map_proj_forward(floc.lon,floc.lat,map_structure=proj_info)',$
      'ctrl_coords_xy = map_proj_forward(ctrl_coords[*,1],ctrl_coords[*,0],map_structure=proj_info)',$
      'c_dist_3D = sqrt( (floc_xy[0,*]-ctrl_coords_xy[0,*])^2 + (floc_xy[1,*]-ctrl_coords_xy[1,*])^2 + (floc.elevation-ctrl_coords[*,2])^2)',$
      'mcp_dist_all_3D[ci,*] = c_dist_3D',$
      'mcp_dist_3D[ci] = mean(c_dist_3D)',$
      'data_index[ci] = 1B'],$
      varnames = ['pitch_vals','roll_vals','heading_vals','img_ncol','img_nrow','ctrl_pixels','ctrl_pixels_elev','ac_lat','ac_lon','ac_elev','ac_pitch','ac_roll','ac_azim','focal_length','sensor_x','sensor_y','camera_offset','xcorr','ycorr','mcp_lats','mcp_lons','npt','ctrl_coords','mcp_dist_all','mcp_dist','data_index','mcp_dist_3D','mcp_dist_all_3D','proj_file'],$
      outvar = ['data_index','mcp_dist','mcp_dist_all','mcp_dist_3D','mcp_dist_all_3D','mcp_lats','mcp_lons'],$
      nsplit=nthreads
      
    ;combine thread output together into single threads
    FOR t=0,nthreads-1 DO BEGIN
      st = strcompress(t,/remove_all)
      tmp1 = execute('sel = where(data_index'+st+' eq 1B)')
      tmp2 = execute('mcp_dist[sel] = mcp_dist'+st+'[sel]')
      tmp3 = execute('mcp_dist_all[sel,*] = mcp_dist_all'+st+'[sel,*]')
      tmp4 = execute('mcp_lats[sel,*] = mcp_lats'+st+'[sel,*]')
      tmp5 = execute('mcp_lons[sel,*] = mcp_lons'+st+'[sel,*]')
      tmp6 = execute('mcp_dist_3D[sel] = mcp_dist_3D'+st+'[sel]')
      tmp7 = execute('mcp_dist_all_3D[sel,*] = mcp_dist_all_3D'+st+'[sel,*]')
    ENDFOR
      
    toc,parallelsection
  
  ENDIF ELSE BEGIN
  
    ;not parallel section (for debugging)
    FOR ci=0,ntodo DO BEGIN
      pitch_bias_in = pitch_vals[ci]
      roll_bias_in = roll_vals[ci]
      heading_bias_in = heading_vals[ci]
      floc = geolocate_ControlPts_v6p0(img_ncol,img_nrow,ctrl_pixels,ctrl_pixels_elev,ac_lat,ac_lon,ac_elev,ac_pitch,ac_roll,ac_azim,focal_length,sensor_x,sensor_y,camera_offset,pitch_bias=pitch_bias_in,roll_bias=roll_bias_in,heading_bias=heading_bias_in,xcorr=xcorr,ycorr=ycorr)
      szl = size(floc.lat,/dim)
      mcp_lats[ci,*] = floc.lat
      mcp_lons[ci,*] = floc.lon
      c_dist = fltarr(npt)
      FOR pt=0,npt-1 DO c_dist[pt] = map_2points(floc.lon[pt],floc.lat[pt],ctrl_coords[pt,1],ctrl_coords[pt,0],/meters)
      mcp_dist_all[ci,*] = c_dist
      mcp_dist[ci] = mean(c_dist)
      restore,proj_file
      floc_xy = map_proj_forward(floc.lon,floc.lat,map_structure=proj_info)
      ctrl_coords_xy = map_proj_forward(ctrl_coords[*,1],ctrl_coords[*,0],map_structure=proj_info)
      c_dist_3D = sqrt( (floc_xy[0,*]-ctrl_coords_xy[0,*])^2 + (floc_xy[1,*]-ctrl_coords_xy[1,*])^2 + (floc.elevation-ctrl_coords[*,2])^2)
      mcp_dist_all_3D[ci,*] = c_dist_3D
      mcp_dist[ci] = mean(c_dist_3D)
      data_index[ci] = 1B
    ENDFOR
  
  ENDELSE
  
  if file_test(proj_file) eq 1 then file_delete,proj_file
  
  ;---------------------------------------------------------------------------------------------------
  ;                              export mounting bias runs stats
  ;---------------------------------------------------------------------------------------------------
  out_file = out_dir+file_basename(file_in,'.jpg')+'_mounting_bias_3D.sav'
  if first_run eq 0 then save,filename=out_file,mcp_dist,mcp_dist_all,mcp_dist_3D,mcp_dist_all_3D,pitch_vals,roll_vals,heading_vals,timing_vals,mcp_lats,mcp_lons,npt,this_GCP_names
  
  ;---------------------------------------------------------------------------------------------------
  ;                                 find the best mounting bias
  ;---------------------------------------------------------------------------------------------------
  FOR b=0,npt-1 DO BEGIN
    best = min(mcp_dist_all[*,b],best_idx)
    Print,'GCP: '+strcompress(b)
    ;print,'Best value for pitch, roll, heading, time offset, alt offset: ',pitch_vals[best_idx],roll_vals[best_idx],heading_vals[best_idx],timing_vals[best_idx],altitude_vals[best_idx]
    print,'Best value for pitch, roll, heading, time offset: ',pitch_vals[best_idx],roll_vals[best_idx],heading_vals[best_idx],timing_vals[best_idx]
    print,'Best value for distance: ',mcp_dist_all[best_idx,b]
    print,'---------------------------------------------------'
    ;print,'Distance between control points: ',mcp_dist_all[best_idx,*]
  ENDFOR
  
  best = min(mcp_dist,best_idx) ;all-points minimum distance method
  Print,'Best fit overall for all points'
  ;print,'Best value for pitch, roll, heading, time offset, alt offset: ',pitch_vals[best_idx],roll_vals[best_idx],heading_vals[best_idx],timing_vals[best_idx],altitude_vals[best_idx]
  print,'Best value for pitch, roll, heading, time offset: ',pitch_vals[best_idx],roll_vals[best_idx],heading_vals[best_idx],timing_vals[best_idx]
  print,'Best value for distance: ',mcp_dist[best_idx]
  print,'==================================================='
  toc,allgeo

ENDFOR ;end FOR loop over all images here

Print,'Finished with current run of mounting bias code.'

end