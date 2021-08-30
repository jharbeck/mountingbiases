PRO mounting_bias_best_fit_v1p4
;---------------------------------------------------------------
;Program to read in mounting bias files from individual images
;and find the best fit to all of them.
;
;Version 1p0: 20181010 by Jeremy Harbeck
;Version 1p1: 20181130, adjusted save file in values to account for
;                       timing adjustments & removal of init values
;Version 1p2: 20190312 - updated to remove timing offsets
;Version 1p3: 20200225 - updated to include 3D distances
;Version 1p4: 20200226 - a temp version to look at focal lengths
;--------------------------------------------------------------- 

;date = '20181002'
date = '20190327'
;date = '20190819'

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

Print,'Processing date: '+date

;setup directories

;Arctic 2018
;in_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20180314'+sl+'ramp_passes'+sl+'mounting_bias_runs'+sl

;Antarctic 2018
;in_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF2-20181002'+sl+'ramp_passes'+sl+'mounting_bias_runs'+sl

;Arctic 2019
in_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20190327'+sl+'ramp_passes'+sl+'mounting_bias_runs'+sl

;Arctic Summer 2019
;in_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20190819'+sl+'ramp_passes'+sl+'mounting_bias_runs'+sl

;Antarctic 2019
;in_dir = base_dir1+sl+'jharbeck'+sl+'Working'+sl+'Cambot'+sl+'PCF1-20191017'+sl+'ramp_passes'+sl+'mounting_bias_runs'+sl

if file_test(in_dir,/dir) eq 0 then STOP,'cannot find primary directory'

files = file_search(in_dir+'*.sav')
if n_elements(files) lt 1 then STOP,'did not find any mounting bias files'

;setup structure to read data from all files into
timing_bias_init = 0.0 ;-0.4338

FOR f=0,n_elements(files)-1 do begin
  file = files[f]
  restore,file,/verbose
  ;restores:mcp_dist,mcp_dist_all,pitch_vals,roll_vals,heading_vals,timing_vals,npt
  
  if f eq 0 then begin
    pitch_bias_in = pitch_vals
    roll_bias_in = roll_vals
    heading_bias_in = heading_vals
    timing_vals_in = timing_vals
    focal_vals_in = focal_vals
    
    names = replicate('',n_elements(files))
    stats_median = replicate(-1.0,n_elements(files),n_elements(pitch_vals))
    stats_mean = replicate(-1.0,n_elements(files),n_elements(pitch_vals))
    stats_median_3D = stats_median
    stats_mean_3D = stats_mean
    
  endif else begin
    ;ensure the roll,pitch,heading and timing values from this file match the initial one so we're comparing apples to apples
    if max(abs(pitch_bias_in-pitch_vals)) gt 0 then STOP
    if max(abs(roll_bias_in-roll_vals)) gt 0 then STOP
    if max(abs(heading_bias_in-heading_vals)) gt 0 then STOP
    if max(abs(timing_vals_in-timing_vals)) gt 0 then STOP
    if max(abs(focal_vals_in-focal_vals)) gt 0 then STOP
  endelse
  
  ;calculate mean mcp_dist based on npt and save to stats variable; assign name of file to varibale
  names = files[f]
  ;stats[f,*] = mcp_dist/float(npt)
  stats_median[f,*] = median(mcp_dist_all,dimension=2) ;calculate median distance from all GCPs to their respective point for each roll, pitch & heading permutation
  stats_mean[f,*] = mean(mcp_dist_all,dimension=2) ;calculate median distance from all GCPs to their respective point for each roll, pitch & heading permutation
  stats_median_3D[f,*] = median(mcp_dist_all_3D,dimension=2)
  stats_mean_3D[f,*] = mean(mcp_dist_all_3D,dimension=2)
  
ENDFOR

;decide which mounting bias combo is best. Take the smallest distance mean.
stats_median = median(stats_median,dimension=1)
stats_mean = mean(stats_mean,dimension=1)
stats_median_3D = median(stats_median_3D,dimension=1)
stats_mean_3D = mean(stats_mean_3D,dimension=1)

best_median = min(stats_median,best_median_idx)
best_mean = min(stats_mean,best_mean_idx)
best_median_3D = min(stats_median_3D,best_median_idx_3D)
best_mean_3D = min(stats_mean_3D,best_mean_idx_3D)

print,'------------------------------------------------------------------------'
print,'Best value for lowest median distance: ',stats_median[best_median_idx]
print,'Best median value for pitch, roll, heading, time offset, focal length: ',pitch_bias_in[best_median_idx],roll_bias_in[best_median_idx],heading_bias_in[best_median_idx],timing_vals_in[best_median_idx];,focal_vals_in[best_mean_idx_3D]
;print,'Best value for mean summed distance: ',mcp_dist[best_idx]
print,'------------------------------------------------------------------------'
print,'Best value for lowest mean distance: ',stats_mean[best_mean_idx]
print,'Best mean value for pitch, roll, heading, time offset, focal length: ',pitch_bias_in[best_mean_idx],roll_bias_in[best_mean_idx],heading_bias_in[best_mean_idx],timing_vals_in[best_mean_idx];,focal_vals_in[best_mean_idx_3D]
print,'------------------------------------------------------------------------'
print,''
print,'------------------------------------------------------------------------'
print,'Best value for lowest 3D median distance: ',stats_median_3D[best_median_idx_3D]
print,'Best median value for 3D pitch, roll, heading, time offset, focal length: ',pitch_bias_in[best_median_idx_3D],roll_bias_in[best_median_idx_3D],heading_bias_in[best_median_idx_3D],timing_vals_in[best_median_idx_3D];,focal_vals_in[best_mean_idx_3D]
;print,'Best value for mean summed distance: ',mcp_dist[best_idx]
print,'------------------------------------------------------------------------'
print,'Best value for lowest 3D mean distance: ',stats_mean_3D[best_mean_idx_3D]
print,'Best mean value for 3D pitch, roll, heading, time offset, focal length: ',pitch_bias_in[best_mean_idx_3D],roll_bias_in[best_mean_idx_3D],heading_bias_in[best_mean_idx_3D],timing_vals_in[best_mean_idx_3D];,focal_vals_in[best_mean_idx_3D]
print,'------------------------------------------------------------------------'

STOP

END