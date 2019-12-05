;*********************************************************************
; Calculate trend in interhemispheric asymmetry in precip and SST 
; for HIST or RCP runs and plot as scatter plot 
; (Rotstayn et al. 2015 Fig 4)
;
; JRB 05.04.18
;*********************************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  modall = (/ "bcc-csm1-1","bcc-csm1-1-m","BNU-ESM","CCSM4","CESM1-BGC",\
              "CESM1-WACCM","EC-EARTH","FGOALS-s2","FIO-ESM",\
              "GFDL-ESM2G","GFDL-ESM2M","MPI-ESM-LR","MPI-ESM-MR", \
              "ACCESS1-0","ACCESS1-3","CESM1-CAM5","CSIRO-Mk3-6-0",\
              "FGOALS-g2","GFDL-CM3","HadGEM2-CC","HadGEM2-ES",\
              "MIROC-ESM","MIROC-ESM-CHEM","MIROC5","MRI-CGCM3",\
              "NorESM1-M","NorESM1-ME" /)


 model = modall

 nmod = 27
 
 modnum_lo = (/0,1,2,3,4,5,6,7,8,9,10,11,12/)
 modnum_hi = (/13,14,15,16,17,18,19,20,21,22,23,24,25,26/)

 nlat = 120
 nlon = 240

; Domains

 latmin1 =  -20.
 latmax1 =    0.
 
 latmin2 =    0.
 latmax2 =   20.

 latmin3 =   -90.
 latmax3 =    90.

 lonmin =   0.  
 lonmax = 360.

 nyr = 106  ; 1900-2005
 ntime = nyr*12
 
 tempseas1 = new((/nmod,nyr,nlat,nlon/),float)

 nyr2 = 95  ; 2006-2100
 ntime2 = nyr2*12

 tempseas2 = new((/nmod,nyr2,nlat,nlon/),float)

; combined variables

 nyr_all = nyr + nyr2
 ntime_all = nyr_all*12

 xtime1 = ispan(0,nyr-1,1)
 xtime2 = ispan(0,nyr2-1,1)

 temp = new((/nmod,ntime_all,nlat,nlon/),float)
 temp_seas = new((/nmod,nyr_all,nlat,nlon/),float)

 temp_djf = new((/nmod,nyr_all,nlat,nlon/),float) 
 temp_djf_ocean = new((/nmod,nyr_all,nlat,nlon/),float)
 temp_areaave1 = new((/nmod,nyr_all/),float)
 temp_areaave2 = new((/nmod,nyr_all/),float)

 sst_asym = new((/nmod,nyr_all/),float)
 sst_trend = new((/nmod/),float)

 prec = new((/nmod,ntime_all,nlat,nlon/),float)
 prec_seas = new((/nmod,nyr_all,nlat,nlon/),float)

 prec_djf = new((/nmod,nyr_all,nlat,nlon/),float) 
 prec_djf_ocean = new((/nmod,nyr_all,nlat,nlon/),float)
 prec_areaave1 = new((/nmod,nyr_all/),float)
 prec_areaave2 = new((/nmod,nyr_all/),float)

 prec_asym = new((/nmod,nyr_all/),float)
 prec_trend = new((/nmod/),float)

 path1 = "/g/data1/eg3/pfd548/IPCC/CMIP5/historical/ts_processed/grid1.5"
 path2 = "/g/data1/eg3/pfd548/IPCC/CMIP5/rcp85/ts_processed/grid1.5"

 path3 = "/g/data1/eg3/pfd548/IPCC/CMIP5/historical/pr_processed/grid1.5"
 path4 = "/g/data1/eg3/pfd548/IPCC/CMIP5/rcp85/pr_processed/grid1.5"

 maskfile = addfile("/g/data/eg3/jxb548/CMIP5/z7000.nc","r")
 mskdats = maskfile->data(0,:,:)
 mskdat = short2flt(mskdats)

 do ii = 0,nmod-1 ; loop over models

 print(model(ii))

 ; a1 = addfile(path1+"/ts_"+model(ii)+"_r1i1p1_1900_2005_r240x120.nc","r") 
 ; ts1 = a1->ts(0:ntime-1,:,:)
 ; ts1 = ts1 - 273.15
 
 ; a2 = addfile(path2+"/ts_"+model(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 
 ; ts2 = a2->ts(0:ntime2-1,:,:)
 ; ts2 = ts2 - 273.15

  if (ii .ne. 6) then
   a1 = addfile(path1+"/ts_"+model(ii)+"_r1i1p1_1900_2005_r240x120.nc","r") 
  end if
  if (ii .eq. 6) then
   a1 = addfile(path1+"/ts_"+model(ii)+"_r8i1p1_1900_2005_r240x120.nc","r") 
  end if
  ts1 = a1->ts(0:ntime-1,:,:)
  ts1 = ts1 - 273.15

  if (ii .ne. 5 .and. ii .ne. 6) then
   a2 = addfile(path2+"/ts_"+model(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 
   end if
  if (ii .eq. 5 ) then
   a2 = addfile(path2+"/ts_"+model(ii)+"_r2i1p1_2006_2100_r240x120.nc","r") 
  end if
  if (ii .eq. 6 ) then
   a2 = addfile(path2+"/ts_"+model(ii)+"_r8i1p1_2006_2100_r240x120.nc","r") 
  end if
  ts2 = a2->ts(0:ntime2-1,:,:)
  ts2 = ts2 - 273.15

  ntime_end2 = ntime + ntime2 - 1

  temp(ii,0:ntime-1,:,:) = ts1
  temp(ii,ntime:ntime_end2,:,:) = ts2

 ; temp_seas(ii,:,:,:) = month_to_season(temp(ii,:,:,:),"DJF")
  temp_seas(ii,:,:,:) = month_to_annual(temp(ii,:,:,:),1)

;  temp_djf(ii,:,:,:) = temp_seas(ii,1:nyr_all-2,:,:)
  temp_djf(ii,:,:,:) = temp_seas(ii,:,:,:)

  temp_djf_ocean(ii,:,:,:) = temp_djf(ii,:,:,:)
  temp_djf_ocean(ii,:,:,:) = mask(temp_djf(ii,:,:,:),mskdat,0)  ; 0 is ocean

; select LAND, OCEAN or LAND+OCEAN:

  temp_areaave1(ii,:) = wgt_areaave(temp_djf(ii,:,{latmin1:latmax1},{lonmin:lonmax}), 1, 1.0, 0) 
  temp_areaave2(ii,:) = wgt_areaave(temp_djf(ii,:,{latmin3:latmax3},{lonmin:lonmax}), 1, 1.0, 0) 

;  sst_asym(ii,:) = temp_areaave2(ii,:) - temp_areaave1(ii,:)
 sst_asym(ii,:) = temp_areaave2(ii,:)

; calculate trend over period of interest

 ; sst_trend(ii) = regCoef(xtime1, sst_asym(ii,0:105))
  sst_trend(ii) = regCoef(xtime2, sst_asym(ii,106:200))

; repeat for precip

  b1 = addfile(path3+"/pr_"+model(ii)+"_r1i1p1_1900_2005_r240x120.nc","r") 
  pr1 = b1->pr(0:ntime-1,:,:)
  pr1 = pr1*24*60*60
 
 ; b2 = addfile(path4+"/pr_"+model(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 
 ; pr2 = b2->pr(0:ntime2-1,:,:)
 ; pr2 = pr2*24*60*60

  if (ii .ne. 5 .and. ii .ne. 6) then
   b2 = addfile(path4+"/pr_"+model(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 
   end if
  if (ii .eq. 5 .or. ii .eq. 6) then
   b2 = addfile(path4+"/pr_"+model(ii)+"_r2i1p1_2006_2100_r240x120.nc","r") 
  end if
  pr2 = b2->pr(0:ntime2-1,:,:)
  pr2 = pr2*24*60*60

  prec(ii,0:ntime-1,:,:) = pr1
  prec(ii,ntime:ntime_end2,:,:) = pr2

 ; prec_seas(ii,:,:,:) = month_to_season(prec(ii,:,:,:),"DJF")
  prec_seas(ii,:,:,:) = month_to_annual(prec(ii,:,:,:),1)

;  prec_djf(ii,:,:,:) = prec_seas(ii,1:nyr_all-2,:,:)
  prec_djf(ii,:,:,:) = prec_seas(ii,:,:,:)

  prec_areaave1(ii,:) = wgt_areaave(prec_djf(ii,:,{latmin1:latmax1},{lonmin:lonmax}), 1, 1.0, 0) 
  prec_areaave2(ii,:) = wgt_areaave(prec_djf(ii,:,{latmin2:latmax2},{lonmin:lonmax}), 1, 1.0, 0) 

  prec_asym(ii,:) = prec_areaave2(ii,:) - prec_areaave1(ii,:)

; calculate trend over period of interest

;  prec_trend(ii) = regCoef(xtime1, prec_asym(ii,0:105))
  prec_trend(ii) = regCoef(xtime2, prec_asym(ii,106:200))

 end do ; end loop over models

;************************************************
; set up variables for scatter plot
;************************************************

  tts1 = sst_trend*100.
  tts2 = prec_trend*100.

  ccr = escorc(tts1,tts2)
  prob  = rtest(ccr, nmod, 0)

;  print( sprintf("%6.2f",ccr) )

;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("eps","scatter")      

  res             = True 
  
  res@gsnMaximize = True

  res@vpHeightF= 0.5                           
  res@vpWidthF = 0.55 
  
  res@tiMainFontHeightF = 0.018
  res@tiXAxisFontHeightF = 0.018
  res@tiYAxisFontHeightF = 0.018

  res@gsnLeftStringFontHeightF = 0.018
  res@gsnRightStringFontHeightF = 0.018

  rsq = ccr * ccr

  rc  = regline(tts1,tts2)

  res@gsnLeftString = "(d) RCP8.5 (n=27)"
  res@gsnRightString = "r="+sprintf("%5.2f",ccr)+", p="+sprintf("%5.2f",prob)
;  res@gsnRightString = "r~S~2~N~="+sprintf("%5.2f",rsq)

  res@tiYAxisString = "Pr trend (mm/day/century)"  
  res@tiXAxisString = "SST trend (K/century)"
 
  res@trXMinF =  -0.5
  res@trXMaxF =   0.5

  res@trXMinF =   0
  res@trXMaxF =   8.0

  res@trYMinF =  -1.0
  res@trYMaxF =   1.0

  res@gsnYRefLine           = 0.0              
  res@gsnYRefLineDashPattern = 1

  res@gsnXRefLine           = 0.0                
  res@gsnXRefLineDashPattern = 1

;  res@xyMarkLineMode      = "Markers"             
 
  res@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
  res@xyDashPatterns      = 1                      ; solid line 
  res@xyLineThicknesses   = (/1,2/)                ; set second line to 2

  res@xyMarkers            =  16                    
  res@xyMarkerColor        =  "black"
  res@xyMarkerSizeF       = 0.001  
  
  res@gsnFrame               = False  
  res@gsnDraw                = False  
   
  data_line = new ( (/2,dimsizes(tts2)/), typeof(tts2))
  data_line(0,:) = tts2
  data_line(1,:) = rc*(tts1-rc@xave) + rc@yave

  plot  = gsn_csm_xy (wks,tts1,tts2,res)  
;  plot  = gsn_csm_xy (wks,tts1,data_line,res)        ; create plot

 ; plot DRY and WET

  polyres                     = True   

 ; polyres@gsMarkerThicknessF  = 2.5

  polyres@gsMarkerIndex       = 16         
  polyres@gsMarkerSizeF       = 0.02  
     
  polyres@gsMarkerColor       = "blue"

  dum1 = gsn_add_polymarker(wks,plot,tts1(modnum_lo),tts2(modnum_lo),polyres)   

  polyres@gsMarkerColor       = "red"

  dum2 = gsn_add_polymarker(wks,plot,tts1(modnum_hi),tts2(modnum_hi),polyres) 

  draw(plot)
  frame(wks)
     
end
