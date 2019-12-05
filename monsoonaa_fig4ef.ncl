;***********************************************
; Plot prec asymmetry in RCP8.5
; For NESP aerosol work
; MONSOON_AA FIGURE 4e,f
;
; JRB 29.03.18
;***********************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

 begin

; read in model files

 mod_loforc_all = (/ "bcc-csm1-1","bcc-csm1-1-m","BNU-ESM","CCSM4","CESM1-BGC",\
                     "CESM1-FASTCHEM","CESM1-WACCM","EC-EARTH","FGOALS-s2","FIO-ESM",\
                     "GFDL-ESM2G","GFDL-ESM2M","MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P"/)

  mod_loforc_rcp85 = (/ "bcc-csm1-1","bcc-csm1-1-m","BNU-ESM","CCSM4","CESM1-BGC",\
                       "CESM1-WACCM","EC-EARTH","FGOALS-s2","FIO-ESM",\
                       "GFDL-ESM2G","GFDL-ESM2M","MPI-ESM-LR","MPI-ESM-MR"/)

 mod_hiforc_all = (/ "ACCESS1-0","ACCESS1-3","CESM1-CAM5","CESM1-CAM5-1-FV2","CSIRO-Mk3-6-0",\
                     "FGOALS-g2","GFDL-CM3","HadGEM2-AO","HadGEM2-CC","HadGEM2-ES",\
                     "MIROC-ESM","MIROC-ESM-CHEM","MIROC5","MRI-CGCM3","MRI-ESM1",\
                     "NorESM1-M","NorESM1-ME" /)

;; exclude HadGEM2-AO due to strange pr fields.

 
 mod_hiforc_rcp85 = (/ "ACCESS1-0","ACCESS1-3","CESM1-CAM5","CSIRO-Mk3-6-0",\
                     "FGOALS-g2","GFDL-CM3","HadGEM2-CC","HadGEM2-ES",\
                     "MIROC-ESM","MIROC-ESM-CHEM","MIROC5","MRI-CGCM3",\
                     "NorESM1-M","NorESM1-ME" /)

 mod_all_all = (/ "ACCESS1-0","ACCESS1-3","bcc-csm1-1","bcc-csm1-1-m","BNU-ESM",\
                 "CCSM4","CESM1-BGC","CESM1-CAM5","CESM1-CAM5-1-FV2","CESM1-FASTCHEM",\
                 "CESM1-WACCM","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FGOALS-s2",\
                 "FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","HadGEM2-AO",\
                 "HadGEM2-CC","HadGEM2-ES","MIROC-ESM","MIROC-ESM-CHEM","MIROC5",\
                 "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P","MRI-CGCM3","MRI-ESM1",\
                 "NorESM1-M","NorESM1-ME"/)

 model = mod_loforc_rcp85
 nmodel = model
 
 nmod = 13 ; 13 or 14 RCP8.5

 nlat = 120
 nlon = 240

; NH and SH domains 

 latmin1 =  -20.
 latmax1 =    0.
 
 latmin2 =    0.
 latmax2 =   20.

 lonmin =   0.  
 lonmax = 360.

 nyr = 106  ; 1900-2005
 ntime = nyr*12
 
 prcseas1 = new((/nmod,nyr,nlat,nlon/),float)

 nyr2 = 95  ; 2006-2100
 ntime2 = nyr2*12

 prcseas2 = new((/nmod,nyr2,nlat,nlon/),float)

; combined variables

 nyr_all = nyr + nyr2
 ntime_all = nyr_all*12

 precip = new((/nmod,ntime_all,nlat,nlon/),float)
 precip_seas = new((/nmod,nyr_all,nlat,nlon/),float)

 precip_djf = new((/nmod,nyr_all-2,nlat,nlon/),float) 
 precip_djf_land = new((/nmod,nyr_all-2,nlat,nlon/),float)
 areaave1 = new((/nmod,nyr_all-2/),float)
 areaave2 = new((/nmod,nyr_all-2/),float)

 pr_asym = new((/nmod,nyr_all-2/),float)

 path1 = "/g/data1/eg3/pfd548/IPCC/CMIP5/historical/pr_processed/grid1.5"
 path2 = "/g/data1/eg3/pfd548/IPCC/CMIP5/rcp85/pr_processed/grid1.5"

 maskfile = addfile("/g/data/eg3/jxb548/CMIP5/z7000.nc","r")
 mskdats = maskfile->data(0,:,:)
 mskdat = short2flt(mskdats)

 do ii = 0,nmod-1 ; loop over models

 print(model(ii))

  a1 = addfile(path1+"/pr_"+model(ii)+"_r1i1p1_1900_2005_r240x120.nc","r") 
  prc1 = a1->pr(0:ntime-1,:,:)
  prc1 = prc1*24*60*60

; for loforce
  if (ii .ne. 5 .and. ii .ne. 6) then
   a2 = addfile(path2+"/pr_"+model(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 
   end if
  if (ii .eq. 5 .or. ii .eq. 6) then
   a2 = addfile(path2+"/pr_"+model(ii)+"_r2i1p1_2006_2100_r240x120.nc","r") 
  end if

; for hiforce
;  a2 = addfile(path2+"/pr_"+model(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 

  prc2 = a2->pr(0:ntime2-1,:,:)
  prc2 = prc2*24*60*60

  ntime_end2 = ntime + ntime2 - 1

  precip(ii,0:ntime-1,:,:) = prc1
  precip(ii,ntime:ntime_end2,:,:) = prc2

 ; precip_seas(ii,:,:,:) = month_to_season(precip(ii,:,:,:),"DJF")
  precip_seas(ii,:,:,:) = month_to_annual(precip(ii,:,:,:),1)

  precip_djf(ii,:,:,:) = precip_seas(ii,1:nyr_all-2,:,:)

  precip_djf_land(ii,:,:,:) = precip_djf(ii,:,:,:)
  precip_djf_land(ii,:,:,:) = mask(precip_djf(ii,:,:,:),mskdat,1)  ; 0 is ocean

; select LAND or LAND+OCEAN below:

  areaave1(ii,:) = wgt_areaave(precip_djf(ii,:,{latmin1:latmax1},{lonmin:lonmax}), 1, 1.0, 0) 
  areaave2(ii,:) = wgt_areaave(precip_djf(ii,:,{latmin2:latmax2},{lonmin:lonmax}), 1, 1.0, 0) 

  pr_asym(ii,:) = areaave2(ii,:) - areaave1(ii,:)

 end do ; end loop over models

; now we have annual area ave value for each model, need to plot these

  x   = ispan(0,nyr_all-3,1) 
  xf = x + 1900 + 1
 
 ; calculate anomalies relative to 1986-2005 period ??

  do ii = 0, nmod-1 
   pr_asym(ii,:) =  pr_asym(ii,:) - dim_avg(pr_asym(ii,85:104))
   pr_asym(ii,:) = runave_Wrap(pr_asym(ii,:),11,0)
  end do

; set up variable for plotting

  pr_asym!0    = "models"
  pr_asym!1    = "time"
  pr_asym&time = xf
  pr_asym&models = ispan(0,nmod-1,1) 
  
  printVarSummary(pr_asym)  

  mnmx      = new ((/3,nyr_all-2/), float)
;  mnmx(0,:) = dim_min( pr_asym(time|:,models|:) )
;  mnmx(1,:) = dim_max( pr_asym(time|:,models|:) )
  mnmx(0,:) = dim_avg( pr_asym(time|:,models|:) ) - dim_stddev( pr_asym(time|:,models|:) )
  mnmx(1,:) = dim_avg( pr_asym(time|:,models|:) ) + dim_stddev( pr_asym(time|:,models|:) )
  mnmx(2,:) = dim_avg( pr_asym(time|:,models|:) )    

;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("oldeps","time") 
 gsn_define_colormap(wks,"cosam12") 

 res = True                                    

 res@gsnDraw            = False             ; don't draw yet
 res@gsnFrame           = False             ; don't advance frame yet

 res@tiMainString = "(e) RCP8.5 LoForce: Precipitation Asymmetry"
 res@txFontHeightF = 0.020                     
 res@tiMainFontHeightF = 0.020

 res@trXMinF = 1900.0
 res@trXMaxF = 2100.0
; res@trYMinF = 0
; res@trYMaxF = 16
 res@trYMinF = -0.5
 res@trYMaxF =  0.5

 res@vpHeightF= 0.3                             ; change aspect ratio of plot
 res@vpWidthF = 0.8
 res@vpXF     = 0.15                            ; start plot at x ndc coord

 res@gsnYRefLine           = 0.0         
 res@gsnYRefLineDashPattern = 1

; res@gsnAboveYRefLineColor = "red"              ; above ref line fill red
; res@gsnBelowYRefLineColor = "blue"             ; below ref line fill blue                 
 res@tiYAxisString = "Precipitation (mm/day)"  
; res@tiXAxisString = "Year"

 res@tmXBMode = "Manual"
 res@tmXBTickStartF = 1900
 res@tmXBTickSpacingF = 20
 res@tmXBTickEndF = 2100

 res@tmXBLabelFontHeightF = 0.018
 res@tmYLLabelFontHeightF = 0.018
 res@tiYAxisFontHeightF = 0.018
 res@tiXAxisFontHeightF = 0.018 
 res@xyDashPattern = 0. 

; res@xyLineThicknessF = 1.
; bot_plot  = gsn_csm_xy (wks,xf,pr_asym(:,:),res)   ; create plot
 
 res@xyLineColor = "Black"
 res@xyLineThicknessF = 4.
 top_plot = gsn_csm_xy (wks,xf,mnmx(2,:),res)       ; create line plot

; Create a plot with the area between both curves filled in pink.

 res@gsnXYFillColors = "gray80"
 res@xyLineColor     = -1                        

 mid_plot  = gsn_csm_xy (wks,xf,mnmx(0:1,:),res)    ; create plot 

 overlay(top_plot,mid_plot) 
;  overlay(top_plot,bot_plot) 

 pres = True
 maximize_output(wks,pres)

end
