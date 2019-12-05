;***********************************************
; Calculate temp in ausm region 
; from CESM large ensemble 
; MONSOON_AA FIGURE 10a
; JRB 04.12.18
;***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

 begin

; read in model files

; nmod = 33
; nrun = (/"001","002","003","004","005","006","007","008","009","010",\ 
;          "011","012","013","014","015","016","017","018","019","020",\
;          "021","022","023","024","025","026","027","028","029","030",\
;          "031","032","033","034" /)

 nmod = 15
 nrun = (/"016","017","018","019","020",\
          "021","022","023","024","025","026","027","028","029","030" /)

 nlat = 192
 nlon = 288

; ave domains

 latmin1 =  -20.
 latmax1 =    0.
 
 latmin2 =    0.
 latmax2 =   20.

 lonmin =   0.  
 lonmax = 360.

 nyr1 = 86  ; 1920-2005
 ntime1 = nyr1*12
 nyr2 = 75  ; 2006-2080
 ntime2 = nyr2*12
 nyr3 = 20  ; 2081-2100
 ntime3 = nyr3*12
  
 nyr = nyr1 + nyr2 + nyr3
 ntime = nyr*12

 temp = new((/nmod,ntime,nlat,nlon/),float)
 temp_seas = new((/nmod,nyr,nlat,nlon/),float)

 temp_djf = new((/nmod,nyr-1,nlat,nlon/),float) 
 temp_djf_land = new((/nmod,nyr-1,nlat,nlon/),float)
 areaave = new((/nmod,nyr-1/),float)
 areaave1 = new((/nmod,nyr-1/),float)
 areaave2 = new((/nmod,nyr-1/),float)
 ts_asym = new((/nmod,nyr-1/),float)

 path1 = "/g/data/eg3/jxb548/CESM/HIST/"
 path2 = "/g/data/eg3/jxb548/CESM/RCP85/"
 path4 = "/g/data1/eg3/afm599/lih-legacy/afmDat/afm/p37/landMask0p3/"

 maskfile = addfile(path4+"sftlf_CESM1-CAM5_landMask.nc","r")
 mskdat = maskfile->sftlf(0,:,:)

 do ii = 0,nmod-1 ; loop over models

  print(nrun(ii))
  
  a1 = addfile(path1+"b.e11.B20TRC5CNBDRD.f09_g16."+nrun(ii)+".cam.h0.TS.192001-200512.nc","r") 
  ts1 = a1->TS(0:ntime1-1,:,:)
  ts1 = ts1 - 273.15

  a2 = addfile(path2+"b.e11.BRCP85C5CNBDRD.f09_g16."+nrun(ii)+".cam.h0.TS.200601-208012.nc.nc","r") 
  ts2 = a2->TS(0:ntime2-1,:,:)
  ts2 = ts2 - 273.15

  a3 = addfile(path2+"b.e11.BRCP85C5CNBDRD.f09_g16."+nrun(ii)+".cam.h0.TS.208101-210012.nc.nc","r") 
  ts3 = a3->TS(0:ntime3-1,:,:)
  ts3 = ts3 - 273.15

  temp(ii,0:ntime1-1,:,:) = ts1
  temp(ii,ntime1:ntime1+ntime2-1,:,:) = ts2
  temp(ii,ntime1+ntime2:ntime1+ntime2+ntime3-1,:,:) = ts3

 ; temp_seas(ii,:,:,:) = month_to_season(temp(ii,:,:,:),"DJF")
  temp_seas(ii,:,:,:) = month_to_annual(temp(ii,:,:,:),1)
  temp_djf(ii,:,:,:) = temp_seas(ii,1:nyr-1,:,:)

  temp_djf_land(ii,:,:,:) = temp_djf(ii,:,:,:)
  temp_djf_land(ii,:,:,:) = mask(temp_djf(ii,:,:,:),mskdat,0)  ; 0 is ocean

; select land or all below:
  
  areaave1(ii,:) = wgt_areaave(temp_djf(ii,:,{latmin1:latmax1},{lonmin:lonmax}), 1, 1.0, 0) 
  areaave2(ii,:) = wgt_areaave(temp_djf(ii,:,{latmin2:latmax2},{lonmin:lonmax}), 1, 1.0, 0) 

  ts_asym(ii,:) = areaave2(ii,:) - areaave1(ii,:)
   
 areaave(ii,:) =  ts_asym(ii,:) ; borrow name to save time

 end do ; end loop over models

; now we have annual area ave value for each model, need to plot these

  x   = ispan(0,nyr-2,1) 
  xf = x + 1920 + 1

; calculate anomalies relative to 1960-1990 period

  do ii = 0, nmod-1 
   areaave(ii,:) =  areaave(ii,:) - dim_avg(areaave(ii,40:69))
   areaave(ii,:) = runave_Wrap(areaave(ii,:),11,0)
  end do

; set up variable for plotting

  areaave!0    = "runs"
  areaave!1    = "time"
  areaave&time = xf
  areaave&runs = ispan(0,nmod-1,1) 
  
  printVarSummary(areaave)  

  mnmx      = new ((/3,nyr-1/), float)
;  mnmx(0,:) = dim_min( areaave(time|:,runs|:) )
;  mnmx(1,:) = dim_max( areaave(time|:,runs|:) )
  mnmx(0,:) = dim_avg( areaave(time|:,runs|:) ) - dim_stddev( areaave(time|:,runs|:) )
  mnmx(1,:) = dim_avg( areaave(time|:,runs|:) ) + dim_stddev( areaave(time|:,runs|:) )
  mnmx(2,:) = dim_avg( areaave(time|:,runs|:) )    

;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("oldeps","time") 
 gsn_define_colormap(wks,"cosam12") 

 res = True                                    

 res@gsnDraw            = False             ; don't draw yet
 res@gsnFrame           = False             ; don't advance frame yet

 res@tiMainString = "(a) CESM RCP8.5 Standard: SST asymmetry"
 res@txFontHeightF = 0.020                    
 res@tiMainFontHeightF = 0.020

 res@trXMinF = 1920.0
 res@trXMaxF = 2100.0
; res@trYMinF = 0
; res@trYMaxF = 16
 res@trYMinF = -0.3
 res@trYMaxF =  0.3

 res@vpHeightF= 0.3                             ; change aspect ratio of plot
 res@vpWidthF = 0.8 
 res@vpXF     = 0.15                            ; start plot at x ndc coord

 res@gsnYRefLine           = 0.0         
 res@gsnYRefLineDashPattern = 1

; res@gsnAboveYRefLineColor = "red"              ; above ref line fill red
; res@gsnBelowYRefLineColor = "blue"             ; below ref line fill blue                 
 res@tiYAxisString = "Temperature (K)"
; res@tiXAxisString = "Year"

; res@tmXBMode = "Manual"
; res@tmXBTickStartF = 1920
; res@tmXBTickSpacingF = 20
; res@tmXBTickEndF = 2006

 res@tmXBLabelFontHeightF = 0.018
 res@tmYLLabelFontHeightF = 0.018
 res@tiYAxisFontHeightF = 0.018
 res@tiXAxisFontHeightF = 0.018 
 res@xyDashPattern = 0. 

; res@xyLineThicknessF = 1.
; bot_plot  = gsn_csm_xy (wks,xf,areaave(:,:),res)   ; create plot
 
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
