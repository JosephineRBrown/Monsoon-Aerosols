;**********************************************
; Calculate trend in PR for ECP4.5 and plot MMM
;
; JRB 19.06.17
;***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

 begin

; read in model files

  modall = (/ "ACCESS1-0","ACCESS1-3","bcc-csm1-1","BNU-ESM","CCSM4",\
              "CESM1-WACCM","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FGOALS-s2",\
              "FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","HadGEM2-ES",\
              "MIROC5","MPI-ESM-LR","MRI-CGCM3","NorESM1-M"/)

 ; mod_loforc = (/ "bcc-csm1-1","bcc-csm1-1-m","BNU-ESM","CCSM4","CESM1-BGC",\
 ;                 "CESM1-WACCM","EC-EARTH","FGOALS-s2","FIO-ESM",\
 ;                 "GFDL-ESM2G","GFDL-ESM2M","MPI-ESM-LR","MPI-ESM-MR"/)  

 ; mod_hiforc = (/ "ACCESS1-0","ACCESS1-3","CESM1-CAM5","CSIRO-Mk3-6-0",\
 ;                 "FGOALS-g2","GFDL-CM3","HadGEM2-CC","HadGEM2-ES",\
 ;                 "MIROC-ESM","MIROC-ESM-CHEM","MIROC5","MRI-CGCM3",\
 ;                 "NorESM1-M","NorESM1-ME" /)

; for RCP2.6

  mod_loforc = (/ "bcc-csm1-1","bcc-csm1-1-m","BNU-ESM","CCSM4",\
                  "FGOALS-s2","FIO-ESM",\
                  "GFDL-ESM2G","GFDL-ESM2M","MPI-ESM-LR","MPI-ESM-MR"/)  

  mod_hiforc = (/ "CESM1-CAM5","CSIRO-Mk3-6-0",\
                  "FGOALS-g2","GFDL-CM3","HadGEM2-ES",\
                  "MIROC-ESM","MIROC-ESM-CHEM","MIROC5","MRI-CGCM3",\
                  "NorESM1-M","NorESM1-ME" /)

 

; model1 = modall
 model1 = mod_loforc
 model2 = mod_hiforc

; nmodel = model
            
 nmod1 = 10  ; LoForce = 10 or 13
 nmod2 = 11   ; HiForce = 11 or 14

 nlat = 120
 nlon = 240

 nyr = 30
 ntime = nyr*12

 lprec = new((/nmod1,ntime,nlat,nlon/),float)

 lprec1 = new((/nmod1,ntime,nlat,nlon/),float)
 lprec2 = new((/nmod1,ntime,nlat,nlon/),float)

 lprec_ave1 = new((/nmod1,nlat,nlon/),float)
 lprec_ave2 = new((/nmod1,nlat,nlon/),float)

 lprec_diff = new((/nmod1,nlat,nlon/),float)
 lprec_diffmean = new((/nlat,nlon/),float)

 hprec = new((/nmod2,ntime,nlat,nlon/),float)

 hprec1 = new((/nmod2,ntime,nlat,nlon/),float)
 hprec2 = new((/nmod2,ntime,nlat,nlon/),float)

 hprec_ave1 = new((/nmod2,nlat,nlon/),float)
 hprec_ave2 = new((/nmod2,nlat,nlon/),float)

 hprec_diff = new((/nmod2,nlat,nlon/),float)
 hprec_diffmean = new((/nlat,nlon/),float) 

 prec_diffmean = new((/nlat,nlon/),float)

 path1 = "/g/data1/eg3/pfd548/IPCC/CMIP5/historical/pr_processed/grid1.5"
 path2 = "/g/data1/eg3/pfd548/IPCC/CMIP5/rcp26/pr_processed/grid1.5"

 do ii = 0,nmod1-1 ; loop over models

 print(model1(ii))

  a1 = addfile(path1+"/pr_"+model1(ii)+"_r1i1p1_1900_2005_r240x120.nc","r") 
  lprc1 = a1->pr(840:1199,:,:)  ; 1970-1999
  lprc1 = lprc1*24*60*60

;  lprc2 = a1->pr(0:371,:,:)  ; 1900-1930
;  lprc2 = lprc2*24*60*60

; below is for LoForc RCP4.5 and RCP8.5 only 
;  if (ii .ne. 5 .and. ii .ne. 6) then
;   a2 = addfile(path2+"/pr_"+model1(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 
;   end if
;  if (ii .eq. 5 .or. ii .eq. 6) then
;   a2 = addfile(path2+"/pr_"+model1(ii)+"_r2i1p1_2006_2100_r240x120.nc","r") 
;  end if
; below is for LoForc RCP2.6
  a2 = addfile(path2+"/pr_"+model1(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 

  lprc2 = a2->pr(768:1127,:,:)  ; 2070-2099
  lprc2 = lprc2*24*60*60

  lprec1(ii,:,:,:) = lprc1
  lprec2(ii,:,:,:) = lprc2

  lat = a1->lat
  lon = a1->lon

  lprec_ann1 = month_to_annual(lprec1(ii,:,:,:),1)
  lprec_ann2 = month_to_annual(lprec2(ii,:,:,:),1)

;  lprec_ann1 = month_to_season(lprec1(ii,:,:,:), "DJF" )
;  lprec_ann2 = month_to_season(lprec2(ii,:,:,:), "DJF" )

; average over time

  lprec_ave1(ii,:,:) = dim_avg_n_Wrap(lprec_ann1,0) 
  lprec_ave2(ii,:,:) = dim_avg_n_Wrap(lprec_ann2,0)

  lprec_diff(ii,:,:) =  lprec_ave2(ii,:,:) - lprec_ave1(ii,:,:)
 ; lprec_diff(ii,:,:) =  lprec_ave1(ii,:,:) - lprec_ave2(ii,:,:)

 end do

 lprec_ave1!0 = "modvar"
 lprec_ave1&modvar = ispan(1,nmod1,1)
 lprec_ave1!1 = "lat"
 lprec_ave1&lat = lat
 lprec_ave1!2 = "lon"
 lprec_ave1&lon = lon

 lprec_diff!0 = "modvar"
 lprec_diff&modvar = ispan(1,nmod1,1)
 lprec_diff!1 = "lat"
 lprec_diff&lat = lat
 lprec_diff!2 = "lon"
 lprec_diff&lon = lon

; lprec_diffmean(:,:) = dim_avg_Wrap( lprec_ave1(lat|:,lon|:,modvar|:) )   ; average over models
 lprec_diffmean(:,:) = dim_avg_Wrap( lprec_diff(lat|:,lon|:,modvar|:) )   ; average over models

 do ii = 0,nmod2-1 ; loop over models

 print(model2(ii))

  a1 = addfile(path1+"/pr_"+model2(ii)+"_r1i1p1_1900_2005_r240x120.nc","r") 
  hprc1 = a1->pr(840:1199,:,:)  ; 1970-1999
  hprc1 = hprc1*24*60*60

 ; hprc2 = a1->pr(0:371,:,:)     ; 1900-1930
 ; hprc2 = hprc2*24*60*60

; below is for HiForc
  a2 = addfile(path2+"/pr_"+model2(ii)+"_r1i1p1_2006_2100_r240x120.nc","r") 

  hprc2 = a2->pr(768:1127,:,:)  ; 2070-2099 
  hprc2 = hprc2*24*60*60

  hprec1(ii,:,:,:) = hprc1
  hprec2(ii,:,:,:) = hprc2
 
  hprec_ann1 = month_to_annual(hprec1(ii,:,:,:),1)
  hprec_ann2 = month_to_annual(hprec2(ii,:,:,:),1)

;  hprec_ann1 = month_to_season(hprec1(ii,:,:,:), "DJF" )
;  hprec_ann2 = month_to_season(hprec2(ii,:,:,:), "DJF" )

; average over time

  hprec_ave1(ii,:,:) = dim_avg_n_Wrap(hprec_ann1,0) 
  hprec_ave2(ii,:,:) = dim_avg_n_Wrap(hprec_ann2,0)

  hprec_diff(ii,:,:) =  hprec_ave2(ii,:,:) - hprec_ave1(ii,:,:)
 ; hprec_diff(ii,:,:) =  hprec_ave1(ii,:,:) -  hprec_ave2(ii,:,:)

end do

 hprec_ave1!0 = "modvar"
 hprec_ave1&modvar = ispan(1,nmod2,1)
 hprec_ave1!1 = "lat"
 hprec_ave1&lat = lat
 hprec_ave1!2 = "lon"
 hprec_ave1&lon = lon

 hprec_diff!0 = "modvar"
 hprec_diff&modvar = ispan(1,nmod2,1)
 hprec_diff!1 = "lat"
 hprec_diff&lat = lat
 hprec_diff!2 = "lon"
 hprec_diff&lon = lon

; hprec_diffmean(:,:) = dim_avg_Wrap( hprec_ave1(lat|:,lon|:,modvar|:) )   ; average over models
 hprec_diffmean(:,:) = dim_avg_Wrap( hprec_diff(lat|:,lon|:,modvar|:) )   ; average over models

 prec_diffmean = lprec_diffmean
 prec_diffmean(:,:) = hprec_diffmean(:,:) - lprec_diffmean(:,:)

;************************************************
;set up colors
;*************************************************
  wks = gsn_open_wks("oldeps","climo")            
  gsn_define_colormap(wks,"CBR_drywet")  
 ; gsn_define_colormap(wks,"BlueDarkRed18")       
;************************************************
; create panel plots
;*************************************************

  nplot = 1

  plot   = new (nplot , graphic)                    ; create graphical array
  
  res                      = True               ; plot options desired
  
;  res@gsnMaximize = True
 
;  res@cnFillMode           = "RasterFill"

  res@gsnSpreadColors      = True    		; use full range of colors
 ; res@gsnSpreadColorStart  = -1	  	        ; to reverse color table
 ; res@gsnSpreadColorEnd    =  2     		; to reverse color table
  res@cnFillOn             = True               ; turn on color fill
  res@cnInfoLabelOn        = False              ; turn off contour info label
  res@cnLinesOn            = False              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off line labels
 
  res@cnLevelSelectionMode = "ManualLevels"    ; set explicit contour levels

  res@cnMinLevelValF  =  -1.0			; set the minimum contour level
  res@cnMaxLevelValF  =   1.0			; set the maximum contour level
  res@cnLevelSpacingF =   0.1			; set the interval between contours 

  res@mpFillOn             = True            
  res@mpCenterLonF         = 180                
  
  res@mpMinLonF            =    0.           
  res@mpMaxLonF            =  360.
  res@mpMinLatF            =  -50. 
  res@mpMaxLatF            =   50.
 
  res@txFontHeightF   	   = 0.024          
  res@tiMainFontHeightF    = 0.024             
  res@tiMainString = ""

  res@tmXBMajorLengthF	   = 0.014	       
  res@tmXBMinorLengthF 	   = 0.007
  res@tmYLMajorLengthF 	   = 0.014
  res@tmYLMinorLengthF 	   = 0.007

  res@tmXBLabelFontHeightF = 0.018	       
  res@tmYLLabelFontHeightF = 0.018

  res@lbLabelBarOn         = True

 ; res@lbOrientation        = "vertical" 
  res@lbLabelFontHeightF = 0.016
 ; res@lbLabelStride        = 2
  res@pmLabelBarOrthogonalPosF = .30           ; move whole thing down

  res@gsnDraw              = False
  res@gsnFrame             = False

  res@mpFillOn             = False           

  res@gsnStringFontHeightF = 0.026

  res@gsnLeftStringFontHeightF = 0.020
  res@gsnRightStringFontHeightF = 0.020

  res@gsnLeftString = "(a) RCP2.6 - Historical (HiForce - LoForce)"
  res@gsnRightString = "mm/day" 

  plot = gsn_csm_contour_map(wks,prec_diffmean(:,:),res)  
 
 draw (plot)
 frame(wks)

end
