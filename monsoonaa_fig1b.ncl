;**********************************************
; Calculate trend in PR for ECP4.5 and plot MMM
;
; MONSOON_AA FIGURE 1B
;
; JRB 19.06.17
;***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

 begin

 model = (/ "CanESM2","CCSM4","CSIRO-Mk3-6-0","FGOALS-g2","GFDL-CM3" ,\
            "GFDL-ESM2M","GISS-E2-H","GISS-E2-R","IPSL-CM5A-LR","NorESM1-M"/)

 nmodel = model

 nmod = 10
 
 nlat = 120
 nlon = 240

 nyr = 106
 ntime = nyr*12
 nyr1 = 31
 ntime1 = nyr1*12

 prec = new((/nmod,ntime,nlat,nlon/),float)

 prec1 = new((/nmod,ntime1,nlat,nlon/),float)
 prec2 = new((/nmod,ntime1,nlat,nlon/),float)

 prec_ave1 = new((/nmod,nlat,nlon/),float)
 prec_ave2 = new((/nmod,nlat,nlon/),float)

 prec_dy= new((/nmod,nlat,nlon/),float)
 slopemean = new((/nlat,nlon/),float)

 path = "/g/data/eg3/jxb548/CMIP5/historicalAA/PR/PROCNEW/"

 do ii = 0, nmod-1
 print(nmodel(ii))
 
; a1 = addfile("/g/data/eg3/jxb548/ECP45/rm_pr_mon_"+model(ii)+"_2006-2300.nc","r") 

 if (ii .eq. 0) then ; CanESM2
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p4_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 1) then ; CCSM4
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p10_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 2) then ; CSIRO
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p4_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 3) then ; FGOALS
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r2i1p1_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 4) then ; GFDL-CM3
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p1_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 5) then; GFDL-ESM2M
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p5_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 6) then ; GISS-E2-H
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p107_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 7) then ; GISS-E2-R
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p107_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 8) then ; IPSL-CM5A-LR
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p3_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 9) then ; NorESM1-M
  a1 = addfile(path+"pr_Amon_"+model(ii)+"_historicalMisc_r1i1p1_1900-2005_r240x120.nc","r") 
 end if


 prec(ii,:,:,:) = a1->pr(:,:,:)
 prec(ii,:,:,:) = 86400*prec(ii,:,:,:)

 prec1(ii,:,:,:) = prec(ii,0:371,:,:)     ; 1900-1930
 prec2(ii,:,:,:) = prec(ii,900:1271,:,:)  ; 1975-2005

 lat = a1->lat
 lon = a1->lon

  prec_ann1 = month_to_annual(prec1(ii,:,:,:),1)
  prec_ann2 = month_to_annual(prec2(ii,:,:,:),1)

;  prec_ann1 = month_to_season( prec1(ii,:,:,:), "DJF" )
;  prec_ann2 = month_to_season( prec2(ii,:,:,:), "DJF" )

; average over time

  prec_ave1(ii,:,:) = dim_avg_n_Wrap(prec_ann1,0) 
  prec_ave2(ii,:,:) = dim_avg_n_Wrap(prec_ann2,0)

  prec_dy(ii,:,:) = prec_ave2(ii,:,:) - prec_ave1(ii,:,:)

 end do ; end loop over models

 prec_dy!0 = "modvar"
 prec_dy&modvar = ispan(1,nmod,1)
 prec_dy!1 = "lat"
 prec_dy&lat = lat
 prec_dy!2 = "lon"
 prec_dy&lon = lon

 slopemean(:,:) = dim_avg_Wrap( prec_dy(lat|:,lon|:,modvar|:) )   ; average over models

; calculate fraction of model agreement for stippling

 pct = new((/nlat,nlon/),float)
 nct = new((/nlat,nlon/),float)
 cct = new((/nlat,nlon/),float)
 cct2 = new((/nlat,nlon/),float)
 
 pct(:,:) = 0
 nct(:,:) = 0
 cct(:,:) = 0
 cct2(:,:) = 0

 do ii = 0, nmod - 1
  do i = 0, nlat - 1 
   do j = 0, nlon - 1

   if (.not. ismissing(prec_dy(ii,i,j)) .and. prec_dy(ii,i,j) .gt. 0.0 ) then
    pct(i,j) = pct(i,j) + 1
   end if
   if (.not. ismissing(prec_dy(ii,i,j)) .and. prec_dy(ii,i,j) .lt. 0.0 ) then
    nct(i,j) = nct(i,j) + 1
   end if
 
   end do
  end do
 end do 

 do i = 0, nlat - 1 
 do j = 0, nlon - 1
  if ( pct(i,j) .gt. 0.70*nmod .or. nct(i,j) .gt. 0.70*nmod ) then ; for stippling
   cct(i,j) = 1
  end if
 end do
 end do

 copy_VarAtts(slopemean,cct)                    
 copy_VarCoords(slopemean,cct) 

;************************************************
;set up colors
;*************************************************
  wks = gsn_open_wks("oldeps","climo")            
  gsn_define_colormap(wks,"CBR_drywet")     
;************************************************
; create panel plots
;*************************************************

  nplot = 1

  plot   = new (nplot , graphic)                    ; create graphical array
  
  res                      = True               ; plot options desired
  
;  res@gsnMaximize = True
 
;  res@cnFillMode           = "RasterFill"

  res@gsnSpreadColors      = True    		; use full range of colors
;  res@gsnSpreadColorStart  = -1	  	        ; to reverse color table
;  res@gsnSpreadColorEnd    =  2     		; to reverse color table
  res@cnFillOn             = True               ; turn on color fill
  res@cnInfoLabelOn        = False              ; turn off contour info label
  res@cnLinesOn            = False              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off line labels
 
  res@cnLevelSelectionMode = "ManualLevels"    ; set explicit contour levels

  res@cnMinLevelValF  =  -0.4			; set the minimum contour level
  res@cnMaxLevelValF  =   0.4			; set the maximum contour level
  res@cnLevelSpacingF =   0.04			; set the interval between contours 

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

  res@mpFillOn             = False               ; turn off gray continents

  res@gsnStringFontHeightF = 0.026

  res@gsnLeftStringFontHeightF = 0.020
  res@gsnRightStringFontHeightF = 0.020

  res@gsnLeftString = "(b) HistoricalAA precipitation difference"
  res@gsnRightString = "mm/day" 

  plot = gsn_csm_contour_map(wks,slopemean(:,:),res)  
 
 ;;---- overplot stippling where X% of models agree on sign of change
 res2 = True

 res2@gsnDraw              = False
 res2@gsnFrame             = False

 res2@gsnLeftString = "" 
 res2@gsnRightString = ""

 res2@gsnShadeFillType = "pattern"      ; pattern fill
 res2@gsnShadeHigh = 17
 res2@cnFillScaleF = 1.0
 res2@cnFillDotSizeF = 0.003
 
 res2@cnLinesOn           = False       ; do not draw contour lines
 res2@cnLineLabelsOn      = False       ; do not draw contour labels
 res2@cnInfoLabelOn        = False      ; turn off contour info label

 plot2 = gsn_csm_contour(wks,cct(:,:),res2)

 plot2 = gsn_contour_shade(plot2,-999,0,res2)        ; Pattern shade contours above 

 overlay(plot, plot2)

 draw (plot)
 frame(wks)

end
