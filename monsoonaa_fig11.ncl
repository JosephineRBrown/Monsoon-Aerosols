;**********************************************
; Calculate trend in PR for ECP4.5 and plot MMM
; MONSOON_AA FIGURE 11
;
; JRB 19.06.17
;***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

 begin

; read in model files

 nmod = 15
 nrun = (/"016","017","018","019","020","021","022",\
          "023","024","025","026","027","028","029","030" /)

 nlat = 192
 nlon = 288

 nyr = 30
 ntime = nyr*12

 nyr1 = 86  ; 1920-2005
 ntime1 = nyr1*12
 nyr2 = 75  ; 2006-2080
 ntime2 = nyr2*12
 nyr3 = 20  ; 2081-2100
 ntime3 = nyr3*12
  
 nyr_all = nyr1 + nyr2 + nyr3
 ntime_all = nyr*12

 re   = 6.37122e06
 rad  = 4.0 * atan(1.0) / 180.0
 con  = re * rad

 latmin = -90.
 latmax = 90.  

 prec = new((/nmod,ntime,nlat,nlon/),float)

 prec1 = new((/nmod,ntime,nlat,nlon/),float)
 prec2 = new((/nmod,ntime,nlat,nlon/),float)
 prec3 = new((/nmod,ntime,nlat,nlon/),float)
 prec3_all = new((/nmod,ntime2+ntime3,nlat,nlon/),float)

 prec_ave1 = new((/nmod,nlat,nlon/),float)
 prec_ave2 = new((/nmod,nlat,nlon/),float)
 prec_ave3 = new((/nmod,nlat,nlon/),float)

 prec_diff = new((/nmod,nlat,nlon/),float)
 prec_diffmean = new((/nlat,nlon/),float)

 path1 = "/g/data/eg3/jxb548/CESM/HIST/"
 path2 = "/g/data/eg3/jxb548/CESM/RCP85/"
 path3 = "/g/data/eg3/jxb548/CESM/FXARCP85/"

 do ii = 0,nmod-1 ; loop over models

; print(model(ii))

  a1 = addfile(path1+"b.e11.B20TRC5CNBDRD.f09_g16."+nrun(ii)+".cam.h0.PRECTOT.192001-200512.nc","r") 
  prc1 = a1->PRECC(600:959,:,:) ; 1970-1999
  prc1 = prc1*24*60*60*1000

  a2 = addfile(path3+"b.e11.BRCP85C5CNBDRD.f09_g16."+nrun(ii)+".b.cam.h0.PRECTOT.200601-210012.nc.nc","r") 
  prc2 = a2->PRECC(768:1127,:,:) ; 2070-2099
  prc2 = prc2*24*60*60*1000

  lat = a1->lat
  lon = a1->lon

  prec1(ii,:,:,:) = prc1
  prec2(ii,:,:,:) = prc2

  a3 = addfile(path2+"b.e11.BRCP85C5CNBDRD.f09_g16."+nrun(ii)+".cam.h0.PRECTOT.200601-208012.nc.nc","r") 
  prc3 = a3->PRECC(0:ntime2-1,:,:)
  prc3 = prc3*24*60*60*1000

  a4 = addfile(path2+"b.e11.BRCP85C5CNBDRD.f09_g16."+nrun(ii)+".cam.h0.PRECTOT.208101-210012.nc.nc","r") 
  prc4 = a4->PRECC(0:ntime3-1,:,:)
  prc4 = prc4*24*60*60*1000

  prec3_all(ii,0:ntime2-1,:,:) = prc3
  prec3_all(ii,ntime2:ntime2+ntime3-1,:,:) = prc4

  prec3(ii,:,:,:) = prec3_all(ii,768:1127,:,:)

  prec_ann1 = month_to_annual(prec1(ii,:,:,:),1)
  prec_ann2 = month_to_annual(prec2(ii,:,:,:),1)
  prec_ann3 = month_to_annual(prec3(ii,:,:,:),1)

;  prec_ann1 = month_to_season( prec1(ii,:,:,:), "DJF" )
;  prec_ann2 = month_to_season( prec2(ii,:,:,:), "DJF" )

; average over time

  prec_ave1(ii,:,:) = dim_avg_n_Wrap(prec_ann1,0) 
  prec_ave2(ii,:,:) = dim_avg_n_Wrap(prec_ann2,0)
  prec_ave3(ii,:,:) = dim_avg_n_Wrap(prec_ann3,0)

  prec_diff(ii,:,:) =  prec_ave3(ii,:,:) - prec_ave2(ii,:,:)

 end do ; end loop over models

 prec_diff!0 = "modvar"
 prec_diff&modvar = ispan(1,nmod,1)
 prec_diff!1 = "lat"
 prec_diff&lat = lat
 prec_diff!2 = "lon"
 prec_diff&lon = lon

 prec_diffmean(:,:) = dim_avg_Wrap( prec_diff(lat|:,lon|:,modvar|:) )   ; average over models

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

   if (.not. ismissing(prec_diff(ii,i,j)) .and. prec_diff(ii,i,j) .gt. 0.0 ) then
    pct(i,j) = pct(i,j) + 1
   end if
   if (.not. ismissing(prec_diff(ii,i,j)) .and. prec_diff(ii,i,j) .lt. 0.0 ) then
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

 copy_VarAtts(prec_diffmean,cct)                    
 copy_VarCoords(prec_diffmean,cct) 

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
;  res@gsnSpreadColorStart  = -1	  	        ; to reverse color table
;  res@gsnSpreadColorEnd    =  2     		; to reverse color table
  res@cnFillOn             = True               ; turn on color fill
  res@cnInfoLabelOn        = False              ; turn off contour info label
  res@cnLinesOn            = False              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off line labels
 
  res@cnLevelSelectionMode = "ManualLevels"    ; set explicit contour levels

  res@cnMinLevelValF  =   -0.5			; set the minimum contour level
  res@cnMaxLevelValF  =    0.5			; set the maximum contour level
  res@cnLevelSpacingF =    0.1			; set the interval between contours 

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

  res@gsnLeftString = "Standard RCP8.5 - FixedAA RCP8.5"
  res@gsnRightString = "mm/day" 

  plot = gsn_csm_contour_map(wks,prec_diffmean(:,:),res)  
 
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
