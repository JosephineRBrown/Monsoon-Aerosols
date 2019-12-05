;**********************************************
; Calculate differencein TAS for two periods in 
; HISTAA runs and subtract ocean 30S-30N
; Following Chung and Soden (2017) Figure 1b.
;
; MONSOON_AA FIGURE 1A
;
; JRB 14.12.17
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

 latmin = -30.0
 latmax = 30.0

 lonmin = 0.
 lonmax = 360. 

 temp = new((/nmod,ntime,nlat,nlon/),float)

 temp1 = new((/nmod,ntime1,nlat,nlon/),float)
 temp2 = new((/nmod,ntime1,nlat,nlon/),float)

 temp_ave1 = new((/nmod,nlat,nlon/),float)
 temp_ave2 = new((/nmod,nlat,nlon/),float)

 temp_dy= new((/nmod,nlat,nlon/),float) 
 temp_dy_land= new((/nmod,nlat,nlon/),float)
 slopemean = new((/nlat,nlon/),float)

 tropdiff = new((/nmod/),float)

 path = "/g/data/eg3/jxb548/CMIP5/historicalAA/TAS/PROCNEW/"

 do ii = 0, nmod-1
 print(nmodel(ii))
 
; a1 = addfile("/g/data/eg3/jxb548/ECP45/rm_pr_mon_"+model(ii)+"_2006-2300.nc","r") 

 if (ii .eq. 0) then ; CanESM2
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p4_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 1) then ; CCSM4
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p10_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 2) then ; CSIRO
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p4_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 3) then ; FGOALS
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r2i1p1_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 4) then ; GFDL-CM3
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p1_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 5) then; GFDL-ESM2M
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p5_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 6) then ; GISS-E2-H
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p107_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 7) then ; GISS-E2-R
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p107_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 8) then ; IPSL-CM5A-LR
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p3_1900-2005_r240x120.nc","r") 
 end if
 if (ii .eq. 9) then ; NorESM1-M
  a1 = addfile(path+"tas_Amon_"+model(ii)+"_historicalMisc_r1i1p1_1900-2005_r240x120.nc","r") 
 end if

 lat = a1->lat
 lon = a1->lon

 temp(ii,:,:,:) = a1->tas(:,:,:)
 temp(ii,:,:,:) = temp(ii,:,:,:) - 273.15

 temp1(ii,:,:,:) = temp(ii,0:371,:,:)     ; 1900-1930
 temp2(ii,:,:,:) = temp(ii,900:1271,:,:)  ; 1975-2005

  temp_ann1 = month_to_annual(temp1(ii,:,:,:),1)
  temp_ann2 = month_to_annual(temp2(ii,:,:,:),1)

;  temp_ann1 = month_to_season( temp1(ii,:,:,:), "JJA" )
;  temp_ann2 = month_to_season( temp2(ii,:,:,:), "JJA" )

; average over time

  temp_ave1(ii,:,:) = dim_avg_n_Wrap(temp_ann1,0) 
  temp_ave2(ii,:,:) = dim_avg_n_Wrap(temp_ann2,0)

 end do ; end loop over models
 
 temp_dy!0 = "modvar"
 temp_dy&modvar = ispan(1,nmod,1)
 temp_dy!1 = "lat"
 temp_dy&lat = lat
 temp_dy!2 = "lon"
 temp_dy&lon = lon

 maskfile = addfile("/g/data/eg3/jxb548/CMIP5/z7000.nc","r")
 mskdats = maskfile->data(0,:,:)
 mskdat = short2flt(mskdats)

 temp_dy_land = temp_dy

 do ii = 0, nmod-1

  temp_dy(ii,:,:) = temp_ave2(ii,:,:) - temp_ave1(ii,:,:)

  temp_dy_land(ii,:,:) = mask(temp_dy(ii,:,:),mskdat,0)

  tropdiff(ii) = wgt_areaave(temp_dy_land(ii,{-30.:30.},:),1.0,1.0,0)
   
; print(nmodel(ii))  
 print(tropdiff(ii))

   temp_dy(ii,:,:) = temp_dy(ii,:,:) - tropdiff(ii)

 end do ; end loop over models

 slopemean(:,:) = dim_avg_Wrap( temp_dy(lat|:,lon|:,modvar|:) )   ; average over models

;************************************************
;set up colors
;*************************************************
  wks = gsn_open_wks("oldeps","climo")            
  gsn_define_colormap(wks,"CBR_coldhot")     
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

  res@cnMinLevelValF  =  -0.5			; set the minimum contour level
  res@cnMaxLevelValF  =   0.5			; set the maximum contour level
  res@cnLevelSpacingF =   0.05			; set the interval between contours 

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

  res@gsnLeftString = "(a) HistoricalAA temperature difference"
  res@gsnRightString = "K" 

  plot = gsn_csm_contour_map(wks,slopemean(:,:),res)  
 

 draw (plot)
 frame(wks)

end
