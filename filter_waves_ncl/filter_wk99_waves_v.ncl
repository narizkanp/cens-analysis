;=====================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
; User Contributed routines follow. They are NOT supported by the NCL team. 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/kf_filter.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
;http://www.atmos.albany.edu/student/ventrice/documents/filter_rt.ncl
;scipted by 31Okt2018
;=====================================================================

begin
;-----------------reading data------------------------
  inFile = addfile("vrm3harmean.nc","r")
  olr = short2flt(inFile->vrm3harmean(:,{-25:25},:))  ;[TimexLatxLon] 1979-2019
  printVarSummary(olr)
  dimensi = dimsizes(olr)
  nt = dimensi(0)
  nlat = dimensi(1)
  nlon = dimensi(2)
  lat = inFile->lat({-25:25})
;-------------------------------------------------------
;------define symmetri & anti-symmetric compenents------
  olrAS = olr
  olrAS  = decompose2SymAsym( olr , 0 ) 
  olrsymm = olrAS
  olrasym = olrAS
  printVarSummary(olrsymm)
  printVarSummary(olrasym)
  do yy = 0,nlat/2-1
      olrsymm(:,20-yy,:) = olrAS(:,yy,:)
      olrasym(:,yy,:) = -olrAS(:,20-yy,:)
  end do
  olrasym(:,{0},:) = 0
;--------------------------------------------------------
;--------set the wave parameters--------------
; KELVIN
  tMinK = 2.5
  tMaxK = 30
  kMinK = 1
  kMaxK = 14
  hMinK = 8
  hMaxK = 90
  waveNameK = "Kelvin"

; EQUATORIAL ROSSBY
  tMinER = 9.7
  tMaxER = 48
  kMinER = -10
  kMaxER = -1
  hMinER = 8
  hMaxER = 90
  waveNameER = "ER"


; MRG
  tMinMRG = 3
  tMaxMRG = 10
  kMinMRG = -10
  kMaxMRG = -1
  hMinMRG = 8
  hMaxMRG = 90
  waveNameMRG = "MRG"


; MJO
  tMinMJO = 30
  tMaxMJO = 96
  kMinMJO = 1
  kMaxMJO = 5
  waveNameMJO = "MJO"
;--------------------------------------------

;---create new variable for filtered data---
  kelvin = olr
  kelvin@wave_name = "Kelvin Wave"
  kelvin@filter = "Wheeler & Kiladis (1999)"
  kelvin@wavenumber = (/ kMinK, kMaxK /)
  kelvin@period = (/ tMinK, tMaxK /)

  er = olr
  er@wave_name = "Equatorial Rossby Wave"
  er@filter = "Wheeler & Kiladis (1999)"
  er@wavenumber = (/ kMinER, kMaxER /)
  er@period = (/ tMinER, tMaxER /)
  
  mrg = olr
  mrg@wave_name = "Mixed-Rossby Gravity"
  mrg@filter = "Wheeler & Kiladis (1999)"
  mrg@wavenumber = (/ kMinMRG, kMaxMRG /)
  mrg@period = (/ tMinMRG, tMaxMRG /)

  mjo = olr
  mjo@wave_name = "Madden-Julian Oscillation"
  mjo@filter = "Wheeler & Kiladis (1999)"
  mjo@wavenumber = (/ kMinMJO, kMaxMJO /)
  mjo@period = (/ tMinMJO, tMaxMJO /)

;-------------------------------------------
;-------------------filtering proccess-----------
; olr : original data
; olrsymm : symmteric compenent
; olrasym : anti-symmetric component
  mis = -999
  mis@_FillValue = -999
  obsPerDay = 1

  print( "Filtering..." )
  do y = 0,nlat-1
    print( y )
    kelvin(time|:,lat|y,lon|:) = kf_filter( olr(time|:,lat|y,lon|:), obsPerDay \
                      , tMinK, tMaxK, kMinK, kMaxK, hMinK, hMaxK, waveNameK )

    er(time|:,lat|y,lon|:) = kf_filter( olr(time|:,lat|y,lon|:), obsPerDay \
                      , tMinER, tMaxER, kMinER, kMaxER, hMinER, hMaxER, waveNameER )

    mrg(time|:,lat|y,lon|:) = kf_filter( olr(time|:,lat|y,lon|:), obsPerDay \
                      , tMinMRG, tMaxMRG, kMinMRG, kMaxMRG, hMinMRG, hMaxMRG, waveNameMRG )

    mjo(time|:,lat|y,lon|:) = kf_filter( olr(time|:,lat|y,lon|:), obsPerDay \
                      , tMinMJO, tMaxMJO, kMinMJO, kMaxMJO, mis, mis, "None" )
  end do

;-----------------------------------------------------
;-------------------saving data to .nc--------------

  f2 = addfile("wk_kelvin_v.nc","c")
  f2->v_kelvin = kelvin

  f3 = addfile("wk_er_v.nc","c")
  f3->v_er = er

  f4 = addfile("wk_mrg_v.nc","c")
  f4->v_mrg = mrg

  f5 = addfile("wk_mjo_v.nc","c")
  f5->v_mjo = mjo

;---------------------------------------------------

end
