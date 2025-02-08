;===========================================================
;THIS PROGRAM CALCULATES THE FIRST THREE HARMONICS OF ANNUAL
;CYCLE ALL YEAR
;Coded by : faiz@ganteng
;===========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin
;---------------------------user's specification-------------------
 fa  = addfile("u925edit.nc","r")
 yST   = 1979      ; year start
 yET   = 2021      ; year end
 LAT   = 80        ; used in plotting
 LON   = 95        ; used in plotting
 nhar  = 3         ; the first nhar harmonics of annual cycle
;------------------------------------------------------------------
 ymdST = yST*10000+1*100+1
 ymdET = yET*10000+12*100+31
 time = fa->time
 ymd = cd_calendar(time,-2)
 iST = toint(ind(ymd.eq.ymdST))
 iET = toint(ind(ymd.eq.ymdET))
 da  = short2flt(fa->u(iST:iET,:,:))
 printVarSummary(da)
 
 nt = dimsizes (da&time)
 nlat = dimsizes (da&lat) 
 nlon = dimsizes (da&lon)
 nY = yET-yST+1
 h = new ((nhar), integer)
 do i = 0,nhar-1,1
   h(i)   = nY*(i+1)-1     ; assuming that number of day in a year is 366
 end do

 db = da
 copy_VarCoords(da,db)
 yrata = dim_avg_n(db,0)
 copy_VarCoords(db(0,:,:),yrata(:,:))
 finfo_db = fourier_info(db(lat|:,lon|:,time|:),0,1.0)
 copy_VarCoords(db(0,:,:),finfo_db(0,:,:,0))
 finfo_db!0 = "ftype"
 finfo_db!3 = "harmonics" 
 y0 = new ((/nt,nlat,nlon/), float)
 y1 = new ((/nt,nlat,nlon/), float)
 y2 = new ((/nt,nlat,nlon/), float)
 yr = new ((/nt,nlat,nlon/), float)
 yt = new ((/nt,nlat,nlon/), float)
 pi = 4*atan(1.0)
 tpN = 2.0*pi/nt
 phi0 = finfo_db(1,:,:,h(0))
 phi1 = finfo_db(1,:,:,h(1))
 phi2 = finfo_db(1,:,:,h(2))
 A0 = finfo_db(0,:,:,h(0))
 A1 = finfo_db(0,:,:,h(1))
 A2 = finfo_db(0,:,:,h(2))
 printVarSummary(y0)
 printVarSummary(A0)
 do t1 = 0,nt-1,1
    y0(t1,:,:) = A0(:,:)*cos(tpN*1*nY*(t1-phi0(:,:)))
    y1(t1,:,:) = A1(:,:)*cos(tpN*2*nY*(t1-phi1(:,:)))
    y2(t1,:,:) = A2(:,:)*cos(tpN*3*nY*(t1-phi2(:,:)))
    yr(t1,:,:) = yrata(:,:)
	print(t1)
 end do   
 
 
 yt(:,:,:) = yr(:,:,:)+y0(:,:,:)+y1(:,:,:)+y2(:,:,:)
 ym = da-yt
 copy_VarCoords(da,yt)   ; yt = the 1st 3 har of seasonal cycle + mean  of data 
 copy_VarCoords(da,yr)   ; yr = mean of data
 copy_VarCoords(da,y0)   ; y0 = the 1st har of seasonal cycle of data
 copy_VarCoords(da,y1)   ; y1 = the 2nd har of seasonal cycle of data   
 copy_VarCoords(da,y2)   ; y2 = the 3rd har of seasonal cycle of data
 copy_VarCoords(da,ym)   ; ym = data - yt

;-------------------saving files to .nc ------------------------
; f1=addfile("olr.3harmean.nc","c")
; f1->olr3harmean = yt

; f2=addfile("olr.mean.nc","c")
; f2->olrmean = yr

; f3=addfile("olr.har1.nc","c")
; f3->olrhar1 = y0

; f4=addfile("olr.har2.nc","c")
; f4->olrhar2 = y1

; f5=addfile("olr.har3.nc","c")
; f5->olrhar3 = y2

 f6=addfile("urm3harmean.nc","c")
 f6->urm3harmean = ym

end