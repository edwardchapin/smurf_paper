; convert between spherical and pixel coordaintes

pro cat_pix, lon, $           ; Lon. in degrees
             lat, $           ; Lat. in degrees
             x, $             ; x pixel coordinate
             y, $             ; y pixel coordinate
             header, $        ; FITS header for transformation
             pix2cat=pix2cat,$; if set, convert x,y -> lon,lat, otherwise does
                              ; the inverse by default
             pixres=pixres, $ ; return projection information if requested
             cdelt1=cdelt1, $
             cdelt2=cdelt2, $
             crpix1=crpix1, $
             crpix2=crpix2, $
             crval1=crval1, $
             crval2=crval2

  ; Get projection from the header
  CDELT1 = sxpar(header,"CD1_1")
  CDELT2 = sxpar(header,"CD2_2")
  if CDELT1 eq 0 then begin
    CDELT1 = sxpar(header,"CDELT1")
    CDELT2 = sxpar(header,"CDELT2")
  endif
  CRPIX1 = sxpar(header,"CRPIX1")
  CRPIX2 = sxpar(header,"CRPIX2")
  CRVAL1 = sxpar(header,"CRVAL1")
  CRVAL2 = sxpar(header,"CRVAL2")

  if keyword_set(pix2cat) then begin
    ; convert x,y --> lon, lat
    lon_off  = (x - (CRPIX1-1))*CDELT1
    lat_off = (y - (CRPIX2-1))*CDELT2
    tan2radec, CRVAL1, CRVAL2, lon_off, lat_off, lon, lat
  endif else begin
    ; convert lon, lat --> x,y
    radec2tan, CRVAL1, CRVAL2, lon, lat, lon_off, lat_off
    x = lon_off/CDELT1 + CRPIX1 - 1
    y = lat_off/CDELT2 + CRPIX2 - 1
  endelse

  pixres = abs(CDELT1)

end
