; identify list of local maxima in a 2-d map
; optionally choose threshold, and a minimum distance to reject many
; sources that are close to eachother
;
; In:  map     = 2d array
;      thresh  = threshold in same units as map
;      mindist = minimum distance between sources
;      header  = FITS header for astrometry
;
; Out: x,y,z   = x, y coordinate and brightness z
;      ra, dec = decimal degrees if header specified

pro findpeaks, map, x, y, z, mindist=mindist, thresh=thresh, header=header, $
               ra=ra, dec=dec

  xpix = n_elements(map[*,0])
  ypix = n_elements(map[0,*])

  if keyword_set(header) then begin
    cdelt1 = sxpar(header,"CDELT1")
    crpix1 = sxpar(header,"CRPIX1")
    crval1 = sxpar(header,"CRVAL1")
    cdelt2 = sxpar(header,"CDELT2")
    crpix2 = sxpar(header,"CRPIX2")
    crval2 = sxpar(header,"CRVAL2")

    if( cdelt1 eq 0 ) then begin
      cdelt1 = sxpar(header,"CD1_1")
      cdelt2 = sxpar(header,"CD2_2")
    endif

    pixres = abs( cdelt1 )*3600.

  endif

  ; x- and y- offsets to all of the pixels in the perimeter of a 3x3 box
  ; centered over the pixel we're examining
  per_xoff = [-1, 0, 1, 1, 1, 0,-1,-1]
  per_yoff = [ 1, 1, 1, 0,-1,-1,-1, 0]

  nsource = 0l

  ; loop over all pixels except for the edges

  srcflag = bytarr(xpix,ypix)

  for i=1l, xpix-2 do begin
    for j=1l, ypix-2 do begin
      
      ; identify indices for the perimeter in the map centered over this
      ; pixel
      per_xind = i+per_xoff
      per_yind = j+per_yoff
      perimeter = map[per_xind, per_yind]

      ; the pixel is a source if it is the local maximum
      if( map[i,j] gt max(perimeter)) then begin
        srcflag[i,j] = 1b
      endif
    endfor
  endfor

  srcind = where(srcflag eq 1b, nsource)
  if srcind[0] ne -1 then begin
    x = srcind mod xpix
    y = srcind / xpix
  endif

  z = map[x,y]

  ; threshold the sources
  if keyword_set(thresh) then begin
    ind = where(z ge thresh)
    if ind[0] ne -1 then begin
      x = x[ind]
      y = y[ind]
      z = z[ind]
    endif
  endif

  ; multiple source rejection

  if keyword_set(mindist) then begin
    goodflag = intarr(n_elements(x))+1

    ; loop over all sources that we found
    for i=0, n_elements(x)-1 do if(goodflag[i]) then begin

      ; calculate dist^2 to all other sources from this source
      d = (x[i]-x)^2d + (y[i]-y)^2d

      ; flag all sources that are closer than the minimum distance, that
      ; haven't been rejected already
      ind = where( (d le mindist^2d) and (goodflag) )

      if ind[0] ne -1 then begin
        ; if we found some other sources nearby, reject all but the brightest
        badind = where( z[ind] lt max(z[ind]) )

        ;if n_elements(badind) gt 1 then badind = badin[1:n_elements(badins)-1]
        if badind[0] ne -1 then goodflag[ind[badind]] = 0
      endif
    endif
    
    ; update the list to contain only the good ones
    ind = where(goodflag)
    if ind[0] ne -1 then begin
      x = x[ind]
      y = y[ind]
      z = z[ind]
    endif
  endif

  ; calculate RA and Dec for each source

  if keyword_set(header) then begin
    ra_off = (x - (crpix1-1))*cdelt1 
    dec_off = (y - (crpix2-1))*cdelt2 

    tan2radec, crval1, crval2, ra_off, dec_off, ra, dec
  endif
      
end
