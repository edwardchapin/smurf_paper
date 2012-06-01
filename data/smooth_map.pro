;------------------------------------------------------------------------------
; smooth_map.pro:
;
;  This function takes a pixel map of non-smoothed fluxes (synthetic survey),
; and smoothes it with a specified gaussian beam (normalized, so spatially
; integrated flux is conserved).
;
; At present, the image is convolved with a gaussian or a disc
;------------------------------------------------------------------------------

function smooth_map, fluxmap, $    ; The unsmoothed bitmap (y=dec,x=ra)
                pixres, $          ; pixel resolution (arcsec/pixel)
                bmwd, $            ; beamwidth (arcsec)
                beammap=beammap, $ ; return the beam in map of same size
                sum=sum, $         ; makes Gaussian peak=1
                squared=squared, $ ; smooth by the beam^2
                areanorm=areanorm,$; if set, normalize to 1/sr
                disc=disc          ; if set, smooth with disc of diameter
                                   ; bdpix, and height 1 (i.e. add up
                                   ; all the pixels in a circular aperture)

  ; allocate space

  xpix = n_elements( fluxmap(*,0) )
  ypix = n_elements( fluxmap(0,*) )
  bdpix = float(bmwd) / float(pixres)

  beam = fltarr(xpix, ypix)
  image = fltarr(xpix, ypix)

  ; create the beam profile

  min_dim = min( [round(bdpix*10 < xpix),round(bdpix*10 < ypix)] )
  if (min_dim le 1) and (xpix ge 3) and (ypix ge 3) then min_dim=3


  xodd = xpix mod 2
  yodd = ypix mod 2

  beamprofile = shift( dist(min_dim), min_dim/2 + xodd, min_dim/2 + yodd )

  ; Disc convolution kernel for summing over apertures

  if keyword_set(disc) then begin
    mask = where(beamprofile le bdpix/2.)
    beamprofile(*) = 0.
    beamprofile(mask) = 1.0
  endif $

  ; Gaussian Beam

  else beamprofile = exp( -beamprofile^2/(2 * (bdpix/2.35)^2) )

  ; insert the profile into a map of the same size as the image so that
  ; we may convolve the two.

  beam(xpix/2 - min_dim/2 : xpix/2 + min_dim/2 - (1 - (min_dim mod 2)), $
       ypix/2 - min_dim/2 : ypix/2 + min_dim/2 - (1 - (min_dim mod 2)) ) = $
       beamprofile

  if not( keyword_set(disc) )  and not(keyword_set(sum))then $
    beam = beam/total(beam)  ; Normalize the beam

  if keyword_set(areanorm) then begin
    ; normalize such that integral over sr subtended by each pixel eq 1
    pixarea = ((pixres/3600d)*!DPI/180d)^2d
    beam = beam/total(beam*pixarea)
  endif

  ; square the beam
  if keyword_set(squared) then beam = beam^2.

  beam = shift(beam,xpix/2,ypix/2)

  ; now convolve the input map with the beam

  image = fft( fft(fluxmap,1) * fft(beam,1), -1 )

  beammap = shift(beam,xpix/2,ypix/2)

  return, image

end













