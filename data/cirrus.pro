;------------------------------------------------------------------------------
; Definition of cirrus power spectrum: (from Gautier et al. 1992)
;
;  P(sigma) = P_0(sigma/sigma_0)^alpha
;
;  where sigma is the spatial frequency in 1/arcminutes
;
; Inputs:
;
;  pixels   = # pixels on side in the cirrus map
;  pixres   = pixel resolution in arcmin/pixel
;  B_0      = mean surface brightness in Jy/sr at 100 um
;  alpha    = index of the radial spectrum fall-off (about -3.0)
;  lambda   = wavelength of observation (optically thin greybody at
;             18k used to extrapolate) in microns 
;
; Returned map is in Jy/sr  
;
; (note: IRAS had 0.57m primary)
;
;------------------------------------------------------------------------------

function cirrus, pixels, pixres, B_0, alpha, lambda, seed=seed

  PI = 3.14159

  ; Calculate power at scale 0.01 arcmin^-1 from surface brightness

  P_0     = 1.4e-12*B_0^3 ; Power at scale 0.01 1/arcmin in Jy^2/sr

  ; Wavelength and cirrus spectra parameters (to extrapolate flux)

  temp = 18.     ; temperature of dust in Kelvin

  scale_flux = (100./lambda)^2 * $
               black(3.e8/(lambda*1.e-6), temp) / black(3.e8/100.e-6, temp)

  ; calculate the spatial frequency at each point in FFT of power
  ; spectrum

  sfreq = findgen(sqrt(2)*pixels)/(pixres*pixels*sqrt(2)) 

  ; Scale P_0 to be in (Jy*arcmin/sr)^2

  SP_0 = P_0 * 3438.^2    ; 3438 arcmin / rad

  ; Calculate the stepsize of each pixel in spatial frequency

  d_sfreq = (1/pixres)/(pixels*sqrt(2))

  ; Generate the filter

  filt = sqrt(SP_0) * $
         ((dist(pixels)/(pixres*pixels*sqrt(2)))/0.01)^(alpha/2.) * $
         d_sfreq

  filt(0) = B_0

  ; Filter some noise

  ;phase = randomu(seed,pixels,pixels)*2.*3.141592
  ;spec = filt * complex(cos(phase),sin(phase))

  fnoise = fft(randomn(seed,pixels,pixels))
  pnoise = sqrt( float(fnoise)^2 + imaginary(fnoise)^2 )
  fnoise = fnoise/mean(pnoise)
  spec = filt * fnoise

  spec(0) = B_0  ; the mean surface brightness level

  cirrus = real( fft(spec, 1) )

  ; Scale the cirrus to the brightness at another wavelength

  cirrus = cirrus * scale_flux

  ;window,0
  ;tvscale,cirrus,/nointerpolation,/keep_aspect_ratio

  return, cirrus

end
