;------------------------------------------------------------------------------
; Simulate some cirrus to add to an M17 map
;------------------------------------------------------------------------------

showplot = 0

; load in the map to which we will be adding the fake signal
; to.

fxread, "m17_bright_extended.fits", data, header
fxread, "m17_bright_extended.fits", qual, qheader, extension=2
bad = where( finite(data) eq 0 )
if bad[0] ne -1 then data[bad]=0

pixres = abs( sxpar(header,"CDELT1") )*3600.
nx = n_elements( data[*,0] )
ny = n_elements( data[0,*] )

;if showplot then begin
;  window,0
;  tvscale, data, /keep, /noint, maxval=0.002, minval=-0.001
;endif


; create some cirrus

openr,1,"cirrus_position"
xl=0 & yl=0 & n=0
readf,1,xl,yl,n
close,1

c = cirrus( n, pixres/60., 1, -3, 100, seed=3 ) ; scale is meaningless
c = c - min(c)
c = 0.002 * c / stdev(c)
c_sm = smooth_map( c, pixres, 14.5 )

r = shift(dist(n),n/2,n/2)
apod = r*0
apod = (1+cos((r-75)*!DPI/25.))/2d
ind = where(r lt 75)
apod[ind]=1
ind = where(r ge 100)
apod[ind]=0

c_sm_apod = c_sm * apod

;if showplot then begin
;  window,1
;  tvscale, c_sm_apod, /keep, /noint, maxval=0.002, minval=-0.001
;endif

;smooth_data = smooth_map(data,pixres,14.5)

fakemap = data*0
fakemap[xl:xl+n-1, yl:yl+n-1] = c_sm_apod

data = data + fakemap

;if showplot then begin
;  window,2
;  tvscale, data, /keep, /noint, maxval=0.002, minval=-0.001
;endif


mn = -0.001
mx = 0.01

map = data
map = map < mx
map = ((map - mn) > 0)-mn
map = alog10(map)

;if showplot then begin
;  window,3
;  tvscale, -map, /noint, /keep, pos=pos
;endif

fxwrite, "m17_fakemap.fits", header, fakemap

exit

end
