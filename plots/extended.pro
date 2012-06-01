;------------------------------------------------------------------------------
; M17 extended data reduction example
;------------------------------------------------------------------------------

datadir = '../data/'

black = 0
red = 1
green = 2
blue = 3
yellow = 4
grey = 5
orange=9
purple=10
white=11
lightblue=12

; default reduction showing ripples at different iteration numbers

iter = ['02', '30', '02', '84']
labels = ["!6(a) default", "(b) default", "(c) bright extended", $
          "(d) bright extended"]

badval=-1000

fxread, datadir+'m17_resp_fplane_align.fits', resp, header
bad = where( finite(resp) eq 0, complement=good )
if bad[0] ne -1 then resp[bad] = 0
if good[0] ne -1 then resp[good] = 1


nx = n_elements(resp[*,0])
ny = n_elements(resp[0,*])

pixres = abs(sxpar( header, "CDELT1" ))

aspect = double(ny)/double(nx)


map_iter = dblarr( nx, ny, 4 )
for i=0, 1 do begin
  fxread, datadir+'m17_default_'+iter[i]+'.fits', map, header
  bad = where( finite(map) eq 0 )
  if bad[0] ne -1 then map[bad] = badval
  map_iter[*,*,i] = map
endfor






set_plot, 'ps'

!p.thick=3.
!x.thick = !p.thick
!y.thick = !p.thick

device, filename="m17.eps", /encapsulated, bits_per_pixel=8, $
        xsize=20, ysize=20*aspect, /color

loadct,0

; clip ranges for intensity maps
mn = -0.0003
mx = 0.01

thick =3.

cs = 1.0
ys = 0.027
xm = 0.00


xl = [0.0, 0.5, 0.0, 0.5]
yl = [0.5, 0.5, 0.0, 0.0]


for i=0, 1 do begin
  pos = [xl[i], yl[i], xl[i]+0.5, yl[i]+0.5]

  ; log-scaled intensity map
  map = map_iter[*,*,i]
  map = map < mx
  map = ((map - mn) > 0)-mn
  map = alog10(map)

  tvscale, -map, /noint, pos=pos
  xyouts, pos[0]+xm, pos[3]-ys, labels[i], /normal, charsize=cs, $
    charthick=thick
  xyouts, pos[0]+xm, pos[3]-2*ys, "n="+strcompress(fix(iter[i]),/remove_all), $
    /normal, charsize=cs, charthick=thick


  if i eq 0 then begin
    ; array footprint
    resp = shift(resp,-250,+25)
    contour, resp, xmargin=[0,0], ymargin=[0,0], $
      xstyle=5, ystyle=5, pos=pos, /noerase, $
      levels=0.5, c_thick=thick*3., c_color=255

    contour, resp, xmargin=[0,0], ymargin=[0,0], $
      xstyle=5, ystyle=5, pos=pos, /noerase, $
      levels=0.5, c_thick=thick, c_color=0

    ; filter scale
    len = (300./3600.)/pixres
    len_norm = (len/double(nx)) * (pos[2]-pos[0])
    xcen = (pos[2]+pos[0])/2. - 0.17
    ycen = (pos[3]+pos[1])/2. - 0.07
    plots, xcen+len_norm*[-0.5,0.5], ycen*[1., 1.], /normal, thick=thick*3, $
      color=255
    plots, xcen+len_norm*[-0.5,0.5], ycen*[1., 1.], /normal, thick=thick, $
      color=0
    xyouts, xcen, ycen-0.02, '300 arcsec', /normal, charsize=cs, $
      charthick=thick*3., align=0.5, color=255
    xyouts, xcen, ycen-0.02, '300 arcsec', /normal, charsize=cs, $
      charthick=thick, align=0.5, color=0
  endif
endfor



; --- M17 bright_extended ---

fname = datadir+'m17_bright_extended.fits'
fxread, fname, map_bright_final, header
bad = where( finite(map_bright_final) eq 0 )
if bad[0] ne -1 then map_bright_final[bad] = badval
fxread, fname, qual_final, header, extension=2

fname = datadir+'m17_bright_extended_2.fits'
fxread, fname, map_bright_2, header
bad = where( finite(map_bright_2) eq 0 )
if bad[0] ne -1 then map_bright_2[bad] = badval
fxread, fname, qual_2, header, extension=2


for i=2, 3 do begin
  pos = [xl[i], yl[i], xl[i]+0.5, yl[i]+0.5]

  if i eq 2 then begin
    map = map_bright_2 < mx
    qual = qual_2
  endif else begin
    map = map_bright_final < mx
    qual = qual_final
  endelse

  map = ((map - mn) > 0)-mn
  map = alog10(map)

  tvscale, -map, /noint, pos=pos
  xyouts, pos[0]+xm, pos[3]-ys, labels[i], /normal, charsize=cs, $
    charthick=thick
  xyouts, pos[0]+xm, pos[3]-2*ys, "n="+strcompress(fix(iter[i]),/remove_all), $
    /normal, charsize=cs, charthick=thick

  mycolour

  contour, qual, xmargin=[0,0], ymargin=[0,0], $
    xstyle=5, ystyle=5, pos=pos, /noerase, $
    levels=0.5, c_thick=1., c_color=red

  loadct,0

endfor


device, /close



set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1

end
