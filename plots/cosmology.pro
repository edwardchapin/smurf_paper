;------------------------------------------------------------------------------
; Plots for lockman data reduction
;------------------------------------------------------------------------------

datadir = '../data/'

jacksize = 110 ; needs to match value in cosmology script

; load in the data

badval=-1000
badvar=1000

fxread, datadir+'lockman_s8d_blank_field.fits', rawmap, mapheader
bad = where( finite(rawmap) eq 0 )
if bad[0] ne -1 then rawmap[bad] = badval

fxread, datadir+'lockman_s8d_blank_field_whitened.fits', whitemap, mapheader
bad = where( finite(whitemap) eq 0 )
if bad[0] ne -1 then whitemap[bad] = badval

fxread, datadir+'lockman_s8d_blank_field.fits', rawvar, header, extension=1
bad = where( finite(rawvar) eq 0 )
if bad[0] ne -1 then rawvar[bad] = badvar
rawerr = sqrt(rawvar)

fxread, datadir+'jackknife.fits', jkmap, header
bad = where( finite(jkmap) eq 0 )
if bad[0] ne -1 then jkmap[bad] = badval

fxread, datadir+'jackknife.fits', jkvar, header, extension=1
bad = where( finite(jkvar) eq 0 )
if bad[0] ne -1 then jkvar[bad] = badvar
jkerr = sqrt(jkvar)

fxread, datadir+'pspec0.fits', rawps, pheader
fxread, datadir+'pspec1.fits', jkps, pheader
fxread, datadir+'pspec2.fits', whiteps, pheader


nx = n_elements(rawmap[*,0])
ny = n_elements(rawmap[0,*])
xcen = sxpar(mapheader, "CRPIX1")-1
ycen = sxpar(mapheader, "CRPIX2")-1

pspec_xunit = sxpar(pheader,"CUNIT1")
pspec_df = sxpar(pheader,"CDELT1")
pspec_nf = n_elements(rawps)
pspec_f = (dindgen(pspec_nf)+0.5)*pspec_df

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

pixres = sxpar( mapheader, "CDELT1" )

aspect = double(ny)/double(nx)


;window, 0, xsize=600, ysize=600*aspect

set_plot, 'ps'

!p.thick=3.
!x.thick = !p.thick
!y.thick = !p.thick

device, filename="lockman_rawmaps.eps", /encapsulated, bits_per_pixel=8, $
        xsize=20, ysize=20*aspect

loadct,0

m = 0.00003
thick =3.

cs = 1.5
ys = 0.027
xm = 0.00

pos = [0, 0.5, 0.5, 1.]
tvscale, -rawmap, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(a) raw", /normal, charsize=cs, $
        charthick=thick
contour, rawerr, xmargin=[0,0], ymargin=[0,0], $
         xstyle=5, ystyle=5, pos=pos, /noerase, $
         levels=[1.25,2.5,5.0]*min(rawerr), c_thick=thick, c_color=255


pos = [0.5, 0.5, 1., 1.]
ind = where( rawerr le 1.5*min(rawerr) )
scalej = 2

tvscale, -jkmap, minval=-m*scalej, maxval=m*scalej, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(b) jackknife", /normal, charsize=cs, $
        charthick=thick

; draw box where we measure PSPEC in normal coordinates
xbox = [-jacksize, -jacksize, +jacksize, +jacksize] + xcen
ybox = [-jacksize, +jacksize, +jacksize, -jacksize] + ycen
xbox = (xbox/double(nx))/2. + pos[0]
ybox = (ybox/double(ny))/2. + pos[1]

for i=0, 3 do begin
  i1 = i
  i2 = (i+1) mod 4

  plots, [xbox[i1],xbox[i2]], [ybox[i1],ybox[i2]], /normal, color=244, $
         thick=thick*1.5
endfor


pos = [0, 0, 0.5, 0.5]

ppos = pos + [0.09, 0.06, -0.02, -0.05]


plot, pspec_f, rawps, xtitle="!6spatial frequency ("+pspec_xunit+")", /ylog, $
      /xlog, charsize=1, /noerase, pos=ppos, charthick=!p.thick, $
      ytitle="PSD (pW!u2!n arcsec!u2!n)", xstyle=1, linestyle=2, $
      yrange=[1d-10,3d-8], ystyle=1

oplot, pspec_f, jkps/4, color=128

oplot, pspec_f, whiteps

xyouts, pos[0]+xm, pos[3]-ys, "(c) radial power spectra", $
        /normal, charsize=cs, charthick=thick

y = [1, 1, 1d-8]
y[1] = y[2]*0.7
y[0] = y[1]*0.7

oplot, [0.02, 0.04], y[0]*[1.,1.]
oplot, [0.02, 0.04], y[1]*[1.,1.], color=128
oplot, [0.02, 0.04], y[2]*[1.,1.], linestyle=2, thick=!p.thick

xyouts, 0.05, y[0], "whitened", charsize=1
xyouts, 0.05, y[1], "jackknife", charsize=1
xyouts, 0.05, y[2], "raw", charsize=1

pos = [0.5, 0, 1.0, 0.5]
tvscale, -whitemap, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(d) whitened", /normal, charsize=cs, $
        charthick=thick

device, /close


set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1


end

