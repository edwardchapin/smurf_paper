;------------------------------------------------------------------------------
; Plots for lockman data reduction
;------------------------------------------------------------------------------

datadir = '../data/'

jacksize = 110 ; needs to match value in cosmology script
sncut = 3.8
subrad = 8     ; radius of submm blob circles

; load in the data

badval=-1000
badvar=1000

readcol, "radio_pos.tst", rid, rra, rdec, format="L,D,D"

fxread, datadir+'lockman_s8d_blank_field.fits', rawmap, mapheader
bad = where( finite(rawmap) eq 0 )
if bad[0] ne -1 then rawmap[bad] = badval

fxread, datadir+'lockman_s8d_blank_field.fits', rawvar, header, extension=1
bad = where( finite(rawvar) eq 0 )
if bad[0] ne -1 then rawvar[bad] = badvar
rawerr = sqrt(rawvar)

fxread, datadir+'lockman_s8d_blank_field_whitened.fits', whitemap, mapheader
bad = where( finite(whitemap) eq 0 )
if bad[0] ne -1 then whitemap[bad] = badval

fxread, datadir+'lockman_s8d_blank_field_filtered.fits', filtered, header
bad = where( finite(filtered) eq 0 )
if bad[0] ne -1 then filtered[bad] = badval

fxread, datadir+'lockman_s8d_blank_field_filtered.fits', filteredvar, header, $
        extension=1
bad = where( finite(filteredvar) eq 0 )
if bad[0] ne -1 then filteredvar[bad] = badvar
filterederr = sqrt(filteredvar)

filteredsnr = filtered/filterederr

fxread, datadir+'jackknife.fits', jkmap, header
bad = where( finite(jkmap) eq 0 )
if bad[0] ne -1 then jkmap[bad] = badval
jkmap = jkmap/2

fxread, datadir+'jackknife.fits', jkvar, header, extension=1
bad = where( finite(jkvar) eq 0 )
if bad[0] ne -1 then jkvar[bad] = badvar
jkerr = sqrt(jkvar)/2

fxread, datadir+'jackknife_whitened.fits', jkwhite, header
bad = where( finite(jkwhite) eq 0 )
if bad[0] ne -1 then jkwhite[bad] = badval
jkwhite = jkwhite/2.

fxread, datadir+'jackknife_filtered.fits', jkfilt, header
bad = where( finite(jkfilt) eq 0 )
if bad[0] ne -1 then jkfilt[bad] = badval
jkfilt = jkfilt/2.

fxread, datadir+'jackknife_filtered.fits', jkfiltvar, header, $
        extension=1
bad = where( finite(jkfiltvar) eq 0 )
if bad[0] ne -1 then jkfiltvar[bad] = badvar
jkfilterr = sqrt(jkfiltvar)
jkfilterr = jkfilterr/2

jkfiltsnr = jkfilt/jkfilterr

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

pixres = abs(sxpar( mapheader, "CDELT1" ))

aspect = double(ny)/double(nx)

;window, 0, xsize=600, ysize=600*aspect

set_plot, 'ps'

!p.thick=3.
!x.thick = !p.thick
!y.thick = !p.thick

; --- MAPS -------------------------------------------------------------------

device, filename="lockman_maps.eps", /encapsulated, bits_per_pixel=8, $
        xsize=3./2.*20, ysize=20*aspect, /color

loadct,0

m = 0.00003
thick =3.

cs = 1.5
ys = 0.027
xm = 0.00

pos = [0, 0.5, 0.333, 1.]
tvscale, -rawmap, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "!6(a) raw", /normal, charsize=cs, $
        charthick=thick

sig_levels = [1.25,2.5,5.0]*min(rawerr)
sig_thresh = sig_levels[1]

mask = rawmap*0
ind = where(rawerr le sig_thresh, complement=bad)
mask[ind] = 1
filteredsnr[bad] = -1000
jkfiltsnr[bad] = -1000

contour, rawerr, xmargin=[0,0], ymargin=[0,0], $
         xstyle=5, ystyle=5, pos=pos, /noerase, $
         levels=sig_levels, c_thick=thick, c_color=255

pos = [0.333, 0.5, 0.666, 1.0]
tvscale, -whitemap, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(b) whitened", /normal, charsize=cs, $
        charthick=thick*5.,color=255
xyouts, pos[0]+xm, pos[3]-ys, "(b) whitened", /normal, charsize=cs, $
        charthick=thick

pos = [0.666, 0.5, 1, 1.0]
tvscale, -filteredsnr, minval=-sncut, maxval=sncut, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(c) smoothed SNR", /normal, charsize=cs, $
        charthick=thick

; overplot radio catalogue: only positions that land on map and are
; within the mask

cat_pix, rra, rdec, rx, ry, mapheader
ind = where( (rx ge 0) and (rx lt nx) and (ry ge 0) and (ry lt ny) )
rid = rid[ind] & rx = rx[ind] & ry=ry[ind]
good = where(mask[rx,ry])
rid = rid[good] & rx = rx[good] & ry=ry[good]

xrpix = (rx/double(nx))/3. + pos[0] ; normal coords
yrpix = (ry/double(ny))/2. + pos[1]

mycolour
plots, xrpix, yrpix, psym=1, symsize=0.4, /normal, color=orange
loadct,0

; overplot peaks
findpeaks, filteredsnr, xmap, ymap, junk, thresh=sncut
ind = where( mask[xmap,ymap] )
xmap = xmap[ind] & ymap = ymap[ind]

xmpix = (xmap/double(nx))/3. + pos[0] ; normal coords
ympix = (ymap/double(ny))/2. + pos[1]

rad =  subrad / (nx*pixres*3600d)

mycolour
for i=0, n_elements(xmpix)-1 do begin
  circle, xmpix[i], ympix[i], rad, aspect=aspect*2./3., /oplot, $
          color=lightblue, thick=thick, /normal
endfor
loadct,0






pos = [0, 0, 0.333, 0.5]
ind = where( rawerr le 1.5*min(rawerr) )

tvscale, -jkmap, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(d) JK", /normal, charsize=cs, $
        charthick=thick

; draw box where we measure PSPEC in normal coordinates
xbox = [-jacksize, -jacksize, +jacksize, +jacksize] + xcen
ybox = [-jacksize, +jacksize, +jacksize, -jacksize] + ycen
xbox = (xbox/double(nx))/3. + pos[0]
ybox = (ybox/double(ny))/2. + pos[1]

for i=0, 3 do begin
  i1 = i
  i2 = (i+1) mod 4

  plots, [xbox[i1],xbox[i2]], [ybox[i1],ybox[i2]], /normal, color=244, $
         thick=thick*1.5
endfor

pos = [0.333, 0, 0.666, 0.5]
tvscale, -jkwhite, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(e) JK whitened", /normal, charsize=cs, $
        charthick=thick*5., color=255
xyouts, pos[0]+xm, pos[3]-ys, "(e) JK whitened", /normal, charsize=cs, $
        charthick=thick

pos = [0.666, 0.0, 1, 0.5]
tvscale, -jkfiltsnr, minval=-sncut, maxval=sncut, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(f) JK smoothed SNR", /normal, charsize=cs, $
        charthick=thick

; overplot radio catalogue: only positions that land on map and are
; within the mask

xrpix = (rx/double(nx))/3. + pos[0] ; normal coords
yrpix = (ry/double(ny))/2. + pos[1]

mycolour
plots, xrpix, yrpix, psym=1, symsize=0.4, /normal, color=orange
loadct,0

; overplot peaks
findpeaks, jkfiltsnr, xjk, yjk, junk, thresh=sncut
ind = where( mask[xjk,yjk] )
xjk = xjk[ind] & yjk = yjk[ind]

xjpix = (xjk/double(nx))/3. + pos[0] ; normal coords
yjpix = (yjk/double(ny))/2. + pos[1]

rad =  subrad / (nx*pixres*3600d)

mycolour
for i=0, n_elements(xjpix)-1 do begin
  circle, xjpix[i], yjpix[i], rad, aspect=aspect*2./3., /oplot, $
          color=lightblue, thick=thick, /normal
endfor
loadct,0

device, /close

; --- FILTERED PDF -------------------------------------------------------------

fxread, datadir+"psf.fits", rawpsf,psfhead
fxread, datadir+"lockman_psf_mapfilt.fits", mapfiltpsf,psfhead
fxread, datadir+"lockman_psf_mapfilt_whitened.fits", whitenedpsf,psfhead

ind = where(rawpsf eq max(rawpsf))
xc = (ind mod nx)[0]
yc = (ind / nx)[0]

x = (dindgen(nx)-xc)*pixres*3600d

pos = [0.1,0.14,0.98,0.99]

device, filename="lockman_psf.eps", /encapsulated, xsize=20, ysize=15

plot, x, rawpsf[*,yc], xtitle="!6R.A. offset (arcsec)", $
      ytitle="Relative response", xstyle=1, $
      yrange=[-0.1,1], ystyle=1, charsize=cs, charthick=thick, $
      xrange=[-60,60], pos=pos

oplot, x, mapfiltpsf[*,yc], linestyle=1
oplot, x, whitenedpsf[*,yc], linestyle=2

plots, [0.12, 0.17], 0.9*[1.,1.], /normal, thick=thick
xyouts, 0.18, 0.89, 'Gaussian FWHM=14.5"', charsize=cs, charthick=thick, /normal

plots, [0.12, 0.17], 0.8*[1.,1.], /normal, thick=thick, linestyle=1
xyouts, 0.18, 0.79, 'map filtered', charsize=cs, charthick=thick, /normal

plots, [0.12, 0.17], 0.7*[1.,1.], /normal, thick=thick, linestyle=2
xyouts, 0.18, 0.69, 'whitened', charsize=cs, charthick=thick, /normal

device, /close

; --- POWER SPECTRA ------------------------------------------------------------

device, filename="lockman_pspec.eps", /encapsulated, xsize=20, ysize=20, $
        bits_per_pixel=8, /color

pos = [0.13,0.09,0.99,0.99]

plot, pspec_f, whiteps, xtitle="!6spatial frequency ("+pspec_xunit+")", /ylog, $
      /xlog, /noerase, $
      ytitle="PSD (pW!u2!n arcsec!u2!n)", xstyle=1, $
      yrange=[1d-10,3d-8], ystyle=1, charsize=cs, charthick=thick, pos=pos

mycolour
oplot, pspec_f, jkps/4, color=orange
loadct,0

oplot, pspec_f, rawps, linestyle=2

y = [1, 1, 1d-8]
y[1] = y[2]*0.7
y[0] = y[1]*0.7

mycolour
oplot, [0.02, 0.037], y[1]*[1.,1.], color=orange
loadct,0
xyouts, 0.045, y[1]*0.97, "jackknife", charsize=cs, charthick=thick


oplot, [0.02, 0.037], y[2]*[1.,1.], linestyle=2, thick=!p.thick
xyouts, 0.045, y[2]*0.97, "raw", charsize=cs, charthick=thick

oplot, [0.02, 0.037], y[0]*[1.,1.]
xyouts, 0.045, y[0]*0.97, "whitened", charsize=cs, charthick=thick


; --- FILTERED PDF -------------------------------------------------------------

ind = where(mask)

n = n_elements(ind)

jksnr = (jkmap/jkerr)[ind]
rawsnr = (rawmap/rawerr)[ind]

minsnr=-5
maxsnr=8
nsnr = (maxsnr-minsnr)*10 + 1


rawhist = histogram( rawsnr, min=minsnr, max=maxsnr, nbins=nsnr, $
                     locations=bins )
jkhist = histogram( jksnr, min=minsnr, max=maxsnr, nbins=nsnr, locations=bins )

filteredhist = histogram( filteredsnr[ind], min=minsnr, max=maxsnr, nbins=nsnr,$
                          locations=bins )
jkfilthist = histogram( jkfiltsnr[ind], min=minsnr, max=maxsnr, nbins=nsnr, $
                        locations=bins )

print, "Raw: mean=", mean(rawsnr), " sig=", stdev(rawsnr)
print, " JK: mean=", mean(jksnr), " sig=", stdev(jksnr)
print, "-------------------------------------------------"
print, "Flt: mean=", mean(filteredsnr[ind]), " sig=", stdev(filteredsnr[ind])
print, "JFl: mean=", mean(jkfiltsnr[ind]), " sig=", stdev(jkfiltsnr[ind])

dsnr = bins[1]-bins[0]
snr = bins+dsnr/2.

g = exp(-snr^2d/2d)
g = n * g / total(g)

device, filename="lockman_hist.eps", /encapsulated, xsize=20, ysize=20


pos = [0.1,0.54,0.99,0.99]

plot, snr, g, xtitle="S/N", ytitle="N pixels", $
      charsize=1.5, charthick=thick, yrange=[1,max(jkhist)], $
      ytickformat="exponent", linestyle=2, /ylog, xstyle=5, pos=pos
xr = [min(snr),max(snr)]

xyouts, 0.12, pos[3]-0.04, "whiten", charsize=cs, charthick=thick, /normal

mycolour
oplot, snr, jkhist, psym=10, color=orange
loadct,0

mycolour
oplot, snr, rawhist, psym=10, color=blue
loadct,0

axis, xaxis=1, xrange=xr, xstyle=1, xtickformat="(A1)"
axis, xaxis=0, xrange=xr, xstyle=1, xtickformat="(A1)"


pos = [0.1,0.09,0.99,0.54]

plot, snr, g, xtitle="S/N", ytitle="N pixels", $
      charsize=1.5, charthick=thick, yrange=[1,max(jkhist)], $
      ytickformat="exponent", linestyle=2, xstyle=5, pos=pos, /noerase, /ylog

xyouts, 0.12, pos[3]-0.04, "whiten + matched-filter", charsize=cs, $
        charthick=thick, /normal

mycolour
oplot, snr, jkfilthist, psym=10, color=orange
loadct,0

mycolour
oplot, snr, filteredhist, psym=10, color=blue
loadct,0

axis, xaxis=1, xrange=xr, xstyle=1, xtickformat="(A1)"
axis, xaxis=0, xrange=xr, xstyle=1, xtitle="S/N", charsize=cs, charthick=thick

mycolour
plots, [0.65,0.7], (pos[3]-0.05)*[1.,1.], color=blue, /normal
plots, [0.65,0.7], (pos[3]-0.1)*[1.,1.], color=orange, /normal
plots, [0.65,0.7], (pos[3]-0.15)*[1.,1.], linestyle=2, /normal
loadct,0

xyouts, 0.72, pos[3]-0.055, "signal map", /normal, charthick=thick, charsize=cs
xyouts, 0.72, pos[3]-0.105, "JK map", /normal, charthick=thick, charsize=cs
xyouts, 0.72, pos[3]-0.155, "ideal Gaussian", /normal, charthick=thick, $
        charsize=cs

device, /close

set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1


end

