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

; --- show jackknife and sum images for 3 different filter scales ---

scale = ['150','300','600']
nscale = n_elements(scale)

device, filename="m17_jk.eps", /encapsulated, bits_per_pixel=8, $
        xsize=20*nscale/4., ysize=20*aspect, /color

dxplot = 1./nscale
dyplot = 1./4.

for i=0, nscale-1 do begin

  ; default
  fxread, datadir+'m17_default_'+scale[i]+'_sum.fits', $
    map_de_sum, header
  bad = where( finite(map_de_sum) eq 0 )
  if bad[0] ne -1 then map_de_sum[bad] = badval

  fxread, datadir+'m17_default_'+scale[i]+'_jk.fits', $
    map_de_jk, header
  bad = where( finite(map_de_jk) eq 0 )
  if bad[0] ne -1 then map_de_jk[bad] = badval

  ; bright extended
  fxread, datadir+'m17_bright_extended_'+scale[i]+'_sum.fits', $
    map_be_sum, header
  bad = where( finite(map_be_sum) eq 0 )
  if bad[0] ne -1 then map_be_sum[bad] = badval

  fxread, datadir+'m17_bright_extended_'+scale[i]+'_jk.fits', $
    map_be_jk, header
  bad = where( finite(map_be_jk) eq 0 )
  if bad[0] ne -1 then map_be_jk[bad] = badval

  fxread, datadir+'m17_bright_extended_half1_'+scale[i]+'_qmask.fits', $
    qual_be_1, header

  fxread, datadir+'m17_bright_extended_half2_'+scale[i]+'_qmask.fits', $
    qual_be_2, header


  ; default sum
  pos = [dxplot*i,1-dyplot,dxplot*(i+1),1.]
  map = map_de_sum
  map = map < mx
  map = ((map - mn) > 0)-mn
  map = alog10(map)
  tvscale, -map, /noint, pos=pos

 ; filter scale
  len = (scale[i]/3600.)/pixres
  len_norm = (len/double(nx)) * (pos[2]-pos[0])
  xcen = (pos[2]+pos[0])/2. + 0.26*dxplot
  ycen = (pos[3]+pos[1])/2. + 0.06*dyplot
  plots, xcen+len_norm*[-0.5,0.5], ycen*[1., 1.], /normal, thick=thick*3, $
    color=255
  plots, xcen+len_norm*[-0.5,0.5], ycen*[1., 1.], /normal, thick=thick, $
    color=0
  xyouts, xcen, ycen-0.015, scale[i]+' arcsec', /normal, charsize=2.*cs/3., $
    charthick=thick*3., align=0.5, color=255
  xyouts, xcen, ycen-0.015, scale[i]+' arcsec', /normal, charsize=2.*cs/3., $
    charthick=thick, align=0.5, color=0

  ; default jk
  pos = [dxplot*i,1-2*dyplot,dxplot*(i+1),1.-dyplot]
  map = map_de_jk
  map = map < mx
  map = ((map - mn) > 0)-mn
  map = alog10(map)
  tvscale, -map, /noint, pos=pos

  ; bright extended sum
  pos = [dxplot*i,1-3*dyplot,dxplot*(i+1),1.-2*dyplot]
  map = map_be_sum
  map = map < mx
  map = ((map - mn) > 0)-mn
  map = alog10(map)
  tvscale, -map, /noint, pos=pos

  mycolour

  contour, qual_be_1, xmargin=[0,0], ymargin=[0,0], $
    xstyle=5, ystyle=5, pos=pos, /noerase, $
    levels=0.5, c_thick=1., c_color=blue

  contour, qual_be_2, xmargin=[0,0], ymargin=[0,0], $
    xstyle=5, ystyle=5, pos=pos, /noerase, $
    levels=0.5, c_thick=1., c_color=red

  loadct,0


  ; bright extended jk
  pos = [dxplot*i,1-4*dyplot,dxplot*(i+1),1.-3*dyplot]
  map = map_be_jk
  map = map < mx
  map = ((map - mn) > 0)-mn
  map = alog10(map)
  tvscale, -map, /noint, pos=pos
endfor


labels = ['Default Average','Default Jackknife','Bright Extended Average', $
          'Bright Extended Jackknife']
cs = 0.666
ys = 0.01
xm = 0.00

for i=0, 3 do begin
  xyouts, xm, 1-i*dyplot-ys, labels[i], /normal, charsize=cs, $
    charthick=thick
endfor

device, /close

; --- power spectra plots -----------------------------------------------------

col = [red,green,blue]
scale = ['150','300','600']
;scale = ['600']
type = 'pspec_m17_'+['bright_extended','default']
nscale = n_elements(scale)

fxread, datadir+'pspec_m17_fakemap.fits', fakeps, pheader
pspec_xunit = sxpar(pheader,"CUNIT1")
pspec_df = sxpar(pheader,"CDELT1")
pspec_nf = n_elements(fakeps)
pspec_f = (dindgen(pspec_nf)+0.5)*pspec_df

prange = [1d-8, max(fakeps)]

pos = [0.13,0.09,0.9,0.94]

for i=0, n_elements(type)-1 do begin

  ; --- raw power spectra ---

  device, filename=type[i]+'.eps', /encapsulated, bits_per_pixel=8, /color, $
    xsize=20, ysize=20

  plot, pspec_f, fakeps, /xlog, /ylog, xstyle=5, $
    ytitle="Raw PSD (pW!u2!n arcsec!u2!n)", ystyle=1, $
    charthick=thick, charsize=cs*2., $
    yrange=prange, thick=thick*4., pos=pos

  for j=0, n_elements(scale)-1 do begin
    mycolour
    fxread, datadir+type[i]+'_'+scale[j]+'_jk.fits', jkps, pheader
    fxread, datadir+type[i]+'_'+scale[j]+'_sum.fits', sumps, pheader

    oplot, pspec_f, jkps, color=col[j], thick=thick*2, linestyle=2
    oplot, 1./scale[j]*[1.,1.],[1d-10,1d10], linestyle=1, col=col[j]

    oplot, pspec_f, sumps, color=col[j], thick=thick*2,linestyle=0

    loadct,0
  endfor

  axis, xaxis=0, xstyle=1, xrange=[min(pspec_f),max(pspec_f)], $
    xtitle="!6spatial frequency ("+pspec_xunit+")", charsize=cs*2, $
    /xlog, charthick=thick

  axis, xaxis=1, xstyle=1, xrange=1d/[min(pspec_f),max(pspec_f)], $
    xtitle="!6angular scale (arcsec)", charsize=cs*2, $
    /xlog, charthick=thick


  xl = 0.6
  yl = 0.8
  ys = 0.04
  xt = 0.06

  plots, xl+[0.,0.05], yl*[1.,1.], thick=thick*4, /normal
  xyouts, xl+xt, yl-0.002, "input signal PSD", /normal, charsize=cs*1.5, $
    charthick=thick

  plots, xl+[0.,0.05], (yl-1*ys)*[1.,1.], thick=thick*2, /normal
  xyouts, xl+xt, yl-1*ys-0.002, "map PSDs", /normal, $
    charsize=cs*1.5, charthick=thick

  plots, xl+[0.,0.05], (yl-2*ys)*[1.,1.], thick=thick*2, linestyle=2, /normal
  xyouts, xl+xt, yl-2*ys-0.002, "jackknife PSDs", /normal, $
    charsize=cs*1.5, charthick=thick

  plots, xl+[0.,0.05], (yl-3*ys)*[1.,1.], thick=thick, linestyle=1, /normal
  xyouts, xl+xt, yl-3*ys-0.002, "filter edges", /normal, charsize=cs*1.5, $
    charthick=thick

  device,/close


  ; --- power spectra corrected for transfer function ---

  device, filename='cor_'+type[i]+'.eps', /encapsulated, bits_per_pixel=8, /color, $
    xsize=20, ysize=20


  plot, [1d-7,1d-7], fakeps, /xlog, /ylog, xstyle=5, $
    ystyle=5, charthick=thick, charsize=cs*2., $
    yrange=prange, thick=thick*2., pos=pos, $
    xrange=[min(pspec_f),max(pspec_f)]

  for j=0, n_elements(scale)-1 do begin
    mycolour
    fxread, datadir+type[i]+'_'+scale[j]+'_jk.fits', jkps, pheader
    fxread, datadir+type[i]+'_'+scale[j]+'_sum.fits', sumps, pheader

    ; work out transfer function
    xfer = jkps*0 + 1.
    ind = where(fakeps gt sumps)
    xfer[ind] = sumps[ind]/fakeps[ind]
    jkps_cor = jkps/xfer
    sumps_cor = sumps/xfer

    oplot, pspec_f, jkps_cor, color=col[j], thick=thick*2, linestyle=2
    oplot, 1./scale[j]*[1.,1.],[1d-10,1d10], linestyle=1, col=col[j]

    oplot, pspec_f, sumps_cor, thick=thick*2, color=col[j], linestyle=0

    loadct,0
  endfor

  axis, xaxis=0, xstyle=1, xrange=[min(pspec_f),max(pspec_f)], $
    xtitle="!6spatial frequency ("+pspec_xunit+")", charsize=cs*2, $
    /xlog, charthick=thick

  axis, xaxis=1, xstyle=1, xrange=1d/[min(pspec_f),max(pspec_f)], $
    xtitle="!6angular scale (arcsec)", charsize=cs*2, $
    /xlog, charthick=thick

  axis, yaxis=0, ystyle=1, yrange=prange, charsize=cs*2, $
    charthick=thick, $
    ytitle="Corrected PSD (pW!u2!n arcsec!u2!n)"

  axis, yaxis=1, ystyle=1, yrange=[0,1], charsize=cs*2, $
    charthick=thick, ytitle="Transfer function", /save, ylog=0

  for j=0, n_elements(scale)-1 do begin
    mycolour
    fxread, datadir+type[i]+'_'+scale[j]+'_sum.fits', sumps, pheader

    ; work out transfer function again...
    xfer = jkps*0 + 1.
    ind = where(fakeps gt sumps)
    xfer[ind] = sumps[ind]/fakeps[ind]

    oplot, pspec_f, xfer, color=col[j], linestyle=0, thick=0

    loadct,0
  endfor


  xl = 0.57
  yl = 0.7
  ys = 0.04
  xt = 0.06

  plots, xl+[0.,0.05], yl*[1.,1.], thick=thick*2, /normal
  xyouts, xl+xt, yl-0.002, "corrected map PSDs", /normal, charsize=cs*1.5, $
    charthick=thick

  plots, xl+[0.,0.05], (yl-ys)*[1.,1.], thick=thick*2, linestyle=2, /normal
  xyouts, xl+xt, yl-ys-0.002, "corrected jackknife PSDs", /normal, $
    charsize=cs*1.5, charthick=thick

  plots, xl+[0.,0.05], (yl-2*ys)*[1.,1.], thick=thick, linestyle=1, /normal
  xyouts, xl+xt, yl-2*ys-0.002, "filter edges", /normal, charsize=cs*1.5, $
    charthick=thick

  plots, xl+[0.,0.05], (yl-3*ys)*[1.,1.], thick=1, /normal
  xyouts, xl+xt, yl-3*ys-0.002, "transfer functions", /normal, $
    charsize=cs*1.5, charthick=thick


  device, /close

endfor


set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1

end
