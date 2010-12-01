;-------------------------------------------------------------------------------
; Some basic plots showing bolometer time series and correlation with
; the fridge oscillations. Also show some basic power spectra. To
; generate the input files, run prepdata.sh in ../data
;-------------------------------------------------------------------------------

datadir = '../data/'

day2sec = 24d*3600d
r2a = (180d/!DPI)*3600d

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

micron = '!7'+!gr.mu+'!6m'

; first, compare s4a, s8d, and mixtemp, and scan position ----------------------

cs = 1.5

loadct, 0

if 0 then begin
  ; cookbook uranus
  ;obs = '20091214_00015'

  ; cookbook debris disk
  obs = '20100313_00029'

  ; lockman
  ;obs = '20100311_00065'

  fxread, datadir+'s4a'+obs+'_con_clean.fits', data450, header
  fxread, datadir+'s8d'+obs+'_con_clean.fits', data850, header

  fxread, datadir+'s4a'+obs+'_con_nocom.fits', nocom450, header
  fxread, datadir+'s8d'+obs+'_con_nocom.fits', nocom850, header

  state = scuba2_readstate( datadir+'state_'+obs+'.tst' )
endif

t = state.rts_end
t = (t - min(t))*day2sec
nt = n_elements(t)

set_plot,'ps'

!p.thick=3.
!x.thick=!p.thick
!y.thick=!p.thick

device, filename='bolos_point_mix.eps', /encapsulated, xsize=20., $
        ysize=30.

;thex = 16 & they = 15
thex = 18 & they = 22

b450 = data450[thex,they,*]
b850 = data850[thex,they,*]

nc450 = nocom450[thex,they,*]
nc850 = nocom850[thex,they,*]

; common mode

com450 = b450 - nc450
com850 = b850 - nc850

com450[nt-1] = com450[nt-2]
com850[nt-1] = com850[nt-2]


xl = 0.15
xr = 0.99

yscl = 0.93
yoff = 0.06

yt = 0.2
xt = 0.03
xrange = [min(t),max(t)]

label = strarr(30)+' '

pos = [xl,0.75*yscl+yoff,xr, 1.0*yscl+yoff]
plot, [0], [0], xstyle=5, charsize=cs, pos=pos, ytitle='Power (pW)', ystyle=1, $
      charthick=!p.thick,xrange=xrange,yrange=[min(b450),max(b450)]
oplot, t, nc450, color=128
oplot, t, b450
;oplot, t, com450, color=128, linestyle=2
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, '450'+micron, charsize=cs, $
        charthick=!p.thick, /normal

pos = [xl,0.5*yscl+yoff,xr, 0.75*yscl+yoff]
plot, [0], [0], xstyle=5, charsize=cs, pos=pos, ytitle='Power (pW)', /noerase, $
      ystyle=1, charthick=!p.thick,xrange=xrange,yrange=[min(b850),max(b850)]
oplot, t, nc850, color=128
oplot, t, b850
;oplot, t, com850, color=128, linestyle=2
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, '850'+micron, charsize=cs, $
        charthick=!p.thick, /normal

pos = [xl,0.25*yscl+yoff,xr, 0.5*yscl+yoff]
m = state.sc2_mixtemp
m = m-mean(m)
plot, t, m*1d3, xstyle=5, charsize=cs, pos=pos, ystyle=1,$
      ytitle='Temperature (mK)', /noerase, charthick=!p.thick, $
      yrange=[-0.07 ,0.07 ]
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, 'Mixing Chamber', charsize=cs,$
        charthick=!p.thick, /normal

pos = [xl,0.+yoff,xr, 0.25*yscl+yoff]
plot, t, state.daz, xstyle=5, charsize=cs, pos=pos, $
      ytitle='Offset (arcsec)', /noerase, charthick=!p.thick, $
      yrange=[-300,300]
oplot, t, state.del, linestyle=2
axis, xaxis=0, xstyle=1, xrange=xrange, xtitle='Time (s)', charthick=!p.thick, $
      charsize=cs
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label

plots, pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal
plots, 0.5+pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal, linestyle=2

xyouts, pos[0]+xt+0.06, pos[1]+yt, '!7'+!gr.delta+'!6Azimuth', charsize=cs, $
        charthick=!p.thick, /normal

xyouts, pos[0]+xt+0.56, pos[1]+yt, '!7'+!gr.delta+'!6Elevation', charsize=cs, $
        charthick=!p.thick, /normal


device, /close

; power spectra ----------------------------------------------------------------

col = [grey,red,green,blue]

cs = 1.5

srate = 1d / (max(t)/double(nt))
nf = nt / 2
freq = (srate/2d)*dindgen(nf)/double(nf)
df = srate / nt


x450 = [ 3,10,16,thex];,28]
y450 = [30,13,16,they];,10]

;x850 = [ 5,13,16,thex]
;y850 = [10,10,16,they]

x850 = [ 4,31,25,thex]
y850 = [16,20,32,they]

box = round(0.1/df) ; width of boxcar in Hz / freq. step size

device, filename='pspec.eps', xsize=20, ysize=20, /color, $
        bits_per_pixel=24

xl = 0.14
xr = 0.99

yscl = 0.845
yoff = 0.09

vel = 200.  ; arcsec/sec

xrange = [0.1,max(freq)]
yrange = [2d-10,1.5d-4]

pos = [xl,0.5*yscl+yoff,xr, 1.0*yscl+yoff]

loadct,0

plot, [0], [0], /xlog, /ylog, xrange=xrange, $
      yrange=yrange, xstyle=5, ystyle=1, $
      ytitle='PSD (pW!u2!n Hz!u-1!n)', $
      charsize=cs, charthick=!p.thick, pos=pos, xtickname=label

thetarange = (1/xrange)*vel/60.

axis, xaxis=1, xrange=thetarange, xstyle=1, xtitle='Angular scale (arcmin)', $
      charsize=cs, charthick=!p.thick

axis, xaxis=0, xrange=xrange, xstyle=1, xtickname=label


mycolour

; legend
xyouts, 0.85, pos[3]-0.04, '450'+micron, charsize=cs, charthick=!p.thick, $
        /normal

xt = 0.6
yt = 0.43 ;0.90

dx = 0.025
dy = 0.03

for i=0, 3 do begin
  plots, xt+[i*dx,(i+1)*dx],      yt*[1.,1.], linestyle=1.,color=col[i], /normal
  plots, xt+[i*dx,(i+1)*dx], (yt-dy)*[1.,1.], linestyle=0.,color=col[i], /normal
endfor
plots, xt+[0,4*dx], (yt-2*dy)*[1.,1.], color=black, /normal
plots, xt+[0,4*dx], (yt-3*dy)*[1.,1.], color=black, linestyle=2, /normal

xyouts, xt+5*dx, yt-0.005, 'raw bolos', charsize=cs, charthick=!p.thick, $
        /normal
xyouts, xt+5*dx, (yt-dy)-0.005, 'com cleaned', charsize=cs, $
        charthick=!p.thick, /normal
xyouts, xt+5*dx, (yt-2*dy)-0.005, 'common-mode', charsize=cs, $
        charthick=!p.thick, /normal
xyouts, xt+5*dx, (yt-3*dy)-0.005, 'PSF', charsize=cs, $
        charthick=!p.thick, /normal


; 450 power spectra
for i=0, n_elements(x450)-1 do begin
  f = fft(data450[x450[i],y450[i],*])
  p = (abs(f)^2d)/df
  oplot, freq, smooth(p,box), color=col[i], linestyle=1

  f = fft(nocom450[x450[i],y450[i],*])
  p = (abs(f)^2d)/df
  oplot, freq, smooth(p,box), color=col[i]
endfor

; 450 PSF
fwhm_f = vel / 7.5
psf = 1d-6 * exp( -freq^2d/(2d*(2.35*fwhm_f)^2d) )
oplot, freq, psf, linestyle=2, color=black

; 450 common-mode power spectrum
f = fft(com450)
p = (abs(f)^2d)/df
oplot, freq, smooth(p,box), color=white, thick=!p.thick*3.
oplot, freq, smooth(p,box), color=black, thick=!p.thick*1.

; 450 reference noise value
oplot, [1d-10,1d10], 5d-8*[1.,1.], linestyle=1, thick=!p.thick*3.,color=white
oplot, [1d-10,1d10], 5d-8*[1.,1.], linestyle=1, color=black

axis, yaxis=1, yrange=yrange, ytickname=label, ystyle=1

pos = [xl,0.0*yscl+yoff,xr, 0.5*yscl+yoff]

loadct,0
plot, [0], [0], /xlog, /ylog, xrange=xrange, $
      yrange=yrange, xstyle=1, ystyle=1, $
      ytitle='PSD (pW!u2!n Hz!u-1!n)', $
      charsize=cs, charthick=!p.thick, pos=pos, xtickname=label, /noerase

axis, xaxis=0, xstyle=1, xrange=xrange, xtitle='Frequency (Hz)', $
      charsize=cs, charthick=!p.thick

mycolour

; legend
xyouts, 0.85, pos[3]-0.04, '850'+micron, charsize=cs, charthick=!p.thick, $
        /normal

; 850 power spectra
for i=0, n_elements(x850)-1 do begin
  f = fft(data850[x850[i],y850[i],*])
  p = (abs(f)^2d)/df
  oplot, freq, smooth(p,box), color=col[i], linestyle=1

  f = fft(nocom850[x850[i],y850[i],*])
  p = (abs(f)^2d)/df
  oplot, freq, smooth(p,box), color=col[i]
endfor

; 850 common-mode power spectrum
f = fft(com850)
p = (abs(f)^2d)/df
oplot, freq, smooth(p,box), color=white, thick=!p.thick*3.
oplot, freq, smooth(p,box), color=black, thick=!p.thick*1.

; 850 PSF
fwhm_f = vel / 14.5
psf = 1d-7 * exp( -freq^2d/(2d*(2.35*fwhm_f)^2d) )
oplot, freq, psf, linestyle=2, color=black

; 850 reference noise value
oplot, [1d-10,1d10], 6d-9*[1.,1.], linestyle=1, thick=!p.thick*3.,color=white
oplot, [1d-10,1d10], 6d-9*[1.,1.], linestyle=1, color=black

axis, yaxis=1, yrange=yrange, ytickname=label, ystyle=1

device, /close

set_plot,'x'

!p.thick=1.
!x.thick=!p.thick
!y.thick=!p.thick

end
