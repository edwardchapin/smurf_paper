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


; first, compare s4a, s8d, and mixtemp, and scan position ----------------------

cs = 1.5

if 0 then begin
  ; cookbook uranus
  ;obs = '20091214_00015'

  ; cookbook debris disk
  obs = '20100313_00029'

  fxread, datadir+'s4a'+obs+'_con_clean.fits', data450, header
  fxread, datadir+'s8d'+obs+'_con_clean.fits', data850, header

  fxread, datadir+'s4a'+obs+'_con_nocom.fits', nocom450, header
  fxread, datadir+'s8d'+obs+'_con_nocom.fits', nocom850, header

  state = scuba2_readstate( datadir+'state_'+obs+'.tst' )
endif

t = state.rts_end
t = (t - min(t))*day2sec

set_plot,'ps'

!p.thick=3.
!x.thick=!p.thick
!y.thick=!p.thick

device, filename='bolos_point_mix.eps', /encapsulated, xsize=20., $
        ysize=30.

thex = 16 & they = 15

b450 = data450[thex,they,*]
b850 = data850[thex,they,*]

nc450 = nocom450[thex,they,*]
nc850 = nocom850[thex,they,*]

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
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, '450 Bolo', charsize=cs, charthick=!p.thick,$
        /normal

pos = [xl,0.5*yscl+yoff,xr, 0.75*yscl+yoff]
plot, [0], [0], xstyle=5, charsize=cs, pos=pos, ytitle='Power (pW)', /noerase, $
      ystyle=1, charthick=!p.thick,xrange=xrange,yrange=[min(b850),max(b850)]
oplot, t, nc850, color=128
oplot, t, b850
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, '850 Bolo', charsize=cs, charthick=!p.thick,$
        /normal

pos = [xl,0.25*yscl+yoff,xr, 0.5*yscl+yoff]
m = state.sc2_mixtemp
m = m-mean(m)
plot, t, m*1d3, xstyle=5, charsize=cs, pos=pos, ystyle=1,$
      ytitle='Temperature (DAC)', /noerase, charthick=!p.thick, $
      yrange=[-0.09,0.13]
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

cs = 1.5

nt = n_elements(t)
df = 1d / (max(t)/double(nt))
nf = nt / 2
freq = (df/2d)*dindgen(nf)/double(nf)

x450 = [ 3,10,16,22];,28]
y450 = [30,13,16,31];,10]

x850 = [ 5,13,16,25]
y850 = [10,10,16,20]

box = 10

device, filename='pspec.eps', xsize=20, ysize=20, /color, $
        bits_per_pixel=24

xl = 0.14
xr = 0.99

yscl = 0.89
yoff = 0.09

xrange = [0.5,max(freq)]
yrange = [1d-7,1d-3]

pos = [xl,0.5*yscl+yoff,xr, 1.0*yscl+yoff]

loadct,0

plot, [0], [0], /xlog, /ylog, xrange=xrange, $
      yrange=yrange, xstyle=1, ystyle=1, $
      ytitle='PSD (pW!u2!n Hz!u-1!n)', $
      charsize=cs, charthick=!p.thick, pos=pos, xtickname=label

mycolour

for i=0, n_elements(x450)-1 do begin
  f = fft(data450[x450[i],y450[i],*])
  p = nt*abs(f)^2d
  oplot, freq, smooth(p,box), color=3-i, thick=1., linestyle=1

  f = fft(nocom450[x450[i],y450[i],*])
  p = nt*abs(f)^2d
  oplot, freq, smooth(p,box), color=3-i
endfor

pos = [xl,0.0*yscl+yoff,xr, 0.5*yscl+yoff]

loadct,0
plot, [0], [0], /xlog, /ylog, xrange=xrange, $
      yrange=yrange, xstyle=1, ystyle=1, $
      ytitle='PSD (pW!u2!n Hz!u-1!n)', $
      charsize=cs, charthick=!p.thick, pos=pos, xtickname=label, /noerase

axis, xaxis=0, xstyle=1, xrange=xrange, xtitle='Frequency (Hz)', $
      charsize=cs, charthick=!p.thick

mycolour
for i=0, n_elements(x850)-1 do begin
  f = fft(data850[x850[i],y850[i],*])
  p = nt*abs(f)^2d
  oplot, freq, smooth(p,box), color=3-i, thick=1, linestyle=1

  f = fft(nocom850[x450[i],y450[i],*])
  p = nt*abs(f)^2d
  oplot, freq, smooth(p,box), color=3-i
endfor

device, /close

set_plot,'x'

!p.thick=1.
!x.thick=!p.thick
!y.thick=!p.thick

end
