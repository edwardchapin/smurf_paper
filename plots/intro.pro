;-------------------------------------------------------------------------------
; some basic plots showing bolometer time series and correlation with
; the fridge oscialltions
;-------------------------------------------------------------------------------

datadir = '../data/'

day2sec = 24d*3600d
r2a = (180d/!DPI)*3600d

; first, compare s4a, s8d, and mixtemp, and scan position

if 0 then begin
  ; these are concatenations of files 2 & 3 and then run through
  ; sc2clean with dimmconfig.lis
  fxread, datadir+'s4a20091214_00015_con_clean.fits', data450, header
  fxread, datadir+'s8d20091214_00015_con_clean.fits', data850, header
  state = scuba2_readstate( datadir+'state_20091214_00015.tst' )
endif

t = state.rts_end
t = (t - min(t))*day2sec

set_plot,'ps'

!p.thick=3.
!x.thick=!p.thick
!y.thick=!p.thick

device, filename='bolos_point_mix.eps', /encapsulated, xsize=20., $
        ysize=30.

b450 = data450[16,16,*]
b850 = data850[16,16,*]

xl = 0.15
xr = 0.99

yscl = 0.93
yoff = 0.06

yt = 0.2
xt = 0.03
xrange = [min(t),max(t)]

label = strarr(30)+' '

pos = [xl,0.75*yscl+yoff,xr, 1.0*yscl+yoff]
plot, t, b450, xstyle=5, charsize=1.5, pos=pos, ytitle='Power (pW)', ystyle=1, $
      charthick=!p.thick
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, '450 Bolo', charsize=1.5, charthick=!p.thick,$
        /normal

pos = [xl,0.5*yscl+yoff,xr, 0.75*yscl+yoff]
plot, t, b850, xstyle=5, charsize=1.5, pos=pos, ytitle='Power (pW)', /noerase, $
      ystyle=1, charthick=!p.thick
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, '850 Bolo', charsize=1.5, charthick=!p.thick,$
        /normal

pos = [xl,0.25*yscl+yoff,xr, 0.5*yscl+yoff]
m = state.sc2_mixtemp
m = m-mean(m)
plot, t, m*1d3, xstyle=5, charsize=1.5, pos=pos, ystyle=1,$
      ytitle='Temperature (DAC)', /noerase, charthick=!p.thick, $
      yrange=[-0.09,0.13]
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, 'Mixing Chamber', charsize=1.5,$
        charthick=!p.thick, /normal

pos = [xl,0.+yoff,xr, 0.25*yscl+yoff]
plot, t, state.daz, xstyle=5, charsize=1.5, pos=pos, $
      ytitle='Offset (arcsec)', /noerase, charthick=!p.thick, $
      yrange=[-200,200]
oplot, t, state.del, linestyle=2
axis, xaxis=0, xstyle=1, xrange=xrange, xtitle='Time (s)', charthick=!p.thick, $
      charsize=1.5
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label

plots, pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal
plots, 0.5+pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal, linestyle=2

xyouts, pos[0]+xt+0.06, pos[1]+yt, '!7'+!gr.delta+'!6Azimuth', charsize=1.5, $
        charthick=!p.thick, /normal

xyouts, pos[0]+xt+0.56, pos[1]+yt, '!7'+!gr.delta+'!6Elevation', charsize=1.5, $
        charthick=!p.thick, /normal


device, /close

set_plot,'x'

!p.thick=1.
!x.thick=!p.thick
!y.thick=!p.thick

end
