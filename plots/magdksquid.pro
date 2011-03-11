;-------------------------------------------------------------------------------
; Some basic plots for magnetic field pickup
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

; compare bolos, dksquids and pointing ----------------------------------------

cs = 1.5

loadct, 0

if 1 then begin
  subarray = 's4a'
  obs = '20100228_00016'

  fxread, datadir+subarray+obs+'_con_clean.fits', bolos, header
  fxread, datadir+subarray+obs+'_con_clean_dksquid.fits', dksquid, header2

  state = scuba2_readstate( datadir+'state_'+obs+'.tst' )
endif

t = state.rts_end
t = (t - min(t))*day2sec
nt = n_elements(t)

; 2 bolometers with different phases for the pickup, and DKS same column
b1 = bolos[5,11,*]
b2 = bolos[5,12,*]
dks = dksquid[5,*]

; convert to kDAC units
b1 = b1 / 1000d
b2 = b2 / 1000d
dks = (dks - mean(dks))/1000d

; smoothing box size
box = 200

set_plot,'ps'

!p.thick=3.
!x.thick=!p.thick
!y.thick=!p.thick

device, filename='magpickup.eps', /encapsulated, xsize=20., $
        ysize=20.

thex = 18 & they = 22

xl = 0.13
xr = 0.99

yscl = 0.91
yoff = 0.08

yt = 0.275
xt = 0.495
xrange = [min(t),max(t)]
yrange = [-9.5,9.5]

label = strarr(30)+' '

; bolos
pos = [xl,0.666*yscl+yoff,xr, 1.0*yscl+yoff]
plot, [0], [0], xstyle=5, charsize=cs, pos=pos, ytitle='Signal (kDAC)', $
      ystyle=1, charthick=!p.thick,xrange=xrange,yrange=yrange
oplot, t, smooth(b2,box,/edge_truncate)
oplot, t, smooth(b1,box,/edge_truncate), linestyle=1, thick=!p.thick/2.

axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, 'smoothed 450'+micron+' bolos', charsize=cs, $
        charthick=!p.thick, /normal

; dksquid
pos = [xl,0.333*yscl+yoff,xr, 0.666*yscl+yoff]
plot, [0], [0], xstyle=5, charsize=cs, pos=pos, ytitle='Signal (kDAC)', $
      ystyle=1, charthick=!p.thick,xrange=xrange,yrange=yrange,/noerase
oplot, t, smooth(-dks,box,/edge_truncate)

axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label
xyouts, pos[0]+xt, pos[1]+yt, 'smoothed dark squid', charsize=cs, $
        charthick=!p.thick, /normal

;pointing
pos = [xl,0.+yoff,xr, 0.333*yscl+yoff]
plot, t, state.daz/3600., xstyle=5, charsize=cs, pos=pos, $
      ytitle='Offset (deg)', /noerase, charthick=!p.thick, $
      yrange=[-0.4,0.2], ystyle=1
oplot, t, state.del/3600., linestyle=2
axis, xaxis=0, xstyle=1, xrange=xrange, xtitle='Time (s)', charthick=!p.thick, $
      charsize=cs
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label

xt = 0.37
loff = 0.25

plots, pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal
plots, loff+pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal, linestyle=2

xyouts, pos[0]+xt+0.06, pos[1]+yt, '!7'+!gr.delta+'!6Azimuth', $
        charsize=cs, charthick=!p.thick, /normal

xyouts, pos[0]+xt+loff+0.06, pos[1]+yt, '!7'+!gr.delta+'!6Elevation', $
        charsize=cs, charthick=!p.thick, /normal


device, /close



set_plot,'x'

!p.thick=1.
!x.thick=!p.thick
!y.thick=!p.thick

end
