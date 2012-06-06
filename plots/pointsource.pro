;------------------------------------------------------------------------------
; Use some point source data to illustrate: simple common-mode removal
; and large-scale degeneracy; effect of adding high-pass filtering;
; and adding a zero boundary condition. First you need to run the
; "pointmap" script in ../data.
;------------------------------------------------------------------------------

datadir = '../data/'

; load in the maps, ast and com models for 2 and 50 iterations without
; any additional filtering

badval=-1000

fxread, datadir+'uranus_com_2.fits', map2, mapheader
bad = where( finite(map2) eq 0 )
if bad[0] ne -1 then map2[bad] = badval

fxread, datadir+'uranus_com_100.fits', map100, header
bad = where( finite(map100) eq 0 )
if bad[0] ne -1 then map100[bad] = badval

fxread, datadir+'uranus_default.fits', map_flt, header
bad = where( finite(map_flt) eq 0 )
if bad[0] ne -1 then map_flt[bad] = badval

fxread, datadir+'uranus_brightcompact.fits', map_compact, header
fxread, datadir+'uranus_brightcompact.fits', map_mask, qualhead, extension=2
bad = where( finite(map_compact) eq 0 )
if bad[0] ne -1 then map_compact[bad] = badval

fxread, datadir+'resp_fplane_align.fits', resp, header
bad = where( finite(resp) eq 0, complement=good )
if bad[0] ne -1 then resp[bad] = 0
if good[0] ne -1 then resp[good] = 1


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


nx = n_elements(map2[*,0])
ny = n_elements(map2[0,*])

pixres = sxpar( mapheader, "CDELT1" )

aspect = double(ny)/double(nx)


;window, 0, xsize=600, ysize=600*aspect

set_plot, 'ps'

!p.thick=3.
!x.thick = !p.thick
!y.thick = !p.thick

device, filename="pointmaps.eps", /encapsulated, bits_per_pixel=8, $
        xsize=20, ysize=20*aspect

loadct,0

m = 0.002
thick =3.

cs = 1
ys = 0.027
xm = 0.00

pos = [0, 0.5, 0.5, 1.]
tvscale, -map2, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(a) COM", /normal, charsize=cs, charthick=thick
xyouts, pos[0]+xm, pos[3]-2*ys, '    n=2', /normal, charsize=cs, charthick=thick

pos = [0.5, 0.5, 1., 1.]
tvscale, -map100, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(b) COM", /normal, charsize=cs, charthick=thick
xyouts, pos[0]+xm, pos[3]-2*ys, '    n=100', /normal, charsize=cs,charthick=thick

resp = shift(resp,-12,2)
contour, resp, xmargin=[0,0], ymargin=[0,0], $
         xstyle=5, ystyle=5, pos=pos, /noerase, $
         levels=0.5, c_thick=thick*2., c_color=255

contour, resp, xmargin=[0,0], ymargin=[0,0], $
         xstyle=5, ystyle=5, pos=pos, /noerase, $
         levels=0.5, c_thick=thick, c_color=0


pos = [0, 0.0, 0.5, 0.5]
tvscale, -map_flt, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(c) COM+FLT", /normal, charsize=cs, $
        charthick=thick
xyouts, pos[0]+xm, pos[3]-2*ys, '    n=10', /normal, charsize=cs,charthick=thick
len = (200./3600.)/pixres
len_norm = (len/double(nx)) * (pos[2]-pos[0])
xcen = (pos[2]+pos[0])/2. - 0.01
plots, xcen+len_norm*[-0.5,0.5], 0.16*[1., 1.], /normal, thick=thick, $
       color=255
xyouts, xcen, 0.13, '200 arcsec', /normal, charsize=cs, charthick=thick, $
        align=0.5, color=255


pos = [0.5, 0, 1.0, 0.5]
tvscale, -map_compact, minval=-m, maxval=m, /noint, pos=pos
xyouts, pos[0]+xm, pos[3]-ys, "(d) COM+FLT                          zeromask", $
        /normal, charsize=cs, $
        charthick=thick
xyouts, pos[0]+xm, pos[3]-2*ys, '    n=6', /normal, charsize=cs,charthick=thick

contour, map_mask, xmargin=[0,0], ymargin=[0,0], $
         xstyle=5, ystyle=5, pos=pos, /noerase, $
         levels=0.5, c_linestyle=3, c_thick=thick, c_color=255


; az/el axes

xcomp = 0.02
ycomp = 0.02

arrow, xcomp, ycomp, xcomp, xcomp+0.04, $
       /normal, thick=!p.thick, hsize=!D.X_SIZE/64.
arrow, xcomp, ycomp, xcomp+0.04, ycomp, $
       /normal, thick=!p.thick, hsize=!D.X_SIZE/64.

xyouts, xcomp, ycomp+0.05, "El", charsize=1., charthick=!p.thick, /normal, $
        align=0.5

xyouts, xcomp+0.056, ycomp-0.004, "Az", charsize=1., charthick=!p.thick, $
        /normal, align=0.5

device, /close



; plot showing degeneracy between map and COM

cs = 1.5

device, filename='degeneracy.eps', /encapsulated, xsize=20., $
        ysize=15., /color, bits_per_pixel=8

if 1 then begin
  fxread, datadir+'com_2.fits', com2, header
  fxread, datadir+'com_100.fits', com100, header

  fxread, datadir+'s4a_gai_2.fits', gai2, header
  fxread, datadir+'s4a_gai_100.fits', gai100, header

  fxread, datadir+'s4a_res_2.fits', res2, header
  fxread, datadir+'s4a_res_100.fits', res100, header

  fxread, datadir+'s4a_ast_2.fits', ast2, header
  fxread, datadir+'s4a_ast_100.fits', ast100, header

  state = scuba2_readstate( datadir+"state_uranus.tst" )
endif

steptime = (max(state.rts_end) - min(state.rts_end)) / n_elements(state.rts_end)
steptime = steptime*24d*3600d

tstart = 1700;10000;5000
tlen = 10.1/steptime

t = state.rts_end[tstart:tstart+tlen-1]
t = (t - t[0])*24d*3600d

xl = 0.16
xr = 0.98

yscl = 0.87
yoff = 0.12

yt = 0.2
xt = 0.03
xrange = [min(t),max(t)]

label = strarr(30)+' '

bolx = 13;20
boly = 29;15


pos = [xl,0.5*yscl+yoff,xr, 1.0*yscl+yoff]

c2 = (com2*gai2[bolx,boly,0])[tstart:tstart+tlen-1]
c100 = (com100*gai100[bolx,boly,0])[tstart:tstart+tlen-1]

a2 = (ast2[bolx,boly,*])[tstart:tstart+tlen-1]
a100 = (ast100[bolx,boly,*])[tstart:tstart+tlen-1]

r2 = (res2[bolx,boly,*])[tstart:tstart+tlen-1]
r100 = (res100[bolx,boly,*])[tstart:tstart+tlen-1]

s = [c100-c2,a100-a2,r2,r100]

plot, [0], [0], xstyle=5, charsize=cs, pos=pos, ytitle='Signal (pW)', $
      charthick=!p.thick,xrange=xrange,yrange=[min(s),max(s)]*1.1
axis, xaxis=0, xstyle=1, xrange=xrange, xtickname=label
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label

mycolour

oplot, t, r2, thick=!p.thick*2.0
oplot, t, r100, thick=!p.thick*0.6, color=128
oplot, t, a100-a2, color=green
oplot, t, c100-c2, color=red

yt = 0.4

plots, pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal, thick=2.
plots, 0.2+pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal, thick=0.6, $
       color=128
plots, 0.4+pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal, color=green
plots, 0.6+pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal, $
       color=red

xyouts, pos[0]+xt+0.06, pos[1]+yt, '!6RES2', charsize=cs, $
        charthick=!p.thick, /normal
xyouts, pos[0]+xt+0.26, pos[1]+yt, '!6RES100', charsize=cs, $
        charthick=!p.thick, /normal, color=128
xyouts, pos[0]+xt+0.46, pos[1]+yt, '!7'+!gr.delta+'!6COM', charsize=cs, $
        charthick=!p.thick, /normal, color=green
xyouts, pos[0]+xt+0.66, pos[1]+yt, '!7'+!gr.delta+'!6AST', charsize=cs, $
        charthick=!p.thick, /normal, color=red

;plot, t, (c100-c2)[tstart:tstart+tlen-1], xstyle=1

;oplot, t, (a100-a2)[tstart:tstart+tlen-1], linestyle=2

loadct,0

pos = [xl,0.0*yscl+yoff,xr, 0.5*yscl+yoff]
d = [state.daz[tstart:tstart+tlen-1],state.del[tstart:tstart+tlen-1]]

plot, [0], [0], xstyle=5, charsize=cs, pos=pos, ytitle='Offset (arcsec)', $
      charthick=!p.thick,xrange=xrange,yrange=[min(d),max(d)], /noerase

axis, xaxis=0, xstyle=1, xrange=xrange, xtitle='Time (s)', charthick=!p.thick, $
      charsize=cs
axis, xaxis=1, xstyle=1, xrange=xrange, xtickname=label


yt = 0.4

plots, pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal
plots, 0.5+pos[0]+xt+[0,0.05], pos[1]+yt*[1.0,1.0]+0.006, /normal, linestyle=2
xyouts, pos[0]+xt+0.06, pos[1]+yt, '!7'+!gr.delta+'!6Azimuth', charsize=cs, $
        charthick=!p.thick, /normal
xyouts, pos[0]+xt+0.56, pos[1]+yt, '!7'+!gr.delta+'!6Elevation', charsize=cs, $
        charthick=!p.thick, /normal

oplot, t, state.daz[tstart:tstart+tlen-1]
oplot, t, state.del[tstart:tstart+tlen-1], linestyle=2

device, /close


set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1


end

