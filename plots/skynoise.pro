; Plot showing bolo common-mode and WVM Tau
;
; data produced by Tim Jenness from 20110324_00018

mycolour

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

if 0 then begin
  readcol, "../data/wvm_and_com8b.asc", idx, tau, rts_num, rts_end, wvm_time, $
    com8b, format="L,D,L,D,D,D"

  readcol, "../data/s4dcom.dat", com4draw, format='D'

  i = where(finite(tau) and finite(com8b))

  idx = idx[i]
  tau = tau[i]
  rts_num = rts_num[i]
  rts_end = rts_end[i]*3600d*24d
  wvm_time = wvm_time[i]*3600d*24d
  com8b = com8b[i]

  ; different sampling for com4d, but hack assuming it happens at the same time

  j = where(finite(com4draw))
  com4draw=com4draw[j]
  t_interpol=loggen(rts_end[0],rts_end[n_elements(rts_end)-1], $
                    n_elements(com4d), /linear);
  com4d = interpol( com4draw, t_interpol, rts_end )

  u = uniq(wvm_time)
  tau = tau[u]
  wvm_time = wvm_time[u]
end

t0 = rts_end[0]


cs = 1.5

set_plot,'ps'

!p.thick=3.
!x.thick=!p.thick
!y.thick=!p.thick

xl = 0.12
xr = 0.87

yscl = 0.85
yoff = 0.14

pos = [xl,0*yscl+yoff,xr, 1.0*yscl+yoff]

device, filename='skynoise.eps', /color, bits_per_pixel=24

plot, [0], [0], xstyle=1, ystyle=5, xtitle='!6Time (s)', $
  pos=pos, charthick=!p.thick, charsize=cs, $
  xrange=[min(rts_end-t0),max(rts_end-t0)], $
  yrange=[min(com8b),max(com8b)]

oplot, rts_end-t0, com8b, color=red

axis, yaxis=0, ystyle=1, ytitle='!6COM!d850!n (pW)', $
  charthick=!p.thick, charsize=cs, color=red

axis, yaxis=1, ystyle=1, yrange=[min(com4d),max(com4d)], $
  ytitle='!6COM!d450!n (pW)', charthick=!p.thick, charsize=cs, /save, $
  color=blue

oplot, rts_end-t0, com4d, color=blue

axis, yaxis=1, ystyle=5, yrange=[min(tau),max(tau)], /save

oplot, wvm_time-t0, tau

device, /close

set_plot,'x'

!p.thick=1.
!x.thick=!p.thick
!y.thick=!p.thick

loadct,0


end
