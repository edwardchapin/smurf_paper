; Plot showing bolo common-mode and WVM Tau
;
; data produced by Tim from 20110324_00018 at 850um

if 0 then begin
  readcol, "../data/wvm_and_com8b.asc", idx, tau, rts_num, rts_end, wvm_time, $
    com, format="L,D,L,D,D,D"

  i = where(finite(tau) and finite(com))

  idx = idx[i]
  tau = tau[i]
  rts_num = rts_num[i]
  rts_end = rts_end[i]*3600d*24d
  wvm_time = wvm_time[i]*3600d*24d

  u = uniq(wvm_time)
  tau = tau[u]
  wvm_time = wvm_time[u]
end

t0 = rts_end[0]


cs = 1.5
loadct, 0

set_plot,'ps'

!p.thick=3.
!x.thick=!p.thick
!y.thick=!p.thick

xl = 0.12
xr = 0.88

yscl = 0.86
yoff = 0.13

pos = [xl,0*yscl+yoff,xr, 1.0*yscl+yoff]

device, filename='skynoise.eps', /color, bits_per_pixel=24

plot, rts_end-t0, com, xstyle=1, ystyle=5, xtitle='!6Time (s)', $
  pos=pos, charthick=!p.thick, charsize=cs

axis, yaxis=0, ystyle=1, ytitle='Power (pW)', charthick=!p.thick, charsize=cs

axis, yaxis=1, ystyle=1, yrange=[min(tau),max(tau)], ytitle='Opacity', $
  charthick=!p.thick, charsize=cs, /save

oplot, wvm_time-t0, tau, linestyle=2


device, /close

set_plot,'x'

!p.thick=1.
!x.thick=!p.thick
!y.thick=!p.thick

end
