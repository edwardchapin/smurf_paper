; Procedure for over-plotting a circe

pro circle, x, y, radius, aspect=aspect, oplot=oplot, rot=rot, _EXTRA=extra

  npoints = 100.0

  if keyword_set(aspect) eq 0 then aspect = 1.

  theta = findgen(npoints)*(2*3.14159/(npoints-1.))
  theta(npoints-1) = 0.0

  xoff = radius*cos(theta)*aspect
  yoff = radius*sin(theta)

  if keyword_set( rot ) then begin
    xn = xoff*cos(rot) + yoff*sin(rot)
    yn = -xoff*sin(rot) + yoff*cos(rot)

    xoff = xn
    yoff = yn
  endif

  xp = x + xoff 
  yp = y + yoff 



  ; previously used plot and oplot

  if keyword_set(oplot) then plots,xp,yp, _EXTRA=extra $
  else plots,xp,yp, _EXTRA=extra

end
