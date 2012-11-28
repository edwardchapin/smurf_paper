readcol,"delayfiles.txt",filenames,format="A"

n = n_elements(filenames)

loadct,13

all_delay=dblarr(n)
all_max=dblarr(n)

window, 0

for i=0, n-1 do begin
  
  fxread, filenames[i], map, header
  delay = strjoin( strsplit( strmid( filenames[i], 27,6 ), "_", /extract ), $
                   "." )

  ind=where(finite(map) eq 0)
  map[ind] = 0

  m = map[500:1000,650:1150]
  tvscale, m, maxval=0.2, minval=-0.01, /keep, /noint

  all_delay[i] = delay
  all_max[i] = max(m)

  print, all_delay[i], all_max[i]

  wait,0.5
endfor

loadct,0

window, 1

fit = gaussfit( all_delay, all_max, A, nterms=4, $
                estimates=[0.22,-0.002,0.01,0] )
d_m = loggen(-0.02, 0.02, 100, /linear)
y_m = A[0]*exp( -((d_m-A[1])/A[2])^2d/2d) + A[3]

plot, all_delay, all_max, psym=1, xrange=[-0.015, 0.015], $
  xtitle='Delay (s)', ytitle='Peak', charsize=2, $
  title="Peak at delay="+string(A[1])+" s", symsize=1.5


oplot, d_m, y_m

plots, A[1]*[1,1], [0,0.25], linestyle=1


end
