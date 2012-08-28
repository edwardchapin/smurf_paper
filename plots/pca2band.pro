;------------------------------------------------------------------------------
; PCA some bolo data at 450 and 850 simultaneously
;------------------------------------------------------------------------------

pre='./2band_'
datadir='../data/'

; debris disk (paper!)
;obs = '20100313_00029'

; uranus
;obs = '20091214_00015'

; lockman
;obs = '20100311_00065'

; 2011 data of NGC207IR post-upgrade (s4a,s4c,s4d,s8a,s8b,s8d)
;obs = '20110203_00016'

; 2011 obs recommended by Harriet. 900" pong.
obs = '20111112_00038'

subarr850 = 's8b'
subarr450 = 's4a'

ncomp = 32    ; single component plots (up to this value)
nrem = 0      ; number of components to remove
dobigplot = 0 ; do the big plots?

micron = '!7'+!gr.mu+'!6m'
day2sec = 24d*3600d

if 0 then begin

  ; load bolo data and mask
  fxread, datadir+subarr850+obs+'_con_clean.fits', data8, head8
  fxread, datadir+subarr450+obs+'_con_clean.fits', data4, head4

  fxread, datadir+'mask_'+subarr850+'.fits', mask8, mhead4
  fxread, datadir+'mask_'+subarr450+'.fits', mask4, mhead4

  state = scuba2_readstate( datadir+'state_'+obs+'.tst' )
  t = state.rts_end
  t = (t - min(t))*day2sec

  ; check that good masked data don't have NaNs
  nx = n_elements( data4[*,0,0] )
  ny = n_elements( data4[0,*,0] )
  nbolo = nx*ny

  ; figure out masks
  bad = where(mask8 ne 0, complement=good)
  mask8[good] = 1
  mask8[bad] = 0

  bad = where(mask4 ne 0, complement=good)
  mask4[good] = 1
  mask4[bad] = 0

  good4 = where(mask4)
  good8 = where(mask8)

  valbad = data4[0]

  nt = n_elements(data4)/nbolo

  df = 1d/(nt*0.005d)
  freq = (dindgen(nt))*df

  ; use only a subset of the detectors (1/2 say)
  fact = 1 ; 8

  ind4 = where(mask4)
  ind8 = where(mask8)

  if fact ne 1 then begin
    themin = min([n_elements(ind4),n_elements(ind8)])

    m = mask4*0
    for i=0, themin/fact-1 do m[ind4[i*fact]] = 1
    mask4 = m

    m = mask8*0
    for i=0, themin/fact-1 do m[ind8[i*fact]] = 1
    mask8 = m
  endif

  ; combine data from each band into one set and merge masks

  newdata = dblarr(nx*2,ny,nt)
  newdata[0:nx-1,*,*] = data4
  newdata[nx:2*nx-1,*,*] = data8

  ind4 = where(mask4)
  xg4 = ind4 mod nx
  yg4 = ind4 / nx

  ind8 = where(mask8)
  xg8 = (ind8 mod nx) + nx
  yg8 = ind8 / nx

  xg = [xg4,xg8]
  yg = [yg4,yg8]

  n = n_elements(xg)

  ; subtract off the mean signal of each detector
  for i=0, n-1 do begin
      newdata[xg[i],yg[i],*] = newdata[xg[i],yg[i],*] - $
        mean(newdata[xg[i],yg[i],*])
  endfor

  ; measure the covariance matrix and do PCA

  print, "Measuring covariance matrix"
  c = dblarr(n,n)

  for i=0, n-1 do begin
    print, i, ' / ', n

    for j=i, n-1 do begin
      c[i,j] = correlate(newdata[xg[i],yg[i],*], newdata[xg[j],yg[j],*], $
                         /covariance)
    endfor

    for j=0, i-1 do begin
      c[i,j] = c[j,i]
    endfor
  endfor

  la_svd, c, lambda, u, v ; eigenval, eigenvect (cols), right-hand matrix

  ; calculate eigenvectors

  print, "Calculating eigenvectors"
  eof = dblarr( n, nt )         ; component, time
  for i=0, n-1 do begin
      print, i, ' / ', n
      for j=0, n-1 do begin
          eof[i,*] = eof[i,*] + newdata[xg[j],yg[j],*] * u[i,j]
      endfor
      eof[i,*] = eof[i,*] / sqrt(total(eof[i,*]^2d))
  endfor

  ; project data along eigenvectors to get principal components

  print, "Projecting data along eigenvectors"
  pc = dblarr( n, n )           ; bolometer, component

  for i=0, n-1 do begin
      print, i, ' / ', n
      for j=0, n-1 do begin
          pc[j,i] = total(newdata[xg[j],yg[j],*]*eof[i,*])
      endfor
  endfor

  ; save for future use
  print, 'Saving PCA data...'
  save, data4, data8, mask4, mask8, head4, head8, $
        valbad, nx, ny, nbolo, nt, df, freq, newdata, nt, xg, yg, n, c, eof, $
        pc, t, xg4, yg4, $
        filename=pre+'data.sav'
endif else begin
  print, 'Restoring PCA data...'
  restore,pre+'data.sav'
endelse

; calculate correct frequency
srate = 1d / (max(t)/double(nt))
nf = nt / 2
freq = (srate/2d)*dindgen(nf)/double(nf)
df = srate / nt

; flip sign of eigenvectors so that eigenvalues are positive
for i=0, n-1 do begin
  lambda = mean( pc[*,i] )

  if( lambda lt 0 ) then begin
    pc[*,i] = -pc[*,i]
    eof[i,*] = -eof[i,*]
  endif

endfor

; dimensions for the big combined array
nxb = nx*2
nyb = ny

; remove first n-components from the data

if nrem gt 0 then begin
  print, "Removing first ",nrem," components"
  newdata_pca = newdata*0
  for i=0, n-1 do begin
    newdata_pca[xg[i],yg[i],*] = newdata[xg[i],yg[i],*]
    for j=0, nrem-1 do begin
      newdata_pca[xg[i],yg[i],*] = newdata_pca[xg[i],yg[i],*] - $
                                   eof[j,*]*pc[i,j]
    endfor
  endfor
endif

; plots of the bolo data, principal components, and residual after
; removing first mode
set_plot,'ps'

if dobigplot then begin
  print, "Doing large plots:"

  loadct, 0

  !p.thick=1.
  !x.thick=1.
  !y.thick=1.

  ;sx = 3
  sx = 4
  sy = 4

  print, 'Bolos...'

  !p.multi=[0,sx,sy]
  device, filename=pre+'bolos.ps', xsize=27.94, ysize=21.59

  for i=0, n-1 do begin
    if xg[i] ge nx then array=subarr850 $
    else array=subarr450

    plot, newdata[xg[i],yg[i],*], xtitle="Sample #", ytitle="Raw Data", $
          title=array+strcompress(xg[i]+1)+strcompress(yg[i]+1), $
          psym=5, ystyle=1, symsize=0.01, yrange=[-0.2,0.2], xstyle=1
  endfor

  device,/close

  print, 'Modes...'

  !p.multi=[0,1,1]
  !p.multi=[0,sx,sy]
  device, filename=pre+'modes.ps'

  for i=0, n-1 do $
    plot, eof[i,*], xtitle="Sample #", ytitle="Normalized Signal", $
          title='Component '+strcompress(i+1), $
          psym=5, ystyle=1, symsize=0.01, yrange=[-0.05,0.05], xstyle=1

  device,/close

  if nrem gt 0 then begin
    print, 'Bolos_clean...'

    !p.multi=[0,1,1]
    !p.multi=[0,sx,sy]
    device, filename=pre+'bolos_clean.ps'

    for i=0, n-1 do begin
      if xg[i] ge nx then array=subarr850 $
      else array=subarr450

      plot, newdata_pca[xg[i],yg[i],*], xtitle="Sample #", $
            ytitle="Raw Cleaned Data", $
            title=array+strcompress(xg[i]+1)+strcompress(yg[i]+1), $
            psym=5, ystyle=1, symsize=0.01, yrange=[-0.05,0.05], xstyle=1
    endfor

    device,/close
  endif

endif

print, "Doing fast maps up to component", ncomp

!p.multi=0

!p.thick=3.
!x.thick=3.
!y.thick=3.

cs = 1.5

for thecomp=0, ncomp-1 do begin

  device, filename=pre+'eigenmap'+strcompress(thecomp+1,/remove_all)+'.eps', $
          xsize=20, ysize=27, bits_per_pixel=8, /color, $
          /encapsulated

  print, thecomp, ' / ', ncomp
  compstr = strcompress(thecomp+1,/remove_all)

  ; write out FITS files
  ;decomp = newdata_pca*0

  ; project out single  components
  ;proj = newdata_pca*0 + valbad
  ;for i=0, n-1 do $
  ;  proj[xg[i],yg[i],*] = eof[thecomp,*]*pc[i,thecomp]

  ;fxwrite, pre+'proj'+compstr+'4.fits', head4, proj[0:nx-1,*,*]
  ;fxwrite, pre+'proj'+compstr+'8.fits', head8, proj[nx:2*nx-1,*,*]

  ; zero unused parts of cleaned data
  ;for i=0, nx-1 do begin
  ;  for j=0, ny-1 do begin
  ;    if mask4[i,j] ne 1 then begin
  ;      newdata_pca[i,j,*] = valbad
  ;    endif

  ;    if mask8[i,j] ne 1 then begin
  ;      newdata_pca[i+nx,j,*] = valbad
  ;    endif
  ;  endfor
  ;endfor

  ;fxwrite, pre+'cleaned'+strcompress(nrem,/remove_all)+'4.fits', head4, $
  ;         newdata_pca[0:nx-1,*,*]
  ;fxwrite, pre+'cleaned'+strcompress(nrem,/remove_all)+'8.fits', head8, $
  ;         newdata_pca[nx:2*nx-1,*,*]

  ; map eigenvalues for the given mode
  eigenmap = dblarr(nxb,nyb)
  eigenmap[xg,yg] = pc[*,thecomp];/pc[*,0] ; normalized by main mode
  eigenmap = eigenmap / 1d6 ; convert to uW from pW

  i = where(eigenmap ne 0, complement=j)

  m = mean(eigenmap[i])
  sig = stdev(eigenmap[i])
  maxval = m + 2*sig
  minval = m - 2*sig

  eigenmap[j] = 0 ; 1d30

  loadct, 0

  lambda = ['450','850']+micron

  cb = 0.595      ; bottom of component line plots
  ch = 0.145      ; height of component line plots
  csp = 0.065     ; space between component line plots
  width = 0.355
  xoff = 0.14


  labels=strarr(30)+' '

  loadct, 0

  pos = [xoff,cb+ch+csp,xoff+2*width,cb+csp+2*ch]

  ; time series
  plot, t, eof[thecomp,*], xstyle=1, ystyle=1, pos=pos, $
        ytitle='Eigenvector', charsize=cs, charthick=!p.thick, $
        xtitle='Time (s)', yrange=[-0.025,0.025]

  xyouts, (pos[0]+pos[2])/2., 0.965, 'Component'+strcompress(thecomp+1), $
    charsize=cs*2, charthick=!p.thick, /normal, alignment=0.5

  ; power spectra
  box = round(0.1/df)        ; width of boxcar in Hz / freq. step size
  ft = fft(eof[thecomp,*])
  p = (abs(ft)^2d)/df

  pos = [xoff,cb,xoff+2*width,cb+ch]

  plot, freq[1:nf-1], smooth(p[1:nf-1],box), xstyle=1, ystyle=1, pos=pos, $
    ytitle="PSD (Hz!u-1!N)", charsize=cs, charthick=!p.thick, $
    xtitle="Frequency (Hz)", /noerase, /ylog, /xlog, ytickformat='exponent', $
    xrange=[0.1,max(freq)]

  ; eigenmap
  for i=0, 1 do begin

    pos = [xoff+width*i,0.06,xoff+width*[i+1],0.53]

    plot,[0],[0], xrange=[0,nx], yrange=[0,ny], pos=pos, $
         xstyle=5, ystyle=5, $; title=lambda[i], $
         charthick=!p.thick, /noerase, charsize=cs

    loadct, 13
    tvlct, r, g, b, /get
    r[255] = 255 & g[255] = 255 & b[255] = 255 ; set 255th element to white
    tvlct, r, g, b

    map = eigenmap[i*nx:(i+1)*nx-1,*]
    ind = where(map eq 0)
    ;map[ind] = minval

    themap = bytscl( map, min=minval, max=maxval, top=254 )
    themap[ind] = 255

    tvscale, themap, /noint, minval=0, maxval=255, pos=pos

    cbar, minval, maxval, pos, 0.01, 0.02, $
          title='!7'+!gr.lambda+'!6'+' (!7'+!gr.mu+'!6W)', cs=cs

    loadct, 0
    axis, xaxis=0, xrange=[0,nx], xstyle=1, charthick=3., xtitle = 'Column', $
          charsize=cs
    axis, xaxis=1, xrange=[0,nx], xstyle=1, charthick=3., xtickname=labels
    if i eq 0 then begin
      axis, yaxis=0, yrange=[0,ny], xstyle=1, charthick=3., ytitle='Row', $
            charsize=cs
    endif else begin
      axis, yaxis=0, yrange=[0,ny], xstyle=1, charthick=3., ytickname=labels
    endelse

    axis, yaxis=1, yrange=[0,ny], ystyle=1, charthick=3., ytickname=labels

    xyouts, 1, 37, lambda[i], charsize=cs, charthick=!p.thick*5., color=255
    xyouts, 1, 37, lambda[i], charsize=cs, charthick=!p.thick, color=0

    if( i eq 0 ) then begin
      d = pc[0:n_elements(xg4)-1,thecomp]
    endif else begin
      d = pc[n_elements(xg4):n_elements(xg)-1,thecomp]
    endelse

    m = mean(d) / 1d6                 ; convert to uW
    sig = sqrt( mean(d^2d) ) / 1d6    ; convert to uW

    xyouts, 1,33, '!7'+!gr.lambda+'!6 = '+string(m,format='(F6.3)')+ ' !7' + $
      !gr.mu + '!6W', $
      charsize=cs, charthick=!p.thick*5., color=255

    plots, [1,2], 34.5*[1,1], thick=!p.thick*5, color=255

    xyouts, 1,31, '!7'+!gr.lambda+'!6!drms!n = '+ $
      string(sig,format='(F6.3)')+ ' !7' + $
      !gr.mu + '!6W', $
      charsize=cs, charthick=!p.thick*5., color=255

    xyouts, 1,33, '!7'+!gr.lambda+'!6 = '+string(m,format='(F6.3)')+ ' !7' + $
      !gr.mu + '!6W', $
      charsize=cs, charthick=!p.thick, color=0
    plots, [1,2], 34.5*[1,1]

    xyouts, 1,31, '!7'+!gr.lambda+'!6!drms!n = '+ $
      string(sig,format='(F6.3)')+ ' !7' + $
      !gr.mu + '!6W', $
      charsize=cs, charthick=!p.thick, color=0
  endfor

  device, /close

endfor



loadct,3

!p.multi=0
set_plot,'x'

!p.thick=1.
!x.thick=1.
!y.thick=1.

end
