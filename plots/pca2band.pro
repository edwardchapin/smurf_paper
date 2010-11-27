;------------------------------------------------------------------------------
; PCA some bolo data at 450 and 850 simultaneously
;------------------------------------------------------------------------------

pre='./2band_'
datadir='../data/'
obs = '20100313_00029'

ncomp = 32    ; single component plots (up to this value)
nrem = 1      ; number of components to remove
dobigplot = 1 ; do the big plots?

if 1 then begin

  ; load bolo data and mask
  fxread, datadir+'s8d'+obs+'_con_clean.fits', data8, head8
  fxread, datadir+'s4a'+obs+'_con_clean.fits', data4, head4

  fxread, datadir+'mask_s8d.fits', mask8, mhead4
  fxread, datadir+'mask_s4a.fits', mask4, mhead4

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
  fact = 1;8

  ind4 = where(mask4)
  ind8 = where(mask8)

  themin = min([n_elements(ind4),n_elements(ind8)])

  m = mask4*0
  for i=0, themin/fact-1 do m[ind4[i*fact]] = 1
  mask4 = m

  m = mask8*0
  for i=0, themin/fact-1 do m[ind8[i*fact]] = 1
  mask8 = m

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
        pc, $
        filename=pre+'data.sav'
endif else begin
  print, 'Restoring PCA data...'
  restore,pre+'data.sav'
endelse

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
  print, "Doing large plots..."

  loadct, 0

  !p.thick=1.
  !x.thick=1.
  !y.thick=1.

  ;sx = 3
  sx = 4
  sy = 4

  !p.multi=[0,sx,sy]
  device, filename=pre+'bolos.ps', xsize=27.94, ysize=21.59

  for i=0, n-1 do begin
    if xg[i] ge nx then array='s8d' $
    else array='s4a'

    plot, newdata[xg[i],yg[i],*], xtitle="Sample #", ytitle="Raw Data", $
          title=array+strcompress(xg[i]+1)+strcompress(yg[i]+1), $
          psym=5, ystyle=1, symsize=0.01, yrange=[-0.2,0.2], xstyle=1
  endfor

  device,/close

  !p.multi=[0,1,1]
  !p.multi=[0,sx,sy]
  device, filename=pre+'modes.ps'

  for i=0, n-1 do $
    plot, eof[i,*], xtitle="Sample #", ytitle="Normalized Signal", $
          title='Component '+strcompress(i+1), $
          psym=5, ystyle=1, symsize=0.01, yrange=[-0.05,0.05], xstyle=1

  device,/close

  !p.multi=[0,1,1]
  !p.multi=[0,sx,sy]
  device, filename=pre+'bolos_clean.ps'

  for i=0, n-1 do begin
    if xg[i] ge nx then array='s8d' $
    else array='s4a'

    plot, newdata_pca[xg[i],yg[i],*], xtitle="Sample #", $
          ytitle="Raw Cleaned Data", $
          title=array+strcompress(xg[i]+1)+strcompress(yg[i]+1), $
          psym=5, ystyle=1, symsize=0.01, yrange=[-0.05,0.05], xstyle=1
  endfor

  device,/close

endif

print, "Doing fast maps up to component", ncomp

!p.multi=0

device, filename=pre+'eigenmap.ps', $
        xsize=2.*38./3., ysize=40./3., bits_per_pixel=8, $
        /color

!p.thick=3.
!x.thick=3.
!y.thick=3.

for thecomp=0, ncomp-1 do begin
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
  i = where(eigenmap ne 0, complement=j)

  m = mean(eigenmap[i])
  sig = stdev(eigenmap[i])
  maxval = m + 2*sig
  minval = m - 2*sig

  ;maxval = max(eigenmap[i])
  ;minval = min(eigenmap[i])

  eigenmap[j] = 0 ; 1d30

  pos =[0.1,0.1,0.75,0.95]
  plot,[0],[0], xrange=[0,nxb], yrange=[0,nyb], pos=pos, $
       xstyle=5, ystyle=5, title='Eigenmap Mode = '+compstr, $
       charthick=!p.thick

  loadct,13
  tvscale, eigenmap, /noint, minval=minval, maxval=maxval, pos=pos

  labels=strarr(30)+' '

  axis, xaxis=0, xrange=[0,nxb], xstyle=1, charthick=3.
  axis, xaxis=1, xrange=[0,nxb], xstyle=1, charthick=3., xtickname=labels
  axis, yaxis=0, yrange=[0,nyb], xstyle=1, charthick=3.
  axis, yaxis=1, yrange=[0,nyb], ystyle=1, charthick=3., ytickname=labels

  cbar, minval, maxval, pos, 0.05, 0.05

endfor

device, /close

loadct,3

!p.multi=0
set_plot,'x'

!p.thick=1.
!x.thick=1.
!y.thick=1.

end
