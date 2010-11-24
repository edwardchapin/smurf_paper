;------------------------------------------------------------------------------
; Parse a .tst file produced by JCMTSTATE2CAT (JCMT state information
; specifically for SCUBA-2 data files
;------------------------------------------------------------------------------

function scuba2_readstate, filename

  tab = string(9b)

  ; skip over header
  openr, 1, filename
  ln = ''
  for i=0, 2 do readf, 1, ln

  ; counts lines
  n = 0l
  while eof(1) eq 0 do begin
    n = n + 1
    readf, 1, ln
  endwhile
  close,1

  ; read in data

  Id = lonarr(n)
  RA = dblarr(n)
  DEC = dblarr(n)
  DRA = dblarr(n)
  DDEC = dblarr(n)
  AZ = dblarr(n)
  EL = dblarr(n)
  DAZ = dblarr(n)
  DEL = dblarr(n)
  WVMTAU = dblarr(n)
  FTS_POS = dblarr(n)
  JOS_DRCONTROL = lonarr(n)
  POL_ANG = dblarr(n)
  RTS_END = dblarr(n)
  RTS_NUM = lonarr(n)
  SC2_BIAS = dblarr(n)
  SC2_FPUTEMP = dblarr(n)
  SC2_HEAT = dblarr(n)
  SC2_MIXTEMP = dblarr(n)
  SMU_AZ_CHOP_X = strarr(n)
  SMU_AZ_CHOP_Y = strarr(n)
  SMU_AZ_JIG_X = strarr(n)
  SMU_AZ_JIG_Y = strarr(n)
  SMU_CHOP_PHASE = strarr(n)
  SMU_JIG_INDEX = strarr(n)
  SMU_TR_CHOP_X = strarr(n)
  SMU_TR_CHOP_Y = strarr(n)
  SMU_TR_JIG_X = strarr(n)
  SMU_TR_JIG_Y = strarr(n)
  TCS_AIRMASS = dblarr(n)
  TCS_AZ_AC1 = dblarr(n)
  TCS_AZ_AC2 = dblarr(n)
  TCS_AZ_ANG = dblarr(n)
  TCS_AZ_BC1 = dblarr(n)
  TCS_AZ_BC2 = dblarr(n)
  TCS_AZ_DC1 = dblarr(n)
  TCS_AZ_DC2 = dblarr(n)
  TCS_BEAM = strarr(n)
  TCS_DM_ABS = dblarr(n)
  TCS_DM_REL = dblarr(n)
  TCS_EN_DC1 = dblarr(n)
  TCS_EN_DC2 = dblarr(n)
  TCS_INDEX = strarr(n)
  TCS_PERCENT_CMP = strarr(n)
  TCS_SOURCE = strarr(n)
  TCS_TAI = dblarr(n)
  TCS_TR_AC1 = dblarr(n)
  TCS_TR_AC2 = dblarr(n)
  TCS_TR_ANG = dblarr(n)
  TCS_TR_BC1 = dblarr(n)
  TCS_TR_BC2 = dblarr(n)
  TCS_TR_DC1 = dblarr(n)
  TCS_TR_DC2 = dblarr(n)
  TCS_TR_SYS = strarr(n)
  WVM_T12 = dblarr(n)
  WVM_T42 = dblarr(n)
  WVM_T78 = dblarr(n)
  WVM_TIME = dblarr(n)

  ; skip over header
  openr, 1, filename
  ln = ''
  for i=0, 2 do readf, 1, ln

  for i=0l, n-1 do begin
    readf, 1, ln
    words = strsplit(ln, tab, /extract)

    if n_elements(words) eq 58 then begin
      Id[i] = words[0]
      RA[i] = words[1]
      DEC[i] = words[2]
      DRA[i] = words[3]
      DDEC[i] = words[4]
      AZ[i] = words[5]
      EL[i] = words[6]
      DAZ[i] = words[7]
      DEL[i] = words[8]
      WVMTAU[i] = words[9]
      FTS_POS[i] = words[10]
      JOS_DRCONTROL[i] = words[11]
      POL_ANG[i] = words[12]
      RTS_END[i] = words[13]
      RTS_NUM[i] = words[14]
      SC2_BIAS[i] = words[15]
      SC2_FPUTEMP[i] = words[16]
      SC2_HEAT[i] = words[17]
      SC2_MIXTEMP[i] = words[18]
      SMU_AZ_CHOP_X[i] = words[19]
      SMU_AZ_CHOP_Y[i] = words[20]
      SMU_AZ_JIG_X[i] = words[21]
      SMU_AZ_JIG_Y[i] = words[22]
      SMU_CHOP_PHASE[i] = words[23]
      SMU_JIG_INDEX[i] = words[24]
      SMU_TR_CHOP_X[i] = words[25]
      SMU_TR_CHOP_Y[i] = words[26]
      SMU_TR_JIG_X[i] = words[27]
      SMU_TR_JIG_Y[i] = words[28]
      TCS_AIRMASS[i] = words[29]
      TCS_AZ_AC1[i] = words[30]
      TCS_AZ_AC2[i] = words[31]
      TCS_AZ_ANG[i] = words[32]
      TCS_AZ_BC1[i] = words[33]
      TCS_AZ_BC2[i] = words[34]
      TCS_AZ_DC1[i] = words[35]
      TCS_AZ_DC2[i] = words[36]
      TCS_BEAM[i] = words[37]
      TCS_DM_ABS[i] = words[38]
      TCS_DM_REL[i] = words[39]
      TCS_EN_DC1[i] = words[40]
      TCS_EN_DC2[i] = words[41]
      TCS_INDEX[i] = words[42]
      TCS_PERCENT_CMP[i] = words[43]
      TCS_SOURCE[i] = words[44]
      TCS_TAI[i] = words[45]
      TCS_TR_AC1[i] = words[46]
      TCS_TR_AC2[i] = words[47]
      TCS_TR_ANG[i] = words[48]
      TCS_TR_BC1[i] = words[49]
      TCS_TR_BC2[i] = words[50]
      TCS_TR_DC1[i] = words[51]
      TCS_TR_DC2[i] = words[52]
      TCS_TR_SYS[i] = words[53]
      WVM_T12[i] = words[54]
      WVM_T42[i] = words[55]
      WVM_T78[i] = words[56]
      WVM_TIME[i] = words[57]
    endif else begin
      print, "Incorrect number of values in line", i
      JOS_DRCONTROL[i] = -1
      Id[i] = -1
    endelse
  endfor

  close,1

  retval = {$
           Id:Id, $
           RA:RA, $
           DEC:DEC, $
           DRA:DRA, $
           DDEC:DDEC, $
           AZ:AZ, $
           EL:EL, $
           DAZ:DAZ, $
           DEL:DEL, $
           WVMTAU:WVMTAU, $
           FTS_POS:FTS_POS, $
           JOS_DRCONTROL:JOS_DRCONTROL, $
           POL_ANG:POL_ANG, $
           RTS_END:RTS_END, $
           RTS_NUM:RTS_NUM, $
           SC2_BIAS:SC2_BIAS, $
           SC2_FPUTEMP:SC2_FPUTEMP, $
           SC2_HEAT:SC2_HEAT, $
           SC2_MIXTEMP:SC2_MIXTEMP, $
           SMU_AZ_CHOP_X:SMU_AZ_CHOP_X, $
           SMU_AZ_CHOP_Y:SMU_AZ_CHOP_Y, $
           SMU_AZ_JIG_X:SMU_AZ_JIG_X, $
           SMU_AZ_JIG_Y:SMU_AZ_JIG_Y, $
           SMU_CHOP_PHASE:SMU_CHOP_PHASE, $
           SMU_JIG_INDEX:SMU_JIG_INDEX, $
           SMU_TR_CHOP_X:SMU_TR_CHOP_X, $
           SMU_TR_CHOP_Y:SMU_TR_CHOP_Y, $
           SMU_TR_JIG_X:SMU_TR_JIG_X, $
           SMU_TR_JIG_Y:SMU_TR_JIG_Y, $
           TCS_AIRMASS:TCS_AIRMASS, $
           TCS_AZ_AC1:TCS_AZ_AC1, $
           TCS_AZ_AC2:TCS_AZ_AC2, $
           TCS_AZ_ANG:TCS_AZ_ANG, $
           TCS_AZ_BC1:TCS_AZ_BC1, $
           TCS_AZ_BC2:TCS_AZ_BC2, $
           TCS_AZ_DC1:TCS_AZ_DC1, $
           TCS_AZ_DC2:TCS_AZ_DC2, $
           TCS_BEAM:TCS_BEAM, $
           TCS_DM_ABS:TCS_DM_ABS, $
           TCS_DM_REL:TCS_DM_REL, $
           TCS_EN_DC1:TCS_EN_DC1, $
           TCS_EN_DC2:TCS_EN_DC2, $
           TCS_INDEX:TCS_INDEX, $
           TCS_PERCENT_CMP:TCS_PERCENT_CMP, $
           TCS_SOURCE:TCS_SOURCE, $
           TCS_TAI:TCS_TAI, $
           TCS_TR_AC1:TCS_TR_AC1, $
           TCS_TR_AC2:TCS_TR_AC2, $
           TCS_TR_ANG:TCS_TR_ANG, $
           TCS_TR_BC1:TCS_TR_BC1, $
           TCS_TR_BC2:TCS_TR_BC2, $
           TCS_TR_DC1:TCS_TR_DC1, $
           TCS_TR_DC2:TCS_TR_DC2, $
           TCS_TR_SYS:TCS_TR_SYS, $
           WVM_T12:WVM_T12, $
           WVM_T42:WVM_T42, $
           WVM_T78:WVM_T78, $
           WVM_TIME:WVM_TIME }

  return, retval

end
