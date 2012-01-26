FUNCTION Exponent, axis, index, number

  IF number EQ 0 THEN RETURN, '0' ;; Special case
  
  ex = String(number, '(e7.0)') 
  pt = StrPos(ex, '.')
  
  first = StrMid(ex, pt-1, 1) 
  sign = StrMid(ex, pt+2, 1)
  exponent = StrMid(ex, pt+3, 100)

  ;; Shave off leading zero in exponent

  WHILE StrMid(exponent, 0, 1) EQ '0' DO exponent = StrMid(exponent, 1, 100) 
  
  if exponent eq '' then exponent = '0'

  IF sign EQ '-' THEN   RETURN, '10!U' + sign + exponent $
     ELSE               RETURN, '10!U' + exponent

  ;IF sign EQ '-' THEN   RETURN, first + 'x10!U' + sign + exponent $
  ;   ELSE               RETURN, first + 'x10!U' + exponent
  END
