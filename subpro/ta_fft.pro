pro ta_fft,y,dx,k,a,hwin=hwin

n=(size(y))[1]
if keyword_set(hwin) then hwindow=hanning(n,/double) else hwindow=fltarr(n)>1.<1.
a=fft(hwindow*y)
if max(imaginary(y)) eq 0 and min(imaginary(y)) eq 0 then begin
	if n mod 2 eq 0 then begin
		k=[findgen(n/2+1),((n/2-1)-findgen(n/2-1))]/float(n*dx)
	endif else begin
		k=[findgen((n-1)/2+1),(((n-1)/2)-findgen((n-1)/2))]/float(n*dx)
	endelse
endif else begin
	if n mod 2 eq 0 then begin
		k=[findgen(n/2+1),-((n/2-1)-findgen(n/2-1))]/float(n*dx)
	endif else begin
		k=[findgen((n-1)/2+1),-(((n-1)/2)-findgen((n-1)/2))]/float(n*dx)
	endelse
endelse 

end
