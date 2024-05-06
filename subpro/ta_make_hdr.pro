FUNCTION ta_make_hdr,index0,tmp=tmp,remove_z=remove_z

nh=n_tags(index0)
tag=tag_names(index0)
res=strarr(nh+1)
FOR ih=0,nh-1 DO BEGIN
	;REPLACE STRENGE CHARACTER
	pos=strpos(tag[ih],'_D$')
	IF pos ge 0 THEN BEGIN
		a=tag[ih]
		strput,a,'_  ',pos
		tag[ih]=strcompress(a,/remove_all)
	ENDIF
	;-------------------------
	IF strlen(tag[ih]) eq 8 THEN BEGIN
		tag_str=tag[ih]
	ENDIF ELSE BEGIN
		tag_str=tag[ih]+strjoin(replicate(' ',8-strlen(tag[ih])))
	ENDELSE
	IF tag[ih] ne 'COMMENT' AND tag[ih] ne 'HISTORY' THEN BEGIN
		res[ih]=tag_str+'='+string(index0.(ih))+' /'
	ENDIF ELSE BEGIN
		res[ih]=tag_str+'  /'
	ENDELSE
ENDFOR
res[nh]='END     '

if strmid(res[0],0,6) ne 'SIMPLE' then res = ['SIMPLE  =      T /',res]

if keyword_set(remove_z) then begin
	pos = where(strmid(res,0,1) ne 'Z')
	res = res[pos]
endif

RETURN,res
END
