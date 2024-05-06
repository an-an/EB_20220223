;+
; NAME       : ta_rm_fringe_pca.pro (function)
; PURPOSE :
; 	remove fringe pattern by using PCA 
; CATEGORY :
; CALLING SEQUENCE :
;	 p = ta_rm_fringe_pca(arr,index=index,/show)
; INPUTS :
;	ii	-- time variation of spectra ;(#slit position,#wavelebgth,#time)
; OUTPUT :
;	spectra without fringe and removal fringe
; OPTIONAL INPUT PARAMETERS : 
;	index   -- index of eigenfeatures used to reconstruct to fringe 
; KEYWORD PARAMETERS :
; 	show   -- show the results or not
; PROCEDURE:
; REFERENCE:
;	Casini, R. 2012, ApJ, 756, 194 
; MODIFICATION HISTORY :
;        t.a. '14/01/28

function ta_rm_fringe_pca,ii,index=index,show=show


ss=size(ii)
nx=ss[1]
ny=ss[2]
nt=ss[3]
avg=total(ii,3)/float(nt)
covar=fltarr(ny,ny)
for i=0,nt-1 do covar=covar+(ii[*,*,i]-avg)##transpose(ii[*,*,i]-avg)
covar=covar/float(nt)

svdc,covar,sigma,U,V,/double,itmax=100
;svdc,covar,sigma,U,V,/double;,itmax=100
u=u[where(sigma ne 0),*]
sigma=sigma[where(sigma ne 0,nsigma)]

eval=fltarr(ny);neigenprofs)
eprof=fltarr(ny,ny)
for j=0,nsigma-1 do begin
    eval(j)=sigma(j)
    eprof(*,j)=U(j,*)
endfor

if keyword_set(show) then begin
	window,4
	plot,eval,/ylog
endif
	
window,2,xs=600,ys=700
!p.multi=[0,2,10]
for i=0,19 do plot,eprof[*,i],title=string(i+1),charsize=1.5
!p.multi=0

ans=''
if keyword_set(index) eq 0 then begin
	print,'write index of eigenfeatures used to reconstrcut to fringe'
	print,'ex) 1,2,10,12'
	read,':',ans
	pos=strpos(ans,',')
	key=1	
	while pos ne -1 do begin
		if key eq 1 then begin
			key=0
			jj=fix(strmid(ans,0,pos))
		endif else jj=[jj,fix(strmid(ans,0,pos))]
		ans=strmid(ans,pos+1,strlen(ans)-pos+1)
		pos=strpos(ans,',')
	endwhile 
	if key eq 1 then begin
		key=0
		jj=fix(ans)
	endif else jj=[jj,fix(ans)]
endif else jj=index

jj=jj-1
new=fltarr(nx,ny,nt)
fri=fltarr(nx,ny,nt)
ss=size(jj)
if ss[0] eq 0 then nj=1 else nj=(size(jj))[1]
for i=0,nt-1 do begin
	vv=ii[*,*,i]#eprof
	for j=0,nj-1 do fri[*,*,i]=fri[*,*,i]+$
			vv[*,jj[j]]#transpose(eprof(*,jj[j]))
	new[*,*,i]=ii[*,*,i]-fri[*,*,i]

	if keyword_set(show) and i eq 0 then begin
		window,0,xs=nx,ys=3*ny
		tvscl,ii[*,*,i],0;<2e4>1e3,0
		tvscl,fri[*,*,i],1
		tvscl,new[*,*,i],2;<2e4>1e3,2
	endif
endfor

return,new

END
