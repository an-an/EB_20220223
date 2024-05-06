@./subpro/ta_fft.pro
@./subpro/ta_make_hdr.pro
@./subpro/ta_atlas.pro
@./subpro/ta_rm_fringe_pca.pro
@./subpro/ta_cthph2thph.pro

PRO set_parameters
COMMON ANALYSIS, data_dir, wavelength, index_for_map,  $
		 bomb_data, bomb_x, bomb_y, wvl_cah3, wvl_cah2v, wvl_cah1v, wvl_fev, wvl_cairv, wvl_hep_v, $
		 cal01b_param,cal01_param,hairline,ref_pixel,	$
		 nicole_param,bomb_pixel,imax_param, observation, binx

;dir0 = '/data/tanan/dkist/20220223_ocp1/20230425_release/pid_1_50/'
dir0  = './demodata/'
data_dir={	dir397:	[dir0+'BOLDQ'], $
		dir630:	[dir0+'BPKJO'], $
		dir854: [dir0+'BQKNO']}

hairline={	pos397 : [122.,1943.]-1,	$
		pos630 : [344.,1843.]-1,	$
		pos854 : [111.,2406.]-3}
	
observation0 = {	date: '2022-02-23',	$
			hxy : [-720.,310.]	$
		}
observation = replicate(observation0, 1)
;-------------------	
	
binx = {	bin397: 8,	$
		bin630:	6,	$
		bin854: 10}

;Wavelength in angstrom
wavelength={	wvl397:	3964.18d + dindgen(1000)*7.7d   /1000d	, $
		wvl630:	6294.95d + dindgen(1000)*12.85d /1000d	, $
		wvl854:	8531.82d + dindgen(1000)*18.82d /1000d	, $
		; wavelength for specific purposes
		wvl_cah3: 3968.4535d			, $
		wvl_cah2v:3968.2918d			, $
		wvl_cah1v:3968.0531d				, $
		wvl_cahcon:3965.00d			, $
		wvl_fecon :6297.00d			, $
		wvl_caicon:8535.00d			, $
		wvl_fev:  6301.508d + 0.1d			, $
		wvl_cairv:8542.144d + 0.5d			, $
		wvl_h1e  :3970.05				, $
		;wvl_h1e_r:[3969.85,3970.18,3970.50,3970.9]			, $
		;wvl_h1e_r1:[3969.5,3969.55,3969.8,3970.2,3970.8,3970.9]			, $ ; for make_data_06
		wvl_h1e_r1:[3969.5,3969.55,3969.8,3970.2,3970.58,3970.6]			, $ ; for make_data_06
		wvl_h1e_r:[3969.8,3970.18,3970.45,3970.9]			, $
		wvl_hep_v:[3969.97d, 3970.17d]	}

cal01b_param = { $
		;con397:	[80+indgen(125-80),280+indgen(290-280),360+indgen(400-360),  $
		;	450+indgen(475-450),630+indgen(640-630),695+indgen(700-695),$
		;	740+indgen(750-740),850+indgen(870-850),970+indgen(980-970)],  $
		con397: [65+indgen(80-65), 107+indgen(123-107),269+indgen(291-269),358+indgen(401-358), 436+indgen(470-436),853+indgen(865-853), 951+indgen(961-951)],  $
		con630:	[40+indgen(60-40),150+indgen(190-150),280+indgen(300-250),  $
			390+indgen(400-390),750+indgen(800-750),900+indgen(1000-900)],  $
		con854:	[0+indgen(200-0),270+indgen(300-270),700+indgen(830-700),950+indgen(1000-950)],   $
		smooth_width	: 100,	$
		poly_order	: 1	$
		}

con397 = [75+indgen(80-75), 107+indgen(123-107),275+indgen(291-275),358+indgen(401-358),  $
	688+indgen(7),830+indgen(865-830), 955+indgen(970-955)]
con630 = [10+indgen(20),120+indgen(10),160+indgen(10),260+indgen(30), $
       355+indgen(15),700+indgen(80),900+indgen(50)]	
con854 = [0+indgen(190),260+indgen(300-260),360+indgen(20),430+indgen(20),700+indgen(830-700),940+indgen(980-940)]
cal01_param = { $
		con397: con397,		$
		con630:	con630,		$
		con854:	con854,		$
		err397: con397>1.<1.,	$ 
		err630: con630>1.<1.,	$ 
		err854: con854>1.<1.,	$ 
		polyorder397: 10,	$
		polyorder630: 5,	$
		polyorder854: 5,	$
		smooth_width:50,	$
		plot_idir: 1,	$
		plot_ifile: 111,	$
		plot_ix: 1020,		$
		smooth_offset: 1,	$
		smooth_factor:1e10	$
		}
cal01_param.err397[where(cal01_param.con397 le 500)] = 5.
cal01_param.err397[where(cal01_param.con397 ge 720 and cal01_param.con397 le 900)] = .5





;INDEX_FOR_MAP
index0={        SINPLE:0,       $
                BITPIX:18l,     $
                NAXIS:2l,       $
                NAXIS1:400l,   $
                NAXIS2:2560l,   $
                EXTEND:0,       $
                DATE_OBS:'23-Feb-2022',    $
                FILENAME:' ',    $
                STATUS:'',   $
                INSTRUME:'ViSP',$
                TELESCOP:'DKIST',$
                XCEN:0d,        $
                YCEN:0d,        $
                CDELT1:0.2d,      $
                CDELT2:0.0d,      $
                CUNIT1:'arcsec',        $
                CUNIT2:'arcsec',        $
                CROTA2:0.,      $
                COMMENT:strarr(5), $
                HISTORY:strarr(5) $
        }
cdelt2_0 = 0.0245
for i=0,2 do begin
        case i of
                0:begin
			cdelt2 = cdelt2_0
                        ycen=0.
                        imax=4e5
                end
                1:begin
			cdelt2 = cdelt2_0*(hairline.pos397[1]-hairline.pos397[0])/(hairline.pos630[1]-hairline.pos630[0]);   0.0298
                        ycen=-0.40
                        imax=1e6
                end
                2:begin
			cdelt2 = cdelt2_0*(hairline.pos397[1]-hairline.pos397[0])/(hairline.pos854[1]-hairline.pos854[0]);  0.0194
                        ycen=-5.63
                        imax=5e5
                end
        endcase
        index=index0

        index.cdelt2=cdelt2
        index.ycen=ycen

        case i of
                0:index_nosun_397 = index
                1:index_nosun_630 = index
                2:index_nosun_854 = index
        endcase
endfor
index_for_map = {	index_nosun_397:	index_nosun_397, $
			index_nosun_630:	index_nosun_630, $
			index_nosun_854:	index_nosun_854 $
			}





;The position of Bombs in no sun map cordination in arcsec
bomb_data = [1,1]
bomb_x	  = [  $
		[-6.5+[-2.,2.]],  $
		[-5.3+[-2.,2.]]	  $
		]
bomb_y	  = [  $
		[-17.7+[-2.,2.]],  $
		[-21.7+[-2.,2.]]  $
		]

ref_pixel = {  	data    : 0,    $
                ifile   : 130,   $
                ;ix      : [500,550],        $
                ix      : [300,310],        $
                comment : ''           $
                }


;The position of Bombs in pixel
bomb_pixel = {	data	: 1,	$
		ifile	: 64,	$
		ix	: 1870+[-20,20],	$
		ix_ref	: 1750+[-20,20],	$
		ix_ref1	: 1750+[-20,20],	$   ; for reference of fitting
		wvl_lc  : 3970.15,		$
		comment : 'Zeeman V?'		$
		}
bomb_pixel = replicate(bomb_pixel,1)
bomb_pixel[0].data  = 0
bomb_pixel[0].ifile = 111
bomb_pixel[0].ix    = 1020+[-20,20]
bomb_pixel[0].ix_ref= 1100+[-20,20]
bomb_pixel[0].ix_ref1= [50,60]
bomb_pixel[0].wvl_lc= 3970.05
bomb_pixel[0].comment= 'Good'
;---

nicole_param = {	$
	wvl397_range	: [3966.,3971.],	$ ; middle range
	wvl630_range	: [6300.,6304.],	$ ; middle
	wvl854_range	: [8539.,8545.],	$ ; middle range
; Weak Field Approximation
	wfa397_range	: [6301.3,6301.6],	$ ; wavelength range for WFA
	wfa630_range	: [6301.3,6301.6],	$ ; wavelength range for WFA
	wfa854_range	: [8541.3,8542.3],	$ ; wavelength range for WFA
	mu		: cos(asin(sqrt(observation.hxy[0]^2 + observation.hxy[1]^2)/get_rb0p(observation.date,/radius,/quiet)))	$
	}

imax_param={	qs_ifile_range	: [390,399],	$
		qs_ix_range	: [0,200],	$
		pos397		: [70+indgen(120-70),260+indgen(285-260),350+indgen(400-350),  $
					440+indgen(640-440),675+indgen(695-675),		$
					850+indgen(870-850),945+indgen(965-945)],		$
		pos630          : [350+indgen(380-350),435+indgen(450-435),600+indgen(630-600)],  $ ; con630[70:150]
		;pos854		: [findgen(80),findgen(80)+211-80],	$
		pos854		: [findgen(180),255+indgen(50),350+indgen(600)],	$
		;wvl630_wr	: [6301.,6301.8],  $
		wvl630_wr1	: [6301.,6301.9],  $
		wvl630_wr2	: [6301.15,6302.05],  $
		imax397		: 0.,		$
		imax630		: 0.,		$
		imax854		: 0.		$
	}


index0={        SINPLE:0,       $
                BITPIX:18l,     $
                NAXIS:2l,       $
                NAXIS1:2560l,   $
                NAXIS2:400l,   $
                EXTEND:0,       $
                DATE_OBS:'23-Feb-2022',    $
                FILENAME:' ',    $
                STATUS:'',   $
                INSTRUME:'ViSP',$
                TELESCOP:'DKIST',$
                XCEN:0d,        $
                YCEN:0d,        $
                CDELT1:0d,      $
                CDELT2:0.2,      $
                CUNIT1:'arcsec',        $
                CUNIT2:'arcsec',        $
                CROTA2:0.,      $
                COMMENT:strarr(5), $
                HISTORY:strarr(5) $
        }
for i=0,2 do begin
	case i of
                0:begin
			factor=1.21859
                        xcen=0.
                        imax=4e5
                end
                1:begin
                        factor=1.
                        xcen=0.
                        imax=1e6
                end
                2:begin
                        factor=1.53357
                        xcen=-5.
                        imax=5e5
                end
        endcase
	index=index0

        index.cdelt1=0.03/factor
        index.xcen=xcen

	case i of
		0:index_nosun_397 = index
		1:index_nosun_630 = index
		2:index_nosun_854 = index
	endcase
endfor

END
;-----------------------------------------------
;-----------------------------------------------
PRO  WVL_CAL, data_dir=data_dir, wavelength0=wavelength0, image_file=image_file, ref_wvl=ref_wvl;, ref_pix=ref_pix
;COMMON ANALYSIS, dir397, dir630, dir854, imgdir, wvl397, wvl630, wvl854

serch_width = 10 ; pix


if 1 then begin
file1=file_search(data_dir,'*_I_*.fits',count=nf)
for i=0,nf-1 do begin
        read_sdo,file1[i],index,data1, /use_shared_lib
	if i eq 0 then begin
		sp1=data1
		sp2=fltarr(index[0].naxis2,nf)
	endif
	sp2[*,i]=rebin(data1,1,index[0].naxis2)
endfor
sp3=rebin(sp2,index[0].naxis2,1)
save,sp1,sp2,sp3,index,file='/data/tanan/idl/tmp.sav'
endif else restore,'/data/tanan/idl/tmp.sav'


nref=n_elements(ref_wvl)
pix = fltarr(nref)
for iref=0,nref-1 do begin
;	tmp=min(sp3[ref_pix[iref]-serch_width:ref_pix[iref]+serch_width],pos)
;	pix[iref] = ref_pix[iref]-serch_width + pos
	tmp=min(abs(wavelength0-ref_wvl[iref]),pos)
	pix[iref] = pos
endfor
;coe=poly_fit(pix,ref_wvl,1)
;wvl = coe[1]*findgen(index[0].naxis2)+coe[0]
wvl = wavelength0


xs=0.15
xe=0.95
ys=0.07
ye=0.95
ydd=0.01
chs=3

yd= (ye-ys-2.*ydd)/3.
brank=replicate(' ',10)

wdef,0,600,1000
!p.multi=[0,1,3]
loadct,0

iy=2
loadct,0
plot_image,transpose(sp1), $
  color=0,background=255,  $
  charsize=chs,xtickname=brank,ytitle='Slit (pix)',title='Wavelength calibration',  $
  noerase=0,/norm,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd]
set_line_color
for iref=0,nref-1 do oplot,replicate(pix[iref],2),[0,index[0].naxis1],color=3

iy=1
loadct,0
plot_image,sp2,  $
  color=0,background=255,  $
  charsize=chs,xtickname=brank,ytitle='Scan (pix)',  $
  noerase=0,/norm,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd]
set_line_color
for iref=0,nref-1 do oplot,replicate(pix[iref],2),!y.crange,color=3

iy=0  
set_line_color
plot,wvl,sp3,/xstyle,  $
  color=0,background=1,  $
  xtitle='Wavelength (A)',ytitle='Intensity (ADU)',  $
  charsize=chs,noerase=0,/norm,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd]
for iref=0,nref-1 do oplot,replicate(ref_wvl[iref],2),!y.crange,color=3

write_png,image_file,tvrd(/true)
loadct,0
!p.multi=0


END
;-----------------------------------------------
FUNCTION  WVL_CAL_1, data_dir=data_dir, wavelength0=wavelength0, image_file=image_file, ref_wvl=ref_wvl, ref_pix=ref_pix
;COMMON ANALYSIS, dir397, dir630, dir854, imgdir, wvl397, wvl630, wvl854

serch_width = 10 ; pix

if 1 then begin
file1=file_search(data_dir,'*_I_*.fits',count=nf)
for i=0,nf-1 do begin
        read_sdo,file1[i],index,data1, /use_shared_lib
        if i eq 0 then begin
                sp1=data1
                sp2=fltarr(index[0].naxis2,nf)
        endif
        sp2[*,i]=rebin(data1,1,index[0].naxis2)
endfor
sp3=rebin(sp2,index[0].naxis2,1)
save,sp1,sp2,sp3,index,file='/data/tanan/idl/tmp.sav'
endif else restore,'/data/tanan/idl/tmp.sav',/verb


nref=n_elements(ref_wvl)
pix = fltarr(nref)
for iref=0,nref-1 do begin
       tmp=min(sp3[ref_pix[iref]-serch_width:ref_pix[iref]+serch_width],pos)
       pix[iref] = ref_pix[iref]-serch_width + pos
endfor
coe=poly_fit(pix,ref_wvl,1)
wvl = coe[1]*findgen(index[0].naxis2)+coe[0]

;wvl = 12.85e-3*findgen(index[0].naxis2)+coe[0]

;atlas=ta_atlas(min(wvl),max(wvl),coe[1],wl=wl);   ,liege=liege,spot=spot




xs=0.15
xe=0.95
ys=0.07
ye=0.95
ydd=0.01
chs=3

yd= (ye-ys-2.*ydd)/3.
brank=replicate(' ',10)

wdef,0,600,1000
!p.multi=[0,1,3]
loadct,0

iy=2
loadct,0
plot_image,transpose(sp1), $
  color=0,background=255,  $
  charsize=chs,xtickname=brank,ytitle='Slit (pix)',title='Wavelength calibration',  $
  noerase=0,/norm,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd]
set_line_color
for iref=0,nref-1 do oplot,replicate(pix[iref],2),[0,index[0].naxis1],color=3

iy=1
loadct,0
plot_image,sp2,  $
  color=0,background=255,  $
  charsize=chs,xtickname=brank,ytitle='Scan (pix)',  $
  noerase=0,/norm,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd]
set_line_color
for iref=0,nref-1 do oplot,replicate(pix[iref],2),!y.crange,color=3

iy=0
set_line_color
plot,wvl,sp3,/xstyle,  $
  color=0,background=1,  $
  xtitle='Wavelength (A)',ytitle='Intensity (ADU)',  $
  charsize=chs,noerase=0,/norm,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd]
for iref=0,nref-1 do oplot,replicate(ref_wvl[iref],2),!y.crange,color=3


write_png,image_file,tvrd(/true)
loadct,0
!p.multi=0


RETURN, wvl

END
;-----------------------------------------------
FUNCTION pix2arcsec_nosun, xpix, ypix, index, atm_disp = atm_disp
; xpix in slit direction
; ypix in scan direction
; return[0] in scan direction
; return[1] in slit direction

if not keyword_set(atm_disp) then atm_disp = 0.

yarc = (xpix + atm_disp - index.naxis2*0.5)*index.cdelt2 + index.ycen + 0.2
xarc = (ypix - index.naxis1*0.5)*index.cdelt1 + index.xcen + 0.2

RETURN, [[xarc],[yarc]]

END
;-----------------------------------------------
FUNCTION arcsec2pix_nosun, xarcsec, yarcsec, index, atm_disp = atm_disp
; xarcsec in scan direction
; yarcsec in slit direction



if not keyword_set(atm_disp) then atm_disp = 0.
xpix = (yarcsec - index.ycen - 0.2) / index.cdelt2 + index.naxis2*0.5 - atm_disp
ypix = (xarcsec - index.xcen - 0.2) / index.cdelt1 + index.naxis1*0.5


RETURN, [[xpix],[ypix]]

END
;-----------------------------------;
;-----------------------------------------------
PRO CAL01,image_dir=image_dir, image_file=image_file, name_tail=name_tail, result_check_crosstalk_01=result_check_crosstalk_01, save_file=save_file
COMMON ANALYSIS

; offset_ifile : offset in ifile to make artpol (artification polarization spectra)

plot_cah_special = 0

if not keyword_set(name_tail) then outname_tail='_cal01' else outname_tail=name_tail 
if is_dir(image_dir) eq 0 then spawn,'mkdir '+image_dir

wdef,0,900,900
!p.multi=[0,3,3]
loadct,0


ndir=n_elements(data_dir.dir397)

for idir=0,ndir-1 do begin
;idir = 1
	angle_q = atan(observation[idir].hxy[1], observation[idir].hxy[0])
	for iarm=0,2 do begin
;iarm = 0
		case iarm of
			0:begin
				title='Ca II H'
				dir=data_dir.dir397[idir]
				con=cal01_param.con397
			end
			1:begin
				title='Fe I 630 nm'
				dir=data_dir.dir630[idir]
				con=cal01_param.con630
			end
			2:begin
				title='Ca II 854 nm'
				dir=data_dir.dir854[idir]
				con=cal01_param.con854
			end
		endcase
		outfile0 = image_dir + string(iarm,format='(i1.1)') + '_'+string(idir,format='(i1.1)')+'.png'
		
		outdir=dir+outname_tail ; '_cal01'
		if is_dir(outdir) eq 0 then spawn,'mkdir '+outdir
		if iarm eq 0 then begin
			outdir_err=dir+outname_tail+'_err' ; '_cal01'
			if is_dir(outdir_err) eq 0 then spawn,'mkdir '+outdir_err
		endif

		filei=file_search(dir,'*_I_*.fits')
		fileq=file_search(dir,'*_Q_*.fits')
		fileu=file_search(dir,'*_U_*.fits')
		filev=file_search(dir,'*_V_*.fits')

		for ifile=0,n_elements(filei)-1 do begin
		;for ifile=cal01_param.plot_ifile,cal01_param.plot_ifile do begin
;ifile = 230
;ifile = 132
;ifile = cal01_param.plot_ifile
;ifile = 111
			read_sdo,filei[ifile],indexi,datai0, /use_shared_lib
			read_sdo,fileq[ifile],indexq,dataq0, /use_shared_lib
			read_sdo,fileu[ifile],indexu,datau0, /use_shared_lib
			read_sdo,filev[ifile],indexv,datav0, /use_shared_lib
			tmp=size(datai0) & nx=tmp[1] & ny=tmp[2]
		
			if idir eq 0 and iarm eq 0 and ifile eq 0 then begin
;if idir eq 1 and iarm eq 0 and ifile eq 111 then begin
				coe0_arr = fltarr(2600,n_elements(filei),ndir,3,3)  ; slit, scan, data, quv, lines
				coe1_arr = fltarr(2600,n_elements(filei),ndir,3,3)  ; slit, scan, data, quv, lines
			endif


			data1 = datai0>0.<0.
			for ipol=1,3 do begin
;ipol=2
				case ipol of
					1:begin
						data0=dataq0
						str_stks='Q'
					end
					2:begin
						data0=datau0
						str_stks='U'
					end
					3:begin
						data0=datav0
						str_stks='V'
					end
				endcase
				if iarm eq 0 and (ipol eq 1 or ipol eq 2) then begin
					lcen = [44,171,242,310,653,715,802,923]
					dd = 20
					nlcen = n_elements(lcen)
					coes = fltarr(2,nlcen)
					sigs = fltarr(2,nlcen)
					;if 0 then begin;--------------------------SKIP
					for l = 0,nlcen-1 do begin
						err = datai0[*,lcen[l]-dd:lcen[l]+dd]>1e-3<1e-3
						coes[*,l]=poly_fit(datai0[*,lcen[l]-dd:lcen[l]+dd], data0[*,lcen[l]-dd:lcen[l]+dd],1, measure_errors = err, sigma=sigma)
						sigs[*,l]=sigma
					endfor

					;coe = poly_fit(lcen,coes[0,*],1, measure_errors = sigs[0,*], sigma=sigma)
					coe = poly_fit(lcen,coes[0,*],1, sigma=sigma)
					coe0= coe[0]+coe[1]*findgen(ny)
					dcoe0=sqrt(sigma[0]^2 + (sigma[1]*findgen(ny))^2)

					coe = poly_fit(lcen,coes[1,*],1, sigma=sigma)
					coe1= coe[0]+coe[1]*findgen(ny)
					dcoe1=sqrt(sigma[0]^2 + (sigma[1]*findgen(ny))^2)

					if plot_cah_special then begin;----------------------------			
					if ifile eq cal01_param.plot_ifile then begin
						outfile = image_dir + 'cah_special_'+string(idir,format='(i1.1)')+'_'+  $
							string(ifile,format='(i3.3)')+'_'+string(ipol,format='(i2.2)')+'.png'
						wdef,1,800,800
						!p.multi=[0,2,2]
						chs = 1.5
						loadct,0
						plot_image,data0,/nosquare,title='Stokes '+str_stks, 	$
							xtitle='X (pixel)',ytitle='Y (pixel)',		$
							color=0,background=255,charsize=chs,ticklen=-0.02
						set_line_color
						for l = 0,nlcen-1 do begin
							oplot,!x.crange,replicate(lcen[l]-dd,2),color=3,line=2
							oplot,!x.crange,replicate(lcen[l]+dd,2),color=3,line=2
						endfor

						l=0
						plot,datai0[*,lcen[l]-dd:lcen[l]+dd], data0[*,lcen[l]-dd:lcen[l]+dd], psym=1,  $
							color=0,background=255,charsize=chs,xtitle='Stokes I',ytitle='Stokes '+str_stks
						oplot,[0,2],coes[0,l]+coes[1,l]*[0.,2.],color=3

						plot,lcen,coes[0,*],psym=1,  $
							color=0,background=1,charsize=chs,xtitle='Y (pixel)',ytitle='Coe 0'
						oplot,coe0,color=3
						oplot,coe0+dcoe0,color=3,line=1
						oplot,coe0-dcoe0,color=3,line=1

						plot,lcen,coes[1,*],psym=1,  $
							color=0,background=1,charsize=chs,xtitle='Y (pixel)',ytitle='Coe 1'
						oplot,coe1,color=3
						oplot,coe1+dcoe1,color=3,line=1
						oplot,coe1-dcoe1,color=3,line=1

						write_png,outfile,tvrd(/true)
					endif
					endif;------------------------------------------------
					data00 = data0					
					edata0 = data0
					for iy=0,ny-1 do begin
						data0[*,iy] = data0[*,iy] - (datai0[*,iy]*coe1[iy]+coe0[iy]) 
						edata0[*,iy]= sqrt( (datai0[*,iy]*dcoe1[iy])^2 + dcoe0[iy]^2 ) 
					endfor

					poly_order = 3
				endif else poly_order = 1
				for ix=0,nx-1 do begin
					prof = data0[ix,*]/datai0[ix,*]	
					coe = poly_fit(con,prof[con],poly_order)
					yfit= poly(findgen(ny),coe)
					data1[ix,*] = (prof - yfit)*datai0[ix,*]

					coe0_arr[ix,ifile,idir,ipol-1,iarm] = coe[0]
					coe1_arr[ix,ifile,idir,ipol-1,iarm] = coe[1]

					if ifile eq cal01_param.plot_ifile and ix eq cal01_param.plot_ix then begin
						wset,0
						set_line_color
						x=findgen(ny)
						plot,x,prof,/xstyle,yr=[-0.1,0.1],/ystyle,  $
						  xtitle='wavelength (pix)',ytitle=ytitle,title=title,  $
      						  color=0,background=1,charsize=2.5
  						oplot,con,prof[con],psym=1,color=3
      						oplot,x,yfit,color=5,thick=2
					endif					
				endfor

				case ipol of 
					1: begin
						dataq1=data1
						if iarm eq 0 then edataq1=edata0
					end
					2: begin
						datau1=data1
						if iarm eq 0 then edatau1=edata0
					end
					3: datav1=data1
				endcase
				if ifile eq cal01_param.plot_ifile then begin
					wset,0
					loadct,0
					plot_image,data0,/nosquare,  $
						color=0,background=255,charsize=2.5
					plot_image,data1,/nosquare,  $
						color=0,background=255,charsize=2.5
					if ipol eq 3 then write_png,outfile0,tvrd(/true)
				endif
			endfor

			; rotate +Q so that the +Q directs to the disk center
			dataq2 = cos(2.*angle_q)*dataq1 + sin(2.*angle_q)*datau1
			datau2 = -sin(2.*angle_q)*dataq1 + cos(2.*angle_q)*datau1

			if iarm eq 0 then begin
	                        edataq2 = sqrt( (cos(2.*angle_q)*edataq1)^2 + (sin(2.*angle_q)*edatau1)^2 )
	                        edatau2 = sqrt( (sin(2.*angle_q)*edataq1)^2 + (cos(2.*angle_q)*edatau1)^2 )

				outfile=outdir_err+strmid(fileq[ifile],strlen(dir),strlen(fileq[ifile])-strlen(dir))
				header=ta_make_hdr(indexq, /remove_z)
				writefits,outfile,edataq2,header

				outfile=outdir_err+strmid(fileu[ifile],strlen(dir),strlen(fileu[ifile])-strlen(dir))
				header=ta_make_hdr(indexu, /remove_z)
				writefits,outfile,edatau2,header
			endif
			

			; save data
			outfile=outdir+strmid(filei[ifile],strlen(dir),strlen(filei[ifile])-strlen(dir))
			header=ta_make_hdr(indexi, /remove_z)
			writefits,outfile,datai0,header

			outfile=outdir+strmid(fileq[ifile],strlen(dir),strlen(fileq[ifile])-strlen(dir))
			header=ta_make_hdr(indexq, /remove_z)
			writefits,outfile,dataq2,header

			outfile=outdir+strmid(fileu[ifile],strlen(dir),strlen(fileu[ifile])-strlen(dir))
			header=ta_make_hdr(indexu, /remove_z)
			writefits,outfile,datau2,header

			outfile=outdir+strmid(filev[ifile],strlen(dir),strlen(filev[ifile])-strlen(dir))
			header=ta_make_hdr(indexv, /remove_z)
			writefits,outfile,datav1,header
		endfor;ifile		
	endfor
endfor

!p.multi=0
loadct,0

save,coe0_arr,coe1_arr,file=save_file

END
;-----------------------------------------------
;-----------------------------------------------
PRO WFA_01, save_file=save_file
COMMON ANALYSIS

ix_hair_397 = [122.,1943.]
ix_hair_630 = [344.,1843.]
ix_hair_854 = [111.,2406.]
ibomb=0
ix   = bomb_pixel[ibomb].ix
;ix = [122,1943]
ix397= ix
ix630= [ (ix_hair_630[0] - ix_hair_630[1]) / (ix_hair_397[0] - ix_hair_397[1]) * (ix397[0] - ix_hair_397[1]) + ix_hair_630[1], $ 
	 (ix_hair_630[0] - ix_hair_630[1]) / (ix_hair_397[0] - ix_hair_397[1]) * (ix397[1] - ix_hair_397[1]) + ix_hair_630[1]  ]
ix854= (ix_hair_854[0] - ix_hair_854[1]) / (ix_hair_397[0] - ix_hair_397[1]) * (ix - ix_hair_397[1]) + ix_hair_854[1] 




cc = 2.99792458d+5 ; km/s 
range_630_vlos = [-20,20] 
range_630_width= 20.
range_854_width= 40.


ndir = n_elements(data_dir.dir397)

; for make array
nz=0
nx=0
for idir=0,ndir-1 do begin
	dir=data_dir.dir630[idir]+'_cal01/'
	filei=file_search(dir,'*_I_*.fits',count=nf)
	nz=nz>nf
	mreadfits,filei[0],index,datai,/nodata
	nx=nx>index[0].naxis1
endfor
nf630=nz
nx630=nx

nz=0
nx=0
for idir=0,ndir-1 do begin
        dir=data_dir.dir630[idir]+'_cal01/'
        filei=file_search(dir,'*_I_*.fits',count=nf)
        nz=nz>nf
        mreadfits,filei[0],index,datai,/nodata
        nx=nx>index[0].naxis1
endfor
nf854=nz
nx854=nx




for idir=0,ndir-1 do begin
;goto,jump
	; Fe I 630
	lambda0 = 6302.503 ;A
	geff    = 2.5
	gcap    = 6.25
	wvl=wavelength.wvl630
	dir=data_dir.dir630[idir]+'_cal01/'

	filei=file_search(dir,'*_I_*.fits',count=nf)
        fileq=file_search(dir,'*_Q_*.fits',count=nf)
        fileu=file_search(dir,'*_U_*.fits',count=nf)
        filev=file_search(dir,'*_V_*.fits',count=nf)
	c = 4.66e-13 * geff * lambda0^2 
	ct = (4.6686e-13 * lambda0^2)^2 * gcap 
	dlam = wvl[1] - wvl[0] 

	;ifile0=230
	;ifile0=135
	ifile0=0
	for ifile=ifile0,nf-1 do begin
		mreadfits,filei[ifile],index,datai
        	mreadfits,fileq[ifile],index,dataq
        	mreadfits,fileu[ifile],index,datau
        	mreadfits,filev[ifile],index,datav
		ll = sqrt(dataq^2 + datau^2)
		nx = index.naxis1
		ny = index.naxis2

		if ifile eq ifile0 and idir eq 0 then begin
			bv630={	lambda0: lambda0,		$
				geff: geff,			$
				gcap: gcap,			$
				blos: fltarr(nf630,nx630,ndir),	$
				bpos: fltarr(nf630,nx630,ndir),  $
				cph : fltarr(nf630,nx630,ndir),  $
				vlos: fltarr(nf630,nx630,ndir),  $
				width:fltarr(nf630,nx630,ndir)	}
		endif
		d1i = datai>0.<0.
		d2i = datai>0.<0.
		d1i[*,1:ny-2] = (datai[*,2:*] - datai[*,0:ny-3])/2./dlam
		d2i[*,1:ny-2] = (d1i[*,2:*] - d1i[*,0:ny-3])/2./dlam

		x = findgen(ny)
		wvl_center_pix = interpol(x,wvl[0:ny-1],lambda0)
		for ix=0,nx-1 do begin
;ix=1000
			;VLOS
			y = reform(datai[ix,*])
			estimates=[-0.5,wvl_center_pix,5.,1.]
			yfit = gaussfit(x,y,a,estimates=estimates,nterm=4)
			bv630.width[ifile,ix,idir] = a[2]
			if a[2] ge range_630_width then begin
;wdef,0,1200,500
;loadct,0
;!p.multi=0
;plot,x,y,psym=1	
;oplot,x,yfit
;stop
				a = [a[0], wvl_center_pix,5.,a[3]]
			endif
			wvl_center_pix_fit = interpol(wvl[0:ny-1],x,a[1])  ; A
			vlos = (wvl_center_pix_fit - lambda0)/lambda0*cc
			bv630.vlos[ifile,ix,idir] = vlos

			;BLOS
			wv_range_pix = a[1]+a[2]*[-3.,3.]
			wv_range_pix = wv_range_pix>0<(ny-3)
			bv630.blos[ifile,ix,idir] = - total( d1i[ix,wv_range_pix[0]:wv_range_pix[1]]*datav[ix,wv_range_pix[0]:wv_range_pix[1]] )/c/   $
							total( d1i[ix,wv_range_pix[0]:wv_range_pix[1]]^2)
;wdef,0,500,500
;loadct,0
;!p.multi=0
;plot,datav[ix,wv_range_pix[0]:wv_range_pix[1]],psym=1
;oplot,d1i[ix,wv_range_pix[0]:wv_range_pix[1]]*(-c)*bv630.blos[ix,ifile,idir]

			;BPOS
			wv_range_pix = a[1]+a[2]*[-4.,-1.5]
			wv_range_pix = wv_range_pix>0<(ny-3)
			bv630.bpos[ifile,ix,idir] = sqrt(  $
			       total( 4./3. * ll[ix,wv_range_pix[0]:wv_range_pix[1]] / ct 		$
			       		/ abs(wvl[wv_range_pix[0]:wv_range_pix[1]] - lambda0) * 	$
					abs(d1i[ix,wv_range_pix[0]:wv_range_pix[1]]) 			$
					) / 		$
			       total( abs((1./(wvl[wv_range_pix[0]:wv_range_pix[1]] - lambda0)))^2 * 	$
			       	      abs(d1i[ix,wv_range_pix[0]:wv_range_pix[1]])^2			$
				      ) 		$
				      )

;wdef,1,500,500
;loadct,0
;!p.multi=0
;plot,ll[ix,wv_range_pix[0]:wv_range_pix[1]],psym=1
;oplot,3./4.*bv630.bpos[ix,ifile,idir]^2 *ct /abs(wvl[wv_range_pix[0]:wv_range_pix[1]] - lambda0)*abs(d1i[ix,wv_range_pix[0]:wv_range_pix[1]])
;stop

			;Capital letter phi, Azimuth wrt LOS
			wv_range_pix = a[1]+a[2]*[-4.,-1.5]
			wv_range_pix = wv_range_pix>0<(ny-3)
                        bv630.cph[ifile,ix,idir] = 0.5 * atan(total(datau[ix,wv_range_pix[0]:wv_range_pix[1]]), total(dataq[ix,wv_range_pix[0]:wv_range_pix[1]])) 

		endfor
	endfor

jump:
	; Ca II 854
        lambda0 = 8542. ;A
        geff    = 1.1
        gcap    = 1.21
        wvl=wavelength.wvl854
        dir=data_dir.dir854[idir]+'_cal01/'

        filei=file_search(dir,'*_I_*.fits',count=nf)
        fileq=file_search(dir,'*_Q_*.fits',count=nf)
        fileu=file_search(dir,'*_U_*.fits',count=nf)
        filev=file_search(dir,'*_V_*.fits',count=nf)
        c = 4.66e-13 * geff * lambda0^2
        ct = (4.6686e-13 * lambda0^2)^2 * gcap
        dlam = wvl[1] - wvl[0]

;        ifile0=230
;        ifile0=135
        ifile0=0
;ifile0=120
        for ifile=ifile0,nf-1 do begin
                mreadfits,filei[ifile],index,datai
                mreadfits,fileq[ifile],index,dataq
                mreadfits,fileu[ifile],index,datau
                mreadfits,filev[ifile],index,datav
                ll = sqrt(dataq^2 + datau^2)
                nx = index.naxis1
                ny = index.naxis2

                if ifile eq ifile0 and idir eq 0 then begin
                        bv854={ lambda0: lambda0,               $
                                geff: geff,                     $
                                gcap: gcap,                     $
                                blos: fltarr(nf854,nx854,ndir),       $
                                bpos: fltarr(nf854,nx854,ndir),  $
                                cph : fltarr(nf854,nx854,ndir),  $
				vlos: fltarr(nf854,nx854,ndir),  $
				width:fltarr(nf854,nx854,ndir)	}
                endif
                d1i = datai>0.<0.
                d2i = datai>0.<0.
                d1i[*,1:ny-2] = (datai[*,2:*] - datai[*,0:ny-3])/2./dlam
                d2i[*,1:ny-2] = (d1i[*,2:*] - d1i[*,0:ny-3])/2./dlam

                x = findgen(ny)
                wvl_center_pix = interpol(x,wvl[0:ny-1],lambda0)
                for ix=0,nx-1 do begin
;ix=830
;ix=200
                        ;VLOS
                        y = reform(datai[ix,*])
                        estimates=[-0.9,wvl_center_pix,5.,1.]
			pos=where(wvl ge (lambda0-.5) and wvl le (lambda0+.5),npos)
                        yfit = gaussfit(x[pos],y[pos],a,estimates=estimates,nterm=4)
			bv854.width[ifile,ix,idir] = a[2]
range_854_width = 80 

			if a[2] ge range_854_width then begin
;wdef,0,1200,500
;loadct,0
;!p.multi=0
;set_line_color
;plot,x,y,psym=1	
;oplot,x[pos],yfit,color=3
;loadct,0
;stop
				a = [a[0], wvl_center_pix,10.,a[3]]
			endif
                        wvl_center_pix_fit = interpol(wvl[0:ny-1],x,a[1])  ; A
                        bv854.vlos[ifile,ix,idir] = (wvl_center_pix_fit - lambda0)/lambda0*cc

                        ;BLOS
                        wv_range_pix = a[1]+a[2]*[-2.,2.]
			wv_range_pix = wv_range_pix>0<(ny-3)
                        bv854.blos[ifile,ix,idir] = - total( d1i[ix,wv_range_pix[0]:wv_range_pix[1]]*datav[ix,wv_range_pix[0]:wv_range_pix[1]] )/c/   $
                                                        total( d1i[ix,wv_range_pix[0]:wv_range_pix[1]]^2)


;wdef,0,500,500
;loadct,0
;!p.multi=0
;tmp1=datav[ix,wv_range_pix[0]:wv_range_pix[1]]
;tmp2=d1i[ix,wv_range_pix[0]:wv_range_pix[1]]*(-c)*bv854.blos[ifile,ix,idir]
;plot,tmp1,psym=1,yr=minmax([tmp1,tmp2])
;oplot,tmp2
;stop

                        ;BPOS
                        wv_range_pix = a[1]+a[2]*[-2.,-0.8]
			wv_range_pix = wv_range_pix>0<(ny-3)
                        bv854.bpos[ifile,ix,idir] = sqrt(  $
                               total( 4./3. * ll[ix,wv_range_pix[0]:wv_range_pix[1]] / ct               $
                                        / abs(wvl[wv_range_pix[0]:wv_range_pix[1]] - lambda0) *         $
                                        abs(d1i[ix,wv_range_pix[0]:wv_range_pix[1]])                    $
                                        ) /             $
                               total( abs((1./(wvl[wv_range_pix[0]:wv_range_pix[1]] - lambda0)))^2 *    $
                                      abs(d1i[ix,wv_range_pix[0]:wv_range_pix[1]])^2                    $
                                      )                 $
                                      )

;wdef,1,500,500
;loadct,0
;!p.multi=0
;plot,ll[ix,wv_range_pix[0]:wv_range_pix[1]],psym=1
;oplot,3./4.*bv854.bpos[ix,ifile,idir]^2 *ct /abs(wvl[wv_range_pix[0]:wv_range_pix[1]] - lambda0)*abs(d1i[ix,wv_range_pix[0]:wv_range_pix[1]])

                        ;Capital letter phi, Azimuth wrt LOS
                        wv_range_pix = a[1]+a[2]*[-2.,-0.8]
			wv_range_pix = wv_range_pix>0<(ny-3)
                        bv854.cph[ifile,ix,idir] = 0.5 * atan(total(datau[ix,wv_range_pix[0]:wv_range_pix[1]]), total(dataq[ix,wv_range_pix[0]:wv_range_pix[1]]))

                endfor
;plot,bv854.blos[*,ifile,idir],yr=[-1000,1000] & oplot,replicate(ix854[0],2),!y.crange & oplot,replicate(ix854[1],2),!y.crange & print,mean(bv854.blos[ix854[0]:ix854[1],ifile,idir])
;stop
;plot,bv854.bpos[*,ifile,idir],yr=[0,1000] & oplot,replicate(ix854[0],2),!y.crange & oplot,replicate(ix854[1],2),!y.crange & print,mean(bv854.bpos[ix854[0]:ix854[1],ifile,idir])
;stop
;plot,bv854.cph[*,ifile,idir]*!radeg,psym=1 & oplot,replicate(ix854[0],2),!y.crange & oplot,replicate(ix854[1],2),!y.crange & print,mean(bv854.cph[ix854[0]:ix854[1],ifile,idir])*!radeg
;stop		
        endfor
endfor

save,bv630,bv854,file=save_file

END
;-----------------------------------------------
FUNCTION WFA,stki,stkq,stku,stkv
COMMON ANALYSIS
; ref) Kuridze et al. 2019, ApJ

pos397 = where(wvl397 ge nicole_param.wvl397_range[0] and wvl397 le nicole_param.wvl397_range[1], npos397)
pos630 = where(wvl630 ge nicole_param.wvl630_range[0] and wvl630 le nicole_param.wvl630_range[1], npos630)
pos854 = where(wvl854 ge nicole_param.wvl854_range[0] and wvl854 le nicole_param.wvl854_range[1], npos854)
wave397= wvl397[pos397]
wave630= wvl630[pos630]
wave854= wvl854[pos854]


stki397 = stki[*,*,0:npos397-1]
stkq397 = stkq[*,*,0:npos397-1]
stku397 = stku[*,*,0:npos397-1]
stkv397 = stkv[*,*,0:npos397-1]
stki630 = stki[*,*,npos397:npos397+npos630-1]
stkq630 = stkq[*,*,npos397:npos397+npos630-1]
stku630 = stku[*,*,npos397:npos397+npos630-1]
stkv630 = stkv[*,*,npos397:npos397+npos630-1]
stki854 = stki[*,*,npos397+npos630:*]
stkq854 = stkq[*,*,npos397+npos630:*]
stku854 = stku[*,*,npos397+npos630:*]
stkv854 = stkv[*,*,npos397+npos630:*]



tmp=size(stki) & nx=tmp[1] & ny=tmp[2]
res = fltarr(nx,ny,3,2)

for iline=0,1 do begin
	case iline of
		0:begin
			; Fe 630
			lambda0 = 6301. ;A
			geff    = 1.6667
			gcap    = 2.5167

			wvl = wave630
			wvlr= nicole_param.wfa630_range
			si  = stki630
			sq  = stkq630
			su  = stku630
			sv  = stkv630
		end
		2:begin
			print,'Not yet detrmine nicole_param.wfa397_range'
			; Ca 397
			lambda0 = 3968. ;A
			geff    = 1.3333
			gcap    = 1.3333

			wvl = wave397
			wvlr= nicole_param.wfa397_range
			si  = stki397
			sq  = stkq397
			su  = stku397
			sv  = stkv397
		end
		1:begin
			; Ca 854
			lambda0 = 8542. ;A
			geff    = 1.1
			gcap    = 1.21

			wvl = wave854
			wvlr= nicole_param.wfa854_range
			si  = stki854
			sq  = stkq854
			su  = stku854
			sv  = stkv854
		end
	endcase

	factor_v_di = -4.67e-13 * geff * lambda0^2 ; v/(di/dl) 
	factor_l_di = -5.45e-26 * gcap * lambda0^4 ; v/(di/dl) 
	dlam = wvl[1] - wvl[0] 
	npos = n_elements(wvl)
	pos  = where(wvl ge wvlr[0] and wvl le wvlr[1])
	for ix=0,nx-1 do begin
		for iy=0,ny-1 do begin

			d1prof = fltarr(npos)
			d1prof[1:npos-2] = (si[ix,iy,2:*] - si[ix,iy,0:npos-3]) /2./dlam

			d2prof = fltarr(npos)
			d2prof[2:npos-3] = (d1prof[3:*] - d1prof[1:npos-4]) /2./dlam

			blos0 = reform(sv[ix,iy,1:npos-2]) / factor_v_di / d1prof[1:npos-2]
			qu = reform(sqrt(sq[ix,iy,*]^2 + su[ix,iy,*]^2))
			bpos0 = sqrt(abs(qu[2:npos-3]/factor_l_di/d2prof[2:npos-3]))

			blos  = median(blos0[pos]) 
			bpos  = median(bpos0[pos]) 

			tmp = max([max(abs(sq[ix,iy,*])),max(abs(su[ix,iy,*]))],key)
			if key eq 0 then tmp=max(abs(sq[ix,iy,*]),pos1) else tmp=max(abs(su[ix,iy,*]),pos1)
			chi = atan(su[ix,iy,pos1],sq[ix,iy,pos1])/2.
		
			res[ix,iy,*,iline] = [bpos*cos(chi), bpos*sin(chi), blos]  ; Bx, By, Bz		
		endfor
	endfor
endfor

RETURN, res
END
;-----------------------------------------------
;-----------------------------------------------
PRO MAKE_DATA_01, save_dir = save_dir, atm_disp_file = atm_disp_file, wfa_los = wfa_los, h1e_ncp_01 = h1e_ncp_01, ref = ref, wfa_solarnormal = wfa_solarnormal
COMMON ANALYSIS


if is_dir(save_dir) eq 0 then spawn,'mkdir '+save_dir

restore,h1e_ncp_01  & h1e_ncp0 = h1e_ncp
;save,h1e_ncp,file=save_file

restore,wfa_solarnormal
;save,b630,thb630,phb630,b854,thb854,phb854,file=save_file

; Wavelength
mreadfits,(file_search(data_dir.dir630[0]+'_cal01/','*_I_*.fits'))[0],tmp,da,/nodata
npos630 = tmp.naxis2
pos630 = where(wavelength.wvl630[0:npos630-1] ge nicole_param.wvl630_range[0] and wavelength.wvl630[0:npos630-1] le nicole_param.wvl630_range[1], npos630)
wave_630=wavelength.wvl630[pos630]

mreadfits,(file_search(data_dir.dir397[0]+'_cal01/','*_I_*.fits'))[0],tmp,da,/nodata
npos397 = tmp.naxis2
pos397 = where(wavelength.wvl397[0:npos397-1] ge nicole_param.wvl397_range[0] and wavelength.wvl397[0:npos397-1] le nicole_param.wvl397_range[1], npos397)
wave_397=wavelength.wvl397[pos397]

mreadfits,(file_search(data_dir.dir854[0]+'_cal01/','*_I_*.fits'))[0],tmp,da,/nodata
npos854 = tmp.naxis2
pos854 = where(wavelength.wvl854[0:npos854-1] ge nicole_param.wvl854_range[0] and wavelength.wvl854[0:npos854-1] le nicole_param.wvl854_range[1], npos854)
wave_854=wavelength.wvl854[pos854]


; Stokes profiles 10"x~10", 0.2"x~0.2"
; & wfa results
dd    = 10. ; arcsec
;dd    = 1. ; arcsec for test
restore,atm_disp_file,/verb & atm_disp = res

nbomb = n_elements(bomb_pixel)
for ibomb=0,nbomb-1 do begin
	save_file = save_dir+string(ibomb,format='(i2.2)')+'.sav'

	idata=bomb_pixel[ibomb].data

	filei_630=file_search(data_dir.dir630[idata]+'_cal01/','*_I_*.fits')
	filei_397=file_search(data_dir.dir397[idata]+'_cal01/','*_I_*.fits')
	filei_854=file_search(data_dir.dir854[idata]+'_cal01/','*_I_*.fits')

	fileq_630=file_search(data_dir.dir630[idata]+'_cal01/','*_Q_*.fits')
	fileq_397=file_search(data_dir.dir397[idata]+'_cal01/','*_Q_*.fits')
	fileq_854=file_search(data_dir.dir854[idata]+'_cal01/','*_Q_*.fits')

	fileu_630=file_search(data_dir.dir630[idata]+'_cal01/','*_U_*.fits')
	fileu_397=file_search(data_dir.dir397[idata]+'_cal01/','*_U_*.fits')
	fileu_854=file_search(data_dir.dir854[idata]+'_cal01/','*_U_*.fits')

	filev_630=file_search(data_dir.dir630[idata]+'_cal01/','*_V_*.fits')
	filev_397=file_search(data_dir.dir397[idata]+'_cal01/','*_V_*.fits')
	filev_854=file_search(data_dir.dir854[idata]+'_cal01/','*_V_*.fits')

	if not keyword_set(ref) then begin
		ycen_pix = mean(bomb_pixel[ibomb].ix)
	endif else begin
		ycen_pix = mean(bomb_pixel[ibomb].ix_ref1)
	endelse
	xcen_pix = bomb_pixel[ibomb].ifile

	xycen_arc= pix2arcsec_nosun(ycen_pix, xcen_pix, index_for_map.index_nosun_397)
	xr_arc   = xycen_arc[0] + dd*0.5*[-1.,1.]
	yr_arc   = xycen_arc[1] + dd*0.5*[-1.,1.]

	xr_pix630= arcsec2pix_nosun(xr_arc, yr_arc, index_for_map.index_nosun_630, atm_disp=atm_disp[1,idata])

	xs    = round(xr_pix630[0,0])
	nx    = round(xr_pix630[1,0]) - xs + 1
	nx    = (fix(nx/binx.bin630) + 1)*binx.bin630
	sfile = round(xr_pix630[0,1])
	nfile = round(xr_pix630[1,1]) - sfile + 1

	tmp = fltarr(nfile,nx,npos630) & stki0_630 = tmp & stkq0_630 = tmp & stku0_630 = tmp & stkv0_630 = tmp
	tmp = fltarr(nfile,nx,npos397) & stki0_397 = tmp & stkq0_397 = tmp & stku0_397 = tmp & stkv0_397 = tmp
	tmp = fltarr(nfile,nx,npos854) & stki0_854 = tmp & stkq0_854 = tmp & stku0_854 = tmp & stkv0_854 = tmp
	h1e_ncp1 = fltarr(nfile,nx)
	fe_vlos1  = fltarr(nfile,nx)
	restore,wfa_los,/verb
	b_los_xyz_0 = fltarr(nfile,nx,3,2)
	bv    = bv630
	iii   = 0
	b_los_xyz_0[*,*,0,iii] = bv.blos[sfile:sfile+nfile-1, xs:xs+nx-1 , idata]
	b_los_xyz_0[*,*,1,iii] = (bv.bpos*cos(bv.cph))[sfile:sfile+nfile-1, xs:xs+nx-1, idata]
	b_los_xyz_0[*,*,2,iii] = (bv.bpos*sin(bv.cph))[sfile:sfile+nfile-1, xs:xs+nx-1, idata]
	fe_vlos1[*,*] = bv.vlos[sfile:sfile+nfile-1, xs:xs+nx-1 , idata] - median(bv.vlos[*, *, idata])
;print,sfile,sfile+nfile-1, xs,xs+nx-1 , idata
;print,median(bv.vlos[*, *, idata])
;stop
	wfa630_los_blos1    = bv.blos[sfile:sfile+nfile-1, xs:xs+nx-1 , idata]
       	wfa630_los_bpos1    = bv.bpos[sfile:sfile+nfile-1, xs:xs+nx-1 , idata] 	
       	wfa630_los_bph1     = bv.cph[sfile:sfile+nfile-1, xs:xs+nx-1 , idata] 	


	bv    = bv854
	iii   = 1

	wfa630_local_b1    = b630[sfile:sfile+nfile-1, xs:xs+nx-1 , idata]
	wfa630_local_thb1  = thb630[sfile:sfile+nfile-1, xs:xs+nx-1 , 1, idata]*!radeg
	wfa630_local_phb1  = phb630[sfile:sfile+nfile-1, xs:xs+nx-1 , 1, idata]*!radeg

	if ibomb eq 7 then begin ; Manuaaly handle 180 degree ambiguity
		wfa630_local_thb0 = wfa630_local_thb1
		tmp = thb630[sfile:sfile+nfile-1, xs:xs+nx-1 , 2, idata]*!radeg
		pos = where(wfa630_local_thb1 ge 50)
		wfa630_local_thb1[pos] = tmp[pos]
		
		tmp = phb630[sfile:sfile+nfile-1, xs:xs+nx-1 , 2, idata]*!radeg
		wfa630_local_phb1[pos] = tmp[pos]
	endif

	tmp = wfa630_local_phb1
	pos = where(tmp lt 0)
	wfa630_local_phb1[pos] = tmp[pos] + 360.

	for ifile=sfile,sfile+nfile-1 do begin
	        mreadfits,filei_630[ifile],index,data  &  stki0_630[ifile-sfile,*,*] = data[xs:xs+nx-1,pos630]
	        mreadfits,fileq_630[ifile],index,data  &  stkq0_630[ifile-sfile,*,*] = data[xs:xs+nx-1,pos630]
	        mreadfits,fileu_630[ifile],index,data  &  stku0_630[ifile-sfile,*,*] = data[xs:xs+nx-1,pos630]
	        mreadfits,filev_630[ifile],index,data  &  stkv0_630[ifile-sfile,*,*] = data[xs:xs+nx-1,pos630]

		mreadfits,filei_397[ifile],index,datai397
		mreadfits,fileq_397[ifile],index,dataq397
		mreadfits,fileu_397[ifile],index,datau397
		mreadfits,filev_397[ifile],index,datav397
		nx397 = index[0].naxis1

		mreadfits,filei_854[ifile],index,datai854
	        mreadfits,fileq_854[ifile],index,dataq854
	        mreadfits,fileu_854[ifile],index,datau854
	        mreadfits,filev_854[ifile],index,datav854
		nx854 = index[0].naxis1
		for iix=xs,xs+nx-1 do begin
	                xy_arc= pix2arcsec_nosun(iix, ifile, index_for_map.index_nosun_630, atm_disp=atm_disp[1,idata])
        	        iix1  = (arcsec2pix_nosun(xy_arc[0], xy_arc[1], index_for_map.index_nosun_397, atm_disp=atm_disp[0,idata]))[0]
			if iix1 ge 0 and iix1 le (nx397-1) then begin
				stki0_397[ifile-sfile,iix-xs,*] = datai397[iix1,pos397]
        		        stkq0_397[ifile-sfile,iix-xs,*] = dataq397[iix1,pos397]
        		        stku0_397[ifile-sfile,iix-xs,*] = datau397[iix1,pos397]
        	        	stkv0_397[ifile-sfile,iix-xs,*] = datav397[iix1,pos397]
				h1e_ncp1[ifile-sfile,iix-xs]     = h1e_ncp0[ifile,iix1,idata]
			endif
			;mreadfits,filev_397[ifile],index,data

        	        iix1  = (arcsec2pix_nosun(xy_arc[0], xy_arc[1], index_for_map.index_nosun_854, atm_disp=atm_disp[2,idata]))[0]
			if iix1 ge 0 and iix1 le (nx854-1) then begin
				stki0_854[ifile-sfile,iix-xs,*] = datai854[iix1,pos854]
                		stkq0_854[ifile-sfile,iix-xs,*] = dataq854[iix1,pos854]
                		stku0_854[ifile-sfile,iix-xs,*] = datau854[iix1,pos854]
                		stkv0_854[ifile-sfile,iix-xs,*] = datav854[iix1,pos854]
;mreadfits,filev_854[ifile],index,data
;stop
				b_los_xyz_0[ifile-sfile,iix-xs,0,iii] = bv.blos[ifile, iix1 , idata]
				b_los_xyz_0[ifile-sfile,iix-xs,1,iii] = (bv.bpos*cos(bv.cph))[ifile, iix1, idata]
				b_los_xyz_0[ifile-sfile,iix-xs,2,iii] = (bv.bpos*sin(bv.cph))[ifile, iix1, idata]
			endif
		endfor
	endfor
	stki_630 = rebin(stki0_630,nfile,nx/binx.bin630,npos630)
	stkq_630 = rebin(stkq0_630,nfile,nx/binx.bin630,npos630)
	stku_630 = rebin(stku0_630,nfile,nx/binx.bin630,npos630)
	stkv_630 = rebin(stkv0_630,nfile,nx/binx.bin630,npos630)
	fe_vlos  = rebin(fe_vlos1,nfile,nx/binx.bin630)
	wfa630_local_b   = rebin(wfa630_local_b1,nfile,nx/binx.bin630)
	wfa630_local_thb = rebin(wfa630_local_thb1,nfile,nx/binx.bin630)
	wfa630_local_phb = rebin(wfa630_local_phb1,nfile,nx/binx.bin630)

	wfa630_los_blos    = rebin(wfa630_los_blos1,nfile,nx/binx.bin630)
       	wfa630_los_bpos    = rebin(wfa630_los_bpos1,nfile,nx/binx.bin630)
       	wfa630_los_bph     = rebin(wfa630_los_bph1,nfile,nx/binx.bin630)

	stki_397 = rebin(stki0_397,nfile,nx/binx.bin630,npos397)
	stkq_397 = rebin(stkq0_397,nfile,nx/binx.bin630,npos397)
	stku_397 = rebin(stku0_397,nfile,nx/binx.bin630,npos397)
	stkv_397 = rebin(stkv0_397,nfile,nx/binx.bin630,npos397)
	h1e_ncp  = rebin(h1e_ncp1,nfile,nx/binx.bin630)

	stki_854 = rebin(stki0_854,nfile,nx/binx.bin630,npos854)
	stkq_854 = rebin(stkq0_854,nfile,nx/binx.bin630,npos854)
	stku_854 = rebin(stku0_854,nfile,nx/binx.bin630,npos854)
	stkv_854 = rebin(stkv0_854,nfile,nx/binx.bin630,npos854)

	b_los_xyz = rebin(b_los_xyz_0,nfile,nx/binx.bin630,3,2)

	; index for map
	index = index_for_map.index_nosun_630
	index.CDELT2 = index.CDELT2*binx.bin630
	index.xcen   = xycen_arc[0]
	index.ycen   = xycen_arc[1]
	index.naxis1 = nfile 
	index.naxis2 = nx/binx.bin630


	save, 	stki_630,stkq_630,stku_630,stkv_630,wave_630,	$
		stki_397,stkq_397,stku_397,stkv_397,wave_397,	$
		stki_854,stkq_854,stku_854,stkv_854,wave_854,	$
		b_los_xyz,h1e_ncp,      $
		ibomb,index,dd,	$
		file=save_file


	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_h1encp.fits'
	hdr = ta_make_hdr(index)
	writefits,fits_file,h1e_ncp,hdr

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_fe_vi.fits'
	hdr = ta_make_hdr(index)
	tmp = min(abs(wave_630 - wavelength.wvl_fev),pos)
	writefits,fits_file,reform((stkv_630/stki_630)[*,*,pos]),hdr
	tmp = reform((stkv_630/stki_630)[*,*,pos])

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_cah3.fits'
	hdr = ta_make_hdr(index)
	tmp = min(abs(wave_397 - wavelength.wvl_cah3),pos)
	writefits,fits_file,reform(stki_397[*,*,pos]),hdr

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_cah1v.fits'
	hdr = ta_make_hdr(index)
	tmp = min(abs(wave_397 - wavelength.wvl_cah1v),pos)
	writefits,fits_file,reform(stki_397[*,*,pos]),hdr

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_fe_vlos.fits'
	hdr = ta_make_hdr(index)
	writefits,fits_file,fe_vlos,hdr

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_wfa630_local_b.fits'
	hdr = ta_make_hdr(index)
	writefits,fits_file,wfa630_local_b,hdr

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_wfa630_local_thb.fits'
	hdr = ta_make_hdr(index)
	writefits,fits_file,wfa630_local_thb,hdr

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_wfa630_local_phb.fits'
	hdr = ta_make_hdr(index)
	writefits,fits_file,wfa630_local_phb,hdr

	; B in LOS frame
	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_wfa630_los_blos.fits'
        hdr = ta_make_hdr(index)
        writefits,fits_file,wfa630_los_blos,hdr

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_wfa630_los_bpos.fits'
        hdr = ta_make_hdr(index)
        writefits,fits_file,wfa630_los_bpos,hdr

	fits_file = save_dir+string(ibomb,format='(i2.2)')+'_wfa630_los_bph.fits'
        hdr = ta_make_hdr(index)
        writefits,fits_file,wfa630_los_bph,hdr

endfor

END
;-----------------------------------------------
;-------------------------------------------------
PRO stray1, x, a, f, pder

xx=findgen(n_elements(x))
f = (x - a[0]) * (a[1]*xx+a[2])

END
;-----------------------------------------------
;=========================================================================
PRO MAKE_H1E_NCP_01,save_file=save_file, image_file=image_file
COMMON ANALYSIS

width = 5
plot_idir = 0
plot_iz   = 135
plot_ix   = 800

key=1
ndir=n_elements(data_dir.dir397)
for idir=0,ndir-1 do begin
;idir=1
;idir = plot_idir[2]
	if 1 then begin
		file=file_search(data_dir.dir397[idir]+'_cal01','*_I_*')
		mreadfits,file,indexi,datai
		file=file_search(data_dir.dir397[idir]+'_cal01','*_V_*')
		mreadfits,file,indexv,datav
		save,indexi,datai,indexv,datav,file='/data/tanan/idl/tmp.sav'
	endif else restore,'/data/tanan/idl/tmp.sav',/verb

	tmp = size(datai) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]
	x = findgen(ny)
	if key then begin
		h1e_ncp=fltarr(nz,2558,ndir)
		key=0
	endif
	for iz=0,nz-1 do begin
		if iz mod 10 eq 0 then print,idir,iz
;	for iz=130,140 do begin
;iz=89
;iz=111  ; target
;iz=80
;iz=plot_iz[2]
		for ix=0,nx-1 do begin
;ix=874
;ix=1031 ; target
;ix=1046
;ix=plot_ix[2]
			profi = reform(datai[ix,*,iz])
			profv = reform(datav[ix,*,iz])
			pos = where((wavelength.wvl397 ge wavelength.wvl_h1e_r[0] and wavelength.wvl397 le wavelength.wvl_h1e_r[1]) or $
				(wavelength.wvl397 ge wavelength.wvl_h1e_r[2] and wavelength.wvl397 le wavelength.wvl_h1e_r[3]) )
			tmp = min(abs(wavelength.wvl397-wavelength.wvl_h1e),pos_center)
			a0 = [profi[pos_center]-max(profi[pos]), pos_center,5.,max(profi[pos]),0.]
			yfit = gaussfit(x[pos],profi[pos],a,nterms=5,estimates=a0)

			if wavelength.wvl397[(a[1]-width)>0<(ny-1)] ge wavelength.wvl_h1e_r[0] and $
			   wavelength.wvl397[(a[1]+width)>0<(ny-1)] le wavelength.wvl_h1e_r[1] then begin
				h1e_ncp[iz,ix,idir]=mean((reform(datav[ix,*,iz]/datai[ix,*,iz]))[a[1]-width:a[1]+width])
			endif else begin
				h1e_ncp[iz,ix,idir]=-100.   ;datav[ix,pos_center,iz]/datai[ix,pos_center,iz]
			endelse
;print,h1e_ncp[iz,ix,idir]

                       	if idir eq plot_idir and iz eq plot_iz and ix eq plot_ix then begin
			;if 1 then begin
				wdef,0,1200,800
				set_line_color
				!p.multi=[0,1,2]
				plot,wavelength.wvl397,profi,psym=1,color=0,background=1,  $
					xr=[3968,3971],  $
					charsize=3,/xstyle,xtitle='Wavelength (A)',ytitle='I'
				oplot,wavelength.wvl397[pos],yfit,color=3,thick=3,psym=1
				oplot,replicate(wavelength.wvl397[a[1]],2),!y.crange,line=0,color=3
				oplot,replicate(wavelength.wvl397[a[1]+width],2),!y.crange,line=2,color=3
				oplot,replicate(wavelength.wvl397[a[1]-width],2),!y.crange,line=2,color=3

				plot,wavelength.wvl397,profv/profi,yr=[-0.01,0.01],/ystyle,  $
					xr=[3968,3971],  $
					psym=-1,color=0,background=1,  $
					charsize=3,/xstyle,xtitle='Wavelength (A)',ytitle='V/I'
				oplot,replicate(wavelength.wvl397[a[1]],2),!y.crange,line=0,color=3
				oplot,replicate(wavelength.wvl397[a[1]+width],2),!y.crange,line=2,color=3
				oplot,replicate(wavelength.wvl397[a[1]-width],2),!y.crange,line=2,color=3
				oplot,!x.crange,[0.,0.],color=0
;stop
			endif
		endfor
	endfor
endfor


save,h1e_ncp,file=save_file

write_png,image_file,tvrd(/true)
!p.multi=0
loadct,0

END
;=========================================================================
PRO fourgaussian, x, a, f

z1 = (x - a[6])/a[7]
z2 = (x - a[8])/a[9]
z3 = (x - a[10])/a[11]
z4 = (x - a[12])/a[13]

f = a[0] + a[1]*x + a[2]*exp(-0.5*z1^2) + a[3]*exp(-0.5*z2^2) + a[4]*exp(-0.5*z3^2) + a[5]*exp(-0.5*z4^2)   

END

;=========================================================================
PRO MAKE_H1E_NCP_02,save_file=save_file, image_dir=image_dir
COMMON ANALYSIS

cc=2.99792458d+10

if is_dir(image_dir) eq 0 then spawn,'mkdir '+image_dir

width = 5
plot_idir = [0   ,0    ,0    ,0   ,0   ,0   ,0    ,0   ,0   ,0   ,1   ,1  ,1   ,1   ,1   ,1   ,1   ,1  ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1   ,1]
plot_iz   = [135 ,71   ,135  ,154 ,112 ,75  ,147  ,296 ,335 ,209 ,111 ,89 ,79  ,140 ,126 ,64  ,74  ,44 ,217 ,121 ,161 ,216 ,250 ,295 ,369 ,357 ,243]
plot_ix   = [800 ,1347 ,552  ,1171,721 ,1078,205  ,2195,2008,1522,1031,867,1143,1054,1428,1809,2252,37 ,1578,2364,2202,2208,2189,2189,1466,1004,180]
threshold = 2e-3

pix_wvl_h1e_r_0 = interpol(findgen(n_elements(wavelength.wvl397)),wavelength.wvl397,wavelength.wvl_h1e_r[0])
pix_wvl_h1e_r_1 = interpol(findgen(n_elements(wavelength.wvl397)),wavelength.wvl397,wavelength.wvl_h1e_r[1])

key=1
ndir=n_elements(data_dir.dir397)
for idir=0,ndir-1 do begin
;idir=1
;idir = plot_idir[2]
        if 1 then begin
                file=file_search(data_dir.dir397[idir]+'_cal01','*_I_*')
                mreadfits,file,indexi,datai
                file=file_search(data_dir.dir397[idir]+'_cal01','*_V_*')
                mreadfits,file,indexv,datav
                save,indexi,datai,indexv,datav,file='/data/tanan/idl/tmp.sav'
        endif else restore,'/data/tanan/idl/tmp.sav',/verb

        tmp = size(datai) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]
        x = findgen(ny)
        if key then begin
                h1e_ncp = fltarr(nz,2558,ndir)
                h1e_vlos= fltarr(nz,2558,ndir)
                key=0
        endif


	for iz=0,nz-1 do begin
;       for iz=130,140 do begin
;iz=89
;iz=111  ; target
;iz=80
;iz=plot_iz[2]
                for ix=0,nx-1 do begin
;ix=2
;ix=1031 ; target
;ix=1046 ; -0.1%
;ix=plot_ix[2]
                        profi = reform(datai[ix,*,iz])
                        profv = reform(datav[ix,*,iz])
			dprofi=profi>0.<0.
			dprofi[1:ny-2] = profi[2:ny-1] - profi[0:ny-3]

			nmargin=0.01
			pos = where(wavelength.wvl397 ge (wavelength.wvl_h1e_r[0]+0.05+nmargin) and $
				wavelength.wvl397 le (wavelength.wvl_h1e_r[1] - nmargin), npos)
			coe = poly_fit(x[pos],dprofi[pos],1, yfit=yfit1)
			rms = sqrt(mean((yfit1-dprofi[pos])^2))
			while (rms ge threshold) and (npos ge 10) do begin
				nmargin=nmargin+0.01
				pos = where(wavelength.wvl397 ge (wavelength.wvl_h1e_r[0]+0.05+nmargin) and $
					wavelength.wvl397 le (wavelength.wvl_h1e_r[1] - nmargin), npos)
				coe = poly_fit(x[pos],dprofi[pos],1, yfit=yfit1)
				rms = sqrt(mean((yfit1-dprofi[pos])^2))
				;print,nmargin, npos, rms
				;plot,x[pos],dprofi[pos],psym=1,charsize=2
				;oplot,x[pos],yfit1
				;stop
			endwhile
;plot,x[pmin:pmax],dprofi[pmin:pmax],psym=1
;oplot,x[pmin:pmax],yfit
			pcen = -coe[0]/coe[1]
			if pcen ge pix_wvl_h1e_r_0 and pcen le pix_wvl_h1e_r_1 then $
				h1e_ncp[iz,ix,idir]=mean((reform(datav[ix,*,iz]/datai[ix,*,iz]))[pcen-width:pcen+width])

				wcen = interpol(wavelength.wvl397,findgen(n_elements(wavelength.wvl397)),pcen)			
				h1e_vlos[iz,ix,idir]=(wcen - wavelength.wvl_h1e)/wavelength.wvl_h1e*cc
;print,h1e_ncp[iz,ix,idir]
			for iplot=0,n_elements(plot_idir)-1 do begin
				if idir eq plot_idir[iplot] and iz eq plot_iz[iplot] and ix eq plot_ix[iplot] then begin
;                        if 1 then begin

					outfile=image_dir+'idir'+string(idir,format='(i2.2)')+'_z'+string(iz,format='(i3.3)')+  $
						'_x'+string(ix,format='(i4.4)')+'.png'

					wdef,0,1200,800
        	                        set_line_color
                	                !p.multi=[0,1,3]
                        	        plot,wavelength.wvl397,profi,psym=1,color=0,background=1,  $
                                	        xr=[3968,3971],  $
						title='Data: '+string(idir,format='(i2.2)')+', scan: '+string(iz,format='(i3.3)')+  $
                                                ', slit: '+string(ix,format='(i4.4)'),    $
                                        	charsize=3,/xstyle,xtitle='Wavelength (A)',ytitle='I'
                                ;oplot,wavelength.wvl397[pos],yfit,color=3,thick=3,psym=1
                                ;oplot,replicate(wavelength.wvl397[b[6]],2),!y.crange,line=0,color=3
                                ;oplot,replicate(wavelength.wvl397[b[6]+b[7]],2),!y.crange,line=2,color=3
                                ;oplot,replicate(wavelength.wvl397[b[6]-b[7]],2),!y.crange,line=2,color=3
                                ;oplot,replicate(wavelength.wvl397[b[8]],2),!y.crange,line=0,color=3
                                ;oplot,replicate(wavelength.wvl397[b[8]+b[9]],2),!y.crange,line=2,color=3
                                ;oplot,replicate(wavelength.wvl397[b[8]-b[9]],2),!y.crange,line=2,color=3
                                ;oplot,replicate(wavelength.wvl397[b[10]],2),!y.crange,line=0,color=3
                                ;oplot,replicate(wavelength.wvl397[b[10]+b[11]],2),!y.crange,line=2,color=3
                                ;oplot,replicate(wavelength.wvl397[b[10]-b[11]],2),!y.crange,line=2,color=3
                                ;oplot,replicate(wavelength.wvl397[b[12]],2),!y.crange,line=0,color=3
                                ;oplot,replicate(wavelength.wvl397[b[12]+b[13]],2),!y.crange,line=2,color=3
                                ;oplot,replicate(wavelength.wvl397[b[12]-b[13]],2),!y.crange,line=2,color=3

					plot,wavelength.wvl397,dprofi,psym=1,color=0,background=1,  $
                                        	xr=[3968,3971],  $
                                        	charsize=3,/xstyle,xtitle='Wavelength (A)',ytitle='dI/dlambda'
					oplot,wavelength.wvl397[pos],yfit1,color=3,thick=2
                                	oplot,replicate(wavelength.wvl397[pcen],2),!y.crange,line=2,color=3
                                	oplot,!x.crange,replicate(0,2),line=2,color=3
                                ;oplot,replicate(wavelength.wvl397[pmin],2),!y.crange,line=2,color=5
                                ;oplot,replicate(wavelength.wvl397[b[6] +.5*b[7]],2),!y.crange,line=2,color=5
                                ;oplot,replicate(wavelength.wvl397[b[8] - 4.*b[9]],2),!y.crange,line=2,color=5
                                ;oplot,replicate(wavelength.wvl397[pix_wvl_h1e_r_1],2),!y.crange,line=2,color=5
				;print,pmax, b[6] +.5*b[7], b[8] - 4.*b[9], pix_wvl_h1e_r_1

                                	plot,wavelength.wvl397,profv/profi,yr=[-0.01,0.01],/ystyle,  $
                                        	xr=[3968,3971],  $
                                        	psym=-1,color=0,background=1,  $
						title='Mean V/I = '+string(h1e_ncp[iz,ix,idir],format='(f7.4)'),  $
                                        	charsize=3,/xstyle,xtitle='Wavelength (A)',ytitle='V/I'
                                	oplot,replicate(wavelength.wvl397[pcen],2),!y.crange,line=0,color=3
                                	oplot,replicate(wavelength.wvl397[pcen+width],2),!y.crange,line=2,color=3
                                	oplot,replicate(wavelength.wvl397[pcen-width],2),!y.crange,line=2,color=3
                                	oplot,!x.crange,[0.,0.],color=0

					write_png,outfile,tvrd(/true)
                        	endif
			endfor ; iplot
                endfor ; ix
        endfor ; iz

	outfile=image_dir+'idir'+string(idir,format='(i2.2)')+'.png'

	wdef,0,500,500
	loadct,0
	!P.multi=0
	plot_image,h1e_ncp[*,*,idir],min=-0.005,max=0.005,  $
		xtitle='Scan (pix)', ytitle='Slit (pix)',title='Mean V/I at He line core, +-0.5%',  $
		color=0,background=255,charsize=1.5,/nosquare,ticklen=-0.02

	set_line_color
	pos=where(plot_idir eq idir, npos)
	if npos ge 1 then begin
		;for ipos=0,npos-1 do oplot,[-100,plot_iz[ipos]],[-100,plot_ix[ipos]],psym=1,color=3,symsize=3
		for ipos=0,npos-1 do oplot,plot_iz[pos[ipos]]+[-5,-5,5,5,-5],plot_ix[pos[ipos]]+[-30,30,30,-30,-30],color=3,symsize=3
	endif

	write_png,outfile,tvrd(/true)
endfor


save,h1e_ncp,h1e_vlos,file=save_file

;write_png,image_file,tvrd(/true)
!p.multi=0
loadct,0

END


;=========================================================================
PRO MAKE_AVERAGE_PROFILE,data=data,save_file=save_file
COMMON ANALYSIS

file=file_search(data,'*_I_*')
mreadfits,file,indexi,datai
prof = reform(rebin(datai,1,indexi[0].naxis2,1))
save,prof,file=save_file

END
;=========================================================================

PRO MAKE_CAH_01, save_file=savefile
COMMON ANALYSIS

tmp=min(abs(wavelength.wvl397 - wavelength.wvl_cah3),pos_h3)
tmp=min(abs(wavelength.wvl397 - wavelength.wvl_cah2v),pos_h2v)
tmp=min(abs(wavelength.wvl397 - wavelength.wvl_cah1v),pos_h1v)
tmp=min(abs(wavelength.wvl397 - wavelength.wvl_cahcon),pos_con)

ndir = n_elements(data_dir.dir397)

; make array
nx=0
nz=0
for idir=0,ndir-1 do begin
	file=file_search(data_dir.dir397[idir]+'_cal01','*_I_*',count=nf)
	mreadfits,file[0],index,data,/nodata
	nx = nx > index[0].naxis1
	nz = nz > nf
endfor
cah3 = fltarr(nz,nx,ndir)
cah2v = fltarr(nz,nx,ndir)
cah1v = fltarr(nz,nx,ndir)
con = fltarr(nz,nx,ndir)

; analysis
for idir=0,ndir-1 do begin
	file=file_search(data_dir.dir397[idir]+'_cal01','*_I_*')
	mreadfits,file,index,data

	tmp=size(data) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]

	cah3[0:nz-1,0:nx-1,idir] = transpose(reform(data[*,pos_h3,*]))
	cah2v[0:nz-1,0:nx-1,idir] = transpose(reform(data[*,pos_h2v,*]))
	cah1v[0:nz-1,0:nx-1,idir] = transpose(reform(data[*,pos_h1v,*]))
	con[0:nz-1,0:nx-1,idir] = transpose(reform(data[*,pos_con,*]))
endfor

save,cah3,cah2v,cah1v,con,file=savefile


END
;=========================================================================

PRO MAKE_H1E_I_01, save_file=savefile
COMMON ANALYSIS

tmp=min(abs(wavelength.wvl397 - wavelength.wvl_h1e),pos_h1e)
tmp=min(abs(wavelength.wvl397 - 3970.7),pos_h1e_con)

ndir = n_elements(data_dir.dir397)

; make array
nx=0
nz=0
for idir=0,ndir-1 do begin
        file=file_search(data_dir.dir397[idir]+'_cal01','*_I_*',count=nf)
        mreadfits,file[0],index,data,/nodata
        nx = nx > index[0].naxis1
        nz = nz > nf
endfor
h1e_i = fltarr(nz,nx,ndir)
h1e_i_ic = fltarr(nz,nx,ndir)

; analysis
for idir=0,ndir-1 do begin
        file=file_search(data_dir.dir397[idir]+'_cal01','*_I_*')
        mreadfits,file,index,data

        tmp=size(data) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]


        h1e_i[0:nz-1,0:nx-1,idir] = transpose(reform(data[*,pos_h1e,*]))
        h1e_i_ic[0:nz-1,0:nx-1,idir] = transpose(reform(data[*,pos_h1e,*]/data[*,pos_h1e_con,*]))
endfor

save,h1e_i,h1e_i_ic,file=savefile


END
;=========================================================================
PRO MAKE_H1E_IQUV_01, save_file=savefile
COMMON ANALYSIS

tmp=min(abs(wavelength.wvl397 - wavelength.wvl_h1e),pos_h1e)
tmp=min(abs(wavelength.wvl397 - 3970.7),pos_h1e_con)

ndir = n_elements(data_dir.dir397)

; make array
nx=0
nz=0
for idir=0,ndir-1 do begin
        file=file_search(data_dir.dir397[idir]+'_cal01','*_I_*',count=nf)
        mreadfits,file[0],index,data,/nodata
        nx = nx > index[0].naxis1
        nz = nz > nf
endfor
h1e_i = fltarr(nz,nx,ndir)
h1e_q = fltarr(nz,nx,ndir)
h1e_u = fltarr(nz,nx,ndir)
h1e_v = fltarr(nz,nx,ndir)

; analysis
for idir=0,ndir-1 do begin
        file=file_search(data_dir.dir397[idir]+'_cal01','*_I_*')
        mreadfits,file,index,data
        file=file_search(data_dir.dir397[idir]+'_cal01','*_Q_*')
        mreadfits,file,index,dataq
        file=file_search(data_dir.dir397[idir]+'_cal01','*_U_*')
        mreadfits,file,index,datau
        file=file_search(data_dir.dir397[idir]+'_cal01','*_V_*')
        mreadfits,file,index,datav

        tmp=size(data) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]


        h1e_i[0:nz-1,0:nx-1,idir] = transpose(reform(data[*,pos_h1e,*]))
        h1e_q[0:nz-1,0:nx-1,idir] = transpose(reform(dataq[*,pos_h1e,*]))
        h1e_u[0:nz-1,0:nx-1,idir] = transpose(reform(datau[*,pos_h1e,*]))
        h1e_v[0:nz-1,0:nx-1,idir] = transpose(reform(datav[*,pos_h1e,*]))
endfor

save,h1e_i,h1e_q,h1e_u,h1e_v,file=savefile


END

;=========================================================================

PRO MAKE_FE_V_01, save_file=savefile
COMMON ANALYSIS

tmp=min(abs(wavelength.wvl630 - wavelength.wvl_fev),pos_fev)

ndir = n_elements(data_dir.dir630)

; make array
nx=0
nz=0
for idir=0,ndir-1 do begin
        file=file_search(data_dir.dir630[idir]+'_cal01','*_I_*',count=nf)
        mreadfits,file[0],index,data,/nodata
        nx = nx > index[0].naxis1
        nz = nz > nf
endfor
fe_vi = fltarr(nz,nx,ndir)

; analysis
for idir=0,ndir-1 do begin
        filei=file_search(data_dir.dir630[idir]+'_cal01','*_I_*')
        filev=file_search(data_dir.dir630[idir]+'_cal01','*_V_*')

        mreadfits,filei,index,datai
        mreadfits,filev,index,datav

        tmp=size(datai) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]


        fe_vi[0:nz-1,0:nx-1,idir] = transpose(reform(datav[*,pos_fev,*]/datai[*,pos_fev,*]))
endfor

save,fe_vi,file=savefile


END

;=========================================================================
PRO BLOS2BNORMAL_01, wfa=wfa, save_file=save_file, con=con
COMMON ANALYSIS

restore,wfa
restore,con

;hmidir = [	'/data/tanan/sdo/hmi/20220223/JSOC_20221207_321/',	$
;		'/data/tanan/sdo/hmi/20220223/JSOC_20221207_321/',	$
;		'/data/tanan/sdo/hmi/20220223/JSOC_20221207_321/',	$ ; need update!
;		'/data/tanan/sdo/hmi/20220223/JSOC_20221207_321/'	$ ; need update!
;	]
hmidir = ['~/lib/EB_20220223/demodata/sharp/']

for ii=0,1 do begin
	case ii of
		0:begin
		       	bv = bv630
			con= con630
		end
		1:begin
		       bv = bv854
		       con= con854
		end
	endcase
	;tmp = size(bv.blos) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]
	tmp = size(bv.blos) & nx=tmp[1] & ny=tmp[2] & nz=1

	b  = sqrt(bv.blos^2 + bv.bpos^2)
	cth1= atan(bv.bpos,bv.blos)
	cph1= bv.cph

	thb=fltarr(nx,ny,3,nz)   ; 3rd element) SHARP, SHARP-like smooth solution, opposite of the smooth solution 
	phb=fltarr(nx,ny,3,nz)
	for iz=0,nz-1 do begin
		;for iz=1,1 do begin

		; ==========================================================
		; == READ SDO/HMI/SHARP DATA, 1st data in the data directory
		; AI: not yet for Dec. 27 data
		filebr = file_search(hmidir[iz],'*.Br.*')
		filebt = file_search(hmidir[iz],'*.Bt.*')
		filebp = file_search(hmidir[iz],'*.Bp.*')
		filecn = file_search(hmidir[iz],'*.continuum.*')

		ifile=0
		read_sdo,filebr[ifile],indexbr,databr
		read_sdo,filebt[ifile],indexbt,databt
		read_sdo,filebp[ifile],indexbp,databp
		read_sdo,filecn[ifile],indexcn,cn0
		bpos=sqrt(databt^2 + databp^2)
		b_ref0=sqrt(databr^2 + databt^2 + databp^2)
		thb_ref0 = atan(bpos,databr)
		phb_ref0 = atan(-databt,databp); + cph_offset
		;phb_ref0[where(phb_ref0*!radeg ge 180.)] = phb_ref0[where(phb_ref0*!radeg ge 180.)] - 2.*!pi

		; ========================
		; == Calculate angles
		;pangle = get_rb0p(observation.date,/pangle,/quiet)
		;cph_offset = !pi/2. - (atan(abs(observation.hxy[0]),abs(observation.hxy[1]))+pangle)    ;4.*!dtor
		cph_offset = atan(observation[iz].hxy[1], observation[iz].hxy[0])
		mu = cos(asin(sqrt(observation[iz].hxy[0]^2 + observation[iz].hxy[1]^2)/get_rb0p(observation[iz].date,/radius,/quiet)))
print,mu
;stop
		; ====================================
		; == Convert LOS to solar normal frame
		if 1 then begin
			thb1=fltarr(nx,ny)
			phb1=fltarr(nx,ny)
			thb2=fltarr(nx,ny)
			phb2=fltarr(nx,ny)
			for ix=0,nx-1 do begin
        			for iy=0,ny-1 do begin
					; ta_cthph2thph need the coordinates of the scatter
        			        ta_cthph2thph,cth1[ix,iy,iz],cph1[ix,iy,iz],acos(mu),thba,phba
        		        	ta_cthph2thph,cth1[ix,iy,iz],(cph1[ix,iy,iz] + !pi),acos(mu),thbb,phbb
                			if thba le thbb then begin
                		        	thb1[ix,iy]=thba
                        			phb1[ix,iy]=phba
                        			thb2[ix,iy]=thbb
                        			phb2[ix,iy]=phbb
                			endif else begin
                		        	thb1[ix,iy]=thbb
                		        	phb1[ix,iy]=phbb
                		        	thb2[ix,iy]=thba
                		        	phb2[ix,iy]=phba
                			endelse
        			endfor
			endfor

			save,b,thb1,phb1,thb2,phb2,file='/data/tanan/idl/tmp.sav'
		endif else restore,'/data/tanan/idl/tmp.sav'
		; thb1 < thb2

		; =================================
		; == MAKE REFERENCE FROM SHARP DATA
		case iz of
			0:begin
				rotate_index=0
				rot_angle=-30.
				xcen = 401
				ycen = 87
				if ii eq 0 then begin
					xr = [250,505]-10
					yr = [100,260]-13
				endif else begin
					xr = [240,495]
		                        yr = [106,208]
				endelse
			end
			1:begin
				rotate_index=0
				rot_angle=20-33.
				xcen = 401
				ycen = 87
				if ii eq 0 then begin
					xr = [250,475]-10+30
					yr = [100,260]-13-22
				endif else begin
					xr = [240,475]+25
					yr = [106,208]-21
				endelse
			end
			else:print,'No SHARP data'
		endcase

		;Flip
		cn1 = rotate(cn0,rotate_index)
		b_ref1 = rotate(b_ref0,rotate_index)
		thb_ref1 = rotate(thb_ref0,rotate_index)
		phb_ref1 = rotate(phb_ref0,rotate_index)
		;rotate
		cn2 = rot(cn1,rot_angle,1.,xcen,ycen)
		b_ref2 = rot(b_ref1,rot_angle,1.,xcen,ycen)
		thb_ref2 = rot(thb_ref1,rot_angle,1.,xcen,ycen)
		phb_ref2 = rot(phb_ref1,rot_angle,1.,xcen,ycen)
		;cut & congrid
		cn3 = cn2[xr[0]:xr[1],yr[0]:yr[1]]
		b_ref3 = b_ref2[xr[0]:xr[1],yr[0]:yr[1]]
		thb_ref3 = thb_ref2[xr[0]:xr[1],yr[0]:yr[1]]
		phb_ref3 = phb_ref2[xr[0]:xr[1],yr[0]:yr[1]]
		cn4 = congrid(cn3,nx,ny)
		b_ref4 = congrid(b_ref3,nx,ny)
		thb_ref4 = congrid(thb_ref3,nx,ny)
		phb_ref4 = congrid(phb_ref3,nx,ny)

		;nregion = 20
		;kernel = [	[1.,1.,1.],	$
		;		[-1.,0.,1.],	$
		;		[-1.,-1.,-1.]	$
		;		]
		;threshold=1.0
		;arr0 = abs(convol(phb1,kernel)) ge threshold
		;arr  = morph_close(arr0,replicate(1,5,5))

		;arr0 = abs(convol(phb1,kernel)) le threshold
		;arr  = morph_open(arr0,replicate(1,5,5))

		;label = label_region(arr,/all_neighbors)
		;h = histogram(label, min=0, max=max(label), nbins=max(label), loc=ilabel)
		;sh=reverse(sort(h))
		;label1=fltarr(nx,ny,nregion)
		;for iregion=0,nregion-1 do begin
		;	tmp=fltarr(nx,ny)
		;	tmp[where(label eq sh[iregion])]=1
		;	label1[*,*,iregion] = tmp
		;endfor
		;label2 = fltarr(nx,ny)
                ;for iregion=0,nregion-1 do begin
		;	label2[where(label1[*,*,iregion] eq 1)]=iregion
		;endfor

		thb_god=thb1
		phb_god=phb1
		thb_bad=thb1
		phb_bad=phb1
		for ix=0,nx-1 do begin
			for iy=0,ny-1 do begin
				if ix le 220. or (thb_ref4[ix,iy]*!radeg le 20.) or (thb_ref4[ix,iy]*!radeg ge 160.)then begin
					tmp=min([abs(thb1[ix,iy]-thb_ref4[ix,iy]), abs(thb2[ix,iy]-thb_ref4[ix,iy])], key)
				endif else begin
					tmp=min([abs(phb1[ix,iy]-phb_ref4[ix,iy]), abs(phb2[ix,iy]-phb_ref4[ix,iy])], key)
				endelse

				;tmp=min([abs(thb1[ix,iy]-thb_ref4[ix,iy]) + abs(phb1[ix,iy]-phb_ref4[ix,iy]), $
				;	abs(thb2[ix,iy]-thb_ref4[ix,iy]) + abs(phb2[ix,iy]-phb_ref4[ix,iy])], key)

				; For HMI SHARP
				if key eq 0 then begin
					thb_god[ix,iy]=thb1[ix,iy]
					phb_god[ix,iy]=phb1[ix,iy]
					thb_bad[ix,iy]=thb2[ix,iy]
					phb_bad[ix,iy]=phb2[ix,iy]
				endif else begin
					thb_god[ix,iy]=thb2[ix,iy] 
					phb_god[ix,iy]=phb2[ix,iy] 
					thb_bad[ix,iy]=thb1[ix,iy]
					phb_bad[ix,iy]=phb1[ix,iy]
				endelse
				thb[ix,iy,0,iz] = thb_god[ix,iy]
				phb[ix,iy,0,iz] = phb_god[ix,iy]

				; Modify SHARP TO MAKE IT SMOOTH
				thb[ix,iy,1,iz] = thb_god[ix,iy] ; thb1[ix,iy]
				phb[ix,iy,1,iz] = phb_god[ix,iy] ; phb1[ix,iy]
				thb[ix,iy,2,iz] = thb_bad[ix,iy] ; thb2[ix,iy]
				phb[ix,iy,2,iz] = phb_bad[ix,iy] ; phb2[ix,iy]


				; This is the case when SDO/HMI/SHARP missolving 180 degree ambiguity 	
				if ii eq 0 and iz eq 0 then begin
					xmax = 150
					ymin = 650
				endif
				if (ii eq 1 and iz eq 0) then begin
					xmax = 400 ;250
					ymin = 650
				endif
				if (ii eq 0 and iz eq 1) then begin
					xmax = 220
					ymin = 900
				endif
				if (ii eq 1 and iz eq 1) then begin
					xmax = 220
					ymin = 900
				endif
				if ix le xmax and iy ge ymin then begin
				       thb[ix,iy,1,iz]=thb1[ix,iy]	
				       phb[ix,iy,1,iz]=phb1[ix,iy]
				       thb[ix,iy,2,iz]=thb2[ix,iy]
				       phb[ix,iy,2,iz]=phb2[ix,iy]
				endif
			endfor
		endfor
	endfor



	case ii of
		0:begin
			b630 = b
			thb630 = thb
			phb630 = phb
		end
		1:begin
			b854 = b
			thb854 = thb
			phb854 = phb
		end
	endcase
endfor

save,b630,thb630,phb630,b854,thb854,phb854,file=save_file

END
;=========================================================================
;========================================================================
;=========================================================================
PRO MAKE_CON_01, save_file=save_file
COMMON ANALYSIS

for iarm=0,2 do begin
	case iarm of
		0:begin
			dir=data_dir.dir397
			tmp=min(abs(wavelength.wvl397 - wavelength.wvl_cahcon),pos)
		end
		1:begin
			dir=data_dir.dir630
			tmp=min(abs(wavelength.wvl630 - wavelength.wvl_fecon),pos)
		end
		2:begin
			dir=data_dir.dir854
			tmp=min(abs(wavelength.wvl854 - wavelength.wvl_caicon),pos)
		end
	endcase

	ndir=n_elements(dir)
	nz=0
	nx=0
	for idir=0,ndir-1 do begin
                file=file_search(dir[idir]+'_cal01','*_I_*',count=nf)
		nz=nz>nf
		mreadfits,file[0],index,data,/nodata
		nx=nx>index[0].naxis1
	endfor
	data1=fltarr(nz,nx,ndir)

	for idir=0,ndir-1 do begin
		file=file_search(dir[idir]+'_cal01','*_I_*',count=nf)
		mreadfits,file,index,data
		data1[0:nf-1,0:index[0].naxis1-1,idir] = transpose(reform(data[*,pos,*]))
	endfor

	case iarm of
		0:con397=data1
		1:con630=data1
		2:con854=data1
	endcase
endfor

save,con397,con630,con854,file=save_file


END
;-----------------------------------------------
PRO CAL_EARTH_ATM_DISPERSION, cont_file=cont_file, savefile=savefile, img_dir=img_dir
COMMON ANALYSIS

restore,cont_file,/verb
if is_dir(img_dir) eq 0 then spawn,'mkdir '+img_dir



wdef,0,1200,800
!p.multi=[0,2,2]
chs=1.5

;nf = (size(con397))[3]
nf = 1
res= fltarr(3,nf)  ; pix 
for i=0,nf-1 do begin
	cn397 = reform(con397[*,hairline.pos397[0]:hairline.pos397[1],i])
	cn630 = reform(con630[*,hairline.pos630[0]:hairline.pos630[1],i])
	cn854 = reform(con854[*,hairline.pos854[0]:hairline.pos854[1],i])


	; 630
	iline = 1
	cn  = cn630
	tmp = size(cn) & nx = tmp[1] & ny=tmp[2]
	ref = congrid(cn397, nx, ny)
	img = congrid(cn, nx, ny)
	;off = get_correl_offsets([[[ref]],[[img]]]) ; does not work

	nlag= 400
	lag = findgen(nlag) - nlag/2.  & lag = lag[where(abs(lag) ge 10)]
	nlag= n_elements(lag)
	cc  = fltarr(nlag)
	for j=0,nlag-1 do begin
		alag=abs(lag[j])
		cc[j] = correlate(ref[*,alag:ny-alag-1], (shift_img(img, [0.,lag[j]]))[*,alag:ny-alag-1] )
	endfor
	tmp = max(cc,pos)
	res[iline,i] = lag[pos]

	levels = [0.5,1.0]
	levels = [0.5,0.8]
	loadct,0
	plot_image,ref>0<1.5,/nosquare,ticklen=-0.02,  $
		color=0,background=255,charsize=chs,  $
		xtitle='X (pix)',ytitle='Y (pix)', title='Reference (397 nm)'
	set_line_color
	contour,img,/over,color=3,levels=levels
	contour,shift_img(img, [0,res[iline,i]]),/over,color=5,levels=levels

        loadct,0
        plot_image,img>0<2,/nosquare,ticklen=-0.02,  $
                color=0,background=255,charsize=chs,  $
                xtitle='X (pix)',ytitle='Y (pix)', title='630 nm'
        set_line_color
        contour,img,/over,color=3,levels=levels
        contour,shift_img(img,[0,res[iline,i]]),/over,color=5,levels=levels

        ; 854
	iline = 2
        cn  = cn854
        tmp = size(cn) & nx = tmp[1] & ny=tmp[2]
        ref = congrid(cn397, nx, ny)
        img = congrid(cn, nx, ny)
        ;off = get_correl_offsets([[[ref]],[[img]]]) ; does not work

        ;nlag= 400
        ;lag = findgen(nlag) - nlag/2.  & lag = lag[where(abs(lag) ge 10)]
	lag = findgen(100)+20
	nlag= n_elements(lag)
        cc  = fltarr(nlag)
        for j=0,nlag-1 do begin
                alag=abs(lag[j])
                cc[j] = correlate(ref[*,alag:ny-alag-1], (shift_img(img, [0.,lag[j]]))[*,alag:ny-alag-1] )
        endfor
        tmp = max(cc,pos)
        res[iline,i] = lag[pos]


	levels=[0.55,0.8]
        loadct,0
        plot_image,ref>0<1.5,/nosquare,ticklen=-0.02,  $
                color=0,background=255,charsize=chs,  $
                xtitle='X (pix)',ytitle='Y (pix)', title='Reference (397 nm)'
        set_line_color
        contour,img,/over,color=3,levels=levels
        contour,shift_img(img, [0.,res[iline,i]]),/over,color=5,levels=levels

        loadct,0
        plot_image,img>0<1.5,/nosquare,ticklen=-0.02,  $
                color=0,background=255,charsize=chs,  $
                xtitle='X (pix)',ytitle='Y (pix)', title='854 nm'
        set_line_color
        contour,img,/over,color=3,levels=levels
        contour,shift_img(img, [0,res[iline,i]]),/over,color=5,levels=levels

	set_line_color
	chs2=2
	iline = 1
	xyouts,/norm,0.25,0.63,charsize=chs2,color=3,'Contour for 630 nm image'
	xyouts,/norm,0.25,0.60,charsize=chs2,color=5,'Shifted by '+string(index_for_map.index_nosun_630.cdelt2*res[iline,i],format='(f4.2)')+' arcsec'

	xyouts,/norm,0.75,0.63,charsize=chs2,color=3,'Contour for 630 nm image'
	xyouts,/norm,0.75,0.60,charsize=chs2,color=5,'Shifted by '+string(index_for_map.index_nosun_630.cdelt2*res[iline,i],format='(f4.2)')+' arcsec'

	iline = 2
	xyouts,/norm,0.25,0.13,charsize=chs2,color=3,'Contour for 854 nm image'
	xyouts,/norm,0.25,0.10,charsize=chs2,color=5,'Shifted by '+string(index_for_map.index_nosun_854.cdelt2*res[iline,i],format='(f4.2)')+' arcsec'

	xyouts,/norm,0.75,0.13,charsize=chs2,color=3,'Contour for 854 nm image'
	xyouts,/norm,0.75,0.10,charsize=chs2,color=5,'Shifted by '+string(index_for_map.index_nosun_854.cdelt2*res[iline,i],format='(f4.2)')+' arcsec'

	write_png,img_dir+'/'+string(i,format='(i2.2)')+'.png',tvrd(/true)
endfor
save,res,file=savefile

print,res

END

;=========================================================================
;=========================================================================
PRO MAKE_DATA_08, obs_dir = obs_dir, save_dir = save_dir, format=format, no_hepsilon=no_hepsilon
COMMON ANALYSIS
; SAVE DATA FOR DESIRE INVERSION
; CALIBRATE STOKES-Q IN 854 nm
; CALIBRATE TILT OF CAIIH
; no_hepsilon True: Fit H epsilon Stokes-I profile with a Gaussian and remove it

if keyword_set(format) eq 0 then format='txt'
if keyword_set(no_hepsilon) eq 0 then no_hepsilon=0
if is_dir(save_dir) eq 0 then spawn,'mkdir '+save_dir
obs_file = file_search(obs_dir,'*.sav', count=nf)

for i=0,nf-1 do begin
;i = 7
    	save_dir1 = save_dir+strmid(obs_file[i],strlen(obs_dir),strlen(obs_file[i])-strlen(obs_dir)-4)+'/'
	if is_dir(save_dir1) eq 0 then spawn,'mkdir '+save_dir1

        restore,obs_file[i]
;save,   stki_630,stkq_630,stku_630,stkv_630,wave_630,   $
;        stki_397,stkq_397,stku_397,stkv_397,wave_397,   $
;        stki_854,stkq_854,stku_854,stkv_854,wave_854,   $
;        b_los_xyz,h1e_ncp,      $
;        ibomb,index,dd, $
;        file=save_file

        tmp=size(stkq_854) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]

        ; calibrate Stokes-Q in 854 nm
        prof = reform(rebin(stkq_854,1,1,nz))
        for ix=0,nx-1 do begin
                for iy=0,ny-1 do begin
                        stkq_854[ix,iy,*] = stkq_854[ix,iy,*] - prof
                endfor
        endfor

	; Fit H epsilon Stokes-I profile with a Gaussian and remove it
        if no_hepsilon then begin
                tmp=size(stkq_397) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]
                pos = where(    (wave_397 ge wavelength.wvl_h1e_r1[0] and wave_397 le wavelength.wvl_h1e_r1[1]) or $
                                (wave_397 ge wavelength.wvl_h1e_r1[2] and wave_397 le wavelength.wvl_h1e_r1[3]) or $
                                (wave_397 ge wavelength.wvl_h1e_r1[4] and wave_397 le wavelength.wvl_h1e_r1[5]) )
                pos1 = where(   (wave_397 ge wavelength.wvl_h1e_r1[0] and wave_397 le wavelength.wvl_h1e_r1[1]) or $
                                (wave_397 ge wavelength.wvl_h1e_r1[2] and wave_397 le (wavelength.wvl_h1e_r1[2]+0.01)) or $
                                (wave_397 ge wavelength.wvl_h1e_r1[4] and wave_397 le wavelength.wvl_h1e_r1[5]) )

                tmp = min(abs(wave_397 - wavelength.wvl_h1e),wvl0_pix)
                stki_397_ori = stki_397
                stkq_397_ori = stkq_397
                stku_397_ori = stku_397
                stkv_397_ori = stkv_397

                qi = stkq_397_ori/stki_397_ori
                ui = stku_397_ori/stki_397_ori
                vi = stkv_397_ori/stki_397_ori
                for ix=0,nx-1 do begin
;ix = 25
                        for iy=0,ny-1 do begin
;iy = 30
                                prof = stki_397[ix,iy,*]

                                y = prof[pos1]
                                coe = poly_fit(pos1,y,2,yfit=yfit)
                                ;plot,y,psym=1
                                ;oplot,yfit
                                ;plot,prof
                                ;oplot,pos1,prof[pos1],psym=1
                                ;tmpx = findgen(n_elements(prof))
                                ;oplot,tmpx,coe[0]+coe[1]*tmpx+coe[2]*tmpx*tmpx

                                y1= prof[pos] - (coe[0]+coe[1]*pos+coe[2]*pos^2)
                                estimates = [-0.1,wvl0_pix,5]
                                yfit = gaussfit(pos,y1,coe1, nterms=3, estimates=estimates)
                                z = (findgen(nz) - coe1[1])/coe1[2]
                                yfit1= coe1[0]*exp(-0.5*z^2) + coe[0]+coe[1]*findgen(nz)+coe[2]*findgen(nz)^2
                                prof1=prof-coe1[0]*exp(-0.5*z^2)

;                               plot,pos,y1,psym=1
;                               oplot,pos,yfit
;                               plot,prof
;                               oplot,prof1
                                ;stop
                                stki_397[ix,iy,*] = prof1
                                stkq_397[ix,iy,*] = qi[ix,iy,*]*stki_397[ix,iy,*]
                                stku_397[ix,iy,*] = ui[ix,iy,*]*stki_397[ix,iy,*]
                                stkv_397[ix,iy,*] = vi[ix,iy,*]*stki_397[ix,iy,*]

;stop
                        endfor
                endfor
        endif

	; CALIBRATION FOR TILT IN CA II H STOKES-I PROFILE, BUT THE TILT IS NOT CLEAR IN AVERAGED PROFILE
	tmp = size(stki_397) & nx=tmp[1] & ny=tmp[2] & nw=tmp[3]
	prof = reform(rebin(stki_397,1,1,nw))
	wvl  = wave_397
	wstep = 0.01 ;wvl[1]-wvl[0]
	wmin  = min(wvl)
	wmax  = max(wvl) ;wmin + wstep*(n_elements(wvl))
	atlas0 = ta_atlas(wmin, wmax, wstep,wl=wvl_atlas)/255.
	atlas  = interpol(atlas0,wvl_atlas,wvl)
	pos    = where(wvl le 3967.5 or wvl ge 3969.5)
	coe    = poly_fit(atlas[pos],prof[pos],1)
	atlas1 = coe[0]+coe[1]*atlas
	

        ; SAVE TXT or FITS
        if format eq 'txt' then begin
                for ix=0,nx-1 do begin
                        print,ix
                        for iy=0,ny-1 do begin
                                outfile = save_dir1+'x'+string(ix,format='(i3.3)')+'_y'+string(iy,format='(i3.3)')+'.txt'
                                openw,lun,outfile,/get_lun
                                for iline = 0,2 do begin
                                        case iline of
                                                0:begin
                                                        wave=wave_397
                                                        stki=stki_397
                                                        stkq=stkq_397
                                                        stku=stku_397
                                                        stkv=stkv_397
                                                end
                                                1:begin
                                                        wave=wave_630
                                                        stki=stki_630
                                                        stkq=stkq_630
                                                        stku=stku_630
                                                        stkv=stkv_630
                                                end
                                                2:begin
                                                        wave=wave_854
                                                        stki=stki_854
                                                        stkq=stkq_854
                                                        stku=stku_854
                                                        stkv=stkv_854
                                                end
                                        endcase
                                        printf,lun,strcompress(string(n_elements(wave)),/remove_all)
                                        for il=0,n_elements(wave)-1 do begin
                                                printf,lun,string(wave[il],format='(f9.4)') + ' ' + $
                                                        string(stki[ix,iy,il],format='(f9.6)') + ' ' + $
                                                        string(stkq[ix,iy,il],format='(f9.6)') + ' ' + $
                                                        string(stku[ix,iy,il],format='(f9.6)') + ' ' + $
                                                        string(stkv[ix,iy,il],format='(f9.6)') + ' '
                                        endfor
                                endfor
                                free_lun,lun
                        endfor
                endfor
        endif


        if format eq 'fits' then begin
                for iline=0,2 do begin
                        case iline of
                                0:begin
                                        outfile = save_dir1+'Ca397nm_noise.fits'
                                        outfile1= save_dir1+'Ca397nm_wave.fits'
                                        stki = stki_397
                                        stkq = stkq_397
                                        stku = stku_397
                                        stkv = stkv_397
                                        wave = wave_397
                                end
                                1:begin
                                        outfile = save_dir1+'Fe630nm_noise.fits'
                                        outfile1 = save_dir1+'Fe630nm_wave.fits'
                                        stki = stki_630
                                        stkq = stkq_630
                                        stku = stku_630
                                        stkv = stkv_630
                                        wave = wave_630
                                end
                                2:begin
                                        outfile = save_dir1+'Ca854nm_noise.fits'
                                        outfile1 = save_dir1+'Ca854nm_wave.fits'
                                        stki = stki_854
                                        stkq = stkq_854
                                        stku = stku_854
                                        stkv = stkv_854
                                        wave = wave_854
                                end
                        endcase

                        tmp = size(stki) & nx=tmp[1] & ny=tmp[2] & nl=tmp[3]
                        sn = fltarr(nx,ny,nl,4)
                        sn[*,*,*,0] = stki
                        sn[*,*,*,1] = stkq
                        sn[*,*,*,2] = stku
                        sn[*,*,*,3] = stkv
                        writefits,outfile,sn
                        writefits,outfile1,wave
                endfor
        endif

ENDFOR ; i

END




;=========================================================================
;=========================================================================
PRO MAKE_NCP_MAP_01, image_dir = image_dir, save_file = save_file
COMMON ANALYSIS

if is_dir(image_dir) eq 0 then spawn,'mkdir '+image_dir

keys = ['cah','hepsilon','fe_a','fe_b','cair']
ncp_parameters0 = {key:'', wavelength_low: 0., wavelength_high: 0., plot_idir:0, plot_iz:135, plot_ix:800}
ncp_parameters  = replicate(ncp_parameters0, 5)
ncp_parameters[*].plot_iz = 222  & ncp_parameters[*].plot_ix = 1350 ; umbra

ncp_parameters[0].key             = 'cah'
ncp_parameters[0].wavelength_low  = 3967.6
ncp_parameters[0].wavelength_high = 3969.
ncp_parameters[1].key             = 'hepsilon'
ncp_parameters[1].wavelength_low  = 3969.85
ncp_parameters[1].wavelength_high = 3970.7
ncp_parameters[2].key             = 'fe_a'
ncp_parameters[2].wavelength_low  = 6300.
ncp_parameters[2].wavelength_high = 6302.
ncp_parameters[3].key             = 'fe_b'
ncp_parameters[3].wavelength_low  = 6302.
ncp_parameters[3].wavelength_high = 6303.
ncp_parameters[4].key             = 'cair'
ncp_parameters[4].wavelength_low  = 8540.
ncp_parameters[4].wavelength_high = 8545.

for iwave = 0,2 do begin
	case iwave of
		0:begin
			dir = data_dir.dir397
			wvl = wavelength.wvl397
			key = ['cah','hepsilon']
			hair = hairline.pos397			
		end
		1:begin
		      	dir = data_dir.dir630
			wvl = wavelength.wvl630
			key = ['fe_a','fe_b']
			hair = hairline.pos630			
		end
		2:begin
			dir = data_dir.dir854
			wvl = wavelength.wvl854
			key = ['cair']
			hair = hairline.pos854			
		end
	endcase
	ndir=n_elements(dir)
	for ikey0=0,n_elements(key)-1 do begin
		ikey=where(keys eq key[ikey0])
		for idir=0,ndir-1 do begin
        		if 1 then begin
                		file=file_search(dir[idir]+'_cal01','*_I_*')
                		mreadfits,file,indexi,datai
                		file=file_search(dir[idir]+'_cal01','*_V_*')
                		mreadfits,file,indexv,datav
                		save,indexi,datai,indexv,datav,file='/data/tanan/idl/tmp.sav'
        		endif else restore,'/data/tanan/idl/tmp.sav'
        		tmp = size(datai) & nx=tmp[1] & ny=tmp[2] & nz=tmp[3]
        		x = findgen(ny)
        		if idir eq 0 then ncp=fltarr(nz,nx,ndir)
        		for iz=0,nz-1 do begin
;       for iz=130,140 do begin
                		for ix=0,nx-1 do begin
                        		profi = reform(datai[ix,*,iz])
                        		profv = reform(datav[ix,*,iz])
                        		pos = where(wvl ge ncp_parameters[ikey].wavelength_low and wvl le ncp_parameters[ikey].wavelength_high)
					ncp[iz,(ix)<((size(ncp))[2]-1),idir] = int_tabulated(wvl[pos],profv[pos]/profi[pos])

	                        	if idir eq ncp_parameters[ikey].plot_idir and $
						iz eq ncp_parameters[ikey].plot_iz and $
						ix eq ncp_parameters[ikey].plot_ix then begin

						outfile = image_dir+key[ikey0]+'_prof.png' 
						wdef,0,1200,800
        	                        	set_line_color
                	                	!p.multi=[0,1,2]
                        	        	plot,wvl,profi,psym=1,color=0,background=1,  $
                                		        charsize=3,/xstyle,xtitle='Wavelength (A)',ytitle='I/Ic'
                                		oplot,replicate(ncp_parameters[ikey].wavelength_low,2),!y.crange,line=2,color=0
                                		oplot,replicate(ncp_parameters[ikey].wavelength_high,2),!y.crange,line=2,color=0

                                		plot,wvl,profv/profi,yr=[-0.1,0.1],/ystyle,  $
                                        		psym=-1,color=0,background=1,  $
                                        		charsize=3,/xstyle,xtitle='Wavelength (A)',ytitle='V/I'
                                		oplot,replicate(ncp_parameters[ikey].wavelength_low,2),!y.crange,line=2,color=0
                                		oplot,replicate(ncp_parameters[ikey].wavelength_high,2),!y.crange,line=2,color=0
                                		oplot,!x.crange,[0.,0.],color=0
						
						write_png,outfile,tvrd(/true)

                        		endif
                		endfor
			endfor ; iz

			outfile = image_dir+key[ikey0]+'_idir'+string(idir,format='(i2.2)')+'.png'
			imean = mean(ncp[*,hair[0]+50:hair[1]-50,idir])
			istd  = stddev(ncp[*,hair[0]+50:hair[1]-50,idir])
			imin  = imean - 3.*istd 
			imax  = imean + 3.*istd 

			wdef,0,600,600
			!p.multi=0
			loadct,0
			plot_image,ncp[*,*,idir],/nosquare,min=imin,max=imax,ticklen=-0.02,  $
				color=0,background=255,title=key[ikey0]
			write_png,outfile,tvrd(/true)

		endfor ; idir
 	
		case key[ikey0] of
			'cah':begin
				ncp_cah = ncp
			end
			'hepsilon':begin	
				ncp_he = ncp
			end
			'fe_a':begin
				ncp_fe6301 = ncp
			end
			'fe_b':begin
				ncp_fe6302 = ncp
			end
			'cair':begin
				ncp_ca854 = ncp
			end
		endcase

	endfor ; ikey
endfor


save,ncp_cah,ncp_he,ncp_fe6301,ncp_fe6302,ncp_ca854,file=save_file



END

;=========================================================================
PRO CHECK_CROSSTALK_01, image_dir=image_dir, savefile=savefile
COMMON ANALYSIS

if is_dir(image_dir) eq 0 then spawn,'mkdir '+image_dir
;local_parameters = {	yl397:770,  $
;			yt397:820,  $
local_parameters = {	yl397:200,  $
			yt397:270,  $
			yl630:530,  $
			yt630:600,  $
			yl854:500,  $
			yt854:600 }
wdef,0,800,600
wdef,1,600,600
for iline=0,2 do begin
	case iline of
		0:begin
			title='397 nm'			
			dir = data_dir.dir397
			yl  = local_parameters.yl397
			yt  = local_parameters.yt397
		end
		1:begin
			title='630 nm'			
			dir = data_dir.dir630		
			yl  = local_parameters.yl630
			yt  = local_parameters.yt630
		end
		2:begin
			title='854 nm'			
			dir = data_dir.dir854		
			yl  = local_parameters.yl854
			yt  = local_parameters.yt854
		end
	endcase
	ndir = n_elements(dir)
	for idir=0,ndir-1 do begin
		file1=file_search(dir[idir],'*_I_*.fits',count=nf)
		file2=file_search(dir[idir],'*_Q_*.fits',count=nf)
		file3=file_search(dir[idir],'*_U_*.fits',count=nf)
		file4=file_search(dir[idir],'*_V_*.fits',count=nf)
		if idir eq 0 then begin
			res = {	r_iq	:fltarr(ndir,nf),	$
				coe0_iq	:fltarr(ndir,nf),	$
				coe1_iq	:fltarr(ndir,nf),	$
				r_iu	:fltarr(ndir,nf),	$
				coe0_iu	:fltarr(ndir,nf),	$
				coe1_iu	:fltarr(ndir,nf),	$
				r_iv	:fltarr(ndir,nf),	$
				coe0_iv	:fltarr(ndir,nf),	$
				coe1_iv	:fltarr(ndir,nf),	$
				r_qv	:fltarr(ndir,nf),	$
				coe0_qv	:fltarr(ndir,nf),	$
				coe1_qv	:fltarr(ndir,nf),	$
				r_uv	:fltarr(ndir,nf),	$
				coe0_uv	:fltarr(ndir,nf),	$
				coe1_uv	:fltarr(ndir,nf)	$
				}
		endif
		for ifile=0,nf-1 do begin
	        	read_sdo,file1[ifile],index,data1, /use_shared_lib
        		read_sdo,file2[ifile],index,data2, /use_shared_lib
        		read_sdo,file3[ifile],index,data3, /use_shared_lib
        		read_sdo,file4[ifile],index,data4, /use_shared_lib

			nmod = 40
			if ifile mod nmod eq 0 then begin; == PLOT ==;			
				outfile = image_dir+'plot_'+string(iline,format='(i1.1)')+'_'+string(idir,format='(i1.1)')+'_'+string(ifile,format='(i3.3)')+'.png'
				chs=2.5
				wset,0
				!p.multi=[0,3,2]
				loadct,0
				plot_image,data1,/nosquare,title=title,  $
					ticklen=-0.02,color=0,background=255,charsize=chs
				set_line_color
				oplot,!x.crange,replicate(yl,2),color=3
				oplot,!x.crange,replicate(yt,2),color=3
			endif

			xtitle='Stokes I'
			x = data1
			ytitle='Stokes Q'
			y = data2
			r = correlate(x[*,yl:yt],y[*,yl:yt])
			coe = poly_fit(x[*,yl:yt],y[*,yl:yt],1)
			if ifile mod nmod eq 0 then begin; == PLOT ==;			
				plot,x[*,yl:yt],y[*,yl:yt],psym=3,  $
					xtitle=xtitle,ytitle=ytitle,title=title+', R='+string(r,format='(f5.2)'),  $
					charsize=chs,color=0,background=1
				oplot,[-2.,2.],coe[0]+coe[1]*[-2.,2.],color=3
			endif
			res.r_iq[idir,ifile] 	= r
			res.coe0_iq[idir,ifile] = coe[0]
			res.coe1_iq[idir,ifile] = coe[1]

			xtitle='Stokes I'
			x = data1
			ytitle='Stokes U'
			y = data3
			r = correlate(x[*,yl:yt],y[*,yl:yt])
			coe = poly_fit(x[*,yl:yt],y[*,yl:yt],1)
			if ifile mod nmod eq 0 then begin; == PLOT ==;			
				plot,x[*,yl:yt],y[*,yl:yt],psym=3,  $
					xtitle=xtitle,ytitle=ytitle,title=title+', R='+string(r,format='(f5.2)'),  $
					charsize=chs,color=0,background=1
				oplot,[-2.,2.],coe[0]+coe[1]*[-2.,2.],color=3
			endif
			res.r_iu[idir,ifile] 	= r
			res.coe0_iu[idir,ifile] = coe[0]
			res.coe1_iu[idir,ifile] = coe[1]

			xtitle='Stokes I'
			x = data1
			ytitle='Stokes V'
			y = data4
			r = correlate(x[*,yl:yt],y[*,yl:yt])
			coe = poly_fit(x[*,yl:yt],y[*,yl:yt],1)
			if ifile mod nmod eq 0 then begin; == PLOT ==;			
				plot,x[*,yl:yt],y[*,yl:yt],psym=3,  $
					xtitle=xtitle,ytitle=ytitle,title=title+', R='+string(r,format='(f5.2)'),  $
					charsize=chs,color=0,background=1
				oplot,[-2.,2.],coe[0]+coe[1]*[-2.,2.],color=3
			endif
			res.r_iv[idir,ifile] 	= r
			res.coe0_iv[idir,ifile] = coe[0]
			res.coe1_iv[idir,ifile] = coe[1]

			xtitle='Stokes Q'
			x = data2
			ytitle='Stokes V'
			y = data3
			r = correlate(x[*,yl:yt],y[*,yl:yt])
			coe = poly_fit(x[*,yl:yt],y[*,yl:yt],1)
			if ifile mod nmod eq 0 then begin; == PLOT ==;			
				plot,x[*,yl:yt],y[*,yl:yt],psym=3,  $
					xtitle=xtitle,ytitle=ytitle,title=title+', R='+string(r,format='(f5.2)'),  $
					charsize=chs,color=0,background=1
				oplot,[-2.,2.],coe[0]+coe[1]*[-2.,2.],color=3
			endif
			res.r_qv[idir,ifile] 	= r
			res.coe0_qv[idir,ifile] = coe[0]
			res.coe1_qv[idir,ifile] = coe[1]

			xtitle='Stokes U'
			x = data3
			ytitle='Stokes V'
			y = data4
			r = correlate(x[*,yl:yt],y[*,yl:yt])
			coe = poly_fit(x[*,yl:yt],y[*,yl:yt],1)
			if ifile mod nmod eq 0 then begin; == PLOT ==;			
				plot,x[*,yl:yt],y[*,yl:yt],psym=3,  $
					xtitle=xtitle,ytitle=ytitle,title=title+', R='+string(r,format='(f5.2)'),  $
					charsize=chs,color=0,background=1
				oplot,[-2.,2.],coe[0]+coe[1]*[-2.,2.],color=3
			endif
			res.r_uv[idir,ifile] 	= r
			res.coe0_uv[idir,ifile] = coe[0]
			res.coe1_uv[idir,ifile] = coe[1]

			if ifile mod nmod eq 0 then begin; == PLOT ==;			
				write_png,outfile,tvrd(/true)
			endif
		endfor

		for iplot=0,2 do begin
			case iplot of
				0:begin
					title1 = 'R'
					outfile = image_dir+'time_r_'+string(iline,format='(i1.1)')+'_'+string(idir,format='(i1.1)')+'.png'
					y1 = res.r_iq[idir,*]
					y2 = res.r_iu[idir,*]
					y3 = res.r_iv[idir,*]
					y4 = res.r_qv[idir,*]
					y5 = res.r_uv[idir,*]
				end
				1:begin
					title1 = 'Coe 0'
					outfile = image_dir+'time_coe0_'+string(iline,format='(i1.1)')+'_'+string(idir,format='(i1.1)')+'.png'
					y1 = res.coe0_iq[idir,*]
					y2 = res.coe0_iu[idir,*]
					y3 = res.coe0_iv[idir,*]
					y4 = res.coe0_qv[idir,*]
					y5 = res.coe0_uv[idir,*]
				end
				2:begin
					title1 = 'Coe 1'
					outfile = image_dir+'time_coe1_'+string(iline,format='(i1.1)')+'_'+string(idir,format='(i1.1)')+'.png'
					y1 = res.coe1_iq[idir,*]
					y2 = res.coe1_iu[idir,*]
					y3 = res.coe1_iv[idir,*]
					y4 = res.coe1_qv[idir,*]
					y5 = res.coe1_uv[idir,*]
				end
			endcase

			xs = 0.2
			xe = 0.95
			ys = 0.1
			ye = 0.95
			ydd= 0.02
			yd = (ye-ys-4*ydd)/5.

			wset,1
			!p.multi=[0,1,5]
			set_line_color
			iy=4
			if iplot eq 0 then yr=[-1,1] else yr=minmax(y1)
			plot,y1,/xstyle,yr=yr,/ystyle,  $
				color=0,background=1,   $
				xtitle=' ',xtickname=replicate(' ',10),  $
				noerase=0,norm=1,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd],  $
				charsize=chs,ytitle=title1+' I-Q',title=title
			iy=3
			if iplot eq 0 then yr=[-1,1] else yr=minmax(y2)
			plot,y2,/xstyle,yr=yr,/ystyle,  $
				color=0,background=1,   $
				xtitle=' ',xtickname=replicate(' ',10),  $
				noerase=1,norm=1,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd],  $
				charsize=chs,ytitle=title1+' I-U',title=' '
			iy=2
			if iplot eq 0 then yr=[-1,1] else yr=minmax(y3)
			plot,y3,/xstyle,yr=yr,/ystyle,  $
				color=0,background=1,   $
				xtitle=' ',xtickname=replicate(' ',10),  $
				noerase=1,norm=1,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd],  $
				charsize=chs,ytitle=title1+' I-V',title=' '
			iy=1
			if iplot eq 0 then yr=[-1,1] else yr=minmax(y4)
			plot,y4,/xstyle,yr=yr,/ystyle,  $
				color=0,background=1,   $
				xtitle=' ',xtickname=replicate(' ',10),  $
				noerase=1,norm=1,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd],  $
				charsize=chs,ytitle=title1+' Q-V',title=' '
			iy=0
			if iplot eq 0 then yr=[-1,1] else yr=minmax(y5)
			plot,y5,/xstyle,yr=yr,/ystyle,  $
				color=0,background=1,   $
				noerase=1,norm=1,pos=[xs,ys+iy*(yd+ydd),xe,ys+iy*(yd+ydd)+yd],  $
				charsize=chs,xtitle='Scan (# of data)',ytitle=title1+' U-V',title=' '
			write_png,outfile,tvrd(/true)
		endfor
	endfor	
	case iline of
		0: res397 = res
		1: res630 = res
		2: res854 = res
	endcase
endfor

save,res397, res630, res854, local_parameters,file=savefile

loadct,0
!p.multi=0


END
;=========================================================================
;=========================================================================
PRO MAKE_H1E_VLOS_01, save_file=save_file
COMMON ANALYSIS

cc = 2.99792458d+5 ; km/s

fit_wvl_pix = [680+findgen(700-680), 730+findgen(780-730), 845+findgen(870-845), 940+findgen(960-940)]

ndir = n_elements(data_dir.dir397)

; for make array
nz=0
nx=0
for idir=0,ndir-1 do begin
	dir=data_dir.dir397[idir]+'_cal01/'
	filei=file_search(dir,'*_I_*.fits',count=nf)
	nz=nz>nf
	mreadfits,filei[0],index,datai,/nodata
	nx=nx>index[0].naxis1
endfor

vlos = fltarr(nz,nx,ndir)

for idir=0,ndir-1 do begin
	wvl=wavelength.wvl397
        dir=data_dir.dir397[idir]+'_cal01/'

        filei=file_search(dir,'*_I_*.fits',count=nf)
        fileq=file_search(dir,'*_Q_*.fits',count=nf)
        fileu=file_search(dir,'*_U_*.fits',count=nf)
        filev=file_search(dir,'*_V_*.fits',count=nf)

        ;ifile0=230
        ;ifile0=135
        ifile0=0
        for ifile=ifile0,nf-1 do begin
		mreadfits,filei[ifile],index,datai
                mreadfits,fileq[ifile],index,dataq
                mreadfits,fileu[ifile],index,datau
                mreadfits,filev[ifile],index,datav
                nx0 = index.naxis1
                ny0 = index.naxis2

                x0 = findgen(ny0)
                wvl_center_pix = interpol(x0,wvl[0:ny0-1],wavelength.wvl_h1e)
                for ix=0,nx0-1 do begin
                        ;VLOS
                        y = reform(datai[ix,fit_wvl_pix])
                        estimates=[-0.5,wvl_center_pix,5.,1.,0.]
                        yfit = gaussfit(fit_wvl_pix,y,a,estimates=estimates,nterm=5)
                        wvl_center_pix_fit = interpol(wvl[0:ny0-1],x0,a[1])  ; A
                        vlos0 = (wvl_center_pix_fit - wavelength.wvl_h1e)/wavelength.wvl_h1e*cc
                        vlos[ifile,ix,idir] = vlos0

;wdef,0,600,600
;!p.multi=0
;plot,fit_wvl_pix,y,psym=1
;oplot,fit_wvl_pix,yfit
;stop
		endfor
	endfor
endfor

save,vlos,file=save_file

END
;=========================================================================
;=========================================================================
PRO CHECK_CROSSTALK_02, cal01_coe = cal01_coe, image_dir=image_dir
COMMON ANALYSIS

if is_dir(image_dir) eq 0 then spawn,'mkdir '+image_dir

restore,cal01_coe
;save,coe0_arr,coe1_arr,file=save_file

tmp = size(coe0_arr) & nx=tmp[1] & ny=tmp[2] & ndir=tmp[3] & npol=tmp[4] & nline=tmp[5]

margin=10

wdef,0,800,800
!p.multi=[0,3,3]
for icoe=0,1 do begin
	case icoe of
		0:arr=coe0_arr
		1:arr=coe1_arr
	endcase
	for idir=0,ndir-1 do begin
		outfile=image_dir+'coe'+string(icoe,format='(i1.1)')+'_'+string(idir,format='(i2.2)')+'.png'
		for iline=0,nline-1 do begin
			case iline of
				0:begin
					line='397 nm'
					hair=hairline.pos397
				end
				1:begin
					line='630 nm'
					hair=hairline.pos630
				end
				2:begin
					line='854 nm'
					hair=hairline.pos854
				end
			endcase
			for ipol=0,npol-1 do begin
				case ipol of
					0:pol='Q/I'
					1:pol='U/I'
					2:pol='V/I'
				endcase

				mm = mean(arr[hair[0]+margin:hair[1]-margin,*,idir,ipol,iline])
				ss = stddev(arr[hair[0]+margin:hair[1]-margin,*,idir,ipol,iline])

				plot_image,arr[*,*,idir,ipol,iline],/nosquare,min=mm-3.*ss,max=mm+3.*ss, $
					charsize=2.,color=0,background=255,	$
					ticklen=-0.02,xtitle='Slit (pix)',ytitle='Scan (steps)',  $
					title=line+', '+pol



			endfor
		endfor

		write_png,outfile,tvrd()
	endfor
endfor


END
;=========================================================================
PRO SYNTH_397_WFA_01, wfa = wfa, atm_disp_file = atm_disp_file
COMMON ANALYSIS

restore,wfa
;save,bv630,bv854,file=save_file

restore,atm_disp_file,/verb & atm_disp = res
;yarc = (xpix + atm_disp - index.naxis2*0.5)*index.cdelt2 + index.ycen + 0.2

dlambda = wavelength.wvl397[1] - wavelength.wvl397[0]

; Parameters for H epsilon
f = 1.     ; filling factor
g = 1.3333 ; effective Lande factor, ref) zline.pro
gg= 1.3333 ; the Landé factor for the transverse magnetic field, ref) zline.pro

ndir=n_elements(data_dir.dir397)
for idir=0,ndir-1 do begin
	dir = data_dir.dir397[idir]+'_cal01/'
	filei=file_search(dir,'*_I_*.fits')

	for ifile=0,n_elements(filei)-1 do begin
		mreadfits,filei[ifile],index,data
		tmp=size(data) & nx=tmp[1] & ny=tmp[2]

		; dI/dlambda
		d1i = data>0.<0.
		d1i[*,1:ny-2] = (data[*,2:*] - data[*,0:ny-3])/2./dlambda

		; d2I/dlambda2
		d2i = data>0.<0.
		d2i[*,1:ny-2] = (d1i[*,2:*] - d1i[*,0:ny-3])/2./dlambda

		for iline=0,1 do begin
			case iline of
				0:begin
					;== WFA 630 nm
					;variables
					bv = bv630
					disp = atm_disp[1,idir]
					index_ref = index_for_map.index_nosun_630
					outdir = data_dir.dir397[idir]+'_cal01_wfa630/'
				end
				1:begin
					;== WFA 854 nm
					;variables
					bv = bv854
					disp = atm_disp[2,idir]
					index_ref = index_for_map.index_nosun_854
					outdir = data_dir.dir397[idir]+'_cal01_wfa854/'
				end
			endcase
			if is_dir(outdir) eq 0 then spawn,'mkdir '+outdir

			;synthesize
			tmp=size(bv.blos) & nybv=tmp[2]
			datal = data>0.<0.
			dataq = data>0.<0.
			datau = data>0.<0.
			datav = data>0.<0.
			for iybv=0,nybv-1 do begin
				bv_arcsec = pix2arcsec_nosun(iybv, ifile, index_ref, atm_disp = disp)
				pix_397   = arcsec2pix_nosun(bv_arcsec[0], bv_arcsec[1], index_for_map.index_nosun_397, atm_disp = 0.)
				ix = round(pix_397[0])

				if ix ge 0 and ix lt nx then begin
					blos = bv.blos[ifile,iybv,idir]
					bpos = bv.bpos[ifile,iybv,idir]
					cph  = bv.cph[ifile,iybv,idir]
			
					dlamb = 4.67e-13 * wavelength.wvl_h1e^2 * blos
					datav[ix,*] = -dlamb*f*g*d1i[ix,*]

					dlamb = 4.67e-13 * wavelength.wvl_h1e^2 * bpos
					;datal[ix,*] = 3./4.* dlamb^2 *f * gg * d1i[ix,*] / (wavelength.wvl397 - wavelength.wvl_h1e)
					datal[ix,*] = -1./4.* dlamb^2 *f * gg * d2i[ix,*]
					dataq[ix,*] = datal[ix,*]*cos(cph*2.)
					datau[ix,*] = datal[ix,*]*sin(cph*2.)
				endif
			endfor

                	; save data
                	outfile0=outdir+strmid(filei[ifile],strlen(dir),strlen(filei[ifile])-strlen(dir))
			tmp = strsplit(outfile0,"_I_",/extract,/regex) & outfile=tmp[0]+'_I_'+tmp[1]
			header=ta_make_hdr(index, /remove_z)
                	writefits,outfile,data,header

			tmp = strsplit(outfile0,"_I_",/extract,/regex) & outfile=tmp[0]+'_Q_'+tmp[1]
                	header=ta_make_hdr(index, /remove_z)
                	writefits,outfile,dataq,header

			tmp = strsplit(outfile0,"_I_",/extract,/regex) & outfile=tmp[0]+'_U_'+tmp[1]
                	header=ta_make_hdr(index, /remove_z)
                	writefits,outfile,datau,header

			tmp = strsplit(outfile0,"_I_",/extract,/regex) & outfile=tmp[0]+'_V_'+tmp[1]
                	header=ta_make_hdr(index, /remove_z)
                	writefits,outfile,datav,header
		endfor
	endfor
endfor

END
;=========================================================================
;=========================================================================
          
PRO SPATIAL_RESOLUTION_01, intensity_map = intensity_map, savefile = savefile
COMMON ANALYSIS

margin = 10

restore, intensity_map
;arr = reform(con397[*,*,1])
arr = con397


ss = size(arr)
nx = ss[1]
ny = ss[2]
ny1= hairline.pos397[1]-margin - (hairline.pos397[0]+margin) +1
res = fltarr(nx, ny1)

for ix=0,nx-1 do begin
	;ta_fft,reform(con397[ix,hairline.pos397[0]+margin:hairline.pos397[1]-margin,1]),index_for_map.index_nosun_397.cdelt2,k,a,hwin=1
	ta_fft,reform(con397[ix,hairline.pos397[0]+margin:hairline.pos397[1]-margin]),index_for_map.index_nosun_397.cdelt2,k,a,hwin=1
	res[ix,*] = abs(a)^2
endfor
;plot,k,abs(a)^2,/ylog

plot_image,res

;stop

END
;=========================================================================

PRO MAKE_DATA_03, obs_dir = obs_dir, save_dir = save_dir
COMMON ANALYSIS

if is_dir(save_dir) eq 0 then spawn,'mkdir '+save_dir
obs_file = file_search(obs_dir,'*.sav', count=nf)

for i=0,nf-1 do begin
        save_dir1 = save_dir+strmid(obs_file[i],strlen(obs_dir),strlen(obs_file[i])-strlen(obs_dir)-4)+'/'
        if is_dir(save_dir1) eq 0 then spawn,'mkdir '+save_dir1

        restore,obs_file[i]
;save,   stki_630,stkq_630,stku_630,stkv_630,wave_630,   $
;        stki_397,stkq_397,stku_397,stkv_397,wave_397,   $
;        stki_854,stkq_854,stku_854,stkv_854,wave_854,   $
;        b_los_xyz,h1e_ncp,      $
;        ibomb,index,dd, $
;        file=save_file

        tmp=size(h1e_ncp) & nx=tmp[1] & ny=tmp[2]
        tmp = min(h1e_ncp, pos)
        posx= pos mod nx
        posy= pos / nx
        print,tmp,h1e_ncp[posx,posy]

        for ix=0,nx-1 do begin
                print,ix
                for iy=0,ny-1 do begin
                        outfile = save_dir1+'x'+string(ix,format='(i3.3)')+'_y'+string(iy,format='(i3.3)')+'.txt'
                        openw,lun,outfile,/get_lun
                        for iline = 0,2 do begin
                                case iline of
                                        0:begin
                                                wave=wave_397
                                                stki=stki_397
                                                stkq=stkq_397
                                                stku=stku_397
                                                stkv=stkv_397
                                        end
                                        1:begin
                                                wave=wave_630
                                                stki=stki_630
                                                stkq=stkq_630
                                                stku=stku_630
                                                stkv=stkv_630
                                        end
                                        2:begin
                                                wave=wave_854
                                                stki=stki_854
                                                stkq=stkq_854
                                                stku=stku_854
                                                stkv=stkv_854
                                        end
                                endcase
                                printf,lun,strcompress(string(n_elements(wave)),/remove_all)
                                for il=0,n_elements(wave)-1 do begin
                                        printf,lun,string(wave[il],format='(f9.4)') + ' ' + $
                                                string(stki[ix,iy,il],format='(f9.6)') + ' ' + $
                                                string(stkq[ix,iy,il],format='(f9.6)') + ' ' + $
                                                string(stku[ix,iy,il],format='(f9.6)') + ' ' + $
                                                string(stkv[ix,iy,il],format='(f9.6)') + ' '
                                endfor
                        endfor
                        free_lun,lun
                        if ix eq posx and iy eq posy then begin
                                spawn,'cp '+outfile+' '+save_dir1+'ncp_x'+string(ix,format='(i3.3)')+'_y'+string(iy,format='(i3.3)')+'.txt'
                        endif
                endfor
        endfor
ENDFOR ; i

END
 


;=========================================================================
;=========================================================================
;==== MAIN ====;


COMMON ANALYSIS

time1=systime(/sec)


set_parameters
saveidl_head = './idl/'
imgdir       = './imgs/' 
name_tail    = '_cal01_tmp'

;=============================================
;== check or modify the wavelength calibration
image_dir = imgdir+'wvl_cal/'
if is_dir(image_dir) eq 0 then spawn,'mkdir '+image_dir
;== Feb/2022
idir = 0
wvl_cal,data_dir=data_dir.dir397[idir],image_file=image_dir+'397_feb2022_0.png',  $
	ref_wvl=[3964.500d,3968.492d,3971.31d],wavelength0=wavelength.wvl397   ; Fe I, Ca II, Fe I
wvl_cal,data_dir=data_dir.dir630[idir],image_file=image_dir+'630_feb2022_0.png',  $
	ref_wvl=[6301.508d,6302.000d,6302.499d,6302.764d],wavelength0=wavelength.wvl630
wvl_cal,data_dir=data_dir.dir854[idir],image_file=image_dir+'854_feb2022_0.png',  $
	ref_wvl=[8536.163d,8538.020d,8542.144d],wavelength0=wavelength.wvl854

wavelength.wvl630 = wvl_cal_1(data_dir=data_dir.dir630[0],image_file=image_dir+'630_feb2022_1.png',  $
;	ref_wvl=[6302.000d,6302.764d],ref_pix=[525,580])   ; O2, Fe I lines
	ref_wvl=[6295.960d,6306.565d],ref_pix=[52,881])   ; O2
save,wavelength,file=saveidl_head+'wvl.sav'

;=============================================
;== Calibration for crosstalks, then rotate +Q so that +Q directs to the disc center, wants more?

cal01,image_dir=imgdir+'cal01/', result_check_crosstalk_01 = saveidl_head+'check_crosstalk_01.sav', save_file=saveidl_head+'cal01.sav'



; Make NCP of H epsilon data 
make_h1e_ncp_01,save_file=saveidl_head+'h1e_ncp_01.sav', $
	image_file=imgdir+'make_h1e_ncp_01.png'

; Make NCP of H epsilon data from dI/dlambda = 0, derive Doppler velocity as well from dI/dlambda=0
make_h1e_ncp_02,save_file=saveidl_head+'h1e_ncp_02.sav', $
	image_dir=imgdir+'make_h1e_ncp_02/'

; Gauss fit H epsilon to get VLOS
make_h1e_vlos_01, save_file=saveidl_head+'make_h1e_vlos_01.sav'

; Make Ca H maps
;make_average_profile,data=data_dir.dir397[0]+'_cal02',save_file=saveidl_head+'397profile.sav'
make_cah_01, save_file=saveidl_head+'cah_01.sav'
;make_cah_01, save_file=saveidl_head+'cah_01.sav'
make_con_01, save_file=saveidl_head+'con_01.sav'

make_h1e_i_01, save_file=saveidl_head+'h1e_i_01.sav'

make_h1e_iquv_01, save_file=saveidl_head+'h1e_iquv_01.sav'

make_fe_v_01, save_file=saveidl_head+'fe_v_01.sav'

; Weak Field Approximation, Doppler velocity at line center
wfa_01,save_file=saveidl_head+'wfa01.sav'


; Change B coordinate from LOS to solar normal
blos2bnormal_01, wfa=saveidl_head+'wfa01.sav', $
	save_file=saveidl_head+'wfa_solarnormal.sav',  $
	con=saveidl_head+'con_01.sav'

;-- Calculate how many pixel 630 nm or 854 nm images should be shifted in Y to match 397 nm
;-- res array [3 wavelengths, 2 data]
cal_earth_atm_dispersion, cont_file = saveidl_head+'con_01.sav', savefile=saveidl_head+'atm_disp.sav', img_dir=imgdir+'atm_dispersion'


synth_397_wfa_01, wfa = saveidl_head+'wfa01.sav', atm_disp_file = saveidl_head+'atm_disp.sav' 

spatial_resolution_01, intensity_map = saveidl_head+'con_01.sav', savefile = saveidl_head+'spatial_resolution_01.sav'


; Make data set 10"x10" FOV with 0.2" samplings
make_data_01, save_dir = saveidl_head+'make_data_01/', atm_disp_file = saveidl_head+'atm_disp.sav', wfa_los=saveidl_head+'wfa01.sav', $
	h1e_ncp_01 = saveidl_head+'h1e_ncp_02.sav', wfa_solarnormal = saveidl_head+'wfa_solarnormal.sav'


; Make text and fits data for DeSIRe, calibrate Stokes-Q in 854 nm, (not yet,,,calibrate tilt of CaIIH)
make_data_08, obs_dir = saveidl_head+'make_data_01/', save_dir = saveidl_head+'make_data_08_txt/', format='txt', no_hepsilon=0
make_data_08, obs_dir = saveidl_head+'make_data_01/', save_dir = saveidl_head+'make_data_08_fts/', format='fits', no_hepsilon=0

make_ncp_map_01, image_dir = imgdir+'make_ncp_map_01/', $
	save_file = saveidl_head+'make_ncp_map_01.sav'


;== Make quiet profiles ==;

; Make data set 10"x10" FOV with 0.2" samplings
;make_data_01, save_dir = saveidl_head+'make_data_01_ref/', atm_disp_file = saveidl_head+'atm_disp.sav', wfa=saveidl_head+'wfa01.sav', $
;	h1e_ncp_01 = saveidl_head+'h1e_ncp_01.sav', /ref
make_data_01, save_dir = saveidl_head+'make_data_01_ref/', atm_disp_file = saveidl_head+'atm_disp.sav', wfa_los=saveidl_head+'wfa01.sav', $
        h1e_ncp_01 = saveidl_head+'h1e_ncp_02.sav', wfa_solarnormal = saveidl_head+'wfa_solarnormal.sav',/ref

;skip:;-----------------------------------
;restore,saveidl_head+'wvl.sav' ; ANAN


; Make text data for DeSIRe
make_data_03, obs_dir = saveidl_head+'make_data_01_ref/', save_dir = saveidl_head+'make_data_03_ref/'

print,'It takes ',systime(/sec) - time1

END
