@./subpro/ta_make_hdr.pro



PRO CAL01, data_dir, observation, cal01_param,    image_file=image_file, result_check_crosstalk_01=result_check_crosstalk_01
;COMMON ANALYSIS

; offset_ifile : offset in ifile to make artpol (artification polarization spectra)

plot_cah_special = 0

outname_tail='_cal01'

wdef,0,900,900
!p.multi=[0,3,3]
loadct,0


ndir=n_elements(data_dir.dir397)

for idir=0,ndir-1 do begin
        angle_q = atan(observation[idir].hxy[1], observation[idir].hxy[0])
        for iarm=0,2 do begin
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

                outdir=dir+outname_tail ; '_cal01'
                if is_dir(outdir) eq 0 then spawn,'mkdir '+outdir

                filei=file_search(dir,'*_I_*.fits')
                fileq=file_search(dir,'*_Q_*.fits')
                fileu=file_search(dir,'*_U_*.fits')
                filev=file_search(dir,'*_V_*.fits')
                for ifile=0,n_elements(filei)-1 do begin
                        read_sdo,filei[ifile],indexi,datai0, /use_shared_lib
                        read_sdo,fileq[ifile],indexq,dataq0, /use_shared_lib
                        read_sdo,fileu[ifile],indexu,datau0, /use_shared_lib
                        read_sdo,filev[ifile],indexv,datav0, /use_shared_lib
                        tmp=size(datai0) & nx=tmp[1] & ny=tmp[2]

                        if idir eq 0 and iarm eq 0 and ifile eq 0 then begin
                                coe0_arr = fltarr(2600,n_elements(filei),ndir,3,3)  ; slit, scan, data, quv, lines
                                coe1_arr = fltarr(2600,n_elements(filei),ndir,3,3)  ; slit, scan, data, quv, lines
                        endif


                        data1 = datai0>0.<0.
                        for ipol=1,3 do begin
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
                                                plot_image,data0,/nosquare,title='Stokes '+str_stks,    $
                                                        xtitle='X (pixel)',ytitle='Y (pixel)',          $
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
                                        ;if ipol eq 3 then write_png,outfile0,tvrd(/true)
                                endif
                        endfor

                        ; rotate +Q so that the +Q directs to the disk center
                        dataq2 = cos(2.*angle_q)*dataq1 + sin(2.*angle_q)*datau1
                        datau2 = -sin(2.*angle_q)*dataq1 + cos(2.*angle_q)*datau1

                        ; save data
                        outfile=outdir+'/'+strmid(filei[ifile],strlen(dir)-1,strlen(filei[ifile])-strlen(dir)+1)
                        header=ta_make_hdr(indexi, /remove_z)
                        writefits,outfile,datai0,header

                        outfile=outdir+'/'+strmid(fileq[ifile],strlen(dir)-1,strlen(fileq[ifile])-strlen(dir)+1)
                        header=ta_make_hdr(indexq, /remove_z)
                        writefits,outfile,dataq2,header

                        outfile=outdir+'/'+strmid(fileu[ifile],strlen(dir)-1,strlen(fileu[ifile])-strlen(dir)+1)
                        header=ta_make_hdr(indexu, /remove_z)
                        writefits,outfile,datau2,header

                        outfile=outdir+'/'+strmid(filev[ifile],strlen(dir)-1,strlen(filev[ifile])-strlen(dir)+1)
                        header=ta_make_hdr(indexv, /remove_z)
                        writefits,outfile,datav1,header
                endfor;ifile
        endfor
endfor

!p.multi=0
loadct,0

END
;-------------------------------------------------------------------------------------------------
; MAIN
;-----


time1 = systime(/sec)

;===============================================
;User defined parameters
data_dir={      dir397: ['./demodata/BOLDQ'], $
                dir630: ['./demodata/BPKJO'], $
                dir854: ['./demodata/BQKNO']}

observation = {        date: '2022-02-23',     $
                        hxy : [-720.,310.]      $
                }

con397 = [75+indgen(80-75), 107+indgen(123-107),275+indgen(291-275),358+indgen(401-358),  $
        688+indgen(7),830+indgen(865-830), 955+indgen(970-955)]
con630 = [10+indgen(20),120+indgen(10),160+indgen(10),260+indgen(30), $
       355+indgen(15),700+indgen(80),900+indgen(50)]
con854 = [0+indgen(190),260+indgen(300-260),360+indgen(20),430+indgen(20),700+indgen(830-700),940+indgen(980-940)]

cal01_param = { $
                con397: con397,         $
                con630: con630,         $
                con854: con854,         $
                err397: con397>1.<1.,   $
                err630: con630>1.<1.,   $
                err854: con854>1.<1.,   $
                polyorder397: 10,       $
                polyorder630: 5,        $
                polyorder854: 5,        $
                smooth_width:50,        $
                plot_idir: 1,   $
                plot_ifile: 111,        $
                plot_ix: 1020,          $
                smooth_offset: 1,       $
                smooth_factor:1e10      $
                }




;===============================================



CAL01, data_dir, observation, cal01_param

print,'It takes ', systime(/sec) - time1


END

