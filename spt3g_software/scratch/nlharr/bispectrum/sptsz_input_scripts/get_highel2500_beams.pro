;;;
; NAME: get_highel2500_beams.pro
; PURPOSE:
;   Return a 2d vector of the beam for a given field
;
; INPUTS: year [int]
;    l_beam,         empty, returned vector of l-values
;    bl_beam,        empty, returned vector of beam values
;
; OPTIONAL INPUTS:
;    beamfiles,      return the beam file if requested
;
; OUTPUTS: 
;    beamfiles,      The string path to the beam files used
;
; NOTES:
; 1) I got these files from (on spt): /home/lizinvt/highel2500/xspec/mk_spectra_run1_2011.pro
;
; MODIFICATION HISTORY:
;  08/01/2013: (KTS) Created, based on /home/lizinvt/highel2500/xspec/mk_spectra_run1_2011.pro
;;;

;...................................................................
; Return the beam
FUNCTION get_highel2500_beams, field, l_beam, bl_beam, sim=sim
beamdir = '/data/lizinvt/highel2500/beams/rev3.2/'
compbeamdir = '/data/lizinvt/highel2500/beams/composites_20130616_good/'
print, '********************'
case field of 
    'ra0h50dec-50': begin
        beamfiles=[beamdir+'blgrid_2010_90.txt', $
                   beamdir+'blgrid_2010_150.txt',$
                   beamdir+'blgrid_2010_220.txt']
        simbeamfiles=beamfiles
    end

    'ra1hdec-42.5': begin
        beamfiles=[compbeamdir+'blgridcomp_8x2010_79x2011_ra1hdec-42.5_90.txt', $
                   compbeamdir+'blgridcomp_8x2010_79x2011_ra1hdec-42.5_150.txt',$
                   compbeamdir+'blgridcomp_8x2010_79x2011_ra1hdec-42.5_220.txt']
        simbeamfiles=beamfiles
    end

    'ra1hdec-60': begin
        beamfiles=[beamdir+'blgrid_2010_90.txt', $
                   beamdir+'blgrid_2010_150.txt',$
                   beamdir+'blgrid_2010_220.txt']
        simbeamfiles=beamfiles
    end

    'ra21hdec-42.5': begin
        beamfiles=[compbeamdir+'blgridcomp_5x2010_83x2011_ra21hdec-42.5_90.txt', $
                   compbeamdir+'blgridcomp_5x2010_83x2011_ra21hdec-42.5_150.txt',$
                   compbeamdir+'blgridcomp_5x2010_83x2011_ra21hdec-42.5_220.txt']
        simbeamfiles=beamfiles
    end

    'ra21hdec-50': begin
        ;ellkmask=1 this has az and el data now!
        beamfiles=[beamdir+'blgrid_2009_90.txt',$
                   beamdir+'blgrid_2009_150.txt',$
                   beamdir+'blgrid_2009_220.txt']
        simbeamfiles=beamfiles
    end

    'ra21hdec-60': begin
        beamfiles=[beamdir+'blgrid_2009_90.txt', $
                   beamdir+'blgrid_2009_150.txt',$
                   beamdir+'blgrid_2009_220.txt']
        simbeamfiles=beamfiles
    end
    'ra22h30dec-55': begin
        beamfiles=[compbeamdir+'blgridcomp_4x2010_72x2011_ra22h30dec-55_90.txt', $
                   compbeamdir+'blgridcomp_4x2010_72x2011_ra22h30dec-55_150.txt',$
                   compbeamdir+'blgridcomp_4x2010_72x2011_ra22h30dec-55_220.txt']
        simbeamfiles=beamfiles
    end
    'ra23h30dec-55': begin
        beamfiles=[beamdir+'blgrid_2010_90.txt', $
                   compbeamdir+'blgridcomp_weighted_2008_2010_ra23h30dec-55_150.txt',$
                   compbeamdir+'blgridcomp_weighted_2008_2010_ra23h30dec-55_220.txt']
        simbeamfiles=beamfiles
    end
    'ra23hdec-45': begin
        beamfiles=[compbeamdir+'blgridcomp_4x2010_80x2011_ra23hdec-45_90.txt', $
                   compbeamdir+'blgridcomp_4x2010_80x2011_ra23hdec-45_150.txt',$
                   compbeamdir+'blgridcomp_4x2010_80x2011_ra23hdec-45_220.txt']
        simbeamfiles=beamfiles
    end
    'ra23hdec-62.5': begin
        beamfiles=[compbeamdir+'blgridcomp_4x2010_88x2011_ra23hdec-62.5_90.txt', $
                   compbeamdir+'blgridcomp_4x2010_88x2011_ra23hdec-62.5_150.txt',$
                   compbeamdir+'blgridcomp_4x2010_88x2011_ra23hdec-62.5_220.txt']
        simbeamfiles=beamfiles
    end
    'ra2h30dec-50': begin
        beamfiles=[beamdir+'blgrid_2010_90.txt', $
                   beamdir+'blgrid_2010_150.txt',$
                   beamdir+'blgrid_2010_220.txt']
        simbeamfiles=beamfiles
    end
    'ra3h30dec-42.5': begin
        beamfiles=[compbeamdir+'blgridcomp_3x2010_82x2011_ra3h30dec-42.5_90.txt', $
                   compbeamdir+'blgridcomp_3x2010_82x2011_ra3h30dec-42.5_150.txt',$
                   compbeamdir+'blgridcomp_3x2010_82x2011_ra3h30dec-42.5_220.txt']
        simbeamfiles=beamfiles
    end
    'ra3h30dec-60': begin
        beamfiles=[beamdir+'blgrid_2009_90.txt', $
                   beamdir+'blgrid_2009_150.txt',$
                   beamdir+'blgrid_2009_220.txt']
        simbeamfiles=beamfiles
    end
    'ra4h10dec-50': begin
        beamfiles=[beamdir+'blgrid_2010_90.txt', $
                   beamdir+'blgrid_2010_150.txt',$
                   beamdir+'blgrid_2010_220.txt']
        simbeamfiles=beamfiles
    end
    'ra5h30dec-45': begin
        beamfiles=[beamdir+'blgrid_2010_90.txt', $
                   beamdir+'blgrid_2010_150.txt',$
                   beamdir+'blgrid_2010_220.txt']
        simbeamfiles=beamfiles
    end
    'ra5h30dec-55': begin
        beamfiles = [beamdir+'blgrid_2011_90.txt', $
                     compbeamdir+'blgridcomp_weighted_2008_2011_ra5h30dec-55_150.txt', $
                     compbeamdir+'blgridcomp_weighted_2008_2011_ra5h30dec-55_220.txt']
        simbeamfiles=beamfiles
    end
    'ra6h30dec-45': begin
        beamfiles=[compbeamdir+'blgridcomp_8x2010_70x2011_ra6h30dec-45_90.txt', $
                   compbeamdir+'blgridcomp_8x2010_70x2011_ra6h30dec-45_150.txt',$
                   compbeamdir+'blgridcomp_8x2010_70x2011_ra6h30dec-45_220.txt']
        simbeamfiles=beamfiles
    end
    'ra6h30dec-55': begin
        beamfiles=[compbeamdir+'blgridcomp_29x2010_53x2011_ra6h30dec-55_90.txt', $
                   compbeamdir+'blgridcomp_29x2010_53x2011_ra6h30dec-55_150.txt',$
                   compbeamdir+'blgridcomp_29x2010_53x2011_ra6h30dec-55_220.txt']
        simbeamfiles=beamfiles
    end
    'ra6hdec-62.5': begin
        beamfiles=[compbeamdir+'blgridcomp_11x2010_74x2011_ra6hdec-62.5_90.txt', $
                   compbeamdir+'blgridcomp_11x2010_74x2011_ra6hdec-62.5_150.txt',$
                   compbeamdir+'blgridcomp_11x2010_74x2011_ra6hdec-62.5_220.txt']
        simbeamfiles=beamfiles
    end

    else: stop
endcase

bf = beamfiles
if keyword_set(sim) then bf = simbeamfiles

readcol,bf[0],l_tmp,bl_tmp,format='d,d'
nfreq=3
l_beam = l_tmp
bl_beam=dblarr(nfreq,n_elements(l_tmp))
for ii=0, nfreq-1 do begin
    ;print, "reading "+bf[ii]
    readcol,bf[ii],l_tmp,bl_tmp,format='d,d'
    bl_beam[ii,*] = bl_tmp
    ;print, n_elements(l_tmp), minmax(l_tmp)
endfor

RETURN, bf
END
