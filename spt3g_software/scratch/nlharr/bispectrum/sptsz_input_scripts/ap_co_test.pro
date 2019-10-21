

field = 'ra0h50dec-50'
band = 150
coadd_file = '/home/nlharr/Bispec_2500/coadds/ra0h50dec-50_150GHz_coadd.fits'
apod_file = '/home/nlharr/Bispec_2500/apod_masks/no_padding_square_def_params/ra0h50dec-50_test.sav'

psd_file = '/home/nlharr/Bispec_2500/noise_maps/ra0h50dec-50_150GHz_psd_file.sav'
jack_file = '/home/nlharr/Bispec_2500/noise_maps/ra0h50dec-50_150GHz_jack_file'


do_coadd = 0
do_apod  = 0
do_noise = 1


clean_files = 0

if (do_coadd) then begin
    if (clean_files) then begin
        file_delete, coadd_file
    endif
    coadd_one_field, field, band, coadd_file
endif

if (do_apod) then begin
    if (clean_files) then begin
        file_delete, apod_file
    endif
    dummy = apod_one_field(field, coadd_file, out_file=apod_file)
endif


if (do_noise) then begin
    if (clean_files) then begin
        file_delete, jack_file
        file_delete, psd_file
    endif
    noise_one_field, field, band, apod_file, psd_file, jack_file
endif




end
