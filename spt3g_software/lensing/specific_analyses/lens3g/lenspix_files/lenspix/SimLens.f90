    !Simple program demonstrating how to generate a simulated lensed map
    !AL, Feb 2004; Updated Oct 2007
    program SimLensCMB
    use HealpixObj
    use HealpixVis
    use Random
    use spinalm_tools
    use IniFile
    use AMLUtils
    implicit none
    Type(HealpixInfo)  :: H
    Type(HealpixMap)   :: M,MI, GradPhi,MM
    Type(HealpixPower) :: P
    Type(HealpixAlm)   :: A,AP

    integer            :: nside, lmax
    integer(I_NPIX)    :: npix
    character(LEN=1024)  :: w8name = '../Healpix_2.00/data/'
    character(LEN=1024)  :: file_stem, cls_file, out_file_root, cls_lensed_file, phi_file,TQU_file
    character(LEN=1024) :: healpixloc
    integer, parameter :: lens_interp =1, lens_exact = 2
    integer :: lens_method = lens_interp
    integer :: mpi_division_method = division_equalrows
    integer ::  interp_method,  rand_seed
    logical :: err, want_pol, set_B_zero, set_E_zero
    real :: interp_factor
    integer status
#ifdef MPIPIX
    integer i

    call mpi_init(i)
#endif

    Ini_Fail_On_Not_Found = .true.
    call Ini_Open(GetParam(1), 3,err)
    if (err) then
#ifdef MPIPIX
        call mpi_finalize(i)
#endif
        stop 'No ini'
    end if
    nside  = Ini_Read_Int('nside')
    npix = nside2npix(nside)

    lmax   = Ini_Read_Int('lmax')  
    cls_file = Ini_Read_String('cls_file')
    out_file_root = Ini_Read_String('out_file_root')
    phi_file  = Ini_Read_String('phi_file')
    TQU_file  = Ini_Read_String('TQU_file')

    lens_method = Ini_Read_Int('lens_method')
    want_pol = Ini_Read_Logical('want_pol')
    set_B_zero = Ini_Read_Logical('set_B_zero')
    set_E_zero = Ini_Read_Logical('set_E_zero')
    rand_seed = Ini_Read_Int('rand_seed')

    interp_method = Ini_read_int('interp_method')

    Ini_Fail_On_Not_Found = .false.


    w8name = Ini_Read_String('w8dir')
    interp_factor=0
    if (lens_method == lens_interp) interp_factor = Ini_Read_Real('interp_factor',3.)
#ifdef MPIPIX
    mpi_division_method = Ini_Read_Int('mpi_division_method',division_balanced);
#endif 

    call Ini_Close

    file_stem =  trim(out_file_root)//'_lmax'//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))// &
    '_interp'//trim(RealToStr(interp_factor,3))//'_method'//trim(IntToStr(interp_method))//'_'

    if (want_pol) file_stem=trim(file_stem)//'pol_'
    file_stem = trim(file_stem)//trim(IntToStr(lens_method)) 

    cls_lensed_file  = trim(file_stem)//'.dat'

    call SetIdlePriority()

    if (w8name=='') then
        call get_environment_variable('HEALPIX', healpixloc, status=status)
        if (status==0) then
            w8name = trim(healpixloc)//'/data/'
        end if
    end if

    if (w8name=='') then
        write (*,*) 'Warning: using unit weights as no w8dir found'
        call HealpixInit(H,nside, lmax,.true., w8dir='', method= mpi_division_method) 
    else
        call HealpixInit(H,nside, lmax,.true., w8dir=w8name,method=mpi_division_method) 
    end if 

    if (H%MpiID ==0) then !if we are main thread
        !All but main thread stay in HealpixInit

        call HealpixPower_nullify(P)
        call HealpixAlm_nullify(A)
        call HealpixAlm_nullify(AP)
        call HealpixMap_nullify(GradPhi)
        call HealpixMap_nullify(M)
        call HealpixMap_nullify(MI)
        call HealpixMap_nullify(MM)

        call HealpixPower_ReadFromTextFile(P,cls_file,lmax,pol=.true.,dolens = .true.)

        call HealpixAlm_Sim(A, P, rand_seed,HasPhi=.true., dopol = want_pol)
        call HealpixAlm2Power(A,P)
        call HealpixPower_Write(P,trim(file_stem)//'_unlensed_simulated.dat')
 
        call HealpixAlm2Map(H,A,M, npix, DoPhi=.true.)

        call HealpixMap_Read(MM,TQU_file)
        M%TQU=MM%TQU
        call HealpixMap_Free(MM)

        call HealpixMap_Read(MI,phi_file,phi_map=.true.)
        M%Phi=MI%Phi
        call HealpixMap_Free(MI)

        call HealpixMap2Alm(H,M,A,lmax)
        if (set_B_zero) then
            A%TEB(3:3,:,:)=0
            call HealpixAlm2Power(A,P)
            call HealpixPower_Write(P,trim(file_stem)//'_unlensed_simulated_zeroB.dat')
        end if

        if (set_E_zero) then
            A%TEB(2:2,:,:)=0
            call HealpixAlm2Power(A,P)
            call HealpixPower_Write(P,trim(file_stem)//'_unlensed_simulated_zeroE.dat')
        end if

        call HealpixMap_Write(M,trim(file_stem)//'_phi.fits',overwrite=.true.,phi_map=.true.)

        call HealpixAlm2GradientMap(H,A, GradPhi,npix,'PHI')

        if (lens_method == lens_exact) then
            call HealpixExactLensedMap_GradPhi(H,A,GradPhi,M)
        else if (lens_method == lens_interp) then
            call HealpixInterpLensedMap_GradPhi(H,A,GradPhi, M, interp_factor, interp_method)
        else
            stop 'unknown lens_method'
        end if

        call HealpixMap2Alm(H,M, A, lmax, dopol = want_pol)

        call HealpixAlm2Power(A,P)
        call HealpixAlm_Free(A)

        call HealpixPower_Write(P,cls_lensed_file)

        call HealpixMap_Write(M,trim(file_stem)//'_lensed_map.fits',overwrite=.true.)

    end if

#ifdef MPIPIX
    call HealpixFree(H)
    call mpi_finalize(i)
#endif

#ifdef DEBUG
    write (*,*) 'End of program'
    pause
#endif
    end program SimLensCMB
