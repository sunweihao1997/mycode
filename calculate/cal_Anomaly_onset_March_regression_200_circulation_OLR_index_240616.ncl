;2023-5-5
;This script calculate the regression of wind to the OLR
begin
    ; ================ Data processing1: abnormal, detrend, calculate monthly mean =======================
    ; Wind data
    fvar = addfile("/home/sun/data/monsoon_onset_anomaly_analysis/ERA5_data_monsoon_onset/monthly/ERA5_1980_2021_monthly_200hpa_UVZ_500hpa_w.nc", "r")


    ; Land-Sea thermal contrast
    f1 = addfile("/home/sun/data/monsoon_onset_anomaly_analysis/process/ERA5_LSTC_1980_2021_t2m_ts_sp_psl_month/ERA5_LSTC_1980_2021_t2m_ts_sp_psl_mar.nc", "r")

    ; Read data
    u  = fvar->Mar_u
    v  = fvar->Mar_v
    psl= fvar->Mar_z
    w  = fvar->Mar_w

    lstc = f1->OLR_mari_Afri


;    print(lstc)

    ; Calculates and removes the mean of the given dimension(s) at all other dimensions.
    lstc_abnormal = dim_rmvmean_n_Wrap(lstc, 0)
    u_abnormal    = dim_rmvmean_n_Wrap(u, 0)
    v_abnormal    = dim_rmvmean_n_Wrap(v, 0)
    w_abnormal    = dim_rmvmean_n_Wrap(w, 0)
    psl_abnormal  = dim_rmvmean_n_Wrap(psl/9.8,0)
    

    ;print(u_abnormal)
    ; Detrend
    lstc_detrend  = lstc_abnormal
    lstc_detrend  = (/dim_standardize_n_Wrap(dtrend_n(lstc_abnormal, True, 0), 1, 0)/)

    u_detrend     = u_abnormal
    u_detrend     = (/dtrend_n(u_abnormal, True, 0)/)
    v_detrend     = v_abnormal
    v_detrend     = (/dtrend_n(v_abnormal, True, 0)/)
    w_detrend     = w_abnormal
    w_detrend     = (/dtrend_n(w_abnormal, True, 0)/)
    p_detrend     = psl_abnormal
    p_detrend     = (/dtrend_n(psl_abnormal, True, 0)/)


    ; Change the dimension sequence
    ;printVarSummary(u_abnormal)
    u_detrend_    = u_detrend
    v_detrend_    = v_detrend
    w_detrend_    = w_detrend
    p_detrend_    = p_detrend
    ;printVarSummary(lstc_detrend)
    ;printVarSummary(u_detrend_)

    ; Delete NAN value
    if (any(isnan_ieee(u_detrend_))) then
        value = 1.e20
        replace_ieeenan(u_detrend_, value, 0)
        replace_ieeenan(v_detrend_, value, 0)
        replace_ieeenan(w_detrend_, value, 0)
        replace_ieeenan(p_detrend_, value, 0)
        replace_ieeenan(lstc_detrend, value, 0)
        u_detrend_@_FillValue = value
        v_detrend_@_FillValue = value
        p_detrend_@_FillValue = value
        w_detrend_@_FillValue = value
        lstc_detrend@_FillValue = value

    end if


    ;=================================================
    ; ================ Data processing1: abnormal, detrend, calculate monthly mean =======================

    ; ================ Data processing2: calculate regression =======================
    ;printVarSummary(u_detrend)
    ;printVarSummary(lstc_detrend_april)
    rc_u = regCoef_n(lstc_detrend, u_detrend_, 0, 0)  
    rc_v = regCoef_n(lstc_detrend, v_detrend_, 0, 0)  
    rc_p = regCoef_n(lstc_detrend, p_detrend_, 0, 0)
    rc_w = regCoef_n(lstc_detrend, w_detrend_, 0, 0)

    ; ================ Write to the nc file ==========================================
    copy_VarCoords(u(0, :, :), rc_u)
    copy_VarCoords(v(0, :, :), rc_v)
    copy_VarCoords(w(0, :, :), rc_w)
    copy_VarCoords(psl(0, :, :), rc_p)

    delete(rc_u@nptxy)
    delete(rc_u@rstd)
    delete(rc_u@yintercept)
    delete(rc_u@tval)
    delete(rc_v@nptxy)
    delete(rc_v@rstd)
    delete(rc_v@yintercept)
    delete(rc_v@tval)
    delete(rc_p@nptxy)
    delete(rc_p@rstd)
    delete(rc_p@yintercept)
    delete(rc_p@tval)
    delete(rc_w@nptxy)
    delete(rc_w@rstd)
    delete(rc_w@yintercept)
    delete(rc_w@tval)

;    printMinMax(lstc_detrend, 0)
;    printMinMax(rc_u, 0)
    system("rm -rf /home/sun/data/monsoon_onset_anomaly_analysis/ERA5_data_monsoon_onset/regression/ERA5_regression_200_uvz_500_w_to_OLR.nc")
    fout = addfile("/home/sun/data/monsoon_onset_anomaly_analysis/ERA5_data_monsoon_onset/regression/ERA5_regression_200_uvz_500_w_to_OLR.nc", "c")
    filedimdef(fout,"time",-1,True)
;
    fout->rc_u = rc_u
    fout->rc_v = rc_v
    fout->rc_p = rc_p
    fout->rc_w = rc_w
;
    ;; Set file attribute
    setfileoption(fout,"DefineMode",True)
    fAtt = True
    fAtt@description = "Create on 17/6/2024. This file include regression from u/v/z at 200 hPa and w at 500 hPa to the OLR index. Created by cal_Anomaly_onset_March_regression_200_circulation_OLR_index_240616.ncl"
    fileattdef(fout,fAtt)

end
