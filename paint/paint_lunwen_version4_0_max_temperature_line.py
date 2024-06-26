'''
2022-10-4
This code plot the max temperature line
'''
import warnings
warnings.filterwarnings("ignore")  # close warning message
path  =  "/home/sun/data/composite/"


def cal_zonal_average(lon_slice,time_slice):
    '''
    Calculate zonal mean for input lon_slice
    '''
    # -------  import  --------------------- 
    import xarray as xr
    import numpy as np
    ##--------------------------------------

    #  ---------files-----------------------
    f0    =    xr.open_dataset(path + "composite3.nc").sel(lon=lon_slice).isel(time=time_slice)
    f1    =    xr.open_dataset(path + "composite-heating-merra.nc").sel(lon=lon_slice).isel(time=time_slice)

    # use variables: vwind omega T turbulence moist physics
    vwind     =    np.nanmean(f0.vwind,axis=3)
    omega     =    np.nanmean(f0.OMEGA,axis=3)

    tem_gradient       =    cal_gradient_meridional(f0.T.data)
    tem_gradient_avg   =    np.nanmean(tem_gradient,axis=3)

    latent    =    np.nanmean(f1.moist,axis=3)
    sensible  =    np.nanmean(f1.turbulence,axis=3)
    physics   =    np.nanmean(f1.physics,axis=3)

    return vwind,omega,tem_gradient_avg,latent,sensible,physics

def cal_gradient_meridional(vars):
    '''
    This function calculate meridional gradient. 
    support arrays: (time,level,lat,lon) and (time,lat,lon)

    process: 1. calculate zonal average  2. calculate meridional gradient
    '''
    # ---------- import -----------------------
    import sys
    sys.path.append("/home/sun/mycode/module/")
    from module_sun_new import cal_xydistance
    import xarray as xr
    import numpy as np
    from geopy.distance import distance
    ##-----------------------------------------

    # ----------------- get disx disy location ---------------------
    f0                   =    xr.open_dataset(path + "composite3.nc")
    disy,disx,location   =    cal_xydistance(f0.lat,f0.lon)

    # ------------------ calculate gradient ------------------------
    dimension  =  len(vars.shape)
    if dimension == 3:
        v_gradient  =  np.gradient(vars,location,axis=1)
    elif dimension == 4:
        v_gradient  =  np.gradient(vars,location,axis=2)

    return v_gradient

def generate_xlabel(array):
    '''This code generate labels for x axis'''
    labels = []
    for i in array:
        if i<0:
            labels.append(str(abs(i))+"S")
        elif i>0:
            labels.append(str(i)+"N")
        else:
            labels.append("EQ")
    return labels

def paint_meridional_circulation():
    '''This script paint meridional circulation and diabatic heating'''
    # -------------- import module ---------------------------------
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import xarray as xr
    import plotly.figure_factory as ff
    import cmasher as cmr
    sys.path.append("/home/sun/mycode_git/paint/")
    from paint_lunwen_version3_0_fig2a_tem_gradient_20220426 import add_text
    #import matplotlib
    #matplotlib.use('Agg')
    
    ## -------------------------------------------------------------

    # -------------- reference -------------------------------------
    f0        =  xr.open_dataset(path + "composite3.nc")
    f1        =  xr.open_dataset(path + "composite-heating-merra.nc")

    # -------------- variables -------------------------------------
    var_list  =  cal_zonal_average(lon_slice=slice(90,100),time_slice=[0,10,20,25,30])

    # circulation variables
    v         =  var_list[0]
    w         =  var_list[1] * -1

    # temperature variables
    tem_gradient  =  var_list[2]

    # diabatic heating variables
    latent        =  var_list[3]
    sensible      =  var_list[4]


    ## --------------  unify v and w  ------------------------------
    multiple  =  np.nanmean(abs(v))/np.nanmean(abs(w))
    w         =  w * 500
    print(multiple)

    ## --------------   vertical interpolate  ----------------------
    ## fuck python stream plot, it need reverse and interpolate vertical axis
    new_level =  np.linspace(1000,100,37)



    new_v     =  np.zeros((v.shape[0],37,v.shape[2]))
    new_w     =  new_v.copy()

    new_t     =  new_v.copy()

    new_s     =  np.zeros((latent.shape[0],37,latent.shape[2]))
    new_l     =  new_s.copy()

    for tt in range(v.shape[0]):
        for yy in range(v.shape[2]):
            old_level =  f0.level.data
            new_v[tt,:,yy]     =  np.interp(new_level[::-1],old_level[::-1],v[tt,::-1,yy])  #Here I did not invert after interpolate
            new_w[tt,:,yy]     =  np.interp(new_level[::-1],old_level[::-1],w[tt,::-1,yy])  #But the result paint is not error

            new_t[tt,:,yy]     =  np.interp(new_level[::-1],old_level[::-1],tem_gradient[tt,::-1,yy]) ; new_t[tt,:,yy]  =  new_t[tt,::-1,yy]

    for tt in range(sensible.shape[0]):
        for yy in range(sensible.shape[2]):
            old_level =  f1.level.data  # Because this variables belongs to two files, the vertical corrdinate is not the same
            new_s[tt,:,yy]     =  np.interp(new_level[::-1],old_level[::-1],sensible[tt,::-1,yy])     ; new_s[tt,:,yy]  =  new_s[tt,::-1,yy]
            new_l[tt,:,yy]     =  np.interp(new_level[::-1],old_level[::-1],latent[tt,::-1,yy])       ; new_l[tt,:,yy]  =  new_l[tt,::-1,yy]

    print("Successful interpolate variables")
    # ------------      date    ----------------------------------
    dates  =  ['-30','-20','-10','-5','0']


    # ------------     paint    ----------------------------------
    ## set figure
    fig1 = plt.figure(figsize=(12, 8))
    ax   = fig1.subplots()

    # set axis ticks and label
    ax.set_xticks(np.linspace(-10, 40, 6, dtype=int))
    ax.set_yticks(np.linspace(1000, 200, 5))

    ax.set_xticklabels(generate_xlabel(np.linspace(-10, 40, 6, dtype=int)))
    #ax.set_yticklabels(np.linspace(100,1000,10,dtype=int))
    ax.tick_params(axis='both', labelsize=22.5)

    # set axis limit
    ax.set_xlim((-10,40))

    j = 0
    im0 = ax.contour(f0.lat.data, new_level, new_t[j]*1e5,levels=[0],colors='gray',linewidths=4,linestyles='--',label='D-30') ; j += 1
    im1 = ax.contour(f0.lat.data, new_level, new_t[j]*1e5,levels=[0],colors='yellow',linewidths=4,linestyles='--',label='D-20') ; j += 1
    im2 = ax.contour(f0.lat.data, new_level, new_t[j]*1e5,levels=[0],colors='green',linewidths=4,linestyles='--',label='D-10') ; j += 1
    im3 = ax.contour(f0.lat.data, new_level, new_t[j]*1e5,levels=[0],colors='black',linewidths=4,linestyles='--',label='D-5') ; j += 1
    im4 = ax.contour(f0.lat.data, new_level, new_t[j]*1e5,levels=[0],colors='red',linewidths=4,linestyles='--',label='D0') ; j += 1

    ax.invert_yaxis()

    ax.legend()


    plt_path  =  "/home/sun/paint/lunwen/version4.0/"
    plt.savefig(plt_path+"composite_tem_max.pdf", dpi=400)

    plt.show()

    



def main():
    paint_meridional_circulation()

if __name__ == "__main__":
    main()

