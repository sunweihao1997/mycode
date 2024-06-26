'''
2022-10-7
本代码绘制论文version4.0中的fig5
内容为沿着10N的环流场
出版标准

吴老师意见: 伸展到100hpa
'''
import sys
import xarray as xr
import numpy as np
module_path = ["/home/sun/mycode/module/","/data5/2019swh/mycode/module/"]
sys.path.append(module_path[0])
from module_sun import *

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter

import warnings
warnings.filterwarnings("ignore")

from paint_lunwen_version3_0_fig2b_2m_tem_wind_20220426 import save_fig
from paint_lunwen_version3_0_fig2a_tem_gradient_20220426 import add_text

class number_date:
    '''定义一个日期序号类'''
    dates  =  [-6,-4,-2,0]
    date   =  [24,26,28,30]
    number =  ["a","b","c","d"]

def set_pic_ticks(
    ax,xticks,yticks,
    xlabels,ylabels,nx=1,ny=1,
    labelsize=20,axis='both',
    ):
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
    # 设置tick
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    # 设置tick_label
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(ylabels)

    # 设置次刻度.
    xlocator = mticker.AutoMinorLocator(nx + 1)
    ylocator = mticker.AutoMinorLocator(ny + 1)
    ax.xaxis.set_minor_locator(xlocator)
    ax.yaxis.set_minor_locator(ylocator)

    # 设置x轴 Formatter.
    ax.xaxis.set_major_formatter(LongitudeFormatter())


    # 设置labelsize大小
    ax.tick_params(axis=axis,labelsize=labelsize)



def paint_pic(lon,level,vwind,uwind,omega,topo_lon,topo):
    # 创建画布
    fig  =  plt.figure(figsize=(20,19))
    spec =  fig.add_gridspec(nrows=2,ncols=2)
    

    j  =  0
    for col in range(2):
        for row in range(2):
            ax   =  fig.add_subplot(spec[row,col])
            ax.invert_yaxis()

            # 设置坐标轴属性
            set_pic_ticks(ax=ax,xticks=np.linspace(50,110,4,dtype=int),yticks=np.linspace(1000,200,5,dtype=int),xlabels=np.linspace(50,110,4,dtype=int),ylabels=np.linspace(1000,200,5,dtype=int),labelsize=25)

            # 等值线
            im1   =  ax.contour(lon,level,vwind[number_date.date[j],:],[0],linewidths=2,colors='k')
            im2   =  ax.contourf(lon,level,vwind[number_date.date[j],:],np.linspace(-6,6,13),cmap=create_ncl_colormap("/home/sun/data/color_rgb/GMT_polar.txt",20),extend='both',alpha=1.0,zorder=0)
            ax.clabel(im1, [0], inline=True, fontsize=12)

            # 矢量图
            q  =  ax.quiver(lon[::4], level[::2], uwind[number_date.date[j],::2,::4], omega[number_date.date[j],::2,::4], 
                            angles='uv',# regrid_shape这个参数越小，是两门就越稀疏
                            scale_units='xy', scale=1.1,        # scale是参考矢量，所以取得越大画出来的箭头就越短
                            units='xy', width=2.2,
                            color='k',zorder=2)

            # 添加地形
            ax2  =  ax.twinx()
            ax2.set_ylim((0,4.5))
            ax2.plot(topo_lon,topo/1000,color='k')
            ax2.fill_between(topo_lon.data,0,topo/1000,where=topo>0,color='k')

            ax2.set_yticklabels([])
            ax2.set_yticks([])



            # 加日期
            if number_date.dates[j]<0:
                add_text(ax=ax,string="D"+str(number_date.dates[j]),location=(0.87,0.92),fontsize=30)
            elif number_date.dates[j]>0:
                 add_text(ax=ax,string="D+"+str(number_date.dates[j]),location=(0.87,0.92),fontsize=30)
            else:
                 add_text(ax=ax,string="D"+str(number_date.dates[j]),location=(0.87,0.92),fontsize=30)
            # 加序号
            add_text(ax=ax,string="("+number_date.number[j]+")",location=(0.015,0.92),fontsize=30)

            ax.set_title("D"+str(number_date.dates[j]),loc='right',fontsize=25)
            ax.set_title("("+number_date.number[j]+")",loc='right',fontsize=25)

            j+=1



    # 加colorbar
    fig.subplots_adjust(top=0.8) 
    cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03]) 
    cb  =  fig.colorbar(im2, cax=cbar_ax, shrink=0.1, pad=0.01, orientation='horizontal')
    cb.ax.tick_params(labelsize=20)
    save_fig(path_out="/home/sun/paint/lunwen/version3.0/",file_out="lunwen_fig5_v3.0_along10N_field.pdf")

def main():
    # 读取范围
    lon_slice  =  slice(45,120)
    lat_slice  =  slice(10,15)
    lev_slice  =  slice(1000,100)

    # 读取数据
    path  =  "/home/sun/qomo-data/"
    f1  =  xr.open_dataset(path+"composite3.nc").sel(lat=lat_slice,lon=lon_slice,level=lev_slice)
    f3  =  xr.open_dataset("/home/sun/data/gebco/bathymetric.nc").sel(lat=lat_slice,lon=lon_slice)

    # 计算经向平均
    uwind  =  np.nanmean(f1.uwind,axis=2)
    vwind  =  np.nanmean(f1.vwind,axis=2)
    omega  =  np.nanmean(f1.OMEGA,axis=2)*-60


    # 计算地形平均
    dixing  =  f3.elevation.data
    dixing[dixing <= 0]  =  0
    topo    =  np.average(dixing,axis=0)

    # 绘图
    paint_pic(lon=f1.lon.data,level=f1.level.data,vwind=vwind,uwind=uwind,omega=omega,topo_lon=f3.lon.data,topo=topo)

if __name__ == "__main__":
    main()
