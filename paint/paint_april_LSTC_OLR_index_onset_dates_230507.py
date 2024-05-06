'''
2023-05-07
This code plot the onset dates with LSTC and OLR index to compare the influence of the two elements on it
'''
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import numpy as np
import xarray as xr
import sys
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


sys.path.append("/home/sun/mycode/module/")
from module_sun import *

def open_onsetdate(file):
    with open(file,'r') as load_f:
        a = json.load(load_f)

    year = np.array(list(a.keys()))    ;  year  =  year.astype(int)
    day  = np.array(list(a.values()))  ;  day   =  day.astype(int)

    return year,day

def select_anomaly_year(day):
    # 筛选 晚年用红色，早年用蓝色
    a  =  np.zeros(40,dtype=int).astype(dtype=str) ; a[:]  =  'grey'
    color_list  =  a.tolist()#  ;  color_list[:]  =  'grey'

    for i in range(0,40):
        if day[i] < np.mean(day) - np.std(day):
            color_list[i]  =  'blue'
        if day[i] > np.mean(day) + np.std(day)-1:
            color_list[i]  =  'red'

    return color_list

def set_pic_ticks(
    ax,xticks,yticks,
    xlabels,ylabels,
    x_minorspace=None,y_minorspace=None,
    labelsize=12,axis='both',
    xaxis_label=None,yaxis_label=None,axis_labelsize=20
    ):
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
    # 设置tick
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    # 设置tick_label
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(ylabels)
    # 设置最小刻度 默认为1
    ax.xaxis.set_minor_locator(MultipleLocator(x_minorspace))
    ax.yaxis.set_minor_locator(MultipleLocator(y_minorspace))

    # 设置labelsize大小
    ax.tick_params(axis=axis,labelsize=labelsize)

    # 设置坐标轴标签
    if xaxis_label != None:
        ax.set_xlabel(xaxis_label,fontsize=axis_labelsize)
    if yaxis_label != None:
        ax.set_ylabel(yaxis_label,fontsize=axis_labelsize)

def plot_baseline(axs,day):
    axs.plot([1979,2020],[0,0],color='black')
    axs.plot([1979,2020],[np.ceil(np.mean(day)-np.std(day))-120-1,np.ceil(np.mean(day)-np.std(day))-120-1],color='k',linestyle='dashed')
    axs.plot([1979,2020],[np.floor(np.mean(day)+np.std(day))-120-1,np.floor(np.mean(day)+np.std(day))-120-1],color='k',linestyle='dashed')

def check_path(path):
    # 检查生成的路径是否存在，如果不存在则创建之
    from pathlib import Path
    import os
    path1 = Path(path)
    if path1.exists():
        print("generation path is exists")
    else:
        os.mkdir(path)

def add_legend():
    # 添加legend
    from matplotlib.patches import Patch
    import matplotlib.pyplot as plt

    legend_elements = [
                        Patch(color='grey',label='normal'),
                        Patch(color='blue',label='early'),
                        Patch(color='red',label='late')]
    plt.legend(handles=legend_elements, loc=3,edgecolor='white',facecolor='white',bbox_to_anchor=(0.2, 0.05),prop={'size': 20})                              
    # Patch这个函数我不是很明白，只调用这个函数不会生成什么
    # 得配合legend一起用，可能就是用来生成legend样式的吧？

def main():
    # Read onset date data
    year,day = open_onsetdate("/home/sun/qomo-data/onsetdate.json")
    print(year)

    # Read OLR index
    olr_file = xr.open_dataset("/home/sun/data/ERA5_data_monsoon_onset/index/OLR_maritime_area_average_monthly.nc").sel(year=slice(1980, 2019))

    # Read LSTC index
    lstc_file = xr.open_dataset("/home/sun/data/ERA5_data_monsoon_onset/index/LSTC_index_April_ERA5_1959to2021.nc") # position is 21

    olr      = olr_file['ttr_avg_mon'].data[:, 3]
    print(olr.shape)
    lstc     = lstc_file['LSTC_april'].data[21:-2]
    print(lstc.shape)

    # standardization
    olr      = olr - np.average(olr)
    lstc     = lstc - np.average(lstc)

    colorlist = select_anomaly_year(day)

    fig,axs  =  plt.subplots(2, figsize = (16, 10))

    y_label  =  ['10Apr','20Apr','1May','10May','20May']
    x_label  =  np.arange(1980,2021,10)

    ytick    =  np.arange(-20,25,10)
    xtick    =  np.arange(1980,2021,10)

    # === First Picture ===

    set_pic_ticks(axs[0],xtick,ytick,x_label,y_label,x_minorspace=5,y_minorspace=5,xaxis_label="Year",yaxis_label="BoBSM Onset Date")
    axs[0].set_xlim(1979,2020)

    axs[0].bar(year,day-120,width=1,color=colorlist,edgecolor='black')
    plot_baseline(axs[0],day)
    #add_legend()

    # Add OLR to first pic
    ax2 = axs[0].twinx()

    ax2.set_ylabel("Maritime Continent OLR", fontsize=20)
    ax2.plot(year, olr, color='black', marker='o', alpha=0.8)

    # === Second Picture ===
    set_pic_ticks(axs[1],xtick,ytick,x_label,y_label,x_minorspace=5,y_minorspace=5,xaxis_label="Year",yaxis_label="BoBSM Onset Date")
    axs[1].set_xlim(1979,2020)

    axs[1].bar(year,day-120,width=1,color=colorlist,edgecolor='black')
    plot_baseline(axs[1],day)
    #add_legend()

    # Add OLR to first pic
    ax2 = axs[1].twinx()

    ax2.set_ylabel("LSTC index", fontsize=20)
    ax2.plot(year, lstc, color='green', marker='*', alpha=0.8)

    path_out = "/home/sun/paint/index_correlation_with_onsetdate/" ; check_path(path_out)
    file_out = "onset_dates_BOB_LSTC_OLR.pdf"
    plt.show()
    plt.savefig(path_out+file_out,dpi=450)

if __name__ == "__main__":
    main()