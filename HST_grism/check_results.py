
import os
import numpy as np
#import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import glob
import astropy.io.fits as pyfits
from astropy.table import Table
from astropy.modeling import models, fitting
import pandas as pd

import tkinter as tk
from PIL import ImageTk, Image
import tkinter.font as fnt

from datetime import datetime
now = datetime.now()

import pyds9

import json

from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.widgets import Slider

import sys, getopt

line_list = np.array(['Lya', 'NV', 'NIV', 'NIII', 'HeII', 'OIII-1663', 'CIII', 'CIII', 'CIV', 'MgII', 'NeV', 'NeVI', 
                         'OII', 'NeIII', 
                          'Hd', 'Hg', 'Hb', 'OIII', 
                          'Ha', 'SII', 'SIII', 'SIII','HeI', 'PaB'])
line_lambda_AA = np.array([1215.7, 1240, 1487, 1750, 1640.4, 1663, 1908, 1906, 1549, 2800, 3346, 3426, 
                            (3727.092+3729.875)/2, 3869.87, 
                            4102, 4340.5, 4861.3, 5006.8, 
                            6562.8, 6716, 9068.6, 9531.1, 10830, 12822])
line_lambda_um = line_lambda_AA * 1e-4 # unit in um

line_list1= np.array(['OII', 'NeIII', 'Hg', 'Hb','OIII', 'Ha', 'PaB', 'PaA'])
line_lambda1_AA = np.array([(3727.092+3729.875)/2, 3869.87, 4340.5, 4861.3, 5006.8, 6562.8, 12822, 18756])
line_list2 = np.array(['Lya', 'NV', 'NIV', 'NIII', 'HeII', 'OIII-1663', 'CIII', 'CIII', 'CIV', 'MgII', 'NeV', 'NeVI',  'Hd', 'SII', 'SIII', 'SIII','HeI'])
line_lambda2_AA = np.array([1215.7, 1240, 1487, 1750, 1640.4, 1663, 1908, 1906, 1549, 2800, 3346, 3426, 4102,  6716, 9068.6, 9531.1, 10830])




DEFAULT_LINE_LIST = ['PaA','PaB','PaG','PaD',
                     'HeI-1083', 'SIII', 'SIII', 'OII-7325', 'ArIII-7138',
                     'SII', 'SII', 'Ha', 'OI-6302','OI-6302', 'HeI-5877', 'OIII', 'OIII', 'Hb', 
                     'OIII-4363', 'Hg', 'Hd', 'H7', 'H8', 'H9', 'H10', 
                     'NeIII-3867', 'OII', 'NeVI-3426', 'NeV-3346', 'MgII', 
                     'CIV-1549', 'CIII-1906', 'CIII-1908', 'OIII-1663', 
                     'HeII-1640', 'NIII-1750', 'NIV-1487', 'NV-1240', 'Lya']


DEFAULT_LINE_WAV = [18756.3, 12821.7, 10941.2, 10052.2,
                    (10832.057+10833.306)/2, 9071.1, 9533.2, (7321.9+7332.21)/2, 7137.77,
                    6718.29, 6732.67, 6564.697, 6302.046, 6365.535, 5877.249, 5007, 4960.295, 4862.738,
                    4364.436,  4341.731, 4102.936, 3971.236, 3890.191, 3836.511, 3799.014, 
                    3869.87, (3727.092+3729.875)/2, 3426.85, 3343.5, 2799.117,
                    1549.48, 1906.683, 1908.734, 1665.85, 1640.4, 1750., 1487., 1240.81, 1215.4] 

def Take_input():
    INPUT = inputtxt.get("1.0", "end-1c")
    return int(INPUT)

def New_Window(oned):
    Window = tk.Toplevel()
    Window.title("1D spectra window")
    Window.grab_set()
    #canvas = tk.Canvas(Window, height=200, width=400)
    fig = Figure() #plt.figure(1)
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    ax.plot(oned[1].data['wave'], oned[1].data['flux']/oned[1].data['flat']/oned[1].data['pscale']*1e19, label='spectra')
    ax.plot(oned[1].data['wave'], oned[1].data['cont']/oned[1].data['flat']*1e19, label='fitted continuum')
    ax.plot(oned[1].data['wave'], oned[1].data['line']/oned[1].data['flat']*1e19, label='fitted emission line')
    ax.plot(oned[1].data['wave'], oned[1].data['err']/oned[1].data['flat']/oned[1].data['pscale']*1e19, label='error')

    ax.set_ylabel('$f_\lambda 10^{-19}$ erg/s/cm$^2/AA$', fontsize=16)
    ax.set_xlabel('$\lambda$ AA', fontsize=16)
    ax.legend(fontsize=12)
    
    canvas = FigureCanvasTkAgg(fig, master=Window)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, Window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    

    button = tk.Button(master=Window, text="Quit", command=Window.destroy)
    button.pack(side=tk.BOTTOM)
    
    #Window.after(100, plt.close)
    

def New_Window2(oned, zgrizli):
    Window = tk.Toplevel()
    Window.title("1D spectra window (smoothed)")
    Window.grab_set()
    #canvas = tk.Canvas(Window, height=200, width=400)
    fig = Figure() #plt.figure(1)
    ax = fig.add_axes([0.15, 0.2, 0.8, 0.7])
    
    kernel_size = 3			# pixel size is ~40A, resolution of the WFC3/G141 grism is ~130 -> 110A, smoothing at 3 pixels to make ~ matched filter 
    kernel = np.ones(kernel_size) / kernel_size
    smooth1d = np.convolve(oned[1].data['flux']/oned[1].data['flat']/oned[1].data['pscale']*1e19, kernel, mode='same')
    smooth1derr = np.convolve(oned[1].data['err']/oned[1].data['flat']/oned[1].data['pscale']*1e19, kernel, mode='same')
    ax.plot(oned[1].data['wave'], smooth1d, label='smoothed spectra (3 pix)')
    ax.plot(oned[1].data['wave'], oned[1].data['cont']/oned[1].data['flat']*1e19, label='fitted continuum')
    ax.plot(oned[1].data['wave'], oned[1].data['line']/oned[1].data['flat']*1e19, label='fitted emission line')
    ax.plot(oned[1].data['wave'], smooth1derr, label='smoothed error (3 pix)')
    ax.set_xlim(oned[1].data['wave'].min(), oned[1].data['wave'].max())
    ax.set_ylim(ax.get_ylim())
    xmin, xmax = ax.get_xlim()
    
    plot_lines, plot_texts = [], []
    for i, l in enumerate(line_list1):
        li, = ax.plot([line_lambda1_AA[i]*(1+zgrizli), line_lambda1_AA[i]*(1+zgrizli)], 
                      [ax.get_ylim()[0], ax.get_ylim()[1]], ls='--', color='grey')
        if (line_lambda1_AA[i]*(1+zgrizli) > xmin) & (line_lambda1_AA[i]*(1+zgrizli) < xmax):
            tx = ax.text(line_lambda1_AA[i]*(1+zgrizli)+2, (ax.get_ylim()[1]-ax.get_ylim()[0])*0.6+ax.get_ylim()[0], l)
        else:
            tx = ax.text(0, (ax.get_ylim()[1]-ax.get_ylim()[0])*0.6+ax.get_ylim()[0], l)
        plot_lines.append(li)
        plot_texts.append(tx)
        
    plot_lines2, plot_texts2 = [], []
    for i, l in enumerate(line_list2):
        li, = ax.plot([line_lambda2_AA[i]*(1+zgrizli), line_lambda2_AA[i]*(1+zgrizli)], 
                      [ax.get_ylim()[0], ax.get_ylim()[1]], ls='--', color='grey', alpha=0.5)
        if (line_lambda2_AA[i]*(1+zgrizli) > xmin) & (line_lambda2_AA[i]*(1+zgrizli) < xmax):
            tx = ax.text(line_lambda2_AA[i]*(1+zgrizli)+2, (ax.get_ylim()[1]-ax.get_ylim()[0])*0.6+ax.get_ylim()[0], l, fontsize=8, alpha=0.5)
        else:
            tx = ax.text(0, (ax.get_ylim()[1]-ax.get_ylim()[0])*0.6+ax.get_ylim()[0], l, fontsize=8, alpha=0.5)
        plot_lines2.append(li)
        plot_texts2.append(tx)
        #print(line_lambda[i]*(1+zgrizli))
        
    ax.set_ylabel('$f_\lambda 10^{-19}$ erg/s/cm$^2/AA$', fontsize=16)
    ax.set_xlabel('$\lambda$ AA', fontsize=16)
    ax.legend(fontsize=12)

    canvas = FigureCanvasTkAgg(fig, master=Window)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(canvas, Window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    #Window.after(100, plt.close)
    # Choose the Slider color
    slider_color = 'White'
 
    # Set the axis and slider position in the plot
    axis_position = fig.add_axes([0.2, 0.05, 0.65, 0.03], facecolor = slider_color)
    slider_position = Slider(axis_position, 'z', 0.1, 6, valinit=zgrizli)
    slider_position.label.set_size(16)
    
    def update(val):
        tmpz = slider_position.val
        for i, li in enumerate(plot_lines):
            li.set_xdata([line_lambda1_AA[i] * (1+tmpz), line_lambda1_AA[i] * (1+tmpz)])
            if (line_lambda1_AA[i]*(1+tmpz) > xmin) & (line_lambda1_AA[i]*(1+tmpz) < xmax):
                plot_texts[i].set_position((line_lambda1_AA[i]*(1+tmpz)+2, (ax.get_ylim()[1]-ax.get_ylim()[0])*0.6+ax.get_ylim()[0]))
            else:
                plot_texts[i].set_position((0, (ax.get_ylim()[1]-ax.get_ylim()[0])*0.6+ax.get_ylim()[0]))
        for i, li in enumerate(plot_lines2):
            li.set_xdata([line_lambda2_AA[i] * (1+tmpz), line_lambda2_AA[i] * (1+tmpz)])
            if (line_lambda2_AA[i]*(1+tmpz) > xmin) & (line_lambda2_AA[i]*(1+tmpz) < xmax):
                plot_texts2[i].set_position((line_lambda2_AA[i]*(1+tmpz)+2, (ax.get_ylim()[1]-ax.get_ylim()[0])*0.6+ax.get_ylim()[0]))
            else:
                plot_texts2[i].set_position((0, (ax.get_ylim()[1]-ax.get_ylim()[0])*0.6+ax.get_ylim()[0]))
        
        fig.canvas.draw_idle()
    # update function called using on_changed() function
    slider_position.on_changed(update)  
    
    button = tk.Button(master=Window, text="Quit", command=Window.destroy)
    button.pack(side=tk.BOTTOM)
   
   
def Make_Gui(grizli_id, exist_flags, progress, oned, mag, z=[0, 0, 0], path="Extractions_peak1"):
    
    if np.isnan(exist_flags[0]):
        exist_flags = [-1,-1,-1,-1, '', 0, '']
        
    # ORIGINAL_DPI = 72.05403377110693
    # def get_dpi():
    #     screen = tk.Tk()
    #     current_dpi = screen.winfo_fpixels('1i')
    #     screen.destroy()
    #     print(current_dpi)
    #     return current_dpi
    # SCALE = get_dpi()/ORIGINAL_DPI
    # def scaled(original_width):
    #     return round(original_width * SCALE)

    root = tk.Tk()
    #root.geometry(f'{scaled(1350)}x{scaled(850)}') 
    root.geometry('1350x850')
    root.title("Main window")
    OPTIONS = ["    ", "Great", "Good", "Unclear", "Bad"]
    
    ### spectra quality and fitting quality
    label = tk.Label(root, text='spectra quality:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.place(x=450, y=10)
    
    var0 = tk.StringVar()
    var0.set(OPTIONS[int(exist_flags[0]+1)])
    #inputtxt = tk.Entry(root, textvariable=var0, width = 5, font=fnt.Font(size = 16))
    w0 = tk.OptionMenu(root, var0, *OPTIONS)
    w0.place(x=580, y=10) 

    label = tk.Label(root, text='spectra fitting:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.place(x=450, y=40)

    var1 = tk.StringVar()
    var1.set(OPTIONS[int(exist_flags[1]+1)])
    inputtxt = tk.Entry(root, textvariable=var1, width = 5, font=fnt.Font(size = 16))
    w1 = tk.OptionMenu(root, var1, *OPTIONS)
    w1.place(x=580, y=40) 
    
    label = tk.Label(root, text='phot fitting:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.place(x=450, y=70)
    
    var2 = tk.StringVar()
    var2.set(OPTIONS[int(exist_flags[2]+1)])
    #inputtxt = tk.Entry(root, textvariable=var3, width = 5, font=fnt.Font(size = 16), )
    w2 = tk.OptionMenu(root, var2, *OPTIONS)
    w2.place(x=580, y=70) 
    
    label = tk.Label(root, text='spec+phot fitting:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.place(x=450, y=100)
    
    var3 = tk.StringVar()
    var3.set(OPTIONS[int(exist_flags[3]+1)])
    #inputtxt = tk.Entry(root, textvariable=var3, width = 5, font=fnt.Font(size = 16), )
    w3 = tk.OptionMenu(root, var3, *OPTIONS)
    w3.place(x=580, y=100) 
    #inputtxt.place(x=580, y=70) 

    label = tk.Label(root, text='Comments:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.place(x=450, y=130)
    
    var4 = tk.StringVar()
    inputtxt = tk.Entry(root, textvariable=var4, width = 15, font=fnt.Font(size=16),)
    inputtxt.insert(0, exist_flags[4]) 
    inputtxt.place(x=540, y=130) 
    
    
    ### refit option
    label = tk.Label(root, text='Refit option [comments format]:', font=fnt.Font(size = 16),  anchor="w", justify="left")
    label.place(x=800, y=10)
    
    varr1 = tk.IntVar()
    chk1 = tk.Checkbutton(root, text='Incorrect scale on spectra vs. photometry',  font=fnt.Font(size = 16), variable=varr1)
    if exist_flags[5]//1000 == 1:
        chk1.select()
        exist_flags[5] = exist_flags[5]-1000
    chk1.place(x=800, y=30)
    
    varr2 = tk.IntVar()
    chk2 = tk.Checkbutton(root, text='clip bad photometric points [obs wavelength in um]',  font=fnt.Font(size = 16), variable=varr2)
    if exist_flags[5]//100 == 1:
        chk2.select()
        exist_flags[5] = exist_flags[5]-100
    chk2.place(x=800, y=53)
    
    varr3 = tk.IntVar()
    chk3 = tk.Checkbutton(root, text='refit with different redshift range [zmin zmax]',  font=fnt.Font(size = 16), variable=varr3)
    if exist_flags[5]//10 == 1:
        chk3.select()
        exist_flags[5] = exist_flags[5]-10
    chk3.place(x=800, y=76)
    
    varr4 = tk.IntVar()
    chk4 = tk.Checkbutton(root, text='refit [write reason in comment]',  font=fnt.Font(size = 16), variable=varr4)
    if exist_flags[5]//1 == 1:
        chk4.select()
    chk4.place(x=800, y=101)
    
    label = tk.Label(root, text='Refit comments:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.place(x=800, y=130)
    
    varr5 = tk.StringVar()
    inputtxt = tk.Entry(root, textvariable=varr5, width = 15, font=fnt.Font(size=16))
    inputtxt.insert(0, exist_flags[6]) 
    inputtxt.place(x=930, y=130) 
    
    
    label = tk.Label(root, text='Checking N = %i (total %i): \n'%(progress[0], progress[1]) 
                     + 'Grizli id: ' + grizli_id 
                     + ' (%.1f mag)'%mag 
                     + '\nGrizli fitted to z = %.2f [green in ds9]'%z[0] 
                     + '\nC20 zphot = %.2f [orange in ds9]'%z[1], 
                     font=fnt.Font(size = 18), anchor="w", justify="left")
    label.place(x=20, y=10)
    

    path_img = path+"%s.%s.png"%(grizli_id, 'full')
    img2 = ImageTk.PhotoImage(Image.open(path_img))
    tk.Label(root, image=img2).place(x=20, y=175)

    path_img = path+"%s.%s.png"%(grizli_id, 'stack')
    img3 = ImageTk.PhotoImage(Image.open(path_img))
    tk.Label(root, image=img3).place(x=840, y=175)

    path_img = path+"%s.%s.png"%(grizli_id, 'sed')
    img1 = ImageTk.PhotoImage(Image.open(path_img))
    tk.Label(root, image=img1).place(x=20, y=530)

    # Add information about wavelengths of strong features 
    line_list_gui = np.array(['OII', 'OII', 'NeIII', 'Hg', 'Hb','OIII', 'OIII', 'Ha', 'PaB', 'PaA'])
    line_lambda_gui = np.array([3726.0, 3728.8, 3868, 4361, 4861.3, 4958.9, 5006.8, 6562.8, 12822, 18756])
    OIIobs = (line_lambda_gui[0] + line_lambda_gui[1])/2*(1+z[0])
    NeIIIobs = line_lambda_gui[2]*(1+z[0])
    Hgobs = line_lambda_gui[3]*(1+z[0])
    Hbobs = line_lambda_gui[4]*(1+z[0])
    OIIIbobs = line_lambda_gui[5]*(1+z[0])
    OIIIrobs = line_lambda_gui[6]*(1+z[0])
    Haobs = line_lambda_gui[7]*(1+z[0])
    Pbobs = line_lambda_gui[8]*(1+z[0])
    Paobs = line_lambda_gui[9]*(1+z[0])

    label = tk.Label(root, text='lambda_obs of emission \nfeatures at z_Grizli:', font=fnt.Font(size = 18), anchor="w", justify="left")
    label.place(x=1130, y=530)

    label = tk.Label(root, text='[OII]: %.1f'%OIIobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#6C2DC7', bg='#FFFFFF')
    label.place(x=1130, y=580)

    label = tk.Label(root, text='[NeIII]: %.1f'%NeIIIobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#5865F2', bg='#FFFFFF')
    label.place(x=1130, y=605)

    label = tk.Label(root, text='Hg: %.1f'%Hgobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#16E2F5', bg='#FFFFFF')
    label.place(x=1130, y=630)

    label = tk.Label(root, text='Hb: %.1f'%Hbobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#1B8A6B',bg='#FFFFFF')
    label.place(x=1130, y=655)

    label = tk.Label(root, text='[OIII]: %.1f'%OIIIbobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#16F529', bg='#FFFFFF')
    label.place(x=1130, y=680)

    label = tk.Label(root, text='[OIII]: %.1f'%OIIIrobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#E9AB17', bg='#FFFFFF')
    label.place(x=1130, y=705)

    label = tk.Label(root, text='Ha: %.1f'%Haobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#FF8C00', bg='#FFFFFF')
    label.place(x=1130, y=730)

    label = tk.Label(root, text='Pb: %.1f'%Pbobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#E42217', bg='#FFFFFF')
    label.place(x=1130, y=755)

    label = tk.Label(root, text='Pa: %.1f'%Paobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#800517', bg='#FFFFFF')
    label.place(x=1130, y=780)
    
    ### Add 1d spectra window 
    button = tk.Button(root, text="plot 1D spectra", bg='White', fg='dark orange', font=fnt.Font(size=20), 
                              command=lambda: New_Window(oned))
    button.place(x=20, y=100)
    
    button = tk.Button(root, text="plot 1D spectra (smoothed)", bg='White', fg='dark orange', font=fnt.Font(size=20), 
                              command=lambda: New_Window2(oned, z[0]))
    button.place(x=20, y=140)
    
    
    ### Add exit and next button 
    def Close(): 
        np.savetxt('.exit_tmp', [1])
        root.destroy() 
        
    def Next(): 
        np.savetxt('.exit_tmp', [0])
        root.destroy() 
        
    exit_button = tk.Button(root, text="Exit", command=Close, font=fnt.Font(size=20))  #command=lambda m="It is an apple": which_button(m)
    exit_button.place(x=1250, y=60)
    
    next_button = tk.Button(root, text="Next", command=Next,  font=fnt.Font(size=20))  #command=lambda m="It is an apple": which_button(m)
    next_button.place(x=1250, y=20)
    ###
    
    root.mainloop()
    
    exit_param = np.loadtxt('.exit_tmp')
    
    spectra_flag = var0.get()
    spec_fit_flag = var1.get()
    phot_fit_flag = var2.get()
    specphot_fit_flag = var3.get()
    
    
    refit_flag = 0
    if varr1.get() == 1:
        refit_flag += 1000
    if varr2.get() == 1:
        refit_flag += 100
    if varr3.get() == 1:
        refit_flag += 10
    if varr4.get() == 1:
        refit_flag += 1
    
    return [spectra_flag, spec_fit_flag, phot_fit_flag, specphot_fit_flag], var4.get(), refit_flag, varr5.get(), exit_param


def Make_Gui_BCL(grizli_id, exist_flags, progress, oned, mag, z=[0, 0, 0], path="Extractions_peak1"):
    
    if np.isnan(exist_flags[0]):
        exist_flags = [-1,-1,-1,-1, '', 0, '']
         
    root = tk.Tk()
    root.geometry('1750x850')
    root.title("Main window")
    OPTIONS = ["    ", "Great", "Good", "Unclear", "Bad"]

    label = tk.Label(root, text='spectra quality:', font=fnt.Font(size = 15), anchor="w", justify="left")
    label.place(x=620, y=10)
 
    var0 = tk.StringVar()
     #inputtxt = tk.Entry(root, textvariable=var0, width = 5, font=fnt.Font(size = 15))
    var0.set(OPTIONS[0])
    w0 = tk.OptionMenu(root, var0, *OPTIONS)
    w0.place(x=800, y=10) 
     
    label = tk.Label(root, text='spectra fitting:', font=fnt.Font(size = 15), anchor="w", justify="left")
    label.place(x=620, y=40)

    var1 = tk.StringVar()
    var1.set(OPTIONS[0])
    w1 = tk.OptionMenu(root, var1, *OPTIONS)
    w1.place(x=800, y=40) 

    label = tk.Label(root, text='phot fitting:', font=fnt.Font(size = 15), anchor="w", justify="left")
    label.place(x=620, y=70)

    var2 = tk.StringVar()
    var2.set(OPTIONS[0])
    w2 = tk.OptionMenu(root, var2, *OPTIONS) #tk.Entry(root, textvariable=var2, width = 5, font=fnt.Font(size = 15))
    w2.place(x=800, y=70)

    label = tk.Label(root, text='phot+spec fitting:', font=fnt.Font(size = 15), anchor="w", justify="left")
    label.place(x=620, y=100)

    var3 = tk.StringVar()
    var3.set(OPTIONS[0])
    w3 = tk.OptionMenu(root, var3, *OPTIONS) #(root, textvariable=var3, width = 5, font=fnt.Font(size = 15))
    w3.place(x=800, y=100) 

    # label = tk.Label(root, text='good=0, bad=1, unclear=2', font=fnt.Font(size = 15), anchor="w", justify="left")
    # label.place(x=620, y=130)

    label = tk.Label(root, text='Comments:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.place(x=620, y=130)

    var4 = tk.StringVar()
    inputtxt = tk.Entry(root, textvariable=var4, width = 24, font=fnt.Font(size=15),)
    inputtxt.place(x=715, y=130) 


    label = tk.Label(root, text='Refit option [comments format]:', font=fnt.Font(size = 15), anchor="w", justify="left")
    label.place(x=1000, y=10)

    varr1 = tk.IntVar()
    chk1 = tk.Checkbutton(root, text='Incorrect scale on spectra vs. photometry',  font=fnt.Font(size = 15), variable=varr1)
    chk1.place(x=1000, y=30)

    varr2 = tk.IntVar()
    chk2 = tk.Checkbutton(root, text='Clip bad photometric points [obs wavelength in um]',  font=fnt.Font(size = 15), variable=varr2)
    chk2.place(x=1000, y=53)

    varr3 = tk.IntVar()
    chk3 = tk.Checkbutton(root, text='Refit with different redshift range [zmin zmax]',  font=fnt.Font(size = 15), variable=varr3)
    chk3.place(x=1000, y=76)

    varr4 = tk.IntVar()
    chk4 = tk.Checkbutton(root, text='Refit [write reason in comment]',  font=fnt.Font(size = 15), variable=varr4)
    chk4.place(x=1000, y=101)

    label = tk.Label(root, text='Comments:', font=fnt.Font(size = 15), anchor="w", justify="left")
    label.place(x=1000, y=130)

    varr5 = tk.StringVar()
    inputtxt = tk.Entry(root, textvariable=varr5, width = 40, font=fnt.Font(size=15))
    inputtxt.place(x=1110, y=130) 

    label = tk.Label(root, text='Checking N = %i (total %i): \n'%(progress[0], progress[1]) 
                        + 'Grizli id: ' + grizli_id  + ' (%.1f mag)'%mag 
                        + '\nGrizli fitted to z = %.2f [green in ds9]'%z[0] 
                        + '\nC20 zphot = %.2f [orange in ds9]'%z[1], font=fnt.Font(size = 16), anchor="w", justify="left")
    label.place(x=20, y=10)

    # Add information about wavelengths of strong features 
    line_list_gui = np.array(['OII', 'OII', 'NeIII', 'Hg', 'Hb','OIII', 'OIII', 'Ha', 'PaB', 'PaA'])
    line_lambda_gui = np.array([3726.0, 3728.8, 3868, 4361, 4861.3, 4958.9, 5006.8, 6562.8, 12822, 18756])
    OIIobs = (line_lambda_gui[0] + line_lambda_gui[1])/2*(1+z[0])
    NeIIIobs = line_lambda_gui[2]*(1+z[0])
    Hgobs = line_lambda_gui[3]*(1+z[0])
    Hbobs = line_lambda_gui[4]*(1+z[0])
    OIIIbobs = line_lambda_gui[5]*(1+z[0])
    OIIIrobs = line_lambda_gui[6]*(1+z[0])
    Haobs = line_lambda_gui[7]*(1+z[0])
    Pbobs = line_lambda_gui[8]*(1+z[0])
    Paobs = line_lambda_gui[9]*(1+z[0])

    label = tk.Label(root, text='lambda_obs of emission features at z_Grizli:', font=fnt.Font(size = 18), anchor="w", justify="left")
    label.place(x=1220, y=540)

    label = tk.Label(root, text='[OII]: %.1f'%OIIobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#6C2DC7', bg='#FFFFFF')
    label.place(x=1230, y=580)

    label = tk.Label(root, text='[NeIII]: %.1f'%NeIIIobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#5865F2', bg='#FFFFFF')
    label.place(x=1230, y=605)

    label = tk.Label(root, text='Hg: %.1f'%Hgobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#16E2F5', bg='#FFFFFF')
    label.place(x=1230, y=630)

    label = tk.Label(root, text='Hb: %.1f'%Hbobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#1B8A6B',bg='#FFFFFF')
    label.place(x=1230, y=655)

    label = tk.Label(root, text='[OIII]: %.1f'%OIIIbobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#16F529', bg='#FFFFFF')
    label.place(x=1230, y=680)

    label = tk.Label(root, text='[OIII]: %.1f'%OIIIrobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#E9AB17', bg='#FFFFFF')
    label.place(x=1230, y=705)

    label = tk.Label(root, text='Ha: %.1f'%Haobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#FF8C00', bg='#FFFFFF')
    label.place(x=1230, y=730)

    label = tk.Label(root, text='Pb: %.1f'%Pbobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#E42217', bg='#FFFFFF')
    label.place(x=1230, y=755)

    label = tk.Label(root, text='Pa: %.1f'%Paobs + ' A', font=fnt.Font(size = 17), anchor="w", justify="left", fg='#800517', bg='#FFFFFF')
    label.place(x=1230, y=780)
    
    ### new added window
    button = tk.Button(root, text="plot 1D spectra", bg='White', fg='dark orange', font=fnt.Font(size=15), 
                                command=lambda: New_Window(oned))
    button.place(x=20, y=100)
        ### Add 1d spectra smoothed window  
    button = tk.Button(root, text="plot 1D spectra (smoothed)", bg='White', fg='dark green', font=fnt.Font(size=15),
                                command=lambda: New_Window2(oned, z[0]))
    button.place(x=20, y=140)


    ###
    ### Add exit and next button 
    def Close(): 
        np.savetxt('.exit_tmp', [1])
        root.destroy() 
        
    def Next(): 
        np.savetxt('.exit_tmp', [0])
        root.destroy() 
        
    exit_button = tk.Button(root, text="Exit", command=Close, font=fnt.Font(size=15))  #command=lambda m="It is an apple": which_button(m)
    exit_button.place(x=1550, y=60)

    next_button = tk.Button(root, text="Next", command=Next,  font=fnt.Font(size=15))  #command=lambda m="It is an apple": which_button(m)
    next_button.place(x=1550, y=20)
    ###

    ### new added

    path_img = path+"%s.%s.png"%(grizli_id, 'full')
    img2r = Image.open(path_img)
    img2rr = img2r.resize((900, 325))
    img2 = ImageTk.PhotoImage(img2rr)
    tk.Label(root, image=img2).place(x=20, y=125+75)

    path_img = path+"%s.%s.png"%(grizli_id, 'stack')
    img3r = Image.open(path_img)
    img3rr = img3r.resize((750, 325))
    img3 = ImageTk.PhotoImage(img3rr)
    tk.Label(root, image=img3).place(x=950, y=125+75)

    path_img = path+"%s.%s.png"%(grizli_id, 'sed')
    img1r = Image.open(path_img)
    img1rr = img1r.resize((1200, 310))
    img1 = ImageTk.PhotoImage(img1rr)
    tk.Label(root, image=img1).place(x=20, y=480+50)


    root.mainloop()

    exit_param = np.loadtxt('.exit_tmp')

    spectra_flag = var0.get()
    spec_fit_flag = var1.get()
    phot_fit_flag = var2.get()
    specphot_fit_flag = var3.get()


    refit_flag = 0
    if varr1.get() == 1:
        refit_flag += 1000
    if varr2.get() == 1:
        refit_flag += 100
    if varr3.get() == 1:
        refit_flag += 10
    if varr4.get() == 1:
        refit_flag += 1

    return [spectra_flag, spec_fit_flag, phot_fit_flag, specphot_fit_flag], var4.get(), refit_flag, varr5.get(), exit_param


def Make_Gui2(grizli_id, exist_flags, progress, oned, mag, z=[0, 0, 0], path="Extractions_peak1"):
    
    if np.isnan(exist_flags[0]):
        exist_flags = [-1,-1,-1, '', 0, '']

    root = tk.Tk()
    root.geometry('1350x850') # width x height of the window
    #root.geometry('1200x630')
    root.title("Main window")
    OPTIONS = ["    ", "Great", "Good", "Unclear", "Bad"]
    
    ### spectra quality and fitting quality
    frame1 = tk.Frame(root)
    frame1.grid(column = 0, row = 0, columnspan = 1)
    frame2 = tk.Frame(root)
    frame2.grid(column = 1, row = 0, columnspan = 1)
    frame3 = tk.Frame(root)
    frame3.grid(column = 2, row = 0, columnspan = 1)
    frame4 = tk.Frame(root)
    frame4.grid(column = 3, row = 0, columnspan = 1)
    frame5 = tk.Frame(root)
    frame5.grid(column = 0, row = 1, columnspan = 4, sticky = tk.W+tk.E)
    frame6 = tk.Frame(root)
    frame6.grid(column = 0, row = 2, columnspan = 3, sticky = tk.W+tk.E)
    frame7 = tk.Frame(root)
    frame7.grid(column = 3, row = 2, columnspan = 1)
    
    label = tk.Label(frame1, text='Checking N = %i (total %i): \n'%(progress[0], progress[1]) 
                     + 'Grizli id: ' + grizli_id 
                     + ' (%.1f mag)'%mag 
                     + '\nGrizli fitted to z = %.2f [green in ds9]'%z[0] 
                     + '\nC20 zphot = %.2f [orange in ds9]'%z[1], 
                     font=fnt.Font(size = 16), anchor="w", justify="left")
    label.pack(anchor=tk.W)
    ### Add 1d spectra window 
    button = tk.Button(frame1, text="plot 1D spectra", bg='White', fg='dark orange', font=fnt.Font(size=20), 
                              command=lambda: New_Window(oned))
    button.pack(anchor=tk.W)
    
    button = tk.Button(frame1, text="plot 1D spectra (smoothed)", bg='White', fg='dark orange', font=fnt.Font(size=20), 
                              command=lambda: New_Window2(oned, z[0]))
    button.pack(anchor=tk.W)
    
    
    label = tk.Label(frame2, text='spectra quality:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.grid(column = 0, row = 0, sticky = tk.W, pady=5)
    
    var0 = tk.StringVar()
    var0.set(OPTIONS[int(exist_flags[0]+1)])
    w0 = tk.OptionMenu(frame2, var0, *OPTIONS)
    w0.grid(column = 1, row = 0, sticky = tk.W)

    label = tk.Label(frame2, text='spectra fitting:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.grid(column = 0, row = 1, sticky = tk.W, pady=5)

    var1 = tk.StringVar()
    var1.set(OPTIONS[int(exist_flags[1]+1)])
    w1 = tk.OptionMenu(frame2, var1, *OPTIONS)
    w1.grid(column = 1, row = 1, sticky = tk.W)
    
    label = tk.Label(frame2, text='phot fitting:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.grid(column = 0, row = 2, sticky = tk.W, pady=5)
    
    var2 = tk.StringVar()
    var2.set(OPTIONS[int(exist_flags[2]+1)])
    w2 = tk.OptionMenu(frame2, var2, *OPTIONS)
    w2.grid(column = 1, row = 2, sticky = tk.W, pady=5)
    
    # label = tk.Label(frame2, text='spec+phot fitting:', font=fnt.Font(size = 16), anchor="w", justify="left")
    # label.pack(anchor=tk.W)

    # var3 = tk.StringVar()
    # var3.set(OPTIONS[int(exist_flags[3]+1)])
    # w3 = tk.OptionMenu(frame22, var3, *OPTIONS)
    # w3.pack(anchor=tk.W)

    label = tk.Label(frame2, text='Comments:', font=fnt.Font(size = 16), anchor="w", justify="left")
    label.grid(column = 0, row = 3, sticky = tk.W+tk.E)
    
    var4 = tk.StringVar()
    inputtxt = tk.Entry(frame2, textvariable=var4, width=20, font=fnt.Font(size=16))
    inputtxt.insert(0, exist_flags[3]) 
    inputtxt.grid(column = 0, row = 4, columnspan=2, sticky = tk.W+tk.E)
    
    ### refit option
    label = tk.Label(frame3, text='Refit option [comments format]:', font=fnt.Font(size = 16),  anchor="w", justify="left")
    label.grid(column =0, row = 0, columnspan=2, sticky = tk.W+tk.E)
    
    varr1 = tk.IntVar()
    chk1 = tk.Checkbutton(frame3, text='refit with different redshift range [zmin zmax]',  font=fnt.Font(size = 16), variable=varr1, anchor="w")
    if exist_flags[4]//10000 == 1:
        chk1.select()
        exist_flags[4] = exist_flags[4]-10000
    chk1.grid(column = 0, row = 1, columnspan=2, sticky = tk.W+tk.E)
    
    varr2 = tk.IntVar()
    chk2 = tk.Checkbutton(frame3, text=u'refit only portion of spectrum [\u03bbobs1 \u03bbobs2 in \u03bcm]',  font=fnt.Font(size = 16), variable=varr2, anchor="w")
    if exist_flags[4]//1000 == 1:
        chk2.select()
        exist_flags[4] = exist_flags[4]-1000
    chk2.grid(column = 0, row = 2, columnspan=2, sticky = tk.W+tk.E)
    
    varr3 = tk.IntVar()
    chk3 = tk.Checkbutton(frame3, text='Incorrect scale on spectra vs. photometry',  font=fnt.Font(size = 16), variable=varr3, anchor="w")
    if exist_flags[4]//100 == 1:
        chk3.select()
        exist_flags[4] = exist_flags[4]-100
    chk3.grid(column = 0, row = 3, columnspan=2, sticky = tk.W+tk.E)
    
    varr4 = tk.IntVar()
    chk4 = tk.Checkbutton(frame3, text='clip bad photometric points [\u03bbobs in \u03bcm]',  font=fnt.Font(size = 16), variable=varr4, anchor="w")
    if exist_flags[4]//10 == 1:
        chk4.select()
        exist_flags[4] = exist_flags[4]-10
    chk4.grid(column = 0, row = 4, columnspan=2, sticky = tk.W+tk.E)
    
    varr5 = tk.IntVar()
    chk5 = tk.Checkbutton(frame3, text='other reason',  font=fnt.Font(size = 16), variable=varr5, anchor="w")
    if exist_flags[4]//1 == 1:
        chk5.select()
    chk5.grid(column = 0, row = 5, columnspan=2, sticky = tk.W+tk.E)
    
    
    label = tk.Label(frame3, text='Refit comments:', font=fnt.Font(size = 16), anchor="w")
    label.grid(column = 0, row = 6, sticky = tk.W+tk.E)
    
    varr6 = tk.StringVar()
    inputtxt = tk.Entry(frame3, textvariable=varr6, width = 15, font=fnt.Font(size=16))
    inputtxt.insert(0, exist_flags[5]) 
    inputtxt.grid(column = 1, row = 6, sticky = tk.W+tk.E)
    
    
    ### Add exit and next button 
    def Close(): 
        np.savetxt('.exit_tmp', [1])
        root.destroy() 
        
    def Next(): 
        np.savetxt('.exit_tmp', [0])
        root.destroy() 
        
    next_button = tk.Button(frame4, text="Next", command=Next,  font=fnt.Font(size=20, weight='bold'), fg='dark orange')  
    next_button.pack(anchor=tk.W, pady=10)
    exit_button = tk.Button(frame4, text="Exit", command=Close, font=fnt.Font(size=20, weight='bold'), fg='dark orange') 
    exit_button.pack(anchor=tk.W, pady=10)
    
    ###display png
    
    path_img = path+"%s.%s.png"%(grizli_id, 'full')
    img2 = ImageTk.PhotoImage(Image.open(path_img))
    tk.Label(frame5, image=img2).pack(side='left')

    path_img = path+"%s.%s.png"%(grizli_id, 'stack')
    img3 = ImageTk.PhotoImage(Image.open(path_img))
    tk.Label(frame5, image=img3).pack(side='right')

    path_img = path+"%s.%s.png"%(grizli_id, 'sed')
    img1 = ImageTk.PhotoImage(Image.open(path_img))
    tk.Label(frame6, image=img1).pack(side='left')
    
    # Add information about wavelengths of strong features 
    line_list_gui = np.array(['OII', 'NeIII', 'Hg', 'Hb','OIII', 'OIII', 'Ha', 'PaB', 'PaA'])
    line_lambda_gui = np.array([3727.4, 3868, 4361, 4861.3, 4958.9, 5006.8, 6562.8, 12822, 18756])
    line_lambda_zobs =  line_lambda_gui * (1+z[0])
    line_colors = ['#6C2DC7', '#5865F2', '#16E2F5', '#1B8A6B', '#16F529', '#E9AB17','#FF8C00', '#E42217', '#800517']
    
    
    label = tk.Label(frame7, text='lambda_obs of emission \nfeatures at z_Grizli:', font=fnt.Font(size = 18), anchor="w", justify="left")
    label.pack()
    
    #labels = []
    for i in range(len(line_list_gui)):
        label = tk.Label(frame7, 
                               text='%s: %.1f'%(line_list_gui[i], line_lambda_zobs[i]) + ' A', 
                               font=fnt.Font(size = 17), anchor="w", justify="left", fg=line_colors[i], bg='#FFFFFF')
        label.pack(anchor=tk.W)

    root.mainloop()
    
    exit_param = np.loadtxt('.exit_tmp')
    
    spectra_flag = var0.get()
    spec_fit_flag = var1.get()
    phot_fit_flag = var2.get()
    #specphot_fit_flag = var3.get()
    
    refit_flag = 0
    if varr1.get() == 1:
        refit_flag += 10000
    if varr2.get() == 1:
        refit_flag += 1000
    if varr3.get() == 1:
        refit_flag += 100
    if varr4.get() == 1:
        refit_flag += 10
    if varr5.get() == 1:
        refit_flag += 1
        
    return [spectra_flag, spec_fit_flag, phot_fit_flag], var4.get(), refit_flag, varr6.get(), exit_param

def toregion(xpix, label, outreg="line.reg", reg_color="green"):
    """
    
    reg_color: "red", "green", "blue", "cyan", "magenta", "pink", "yellow", "white", "black"
    """     

    FILE=open(outreg,'w')
    FILE.write("# Region file format: DS9 version 4.1\n")
    FILE.write('global color=%s dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'%reg_color)
    
    for i in range(0,len(xpix)):
        #Write to ds9 region file
        FILE.write("line(" + str(xpix[i]) + "," + str(10) + "," + str(xpix[i])  + ",25)\n")
        FILE.write("line(" + str(xpix[i]) + "," + str(40) + "," + str(xpix[i])  + ",55) # text={" + label[i] + "}\n")
    FILE.close()
    
def getreg(id, z_list, PATH):
    ### z_list = [z_grizli, z_phot, z_spec]
    
    grism_limits = [1.06, 1.73, 0.0046] #unit in um
    
    #line_list = ['PaA':[18700],'PaB':[12800], 'Ha':[6562.819], 'OIII':[4958.911, 5006.843], 'Hb':[4861.333], 'OII':[3726.032, 3728.815]]
    
    region_list = []
    for j, zj in enumerate(z_list):
        if zj > 0:
            line_obs_pix = (line_lambda_um * (1+zj) - 1.06) / 0.0046 + 1 #
            idx = (line_obs_pix>0) & (line_obs_pix < 146)
            toregion(line_obs_pix[idx], label=line_list[idx], outreg=PATH+id+"_line_%s.reg"%z_list_label[j], reg_color=z_list_color[j])
            region_list = np.append(region_list, PATH+id+"_line_%s.reg"%z_list_label[j])
    return region_list

def getFWHM(stack_path):
    hdu = pyfits.open(stack_path) 
    y = np.sum(hdu[5].data, axis=1)
    x = np.arange(0, 64)

    g_init = models.Gaussian1D(amplitude=hdu[5].data.max(), mean=31.5, stddev=2)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y)
    a, b = (g.mean.value - g.stddev.value*1.175)+1, (g.mean.value + g.stddev.value*1.175)+1
    
    outreg = stack_path[:-5]+"_FWHM.reg"
    FILE=open(outreg,'w')
    FILE.write("# Region file format: DS9 version 4.1\n")
    FILE.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    
    FILE.write("line(" + str(1) + "," + '%.3f'%a + "," + str(146)  + ","+ '%.3f'%a + ")\n")
    FILE.write("line(" + str(1) + "," + '%.3f'%b + "," + str(146)  + ","+ '%.3f'%b + ")\n")
    FILE.close()
    hdu.close()
    #print(g.mean.value, g.stddev.value, g.stddev.value*1.175*2)

def get_modelsub(stack_path):
    hdu = pyfits.open(stack_path) 
    remove_list = ['CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2']
    
    if len(hdu) <= 10:
        fits_frame = 1
        modelsub0 = hdu[0].copy()
        modelsub1 = hdu[1].copy()        
        modelsub1.data = hdu[1].data - hdu[4].data
        for l in remove_list:
            if l in modelsub1.header:
                del modelsub1.header[l]

        modelsub = pyfits.HDUList([modelsub0, modelsub1])
        modelsub.writeto(stack_path[:-5]+"_modelsub.fits", overwrite=True)
    else:
        fits_frame = 2
        modelsub0 = hdu[0].copy()
        modelsub1 = hdu[1].copy()
        modelsub1.data = hdu[1].data - hdu[4].data
        for l in remove_list:
            if l in modelsub1.header:
                del modelsub1.header[l]
                
        modelsub2 = hdu[6].copy()
        modelsub2.data = hdu[6].data - hdu[9].data
        for l in remove_list:
            if l in modelsub2.header:
                del modelsub2.header[l]
        
        modelsub = pyfits.HDUList([modelsub0, modelsub1, modelsub2])
        modelsub.writeto(stack_path[:-5]+"_modelsub.fits", overwrite=True)
    hdu.close()
    return fits_frame

def getparams(argv):
   inputfile = ''
   n0 = 0
   id_list = []
   opts, arg = getopt.getopt(argv,"hi:n:l:",["ifile= ","number=", "idlist"])
   for opt, arg in opts:
        if opt == '-h':
            print ('check_results.py -i <inputfile> -n <startnumber> -l <grilziidlist>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-n", "--number"):
            n0 = int(arg)
        elif opt in ("-l", "--idlist"):
            #id_list.append(arg)
            id_list = json.loads(arg)
   print('Input file: ', inputfile)
   file_path = inputfile.split("/")[0]
   print('file folder: ', file_path)
   
   return inputfile, file_path, n0, id_list

if __name__ == "__main__":
    
    #opts, arg = getopt.getopt(sys.argv[1:],"hl:n:",["help", "lfile=","number="])
    #print(opts)
    inputfile, file_path, n0, id_list = getparams(sys.argv[1:])

#inputfile = 'hyperion_peak1_grizli_result_wphotz_wspecz_z2_z3_wflag.cat'
PATH_table = file_path+'/'
tab1 = pd.read_csv(inputfile, sep=',')
print("total number of objects in the input file is ", len(tab1))

#check exist result


if len(id_list) > 0:
    print('running grizli id ' , id_list) 
else:
    print('start from number ', n0)
       
### check previous flag catalog
if os.path.exists(inputfile[:-4]+'.addflags.cat'):
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    print(inputfile[:-4]+'.addflags.cat' + ' exist, move to ' + inputfile[:-4]+'.addflags_%s.cat'%dt_string)
    os.system('cp '+inputfile[:-4]+'.addflags.cat' + ' ' + inputfile[:-4]+'.addflags_%s.cat'%dt_string)  
    tab1 = pd.read_csv(inputfile[:-4]+'.addflags.cat', sep=',')
    #tab1.loc[tab1.spectra_quality_flag.isna(), 'spectra_quality_flag'] = -1
    #tab1.loc[tab1.spectra_fitting_flag.isna(), 'spectra_fitting_flag'] = -1
    #tab1.loc[tab1.phot_fitting_flag.isna(), 'phot_fitting_flag'] = -1
    #tab1.loc[tab1.specphot_fitting_flag.isna(), 'specphot_fitting_flag'] = -1
    
else:
    tab1['spectra_quality_flag'] = -1
    tab1['spectra_fitting_flag'] = -1
    tab1['phot_fitting_flag'] = -1
    tab1['specphot_fitting_flag'] = -1
    tab1['spectra_comment'] = ''
    tab1['refit_flag'] = 0 
    tab1['refit_comment'] = ''
    
### check all files exist:    
check_files = 0
for i in range(len(tab1)):
    id_i =  tab1.loc[i, 'grizli_ID'] 
    # peak = id_i.split("_")[1]
    files = [".stack.fits", ".1D.fits", ".full.png", ".sed.png", ".stack.png"]
    for f in files:
        if not os.path.exists(PATH_table+ id_i + f):
            print(PATH_table+ id_i + f + ' not exist')
            check_files = 1
if check_files == 1:
    sys.exit()
else:
    print("All files exist")
    
np.savetxt('.exit_tmp', [0]) # init a tmp file to store GUI exit: exit = 1, next = 0

for i in range(n0, len(tab1)):
    
    id_i, mag_i =  tab1.loc[i, 'grizli_ID'], tab1.loc[i, 'ir_mag_tot']#.split("_")
    # peak = id_i.split("_")[1]
    # PATH_table = file_path +'/' 
        
    if len(id_list) > 0:
        print(id_i, id_list)
        if id_i not in id_list:
            continue

    z_list = [tab1.loc[i, 'REDSHIFT'], tab1.loc[i, 'ez_z_phot']] #, tab1.loc[i, 'z_spec']]
    z_list_label = ["zgrism",  'zphot'] #, 'zspec'
    z_list_color = ["green", 'orange'] #, 'red'
    
    grism_limits = [1.06, 1.73, 0.0046] #unit in um
    
    region_list = getreg(id=id_i, z_list=z_list, PATH=PATH_table) 
   
    stack_path = PATH_table+ id_i + ".stack.fits"
    oned_path = PATH_table+ id_i + ".1D.fits"
    #fits_frame = 0
    #hdu = pyfits.open(stack_path) 
    oned = pyfits.open(oned_path)
    
    fits_frame = get_modelsub(stack_path)
    getFWHM(stack_path)
        
    d = pyds9.DS9()
    d.set('single')
    d.set('scale mode zscale')

    d.set('frame 1')
    d.set('file '+ stack_path+'[1]')
    d.set('regions', 'text(30, 5) # text={spectra #1} color=red font="helvetica 12 bold italic" ')
    for f in region_list[:2]:
        d.set("region "+ f)
    d.set("region "+stack_path[:-5]+"_FWHM.reg")
        
    d.set('frame 2')
    d.set('file '+ stack_path[:-5]+"_modelsub.fits[1]")
    d.set('regions', 'text(50, 5) # text={continuum sub #1} color=red font="helvetica 12 bold italic" ')
    for f in region_list[:2]:
        d.set("region "+ f)
    d.set("region "+stack_path[:-5]+"_FWHM.reg")
    
    if fits_frame == 2:
        d.set('frame 3')
        d.set('file '+ stack_path+'[6]')
        d.set('regions', 'text(30, 5) # text={spectra #2} color=red font="helvetica 12 bold italic" ')
        for f in region_list[:2]:
            d.set("region "+ f)
        d.set("region "+stack_path[:-5]+"_FWHM.reg")
        
        d.set('frame 4')
        d.set('file '+ stack_path[:-5]+"_modelsub.fits[2]")
        d.set('regions', 'text(50, 5) # text={continuum sub #2} color=red font="helvetica 12 bold italic" ')
        for f in region_list[:2]:
            d.set("region "+ f)
        d.set("region "+stack_path[:-5]+"_FWHM.reg")
        
            
    d.set('lock frame IMAGE')
    d.set('smooth radius 4')
    d.set('smooth function gaussian')
    d.set('smooth yes')
    d.set('lock smooth')
    d.set('lock colorbar')

    flags_cols = ['spectra_quality_flag', 'spectra_fitting_flag', 'phot_fitting_flag', 'spectra_comment', 'refit_flag', 'refit_comment']
    exist_flags = tab1.loc[i, flags_cols].tolist()
    
    
    #flags, tab1.loc[i, 'spectra_comment'],  tab1.loc[i, 'refit_flag'], tab1.loc[i, 'refit_comment'], f = Make_Gui(grizli_id=id_i, exist_flags=exist_flags, progress=[i, len(tab1)], oned=oned, mag=mag_i, z=z_list, path=PATH_table)
    
    # The code below is for Brian's Linux 
    flags, tab1.loc[i, 'spectra_comment'],  tab1.loc[i, 'refit_flag'], tab1.loc[i, 'refit_comment'],  f = Make_Gui2(grizli_id=id_i, exist_flags=exist_flags, progress=[i, len(tab1)], oned=oned, mag=mag_i, z=z_list, path=PATH_table)
    
    flag_system = {"    ":-1,"Great":0, "Good":1, "Unclear":2, "Bad":3}
    flags_col = ['spectra_quality_flag', 'spectra_fitting_flag', 'phot_fitting_flag']
    for fi, flag_i in enumerate(flags):
        tab1.loc[i, flags_col[fi]] =flag_system[flag_i]
    
    # get spec+phot fitting flag
    # "Great": both “phot fitting” and “spec fitting” flags are ‘great’
    # "Good": both “phot fitting” and “spec fitting” flags are ‘good’ or one of the two is ‘good’ and the other is ‘great’
    # "Unclear": both “phot fitting” and “spec fitting” flags are ‘unclear’ or one of the two is ‘unclear’ and the other is ‘good/great’
    # "Bad": “phot fitting” and/or “spec fitting” flags are ‘bad’
    if (flags[1] == "Great") & (flags[2] == "Great"):
        tab1.loc[i, 'specphot_fitting_flag'] = 0
    elif ((flags[1] == "Good") & (flags[2] == "Good")) | ((flags[1] == "Good")&(flags[2] == "Great")) | ((flags[1] == "Great")&(flags[2] == "Good")):
        tab1.loc[i, 'specphot_fitting_flag'] = 1
    elif ((flags[1] == "Unclear") & (flags[2] == "Unclear")) | ((flags[1] == "Good")&(flags[2] == "Unclear")) | ((flags[1] == "Unclear")&(flags[2] == "Good")) |  ((flags[1] == "Great")&(flags[2] == "Unclear")) | ((flags[1] == "Unclear")&(flags[2] == "Great")):
        tab1.loc[i, 'specphot_fitting_flag'] = 2
    elif (flags[1] == "Bad") | (flags[2] == "Bad"):
        tab1.loc[i, 'specphot_fitting_flag'] = 3
    else:
        print('false combination for spec+phot fitting: ', id_i, ' ', flags)

    #tab1.loc[i, 'spectra_quality_flag'], tab1.loc[i, 'spectra_fitting_flag'], tab1.loc[i, 'phot_fitting_flag'], tab1.loc[i, 'specphot_fitting_flag'] = flags
    
    d.set('frame DELETE all')
    
    if i % 2 == 0:
        tab1.to_csv(inputfile[:-4]+'.addflags.cat', sep=',', index=False)
        
    if f == 1:
        tab1.to_csv(inputfile[:-4]+'.addflags.cat', sep=',', index=False)
        print("Save to " + inputfile[:-4]+'.addflags.cat')
        print("Exit at N=%i"%i+" (%s) "%id_i)
        
        for ff in glob.glob(PATH_table+"*.reg"):
            os.remove(ff)
        for ff in glob.glob(PATH_table+"*_modelsub.fits"):
            os.remove(ff)
        print("remove all middle stage files") 
        break
    
tab1.to_csv(inputfile[:-4]+'.addflags.cat', sep=',', index=False)
print("Save to " + inputfile[:-4]+'.addflags.cat')

for ff in glob.glob(PATH_table+"*.reg"):
    os.remove(ff)
for ff in glob.glob(PATH_table+"*_modelsub.fits"):
    os.remove(ff)
print("remove all middle stage files")
