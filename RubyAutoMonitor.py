import os
#from datetime import *
import time
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import tkinter.scrolledtext as ScrolledText
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (	FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
#from PIL import ImageTk,Image
import numpy as np
import scipy
import scipy.constants
from scipy.optimize import curve_fit
import threading
import re



class Canvas_for_plots():
	def __init__(self,Frame,**kwargs):
		self._frame = Frame
		fig = plt.figure(figsize=(kwargs.get('width',4),kwargs.get('height',4)),constrained_layout=True)
		spec2 = gridspec.GridSpec(ncols=(kwargs.get('ncols',1)), nrows=(kwargs.get('nrows',1)), figure=fig)
		num_subplots = kwargs.get('ncols',1)*kwargs.get('nrows',1)
		axs = []
		[i,j] = [0,0]
		for n in range(num_subplots):
			[i,j] = divmod(n, kwargs.get('ncols',1))
			axs.append(fig.add_subplot(spec2[i, j]))
		self._fig = fig
		self._axs = axs

	def initialise(self):		
		canvas = FigureCanvasTkAgg(self._fig, master=self._frame)  # A tk.DrawingArea.
		#toolbar = NavigationToolbar2Tk(canvas, self._frame)
		#toolbar.update()
		#canvas.draw()
		canvas.get_tk_widget().pack(fill=BOTH, expand=0)	
		
class Curve():
    def __init__(self,Canvas,index,**kwargs):
        self._Canvas = Canvas
        self._ax = Canvas._axs[index]
        self._xData = np.empty(0)
        self._yData = np.empty(0)
        '''
        if kwargs.get('linestyle'):
            line, = self._ax.plot(self._xData, self._yData, kwargs.get('linestyle'))
        else:
            line, = self._ax.plot(self._xData, self._yData)
        ''' 
        line, = self._ax.plot(self._xData, self._yData, **kwargs)
        self._curve = line

    def clear(self):
    #self._ax.cla()
        self._xData = np.empty(0)
        self._yData = np.empty(0)
			
	

def FlushPlot(Curve_object, **kwargs):
    if Curve_object._xData.size > 0:
        Curve_object._curve.set_xdata(Curve_object._xData)
        Curve_object._curve.set_ydata(Curve_object._yData)
        if kwargs.get('x_range') is not None:
            Curve_object._ax.set_xlim(0.995*kwargs.get('x_range')[0],1.005*kwargs.get('x_range')[1])
        else:
            Curve_object._ax.set_xlim(0.99*np.min(Curve_object._xData),1.01*np.max(Curve_object._xData))
        if kwargs.get('y_range') is not None:
            Curve_object._ax.set_ylim(0.995*kwargs.get('y_range')[0],1.005*kwargs.get('y_range')[1])
        else:
            Curve_object._ax.set_ylim(0.99*np.min(Curve_object._yData),1.01*np.max(Curve_object._yData))
        Curve_object._Canvas._fig.canvas.flush_events()
    else:
        pass
				
def FlushMultiPlots(Curve_object_list):
	x_low = []
	x_high = []
	y_low = []
	y_high = []
	for Curve_object in Curve_object_list:
		#print("x data =: ", Curve_object._xData)
		#print("y data =: ", Curve_object._yData)
		Curve_object._curve.set_xdata(Curve_object._xData)
		x_low.append(np.min(Curve_object._xData))
		x_high.append(np.max(Curve_object._xData))
		Curve_object._curve.set_ydata(Curve_object._yData)
		y_low.append(np.min(Curve_object._yData))
		y_high.append(np.max(Curve_object._yData))
		Curve_object._Canvas._fig.canvas.flush_events()
	xmin = np.min(x_low)
	xmax = np.max(x_high)
	ymin = np.min(y_low)
	ymax = np.max(y_high)
	Curve_object_list[0]._ax.set_xlim(0.99*xmin,1.01*xmax)
	Curve_object_list[0]._ax.set_ylim(0.99*ymin,1.01*ymax)

def AddPlot(Curve_object,**kwargs):
    Curve_object._curve.set_xdata(Curve_object._xData)
    Curve_object._curve.set_ydata(Curve_object._yData)
    if kwargs.get('x_range') is not None:
        Curve_object._ax.set_xlim(0.995*kwargs.get('x_range')[0],1.005*kwargs.get('x_range')[1])
    if kwargs.get('y_range') is not None:
        Curve_object._ax.set_ylim(0.995*kwargs.get('y_range')[0],1.005*kwargs.get('y_range')[1])
	
def clear_frame(frame):
   for widgets in frame.winfo_children():
      widgets.destroy()
    
    
    

def RubyR1Equ(R1peak_wavelength, wavelength_amb):
    pressure = 1904/7.665*(((1+(R1peak_wavelength/wavelength_amb-1))**7.665)-1)
    return pressure
    
def GaussianEqu(x, A, x0, w, c):
    return A*np.exp(-((x-x0)/w)**2) + c
        
def GaussianFit(slt_Temp, slt_r1, **kwargs):
    initGuess = kwargs.get('p0')
    if kwargs.get('p0') is not None:
        popt, pcov = curve_fit(GaussianEqu, slt_Temp, slt_r1, 
                                    p0=(initGuess)
                                #, bounds=([0,-1e3,0,0], [1e4, 100, 100, 100])
                                )
    else:
        popt, pcov = curve_fit(GaussianEqu, slt_Temp, slt_r1
                                #, bounds=([0,-1e3,0,0], [1e4, 100, 100, 100])
                                )

    return popt
   

   
def tableFrameInit():
    tableHeader = ['File No.', 'Torque (cNm)', 'Wavelength (nm)', 'Pressure (GPa)']
    for i in range(len(tableHeader)):
            Label(tableFrame, width=20, text=tableHeader[i],
                        #relief=SUNKEN
                        ).grid(row=1, column=i)
    
def readFitPlot(folderDirectory, rubyCanvas, tableFrame, torPCurve, **kwargs):
    global currentFileNo, tblHeader, tbl

    dir_list = os.listdir(folderDirectory)
    files = [os.path.join(folderDirectory, i) for i in dir_list]
    #print(files.sort(key=os.path.getctime))
    dir_list = [x for x in dir_list if re.search('.*(_HRD).*(\.txt)$', x)]
    currentFileNo = len(dir_list)
    
    files = [os.path.join(folderDirectory, x) for x in dir_list]
    files.sort(key=os.path.getctime)
    fileIndex = 0
    rawSptra = []
    R1peakSptra = []
    fitSptra = []
    
    wavelength_amb = kwargs.get('wavelength_amb', 694.25)
    wlamb_e.delete(0, END)
    wlamb_e.insert(END, wavelength_amb)
    
    wavelength_low = []
    wavelength_high = []
    intensity_low = []
    intensity_high = []
    
    #tableFrameInit()
    tblHeader = ['File No.', 'Torque (cNm)', 'Wavelength (nm)', 'Pressure (GPa)']
    #tbl = []
    tbl = np.empty([0,4])

    for filepath in files:
        rawSptra.append(Curve(rubyCanvas,0, marker='.'))
        R1peakSptra.append(Curve(rubyCanvas,0, linestyle='-', linewidth=3))
        fitSptra.append(Curve(rubyCanvas, 0, linestyle='-', color='r', linewidth=2))
        #print(dir)
        
        #torque = float(re.split('_', re.search('(.*)cNm', dir).group(1))[-1])
        torque = float(re.split('_', re.search('(.*)cNm', filepath).group(1))[-1])
        
        #filepath = os.path.join(folderDirectory,dir)
        content = np.genfromtxt(filepath,delimiter='\t',dtype=float,skip_header=14)
        content = np.transpose(content)
        wavelength = content[0]
        intensity = content[1]
        
        
        
        rawSptra[fileIndex]._xData = wavelength
        rawSptra[fileIndex]._yData = intensity


        ###----- fitting ----###
        halfwidth = kwargs.get('halfwidth', 15)
        peak_i = np.argmax(intensity)
        #width_nm = wavelength[peak_i+peakwidth/2] - wavelength[peak_i-peakwidth/2]
        
        R1_peak_x = wavelength[peak_i-halfwidth:peak_i+halfwidth]
        R1_peak_y = intensity[peak_i-halfwidth:peak_i+halfwidth]
        
        R1peakSptra[fileIndex]._xData = R1_peak_x
        R1peakSptra[fileIndex]._yData = R1_peak_y
        
        #plt.plot(R1_peak_x, R1_peak_y, 'r-', label='R1 Peak')
        
        popt = GaussianFit(R1_peak_x, R1_peak_y, p0 = [np.amax(intensity),
                                                        wavelength[peak_i],
                                                        wavelength[peak_i+halfwidth]-wavelength[peak_i-halfwidth],
                                                        np.amin(intensity)
                                                        ])
        
        A = popt[0]
        x0 = popt[1]
        w = popt[2]
        c = popt[3]
        
        
        pressure = RubyR1Equ(x0, wavelength_amb)
        
        
        
        tblRow = [fileIndex, torque, x0, pressure]
        tbl = np.vstack([tbl, tblRow])
 
        
        ####---- Simulate fitting curve----####
        sim_x = np.linspace(np.amin(R1_peak_x),np.amax(R1_peak_x),101)
        sim_y = []
        
        for i in range(len(sim_x)):
            sim_y = np.append(sim_y, GaussianEqu(sim_x[i], A, x0, w, c))
    
        fitSptra[fileIndex]._xData = sim_x
        fitSptra[fileIndex]._yData = sim_y
        
        
        wavelength_low.append(np.min(R1_peak_x))
        wavelength_high.append(np.max(R1_peak_x))
        intensity_low.append(np.min(R1_peak_y))
        intensity_high.append(np.max(R1_peak_y))
        x_min = np.min(wavelength_low)
        x_max = np.max(wavelength_high)
        y_min = np.min(intensity_low)
        y_max = np.max(intensity_high)
        
        AddPlot(rawSptra[fileIndex], x_range=[x_min, x_max])
        AddPlot(R1peakSptra[fileIndex], x_range=[x_min, x_max], y_range=[y_min, y_max])
        AddPlot(fitSptra[fileIndex], x_range=[x_min, x_max])
        #FlushPlot(fitSptra[fileIndex])
        fileIndex += 1
        #print('through')
        
    ####=== Plot pressure table ===####
    for i in range(len(tblHeader)):
        Label(tableFrame,text=tblHeader[i], width=15).grid(row=1,column=i)
    for i in range(len(tblRow)):
        for j in range(fileIndex):
            Label(tableFrame,text=round(tbl[j,i],3), width=15, relief=SUNKEN).grid(row=j+2,column=i)
        
    ####=== Plot Torque vs Pressure ===####
    torArray = np.asarray([float(i) for i in tbl[1:,1]])
    presuArray = np.asarray([float(i) for i in tbl[1:,3]])
    torPCurve._xData = torArray
    torPCurve._yData = presuArray
    FlushPlot(torPCurve, x_range=[0,40])
        
def selectdir():
    #global directory
    filedirectory = filedialog.askdirectory()
    fd_e.delete(0, 'end')
    fd_e.insert(END, filedirectory)
    
def selectRbFilesDir():
    fd_buffer = fd_e.get()
    selectdir()
    if fd_e.get() != fd_buffer:
        #print('change fd passed')
        rubyCurve._ax.cla()
        clear_frame(tableSubFrame)
        readFitPlot(fd_e.get(), rubyCanvas, tableSubFrame, torPCurve)
    
def saveRubyTable(): 
    header = np.expand_dims(tblHeader,axis=0)
    unit_row = ['', 'cNm', 'nm', 'GPa']
    unit_row = np.expand_dims(unit_row,axis=0)
    cm_str = ['%s' %'' for j in range(len(unit_row))]
    cm_row = np.expand_dims(cm_str,axis=0)
    
    
    fname = 'RubyR1Fit'
    fname_prefix = fname
    i = 1
    #fpath = os.path.join(current_path,fname)
    fpath = os.path.join(fd_e.get(),fname)
    while os.path.exists(fpath+'.txt'):
        fname = fname_prefix + '-%03d'%i
        fpath = os.path.join(fd_e.get(),fname)
        i += 1
    f = open(fpath+'.txt','w')
    
    
    np.savetxt(f, header, fmt='%5s' ,delimiter='\t',comments='')
    np.savetxt(f, unit_row, fmt='%5s' ,delimiter='\t',comments='')
    np.savetxt(f, cm_row, fmt='%5s' ,delimiter='\t',comments='')
    np.savetxt(f, tbl, fmt='%5s' ,delimiter='\t',comments='')
    f.close()

        
def RubyFilesDirUpdate():
    ct = threading.currentThread()
    while True:
        time.sleep(1)
        folderDirectory = fd_e.get()
        dir_list = os.listdir(folderDirectory)
        dir_list = [x for x in dir_list if re.search('.*(_HRD).*(\.txt)$', x)]
        nowFileNo = len(dir_list)
        #print(nowFileNo != currentFileNo)
        if nowFileNo != currentFileNo:
            rubyCurve._ax.cla()
            clear_frame(tableSubFrame)
            readFitPlot(fd_e.get(), rubyCanvas, tableSubFrame, torPCurve, wavelength_amb=float(wlamb_e.get()))
    
 
 
def convertorRuby(wavelength_amb, **kwargs):
    if kwargs.get('mode') == 'W2P':
        pressure = RubyR1Equ(float(wl_e.get()), wavelength_amb)
        p_e.delete(0, 'end')
        p_e.insert(END, pressure)
    elif kwargs.get('mode') == 'P2W':
        x0 = ((float(p_e.get())*(7.665/1904)+1)**(1/7.665))*(wavelength_amb)
        wl_e.delete(0, 'end')
        wl_e.insert(END, x0)
    else:
        ValueError('No such mode.')
 
def dummy():
    pass
 
 
if __name__ == "__main__":
    current_path = os.path.dirname(__file__)
    root = Tk()
    root.geometry("1280x720")
    #root.geometry("1440x720")
    root.title('Ruby tracker')
    plt.ion()						# enable automatic update plots from plt 


    DirectoryFrame = LabelFrame(root,text="Ruby files directory", width=100)
    DirectoryFrame.grid(row=0,column=0)

    fd_na = Label(DirectoryFrame, text = 'File directory:')
    fd_na.grid(row=3,column=0)
    fd_e = Entry(DirectoryFrame,width=50,borderwidth=2)
    #fd_e.insert(END, current_path)
    fd_e.insert(END, r'\\S4\Datenpool\Yuk Tai\Data and Analysis\Beryl-D2O_wP\BerD-S3\BerD-S3_Ruby')
    fd_e.grid(row=4,column=0,columnspan=10)
    #setattr(sample,'fd_entry',fd_e)
    fd_bn = Button(DirectoryFrame, text = "Select", command=selectRbFilesDir)
    fd_bn.grid(row=3,column=1)
    

    tableFrame = LabelFrame(root,text="Pressure table")
    tableFrame.grid(row=1,column=0,rowspan=10,columnspan=2)
    wlamb_na = Label(tableFrame,text="Wavelength_amb (nm)")
    wlamb_na.grid(row=0,column=0)
    wlamb_e = Entry(tableFrame,width=20,borderwidth=2)
    wlamb_e.insert(END, 694.25)
    wlamb_e.grid(row=0,column=1)
    recal_btn = Button(tableFrame, text='Recal', 
                        command=lambda: readFitPlot(fd_e.get(), rubyCanvas, 
                                                    tableSubFrame, torPCurve,
                                                    wavelength_amb=float(wlamb_e.get())))
    recal_btn.grid(row=0, column=2)
    save_btn = Button(tableFrame, text='Save', command=lambda: saveRubyTable())
    save_btn.grid(row=0, column=20)
    tableSubFrame = LabelFrame(tableFrame)
    tableSubFrame.grid(row=2,column=0,columnspan=200)



    calcuFrame = LabelFrame(root,text="Wavelength / Pressure convertor")
    calcuFrame.grid(row=0,column=1)
    Label(calcuFrame, text = 'Wavelength (nm)').grid(row=3,column=3)
    Label(calcuFrame, text = 'Pressure (GPa)').grid(row=4,column=3)
    wl_e = Entry(calcuFrame, width=20, borderwidth=2)
    wl_e.grid(row=3,column=4)
    p_e = Entry(calcuFrame, width=20, borderwidth=2)
    p_e.grid(row=4,column=4)
    cal_btn_wl = Button(calcuFrame, text = "Cal W2P", command=lambda: convertorRuby(694.25, mode='W2P'))
    cal_btn_wl.grid(row=3,column=5)
    cal_btn_p = Button(calcuFrame, text = "Cal P2W", command=lambda: convertorRuby(694.25, mode='P2W'))
    cal_btn_p.grid(row=4,column=5)
    
    
    rubySptraFrame = LabelFrame(root,text="Sprectra")
    rubySptraFrame.grid(row=0,column=2,rowspan=2,columnspan=1)
    rubyCanvas = Canvas_for_plots(rubySptraFrame,width=5,height=3,ncols=1,nrows=1)
    rubyCanvas.initialise()
    rubyCurve = Curve(rubyCanvas,0)
    rubyCurve._ax.set_xlabel('Wavelength (nm)')
    rubyCurve._ax.set_ylabel('Intensity ()')
    
    torquePressureRelatFrame = LabelFrame(root,text="Torque-Pressure")
    torquePressureRelatFrame.grid(row=2,column=2,rowspan=1,columnspan=1)
    torPCanvas = Canvas_for_plots(torquePressureRelatFrame,width=3,height=3,ncols=1,nrows=1)
    torPCanvas.initialise()
    torPCurve = Curve(torPCanvas,0, marker='o', markersize=5)
    torPCurve._ax.set_xlabel('Torque (cNm)')
    torPCurve._ax.set_ylabel('Pressure (GPa)')
    
    

    
    
    
    
    
    
    
    
    
    Button(root, text="Quit", command=quit).grid(row=0,column=10)
    Button(root, text="Quit", command=quit).grid(row=10,column=0)


    #########
    #DataArray_initialise()
    start_time = time.time()
    
    readFitPlot(fd_e.get(), rubyCanvas, tableSubFrame, torPCurve, wavelength_amb=float(wlamb_e.get()))
    
    threading.Thread(target=RubyFilesDirUpdate, args=()).start()

    root.mainloop()
    root.destroy()	#optional
