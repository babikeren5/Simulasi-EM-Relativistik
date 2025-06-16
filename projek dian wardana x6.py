# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 12:18:00 2023

@author: dian wardana
"""# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 09:39:35 2023

@author: dian wardana
"""


#lebel sumbu dan ganti x/t jadi x,t
#contoh label persamaan nonlinier E.t
#yang linier gak bagus

from tkinter import *
import tkinter as tk
import tkinter as ttk
import tkinter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
import xlsxwriter as xl
from tkinter import filedialog
from random import choice

class warna_ganti(tk.Button):
    def __init__(self, master, **kwargs):
        super().__init__(master,**kwargs)
        self.colors=['grey']
        self.bind("<Button-1>",self.change_color)
    def change_color(self,event):
        new_color=choice(self.colors)
        self.config(bg=new_color)
        
def save_data():
    #biasa
        q=aq.get() or 1
        B0=aB0.get() or 1
        B1=aB1.get() or 2
        m=am.get() or 4
        
        E0=aE0.get() or 1
        E1=aE1.get() or 1
        
        t0=at0.get() or 0
        tn=atn.get() or 10
        ndata=andata.get() or 1000
        
        x0=ax0.get() or 0
        y0=ay0.get() or 0
        z0=az0.get() or 0

        vx0=avx0.get() or 0
        vy0=avy0.get() or 0
        vz0=avz0.get() or 0
    #10^
        qT=aqT.get() or 0
        B0T=aB0T.get() or 0
        B1T=aB1T.get() or 0
        mT=amT.get() or 0
        
        E0T=aE0T.get() or 0
        E1T=aE1T.get() or 0
        
        t0T=at0T.get() or 0
        tnT=atnT.get() or 0
        ndataT=andataT.get() or 0
        
        x0T=ax0T.get() or 0
        y0T=ay0T.get() or 0
        z0T=az0T.get() or 0

        vx0T=avx0T.get() or 0
        vy0T=avy0T.get() or 0
        vz0T=avz0T.get() or 0
        
        file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
        if file_path:
            with open(file_path, "w") as file:
                file.write(f"{q}\n{B0}\n{B1}\n{m}\n{E0}\n{E1}\n{t0}\n{tn}\n{ndata}\n{x0}\n{y0}\n{z0}\n{vx0}\n{vy0}\n{vz0}\n{qT}\n{B0T}\n{B1T}\n{mT}\n{E0T}\n{E1T}\n{t0T}\n{tnT}\n{ndataT}\n{x0T}\n{y0T}\n{z0T}\n{vx0T}\n{vy0T}\n{vz0T}")
def load_data():
        file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt")])
        if file_path:
            with open(file_path, "r") as file:
                lines = file.readlines()
                if len(lines) >= 29:
                    aq.delete(0, tk.END)
                    aq.insert(0, lines[0].strip())
                    
                    aB0.delete(0, tk.END)
                    aB0.insert(0, lines[1].strip())
                    
                    aB1.delete(0, tk.END)
                    aB1.insert(0, lines[2].strip())
                    
                    am.delete(0, tk.END)
                    am.insert(0, lines[3].strip())
                    
                    aE0.delete(0, tk.END)
                    aE0.insert(0, lines[4].strip())
                    
                    aE1.delete(0, tk.END)
                    aE1.insert(0, lines[5].strip())
                    
                    at0.delete(0, tk.END)
                    at0.insert(0, lines[6].strip())
                    
                    atn.delete(0, tk.END)
                    atn.insert(0, lines[7].strip())
                    
                    andata.delete(0, tk.END)
                    andata.insert(0, lines[8].strip())
                    
                    ax0.delete(0, tk.END)
                    ax0.insert(0, lines[9].strip())
                    
                    ay0.delete(0, tk.END)
                    ay0.insert(0, lines[10].strip())
                    
                    az0.delete(0, tk.END)
                    az0.insert(0, lines[11].strip())

                    avx0.delete(0, tk.END)
                    avx0.insert(0, lines[12].strip())
                    
                    avy0.delete(0, tk.END)
                    avy0.insert(0, lines[13].strip())

                    avz0.delete(0, tk.END)
                    avz0.insert(0, lines[14].strip())
                    
                    #10^
                    
                    aqT.delete(0, tk.END)
                    aqT.insert(0, lines[15].strip())
                    
                    aB0T.delete(0, tk.END)
                    aB0T.insert(0, lines[16].strip())
                    
                    aB1T.delete(0, tk.END)
                    aB1T.insert(0, lines[17].strip())
                    
                    amT.delete(0, tk.END)
                    amT.insert(0, lines[18].strip())
                    
                    aE0T.delete(0, tk.END)
                    aE0T.insert(0, lines[19].strip())
                    
                    aE1T.delete(0, tk.END)
                    aE1T.insert(0, lines[20].strip())
                    
                    at0T.delete(0, tk.END)
                    at0T.insert(0, lines[21].strip())
                    
                    atnT.delete(0, tk.END)
                    atnT.insert(0, lines[22].strip())
                    
                    andataT.delete(0, tk.END)
                    andataT.insert(0, lines[23].strip())
                    
                    ax0T.delete(0, tk.END)
                    ax0T.insert(0, lines[24].strip())
                    
                    ay0T.delete(0, tk.END)
                    ay0T.insert(0, lines[25].strip())
                    
                    az0T.delete(0, tk.END)
                    az0T.insert(0, lines[26].strip())

                    avx0T.delete(0, tk.END)
                    avx0T.insert(0, lines[27].strip())
                    
                    avy0T.delete(0, tk.END)
                    avy0T.insert(0, lines[28].strip())

                    avz0T.delete(0, tk.END)
                    avz0T.insert(0, lines[29].strip())
def nonlinier(): 
    nonlinier = Tk()
    nonlinier.title("nonlinier ")
    nonlinier.geometry("200x200+0+480")

    qp=float(aq.get() or 1)
    mp=float(am.get() or 1)
    qT=float(aqT.get() or 0)
    mT=float(amT.get() or 0)
    q=qp*(10**qT)
    m=mp*(10**mT)
    
    k=q/m
    
    E0p=float(aE0.get() or 1)
    E1p=float(aE1.get() or 1)
    E0T=float(aE0T.get() or 0)
    E1T=float(aE1T.get() or 0)
    E0=E0p*(10**E0T)
    E1=E1p*(10**E1T)
    
    B0p=float(aB0.get() or 1)
    B1p=float(aB1.get() or 1)
    B0T=float(aB0T.get() or 0)
    B1T=float(aB1T.get() or 0)
    B0=B0p*(10**B0T)
    B1=B1p*(10**B1T)
    
    t0p=float(at0.get() or 0)
    tnp=float(atn.get() or 10)
    ndatap=int(andata.get() or 1000)
    t0T=float(at0T.get() or 0)
    tnT=float(atnT.get() or 0)
    ndataT=int(andataT.get() or 0)
    t0=t0p*(10**t0T)
    tn=tnp*(10**tnT)
    ndata=ndatap*(10**ndataT)
    t=np.linspace(t0,tn,ndata)
    
    c=3*(10**9)
    
    x0p=float(ax0.get() or 1)
    x0T=float(ax0T.get() or 0)
    x0=x0p*(10**x0T)
    x=np.zeros(ndata)
    x[0]=x0
    
    
    y0p=float(ay0.get() or 1)
    y0T=float(ay0T.get() or 0)
    y0=y0p*(10**y0T)
    y=np.zeros(ndata)
    y[0]=y0
    
    z0p=float(az0.get() or 1)
    z0T=float(az0T.get() or 0)
    z0=z0p*(10**z0T)
    z=np.zeros(ndata)
    z[0]=z0
    
    vx0p=float(avx0.get() or 1)
    vx0T=float(avx0T.get() or 0)
    vx0=vx0p*(10**vx0T)
    vx=np.zeros(ndata)
    vx[0]=vx0
    
    vy0p=float(avy0.get() or 1)
    vy0T=float(avy0T.get() or 0)
    vy0=vy0p*(10**vy0T)
    vy=np.zeros(ndata)
    vy[0]=vy0
    
    vz0p=float(avz0.get() or 1)
    vz0T=float(avz0T.get() or 0)
    vz0=vz0p*(10**vz0T)
    vz=np.zeros(ndata)
    vz[0]=vz0
    
    def rk4():
        yuhu= Tk()
        yuhu.title("grafik yang diinginkan")
        yuhu.geometry("240x470+1110+0")
        
        def funcy (t , vy , y , k, B0, B1, y0, vx0 ):
             return vy
        def guncy (t , vy , y , k, B0, B1, y0, vx0 ):
             return -(k**2) * ((B0+B1*t)**2) *(y-y0)-((k) * ((B0+B1*t)) *vx0)
        def funcx (t , vx , x , k, B0, B1, x0, vy0  ):
             return vx
        def guncx (t , vx , x , k, B0, B1, x0, vy0  ):
             return ((k**2) * ((B0+B1*t)**2) )*(x0-x)+((k) * ((B0+B1*t)) *vy0)
        def funcz (t , vz , z , k, E0, E1 ):
             return vz
        def guncz (t , vz , z , k, E0, E1 ):
             return k*(E0+E1*t)
        for i in range (ndata-1):
            
            h=t[i+1]-t[i]
            k1y=h*funcy(t[i],vy[i],y[i] ,k, B0, B1, y0, vx0  )
            l1y=h*guncy(t[i],vy[i],y[i] ,k, B0, B1, y0, vx0  )
            k2y=h*funcy(t[i]+0.5*h,vy[i]+(0.5*l1y),y[i]+(0.5*k1y) ,k, B0, B1, y0, vx0  )
            l2y=h*guncy(t[i]+0.5*h,vy[i]+(0.5*l1y),y[i]+(0.5*k1y) ,k, B0, B1, y0, vx0  )
            k3y=h*funcy(t[i]+0.5*h,vy[i]+(0.5*l2y),y[i]+(0.5*k2y) ,k, B0, B1, y0, vx0  )
            l3y=h*guncy(t[i]+0.5*h,vy[i]+(0.5*l2y),y[i]+(0.5*k2y) ,k, B0, B1, y0, vx0  )
            k4y=h*funcy(t[i]+h,vy[i]+l3y,y[i]+k3y ,k, B0, B1, y0, vx0  )
            l4y=h*guncy(t[i]+h,vy[i]+l3y,y[i]+k3y ,k, B0, B1, y0, vx0  )
            
            k1x=h*funcx(t[i],vx[i],x[i] ,k, B0, B1, x0, vy0  )
            l1x=h*guncx(t[i],vx[i],x[i] ,k, B0, B1, x0, vy0  )
            k2x=h*funcx(t[i]+0.5*h,vx[i]+(0.5*l1x),x[i]+(0.5*k1x) ,k, B0, B1, x0, vy0  )
            l2x=h*guncx(t[i]+0.5*h,vx[i]+(0.5*l1x),x[i]+(0.5*k1x) ,k, B0, B1, x0, vy0  )
            k3x=h*funcx(t[i]+0.5*h,vx[i]+(0.5*l2x),x[i]+(0.5*k2x) ,k, B0, B1, x0, vy0  )
            l3x=h*guncx(t[i]+0.5*h,vx[i]+(0.5*l2x),x[i]+(0.5*k2x) ,k, B0, B1, x0, vy0  )
            k4x=h*funcx(t[i]+h,vx[i]+l3x,x[i]+k3x ,k, B0, B1, x0, vy0  )
            l4x=h*guncx(t[i]+h,vx[i]+l3x,x[i]+k3x ,k, B0, B1, x0, vy0  )
            
            k1z=h*funcz(t[i]  , vz[i] , z[i], k, E0, E1) 
            l1z=h*guncz(t[i]  , vz[i] , z[i], k, E0, E1)
            k2z=h*funcz(t[i]+0.5*h , vz[i]+(0.5*l1z) , z[i]+(0.5*k1z), k, E0, E1)
            l2z=h*guncz(t[i]+0.5*h , vz[i]+(0.5*l1z) , z[i]+(0.5*k1z), k, E0, E1)
            k3z=h*funcz(t[i]+0.5*h , vz[i]+(0.5*l2z) , z[i]+(0.5*k2z), k, E0, E1)
            l3z=h*guncz(t[i]+0.5*h , vz[i]+(0.5*l2z) , z[i]+(0.5*k2z), k, E0, E1)
            k4z=h*funcz(t[i]+h , vz[i]+l3z , z[i]+k3z , k, E0, E1)
            l4z=h*guncz(t[i]+h , vz[i]+l3z , z[i]+k3z , k, E0, E1)
            
            kky=(k1y+2*k2y+2*k3y+k4y)/6.0
            ly=(l1y+2*l2y+2*l3y+l4y)/6.0
            
            y[i+1]=y[i]+kky
            vy[i+1]=vy[i]+ly
            
            kkx=(k1x+2*k2x+2*k3x+k4x)/6.0
            lx=(l1x+2*l2x+2*l3x+l4x)/6.0

            x[i+1]=x[i]+kkx
            vx[i+1]=vx[i]+lx

            kkz=(k1z+2*k2z+2*k3z+k4z)/6.0
            lz=(l1z+2*l2z+2*l3z+l4z)/6.0
            
            z[i+1]=z[i]+kkz
            vz[i+1]=vz[i]+lz
            
            v=np.sqrt(vx**2+vy**2+vz**2)
            g=np.sqrt(1-(v**2/c**2))
            
            M=m/g
            E=M*(c**2)
            KI=(c**2)*M-(c**2)*m
            
            x_max=np.max(x)
            x_min=np.min(x)
            
            y_max=np.max(y)
            y_min=np.min(y)
            
            z_max=np.max(z)
            z_min=np.min(z)
            
            v_max=np.max(v)
            v_min=np.min(v)
            
            M_max=np.max(M)
            M_min=np.min(M)
            
            E_max=np.max(E)
            E_min=np.min(E)
            
            KI_max=np.max(KI)
            KI_min=np.min(KI)
            
        def gambarx():
            root = tkinter.Tk()
            root.title('ploting posisi x, t pada rk4 ')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,x,'r',label="posisi x terhadap waktu")
            plot1.set_ylabel('posisi x(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambary():
            root = tkinter.Tk()
            root.title('ploting posisi y, t pada rk4 ')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,y,'b',label="posisi y terhadap waktu")
            plot1.set_ylabel('posisi y(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarz():
            root = tkinter.Tk()
            root.title('ploting posisi z, t pada rk4 ')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,z,'g',label="posisi z terhadap waktu")
            plot1.set_ylabel('posisi z(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarv():
            root = tkinter.Tk()
            root.title('ploting kecepatan mutlak dan waktu pada rk4')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,v,'y',label="kecepatan mutlak v terhadap waktu")
            plot1.set_ylabel('kecepatan V(m/s)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {v_min}', xy=(t[np.argmin(v)], v_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {v_max}', xy=(t[np.argmax(v)], v_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarm():
            root = tkinter.Tk()
            root.title('ploting massa dan waktu pada rk4')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,M,'g',label="massa M terhadap waktu")
            plot1.set_ylabel('massa gerak M(kg)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {M_min}', xy=(t[np.argmin(M)], M_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {M_max}', xy=(t[np.argmax(M)], M_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarKI():
            root = tkinter.Tk()
            root.title('ploting energi kinetik waktu pada rk4')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,KI,'g',label="energi kinetik terhadap waktu")
            plot1.set_ylabel('energi kinetik KI(j)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {KI_min}', xy=(t[np.argmin(KI)], KI_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {KI_max}', xy=(t[np.argmax(KI)], KI_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarE():
            root = tkinter.Tk()
            root.title('ploting energi pada rk4')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,E,'r',label="energi gerak")
            plot1.set_ylabel('energi E(j)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {E_min}', xy=(t[np.argmin(E)], E_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {E_max}', xy=(t[np.argmax(E)], E_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
                canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop()  
        def gambar2d():
            root = tkinter.Tk()
            root.title('ploting medan magnet pada rk4')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(x,y,'y',label="posisi x terhadap y")
            plot1.set_xlabel('posisi x')
            plot1.set_ylabel('posisi y')
            plot1.plot(x[0],y[0],marker='o',markersize=8,color='blue',label="stars")
            plot1.plot(x[-1],y[-1],marker='o',markersize=8,color='green',label="end")
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambar3d():
            root = tkinter.Tk()
            root.title('ploting medan elektromagnet pada rk4')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111,projection='3d')
            plot1.plot(x,y,z,'r',label="posisi")
            plot1.set_xlabel('posisi x(m)')
            plot1.set_ylabel('posisi y(m)')
            plot1.set_zlabel('posisi z(m)')
            plot1.plot(x[0],y[0],z[0],marker='o',markersize=8,color='blue',label="stars")
            plot1.plot(x[-1],y[-1],z[-1],marker='o',markersize=8,color='green',label="end")
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
                canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop()  
        def hasilkanxl():
            printa = Tk()
            printa.title("hasilprint pada rk4")
            printa.geometry("240x100+1110+520")
            lbl=Label(printa, text="nama file    .xlsx",width=15)
            lbl.grid(column=0,row=0)
            namafile=Entry(printa, width=10)
            namafile.grid(column=1,row=0)
            def printnya():
                workbook=xl.Workbook(namafile.get())
                worksheet=workbook.add_worksheet()
                worksheet.write(0,0,"q")
                worksheet.write(1,0,q)
                worksheet.write(0,1,"B0")
                worksheet.write(1,1,B0)
                worksheet.write(0,2,"B1")
                worksheet.write(1,2,B1)
                worksheet.write(0,3,"m")
                worksheet.write(1,3,m)
                worksheet.write(0,4,"E0")
                worksheet.write(1,4,E0)
                worksheet.write(0,5,"E1")
                worksheet.write(1,5,E1)
                worksheet.write(0,6,"t0")
                worksheet.write(1,6,t0)
                worksheet.write(0,7,"tn")
                worksheet.write(1,7,tn)
                worksheet.write(0,8,"ndata")
                worksheet.write(1,8,ndata)
                worksheet.write(0,9,"waktu")
                worksheet.write_column(1,9,t)
                worksheet.write(0,10,"x")
                worksheet.write_column(1,10,x)
                worksheet.write(0,11,"y")
                worksheet.write_column(1,11,y)
                worksheet.write(0,12,"z")
                worksheet.write_column(1,12,z)
                worksheet.write(0,13,"vx")
                worksheet.write_column(1,13,vx)
                worksheet.write(0,14,"vy")
                worksheet.write_column(1,14,vy)
                worksheet.write(0,15,"vz")
                worksheet.write_column(1,15,vz)
                worksheet.write(0,16,"v")
                worksheet.write_column(1,16,v)
                worksheet.write(0,17,"M")
                worksheet.write_column(1,17,M)
                worksheet.write(0,18,"E")
                worksheet.write_column(1,18,E)
                worksheet.write(0,19,"KI")
                worksheet.write_column(1,19,KI)
                
                worksheet.write(0,20,"x")
                worksheet.write(1,20,x_max)
                worksheet.write(2,20,x_min)
                worksheet.write(0,21,"y")
                worksheet.write(1,21,y_max)
                worksheet.write(2,21,y_min)
                worksheet.write(0,22,"z")
                worksheet.write(1,22,z_max)
                worksheet.write(2,22,z_min)
                worksheet.write(0,23,"v")
                worksheet.write(1,23,v_max)
                worksheet.write(2,23,v_min)
                worksheet.write(0,24,"M")
                worksheet.write(1,24,M_max)
                worksheet.write(2,24,M_min)
                worksheet.write(0,25,"E")
                worksheet.write(1,25,E_max)
                worksheet.write(2,25,E_min)
                worksheet.write(0,26,"KI")
                worksheet.write(1,26,KI_max)
                worksheet.write(2,26,KI_min)
                workbook.close()

            btn=warna_ganti(printa,text="printnya",command=printnya).grid(column=0,row=3)    
        
        lbl=Label(yuhu, text="RK4 ",width=15).grid(column=0,row=1)
        lbl=Label(yuhu, text="pilih output ",width=10).grid(column=0,row=2)
        lbl=Label(yuhu, text=" ",width=10).grid(column=0,row=3)
        lbl=Label(yuhu, text="posisi ",width=10).grid(column=0,row=4)

        lbl=Label(yuhu, text="x,t",width=10).grid(column=0,row=5)        
        btn=warna_ganti(yuhu,text="_",command=gambarx,bg="green",fg="white").grid(column=1,row=5)
        lbl=Label(yuhu, text="y,t",width=10).grid(column=0,row=6)
        btn=warna_ganti(yuhu,text="_",command=gambary,bg="green",fg="white").grid(column=1,row=6)
        lbl=Label(yuhu, text="z,t",width=10).grid(column=0,row=7)
        btn=warna_ganti(yuhu,text="_",command=gambarz,bg="green",fg="white").grid(column=1,row=7)
        
        lbl=Label(yuhu, text="kecepatan mutlak ",width=15).grid(column=0,row=8)
        lbl=Label(yuhu, text="v,t ",width=10).grid(column=0,row=9)
        btn=warna_ganti(yuhu,text="_",command=gambarv,bg="green",fg="white").grid(column=1,row=9)
        
        lbl=Label(yuhu, text="massa gerak ",width=10).grid(column=0,row=10)
        lbl=Label(yuhu, text="M,t",width=10).grid(column=0,row=11)
        btn=warna_ganti(yuhu,text="_",command=gambarm,bg="green",fg="white").grid(column=1,row=11)
        
        lbl=Label(yuhu, text="energi gerak ",width=10).grid(column=0,row=12)
        lbl=Label(yuhu, text="E,t",width=10).grid(column=0,row=13)
        btn=warna_ganti(yuhu,text="_",command=gambarE,bg="green",fg="white").grid(column=1,row=13)
        
        lbl=Label(yuhu, text="energi kinetik ",width=10).grid(column=0,row=14)
        lbl=Label(yuhu, text="K,t",width=10).grid(column=0,row=15)
        btn=warna_ganti(yuhu,text="_",command=gambarKI,bg="green",fg="white").grid(column=1,row=15)
        
        lbl=Label(yuhu, text="Simulasi gerak ",width=10).grid(column=0,row=16)
        lbl=Label(yuhu, text="2D ",width=10).grid(column=0,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar2d,bg="green",fg="white").grid(column=1,row=17)
        
        lbl=Label(yuhu, text="3D ",width=10).grid(column=0,row=18)
        btn=warna_ganti(yuhu,text="_",command=gambar3d,bg="green",fg="white").grid(column=1,row=18)
        
        lbl=Label(yuhu, text="print xlsx ",width=10).grid(column=0,row=19)
        lbl=Label(yuhu, text="xlsx ",width=10).grid(column=0,row=20)
        btn=warna_ganti(yuhu,text="_",command=hasilkanxl,bg="green",fg="white").grid(column=1,row=20)
        yuhu.mainloop()
        
    def rk45B():
        yuhu= Tk()
        yuhu.title("grafik yang diinginkan")
        yuhu.geometry("240x470+1110+0")
        
        def funcy (t , vy , y , k, B0, B1, y0, vx0 ):
             return vy
        def guncy (t , vy , y , k, B0, B1, y0, vx0 ):
             return -(k**2) * ((B0+B1*t)**2) *(y-y0)-((k) * ((B0+B1*t)) *vx0)
        def funcx (t , vx , x , k, B0, B1, x0, vy0  ):
             return vx
        def guncx (t , vx , x , k, B0, B1, x0, vy0  ):
             return ((k**2) * ((B0+B1*t)**2) )*(x0-x)+((k) * ((B0+B1*t)) *vy0)
        def funcz (t , vz , z , k, E0, E1 ):
             return vz
        def guncz (t , vz , z , k, E0, E1 ):
             return k*(E0+E1*t)
        for i in range (ndata-1):

           h=t[i+1]-t[i]
            
           k1y=h*funcy(t[i] ,vy[i],y[i],k, B0, B1, y0, vx0 )
           l1y=h*guncy(t[i] ,vy[i],y[i],k, B0, B1, y0, vx0 )
           k2y=h*funcy(t[i]+(1*h/4) , vy[i]+(1*l1y/4),y[i]+(1*k1y/4),k, B0, B1, y0, vx0  )
           l2y=h*guncy(t[i]+(1*h/4) , vy[i]+(1*l1y/4),y[i]+(1*k1y/4),k, B0, B1, y0, vx0  )
           k3y=h*funcy(t[i]+(1*h/4) , vy[i]+(1*l1y/8)+(1*l2y/8) ,y[i]+(1*k1y/8)+(1*k2y/8),k, B0, B1, y0, vx0 )
           l3y=h*guncy(t[i]+(1*h/4) , vy[i]+(1*l1y/8)+(1*l2y/8) ,y[i]+(1*k1y/8)+(1*k2y/8),k, B0, B1, y0, vx0 )
           k4y=h*funcy(t[i]+(1*h/2) ,vy[i]-(1*l2y/2)+(l3y) ,y[i]-(1*k2y/2)+k3y,k, B0, B1, y0, vx0 )
           l4y=h*guncy(t[i]+(1*h/2) ,vy[i]-(1*l2y/2)+(l3y) ,y[i]-(1*k2y/2)+k3y,k, B0, B1, y0, vx0 )
           k5y=h*funcy(t[i]+(3*h/4) ,vy[i]+(3*l1y/16)+(9*l4y/16) ,y[i]+(3*k1y/16)+(9*k4y/16),k, B0, B1, y0, vx0 )
           l5y=h*guncy(t[i]+(3*h/4) ,vy[i]+(3*l1y/16)+(9*l4y/16) ,y[i]+(3*k1y/16)+(9*k4y/16),k, B0, B1, y0, vx0 )
           k6y=h*funcy(t[i]+h ,vy[i]-(3*l1y/7)+(2*l2y/7)+(12*l3y/7)-(12*l4y/7)+(8*l5y/7), y[i]-(3*k1y/7)+(2*k2y/7)+(12*k3y/7)-(12*k4y/7)+(8*k5y/7),k, B0, B1, y0, vx0 )
           l6y=h*guncy(t[i]+h ,vy[i]-(3*l1y/7)+(2*l2y/7)+(12*l3y/7)-(12*l4y/7)+(8*l5y/7), y[i]-(3*k1y/7)+(2*k2y/7)+(12*k3y/7)-(12*k4y/7)+(8*k5y/7),k, B0, B1, y0, vx0 ) 
              
               
           k1x=h*funcx(t[i] ,vx[i] ,x[i] ,k, B0, B1, x0, vy0 )
           l1x=h*guncx(t[i] ,vx[i] ,x[i] ,k, B0, B1, x0, vy0 )
           k2x=h*funcx(t[i]+(1*h/4) , vx[i]+(1*l1x/4) ,x[i]+(1*k1x/4),k, B0, B1, x0, vy0 )
           l2x=h*guncx(t[i]+(1*h/4) , vx[i]+(1*l1x/4) ,x[i]+(1*k1x/4),k, B0, B1, x0, vy0 )
           k3x=h*funcx(t[i]+(1*h/4) , vx[i]+(1*l1x/8)+(1*l2x/8) ,x[i]+(1*k1x/8)+(1*k2x/8),k, B0, B1, x0, vy0)
           l3x=h*guncx(t[i]+(1*h/4) , vx[i]+(1*l1x/8)+(1*l2x/8) ,x[i]+(1*k1x/8)+(1*k2x/8),k, B0, B1, x0, vy0)
           k4x=h*funcx(t[i]+(1*h/2) ,vx[i]-(1*l2x/2)+(l3x) ,x[i]-(1*k2x/2)+k3x,k, B0, B1, x0, vy0)
           l4x=h*guncx(t[i]+(1*h/2) ,vx[i]-(1*l2x/2)+(l3x) ,x[i]-(1*k2x/2)+k3x,k, B0, B1, x0, vy0)
           k5x=h*funcx(t[i]+(3*h/4) ,vx[i]+(3*l1x/16)+(9*l4x/16) ,x[i]+(3*k1x/16)+(9*k4x/16),k, B0, B1, x0, vy0)
           l5x=h*guncx(t[i]+(3*h/4) ,vx[i]+(3*l1x/16)+(9*l4x/16) ,x[i]+(3*k1x/16)+(9*k4x/16),k, B0, B1, x0, vy0)
           k6x=h*funcx(t[i]+h ,vx[i]-(3*l1x/7)+(2*l2x/7)+(12*l3x/7)-(12*l4x/7)+(8*l5x/7), x[i]-(3*k1x/7)+(2*k2x/7)+(12*k3x/7)-(12*k4x/7)+(8*k5x/7),k, B0, B1, x0, vy0)
           l6x=h*guncx(t[i]+h ,vx[i]-(3*l1x/7)+(2*l2x/7)+(12*l3x/7)-(12*l4x/7)+(8*l5x/7), x[i]-(3*k1x/7)+(2*k2x/7)+(12*k3x/7)-(12*k4x/7)+(8*k5x/7),k, B0, B1, x0, vy0) 

           
           k1z=h*funcz(t[i] ,vz[i],z[i],k, E0, E1)
           l1z=h*guncz(t[i] ,vz[i],z[i],k, E0, E1 )
           k2z=h*funcz(t[i]+(1*h/4) , vz[i]+(1*l1z/4),z[i]+(1*k1z/4),k, E0, E1 )
           l2z=h*guncz(t[i]+(1*h/4) , vz[i]+(1*l1z/4),z[i]+(1*k1z/4),k, E0, E1 )
           k3z=h*funcz(t[i]+(1*h/4) , vz[i]+(1*l1z/8)+(1*l2z/8) ,z[i]+(1*k1z/8)+(1*k2z/8),k, E0, E1)
           l3z=h*guncz(t[i]+(1*h/4) , vz[i]+(1*l1z/8)+(1*l2z/8) ,z[i]+(1*k1z/8)+(1*k2z/8),k, E0, E1)
           k4z=h*funcz(t[i]+(1*h/2) ,vz[i]-(1*l2z/2)+(l3z) ,z[i]-(1*k2z/2)+k3z,k, E0, E1)
           l4z=h*guncz(t[i]+(1*h/2) ,vz[i]-(1*l2z/2)+(l3z) ,z[i]-(1*k2z/2)+k3z,k, E0, E1)
           k5z=h*funcz(t[i]+(3*h/4) ,vz[i]+(3*l1z/16)+(9*l4z/16) ,z[i]+(3*k1z/16)+(9*k4z/16),k, E0, E1)
           l5z=h*guncz(t[i]+(3*h/4) ,vz[i]+(3*l1z/16)+(9*l4z/16) ,z[i]+(3*k1z/16)+(9*k4z/16),k, E0, E1)
           k6z=h*funcz(t[i]+h ,vz[i]-(3*l1z/7)+(2*l2z/7)+(12*l3z/7)-(12*l4z/7)+(8*l5z/7), z[i]-(3*k1z/7)+(2*k2z/7)+(12*k3z/7)-(12*k4z/7)+(8*k5z/7),k, E0, E1)
           l6z=h*guncz(t[i]+h ,vz[i]-(3*l1z/7)+(2*l2z/7)+(12*l3z/7)-(12*l4z/7)+(8*l5z/7), z[i]-(3*k1z/7)+(2*k2z/7)+(12*k3z/7)-(12*k4z/7)+(8*k5z/7),k, E0, E1)
             
              
           kky=(7*k1y+32*k3y+12*k4y+32*k5y+7*k6y)/90
           ly=(7*l1y+32*l3y+12*l4y+32*l5y+7*l6y)/90
              
           y[i+1]=y[i]+kky
           vy[i+1]=vy[i]+ly
              
           kkx=(7*k1x+32*k3x+12*k4x+32*k5x+7*k6x)/90
           lx=(7*l1x+32*l3x+12*l4x+32*l5x+7*l6x)/90
          
           x[i+1]=x[i]+kkx
           vx[i+1]=vx[i]+lx
          
           kkz=(7*k1z+32*k3z+12*k4z+32*k5z+7*k6z)/90
           lz=(7*l1z+32*l3z+12*l4z+32*l5z+7*l6z)/90
              
           z[i+1]=z[i]+kkz
           vz[i+1]=vz[i]+lz       
           v=np.sqrt(vx**2+vy**2+vz**2)
           
           g=np.sqrt(1-(v**2/c**2))
           
           M=m/g
           E=M*(c**2)
           KI=(c**2)*M-(c**2)*m        
           
           x_max=np.max(x)
           x_min=np.min(x)
            
           y_max=np.max(y)
           y_min=np.min(y)
            
           z_max=np.max(z)
           z_min=np.min(z)
            
           v_max=np.max(v)
           v_min=np.min(v)
           
           M_max=np.max(M)
           M_min=np.min(M)
           
           E_max=np.max(E)
           E_min=np.min(E)
           
           KI_max=np.max(KI)
           KI_min=np.min(KI)
           
        def gambarx():
            root = tkinter.Tk()
            root.title('ploting posisi x, t pada RK45B ')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,x,'r',label="posisi x terhadap waktu")
            plot1.set_ylabel('posisi x(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambary():
            root = tkinter.Tk()
            root.title('ploting posisi y, t pada RK45B ')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,y,'b',label="posisi y terhadap waktu")
            plot1.set_ylabel('posisi y(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarz():
            root = tkinter.Tk()
            root.title('ploting posisi z, t pada RK45B ')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,z,'g',label="posisi z terhadap waktu")
            plot1.set_ylabel('posisi z(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarv():
            root = tkinter.Tk()
            root.title('ploting kecepatan mutlak dan waktu pada RK45B')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,v,'y',label="kecepatan mutlak v terhadap waktu")
            plot1.set_ylabel('kecepatan V(m/s)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {v_min}', xy=(t[np.argmin(v)], v_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {v_max}', xy=(t[np.argmax(v)], v_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarm():
            root = tkinter.Tk()
            root.title('ploting massa dan waktu pada RK45B')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,M,'g',label="massa M terhadap waktu")
            plot1.set_ylabel('massa gerak M(kg)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {M_min}', xy=(t[np.argmin(M)], M_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {M_max}', xy=(t[np.argmax(M)], M_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarKI():
            root = tkinter.Tk()
            root.title('ploting energi kinetik waktu pada RK45B')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,KI,'g',label="energi kinetik terhadap waktu")
            plot1.set_ylabel('energi kinetik KI(j)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {KI_min}', xy=(t[np.argmin(KI)], KI_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {KI_max}', xy=(t[np.argmax(KI)], KI_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarE():
            root = tkinter.Tk()
            root.title('ploting energi pada RK45B')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,E,'r',label="energi gerak")
            plot1.set_ylabel('energi E(j)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {E_min}', xy=(t[np.argmin(E)], E_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {E_max}', xy=(t[np.argmax(E)], E_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
                canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop()  
        def gambar2d():
            root = tkinter.Tk()
            root.title('ploting medan magnet pada RK45B')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(x,y,'y',label="posisi x terhadap y")
            plot1.set_xlabel('posisi x')
            plot1.set_ylabel('posisi y')
            plot1.plot(x[0],y[0],marker='o',markersize=8,color='blue',label="stars")
            plot1.plot(x[-1],y[-1],marker='o',markersize=8,color='green',label="end")
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambar3d():
            root = tkinter.Tk()
            root.title('ploting medan elektromagnet pada RK445B')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111,projection='3d')
            plot1.plot(x,y,z,'r',label="posisi")
            plot1.set_xlabel('posisi x(m)')
            plot1.set_ylabel('posisi y(m)')
            plot1.set_zlabel('posisi z(m)')
            plot1.plot(x[0],y[0],z[0],marker='o',markersize=8,color='blue',label="stars")
            plot1.plot(x[-1],y[-1],z[-1],marker='o',markersize=8,color='green',label="end")
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
                canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop()  
        def hasilkanxl():
            printa = Tk()
            printa.title("hasilprint pada RK45B")
            printa.geometry("240x100+1110+520")
            lbl=Label(printa, text="nama file    .xlsx",width=15)
            lbl.grid(column=0,row=0)
            namafile=Entry(printa, width=10)
            namafile.grid(column=1,row=0)
            def printnya():
                workbook=xl.Workbook(namafile.get())
                worksheet=workbook.add_worksheet()
                worksheet.write(0,0,"q")
                worksheet.write(1,0,q)
                worksheet.write(0,1,"B0")
                worksheet.write(1,1,B0)
                worksheet.write(0,2,"B1")
                worksheet.write(1,2,B1)
                worksheet.write(0,3,"m")
                worksheet.write(1,3,m)
                worksheet.write(0,4,"E0")
                worksheet.write(1,4,E0)
                worksheet.write(0,5,"E1")
                worksheet.write(1,5,E1)
                worksheet.write(0,6,"t0")
                worksheet.write(1,6,t0)
                worksheet.write(0,7,"tn")
                worksheet.write(1,7,tn)
                worksheet.write(0,8,"ndata")
                worksheet.write(1,8,ndata)
                worksheet.write(0,9,"waktu")
                worksheet.write_column(1,9,t)
                worksheet.write(0,10,"x")
                worksheet.write_column(1,10,x)
                worksheet.write(0,11,"y")
                worksheet.write_column(1,11,y)
                worksheet.write(0,12,"z")
                worksheet.write_column(1,12,z)
                worksheet.write(0,13,"vx")
                worksheet.write_column(1,13,vx)
                worksheet.write(0,14,"vy")
                worksheet.write_column(1,14,vy)
                worksheet.write(0,15,"vz")
                worksheet.write_column(1,15,vz)
                worksheet.write(0,16,"v")
                worksheet.write_column(1,16,v)
                worksheet.write(0,17,"M")
                worksheet.write_column(1,17,M)
                worksheet.write(0,18,"E")
                worksheet.write_column(1,18,E)
                worksheet.write(0,19,"KI")
                worksheet.write_column(1,19,KI)
                
                worksheet.write(0,20,"x")
                worksheet.write(1,20,x_max)
                worksheet.write(2,20,x_min)
                worksheet.write(0,21,"y")
                worksheet.write(1,21,y_max)
                worksheet.write(2,21,y_min)
                worksheet.write(0,22,"z")
                worksheet.write(1,22,z_max)
                worksheet.write(2,22,z_min)
                worksheet.write(0,23,"v")
                worksheet.write(1,23,v_max)
                worksheet.write(2,23,v_min)
                worksheet.write(0,24,"M")
                worksheet.write(1,24,M_max)
                worksheet.write(2,24,M_min)
                worksheet.write(0,25,"E")
                worksheet.write(1,25,E_max)
                worksheet.write(2,25,E_min)
                worksheet.write(0,26,"KI")
                worksheet.write(1,26,KI_max)
                worksheet.write(2,26,KI_min)
                workbook.close()

            btn=warna_ganti(printa,text="printnya",command=printnya).grid(column=0,row=3)    
        
        lbl=Label(yuhu, text="RK45B ",width=15).grid(column=0,row=1)
        lbl=Label(yuhu, text="pilih output ",width=10).grid(column=0,row=2)
        lbl=Label(yuhu, text=" ",width=10).grid(column=0,row=3)
        lbl=Label(yuhu, text="posisi ",width=10).grid(column=0,row=4)

        lbl=Label(yuhu, text="x,t",width=10).grid(column=0,row=5)        
        btn=warna_ganti(yuhu,text="_",command=gambarx,bg="green",fg="white").grid(column=1,row=5)
        lbl=Label(yuhu, text="y,t",width=10).grid(column=0,row=6)
        btn=warna_ganti(yuhu,text="_",command=gambary,bg="green",fg="white").grid(column=1,row=6)
        lbl=Label(yuhu, text="z,t",width=10).grid(column=0,row=7)
        btn=warna_ganti(yuhu,text="_",command=gambarz,bg="green",fg="white").grid(column=1,row=7)
        
        lbl=Label(yuhu, text="kecepatan mutlak ",width=15).grid(column=0,row=8)
        lbl=Label(yuhu, text="v,t ",width=10).grid(column=0,row=9)
        btn=warna_ganti(yuhu,text="_",command=gambarv,bg="green",fg="white").grid(column=1,row=9)
        
        lbl=Label(yuhu, text="massa gerak ",width=10).grid(column=0,row=10)
        lbl=Label(yuhu, text="M,t",width=10).grid(column=0,row=11)
        btn=warna_ganti(yuhu,text="_",command=gambarm,bg="green",fg="white").grid(column=1,row=11)
        
        lbl=Label(yuhu, text="energi gerak ",width=10).grid(column=0,row=12)
        lbl=Label(yuhu, text="E,t",width=10).grid(column=0,row=13)
        btn=warna_ganti(yuhu,text="_",command=gambarE,bg="green",fg="white").grid(column=1,row=13)
        
        lbl=Label(yuhu, text="energi kinetik ",width=10).grid(column=0,row=14)
        lbl=Label(yuhu, text="K,t",width=10).grid(column=0,row=15)
        btn=warna_ganti(yuhu,text="_",command=gambarKI,bg="green",fg="white").grid(column=1,row=15)
        
        lbl=Label(yuhu, text="Simulasi gerak ",width=10).grid(column=0,row=16)
        lbl=Label(yuhu, text="2D ",width=10).grid(column=0,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar2d,bg="green",fg="white").grid(column=1,row=17)
        
        lbl=Label(yuhu, text="3D ",width=10).grid(column=0,row=18)
        btn=warna_ganti(yuhu,text="_",command=gambar3d,bg="green",fg="white").grid(column=1,row=18)
        
        lbl=Label(yuhu, text="print xlsx ",width=10).grid(column=0,row=19)
        lbl=Label(yuhu, text="xlsx ",width=10).grid(column=0,row=20)
        btn=warna_ganti(yuhu,text="_",command=hasilkanxl,bg="green",fg="white").grid(column=1,row=20)
        yuhu.mainloop()
                    
    def rk45F():
        yuhu= Tk()
        yuhu.title("agrafik yang diinginkan")
        yuhu.geometry("240x470+1110+0")
        def funcy (t , vy , y , k, B0, B1, y0, vx0 ):
            return vy
        def guncy (t , vy , y , k, B0, B1, y0, vx0 ):
            return -(k**2) * ((B0+B1*t)**2) *(y-y0)-((k) * ((B0+B1*t)) *vx0)
        def funcx (t , vx , x , k, B0, B1, x0, vy0  ):
            return vx
        def guncx (t , vx , x , k, B0, B1, x0, vy0  ):
            return ((k**2) * ((B0+B1*t)**2) )*(x0-x)+((k) * ((B0+B1*t)) *vy0)
        def funcz (t , vz , z , k, E0, E1 ):
            return vz
        def guncz (t , vz , z , k, E0, E1 ):
            return k*(E0+E1*t)
        for i in range (ndata-1):
                
            h=t[i+1]-t[i]
            k1y=h*funcy(t[i]     , vy[i]          ,y[i],k, B0, B1, y0, vx0 ) 
            l1y=h*guncy(t[i]     , vy[i]          ,y[i],k, B0, B1, y0, vx0 ) 
            k2y=h*funcy(t[i]+(1*h/4)      , vy[i]+(1*l1y/4)  ,y[i]+(1*k1y/4),k, B0, B1, y0, vx0 ) 
            l2y=h*guncy(t[i]+(1*h/4)      , vy[i]+(1*l1y/4)         ,y[i]+(1*k1y/4),k, B0, B1, y0, vx0  )
            k3y=h*funcy(t[i]+(3*h/8)      , vy[i]+(3*l1y/32)+(9*l2y/32)    ,y[i]+(3*k1y/32)+(9*k2y/32),k, B0, B1, y0, vx0  )
            l3y=h*guncy(t[i]+(3*h/8)      , vy[i]+(3*l1y/32)+(9*l2y/32)    ,y[i]+(3*k1y/32)+(9*k2y/32),k, B0, B1, y0, vx0  )
            k4y=h*funcy(t[i]+(12*h/13)    , vy[i]+(1932*l1y/2197)-(7200*l2y/2197)+(7296*l3y/2197)      ,y[i]+(1932*k1y/2197)-(7200*k2y/2197)+(7296*k3y/2197),k, B0, B1, y0, vx0  )
            l4y=h*guncy(t[i]+(12*h/13)    , vy[i]+(1932*l1y/2197)-(7200*l2y/2197)+(7296*l3y/2197)      ,y[i]+(1932*k1y/2197)-(7200*k2y/2197)+(7296*k3y/2197),k, B0, B1, y0, vx0  )
            k5y=h*funcy(t[i]+(h)   , vy[i]+(439*l1y/216)-(8*l2y)+(3680*l3y/513)-(845*l4y/4104)  ,y[i]+(439*k1y/216)-(8*k2y)+(3680*k3y/513)-(845*k4y/4104),k, B0, B1, y0, vx0  )
            l5y=h*guncy(t[i]+(h)   , vy[i]+(439*l1y/216)-(8*l2y)+(3680*l3y/513)-(845*l4y/4104)  ,y[i]+(439*k1y/216)-(8*k2y)+(3680*k3y/513)-(845*k4y/4104),k, B0, B1, y0, vx0  )
            k6y=h*funcy(t[i]+(1*h/2)      , vy[i]-(8*l1y/27)+(2*l2y)-(3544*l3y/2565)+(1859*l4y/4104)-(11*l5y/40)  ,y[i]-(8*k1y/27)+(2*k2y)-(3544*k3y/2565)+(1859*k4y/4104)-(11*k5y/40),k, B0, B1, y0, vx0 )
            l6y=h*guncy(t[i]+(1*h/2)      , vy[i]-(8*l1y/27)+(2*l2y)-(3544*l3y/2565)+(1859*l4y/4104)-(11*l5y/40)  ,y[i]-(8*k1y/27)+(2*k2y)-(3544*k3y/2565)+(1859*k4y/4104)-(11*k5y/40),k, B0, B1, y0, vx0 )   
                
            k1x=h*funcx(t[i]       , vx[i]            ,x[i],k, B0, B1, x0, vy0 )
            l1x=h*guncx(t[i]       , vx[i]            ,x[i],k, B0, B1, x0, vy0 )
            k2x=h*funcx(t[i]+(1*h/4)      , vx[i]+(1*l1x/4)         ,x[i]+(1*k1x/4),k, B0, B1, x0, vy0 )
            l2x=h*guncx(t[i]+(1*h/4)      , vx[i]+(1*l1x/4)         ,x[i]+(1*k1x/4),k, B0, B1, x0, vy0 )
            k3x=h*funcx(t[i]+(3*h/8)      , vx[i]+(3*l1x/32)+(9*l2x/32)    ,x[i]+(3*k1x/32)+(9*k2x/32),k, B0, B1, x0, vy0 )
            l3x=h*guncx(t[i]+(3*h/8)      , vx[i]+(3*l1x/32)+(9*l2x/32)    ,x[i]+(3*k1x/32)+(9*k2x/32),k, B0, B1, x0, vy0 )
            k4x=h*funcx(t[i]+(12*h/13)    , vx[i]+(1932*l1x/2197)-(7200*l2x/2197)+(7296*l3x/2197)      ,x[i]+(1932*k1x/2197)-(7200*k2x/2197)+(7296*k3x/2197),k, B0, B1, x0, vy0 )
            l4x=h*guncx(t[i]+(12*h/13)    , vx[i]+(1932*l1x/2197)-(7200*l2x/2197)+(7296*l3x/2197)      ,x[i]+(1932*k1x/2197)-(7200*k2x/2197)+(7296*k3x/2197),k, B0, B1, x0, vy0 )
            k5x=h*funcx(t[i]+(h)   , vx[i]+(439*l1x/216)-(8*l2x)+(3680*l3x/513)-(845*l4x/4104)  ,x[i]+(439*k1x/216)-(8*k2x)+(3680*k3x/513)-(845*k4x/4104),k, B0, B1, x0, vy0 )
            l5x=h*guncx(t[i]+(h)   , vx[i]+(439*l1x/216)-(8*l2x)+(3680*l3x/513)-(845*l4x/4104)  ,x[i]+(439*k1x/216)-(8*k2x)+(3680*k3x/513)-(845*k4x/4104),k, B0, B1, x0, vy0 )
            k6x=h*funcx(t[i]+(1*h/2)      , vx[i]-(8*l1x/27)+(2*l2x)-(3544*l3x/2565)+(1859*l4x/4104)-(11*l5x/40)  ,x[i]-(8*k1x/27)+(2*k2x)-(3544*k3x/2565)+(1859*k4x/4104)-(11*k5x/40),k, B0, B1, x0, vy0)
            l6x=h*guncx(t[i]+(1*h/2)      , vx[i]-(8*l1x/27)+(2*l2x)-(3544*l3x/2565)+(1859*l4x/4104)-(11*l5x/40)  ,x[i]-(8*k1x/27)+(2*k2x)-(3544*k3x/2565)+(1859*k4x/4104)-(11*k5x/40),k, B0, B1, x0, vy0)   
            
            k1z=h*funcz(t[i]       , vz[i]            ,z[i], k, E0, E1  )
            l1z=h*guncz(t[i]       , vz[i]            ,z[i], k, E0, E1  )
            k2z=h*funcz(t[i]+(1*h/4)      , vz[i]+(1*l1z/4)         ,z[i]+(1*k1z/4), k, E0, E1  )
            l2z=h*guncz(t[i]+(1*h/4)      , vz[i]+(1*l1z/4)         ,z[i]+(1*k1z/4) , k, E0, E1  )
            k3z=h*funcz(t[i]+(3*h/8)      , vz[i]+(3*l1z/32)+(9*l2z/32)    ,z[i]+(3*k1z/32)+(9*k2z/32) , k, E0, E1  )
            l3z=h*guncz(t[i]+(3*h/8)      , vz[i]+(3*l1z/32)+(9*l2z/32)    ,z[i]+(3*k1z/32)+(9*k2z/32) , k, E0, E1  )
            k4z=h*funcz(t[i]+(12*h/13)    , vz[i]+(1932*l1z/2197)-(7200*l2z/2197)+(7296*l3z/2197)      ,z[i]+(1932*k1z/2197)-(7200*k2z/2197)+(7296*k3z/2197) , k, E0, E1  )
            l4z=h*guncz(t[i]+(12*h/13)    , vz[i]+(1932*l1z/2197)-(7200*l2z/2197)+(7296*l3z/2197)      ,z[i]+(1932*k1z/2197)-(7200*k2z/2197)+(7296*k3z/2197) , k, E0, E1  )
            k5z=h*funcz(t[i]+(h)   , vz[i]+(439*l1z/216)-(8*l2z)+(3680*l3z/513)-(845*l4z/4104)  ,z[i]+(439*k1z/216)-(8*k2z)+(3680*k3z/513)-(845*k4z/4104) , k, E0, E1  )
            l5z=h*guncz(t[i]+(h)   , vz[i]+(439*l1z/216)-(8*l2z)+(3680*l3z/513)-(845*l4z/4104)  ,z[i]+(439*k1z/216)-(8*k2z)+(3680*k3z/513)-(845*k4z/4104) , k, E0, E1  )
            k6z=h*funcz(t[i]+(1*h/2)      , vz[i]-(8*l1z/27)+(2*l2z)-(3544*l3z/2565)+(1859*l4z/4104)-(11*l5z/40)  ,z[i]-(8*k1z/27)+(2*k2z)-(3544*k3z/2565)+(1859*k4z/4104)-(11*k5z/40), k, E0, E1  )
            l6z=h*guncz(t[i]+(1*h/2)      , vz[i]-(8*l1z/27)+(2*l2z)-(3544*l3z/2565)+(1859*l4z/4104)-(11*l5z/40)  ,z[i]-(8*k1z/27)+(2*k2z)-(3544*k3z/2565)+(1859*k4z/4104)-(11*k5z/40), k, E0, E1  )   
                     
            kky=(16*k1y/135+6656*k3y/12825+28561*k4y/56430-9*k5y/50+2*k6y/55)
            ly=(16*l1y/135+6656*l3y/12825+28561*l4y/56430-9*l5y/50+2*l6y/55)
                  
            y[i+1]=y[i]+kky
            vy[i+1]=vy[i]+ly
                  
            kkx=(16*k1x/135+6656*k3x/12825+28561*k4x/56430-9*k5x/50+2*k6x/55)
            lx=(16*l1x/135+6656*l3x/12825+28561*l4x/56430-9*l5x/50+2*l6x/55)
              
            x[i+1]=x[i]+kkx
            vx[i+1]=vx[i]+lx
              
            kkz=(16*k1z/135+6656*k3z/12825+28561*k4z/56430-9*k5z/50+2*k6z/55)
            lz=(16*l1z/135+6656*l3z/12825+28561*l4z/56430-9*l5z/50+2*l6z/55)
                  
            z[i+1]=z[i]+kkz
            vz[i+1]=vz[i]+lz
            v=np.sqrt(vx**2+vy**2+vz**2)
            g=np.sqrt(1-(v**2/c**2))
            
            M=m/g
            E=M*(c**2)
            KI=(c**2)*M-(c**2)*m
            x_max=np.max(x)
            x_min=np.min(x)
            
            y_max=np.max(y)
            y_min=np.min(y)
            
            z_max=np.max(z)
            z_min=np.min(z)
            
            v_max=np.max(v)
            v_min=np.min(v)
            
            M_max=np.max(M)
            M_min=np.min(M)
            
            E_max=np.max(E)
            E_min=np.min(E)
            
            KI_max=np.max(KI)
            KI_min=np.min(KI)
            
        def gambarx():
            root = tkinter.Tk()
            root.title('ploting posisi x, t pada RK45F')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,x,'r',label="posisi x terhadap waktu")
            plot1.set_ylabel('posisi x(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambary():
            root = tkinter.Tk()
            root.title('ploting posisi y, t pada RK45F ')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,y,'b',label="posisi y terhadap waktu")
            plot1.set_ylabel('posisi z(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarz():
            root = tkinter.Tk()
            root.title('ploting posisi z, t pada RK45F ')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,z,'g',label="posisi z terhadap waktu")
            plot1.set_ylabel('posisi z(m)')
            plot1.set_xlabel('waktu t(s)')
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarv():
            root = tkinter.Tk()
            root.title('ploting kecepatan mutlak dan waktu pada RK45F')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,v,'y',label="kecepatan mutlak v terhadap waktu")
            plot1.set_ylabel('kecepatan V(m/s)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {v_min}', xy=(t[np.argmin(v)], v_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {v_max}', xy=(t[np.argmax(v)], v_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarm():
            root = tkinter.Tk()
            root.title('ploting massa dan waktu pada RK45F')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,M,'g',label="massa M terhadap waktu")
            plot1.set_ylabel('massa gerak M(kg)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {M_min}', xy=(t[np.argmin(M)], M_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {M_max}', xy=(t[np.argmax(M)], M_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarKI():
            root = tkinter.Tk()
            root.title('ploting energi kinetik waktu pada RK45F')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,KI,'g',label="energi kinetik terhadap waktu")
            plot1.set_ylabel('energi kinetik KI(j)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {KI_min}', xy=(t[np.argmin(KI)], KI_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {KI_max}', xy=(t[np.argmax(KI)], KI_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambarE():
            root = tkinter.Tk()
            root.title('ploting energi pada Rk45F')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(t,E,'r',label="energi gerak")
            plot1.set_ylabel('energi E(j)')
            plot1.set_xlabel('waktu t(s)')
            plot1.annotate(f'Min: {E_min}', xy=(t[np.argmin(E)], E_min), xytext=(10, 30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.annotate(f'Max: {E_max}', xy=(t[np.argmax(E)], E_max), xytext=(-120, -30),
                           textcoords='offset points', arrowprops=dict(arrowstyle="->"))
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
                canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop()  
        def gambar2d():
            root = tkinter.Tk()
            root.title('ploting medan magnet pada RK45F')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111)
            plot1.plot(x,y,'y',label="posisi x terhadap y")
            plot1.set_xlabel('posisi x')
            plot1.set_ylabel('posisi y')
            plot1.plot(x[0],y[0],marker='o',markersize=8,color='blue',label="stars")
            plot1.plot(x[-1],y[-1],marker='o',markersize=8,color='green',label="end")
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
            canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop() 
        def gambar3d():
            root = tkinter.Tk()
            root.title('ploting medan elektromagnet pada RK45F')
            root.geometry("850x700+220+0")
            fig = Figure(figsize = (6,3),dpi = 100)
            plot1 = fig.add_subplot(111,projection='3d')
            plot1.plot(x,y,z,'r',label="posisi")
            plot1.set_xlabel('posisi x(m)')
            plot1.set_ylabel('posisi y(m)')
            plot1.set_zlabel('posisi z(m)')
            plot1.plot(x[0],y[0],z[0],marker='o',markersize=8,color='blue',label="stars")
            plot1.plot(x[-1],y[-1],z[-1],marker='o',markersize=8,color='green',label="end")
            plot1.legend(loc="upper left")
            canvas = FigureCanvasTkAgg(fig,master = root)
            canvas.draw()
            canvas.get_tk_widget().place(x=400,y=400)
            canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas,root)
            toolbar.update()
            def on_key_press(event):
                print("you pressed {}".format(event.key))
                key_press_handler(event, canvas, toolbar)
                canvas.mpl_connect("key_press_event", on_key_press)
            def _quit():
                root.quit()   
                root.destroy()
            button = tkinter.Button(master=root, text="Quit", command=_quit)
            button.pack(side=tkinter.BOTTOM)
            tkinter.mainloop()  
        def hasilkanxl():
            printa = Tk()
            printa.title("hasilprint pada RK45F")
            printa.geometry("240x100+1110+520")
            lbl=Label(printa, text="nama file    .xlsx",width=15)
            lbl.grid(column=0,row=0)
            namafile=Entry(printa, width=10)
            namafile.grid(column=1,row=0)
            def printnya():
                workbook=xl.Workbook(namafile.get())
                worksheet=workbook.add_worksheet()
                worksheet.write(0,0,"q")
                worksheet.write(1,0,q)
                worksheet.write(0,1,"B0")
                worksheet.write(1,1,B0)
                worksheet.write(0,2,"B1")
                worksheet.write(1,2,B1)
                worksheet.write(0,3,"m")
                worksheet.write(1,3,m)
                worksheet.write(0,4,"E0")
                worksheet.write(1,4,E0)
                worksheet.write(0,5,"E1")
                worksheet.write(1,5,E1)
                worksheet.write(0,6,"t0")
                worksheet.write(1,6,t0)
                worksheet.write(0,7,"tn")
                worksheet.write(1,7,tn)
                worksheet.write(0,8,"ndata")
                worksheet.write(1,8,ndata)
                worksheet.write(0,9,"waktu")
                worksheet.write_column(1,9,t)
                worksheet.write(0,10,"x")
                worksheet.write_column(1,10,x)
                worksheet.write(0,11,"y")
                worksheet.write_column(1,11,y)
                worksheet.write(0,12,"z")
                worksheet.write_column(1,12,z)
                worksheet.write(0,13,"vx")
                worksheet.write_column(1,13,vx)
                worksheet.write(0,14,"vy")
                worksheet.write_column(1,14,vy)
                worksheet.write(0,15,"vz")
                worksheet.write_column(1,15,vz)
                worksheet.write(0,16,"v")
                worksheet.write_column(1,16,v)
                worksheet.write(0,17,"M")
                worksheet.write_column(1,17,M)
                worksheet.write(0,18,"E")
                worksheet.write_column(1,18,E)
                worksheet.write(0,19,"KI")
                worksheet.write_column(1,19,KI)
                
                worksheet.write(0,20,"x")
                worksheet.write(1,20,x_max)
                worksheet.write(2,20,x_min)
                worksheet.write(0,21,"y")
                worksheet.write(1,21,y_max)
                worksheet.write(2,21,y_min)
                worksheet.write(0,22,"z")
                worksheet.write(1,22,z_max)
                worksheet.write(2,22,z_min)
                worksheet.write(0,23,"v")
                worksheet.write(1,23,v_max)
                worksheet.write(2,23,v_min)
                worksheet.write(0,24,"M")
                worksheet.write(1,24,M_max)
                worksheet.write(2,24,M_min)
                worksheet.write(0,25,"E")
                worksheet.write(1,25,E_max)
                worksheet.write(2,25,E_min)
                worksheet.write(0,26,"KI")
                worksheet.write(1,26,KI_max)
                worksheet.write(2,26,KI_min)
                workbook.close()

            btn=warna_ganti(printa,text="printnya",command=printnya).grid(column=0,row=3)    
        
        lbl=Label(yuhu, text="RK45F ",width=15).grid(column=0,row=1)
        lbl=Label(yuhu, text="pilih output ",width=10).grid(column=0,row=2)
        lbl=Label(yuhu, text=" ",width=10).grid(column=0,row=3)
        lbl=Label(yuhu, text="posisi ",width=10).grid(column=0,row=4)

        lbl=Label(yuhu, text="x,t",width=10).grid(column=0,row=5)        
        btn=warna_ganti(yuhu,text="_",command=gambarx,bg="green",fg="white").grid(column=1,row=5)
        lbl=Label(yuhu, text="y,t",width=10).grid(column=0,row=6)
        btn=warna_ganti(yuhu,text="_",command=gambary,bg="green",fg="white").grid(column=1,row=6)
        lbl=Label(yuhu, text="z,t",width=10).grid(column=0,row=7)
        btn=warna_ganti(yuhu,text="_",command=gambarz,bg="green",fg="white").grid(column=1,row=7)
        
        lbl=Label(yuhu, text="kecepatan mutlak ",width=15).grid(column=0,row=8)
        lbl=Label(yuhu, text="v,t ",width=10).grid(column=0,row=9)
        btn=warna_ganti(yuhu,text="_",command=gambarv,bg="green",fg="white").grid(column=1,row=9)
        
        lbl=Label(yuhu, text="massa gerak ",width=10).grid(column=0,row=10)
        lbl=Label(yuhu, text="M,t",width=10).grid(column=0,row=11)
        btn=warna_ganti(yuhu,text="_",command=gambarm,bg="green",fg="white").grid(column=1,row=11)
        
        lbl=Label(yuhu, text="energi gerak ",width=10).grid(column=0,row=12)
        lbl=Label(yuhu, text="E,t",width=10).grid(column=0,row=13)
        btn=warna_ganti(yuhu,text="_",command=gambarE,bg="green",fg="white").grid(column=1,row=13)
        
        lbl=Label(yuhu, text="energi kinetik ",width=10).grid(column=0,row=14)
        lbl=Label(yuhu, text="K,t",width=10).grid(column=0,row=15)
        btn=warna_ganti(yuhu,text="_",command=gambarKI,bg="green",fg="white").grid(column=1,row=15)
        
        lbl=Label(yuhu, text="Simulasi gerak ",width=10).grid(column=0,row=16)
        lbl=Label(yuhu, text="2D ",width=10).grid(column=0,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar2d,bg="green",fg="white").grid(column=1,row=17)
        
        lbl=Label(yuhu, text="3D ",width=10).grid(column=0,row=18)
        btn=warna_ganti(yuhu,text="_",command=gambar3d,bg="green",fg="white").grid(column=1,row=18)
        
        lbl=Label(yuhu, text="print xlsx ",width=10).grid(column=0,row=19)
        lbl=Label(yuhu, text="xlsx ",width=10).grid(column=0,row=20)
        btn=warna_ganti(yuhu,text="_",command=hasilkanxl,bg="green",fg="white").grid(column=1,row=20)
        yuhu.mainloop()
        
    
    lbl22=Label(nonlinier, text="per.nonlinier ",width=15).grid(column=0,row=1)
    lbl23=Label(nonlinier, text="pilih metode ",width=10).grid(column=0,row=2)
    lbl24=Label(nonlinier, text=" ",width=10).grid(column=0,row=3)
    lbl25=Label(nonlinier, text="rk4 ",width=10).grid(column=0,row=4)
    lbl26=Label(nonlinier, text="rk45B ",width=10).grid(column=0,row=5)
    lbl27=Label(nonlinier, text="rk45F ",width=10).grid(column=0,row=6)
    
    btn66=warna_ganti(nonlinier,text="_",command=rk4).grid(column=1,row=4)
    btn77=warna_ganti(nonlinier,text="_",command=rk45B).grid(column=1,row=5)
    btn88=warna_ganti(nonlinier,text="_",command=rk45F).grid(column=1,row=6)
    nonlinier.mainloop()   
    
def linier():
    linier = Tk()
    linier.title("linier ")
    linier.geometry("200x200+0+480")

    qp=float(aq.get() or 1)
    mp=float(am.get() or 1)
    qT=float(aqT.get() or 0)
    mT=float(amT.get() or 0)
    q=qp*(10**qT)
    m=mp*(10**mT)
    
    k=q/m
    
    E0p=float(aE0.get() or 1)
    E1p=float(aE1.get() or 1)
    E0T=float(aE0T.get() or 0)
    E1T=float(aE1T.get() or 0)
    E0=E0p*(10**E0T)
    E1=E1p*(10**E1T)
    
    B0p=float(aB0.get() or 1)
    B1p=float(aB1.get() or 1)
    B0T=float(aB0T.get() or 0)
    B1T=float(aB1T.get() or 0)
    B0=B0p*(10**B0T)
    B1=B1p*(10**B1T)
    
    t0p=float(at0.get() or 0)
    tnp=float(atn.get() or 10)
    ndatap=int(andata.get() or 1000)
    t0T=float(at0T.get() or 0)
    tnT=float(atnT.get() or 0)
    ndataT=int(andataT.get() or 0)
    t0=t0p*(10**t0T)
    tn=tnp*(10**tnT)
    ndata=ndatap*(10**ndataT)
    t=np.linspace(t0,tn,ndata)
    
    c=3*(10**9)
    
    x0p=float(ax0.get() or 1)
    x0T=float(ax0T.get() or 0)
    x0=x0p*(10**x0T)
    x=np.zeros(ndata)
    x[0]=x0
    
    
    y0p=float(ay0.get() or 1)
    y0T=float(ay0T.get() or 0)
    y0=y0p*(10**y0T)
    y=np.zeros(ndata)
    y[0]=y0
    
    z0p=float(az0.get() or 1)
    z0T=float(az0T.get() or 0)
    z0=z0p*(10**z0T)
    z=np.zeros(ndata)
    z[0]=y0
    
    vx0p=float(avx0.get() or 1)
    vx0T=float(avx0T.get() or 0)
    vx0=vx0p*(10**vx0T)
    vx=np.zeros(ndata)
    vx[0]=vx0
    
    vy0p=float(avy0.get() or 1)
    vy0T=float(avy0T.get() or 0)
    vy0=vy0p*(10**vy0T)
    vy=np.zeros(ndata)
    vy[0]=vy0
    
    vz0p=float(avz0.get() or 1)
    vz0T=float(avz0T.get() or 0)
    vz0=vz0p*(10**vz0T)
    vz=np.zeros(ndata)
    vz[0]=vz0
      

    def rk4():
        yuhu= Tk()
        yuhu.title("grafik yang diinginkan")
        yuhu.geometry("240x475+1110+0")
        
        def funcy (t , vy , y , k, B0, B1,  y0, vx0 ):
             return vy
        def guncy (t , vy , y , k, B0, B1,  y0, vx0 ):
             return -(k**2) * ((B0+B1*t)**2) *(y-y0)-((k) * ((B0+B1*t)) *vx0)
        def funcx (t , vx , y , k, B0, B1, x0, vy0  ):
             return vx
        def guncx (t , vx , x , k, B0, B1, x0, vy0  ):
             return ((k**2) * ((B0+B1*t)**2) )*(x0-x)+((k) * ((B0+B1*t)) *vy0)
        def funcz (t , vz , z , k, E0, E1):
             return vz
        def guncz (t , vz , z , k, E0, E1):
             return k*(E0+E1*t)
        for i in range (ndata-1):
            
            h=t[i+1]-t[i]
            k1y=h*funcy(t[i],vy[i],y[i] ,k, B0, B1, y0, vx0  )
            l1y=h*guncy(t[i],vy[i],y[i] ,k, B0, B1, y0, vx0  )
            k2y=h*funcy(t[i]+0.5*h,vy[i]+(0.5*l1y),y[i]+(0.5*k1y) ,k, B0, B1, y0, vx0  )
            l2y=h*guncy(t[i]+0.5*h,vy[i]+(0.5*l1y),y[i]+(0.5*k1y) ,k, B0, B1, y0, vx0  )
            k3y=h*funcy(t[i]+0.5*h,vy[i]+(0.5*l2y),y[i]+(0.5*k2y) ,k, B0, B1, y0, vx0  )
            l3y=h*guncy(t[i]+0.5*h,vy[i]+(0.5*l2y),y[i]+(0.5*k2y) ,k, B0, B1, y0, vx0  )
            k4y=h*funcy(t[i]+h,vy[i]+l3y,y[i]+k3y ,k, B0, B1, y0, vx0  )
            l4y=h*guncy(t[i]+h,vy[i]+l3y,y[i]+k3y ,k, B0, B1, y0, vx0  )
            
            k1x=h*funcx(t[i],vx[i],x[i] ,k, B0, B1, x0, vy0  )
            l1x=h*guncx(t[i],vx[i],x[i] ,k, B0, B1, x0, vy0  )
            k2x=h*funcx(t[i]+0.5*h,vx[i]+(0.5*l1x),x[i]+(0.5*k1x) ,k, B0, B1, x0, vy0  )
            l2x=h*guncx(t[i]+0.5*h,vx[i]+(0.5*l1x),x[i]+(0.5*k1x) ,k, B0, B1, x0, vy0  )
            k3x=h*funcx(t[i]+0.5*h,vx[i]+(0.5*l2x),x[i]+(0.5*k2x) ,k, B0, B1, x0, vy0  )
            l3x=h*guncx(t[i]+0.5*h,vx[i]+(0.5*l2x),x[i]+(0.5*k2x) ,k, B0, B1, x0, vy0  )
            k4x=h*funcx(t[i]+h,vx[i]+l3x,x[i]+k3x ,k, B0, B1, x0, vy0  )
            l4x=h*guncx(t[i]+h,vx[i]+l3x,x[i]+k3x ,k, B0, B1, x0, vy0  )
            
            k1z=h*funcz(t[i]  , vz[i] , z[i], k, E0, E1) 
            l1z=h*guncz(t[i]  , vz[i] , z[i], k, E0, E1)
            k2z=h*funcz(t[i]+0.5*h , vz[i]+(0.5*l1z) , z[i]+(0.5*k1z), k, E0, E1)
            l2z=h*guncz(t[i]+0.5*h , vz[i]+(0.5*l1z) , z[i]+(0.5*k1z), k, E0, E1)
            k3z=h*funcz(t[i]+0.5*h , vz[i]+(0.5*l2z) , z[i]+(0.5*k2z), k, E0, E1)
            l3z=h*guncz(t[i]+0.5*h , vz[i]+(0.5*l2z) , z[i]+(0.5*k2z), k, E0, E1)
            k4z=h*funcz(t[i]+h , vz[i]+l3z , z[i]+k3z , k, E0, E1)
            l4z=h*guncz(t[i]+h , vz[i]+l3z , z[i]+k3z , k, E0, E1)
            
            kky=(k1y+2*k2y+2*k3y+k4y)/6.0
            ly=(l1y+2*l2y+2*l3y+l4y)/6.0
            
            y[i+1]=y[i]+kky
            vy[i+1]=vy[i]+ly
            
            kkx=(k1x+2*k2x+2*k3x+k4x)/6.0
            lx=(l1x+2*l2x+2*l3x+l4x)/6.0

            x[i+1]=x[i]+kkx
            vx[i+1]=vx[i]+lx

            kkz=(k1z+2*k2z+2*k3z+k4z)/6.0
            lz=(l1z+2*l2z+2*l3z+l4z)/6.0
            
            z[i+1]=z[i]+kkz
            vz[i+1]=vz[i]+lz
            
            vya=-(vx0)*np.sin(k*B0*t)+(vy0)*np.cos(k*B0*t)
            ya=(vx0/(k*B0))*np.cos(k*B0*t)+(vy0/k)*np.sin(k*B0*t)+y0-(vx0/(k*B0))
            kesy=abs((ya-y)/ya)
            ksy=kesy*100
              
            vxa=(vy0)*np.sin((k*B0)*t)+(vx0)*np.cos((k*B0)*t)
            xa=(-vy0/(k*B0))*np.cos((k*B0)*t)+(vx0/(k*B0))*np.sin((k*B0)*t)+x0+vy0/(k*B0)
            kesx=abs((xa-x)/xa)
            ksx=kesx*100
            
            vza= k*(E0*t+E1*1/2* t**2)+vz0
            za= k*(E0*1/2*t**2+E1*1/6* t**3)+vz0*t+z0
            kesz=abs((za-z)/za)
            ksz=kesz*100
             
            v=np.sqrt(vx**2+vy**2+vz**2)
            va=np.sqrt(vxa**2+vya**2+vza**2)
            kesv=abs((va-v)/va)
            ksv=kesv*100
             
            g=np.sqrt(1-(v**2/c**2))
            ga=np.sqrt(1-(va**2/c**2))
             
            M=m/g
            Ma=m/ga
            kem=abs((M-Ma)/Ma)
            km=kem*100
             
            E=M*(c**2)
            Ea=Ma*(c**2)
            keE=abs((E-Ea)/Ea)
            kE=keE*100
             
            KI=(c**2)*M-(c**2)*m
            KIa=(c**2)*Ma-(c**2)*m
            keKI=abs((KI-KIa)/KIa)
            kKI=keKI*100
            
            x_max=np.max(x)
            x_min=np.min(x)
            
            y_max=np.max(y)
            y_min=np.min(y)
            
            z_max=np.max(z)
            z_min=np.min(z)
            
            xa_max=np.max(xa)
            xa_min=np.min(xa)
            
            ya_max=np.max(ya)
            ya_min=np.min(ya)
            
            za_max=np.max(za)
            za_min=np.min(za)
            
            v_max=np.max(v)
            v_min=np.min(v)
            
            M_max=np.max(M)
            M_min=np.min(M)
            
            E_max=np.max(E)
            E_min=np.min(E)
            
            KI_max=np.max(KI)
            KI_min=np.min(KI)
            
            va_max=np.max(va)
            va_min=np.min(va)
            
            Ma_max=np.max(Ma)
            Ma_min=np.min(Ma)
            
            Ea_max=np.max(Ea)
            Ea_min=np.min(Ea)
            
            KIa_max=np.max(KIa)
            KIa_min=np.min(KIa)
            
            ksy_max=np.max(ksy)
            ksy_min=np.min(ksy)
            ksy_mean=np.mean(ksy)
            
            ksx_max=np.max(ksx)
            ksx_min=np.min(ksx)
            ksx_mean=np.mean(ksx)
            
            ksz_max=np.max(ksz)
            ksz_min=np.min(ksz)
            ksz_mean=np.mean(ksz)
            
            ksv_max=np.max(ksv)
            ksv_min=np.min(ksv)
            ksv_mean=np.mean(ksv)
            
            km_max=np.max(km)
            km_min=np.min(km)
            km_mean=np.mean(km)
            
            kE_max=np.max(kE)
            kE_min=np.min(kE)
            kE_mean=np.mean(kE)
            
            kKI_max=np.max(kKI)
            kKI_min=np.min(kKI)
            kKI_mean=np.mean(kKI)
            
        def gambarx():
             root = tkinter.Tk()
             root.title('ploting posisi x, t pada RK4 ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,x,'g',label="posisi x terhadap waktu")
             plot1.set_ylabel('posisi x(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,xa,'y',label="posisi xa terhadap waktu")
             plot2.set_ylabel('posisi xa(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambary():
             root = tkinter.Tk()
             root.title('ploting posisi y, t pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,y,'g',label="posisi y terhadap waktu")
             plot1.set_ylabel('posisi y(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,ya,'y',label="posisi ya terhadap waktu")
             plot2.set_ylabel('posisi ya(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
        def gambarz():
             root = tkinter.Tk()
             root.title('ploting posisi z, t pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,z,'g',label="posisi z terhadap waktu")
             plot1.set_ylabel('posisi z(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,za,'y',label="posisi za terhadap waktu")
             plot2.set_ylabel('posisi za(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambarkx():
             root = tkinter.Tk()
             root.title('ploting kes.posisi x pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksx,'r',label="posisi kes.x terhadap waktu")
             plot1.set_ylabel('kesalahan kx(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot1.annotate(f'Min: {ksx_min}', xy=(t[np.argmin(ksx)], ksx_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksx_max}', xy=(t[np.argmax(ksx)], ksx_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksx_mean}', xy=(t[len(t) // 2], ksx_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))

             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambarky():
             root = tkinter.Tk()
             root.title('ploting kes.posisi y pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksy,'r',label="posisi kes.y terhadap waktu")
             plot1.set_ylabel('kesalahan ky(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot1.annotate(f'Min: {ksy_min}', xy=(t[np.argmin(ksy)], ksy_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksy_max}', xy=(t[np.argmax(ksy)], ksy_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksy_mean}', xy=(t[len(t) // 2], ksy_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
        def gambarkz():
             root = tkinter.Tk()
             root.title('ploting kes. posisi z pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksz,'r',label="posisi kes.z terhadap waktu")
             plot1.set_ylabel('kesalahan kz(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {ksz_min}', xy=(t[np.argmin(ksz)], ksz_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksz_max}', xy=(t[np.argmax(ksz)], ksz_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksz_mean}', xy=(t[len(t) // 2], ksz_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
              
        def gambarv():
             root = tkinter.Tk()
             root.title('ploting kecepatan v,T pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,v,'g',label="kecepatan mutlak terhadap waktu")
             plot1.set_ylabel('kecepata v (m,s)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {v_min}', xy=(t[np.argmin(v)], v_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {v_max}', xy=(t[np.argmax(v)], v_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,va,'y',label="kecepatan mutlak A terhadap waktu")
             plot2.set_ylabel('kecepata va (m,s)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {va_min}', xy=(t[np.argmin(va)], va_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {va_max}', xy=(t[np.argmax(va)], va_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarkv():
             root = tkinter.Tk()
             root.title('ploting kes.kecepatan v pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.set_ylabel('kesalahan kv(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.plot(t,ksv,'r',label="posis kes.v terhadap waktu")
             plot1.annotate(f'Min: {ksv_min}', xy=(t[np.argmin(ksv)], ksv_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksv_max}', xy=(t[np.argmax(ksv)], ksv_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksv_mean}', xy=(t[len(t) // 2], ksv_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarm():
             root = tkinter.Tk()
             root.title('ploting posisi m, t pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,M,'g',label="massa M terhadap waktu")
             plot1.set_ylabel('massa M(kg)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {M_min}', xy=(t[np.argmin(M)], M_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {M_max}', xy=(t[np.argmax(M)], M_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,Ma,'y',label="massa Ma terhadap waktu")
             plot2.set_ylabel('massa Ma(kg)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {Ma_min}', xy=(t[np.argmin(Ma)], Ma_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {Ma_max}', xy=(t[np.argmax(Ma)], Ma_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarkm():
             root = tkinter.Tk()
             root.title('ploting kes. m pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,km,'r',label="kesalaham massa terhadap waktu")
             plot1.set_ylabel('kesalahan km(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {km_min}', xy=(t[np.argmin(km)], km_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {km_max}', xy=(t[np.argmax(km)], km_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {km_mean}', xy=(t[len(t) // 2], km_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarE():
             root = tkinter.Tk()
             root.title('ploting energi E /t pada RK4 ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,E,'g',label="energi gerak")
             plot1.set_ylabel('energi E(j)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {E_min}', xy=(t[np.argmin(E)], E_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {E_max}', xy=(t[np.argmax(E)], E_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,Ea,'y',label="energi gerak analitik")
             plot2.set_ylabel('energi Ea(j)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {Ea_min}', xy=(t[np.argmin(Ea)], Ea_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {Ea_max}', xy=(t[np.argmax(Ea)], Ea_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def gambarkE():
             root = tkinter.Tk()
             root.title('ploting kes.E pada RK4 ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,kE,'r',label="posisi kes.E terhadap waktu")
             plot1.set_ylabel('kesalahan kE(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {kE_min}', xy=(t[np.argmin(kE)], kE_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {kE_max}', xy=(t[np.argmax(kE)], kE_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {kE_mean}', xy=(t[len(t) // 2], kE_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarKI():
             root = tkinter.Tk()
             root.title('ploting energi kinetik pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,KI,'g',label="energi kinetik ")
             plot1.set_ylabel('energi kinetik KI(j)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {KI_min}', xy=(t[np.argmin(KI)], KI_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {KI_max}', xy=(t[np.argmax(KI)], KI_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,KIa,'y',label="energi kinetik analitik")
             plot2.set_ylabel('energi kinetik KIa(j)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {KIa_min}', xy=(t[np.argmin(KIa)], KIa_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {KIa_max}', xy=(t[np.argmax(KIa)], KIa_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def gambarkKI():
             root = tkinter.Tk()
             root.title('ploting energi kinetik pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,kKI,'r',label="posisi kes energi kinetik terhadap waktu")
             plot1.set_ylabel('kesalahan KI(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {kKI_min}', xy=(t[np.argmin(kKI)], kKI_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {kKI_max}', xy=(t[np.argmax(kKI)], kKI_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {kKI_mean}', xy=(t[len(t) // 2], kKI_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambar2d():
             root = tkinter.Tk()
             root.title('ploting simulasi 2D pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(x,y,'r',label="posisi")
             plot1.set_xlabel('posisi x(m)')
             plot1.set_ylabel('posisi y(m)')
             plot1.plot(x[0],y[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(x[-1],y[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
         
        def gambar2da():
             root = tkinter.Tk()
             root.title('ploting simulasi 2D analitik pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(xa,ya,'b',label="posisi")
             plot1.set_xlabel('posisi xa(m)')
             plot1.set_ylabel('posisi ya(m)')
             plot1.plot(xa[0],ya[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(xa[-1],ya[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
             
        def gambar3d():
             root = tkinter.Tk()
             root.title('ploting simulasi 3d pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111,projection='3d')
             plot1.plot(x,y,z,'r',label="posisi")
             plot1.set_xlabel('posisi x(m)')
             plot1.set_ylabel('posisi y(m)')
             plot1.set_zlabel('posisi z(m)')
             plot1.plot(x[0],y[0],z[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(x[-1],y[-1],z[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
         
        def gambar3da():
             root = tkinter.Tk()
             root.title('ploting simulasi 3d analitik pada RK4')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111,projection='3d')
             plot1.plot(xa,ya,za,'b',label="posisi")
             plot1.set_xlabel('posisi xa(m)')
             plot1.set_ylabel('posisi ya(m)')
             plot1.set_zlabel('posisi za(m)')
             plot1.plot(xa[0],ya[0],za[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(xa[-1],ya[-1],za[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def hasilkanxl():
             printa = Tk()
             printa.title("hasilprint pada RK4")
             printa.geometry("240x100+1110+520")
             lbl=Label(printa, text="nama file    .xlsx",width=15)
             lbl.grid(column=0,row=0)
             namafile=Entry(printa, width=10)
             namafile.grid(column=1,row=0)
             def printnya():
                 workbook=xl.Workbook(namafile.get())
                 worksheet=workbook.add_worksheet()
                 worksheet.write(0,0,"q")
                 worksheet.write(1,0,q)
                 worksheet.write(0,1,"B0")
                 worksheet.write(1,1,B0)
                 worksheet.write(0,2,"B1")
                 worksheet.write(1,2,B1)
                 worksheet.write(0,3,"m")
                 worksheet.write(1,3,m)
                 worksheet.write(0,4,"E0")
                 worksheet.write(1,4,E0)
                 worksheet.write(0,5,"E1")
                 worksheet.write(1,5,E1)
                 worksheet.write(0,6,"t0")
                 worksheet.write(1,6,t0)
                 worksheet.write(0,7,"tn")
                 worksheet.write(1,7,tn)
                 worksheet.write(0,8,"ndata")
                 worksheet.write(1,8,ndata)
                 worksheet.write(0,9,"waktu")
                 worksheet.write_column(1,9,t)
                 worksheet.write(0,10,"x")
                 worksheet.write_column(1,10,x)
                 worksheet.write(0,11,"y")
                 worksheet.write_column(1,11,y)
                 worksheet.write(0,12,"z")
                 worksheet.write_column(1,12,z)
                 worksheet.write(0,13,"vx")
                 worksheet.write_column(1,13,vx)
                 worksheet.write(0,14,"vy")
                 worksheet.write_column(1,14,vy)
                 worksheet.write(0,15,"vz")
                 worksheet.write_column(1,15,vz)
                 worksheet.write(0,16,"v")
                 worksheet.write_column(1,16,v)
                 worksheet.write(0,17,"M")
                 worksheet.write_column(1,17,M)
                 worksheet.write(0,18,"E")
                 worksheet.write_column(1,18,E)
                 worksheet.write(0,19,"KI")
                 worksheet.write_column(1,19,KI)
                 worksheet.write(0,20,"xa")
                 worksheet.write_column(1,20,xa)
                 worksheet.write(0,21,"ya")
                 worksheet.write_column(1,21,ya)
                 worksheet.write(0,22,"za")
                 worksheet.write_column(1,22,za)
                 worksheet.write(0,23,"vxa")
                 worksheet.write_column(1,23,vxa)
                 worksheet.write(0,24,"vya")
                 worksheet.write_column(1,24,vya)
                 worksheet.write(0,25,"vza")
                 worksheet.write_column(1,25,vza)
                 worksheet.write(0,26,"va")
                 worksheet.write_column(1,26,va)
                 worksheet.write(0,27,"Ma")
                 worksheet.write_column(1,27,Ma)
                 worksheet.write(0,28,"Ea")
                 worksheet.write_column(1,28,Ea)
                 worksheet.write(0,29,"KIa")
                 worksheet.write_column(1,29,KIa)
                 
                 worksheet.write(0,30,"xRK4")
                 worksheet.write(1,30,x_max)
                 worksheet.write(2,30,x_min)
                 
                 worksheet.write(3,30,"xA")
                 worksheet.write(4,30,xa_max)
                 worksheet.write(5,30,xa_min)
                 
                 worksheet.write(0,31,"yRK4")
                 worksheet.write(1,31,y_max)
                 worksheet.write(2,31,y_min)
                 
                 worksheet.write(3,31,"yA")
                 worksheet.write(4,31,ya_max)
                 worksheet.write(5,31,ya_min)
                 
                 worksheet.write(0,32,"zRK4")
                 worksheet.write(1,32,z_max)
                 worksheet.write(2,32,z_min)
                 
                 worksheet.write(3,32,"za")
                 worksheet.write(4,32,za_max)
                 worksheet.write(5,32,za_min)
                 
                 worksheet.write(0,33,"vRK4")
                 worksheet.write(1,33,v_max)
                 worksheet.write(2,33,v_min)
                 
                 worksheet.write(3,33,"va")
                 worksheet.write(4,33,va_max)
                 worksheet.write(5,33,va_min)
                 
                 worksheet.write(0,34,"MRK4")
                 worksheet.write(1,34,M_max)
                 worksheet.write(2,34,M_min)
                 
                 worksheet.write(3,34,"Ma")
                 worksheet.write(4,34,Ma_max)
                 worksheet.write(5,34,Ma_min)
                 
                 worksheet.write(0,35,"ERK4")
                 worksheet.write(1,35,E_max)
                 worksheet.write(2,35,E_min)
                 
                 worksheet.write(3,35,"Ea")
                 worksheet.write(4,35,Ea_max)
                 worksheet.write(5,35,Ea_min)
                 
                 worksheet.write(0,36,"KIRK4")
                 worksheet.write(1,36,KI_max)
                 worksheet.write(2,36,KI_min)
                 
                 worksheet.write(3,36,"Kia")
                 worksheet.write(4,36,KIa_max)
                 worksheet.write(5,36,KIa_min)
                 
                 worksheet.write(0,37,"kx")
                 worksheet.write(1,37,ksx_max)
                 worksheet.write(2,37,ksx_min)
                 worksheet.write(3,37,ksx_mean)
                 
                 worksheet.write(0,38,"ky")
                 worksheet.write(1,38,ksy_max)
                 worksheet.write(2,38,ksy_min)
                 worksheet.write(3,38,ksy_mean)
                 
                 worksheet.write(0,39,"kz")
                 worksheet.write(1,39,ksz_max)
                 worksheet.write(2,39,ksz_min)
                 worksheet.write(3,39,ksz_mean)
                 
                 worksheet.write(0,40,"kv")
                 worksheet.write(1,40,ksv_max)
                 worksheet.write(2,40,ksv_min)
                 worksheet.write(3,40,ksv_mean)
                 
                 worksheet.write(0,41,"km")
                 worksheet.write(1,41,km_max)
                 worksheet.write(2,41,km_min)
                 worksheet.write(3,41,km_mean)
                 
                 worksheet.write(0,42,"kE")
                 worksheet.write(1,42,kE_max)
                 worksheet.write(2,42,kE_min)
                 worksheet.write(3,42,kE_mean)
                 
                 worksheet.write(0,43,"kKI")
                 worksheet.write(1,43,kKI_max)
                 worksheet.write(2,43,kKI_min)
                 worksheet.write(3,43,kKI_mean)
                 
                 workbook.close()

             btn=warna_ganti(printa,text="printnya",command=printnya).grid(column=0,row=3)    
         
        lbl=Label(yuhu, text="RK4",width=15).grid(column=0,row=0)
        lbl=Label(yuhu, text="pilih output",width=15).grid(column=0,row=1)
        lbl=Label(yuhu, text="",width=15).grid(column=0,row=2)
         
        lbl=Label(yuhu, text="posisi",width=15).grid(column=0,row=3)
        lbl=Label(yuhu, text="x,t",width=15).grid(column=0,row=4)
        btn=warna_ganti(yuhu,text="_",command=gambarx,bg="yellow",fg="white").grid(column=1,row=4)
        lbl=Label(yuhu, text="y,t",width=15).grid(column=0,row=5)
        btn=warna_ganti(yuhu,text="_",command=gambary,bg="yellow",fg="white").grid(column=1,row=5)
        lbl=Label(yuhu, text="z,t",width=15).grid(column=0,row=6)
        btn=warna_ganti(yuhu,text="_",command=gambarz,bg="yellow",fg="white").grid(column=1,row=6)
         
        lbl=Label(yuhu, text="error xyz",width=15).grid(column=2,row=3)
        btn=warna_ganti(yuhu,text="_",command=gambarkx,bg="red",fg="white").grid(column=2,row=4)
        btn=warna_ganti(yuhu,text="_",command=gambarky,bg="red",fg="white").grid(column=2,row=5)
        btn=warna_ganti(yuhu,text="_",command=gambarkz,bg="red",fg="white").grid(column=2,row=6)
         
        lbl=Label(yuhu, text="kecepatan mutlak v",width=15).grid(column=0,row=7)
        lbl=Label(yuhu, text="v,t",width=15).grid(column=0,row=8)
        btn=warna_ganti(yuhu,text="_",command=gambarv,bg="yellow",fg="white").grid(column=1,row=8)
         
        lbl=Label(yuhu, text="error v",width=15).grid(column=2,row=7)
        btn=warna_ganti(yuhu,text="_",command=gambarkv,bg="red",fg="white").grid(column=2,row=8)
         
        lbl=Label(yuhu, text="massa gerak",width=15).grid(column=0,row=9)
        lbl=Label(yuhu, text="M,t",width=15).grid(column=0,row=10)
        btn=warna_ganti(yuhu,text="_",command=gambarm,bg="yellow",fg="white").grid(column=1,row=10)
         
        lbl=Label(yuhu, text="error m",width=15).grid(column=2,row=9)
        btn=warna_ganti(yuhu,text="_",command=gambarkm,bg="red",fg="white").grid(column=2,row=10)
         
        lbl=Label(yuhu, text="energi gerak ",width=15).grid(column=0,row=11)
        lbl=Label(yuhu, text="E,t",width=15).grid(column=0,row=12)
        btn=warna_ganti(yuhu,text="_",command=gambarE,bg="yellow",fg="white").grid(column=1,row=12)
         
        lbl=Label(yuhu, text="error E ",width=15).grid(column=2,row=11)
        btn=warna_ganti(yuhu,text="_",command=gambarkE,bg="red",fg="white").grid(column=2,row=12)

        lbl=Label(yuhu, text="energi kinetik ",width=15).grid(column=0,row=13)
        lbl=Label(yuhu, text="K,t",width=15).grid(column=0,row=14)
        btn=warna_ganti(yuhu,text="_",command=gambarKI,bg="yellow",fg="white").grid(column=1,row=14)
         
        lbl=Label(yuhu, text="error K ",width=15).grid(column=2,row=13)
        btn=warna_ganti(yuhu,text="_",command=gambarkKI,bg="red",fg="white").grid(column=2,row=14)
         
        lbl=Label(yuhu, text="Simulasi gerak ",width=15).grid(column=0,row=15)
        lbl=Label(yuhu, text="analitik ",width=15).grid(column=2,row=15)
        lbl=Label(yuhu, text="2D ",width=15).grid(column=0,row=16)
        btn=warna_ganti(yuhu,text="_",command=gambar2d,bg="yellow",fg="white").grid(column=1,row=16)
        btn=warna_ganti(yuhu,text="_",command=gambar2da,bg="yellow",fg="white").grid(column=2,row=16)
         
        lbl=Label(yuhu, text="3D ",width=15).grid(column=0,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar3d,bg="yellow",fg="white").grid(column=1,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar3da,bg="yellow",fg="white").grid(column=2,row=17)
         
        lbl=Label(yuhu, text="print xlsx ",width=15).grid(column=0,row=18)
        lbl=Label(yuhu, text="xlsx ",width=15).grid(column=0,row=19)
        btn=warna_ganti(yuhu,text="_",command=hasilkanxl,bg="yellow",fg="white").grid(column=1,row=19)
        yuhu.mainloop()      
        
    def rk45B():
        yuhu= Tk()
        yuhu.title("grafik yang diinginkan")
        yuhu.geometry("240x475+1110+0")
        
        def funcy (t , vy , y , k, B0, B1, y0, vx0 ):
             return vy
        def guncy (t , vy , y , k, B0, B1, y0, vx0 ):
             return -(k**2) * ((B0+B1*t)**2) *(y-y0)-((k) * ((B0+B1*t)) *vx0)
        def funcx (t , vx , x , k, B0, B1, x0, vy0  ):
             return vx
        def guncx (t , vx , x , k, B0, B1, x0, vy0  ):
             return ((k**2) * ((B0+B1*t)**2) )*(x0-x)+((k) * ((B0+B1*t)) *vy0)
        def funcz (t , vz , z , k, E0, E1):
             return vz
        def guncz (t , vz , z , k, E0, E1):
             return k*(E0+E1*t)
        for i in range (ndata-1):

           h=t[i+1]-t[i]
           k1y=h*funcy(t[i] ,vy[i],y[i],k, B0, B1, y0, vx0 )
           l1y=h*guncy(t[i] ,vy[i],y[i],k, B0, B1, y0, vx0 )
           k2y=h*funcy(t[i]+(1*h/4) , vy[i]+(1*l1y/4),y[i]+(1*k1y/4),k, B0, B1, y0, vx0  )
           l2y=h*guncy(t[i]+(1*h/4) , vy[i]+(1*l1y/4),y[i]+(1*k1y/4),k, B0, B1, y0, vx0  )               
           k3y=h*funcy(t[i]+(1*h/4) , vy[i]+(1*l1y/8)+(1*l2y/8) ,y[i]+(1*k1y/8)+(1*k2y/8),k, B0, B1, y0, vx0 )
           l3y=h*guncy(t[i]+(1*h/4) , vy[i]+(1*l1y/8)+(1*l2y/8) ,y[i]+(1*k1y/8)+(1*k2y/8),k, B0, B1, y0, vx0 )
           k4y=h*funcy(t[i]+(1*h/2) ,vy[i]-(1*l2y/2)+(l3y) ,y[i]-(1*k2y/2)+k3y,k, B0, B1, y0, vx0 )
           l4y=h*guncy(t[i]+(1*h/2) ,vy[i]-(1*l2y/2)+(l3y) ,y[i]-(1*k2y/2)+k3y,k, B0, B1, y0, vx0 )
           k5y=h*funcy(t[i]+(3*h/4) ,vy[i]+(3*l1y/16)+(9*l4y/16) ,y[i]+(3*k1y/16)+(9*k4y/16),k, B0, B1, y0, vx0 )
           l5y=h*guncy(t[i]+(3*h/4) ,vy[i]+(3*l1y/16)+(9*l4y/16) ,y[i]+(3*k1y/16)+(9*k4y/16),k, B0, B1, y0, vx0 )
           k6y=h*funcy(t[i]+h ,vy[i]-(3*l1y/7)+(2*l2y/7)+(12*l3y/7)-(12*l4y/7)+(8*l5y/7), y[i]-(3*k1y/7)+(2*k2y/7)+(12*k3y/7)-(12*k4y/7)+(8*k5y/7),k, B0, B1, y0, vx0 )
           l6y=h*guncy(t[i]+h ,vy[i]-(3*l1y/7)+(2*l2y/7)+(12*l3y/7)-(12*l4y/7)+(8*l5y/7), y[i]-(3*k1y/7)+(2*k2y/7)+(12*k3y/7)-(12*k4y/7)+(8*k5y/7),k, B0, B1, y0, vx0 ) 
              
               
           k1x=h*funcx(t[i] ,vx[i] ,x[i] ,k, B0, B1, x0, vy0 )
           l1x=h*guncx(t[i] ,vx[i] ,x[i] ,k, B0, B1, x0, vy0 )
           k2x=h*funcx(t[i]+(1*h/4) , vx[i]+(1*l1x/4) ,x[i]+(1*k1x/4),k, B0, B1, x0, vy0 )
           l2x=h*guncx(t[i]+(1*h/4) , vx[i]+(1*l1x/4) ,x[i]+(1*k1x/4),k, B0, B1, x0, vy0 )
           k3x=h*funcx(t[i]+(1*h/4) , vx[i]+(1*l1x/8)+(1*l2x/8) ,x[i]+(1*k1x/8)+(1*k2x/8),k, B0, B1, x0, vy0)
           l3x=h*guncx(t[i]+(1*h/4) , vx[i]+(1*l1x/8)+(1*l2x/8) ,x[i]+(1*k1x/8)+(1*k2x/8),k, B0, B1, x0, vy0)
           k4x=h*funcx(t[i]+(1*h/2) ,vx[i]-(1*l2x/2)+(l3x) ,x[i]-(1*k2x/2)+k3x,k, B0, B1, x0, vy0)
           l4x=h*guncx(t[i]+(1*h/2) ,vx[i]-(1*l2x/2)+(l3x) ,x[i]-(1*k2x/2)+k3x,k, B0, B1, x0, vy0)
           k5x=h*funcx(t[i]+(3*h/4) ,vx[i]+(3*l1x/16)+(9*l4x/16) ,x[i]+(3*k1x/16)+(9*k4x/16),k, B0, B1, x0, vy0)
           l5x=h*guncx(t[i]+(3*h/4) ,vx[i]+(3*l1x/16)+(9*l4x/16) ,x[i]+(3*k1x/16)+(9*k4x/16),k, B0, B1, x0, vy0)
           k6x=h*funcx(t[i]+h ,vx[i]-(3*l1x/7)+(2*l2x/7)+(12*l3x/7)-(12*l4x/7)+(8*l5x/7), x[i]-(3*k1x/7)+(2*k2x/7)+(12*k3x/7)-(12*k4x/7)+(8*k5x/7),k, B0, B1, x0, vy0)
           l6x=h*guncx(t[i]+h ,vx[i]-(3*l1x/7)+(2*l2x/7)+(12*l3x/7)-(12*l4x/7)+(8*l5x/7), x[i]-(3*k1x/7)+(2*k2x/7)+(12*k3x/7)-(12*k4x/7)+(8*k5x/7),k, B0, B1, x0, vy0) 

           
           k1z=h*funcz(t[i] ,vz[i],z[i],k, E0,E1)
           l1z=h*guncz(t[i] ,vz[i],z[i],k, E0,E1)
           k2z=h*funcz(t[i]+(1*h/4) , vz[i]+(1*l1z/4),z[i]+(1*k1z/4),k, E0,E1)
           l2z=h*guncz(t[i]+(1*h/4) , vz[i]+(1*l1z/4),z[i]+(1*k1z/4),k, E0,E1)
           k3z=h*funcz(t[i]+(1*h/4) , vz[i]+(1*l1z/8)+(1*l2z/8) ,z[i]+(1*k1z/8)+(1*k2z/8),k, E0,E1)
           l3z=h*guncz(t[i]+(1*h/4) , vz[i]+(1*l1z/8)+(1*l2z/8) ,z[i]+(1*k1z/8)+(1*k2z/8),k, E0,E1)
           k4z=h*funcz(t[i]+(1*h/2) ,vz[i]-(1*l2z/2)+(l3z) ,z[i]-(1*k2z/2)+k3z,k, E0,E1)
           l4z=h*guncz(t[i]+(1*h/2) ,vz[i]-(1*l2z/2)+(l3z) ,z[i]-(1*k2z/2)+k3z,k, E0,E1)
           k5z=h*funcz(t[i]+(3*h/4) ,vz[i]+(3*l1z/16)+(9*l4z/16) ,z[i]+(3*k1z/16)+(9*k4z/16),k, E0,E1)
           l5z=h*guncz(t[i]+(3*h/4) ,vz[i]+(3*l1z/16)+(9*l4z/16) ,z[i]+(3*k1z/16)+(9*k4z/16),k, E0,E1)
           k6z=h*funcz(t[i]+h ,vz[i]-(3*l1z/7)+(2*l2z/7)+(12*l3z/7)-(12*l4z/7)+(8*l5z/7), z[i]-(3*k1z/7)+(2*k2z/7)+(12*k3z/7)-(12*k4z/7)+(8*k5z/7),k, E0,E1)
           l6z=h*guncz(t[i]+h ,vz[i]-(3*l1z/7)+(2*l2z/7)+(12*l3z/7)-(12*l4z/7)+(8*l5z/7), z[i]-(3*k1z/7)+(2*k2z/7)+(12*k3z/7)-(12*k4z/7)+(8*k5z/7),k, E0,E1)
              
           kky=(7*k1y+32*k3y+12*k4y+32*k5y+7*k6y)/90
           ly=(7*l1y+32*l3y+12*l4y+32*l5y+7*l6y)/90
              
           y[i+1]=y[i]+kky
           vy[i+1]=vy[i]+ly
              
           kkx=(7*k1x+32*k3x+12*k4x+32*k5x+7*k6x)/90
           lx=(7*l1x+32*l3x+12*l4x+32*l5x+7*l6x)/90
          
           x[i+1]=x[i]+kkx
           vx[i+1]=vx[i]+lx
          
           kkz=(7*k1z+32*k3z+12*k4z+32*k5z+7*k6z)/90
           lz=(7*l1z+32*l3z+12*l4z+32*l5z+7*l6z)/90    
           
           z[i+1]=z[i]+kkz
           vz[i+1]=vz[i]+lz
           
           vya=-(vx0)*np.sin(k*B0*t)+(vy0)*np.cos(k*B0*t)
           ya=(vx0/(k*B0))*np.cos(k*B0*t)+(vy0/k)*np.sin(k*B0*t)+y0-(vx0/(k*B0))
           kesy=abs((ya-y)/ya)
           ksy=kesy*100
             
           vxa=(vy0)*np.sin((k*B0)*t)+(vx0)*np.cos((k*B0)*t)
           xa=(-vy0/(k*B0))*np.cos((k*B0)*t)+(vx0/(k*B0))*np.sin((k*B0)*t)+x0+vy0/(k*B0)
           kesx=abs((xa-x)/xa)
           ksx=kesx*100
             
           vza= k*(E0*t+E1*1/2* t**2)+vz0
           za= k*(E0*1/2*t**2+E1*1/6* t**3)+vz0*t+z0
           kesz=abs((za-z)/za)
           ksz=kesz*100
            
           v=np.sqrt(vx**2+vy**2+vz**2)
           va=np.sqrt(vxa**2+vya**2+vza**2)
           kesv=abs((va-v)/va)
           ksv=kesv*100
            
           g=np.sqrt(1-(v**2/c**2))
           ga=np.sqrt(1-(va**2/c**2))
            
           M=m/g
           Ma=m/ga
           kem=abs((M-Ma)/Ma)
           km=kem*100
            
           E=M*(c**2)
           Ea=Ma*(c**2)
           keE=abs((E-Ea)/Ea)
           kE=keE*100
            
           KI=(c**2)*M-(c**2)*m
           KIa=(c**2)*Ma-(c**2)*m
           keKI=abs((KI-KIa)/KIa)
           kKI=keKI*100
             
           x_max=np.max(x)
           x_min=np.min(x)
           
           y_max=np.max(y)
           y_min=np.min(y)
           
           z_max=np.max(z)
           z_min=np.min(z)
           
           xa_max=np.max(xa)
           xa_min=np.min(xa)
           
           ya_max=np.max(ya)
           ya_min=np.min(ya)
           
           za_max=np.max(za)
           za_min=np.min(za)
          
           v_max=np.max(v)
           v_min=np.min(v)
           
           M_max=np.max(M)
           M_min=np.min(M)
           
           E_max=np.max(E)
           E_min=np.min(E)
           
           KI_max=np.max(KI)
           KI_min=np.min(KI)
           
           va_max=np.max(va)
           va_min=np.min(va)
           
           Ma_max=np.max(Ma)
           Ma_min=np.min(Ma)
           
           Ea_max=np.max(Ea)
           Ea_min=np.min(Ea)
           
           KIa_max=np.max(KIa)
           KIa_min=np.min(KIa)
           
           ksy_max=np.max(ksy)
           ksy_min=np.min(ksy)
           ksy_mean=np.mean(ksy)
           
           ksx_max=np.max(ksx)
           ksx_min=np.min(ksx)
           ksx_mean=np.mean(ksx)
           
           ksz_max=np.max(ksz)
           ksz_min=np.min(ksz)
           ksz_mean=np.mean(ksz)
           
           ksv_max=np.max(ksv)
           ksv_min=np.min(ksv)
           ksv_mean=np.mean(ksv)
           
           km_max=np.max(km)
           km_min=np.min(km)
           km_mean=np.mean(km)
           
           kE_max=np.max(kE)
           kE_min=np.min(kE)
           kE_mean=np.mean(kE)
           
           kKI_max=np.max(kKI)
           kKI_min=np.min(kKI)
           kKI_mean=np.mean(kKI)
           
        def gambarx():
             root = tkinter.Tk()
             root.title('ploting posisi x, t pada RK45B ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,x,'g',label="posisi x terhadap waktu")
             plot1.set_ylabel('posisi x(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,xa,'y',label="posisi xa terhadap waktu")
             plot2.set_ylabel('posisi xa(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambary():
             root = tkinter.Tk()
             root.title('ploting posisi y, t pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,y,'g',label="posisi y terhadap waktu")
             plot1.set_ylabel('posisi y(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,ya,'y',label="posisi ya terhadap waktu")
             plot2.set_ylabel('posisi ya(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
        def gambarz():
             root = tkinter.Tk()
             root.title('ploting posisi z, t pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,z,'g',label="posisi z terhadap waktu")
             plot1.set_ylabel('posisi z(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,za,'y',label="posisi za terhadap waktu")
             plot2.set_ylabel('posisi za(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambarkx():
             root = tkinter.Tk()
             root.title('ploting kes.posisi x pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksx,'r',label="posisi kes.x terhadap waktu")
             plot1.set_ylabel('kesalahan kx(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot1.annotate(f'Min: {ksx_min}', xy=(t[np.argmin(ksx)], ksx_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksx_max}', xy=(t[np.argmax(ksx)], ksx_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksx_mean}', xy=(t[len(t) // 2], ksx_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))

             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambarky():
             root = tkinter.Tk()
             root.title('ploting kes.posisi y pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksy,'r',label="posisi kes.y terhadap waktu")
             plot1.set_ylabel('kesalahan ky(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot1.annotate(f'Min: {ksy_min}', xy=(t[np.argmin(ksy)], ksy_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksy_max}', xy=(t[np.argmax(ksy)], ksy_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksy_mean}', xy=(t[len(t) // 2], ksy_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
        def gambarkz():
             root = tkinter.Tk()
             root.title('ploting kes. posisi z pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksz,'r',label="posisi kes.z terhadap waktu")
             plot1.set_ylabel('kesalahan kz(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {ksz_min}', xy=(t[np.argmin(ksz)], ksz_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksz_max}', xy=(t[np.argmax(ksz)], ksz_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksz_mean}', xy=(t[len(t) // 2], ksz_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
              
        def gambarv():
             root = tkinter.Tk()
             root.title('ploting kecepatan v,T pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,v,'g',label="kecepatan mutlak terhadap waktu")
             plot1.set_ylabel('kecepata v (m,s)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {v_min}', xy=(t[np.argmin(v)], v_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {v_max}', xy=(t[np.argmax(v)], v_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,va,'y',label="kecepatan mutlak A terhadap waktu")
             plot2.set_ylabel('kecepata va (m,s)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {va_min}', xy=(t[np.argmin(va)], va_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {va_max}', xy=(t[np.argmax(va)], va_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarkv():
             root = tkinter.Tk()
             root.title('ploting kes.kecepatan v pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.set_ylabel('kesalahan kv(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.plot(t,ksv,'r',label="posis kes.v terhadap waktu")
             plot1.annotate(f'Min: {ksv_min}', xy=(t[np.argmin(ksv)], ksv_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksv_max}', xy=(t[np.argmax(ksv)], ksv_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksv_mean}', xy=(t[len(t) // 2], ksv_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarm():
             root = tkinter.Tk()
             root.title('ploting posisi m, t pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,M,'g',label="massa M terhadap waktu")
             plot1.set_ylabel('massa M(kg)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {M_min}', xy=(t[np.argmin(M)], M_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {M_max}', xy=(t[np.argmax(M)], M_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,Ma,'y',label="massa Ma terhadap waktu")
             plot2.set_ylabel('massa Ma(kg)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {Ma_min}', xy=(t[np.argmin(Ma)], Ma_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {Ma_max}', xy=(t[np.argmax(Ma)], Ma_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarkm():
             root = tkinter.Tk()
             root.title('ploting kes. m pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,km,'r',label="kesalaham massa terhadap waktu")
             plot1.set_ylabel('kesalahan km(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {km_min}', xy=(t[np.argmin(km)], km_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {km_max}', xy=(t[np.argmax(km)], km_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {km_mean}', xy=(t[len(t) // 2], km_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarE():
             root = tkinter.Tk()
             root.title('ploting energi E, t pada RK45B ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,E,'g',label="energi gerak")
             plot1.set_ylabel('energi E(j)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {E_min}', xy=(t[np.argmin(E)], E_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {E_max}', xy=(t[np.argmax(E)], E_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,Ea,'y',label="energi gerak analitik")
             plot2.set_ylabel('energi Ea(j)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {Ea_min}', xy=(t[np.argmin(Ea)], Ea_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {Ea_max}', xy=(t[np.argmax(Ea)], Ea_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def gambarkE():
             root = tkinter.Tk()
             root.title('ploting kes.E pada RK45B ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,kE,'r',label="posisi kes.E terhadap waktu")
             plot1.set_ylabel('kesalahan kE(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {kE_min}', xy=(t[np.argmin(kE)], kE_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {kE_max}', xy=(t[np.argmax(kE)], kE_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {kE_mean}', xy=(t[len(t) // 2], kE_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarKI():
             root = tkinter.Tk()
             root.title('ploting energi kinetik pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,KI,'g',label="energi kinetik ")
             plot1.set_ylabel('energi kinetik KI(j)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {KI_min}', xy=(t[np.argmin(KI)], KI_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {KI_max}', xy=(t[np.argmax(KI)], KI_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,KIa,'y',label="energi kinetik analitik")
             plot2.set_ylabel('energi kinetik KIa(j)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {KIa_min}', xy=(t[np.argmin(KIa)], KIa_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {KIa_max}', xy=(t[np.argmax(KIa)], KIa_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def gambarkKI():
             root = tkinter.Tk()
             root.title('ploting energi kinetik pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,kKI,'r',label="posisi kes energi kinetik terhadap waktu")
             plot1.set_ylabel('kesalahan KI(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {kKI_min}', xy=(t[np.argmin(kKI)], kKI_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {kKI_max}', xy=(t[np.argmax(kKI)], kKI_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {kKI_mean}', xy=(t[len(t) // 2], kKI_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambar2d():
             root = tkinter.Tk()
             root.title('ploting simulasi 2D pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(x,y,'r',label="posisi")
             plot1.set_xlabel('posisi x(m)')
             plot1.set_ylabel('posisi y(m)')
             plot1.plot(x[0],y[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(x[-1],y[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
         
        def gambar2da():
             root = tkinter.Tk()
             root.title('ploting simulasi 2D analitik pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(xa,ya,'b',label="posisi")
             plot1.set_xlabel('posisi xa(m)')
             plot1.set_ylabel('posisi ya(m)')
             plot1.plot(xa[0],ya[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(xa[-1],ya[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
             
        def gambar3d():
             root = tkinter.Tk()
             root.title('ploting simulasi 3d pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111,projection='3d')
             plot1.plot(x,y,z,'r',label="posisi")
             plot1.set_xlabel('posisi x(m)')
             plot1.set_ylabel('posisi y(m)')
             plot1.set_zlabel('posisi z(m)')
             plot1.plot(x[0],y[0],z[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(x[-1],y[-1],z[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
         
        def gambar3da():
             root = tkinter.Tk()
             root.title('ploting simulasi 3d analitik pada RK45B')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111,projection='3d')
             plot1.plot(xa,ya,za,'b',label="posisi")
             plot1.set_xlabel('posisi xa(m)')
             plot1.set_ylabel('posisi ya(m)')
             plot1.set_zlabel('posisi za(m)')
             plot1.plot(xa[0],ya[0],za[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(xa[-1],ya[-1],za[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def hasilkanxl():
             printa = Tk()
             printa.title("hasilprint pada RK45B")
             printa.geometry("240x100+1110+520")
             lbl=Label(printa, text="nama file    .xlsx",width=15)
             lbl.grid(column=0,row=0)
             namafile=Entry(printa, width=10)
             namafile.grid(column=1,row=0)
             def printnya():
                 workbook=xl.Workbook(namafile.get())
                 worksheet=workbook.add_worksheet()
                 worksheet.write(0,0,"q")
                 worksheet.write(1,0,q)
                 worksheet.write(0,1,"B0")
                 worksheet.write(1,1,B0)
                 worksheet.write(0,2,"B1")
                 worksheet.write(1,2,B1)
                 worksheet.write(0,3,"m")
                 worksheet.write(1,3,m)
                 worksheet.write(0,4,"E0")
                 worksheet.write(1,4,E0)
                 worksheet.write(0,5,"E1")
                 worksheet.write(1,5,E1)
                 worksheet.write(0,6,"t0")
                 worksheet.write(1,6,t0)
                 worksheet.write(0,7,"tn")
                 worksheet.write(1,7,tn)
                 worksheet.write(0,8,"ndata")
                 worksheet.write(1,8,ndata)
                 worksheet.write(0,9,"waktu")
                 worksheet.write_column(1,9,t)
                 worksheet.write(0,10,"x")
                 worksheet.write_column(1,10,x)
                 worksheet.write(0,11,"y")
                 worksheet.write_column(1,11,y)
                 worksheet.write(0,12,"z")
                 worksheet.write_column(1,12,z)
                 worksheet.write(0,13,"vx")
                 worksheet.write_column(1,13,vx)
                 worksheet.write(0,14,"vy")
                 worksheet.write_column(1,14,vy)
                 worksheet.write(0,15,"vz")
                 worksheet.write_column(1,15,vz)
                 worksheet.write(0,16,"v")
                 worksheet.write_column(1,16,v)
                 worksheet.write(0,17,"M")
                 worksheet.write_column(1,17,M)
                 worksheet.write(0,18,"E")
                 worksheet.write_column(1,18,E)
                 worksheet.write(0,19,"KI")
                 worksheet.write_column(1,19,KI)
                 worksheet.write(0,20,"xa")
                 worksheet.write_column(1,20,xa)
                 worksheet.write(0,21,"ya")
                 worksheet.write_column(1,21,ya)
                 worksheet.write(0,22,"za")
                 worksheet.write_column(1,22,za)
                 worksheet.write(0,23,"vxa")
                 worksheet.write_column(1,23,vxa)
                 worksheet.write(0,24,"vya")
                 worksheet.write_column(1,24,vya)
                 worksheet.write(0,25,"vza")
                 worksheet.write_column(1,25,vza)
                 worksheet.write(0,26,"va")
                 worksheet.write_column(1,26,va)
                 worksheet.write(0,27,"Ma")
                 worksheet.write_column(1,27,Ma)
                 worksheet.write(0,28,"Ea")
                 worksheet.write_column(1,28,Ea)
                 worksheet.write(0,29,"KIa")
                 worksheet.write_column(1,29,KIa)
                 
                 worksheet.write(0,30,"xRK45B")
                 worksheet.write(1,30,x_max)
                 worksheet.write(2,30,x_min)
                 
                 worksheet.write(3,30,"xA")
                 worksheet.write(4,30,xa_max)
                 worksheet.write(5,30,xa_min)
                 
                 worksheet.write(0,31,"yRK45B")
                 worksheet.write(1,31,y_max)
                 worksheet.write(2,31,y_min)
                 
                 worksheet.write(3,31,"yA")
                 worksheet.write(4,31,ya_max)
                 worksheet.write(5,31,ya_min)
                 
                 worksheet.write(0,32,"zRK45B")
                 worksheet.write(1,32,z_max)
                 worksheet.write(2,32,z_min)
                 
                 worksheet.write(3,32,"za")
                 worksheet.write(4,32,za_max)
                 worksheet.write(5,32,za_min)
                 
                 worksheet.write(0,33,"vRK45B")
                 worksheet.write(1,33,v_max)
                 worksheet.write(2,33,v_min)
                 
                 worksheet.write(3,33,"va")
                 worksheet.write(4,33,va_max)
                 worksheet.write(5,33,va_min)
                 
                 worksheet.write(0,34,"MRK45B")
                 worksheet.write(1,34,M_max)
                 worksheet.write(2,34,M_min)
                 
                 worksheet.write(3,34,"Ma")
                 worksheet.write(4,34,Ma_max)
                 worksheet.write(5,34,Ma_min)
                 
                 worksheet.write(0,35,"ERK45B")
                 worksheet.write(1,35,E_max)
                 worksheet.write(2,35,E_min)
                 
                 worksheet.write(3,35,"Ea")
                 worksheet.write(4,35,Ea_max)
                 worksheet.write(5,35,Ea_min)
                 
                 worksheet.write(0,36,"KIRK45B")
                 worksheet.write(1,36,KI_max)
                 worksheet.write(2,36,KI_min)
                 
                 worksheet.write(3,36,"Kia")
                 worksheet.write(4,36,KIa_max)
                 worksheet.write(5,36,KIa_min)
                 
                 worksheet.write(0,37,"kx")
                 worksheet.write(1,37,ksx_max)
                 worksheet.write(2,37,ksx_min)
                 worksheet.write(3,37,ksx_mean)
                 
                 worksheet.write(0,38,"ky")
                 worksheet.write(1,38,ksy_max)
                 worksheet.write(2,38,ksy_min)
                 worksheet.write(3,38,ksy_mean)
                 
                 worksheet.write(0,39,"kz")
                 worksheet.write(1,39,ksz_max)
                 worksheet.write(2,39,ksz_min)
                 worksheet.write(3,39,ksz_mean)
                 
                 worksheet.write(0,40,"kv")
                 worksheet.write(1,40,ksv_max)
                 worksheet.write(2,40,ksv_min)
                 worksheet.write(3,40,ksv_mean)
                 
                 worksheet.write(0,41,"km")
                 worksheet.write(1,41,km_max)
                 worksheet.write(2,41,km_min)
                 worksheet.write(3,41,km_mean)
                 
                 worksheet.write(0,42,"kE")
                 worksheet.write(1,42,kE_max)
                 worksheet.write(2,42,kE_min)
                 worksheet.write(3,42,kE_mean)
                 
                 worksheet.write(0,43,"kKI")
                 worksheet.write(1,43,kKI_max)
                 worksheet.write(2,43,kKI_min)
                 worksheet.write(3,43,kKI_mean)
                 
                 workbook.close()

             btn=warna_ganti(printa,text="printnya",command=printnya).grid(column=0,row=3)    
         
        lbl=Label(yuhu, text="RK45B",width=15).grid(column=0,row=0)
        lbl=Label(yuhu, text="pilih output",width=15).grid(column=0,row=1)
        lbl=Label(yuhu, text="",width=15).grid(column=0,row=2)
         
        lbl=Label(yuhu, text="posisi",width=15).grid(column=0,row=3)
        lbl=Label(yuhu, text="x,t",width=15).grid(column=0,row=4)
        btn=warna_ganti(yuhu,text="_",command=gambarx,bg="yellow",fg="white").grid(column=1,row=4)
        lbl=Label(yuhu, text="y,t",width=15).grid(column=0,row=5)
        btn=warna_ganti(yuhu,text="_",command=gambary,bg="yellow",fg="white").grid(column=1,row=5)
        lbl=Label(yuhu, text="z,t",width=15).grid(column=0,row=6)
        btn=warna_ganti(yuhu,text="_",command=gambarz,bg="yellow",fg="white").grid(column=1,row=6)
         
        lbl=Label(yuhu, text="error xyz",width=15).grid(column=2,row=3)
        btn=warna_ganti(yuhu,text="_",command=gambarkx,bg="red",fg="white").grid(column=2,row=4)
        btn=warna_ganti(yuhu,text="_",command=gambarky,bg="red",fg="white").grid(column=2,row=5)
        btn=warna_ganti(yuhu,text="_",command=gambarkz,bg="red",fg="white").grid(column=2,row=6)
         
        lbl=Label(yuhu, text="kecepatan mutlak v",width=15).grid(column=0,row=7)
        lbl=Label(yuhu, text="v,t",width=15).grid(column=0,row=8)
        btn=warna_ganti(yuhu,text="_",command=gambarv,bg="yellow",fg="white").grid(column=1,row=8)
         
        lbl=Label(yuhu, text="error v",width=15).grid(column=2,row=7)
        btn=warna_ganti(yuhu,text="_",command=gambarkv,bg="red",fg="white").grid(column=2,row=8)
         
        lbl=Label(yuhu, text="massa gerak",width=15).grid(column=0,row=9)
        lbl=Label(yuhu, text="M,t",width=15).grid(column=0,row=10)
        btn=warna_ganti(yuhu,text="_",command=gambarm,bg="yellow",fg="white").grid(column=1,row=10)
         
        lbl=Label(yuhu, text="error m",width=15).grid(column=2,row=9)
        btn=warna_ganti(yuhu,text="_",command=gambarkm,bg="red",fg="white").grid(column=2,row=10)
         
        lbl=Label(yuhu, text="energi gerak ",width=15).grid(column=0,row=11)
        lbl=Label(yuhu, text="E,t",width=15).grid(column=0,row=12)
        btn=warna_ganti(yuhu,text="_",command=gambarE,bg="yellow",fg="white").grid(column=1,row=12)
         
        lbl=Label(yuhu, text="error E ",width=15).grid(column=2,row=11)
        btn=warna_ganti(yuhu,text="_",command=gambarkE,bg="red",fg="white").grid(column=2,row=12)

        lbl=Label(yuhu, text="energi kinetik ",width=15).grid(column=0,row=13)
        lbl=Label(yuhu, text="K,t",width=15).grid(column=0,row=14)
        btn=warna_ganti(yuhu,text="_",command=gambarKI,bg="yellow",fg="white").grid(column=1,row=14)
         
        lbl=Label(yuhu, text="error K ",width=15).grid(column=2,row=13)
        btn=warna_ganti(yuhu,text="_",command=gambarkKI,bg="red",fg="white").grid(column=2,row=14)
         
        lbl=Label(yuhu, text="Simulasi gerak ",width=15).grid(column=0,row=15)
        lbl=Label(yuhu, text="analitik ",width=15).grid(column=2,row=15)
        lbl=Label(yuhu, text="2D ",width=15).grid(column=0,row=16)
        btn=warna_ganti(yuhu,text="_",command=gambar2d,bg="yellow",fg="white").grid(column=1,row=16)
        btn=warna_ganti(yuhu,text="_",command=gambar2da,bg="yellow",fg="white").grid(column=2,row=16)
         
        lbl=Label(yuhu, text="3D ",width=15).grid(column=0,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar3d,bg="yellow",fg="white").grid(column=1,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar3da,bg="yellow",fg="white").grid(column=2,row=17)
         
        lbl=Label(yuhu, text="print xlsx ",width=15).grid(column=0,row=18)
        lbl=Label(yuhu, text="xlsx ",width=15).grid(column=0,row=19)
        btn=warna_ganti(yuhu,text="_",command=hasilkanxl,bg="yellow",fg="white").grid(column=1,row=19)
        yuhu.mainloop()     
        
    def rk45F():
        yuhu= Tk()
        yuhu.title("agrafik yang diinginkan")
        yuhu.geometry("240x475+1110+0")
        def funcy (t , vy , y , k, B0, B1,  y0, vx0 ):
            return vy
        def guncy (t , vy , y , k, B0, B1,  y0, vx0 ):
            return -(k**2) * ((B0+B1*t)**2) *(y-y0)-((k) * ((B0+B1*t)) *vx0)
        def funcx (t , vx , x , k, B0, B1, x0, vy0  ):
            return vx
        def guncx (t , vx , x , k, B0, B1, x0, vy0  ):
            return ((k**2) * ((B0+B1*t)**2) )*(x0-x)+((k) * ((B0+B1*t)) *vy0)
        def funcz (t , vz , z , k, E0, E1):
            return vz
        def guncz (t , vz , z , k, E0, E1):
            return k*(E0+E1*t)
        for i in range (ndata-1):
                
            h=t[i+1]-t[i]
            
            k1y=h*funcy(t[i]     , vy[i]          ,y[i],k, B0, B1, y0, vx0 ) 
            l1y=h*guncy(t[i]     , vy[i]          ,y[i],k, B0, B1, y0, vx0 ) 
                
            k2y=h*funcy(t[i]+(1*h/4)      , vy[i]+(1*l1y/4)  ,y[i]+(1*k1y/4),k, B0, B1, y0, vx0 ) 
            l2y=h*guncy(t[i]+(1*h/4)      , vy[i]+(1*l1y/4)  ,y[i]+(1*k1y/4),k, B0, B1, y0, vx0  )
                
            k3y=h*funcy(t[i]+(3*h/8)      , vy[i]+(3*l1y/32)+(9*l2y/32)    ,y[i]+(3*k1y/32)+(9*k2y/32),k, B0, B1, y0, vx0  )
            l3y=h*guncy(t[i]+(3*h/8)      , vy[i]+(3*l1y/32)+(9*l2y/32)    ,y[i]+(3*k1y/32)+(9*k2y/32),k, B0, B1, y0, vx0  )
                
            k4y=h*funcy(t[i]+(12*h/13)    , vy[i]+(1932*l1y/2197)-(7200*l2y/2197)+(7296*l3y/2197)      ,y[i]+(1932*k1y/2197)-(7200*k2y/2197)+(7296*k3y/2197),k, B0, B1, y0, vx0  )
            l4y=h*guncy(t[i]+(12*h/13)    , vy[i]+(1932*l1y/2197)-(7200*l2y/2197)+(7296*l3y/2197)      ,y[i]+(1932*k1y/2197)-(7200*k2y/2197)+(7296*k3y/2197),k, B0, B1, y0, vx0  )
                
            k5y=h*funcy(t[i]+(h)   , vy[i]+(439*l1y/216)-(8*l2y)+(3680*l3y/513)-(845*l4y/4104)  ,y[i]+(439*k1y/216)-(8*k2y)+(3680*k3y/513)-(845*k4y/4104),k, B0, B1, y0, vx0  )
            l5y=h*guncy(t[i]+(h)   , vy[i]+(439*l1y/216)-(8*l2y)+(3680*l3y/513)-(845*l4y/4104)  ,y[i]+(439*k1y/216)-(8*k2y)+(3680*k3y/513)-(845*k4y/4104),k, B0, B1, y0, vx0  )
                
            k6y=h*funcy(t[i]+(1*h/2)      , vy[i]-(8*l1y/27)+(2*l2y)-(3544*l3y/2565)+(1859*l4y/4104)-(11*l5y/40)  ,y[i]-(8*k1y/27)+(2*k2y)-(3544*k3y/2565)+(1859*k4y/4104)-(11*k5y/40),k, B0, B1, y0, vx0 )
            l6y=h*guncy(t[i]+(1*h/2)      , vy[i]-(8*l1y/27)+(2*l2y)-(3544*l3y/2565)+(1859*l4y/4104)-(11*l5y/40)  ,y[i]-(8*k1y/27)+(2*k2y)-(3544*k3y/2565)+(1859*k4y/4104)-(11*k5y/40),k, B0, B1, y0, vx0 )   
                
            
            k1x=h*funcx(t[i]       , vx[i]            ,x[i],k, B0, B1, x0, vy0 )
            l1x=h*guncx(t[i]       , vx[i]            ,x[i],k, B0, B1, x0, vy0 )
                
            k2x=h*funcx(t[i]+(1*h/4)      , vx[i]+(1*l1x/4)         ,x[i]+(1*k1x/4),k, B0, B1, x0, vy0 )
            l2x=h*guncx(t[i]+(1*h/4)      , vx[i]+(1*l1x/4)         ,x[i]+(1*k1x/4),k, B0, B1, x0, vy0 )
                
            k3x=h*funcx(t[i]+(3*h/8)      , vx[i]+(3*l1x/32)+(9*l2x/32)    ,x[i]+(3*k1x/32)+(9*k2x/32),k, B0, B1, x0, vy0 )
            l3x=h*guncx(t[i]+(3*h/8)      , vx[i]+(3*l1x/32)+(9*l2x/32)    ,x[i]+(3*k1x/32)+(9*k2x/32),k, B0, B1, x0, vy0 )
                
            k4x=h*funcx(t[i]+(12*h/13)    , vx[i]+(1932*l1x/2197)-(7200*l2x/2197)+(7296*l3x/2197)      ,x[i]+(1932*k1x/2197)-(7200*k2x/2197)+(7296*k3x/2197),k, B0, B1, x0, vy0 )
            l4x=h*guncx(t[i]+(12*h/13)    , vx[i]+(1932*l1x/2197)-(7200*l2x/2197)+(7296*l3x/2197)      ,x[i]+(1932*k1x/2197)-(7200*k2x/2197)+(7296*k3x/2197),k, B0, B1, x0, vy0 )
                
            k5x=h*funcx(t[i]+(h)   , vx[i]+(439*l1x/216)-(8*l2x)+(3680*l3x/513)-(845*l4x/4104)  ,x[i]+(439*k1x/216)-(8*k2x)+(3680*k3x/513)-(845*k4x/4104),k, B0, B1, x0, vy0 )
            l5x=h*guncx(t[i]+(h)   , vx[i]+(439*l1x/216)-(8*l2x)+(3680*l3x/513)-(845*l4x/4104)  ,x[i]+(439*k1x/216)-(8*k2x)+(3680*k3x/513)-(845*k4x/4104),k, B0, B1, x0, vy0 )
                
            k6x=h*funcx(t[i]+(1*h/2)      , vx[i]-(8*l1x/27)+(2*l2x)-(3544*l3x/2565)+(1859*l4x/4104)-(11*l5x/40)  ,x[i]-(8*k1x/27)+(2*k2x)-(3544*k3x/2565)+(1859*k4x/4104)-(11*k5x/40),k, B0, B1, x0, vy0)
            l6x=h*guncx(t[i]+(1*h/2)      , vx[i]-(8*l1x/27)+(2*l2x)-(3544*l3x/2565)+(1859*l4x/4104)-(11*l5x/40)  ,x[i]-(8*k1x/27)+(2*k2x)-(3544*k3x/2565)+(1859*k4x/4104)-(11*k5x/40),k, B0, B1, x0, vy0)   

            k1z=h*funcz(t[i]       , vz[i]            ,z[i], k, E0, E1 )
            l1z=h*guncz(t[i]       , vz[i]            ,z[i], k, E0, E1 )
            
            
            k2z=h*funcz(t[i]+(1*h/4)      , vz[i]+(1*l1z/4)         ,z[i]+(1*k1z/4), k, E0, E1 )
            l2z=h*guncz(t[i]+(1*h/4)      , vz[i]+(1*l1z/4)         ,z[i]+(1*k1z/4) , k, E0, E1 )
            
            k3z=h*funcz(t[i]+(3*h/8)      , vz[i]+(3*l1z/32)+(9*l2z/32)    ,z[i]+(3*k1z/32)+(9*k2z/32) , k, E0, E1 )
            l3z=h*guncz(t[i]+(3*h/8)      , vz[i]+(3*l1z/32)+(9*l2z/32)    ,z[i]+(3*k1z/32)+(9*k2z/32) , k, E0, E1 )
            
            k4z=h*funcz(t[i]+(12*h/13)    , vz[i]+(1932*l1z/2197)-(7200*l2z/2197)+(7296*l3z/2197)      ,z[i]+(1932*k1z/2197)-(7200*k2z/2197)+(7296*k3z/2197) , k, E0, E1 )
            l4z=h*guncz(t[i]+(12*h/13)    , vz[i]+(1932*l1z/2197)-(7200*l2z/2197)+(7296*l3z/2197)      ,z[i]+(1932*k1z/2197)-(7200*k2z/2197)+(7296*k3z/2197) , k, E0, E1 )
                
            k5z=h*funcz(t[i]+(h)   , vz[i]+(439*l1z/216)-(8*l2z)+(3680*l3z/513)-(845*l4z/4104)  ,z[i]+(439*k1z/216)-(8*k2z)+(3680*k3z/513)-(845*k4z/4104) , k, E0, E1 )
            l5z=h*guncz(t[i]+(h)   , vz[i]+(439*l1z/216)-(8*l2z)+(3680*l3z/513)-(845*l4z/4104)  ,z[i]+(439*k1z/216)-(8*k2z)+(3680*k3z/513)-(845*k4z/4104) , k, E0, E1 )
                
            k6z=h*funcz(t[i]+(1*h/2)      , vz[i]-(8*l1z/27)+(2*l2z)-(3544*l3z/2565)+(1859*l4z/4104)-(11*l5z/40)  ,z[i]-(8*k1z/27)+(2*k2z)-(3544*k3z/2565)+(1859*k4z/4104)-(11*k5z/40), k, E0, E1 )
            l6z=h*guncz(t[i]+(1*h/2)      , vz[i]-(8*l1z/27)+(2*l2z)-(3544*l3z/2565)+(1859*l4z/4104)-(11*l5z/40)  ,z[i]-(8*k1z/27)+(2*k2z)-(3544*k3z/2565)+(1859*k4z/4104)-(11*k5z/40), k, E0, E1 )   
                     
            kky=(16*k1y/135+6656*k3y/12825+28561*k4y/56430-9*k5y/50+2*k6y/55)
            ly=(16*l1y/135+6656*l3y/12825+28561*l4y/56430-9*l5y/50+2*l6y/55)
                  
            y[i+1]=y[i]+kky
            vy[i+1]=vy[i]+ly
                  
            kkx=(16*k1x/135+6656*k3x/12825+28561*k4x/56430-9*k5x/50+2*k6x/55)
            lx=(16*l1x/135+6656*l3x/12825+28561*l4x/56430-9*l5x/50+2*l6x/55)
              
            x[i+1]=x[i]+kkx
            vx[i+1]=vx[i]+lx
              
            kkz=(16*k1z/135+6656*k3z/12825+28561*k4z/56430-9*k5z/50+2*k6z/55)
            lz=(16*l1z/135+6656*l3z/12825+28561*l4z/56430-9*l5z/50+2*l6z/55)
                  
            z[i+1]=z[i]+kkz
            vz[i+1]=vz[i]+lz
            
            vya=-(vx0)*np.sin(k*B0*t)+(vy0)*np.cos(k*B0*t)
            ya=(vx0/(k*B0))*np.cos(k*B0*t)+(vy0/k)*np.sin(k*B0*t)+y0-(vx0/(k*B0))
            kesy=abs((ya-y)/ya)
            ksy=kesy*100
              
            vxa=(vy0)*np.sin((k*B0)*t)+(vx0)*np.cos((k*B0)*t)
            xa=(-vy0/(k*B0))*np.cos((k*B0)*t)+(vx0/(k*B0))*np.sin((k*B0)*t)+x0+vy0/(k*B0)
            kesx=abs((xa-x)/xa)
            ksx=kesx*100
              
            vza=k*(E0*t+E1*1/2*t**2)+vz0
            za= k*(E0*1/2*t**2+E1*1/6* t**3)+vz0*t+z0
            kesz=abs((za-z)/za)
            ksz=kesz*100
             
            v=np.sqrt(vx**2+vy**2+vz**2)
            va=np.sqrt(vxa**2+vya**2+vza**2)
            kesv=abs((va-v)/va)
            ksv=kesv*100
             
            g=np.sqrt(1-(v**2/c**2))
            ga=np.sqrt(1-(va**2/c**2))
             
            M=m/g
            Ma=m/ga
            kem=abs((M-Ma)/Ma)
            km=kem*100
             
            E=M*(c**2)
            Ea=Ma*(c**2)
            keE=abs((E-Ea)/Ea)
            kE=keE*100
             
            KI=(c**2)*M-(c**2)*m
            KIa=(c**2)*Ma-(c**2)*m
            keKI=abs((KI-KIa)/KIa)
            kKI=keKI*100
            
            x_max=np.max(x)
            x_min=np.min(x)
            
            y_max=np.max(y)
            y_min=np.min(y)
            
            z_max=np.max(z)
            z_min=np.min(z)
            
            xa_max=np.max(xa)
            xa_min=np.min(xa)
            
            ya_max=np.max(ya)
            ya_min=np.min(ya)
            
            za_max=np.max(za)
            za_min=np.min(za)
            
            v_max=np.max(v)
            v_min=np.min(v)
            
            M_max=np.max(M)
            M_min=np.min(M)
            
            E_max=np.max(E)
            E_min=np.min(E)
            
            KI_max=np.max(KI)
            KI_min=np.min(KI)
            
            va_max=np.max(va)
            va_min=np.min(va)
            
            Ma_max=np.max(Ma)
            Ma_min=np.min(Ma)
            
            Ea_max=np.max(Ea)
            Ea_min=np.min(Ea)
            
            KIa_max=np.max(KIa)
            KIa_min=np.min(KIa)
            
            ksy_max=np.max(ksy)
            ksy_min=np.min(ksy)
            ksy_mean=np.mean(ksy)
            
            ksx_max=np.max(ksx)
            ksx_min=np.min(ksx)
            ksx_mean=np.mean(ksx)
            
            ksz_max=np.max(ksz)
            ksz_min=np.min(ksz)
            ksz_mean=np.mean(ksz)
            
            ksv_max=np.max(ksv)
            ksv_min=np.min(ksv)
            ksv_mean=np.mean(ksv)
            
            km_max=np.max(km)
            km_min=np.min(km)
            km_mean=np.mean(km)
            
            kE_max=np.max(kE)
            kE_min=np.min(kE)
            kE_mean=np.mean(kE)
            
            kKI_max=np.max(kKI)
            kKI_min=np.min(kKI)
            kKI_mean=np.mean(kKI)
            
        def gambarx():
             root = tkinter.Tk()
             root.title('ploting posisi x, t pada RK45F ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,x,'g',label="posisi x terhadap waktu")
             plot1.set_ylabel('posisi x(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,xa,'y',label="posisi xa terhadap waktu")
             plot2.set_ylabel('posisi xa(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambary():
             root = tkinter.Tk()
             root.title('ploting posisi y, t pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,y,'g',label="posisi y terhadap waktu")
             plot1.set_ylabel('posisi y(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,ya,'y',label="posisi ya terhadap waktu")
             plot2.set_ylabel('posisi ya(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
        def gambarz():
             root = tkinter.Tk()
             root.title('ploting posisi z, t pada RK45C')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,z,'g',label="posisi z terhadap waktu")
             plot1.set_ylabel('posisi z(m)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,za,'y',label="posisi za terhadap waktu")
             plot2.set_ylabel('posisi za(m)')
             plot2.set_xlabel('waktu t(s)')
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambarkx():
             root = tkinter.Tk()
             root.title('ploting kes.posisi x pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksx,'r',label="posisi kes.x terhadap waktu")
             plot1.set_ylabel('kesalahan kx(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot1.annotate(f'Min: {ksx_min}', xy=(t[np.argmin(ksx)], ksx_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksx_max}', xy=(t[np.argmax(ksx)], ksx_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksx_mean}', xy=(t[len(t) // 2], ksx_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))

             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
        def gambarky():
             root = tkinter.Tk()
             root.title('ploting kes.posisi y pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksy,'r',label="posisi kes.y terhadap waktu")
             plot1.set_ylabel('kesalahan ky(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.legend(loc="upper left")
             plot1.annotate(f'Min: {ksy_min}', xy=(t[np.argmin(ksy)], ksy_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksy_max}', xy=(t[np.argmax(ksy)], ksy_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksy_mean}', xy=(t[len(t) // 2], ksy_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
        def gambarkz():
             root = tkinter.Tk()
             root.title('ploting kes. posisi z pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,ksz,'r',label="posisi kes.z terhadap waktu")
             plot1.set_ylabel('kesalahan kz(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {ksz_min}', xy=(t[np.argmin(ksz)], ksz_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksz_max}', xy=(t[np.argmax(ksz)], ksz_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksz_mean}', xy=(t[len(t) // 2], ksz_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
              
        def gambarv():
             root = tkinter.Tk()
             root.title('ploting kecepatan v,T pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,v,'g',label="kecepatan mutlak terhadap waktu")
             plot1.set_ylabel('kecepata v (m,s)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {v_min}', xy=(t[np.argmin(v)], v_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {v_max}', xy=(t[np.argmax(v)], v_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,va,'y',label="kecepatan mutlak A terhadap waktu")
             plot2.set_ylabel('kecepata va (m,s)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {va_min}', xy=(t[np.argmin(va)], va_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {va_max}', xy=(t[np.argmax(va)], va_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarkv():
             root = tkinter.Tk()
             root.title('ploting kes.kecepatan v pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.set_ylabel('kesalahan kv(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.plot(t,ksv,'r',label="posis kes.v terhadap waktu")
             plot1.annotate(f'Min: {ksv_min}', xy=(t[np.argmin(ksv)], ksv_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {ksv_max}', xy=(t[np.argmax(ksv)], ksv_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {ksv_mean}', xy=(t[len(t) // 2], ksv_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarm():
             root = tkinter.Tk()
             root.title('ploting posisi m, t pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,M,'g',label="massa M terhadap waktu")
             plot1.set_ylabel('massa M(kg)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {M_min}', xy=(t[np.argmin(M)], M_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {M_max}', xy=(t[np.argmax(M)], M_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,Ma,'y',label="massa Ma terhadap waktu")
             plot2.set_ylabel('massa Ma(kg)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {Ma_min}', xy=(t[np.argmin(Ma)], Ma_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {Ma_max}', xy=(t[np.argmax(Ma)], Ma_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarkm():
             root = tkinter.Tk()
             root.title('ploting kes. m pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,km,'r',label="kesalaham massa terhadap waktu")
             plot1.set_ylabel('kesalahan km(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {km_min}', xy=(t[np.argmin(km)], km_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {km_max}', xy=(t[np.argmax(km)], km_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {km_mean}', xy=(t[len(t) // 2], km_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarE():
             root = tkinter.Tk()
             root.title('ploting energi E /t pada RK45F ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,E,'g',label="energi gerak")
             plot1.set_ylabel('energi E(j)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {E_min}', xy=(t[np.argmin(E)], E_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {E_max}', xy=(t[np.argmax(E)], E_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,Ea,'y',label="energi gerak analitik")
             plot2.set_ylabel('energi Ea(j)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {Ea_min}', xy=(t[np.argmin(Ea)], Ea_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {Ea_max}', xy=(t[np.argmax(Ea)], Ea_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def gambarkE():
             root = tkinter.Tk()
             root.title('ploting kes.E pada RK45F ')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,kE,'r',label="posisi kes.E terhadap waktu")
             plot1.set_ylabel('kesalahan kE(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {kE_min}', xy=(t[np.argmin(kE)], kE_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {kE_max}', xy=(t[np.argmax(kE)], kE_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {kE_mean}', xy=(t[len(t) // 2], kE_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambarKI():
             root = tkinter.Tk()
             root.title('ploting energi kinetik pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(211)
             plot1.plot(t,KI,'g',label="energi kinetik ")
             plot1.set_ylabel('energi kinetik KI(j)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {KI_min}', xy=(t[np.argmin(KI)], KI_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {KI_max}', xy=(t[np.argmax(KI)], KI_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             plot2 = fig.add_subplot(212)
             plot2.plot(t,KIa,'y',label="energi kinetik analitik")
             plot2.set_ylabel('energi kinetik KIa(j)')
             plot2.set_xlabel('waktu t(s)')
             plot2.annotate(f'Min: {KIa_min}', xy=(t[np.argmin(KIa)], KIa_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.annotate(f'Max: {KIa_max}', xy=(t[np.argmax(KIa)], KIa_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot2.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def gambarkKI():
             root = tkinter.Tk()
             root.title('ploting energi kinetik pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(t,kKI,'r',label="posisi kes energi kinetik terhadap waktu")
             plot1.set_ylabel('kesalahan KI(%)')
             plot1.set_xlabel('waktu t(s)')
             plot1.annotate(f'Min: {kKI_min}', xy=(t[np.argmin(kKI)], kKI_min), xytext=(10, 30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Max: {kKI_max}', xy=(t[np.argmax(kKI)], kKI_max), xytext=(-120, -30),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.annotate(f'Average: {kKI_mean}', xy=(t[len(t) // 2], kKI_mean), xytext=(10, 10),
                            textcoords='offset points', arrowprops=dict(arrowstyle="->"))
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
             canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()
             
        def gambar2d():
             root = tkinter.Tk()
             root.title('ploting simulasi 2D pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(x,y,'r',label="posisi")
             plot1.set_xlabel('posisi x(m)')
             plot1.set_ylabel('posisi y(m)')
             plot1.plot(x[0],y[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(x[-1],y[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
         
        def gambar2da():
             root = tkinter.Tk()
             root.title('ploting simulasi 2D analitik pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111)
             plot1.plot(xa,ya,'b',label="posisi")
             plot1.set_xlabel('posisi xa(m)')
             plot1.set_ylabel('posisi ya(m)')
             plot1.plot(xa[0],ya[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(xa[-1],ya[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
             
        def gambar3d():
             root = tkinter.Tk()
             root.title('ploting simulasi 3d pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111,projection='3d')
             plot1.plot(x,y,z,'r',label="posisi")
             plot1.set_xlabel('posisi x(m)')
             plot1.set_ylabel('posisi y(m)')
             plot1.set_zlabel('posisi z(m)')
             plot1.plot(x[0],y[0],z[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(x[-1],y[-1],z[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop()  
         
        def gambar3da():
             root = tkinter.Tk()
             root.title('ploting simulasi 3d analitik pada RK45F')
             root.geometry("850x700+220+0")
             fig = Figure(figsize = (6,3),dpi = 100)
             plot1 = fig.add_subplot(111,projection='3d')
             plot1.plot(xa,ya,za,'b',label="posisi")
             plot1.set_xlabel('posisi xa(m)')
             plot1.set_ylabel('posisi ya(m)')
             plot1.set_zlabel('posisi za(m)')
             plot1.plot(xa[0],ya[0],za[0],marker='o',markersize=8,color='yellow',label="stars")
             plot1.plot(xa[-1],ya[-1],za[-1],marker='o',markersize=8,color='green',label="end")
             plot1.legend(loc="upper left")
             canvas = FigureCanvasTkAgg(fig,master = root)
             canvas.draw()
             canvas.get_tk_widget().place(x=400,y=400)
             canvas.get_tk_widget().pack(side=tkinter.TOP,padx=10,pady=10, fill=tkinter.BOTH, expand=1)
             toolbar = NavigationToolbar2Tk(canvas,root)
             toolbar.update()
             def on_key_press(event):
                 print("you pressed {}".format(event.key))
                 key_press_handler(event, canvas, toolbar)
                 canvas.mpl_connect("key_press_event", on_key_press)
             def _quit():
                 root.quit()   
                 root.destroy()
             button = tkinter.Button(master=root, text="Quit", command=_quit)
             button.pack(side=tkinter.BOTTOM)
             tkinter.mainloop() 
             
        def hasilkanxl():
             printa = Tk()
             printa.title("hasilprint pada RK4F")
             printa.geometry("240x100+1110+520")
             lbl=Label(printa, text="nama file    .xlsx",width=15)
             lbl.grid(column=0,row=0)
             namafile=Entry(printa, width=10)
             namafile.grid(column=1,row=0)
             def printnya():
                 workbook=xl.Workbook(namafile.get())
                 worksheet=workbook.add_worksheet()
                 worksheet.write(0,0,"q")
                 worksheet.write(1,0,q)
                 worksheet.write(0,1,"B0")
                 worksheet.write(1,1,B0)
                 worksheet.write(0,2,"B1")
                 worksheet.write(1,2,B1)
                 worksheet.write(0,3,"m")
                 worksheet.write(1,3,m)
                 worksheet.write(0,4,"E0")
                 worksheet.write(1,4,E0)
                 worksheet.write(0,5,"E1")
                 worksheet.write(1,5,E1)
                 worksheet.write(0,6,"t0")
                 worksheet.write(1,6,t0)
                 worksheet.write(0,7,"tn")
                 worksheet.write(1,7,tn)
                 worksheet.write(0,8,"ndata")
                 worksheet.write(1,8,ndata)
                 worksheet.write(0,9,"waktu")
                 worksheet.write_column(1,9,t)
                 worksheet.write(0,10,"x")
                 worksheet.write_column(1,10,x)
                 worksheet.write(0,11,"y")
                 worksheet.write_column(1,11,y)
                 worksheet.write(0,12,"z")
                 worksheet.write_column(1,12,z)
                 worksheet.write(0,13,"vx")
                 worksheet.write_column(1,13,vx)
                 worksheet.write(0,14,"vy")
                 worksheet.write_column(1,14,vy)
                 worksheet.write(0,15,"vz")
                 worksheet.write_column(1,15,vz)
                 worksheet.write(0,16,"v")
                 worksheet.write_column(1,16,v)
                 worksheet.write(0,17,"M")
                 worksheet.write_column(1,17,M)
                 worksheet.write(0,18,"E")
                 worksheet.write_column(1,18,E)
                 worksheet.write(0,19,"KI")
                 worksheet.write_column(1,19,KI)
                 worksheet.write(0,20,"xa")
                 worksheet.write_column(1,20,xa)
                 worksheet.write(0,21,"ya")
                 worksheet.write_column(1,21,ya)
                 worksheet.write(0,22,"za")
                 worksheet.write_column(1,22,za)
                 worksheet.write(0,23,"vxa")
                 worksheet.write_column(1,23,vxa)
                 worksheet.write(0,24,"vya")
                 worksheet.write_column(1,24,vya)
                 worksheet.write(0,25,"vza")
                 worksheet.write_column(1,25,vza)
                 worksheet.write(0,26,"va")
                 worksheet.write_column(1,26,va)
                 worksheet.write(0,27,"Ma")
                 worksheet.write_column(1,27,Ma)
                 worksheet.write(0,28,"Ea")
                 worksheet.write_column(1,28,Ea)
                 worksheet.write(0,29,"KIa")
                 worksheet.write_column(1,29,KIa)
                 
                 worksheet.write(0,30,"xRK45F")
                 worksheet.write(1,30,x_max)
                 worksheet.write(2,30,x_min)
                 
                 worksheet.write(3,30,"xA")
                 worksheet.write(4,30,xa_max)
                 worksheet.write(5,30,xa_min)
                 
                 worksheet.write(0,31,"yRK45F")
                 worksheet.write(1,31,y_max)
                 worksheet.write(2,31,y_min)
                 
                 worksheet.write(3,31,"yA")
                 worksheet.write(4,31,ya_max)
                 worksheet.write(5,31,ya_min)
                 
                 worksheet.write(0,32,"zRK45F")
                 worksheet.write(1,32,z_max)
                 worksheet.write(2,32,z_min)
                 
                 worksheet.write(3,32,"za")
                 worksheet.write(4,32,za_max)
                 worksheet.write(5,32,za_min)
                 
                 worksheet.write(0,33,"vRK45F")
                 worksheet.write(1,33,v_max)
                 worksheet.write(2,33,v_min)
                 
                 worksheet.write(3,33,"va")
                 worksheet.write(4,33,va_max)
                 worksheet.write(5,33,va_min)
                 
                 worksheet.write(0,34,"MRK45F")
                 worksheet.write(1,34,M_max)
                 worksheet.write(2,34,M_min)
                 
                 worksheet.write(3,34,"Ma")
                 worksheet.write(4,34,Ma_max)
                 worksheet.write(5,34,Ma_min)
                 
                 worksheet.write(0,35,"ERK45F")
                 worksheet.write(1,35,E_max)
                 worksheet.write(2,35,E_min)
                 
                 worksheet.write(3,35,"Ea")
                 worksheet.write(4,35,Ea_max)
                 worksheet.write(5,35,Ea_min)
                 
                 worksheet.write(0,36,"KIRK45F")
                 worksheet.write(1,36,KI_max)
                 worksheet.write(2,36,KI_min)
                 
                 worksheet.write(3,36,"Kia")
                 worksheet.write(4,36,KIa_max)
                 worksheet.write(5,36,KIa_min)
                 
                 worksheet.write(0,37,"kx")
                 worksheet.write(1,37,ksx_max)
                 worksheet.write(2,37,ksx_min)
                 worksheet.write(3,37,ksx_mean)
                 
                 worksheet.write(0,38,"ky")
                 worksheet.write(1,38,ksy_max)
                 worksheet.write(2,38,ksy_min)
                 worksheet.write(3,38,ksy_mean)
                 
                 worksheet.write(0,39,"kz")
                 worksheet.write(1,39,ksz_max)
                 worksheet.write(2,39,ksz_min)
                 worksheet.write(3,39,ksz_mean)
                 
                 worksheet.write(0,40,"kv")
                 worksheet.write(1,40,ksv_max)
                 worksheet.write(2,40,ksv_min)
                 worksheet.write(3,40,ksv_mean)
                 
                 worksheet.write(0,41,"km")
                 worksheet.write(1,41,km_max)
                 worksheet.write(2,41,km_min)
                 worksheet.write(3,41,km_mean)
                 
                 worksheet.write(0,42,"kE")
                 worksheet.write(1,42,kE_max)
                 worksheet.write(2,42,kE_min)
                 worksheet.write(3,42,kE_mean)
                 
                 worksheet.write(0,43,"kKI")
                 worksheet.write(1,43,kKI_max)
                 worksheet.write(2,43,kKI_min)
                 worksheet.write(3,43,kKI_mean)
                 
                 workbook.close()

             btn=warna_ganti(printa,text="printnya",command=printnya).grid(column=0,row=3)    
         
        lbl=Label(yuhu, text="RK45f",width=15).grid(column=0,row=0)
        lbl=Label(yuhu, text="pilih output",width=15).grid(column=0,row=1)
        lbl=Label(yuhu, text="",width=15).grid(column=0,row=2)
         
        lbl=Label(yuhu, text="posisi",width=15).grid(column=0,row=3)
        lbl=Label(yuhu, text="x,t",width=15).grid(column=0,row=4)
        btn=warna_ganti(yuhu,text="_",command=gambarx,bg="yellow",fg="white").grid(column=1,row=4)
        lbl=Label(yuhu, text="y,t",width=15).grid(column=0,row=5)
        btn=warna_ganti(yuhu,text="_",command=gambary,bg="yellow",fg="white").grid(column=1,row=5)
        lbl=Label(yuhu, text="z,t",width=15).grid(column=0,row=6)
        btn=warna_ganti(yuhu,text="_",command=gambarz,bg="yellow",fg="white").grid(column=1,row=6)
         
        lbl=Label(yuhu, text="error xyz",width=15).grid(column=2,row=3)
        btn=warna_ganti(yuhu,text="_",command=gambarkx,bg="red",fg="white").grid(column=2,row=4)
        btn=warna_ganti(yuhu,text="_",command=gambarky,bg="red",fg="white").grid(column=2,row=5)
        btn=warna_ganti(yuhu,text="_",command=gambarkz,bg="red",fg="white").grid(column=2,row=6)
         
        lbl=Label(yuhu, text="kecepatan mutlak v",width=15).grid(column=0,row=7)
        lbl=Label(yuhu, text="v,t",width=15).grid(column=0,row=8)
        btn=warna_ganti(yuhu,text="_",command=gambarv,bg="yellow",fg="white").grid(column=1,row=8)
         
        lbl=Label(yuhu, text="error v",width=15).grid(column=2,row=7)
        btn=warna_ganti(yuhu,text="_",command=gambarkv,bg="red",fg="white").grid(column=2,row=8)
         
        lbl=Label(yuhu, text="massa gerak",width=15).grid(column=0,row=9)
        lbl=Label(yuhu, text="M,t",width=15).grid(column=0,row=10)
        btn=warna_ganti(yuhu,text="_",command=gambarm,bg="yellow",fg="white").grid(column=1,row=10)
         
        lbl=Label(yuhu, text="error m",width=15).grid(column=2,row=9)
        btn=warna_ganti(yuhu,text="_",command=gambarkm,bg="red",fg="white").grid(column=2,row=10)
         
        lbl=Label(yuhu, text="energi gerak ",width=15).grid(column=0,row=11)
        lbl=Label(yuhu, text="E,t",width=15).grid(column=0,row=12)
        btn=warna_ganti(yuhu,text="_",command=gambarE,bg="yellow",fg="white").grid(column=1,row=12)
         
        lbl=Label(yuhu, text="error E ",width=15).grid(column=2,row=11)
        btn=warna_ganti(yuhu,text="_",command=gambarkE,bg="red",fg="white").grid(column=2,row=12)

        lbl=Label(yuhu, text="energi kinetik ",width=15).grid(column=0,row=13)
        lbl=Label(yuhu, text="K,t",width=15).grid(column=0,row=14)
        btn=warna_ganti(yuhu,text="_",command=gambarKI,bg="yellow",fg="white").grid(column=1,row=14)
         
        lbl=Label(yuhu, text="error K ",width=15).grid(column=2,row=13)
        btn=warna_ganti(yuhu,text="_",command=gambarkKI,bg="red",fg="white").grid(column=2,row=14)
         
        lbl=Label(yuhu, text="Simulasi gerak ",width=15).grid(column=0,row=15)
        lbl=Label(yuhu, text="analitik ",width=15).grid(column=2,row=15)
        lbl=Label(yuhu, text="2D ",width=15).grid(column=0,row=16)
        btn=warna_ganti(yuhu,text="_",command=gambar2d,bg="yellow",fg="white").grid(column=1,row=16)
        btn=warna_ganti(yuhu,text="_",command=gambar2da,bg="yellow",fg="white").grid(column=2,row=16)
         
        lbl=Label(yuhu, text="3D ",width=15).grid(column=0,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar3d,bg="yellow",fg="white").grid(column=1,row=17)
        btn=warna_ganti(yuhu,text="_",command=gambar3da,bg="yellow",fg="white").grid(column=2,row=17)
         
        lbl=Label(yuhu, text="print xlsx ",width=15).grid(column=0,row=18)
        lbl=Label(yuhu, text="xlsx ",width=15).grid(column=0,row=19)
        btn=warna_ganti(yuhu,text="_",command=hasilkanxl,bg="yellow",fg="white").grid(column=1,row=19)
        yuhu.mainloop()           
        
    lbl22=Label(linier, text="per.linier ",width=15).grid(column=0,row=1)
    lbl23=Label(linier, text="pilih metode ",width=10).grid(column=0,row=2)
    lbl24=Label(linier, text=" ",width=10).grid(column=0,row=3)
    lbl25=Label(linier, text="rk4 ",width=10).grid(column=0,row=4)
    lbl26=Label(linier, text="rk45B ",width=10).grid(column=0,row=5)
    lbl27=Label(linier, text="rk45F ",width=10).grid(column=0,row=6)
        
    btn66=warna_ganti(linier,text="_",command=rk4,bg="pink",fg="white").grid(column=1,row=4)
    btn77=warna_ganti(linier,text="_",command=rk45B,bg="pink",fg="white").grid(column=1,row=5)
    btn88=warna_ganti(linier,text="_",command=rk45F,bg="pink",fg="white").grid(column=1,row=6)
    linier.mainloop()    
def imputan():        
    user_input=aB1.get()
    if user_input=="0":
        linier()

    else:
        nonlinier()
        
main = Tk()
main.title("input masukan")
main.geometry("200x430+0+0")

lbl=Label(main,text="masukan input",width=15).grid(column=0,row=0)
lbl=Label(main,text="10^",width=10).grid(column=2,row=0)

lbl=Label(main, text="q",width=5)
lbl.grid(column=0,row=1)

aq=Entry(main, width=5)
aq.grid(column=1,row=1)

aqT=Entry(main, width=5)
aqT.grid(column=2,row=1)

lbl1=Label(main, text="B0",width=5)
lbl1.grid(column=0,row=2)

aB0=Entry(main, width=5)
aB0.grid(column=1,row=2)

aB0T=Entry(main, width=5)
aB0T.grid(column=2,row=2)
    
lbl2=Label(main, text="B1",width=5)
lbl2.grid(column=0,row=3)

aB1=Entry(main, width=5)
aB1.grid(column=1,row=3)

aB1T=Entry(main, width=5)
aB1T.grid(column=2,row=3)

lbl3=Label(main, text="m",width=5)
lbl3.grid(column=0,row=4)

am=Entry(main, width=5)
am.grid(column=1,row=4)

amT=Entry(main, width=5)
amT.grid(column=2,row=4)

lbl4=Label(main, text="E0",width=5)
lbl4.grid(column=0,row=5)

aE0=Entry(main, width=5)
aE0.grid(column=1,row=5)

aE0T=Entry(main, width=5)
aE0T.grid(column=2,row=5)

lbl5=Label(main, text="E1",width=5)
lbl5.grid(column=0,row=6)

aE1=Entry(main, width=5)
aE1.grid(column=1,row=6)

aE1T=Entry(main, width=5)
aE1T.grid(column=2,row=6)

#imputan untuk waktu

lbl6=Label(main, text="t0",width=5)
lbl6.grid(column=0,row=7)

at0=Entry(main, width=5)
at0.grid(column=1,row=7)

at0T=Entry(main, width=5)
at0T.grid(column=2,row=7)

lbl7=Label(main, text="tn",width=5)
lbl7.grid(column=0,row=8)

atn=Entry(main, width=5)
atn.grid(column=1,row=8)

atnT=Entry(main, width=5)
atnT.grid(column=2,row=8)

lbl8=Label(main, text="ndata",width=5)
lbl8.grid(column=0,row=9)

andata=Entry(main, width=5)
andata.grid(column=1,row=9)

andataT=Entry(main, width=5)
andataT.grid(column=2,row=9)

#imputan untuk posisi awal

lbl9=Label(main, text="x0",width=5)
lbl9.grid(column=0,row=10)

ax0=Entry(main, width=5)
ax0.grid(column=1,row=10)

ax0T=Entry(main, width=5)
ax0T.grid(column=2,row=10)

lbl10=Label(main, text="y0",width=5)
lbl10.grid(column=0,row=11)

ay0=Entry(main, width=5)
ay0.grid(column=1,row=11)

ay0T=Entry(main, width=5)
ay0T.grid(column=2,row=11)

lbl11=Label(main, text="z0",width=5)
lbl11.grid(column=0,row=12)

az0=Entry(main, width=5)
az0.grid(column=1,row=12)

az0T=Entry(main, width=5)
az0T.grid(column=2,row=12)

#imputan untuk kecepatan awal
lbl12=Label(main, text="vx0",width=5)
lbl12.grid(column=0,row=13)

avx0=Entry(main, width=5)
avx0.grid(column=1,row=13)

avx0T=Entry(main, width=5)
avx0T.grid(column=2,row=13)

lbl13=Label(main, text="vy0",width=5)
lbl13.grid(column=0,row=14)

avy0=Entry(main, width=5)
avy0.grid(column=1,row=14)

avy0T=Entry(main, width=5)
avy0T.grid(column=2,row=14)

lbl14=Label(main, text="vz0",width=5)
lbl14.grid(column=0,row=15)

avz0=Entry(main, width=5)
avz0.grid(column=1,row=15)

avz0T=Entry(main, width=5)
avz0T.grid(column=2,row=15)

btn=warna_ganti(main,text="input",bg="blue",fg="white",command=imputan,width=10).grid(column=0,row=16)
btn=warna_ganti(main,text="save", bg="blue",fg="white",command=save_data,width=10).grid(column=0,row=17)
btn=warna_ganti(main,text="load", bg="blue",fg="white",command=load_data,width=10).grid(column=0,row=18)
main.mainloop()