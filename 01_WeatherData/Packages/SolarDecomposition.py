# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 19:57:20 2020

@author: KJM
"""

import numpy as np

class CalSolar():
    def __init__(self):
        self.Lo_Longit = 126.26
        self.Lo_Latit = 37.48
        self.St_Longit = 135
    
    def SolarTime(self):
        Correction = 4 * (self.Lo_Longit - self.St_Longit) + self.EOT
        Sol_TimeHour = self.St_TimeHour
        Sol_TimeMin = self.St_TimeMin + Correction
        
        if Sol_TimeMin >= 60:
            Sol_TimeHour += 1
            Sol_TimeMin -= 60
    
        elif Sol_TimeMin < 0:
            Sol_TimeHour -= 1
            Sol_TimeMin += 60
        
        return ((Sol_TimeHour + Sol_TimeMin / 60) - 12) * 15
    
    
    def SolarLocation(self, N_Time):
        
        self.DayOfYear = N_Time[0]
        
        # Solar Time
        self.St_TimeHour = N_Time[1]
        if self.St_TimeHour == 0:
            self.St_TimeHour = 24
    
        self.St_TimeMin = N_Time[2]
    
        # Solar Declination
        Sol_Declin = 23.45 * np.sin(2*np.pi * (self.DayOfYear+284)/365)
        
        # Equation of Time
        B = (self.DayOfYear-1)*np.pi*2/365
        self.EOT = 229.2*(0.000075 + 0.001868*np.cos(B) - 0.032077*np.sin(B)
                     - 0.014615*np.cos(2*B) - 0.04089*np.sin(2*B))
        
        # Hour Angle
        HourAngle = self.SolarTime()
        
        # Solar Altitude
        Sol_Altit = 180/np.pi * np.arcsin(np.cos(np.pi/180 * self.Lo_Latit) * np.cos(np.pi/180 * Sol_Declin)
                    * np.cos(np.pi/180 * HourAngle) + np.sin(np.pi/180 * self.Lo_Latit) * np.sin(np.pi/180 * Sol_Declin))
        
        # Solar Azimuth
        Sol_Azim = 180/np.pi * np.arccos((np.sin(np.pi/180 * Sol_Altit) * np.sin(np.pi/180 * self.Lo_Latit) - np.sin(np.pi/180 * Sol_Declin))
                   / (np.cos(np.pi/180 * Sol_Altit) * np.cos(np.pi/180 * self.Lo_Latit)))
                   
        if self.St_TimeHour < 13:
            Sol_Azim *= -1
        
        return Sol_Altit, Sol_Azim
    
    
    def Watanabe(self, Gl_Radiation, N_Time):
        
        # Solar Location
        Sol_Altit, Sol_Azim = self.SolarLocation(N_Time)
        
        # extraterrestrial irradiation
        #Io = 1382 * (1 + 0.033 * np.cos(2*np.pi *self.DayOfYear /365))
        Io = 1367
        
        # CLearness Coefficient
        Kt = Gl_Radiation / (Io * np.sin(np.pi/180 * Sol_Altit))
        shy_3 = Io * np.sin(np.pi/180 * Sol_Altit)
        
        if shy_3 < 0:
            Kt = 0
        #elif shy_3 < Io * np.sin(np.pi/180 * 10):
        #    Kt = Gl_Radiation / Io * np.sin(np.pi/180 * 10)
            
        Ktc = 0.4268 + 0.1934 * np.sin(np.pi/180 * Sol_Altit)
        
        if Kt >= Ktc:
            Kds = Kt - (1.107 + 0.03569 * np.sin(np.pi/180 * Sol_Altit) + 1.681 * np.power(np.sin(np.pi/180 * Sol_Altit), 2)) * np.power((1-Kt),3)
            
        else:
            Kds = (3.996 - 3.862 * np.sin(np.pi/180 * Sol_Altit) + 1.54 * np.power(np.sin(np.pi/180 * Sol_Altit), 2)) * np.power(Kt, 3)
            
        Ib = Io * np.sin(np.pi/180 * Sol_Altit) * Kds * (1-Kt) / (1 - Kds)
        Id = Io * np.sin(np.pi/180 * Sol_Altit) * (Kt - Kds) /  (1 - Kds)
        
        return Ib, Id, Kt, Sol_Altit, Sol_Azim
        
        
        
        
        
        
        
        

