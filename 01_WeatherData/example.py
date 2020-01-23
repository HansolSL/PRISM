# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 20:23:08 2020

@author: KJM
"""

import numpy as np
from SolarDecomposition import CalSolar
from Packages import SolarDecomposition as SD

# Global Radiation
Gl = 4.166667  # W/m2

# Time
Time = [9, 18, 0]  # Day of Year, Hour, Minute

# Class
sd = SD.CalSolar()

# Decomposition
Ib, Id, Kt, Sol_Altit, Sol_Azim = sd.Watanabe(Gl, Time)
