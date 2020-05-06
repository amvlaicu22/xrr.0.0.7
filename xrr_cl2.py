#!/usr/bin/python3
# -*- coding: utf-8 -*-
import time
from XRR import XRR  #, _XLayer

def fit_layer(Xrr):
    Xrr.fitprof = 0    # 0=layer,  1=profile
    Xrr.fitmode = 0    # 0=LM ,    1=DE
    Xrr.fitmax = 500
    Xrr.update = 100
    Xrr.fit_set_fitkeys()  #   d, rho, sigma, a0b0
    Xrr.fit_start('console start and sleep(4)')
    # ~ time.sleep(10)
    # ~ Xrr.fit_stop('console stop')
    Xrr.fit_wait()

    Xrr.layer_print()
    Xrr.layer_profile()
    Xrr._plot()
    print("finished fit layer ... ")


def fit_profile(Xrr):
    Xrr.fitprof = 1    # 0=layer,  1=profile
    Xrr.fitmode = 0    # 0=LM ,    1=DE
    Xrr.fitmax = 400
    Xrr.update = 50
    Xrr._pdz = 2
    # ~ Xrr.layer_profile()  # >> is called in fit profile
    Xrr.fit_start('fit profile ...')
    # ~ time.sleep(10)
    # ~ Xrr.fit_stop('console stop')
    Xrr.fit_wait()

    Xrr._plot()
    print("finished fit profile ... ")
    
    

if __name__ == "__main__":
    Xrr = XRR()
    # ~ Xrr.data_load_gui(".dat")
    Xrr._plot()
    
    # ~ input("Press [CR] to start fit layer")
    fit_layer(Xrr)
    
    # ~ input("Press [CR] to start fit layer")
    fit_profile(Xrr)
    input("Press [CR] to exit program")
    

