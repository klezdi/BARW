#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BRANCHING AND ANNIHILATING RANDOM WALK SIMULATION MODULES

Script to run simulations of branching and annihilating random walks based on 
    https://www.nature.com/articles/s41467-021-27135-5 

Uses the modules defined in the script "branching_rules.py"

@author: mucar
"""
import numpy as np
import branching_rules as br

#%%
# we first set the elementary step size to be equal to 1:
lstep = 1

class Tissue:
    # these attributes will be defined at the first instance of the class (i.e. whenever we generate an object for it)
    def __init__(self, prob, tip, angle, min_angle = np.pi/10, rad_avoid = 3*lstep, rad_termin = 1.5*lstep):
        """
        Initialize parameters

        Parameters
        ----------
        prob : float
            DESCRIPTION.
        tip : array
            DESCRIPTION.
        angle : array
            DESCRIPTION.
        min_angle : float, optional
            DESCRIPTION. The default is np.pi/10.
        rad_avoid : float, optional
            DESCRIPTION. The default is 3*lstep.
        rad_termin : float, optional
            DESCRIPTION. The default is 1.5*lstep.

        Returns
        -------
        None.

        """
        self.tip = tip
        self.angle = angle
        self.min_angle = min_angle
        self.rad_avoid = rad_avoid
        self.rad_termin = rad_termin
        self.prob = prob
        
    def evolve(self, coords):
        """
        Method for branching + elongation + annihilation of tips.

        Parameters
        ----------
        coords : array
            List of ALL coordinate points, together with branch and parent labels for the entire network.
        Returns
        -------
        None.

        """
        
        updated_info = br.branching(self.prob, self.tip, self.angle, coords, self.min_angle, self.rad_termin)
        
        # update the coordiantes and angles of all active tips:
        self.tip = updated_info['tip']
        self.angle = updated_info['angle']
        
    # this method will incorporate the external guidance, self-avoidance and hetero-interactions between different tissues!
    def guidance(self, coord_last, coords_till, fc, fs):
        """
        Method for the effect of external potential and self-avoidance on the tips.

        Parameters
        ----------
        coord_last : array
            List of all active coordinate points at the last time step of the loop.
        coords_till : array
            List of all coordinate points until a N time steps before (to prevent immediate annihilation events).
        fc : float
            Strength of the external guidance/potential (typically between 0 and 0.3 is a good range).
        fav : float
            Strength of self-avoidance / attraction of active tips of the network to its branches 
            (same range as fs). fs<0 for self-repulsion, fs>0 for self-attraction
        Returns
        -------
        None.

        """
        updated_info = br.guidance_avoidance(self.tip, self.angle, coord_last, coords_till, fc, fs, self.rad_avoid)
        
        # update the coordiantes and angles of all active tips:
        self.tip = updated_info['tip']
        self.angle = updated_info['angle']
        
    def __repr__(self):
        return '[Tissue Tips: %s, Angles: %s]' % (self.tip, self.angle)
    

#%%

def simulation_loop(prob,fc,fs,tmax):
    """
    Simulation loop for the BARWs with external guidance and self-interactions.

    Parameters
    ----------
    prob : float
        Branching probability of the active tips.
    fc : float
        Strength of the external guidance/potential (typically between 0 and 0.3 is a good range).
    fs : float
        Strength of self-avoidance / attraction of active tips of the network to its branches 
        (same range as fs). fs<0 for self-repulsion, fs>0 for self-attraction
    tmax : int
        Maximal simulation time.

    Returns
    -------
    d : dict
        List of all coordinates (coordinates) & angle values (angles) 
        (as well as their parent and branch labels), also a list (evolve) 
        containing info on number of active tips at each simulation step 
        (can be used to analyze network statistics at intermediate time points)

    """
    # starting coordinates of the network:
    start_x, start_y = 100.0, 100.0

    # here we generate an instances of the Tissue class:
    neuron = Tissue(prob, tip=np.array([[start_x,start_y,0,1]]), angle = np.array([np.pi/2]), rad_avoid=3)
    
    # these lists will carry the info that we need to update at every simulation step:
    coordinates = np.array([neuron.tip[0]])    
    coords_last = np.array([])
    coords_until = coordinates
    angle_list = np.array([[np.pi/2,0]])
    last_tiplength_fin = np.array([len(neuron.tip)])
    
    
    for t in range(tmax):
        # evaluate as long as there's an active tip in neuron 1:
        if len(neuron.tip)>0:
            neuron.evolve(coordinates)
            neuron.guidance(coords_last,coords_until,fc,fs)
    
            # set angle values to be within [-pi,pi]
            angle = (neuron.angle+np.pi) % (2*np.pi) - np.pi     
            # save the angles of nodes [in degrees!] including the generation number
            angle_list = np.append(angle_list,np.column_stack((np.degrees(angle),neuron.tip[:,-1])),axis=0)
    
            # save the coordinates of all nodes
            coordinates = np.append(coordinates,np.array(neuron.tip),axis=0)    
            last_tiplength_fin = np.append(last_tiplength_fin,len(neuron.tip))
            coords_last = coordinates[-int(np.sum(last_tiplength_fin[-1:])):]
            if t>2:
                coords_until = coordinates[:-int(np.sum(last_tiplength_fin[-2:]))]
                
        # if there are no active tips overall, stop the simulation:
        if len(neuron.tip)==0:
            break
                    
    d = dict()
    d['coordinates'] = coordinates
    d['angles'] = angle_list
    d['evolve'] = last_tiplength_fin
    
    return d
    
