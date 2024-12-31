# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:02:16 2024

@author: coolm
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from multiprocessing import Pool, freeze_support

class Binary_Planet_System():
    
    def __init__(self,star_masses = [1, 1]):
        """
        Intiallises the motion of the two stars

        Parameters
        ----------
        star_masses : LIST, optional
            The Masses of the two stars, intially set to be the mass of the sun
            in units of Solar Mass.For simplisty the larger mass should be first.
            The default is [1, 1].

        Returns
        -------
        None.

        """
        self.m1 = m1 = star_masses[0]
        self.m2 = m2 =  star_masses[1]
        self.M = m1+m2
        self.R1 = 1
        self.R2 = m1 / m2
        self.omega = 1

        
    def star_position(self,t):
        """
        Determines the position of the 2 stars at a specific time

        Returns
        -------
        X and Y positions of the 2 stars.

        """
        
        R1 = self.R1
        R2 = self.R2
        omega = self.omega
        
        X1 = R1*np.cos(omega*t)
        Y1 = R1*np.sin(omega*t)
        
        X2 = -R2*np.cos(omega*t)
        Y2 = -R2*np.sin(omega*t)
        
        return [X1,X2,Y1,Y2]
        
    
    def dSdt(self,t,S):
        """
        Contains the differential equation to be solved
        
        Parameters
        ----------
        t : FLOAT
            Time at which this iteration of the differential equation is calculated.
        S : LIST
            Contains the X and Y position and Velocity.

        Returns
        -------
        vx: FLOAT
            X velocity to be integrated and add to X position.
        vy: FLOAT
            Y velocity to be integrated and add to Y position
        dvx:FLOAT
            X acceleration to be integrated to find new X velocity
        dvy:FLOAT
            Y acceleration to be integrated to find new Y velocity

        """
        #Calls masses
        m1 = self.m1
        m2 = self.m2
        G = 1
        
        X1,X2,Y1,Y2 = self.star_position(t) #finds the X and Y position of stars 1&2
        
        x,y,vx,vy = S
        #Calculates distance between planets and stars
        r1 = np.sqrt((x-X1)**2 + (y-Y1)**2)
        r2 = np.sqrt((x-X2)**2 + (y-Y2)**2)
        #Calculates acceleration due to gravity of the 2 stars
        dvx = -G*m1*((1/r1**2)*(x-X1)/r1)-G*m2*((1/r2**2)*(x-X2)/r2)
        dvy = -G*m1*((1/r1**2)*(y-Y1)/r1)-G*m2*((1/r2**2)*(y-Y2)/r2)
        
        return [vx,
                vy,
                dvx,
                dvy]
    
    
    
    def Planet_path(self,S00,S01,S02,S03):
        """
        Solves the differential equation that governs the motion of the planet.
        G is normalised to 1
        
        Returns
        -------
        path : ODERESULT
            Contains the time, X&Y postion and velocity coordinates.

        """
        S0 = [S00,S01,S02,S03]
        #Solves differential equation with high accuracy and precision requirements. At times within the
        #array t
        #print(S0)
        

        path = solve_ivp(self.dSdt, [0,self.t.max()], S0, t_eval = self.t,
                         method='DOP853',rtol=1e-10,atol=1e-13)
        

        
        #print(path)
        
        return path
    
    def P_type(self,r):
        """
        Sets initial conditions for P-type orbits

        Parameters
        ----------
        r : float
            radius of orbit.

        Returns
        -------
        orbits : array
            Orbit Path.

        """
        #No. of sets of intail conditions in array
        self.t = np.linspace(0,300*np.pi,5000)
        n = 100
        a = np.linspace(0,2*np.pi,n)
        #Finds the velocity for a circular orbit in a two body system
        v = np.sqrt(self.M/r)
        #Sets intial conditions
        S0 = np.array([r*np.cos(a),r*np.sin(a),v*-np.sin(a),v*np.cos(a)]).transpose()
        #print(S0)
        #Runs the differential equation solver for the different intial conditions
        with Pool() as pool:
            orbits = pool.starmap(self.Planet_path,S0)
        #print(orbits)
        print(r/20)
        #self.motion(orbits,r,n)
        return orbits
    
    def S_type(self,r):
        """
        Sets initial conditions for S-type orbits

        Parameters
        ----------
        r : float
            radius of orbit.

        Returns
        -------
        orbits : array
            Orbit Path.

        """
        self.t = np.linspace(0,2*np.pi,1000)
        #No. of sets of intail conditions in array
        n = 50
        a = np.linspace(0,2*np.pi,n)
        self.m = [self.m1,self.m2]
        for j in range(1,2):
            #Finds the velocity for a circular orbit in a two body system
            v = np.sqrt(self.m[j]/r)
            #Sets intial conditions
            X1,X2,Y1,Y2 = self.star_position(self.t)
            X = [X1[0],X2[0]]
            Y = [Y1[0],Y2[0]]
            S0 = np.array([r*np.cos(a)+X[j],r*np.sin(a)+Y[j],v*-np.sin(a),v*np.cos(a)]).transpose()
            #print(S0)
            #Runs the differential equation solver for the different intial conditions
            with Pool() as pool:
                orbits = pool.starmap(self.Planet_path,S0)
            #print(orbits)
            print(r/0.5)
        #     view_r = 3
        #     fig = plt.figure()
        #     plt.axis("square")
        #     plt.xlabel("X (Orbtial Radii)")
        #     plt.ylabel("Y (Orbtial Radii)")
        #     plt.title("S type orbits at r={:0.2f}R".format(r))
        #     plt.ylim(-view_r,view_r)
        #     plt.xlim(-view_r,view_r)
        #     plt.plot(X1,Y1)
        #     plt.plot(X2,Y2)
        #     plt.scatter(X,Y)
        #     for k in range(0,n):
        #         x = orbits[k].y[0]
        #         y = orbits[k].y[1]
        #         print(x)
        #         plt.plot(x,y)
            
        # plt.show()
        
        return orbits
    
    def Set_ICs(self):
        """
        Runs initial conditions code for P and S type planets
        Then takes final plots 

        Returns
        -------
        None.

        """
        masses = np.linspace(1,10,10)
        P_stab_array = np.array([])
        self.r_P_type = r_P_type = np.linspace(3,20,10)
        r_S_type = np.linspace(0.01,0.5,10)
        for i in range(len(masses)):
            self.m1 = m1 = masses[i]
            self.M = m1+self.m2
            
            print(i/len(masses))
            
            self.R2 = m1 / self.m2
            
            
            
            P_type_vec = np.vectorize(self.P_type)
            P_orbits = P_type_vec(r_P_type)
            
            S_type_vec = np.vectorize(self.S_type)
            S_orbits = S_type_vec(r_S_type)
            
            #print(P_orbits[1])
            P_stab = self.stability(P_orbits,r_P_type)
            
            try:
                P_stab_array = np.vstack([P_stab_array,P_stab])
            except:
                P_stab_array = P_stab
        P_stab_array = np.flip(P_stab_array,axis=0)
        fig, ax = plt.subplots()
        plt.title("Stability of P-type Orbits")
        plt.xlabel("Radius of Planetary Orbit (normalised units)")
        plt.ylabel("Mass of Star 1")
        im = ax.imshow(P_stab_array,norm='log')
        cbar = fig.colorbar(im)
        cbar.set_label("Instability")
        ax.tick_params(which='minor', length=5)
        ax.set_yticks(np.arange(20),np.flip(masses))
        
        S_stab_array = np.array([])
        for i in range(len(masses)):
            self.m1 = m1 = masses[i]
            self.M = m1+self.m2
            
            
            print(i/len(masses))
            
            self.R2 = m1 / self.m2
            
            
            S_type_vec = np.vectorize(self.S_type)
            S_orbits = S_type_vec(r_S_type)
            
            #print(S_orbits)
            S_stab = self.stability_S_type(S_orbits,r_S_type)
            #print(S_stab)
            try:
                S_stab_array = np.vstack([S_stab_array,S_stab])
            except:
                S_stab_array = S_stab
        S_stab_array = np.flip(S_stab_array,axis=0)
        print(S_stab_array)
        #print(r_S_type)
        #print(masses)
        fig, ax = plt.subplots()
        plt.title("Stability of S-type Orbits")
        plt.xlabel("Radius of Planetary Orbit (normalised units)")
        plt.ylabel("Mass of Star 1")
        im = ax.imshow(S_stab_array,norm='log')
        cbar = fig.colorbar(im)
        cbar.set_label("Instability")
        ax.set_xticks([0,2,4,6,8],np.around(r_S_type[0::2],2))
        ax.set_yticks(np.arange(10),np.around(np.flip(masses),2))
    
    def stability(self,orb,r_array):
        #calculates the orbital instability
        instability = np.array([])
        radii = np.array([])
        for i in range(len(orb)):
            del_r_array = np.array([])
            for j in range(len(orb[i])):
                orbit = orb[i][j]
                x = orbit.y[0]
                y = orbit.y[1]
                r = np.sqrt(x**2 + y**2)
                #print(r)
                del_r = abs(np.sum((r-r[0])/r[0])/len(r))
                del_r_array = np.append(del_r_array,del_r)
                # print(del_r)
            av_del_r = np.sum(del_r_array)/len(del_r_array)
            # print(del_r_array)
            instability = np.append(instability,av_del_r)
        # print(instability)
        
        plt.plot(r_array,instability)
        return instability
    
    
    
    def stability_S_type(self,orb,r_array):
        #calculates the orbital instability for S-type orbits
        instability = np.array([])
        radii = np.array([])
        omega = self.omega
        t = self.t
        for i in range(len(orb)):
            del_r_array = np.array([])
            for j in range(len(orb[i])):
                orbit = orb[i][j]
                x = orbit.y[0]
                y = orbit.y[1]
                #zeta = x*np.cos(omega*t) + y*np.sin(omega*t)
                #eta = -x*np.sin(omega*t) + y*np.cos(omega*t)
                r = np.sqrt(x**2 + y**2)
                #print(r)
                del_r = abs(np.sum((r-r[0])/r[0])/len(r))
                del_r_array = np.append(del_r_array,del_r)
                # print(del_r)
            av_del_r = np.sum(del_r_array)/len(del_r_array)
            # print(del_r_array)
            
            instability = np.append(instability,av_del_r)
        print(instability)
        #plt.plot(eta,zeta)
        # plt.plot(r_array,instability)
        plt.show()
        return instability

        
    def motion(self,orb,r,n):
        """
        Determines the position of the 2 stars and plots them

        Returns
        -------
        None.

        """
        v = np.sqrt(self.M/r)
        
        
        S0 = [r,0,0,v]
        
        # stationary frame
        R1 = self.R1
        R2 = self.R2
        omega = self.omega
        t = self.t
        path = self.Planet_path(S0[0],S0[1],S0[2],S0[3])
        x = path.y[0]
        y = path.y[1]
        view_r = 25
        
        X1,X2,Y1,Y2 = self.star_position(t)
        
        fig = plt.figure()
        plt.axis("square")
        plt.xlabel("X (Orbtial Radii)")
        plt.ylabel("Y (Orbtial Radii)")
        plt.title("Single P type orbit at r={:0.2f}R".format(r))
        
        plt.ylim(-view_r,view_r)
        plt.xlim(-view_r,view_r)
        
        plt.plot(X1,Y1)
        plt.plot(X2,Y2)
        plt.plot(x,y)

        k = np.sqrt(self.M/(r**3))
        
        og_path = [r*np.cos(k*t),r*np.sin(k*t)]
        
        plt.plot(og_path[0],og_path[1])
        plt.show()
        
        # rotating frame
        Z1 = 1
        E1 = 0
        Z2 = -self.R2 
        E2 = 0
        fig = plt.figure()
        ax = fig.add_subplot()
 
        # square plot
        ax.set_aspect('equal', adjustable='box')
        plt.xlabel("X (Orbtial Radii)")
        plt.ylabel("Y (Orbtial Radii)")
        plt.title("P type orbits at r={:0.2f}R in rotating frame".format(r))
        zeta = x*np.cos(omega*t) + y*np.sin(omega*t)
        eta = -x*np.sin(omega*t) + y*np.cos(omega*t)

        plt.scatter(E1,Z1)
        plt.scatter(E2,Z2)
        plt.plot(zeta,eta,color="green")
        plt.show()
        
        
        fig = plt.figure()
        plt.axis("square")
        plt.xlabel("X (Orbtial Radii)")
        plt.ylabel("Y (Orbtial Radii)")
        plt.title("P type orbits at r={:0.2f}R".format(r))
        plt.ylim(-view_r,view_r)
        plt.xlim(-view_r,view_r)
        for i in range(0,n):
            plt.plot(orb[i].y[0],orb[i].y[1])
        plt.plot(X1,Y1)
        plt.plot(X2,Y2)
        plt.show()
        





planets = Binary_Planet_System([2,1])

if __name__=="__main__":
    freeze_support()
    planets.Set_ICs()






























