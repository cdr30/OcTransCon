"""
Module containing classes for the application of a
Kalman filter and RTS smoother.

"""

import numpy as np

from numpy import mean, sqrt, square


class KalmanModel(object):
    """
    Class containing the data matrices and methods for implementation of a 
    Kalman filter. Forward model predictions of OHC obey heat conservation.
    Forward model predictions of oht and flx are estimated using climatological
    values.
        
    """
    
    def __init__(self, config, flx_ob, ohc_ob, oht_ob, flx_ob_err, ohc_ob_err):
        """ Initialize Kalman Model object """
        self.config = config
        self.flx_ob = flx_ob
        self.ohc_ob = ohc_ob
        self.oht_ob = oht_ob
        self.flx_ob_err = flx_ob_err
        self.ohc_ob_err = ohc_ob_err
        self.dt = self.get_dt()
        self.ohc_fwd_err = self.get_ohc_fwd_err()
        self.flx_fwd_err = self.get_fwd_err(self.flx_ob)
        self.oht_fwd_err = self.get_fwd_err(self.oht_ob)
        self.flx_clim = self.flx_ob.mean()
        self.oht_clim = self.oht_ob.mean()      
                
    def get_dt(self):
        """ Return number of seconds in month, dt """
        return self.config.getfloat('kfilter', 'dt')
        
    def get_nmax(self):
        """ Return length of observation array, nmax """
        return len(self.ohc_ob)
    
    def get_y(self): 
        """ Return observation matrix, y """
        return [np.matrix([self.ohc_ob[n], self.flx_ob[n]]).T for n in range(self.get_nmax())]

    def get_E(self):
        """ Return observation operator, E """
        return np.matrix([[1, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 0, 0, 0, 0]])

    def get_R(self):
        """ Return observation error matrix, R """
        return np.matrix([[self.ohc_ob_err**2, 0],[0, self.flx_ob_err**2]])

    def get_fwd_err(self, dat):
        """ Return error in forward model predictions """
        return dat.std()

    def get_Bq(self):
        """ Return control matrix inputs, Bq"""
        return np.matrix([0, 0, self.flx_clim, 0, 0, self.oht_clim, 0, 0]).T
    
    def get_A(self):
        """ Return state transition matrix, A """
        return np.matrix([[0, 1, 0.5*self.dt, self.dt, 0.5*self.dt,
                            0.5*self.dt, self.dt, 0.5*self.dt],
                          [1, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 1, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 1, 0]])        
    
    def get_ohc_fwd_err(self):
        """ Return error in heat conservation forward model. """
        return np.std(self.ohc_ob) * 1.e-3
    
    def get_Q(self):
        """ Return forward model error covariance matrix, Q """
        return np.matrix([[self.ohc_fwd_err**2, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, self.flx_fwd_err**2, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, self.oht_fwd_err**2, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0]])
       
    def get_Gamma(self):
        """ Return gamma (identity) matrix """
        return np.matrix(np.identity(8))
    
    def get_P0(self):
        """ Return initial prediction uncertainty covariance """
        return np.matrix([[self.ohc_fwd_err**2, 0, 0, 0, 0, 0, 0, 0], 
                          [0, self.ohc_fwd_err**2, 0, 0, 0, 0, 0, 0],
                          [0, 0, self.flx_fwd_err**2, 0, 0, 0, 0, 0],
                          [0, 0, 0, self.flx_fwd_err**2, 0, 0, 0, 0],
                          [0, 0, 0, 0, self.flx_fwd_err**2, 0, 0, 0],
                          [0, 0, 0, 0, 0, self.oht_fwd_err**2, 0, 0],
                          [0, 0, 0, 0, 0, 0, self.oht_fwd_err**2, 0],
                          [0, 0, 0, 0, 0, 0, 0, self.oht_fwd_err**2]])
    
    def get_x0(self):
        """ Return initial values for state estimate """
        return np.matrix([self.ohc_ob[1], self.ohc_ob[0],
                          self.flx_ob[1], self.flx_ob[0], self.flx_clim,
                          self.oht_ob[1], self.oht_ob[0], self.oht_clim]).T    
                    
    def kalman_smooth(self): 
        """ Apply forward Kalman filter and backward RTS smoother """
        
        # Initialize values
        A = self.get_A()
        Bq = self.get_Bq()
        Q = self.get_Q()
        Gamma = self.get_Gamma()
        y = self.get_y()
        E = self.get_E()
        R = self.get_R()
        nmax = self.get_nmax()
        x = self.get_x0()
        P = self.get_P0()
        x_fwd, P_fwd, K_fwd = [], [], []
        x_fwd_pred, P_fwd_pred = [], []
    
        ### Forward Kalman filter
        for n in range(nmax):
            
            # Predict
            x = A * x + Bq
            P = A*P*A.T + Gamma*Q*Gamma.T
            
            # Save pre-update variables
            x_fwd_pred.append(x)
            P_fwd_pred.append(P)
            
            # Update
            K = P*E.T * np.linalg.inv(E*P*E.T + R)
            x = x + K * (y[n] - E * x)
            P = P - K*E*P
            
            # Save post-update variables
            x_fwd.append(x)
            P_fwd.append(P)
            K_fwd.append(K)  
            
        ### Backwards RTS smoother.
        x_bwd = list(x_fwd)
        P_bwd = list(P_fwd)
        
        for t1 in np.arange(nmax-1)[::-1] + 1:
            t =  t1 - 1
            L = P_fwd[t] * A.T * np.linalg.inv(P_fwd_pred[t1])
            x_bwd[t] = x_fwd[t] + L*(x_bwd[t1] - x_fwd_pred[t1])
            P_bwd[t] = P_fwd[t] + L *(P_bwd[t1] - P_fwd_pred[t1])* L.T
            
        return x_fwd, P_fwd, x_bwd, P_bwd


class KalmanModelRatePersist(KalmanModel):
    """ 
    Overrides KalmanModel so that forward model predictions of oht 
    and flx are estimated using rate-of-change persistence.
    
    """
    def get_fwd_err(self, dat):
        """ Return error in forward model predictions """
        tru = dat[2:]
        rate_per = 2 * dat[1:-1] - dat[0:-2]
        diff = tru - rate_per
        return sqrt(mean(square(diff)))
        
    def get_Bq(self):
        """ Return control matrix inputs, Bq"""
        return np.matrix([0, 0, 0, 0, 0, 0, 0, 0]).T
    
    def get_A(self):
        """ Return state transition matrix, A """
        return np.matrix([[0, 1, 0.5*self.dt, self.dt, 0.5*self.dt,
                            0.5*self.dt, self.dt, 0.5*self.dt],
                          [1, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 2,-1, 0, 0, 0, 0],
                          [0, 0, 1, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 2,-1, 0],
                          [0, 0, 0, 0, 0, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 1, 0]])        


class KalmanModelPersist(KalmanModel):
    """ 
    Overrides KalmanModel so that forward model predictions of oht 
    and flx are estimated using anomaly persistence.
    
    """
    def get_fwd_err(self, dat):
        """ Return error in forward model predictions """
        tru = dat[1:]
        anom_per = dat[:-1]
        diff = tru - anom_per
        return sqrt(mean(square(diff)))

    def get_Bq(self):
        """ Return control matrix inputs, Bq"""
        return np.matrix([0, 0, 0, 0, 0, 0, 0, 0]).T
    
    def get_A(self):
        """ Return state transition matrix, A """
        return np.matrix([[0, 1, 0.5*self.dt, self.dt, 0.5*self.dt,
                            0.5*self.dt, self.dt, 0.5*self.dt],
                          [1, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 1, 0, 0, 0, 0, 0],
                          [0, 0, 1, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 1, 0, 0],
                          [0, 0, 0, 0, 0, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 1, 0]])        

        

def apply_ksmooth(config, flx_ob, ohc_ob, oht_ob, flx_ob_err, ohc_ob_err):
    """ Generate model, apply Kalman filter, and return output as a dictionary. """
    
    fwd_model = config.get('kfilter', 'fwd_model').upper().strip()
    
    if fwd_model == 'CLIMATOLOGY':
        kmodel = KalmanModel(config, flx_ob, ohc_ob, oht_ob, flx_ob_err, ohc_ob_err)
    elif fwd_model == 'PERSISTENCE':
        kmodel = KalmanModelPersist(config, flx_ob, ohc_ob, oht_ob, flx_ob_err, ohc_ob_err)
    elif fwd_model == 'RATE_PERSISTENCE':
        kmodel = KalmanModelRatePersist(config, flx_ob, ohc_ob, oht_ob, flx_ob_err, ohc_ob_err)
    else:
        raise ValueError('fwd_model = %s: method not recognized.' % fwd_model)
            
    x_fwd, P_fwd, x_bwd, P_bwd = kmodel.kalman_smooth()
   
    ### Extract results from state matrices
    kout = {}
    kout['ohc_ksmooth'] = np.array([m[0,0] for m in x_bwd])
    kout['flx_ksmooth'] = np.array([m[3,0] for m in x_bwd])
    kout['oht_ksmooth'] = np.array([m[6,0] for m in x_bwd])
    kout['ohc_kfwd'] = np.array([m[0,0] for m in x_fwd])
    kout['flx_kfwd'] = np.array([m[3,0] for m in x_fwd])
    kout['oht_kfwd'] = np.array([m[6,0] for m in x_fwd])

    ### Extract uncertainties from error covariance matrices
    kout['ohc_ksmooth_err'] = np.array([np.sqrt(m[0,0]) for m in P_bwd])
    kout['flx_ksmooth_err'] = np.array([np.sqrt(m[3,3]) for m in P_bwd])
    kout['oht_ksmooth_err'] = np.array([np.sqrt(m[6,6]) for m in P_bwd])
    kout['ohc_kfwd_err'] = np.array([np.sqrt(m[0,0]) for m in P_fwd])
    kout['flx_kfwd_err'] = np.array([np.sqrt(m[3,3]) for m in P_fwd])
    kout['oht_kfwd_err'] = np.array([np.sqrt(m[6,6]) for m in P_fwd])

    return kout
 


