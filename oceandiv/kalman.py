"""
Module containing routines and model definitions for the application 
of a Kalman smoothing algorithm.

"""

import numpy as np


def create_kalman_model(config, flx_ob, ohc_ob, flx_ob_err, ohc_ob_err):
    """
    Return the matrices and initial state vectors for implementation of a 
    Kalman filter. Forward model predictions of OHC obey heat conservation.
    Initial estimates of oht and flx assume persistence of the previous anomaly.
    
    """
    
    ### Read configuration options.
    dt = config.getfloat('kfilter', 'dt')
    oht_err_scale = config.getfloat('kfilter', 'oht_error_scaling')
    
    ### Estimate forward model errors
    ohc_fwd_err = np.std(ohc_ob) * 1.e-3             # Small, as conservation is accurate approximation.
    flx_fwd_err = np.std(flx_ob[1:] - flx_ob[:-1])   # Error in flux persistence assumption.
    oht_fwd_err = flx_fwd_err * ohts_err_scale       # Scaled estimate of error in ohtsport persistence

    ### Create observation matrices
    nmax = len(ohc_ob)
    y = [np.matrix([ohc_ob[n], flx_ob[n]]).T for n in range(nmax)]
    E = np.matrix([[1, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0, 0, 0]])
    R = np.matrix([[ohc_ob_err**2, 0],[0, flx_ob_err**2]])

    ### Create forward model matrices
    A = np.matrix([[0, 1, 0.5*dt, dt, 0.5*dt, 0.5*dt, dt, 0.5*dt],
                   [1, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0, 0, 1, 0]])

    Q = np.matrix([[ohc_fwd_err**2, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, flx_fwd_err**2, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, oht_fwd_err**2, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0]])

    Gamma = np.matrix(np.identity(8))

    ### Initialize state vectors
    P0 = np.matrix([[ohc_ob_err**2, 0, 0, 0, 0, 0, 0, 0], 
                   [0, ohc_ob_err**2, 0, 0, 0, 0, 0, 0],
                   [0, 0, flx_ob_err**2, 0, 0, 0, 0, 0],
                   [0, 0, 0, flx_ob_err**2, 0, 0, 0, 0],
                   [0, 0, 0, 0, flx_ob_err**2, 0, 0, 0],
                   [0, 0, 0, 0, 0, flx_ob_err**2, 0, 0],
                   [0, 0, 0, 0, 0, 0, flx_ob_err**2, 0],
                   [0, 0, 0, 0, 0, 0, 0, flx_ob_err**2]])
    
    x0 = np.matrix([ohc_ob[0],ohc_ob[0],
                    flx_ob[0], flx_ob[0], flx_ob[0],
                    -flx_ob[0], -flx_ob[0], -flx_ob[0]]).T    

    return A, Q, Gamma, y, E, R, P0, x0


def kalman_smooth(A, Q, Gamma, y, E, R, P0, x0):
    """
    Apply forward kalman filter and backward RTS smoother using
    specified linear model and data matrices.
    
    """
    
    ### Initialize 
    nmax = len(y)
    x = x0
    P = P0
    x_fwd, P_fwd, K_fwd = [], [], []
    x_fwd_pred, P_fwd_pred = [], []

    ### Forward Kalman filter
    for n in range(nmax):
        
        # Predict
        x = A * x
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


def apply_ksmooth(config, flx_ob, ohc_ob, flx_ob_err, ohc_ob_err):
    """ Apply Kalman filter to the specified data. """
    
    A, Q, Gamma, y, E, R, P0, x0 = create_kalman_model(
        config, flx_ob, ohc_ob, flx_ob_err, ohc_ob_err)
    x_fwd, P_fwd, x_bwd, P_bwd = kalman_smooth(A, Q, Gamma, y, E, R, P0, x0)
   
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
 


