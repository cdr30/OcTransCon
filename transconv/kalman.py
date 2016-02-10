"""
Module containing code and matrix definitions to apply forward Kalman filter
and backward RTS smoother.

"""

import numpy as np


def create_kalman_model(config, flx_ob, zsum_ob, flx_ob_err, zsum_ob_err):
    """
    Return the matrices and initial state vectors for application of 
    local conservation with a Kalman smoother. Initial predictions of surfaces 
    fluxes and transport convergence are based on persistence.
    
    """
    dt = config.getfloat('kfilter', 'dt')
    trans_err_scale = config.getfloat('kfilter', 'transport_error_scaling')
    
    ### Estimate forward model errors
    zsum_fwd_err = np.std(zsum_ob) * 1e-3           # Small, as conservation should be accurate.
    flx_fwd_err = np.std(flx_ob[1:] - flx_ob[:-1])  # Estimate of error in flux persistence 
    tran_fwd_err = flx_fwd_err * trans_err_scale    # Scaled estimate of error in transport persistence

    ### Create observation matrices
    nmax = len(zsum_ob)
    y = [np.matrix([zsum_ob[n], flx_ob[n]]).T for n in range(nmax)]
    E = np.matrix([[1, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0, 0, 0]])
    R = np.matrix([[zsum_ob_err**2, 0],[0, flx_ob_err**2]])

    ### Create forward model matrices
    A = np.matrix([[0, 1, 0.5*dt, dt, 0.5*dt, 0.5*dt, dt, 0.5*dt],
                   [1, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0, 0, 1, 0]])
    Q = np.matrix([[zsum_fwd_err**2, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, flx_fwd_err**2, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, tran_fwd_err**2, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0]])
    Gamma = np.matrix(np.identity(8))

    # Initialize state vectors
    P0 = np.matrix([[zsum_ob_err**2, 0, 0, 0, 0, 0, 0, 0], # Large and uncorrelated
                   [0, zsum_ob_err**2, 0, 0, 0, 0, 0, 0],
                   [0, 0, flx_ob_err**2, 0, 0, 0, 0, 0],
                   [0, 0, 0, flx_ob_err**2, 0, 0, 0, 0],
                   [0, 0, 0, 0, flx_ob_err**2, 0, 0, 0],
                   [0, 0, 0, 0, 0, flx_ob_err**2, 0, 0],
                   [0, 0, 0, 0, 0, 0, flx_ob_err**2, 0],
                   [0, 0, 0, 0, 0, 0, 0, flx_ob_err**2]])
    
    x0 = np.matrix([zsum_ob[0],zsum_ob[0],
                    flx_ob[0], flx_ob[0], flx_ob[0],
                    -flx_ob[0], -flx_ob[0], -flx_ob[0]]).T    # Zeros...

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
        
        #update
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


def apply_ksmooth(config, flx_ob, zsum_ob, flx_ob_err, zsum_ob_err):
    """
    Generate model matrices and apply forward Kalman filter and
    backward RTS smoother to observed data.
        
    """
    A, Q, Gamma, y, E, R, P0, x0 = create_kalman_model(
        config, flx_ob, zsum_ob, flx_ob_err, zsum_ob_err)
    x_fwd, P_fwd, x_bwd, P_bwd = kalman_smooth(A, Q, Gamma, y, E, R, P0, x0)
   
    ### Extract results from state matrices
    kout = {}
    kout['zsum_kb'] = np.array([m[0,0] for m in x_bwd])
    kout['flx_kb'] = np.array([m[3,0] for m in x_bwd])
    kout['tran_kb'] = np.array([m[6,0] for m in x_bwd])
    kout['zsum_kf'] = np.array([m[0,0] for m in x_fwd])
    kout['flx_kf'] = np.array([m[3,0] for m in x_fwd])
    kout['tran_kf'] = np.array([m[6,0] for m in x_fwd])

    ### Extract uncertainties from error covariance matrices
    kout['zsum_kb_err'] = np.array([np.sqrt(m[0,0]) for m in P_bwd])
    kout['flx_kb_err'] = np.array([np.sqrt(m[3,3]) for m in P_bwd])
    kout['tran_kb_err'] = np.array([np.sqrt(m[6,6]) for m in P_bwd])
    kout['zsum_kf_err'] = np.array([np.sqrt(m[0,0]) for m in P_fwd])
    kout['flx_kf_err'] = np.array([np.sqrt(m[3,3]) for m in P_fwd])
    kout['tran_kf_err'] = np.array([np.sqrt(m[6,6]) for m in P_fwd])

    return kout
 


