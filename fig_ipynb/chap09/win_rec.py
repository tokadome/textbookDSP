import numpy as np

def win_rec(n, M=16, WB=0, fillzero=False):
  '''Rectangular window.
       n: discrete time. ndarray, int
       M: window size; dafault = 16
       WB: begin window; default = 0
       fillzero: value out of window; default value = nan

       return: rectangular window. ndarray'''
  
  l_n = len(n)
  w = np.ones(l_n)
  if (l_n <= 0):
    print("Argument n seems to be null.\n")
    return(None)

  W_lb = WB
  W_ub = M-1+WB
  if (fillzero == False): s = np.nan
  else: s = 0
  for i in range(l_n):
    if ((n[i] < W_lb) or (n[i] > W_ub)): w[i] = s
  return(w)

def extract_by_win_rec(n, f, M=16, WB=0, fillzero=False):
  '''Signal extraction by Rectangular window.
       n: discrete time. ndarray, int
       f: function of n. ndarray, float
       M: window size; dafault = 16
       WB: begin window; default = 0
       fillzero: value out of window; default value = nan

       return: windowed f. ndarray'''
  
  l_n = len(n)
  l_f = len(f)
  f_tmp = np.array(f)
  if ((l_n <= 0) or (l_f <= 0)):
    print("Arguments n or f seem to be null.\n")
    return(None)
  elif (l_n != l_f):
    print("Length of n and f must be equal.")
    return(None)

  W_lb = WB
  W_ub = M-1+WB
  if (fillzero == False): s = np.nan
  else: s = 0
  for i in range(l_n):
    if ((n[i] < W_lb) or (n[i] > W_ub)): f_tmp[i] = s
  return(f_tmp)


def win_rec_spectral(M=16, omg_shift=0, l_range=-2*np.pi, r_range=2*np.pi, step=0.001):
  '''Rectangular window.
       M: window size; dafault = 16
       omg_shift: shift spectral to frequency by omg_shift; dafault = 0
       l_range= left of frequency range; default = -2*pi
       r_range= right of frequency range; default = 2*pi
       step: default=0.001

       return: Spectral of recutanglar window. ndarray'''

  omg = np.arange(l_range, r_range, step) - omg_shift
  W = (np.sin((M/2)*omg)/np.sin((1/2)*omg))*np.exp(omg*(-(M-1)/2)*1j)
  
  return(W)
    
