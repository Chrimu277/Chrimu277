##########################################################################################
## UE Wahrscheinlichkeitstheorie, Datenanalyse, Statistik
## MCMC fuer Aufgabe 3c - Blatt 6 - SS 2019
##########################################################################################

N       = len(x)          # Laenge des Arrays aus Aufgabe b)  
a       = [0,0]           # Startwerte
da      = 0.34            # Schrittweite
N_mess  = 1000            # Messungen
N_skip  = 20              # ausgelassene Werte
N_loop  = N_mess * N_skip # Schritte ingesamt
sigma   = 0.75
alpha   = N/(2*sigma**2)

def calc_phi(a,b,avg_x,avg_y,avg_x2,avg_xy,avg_y2):
    phi = a**2 + b**2 * avg_x2 + 2 * a * b * avg_x - 2 * a * avg_y - 2* b * avg_xy + avg_y2 # avg_x, avg_y, ... aus Aufgabe b)
    return(phi) 

phi0      = calc_phi(a[0], a[1], avg_x,avg_y,avg_x2,avg_xy,avg_y2)
list0     = np.zeros([N_mess,2])
list_corr = np.zeros([N_loop,2]) # korrelierte Messdaten

for i in range(0,N_loop):
    a_t   = a + (np.random.random(2) - 0.5) * da
    phi_t = calc_phi(a_t[0], a_t[1], avg_x,avg_y,avg_x2,avg_xy,avg_y2)
    
    if np.log(np.random.random(1)) < (alpha * (phi0 - phi_t)):
        a       = a_t
        phi0    = phi_t
    list_corr[i,:] = a
    list0[np.int(np.ceil(i/N_skip)-1),:] = a

min_len = np.int(np.ceil(0.1 * N_mess))-1
max_len = len(list0)-1

avg_a    = np.array([0,0],dtype=float)
avg_a[0] = np.mean(list0.T[0][np.linspace(min_len,max_len,1, dtype=int)])
avg_a[1] = np.mean(list0.T[1][np.linspace(min_len,max_len,1, dtype=int)])

C    = np.array([0,0])
C[0] = np.cov(list0.T[0][np.linspace(min_len,max_len,1, dtype=int)])
C[1] = np.cov(list0.T[1][np.linspace(min_len,max_len,1, dtype=int)])
std  = np.sqrt(np.diag(C))
