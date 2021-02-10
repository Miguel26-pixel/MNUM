import matplotlib.pyplot as plt
import numpy as np
from math import e


#Substance elimination rate taken from project specification
#Ke converted from h^(-1) to min^(-1)
Ke = 0.0693 * (1 / 60)


#Ka - TBD
Ka = 0


#Distibution function slope
#Calculated se
m_d = (2 / 9) / (3 * 60)


#Time in minutes from ingestion of substance until it is completly absorbed
tmax = 3 * 60


#Start time in minutes
t_0 = 0


#Mass of substance in central and plasmatic compartiments at t_0
mi_0 = 0
mp_0 = 0


#Period of ingestion of substance in mins
period = 24 * 60 


#Plasmatic volume in ml
Vap = 3650


#End time in hours
t_f = 500


#Time interval for integration
delta_t = 1


def bissect_method_recur(a, b, func, prec):
    if (func(a) * func(b) > 0):
        return "Invalid arguments\n"
    elif (b - a > prec):
        m = (a + b) / 2
        if (func(a) * func(m) < 0):
            return bissect_method_recur(a, m, func, prec)
        else:
            return bissect_method_recur(m, b, func, prec)
    else:
        return (a + b) / 2
        
        
def Ka_eq(K):
    return K * e ** (- K * tmax) - Ke * e ** (- Ke * tmax)
    

#Find Ka in min^(-1)   
Ka = bissect_method_recur(0, 0.01, Ka_eq, 0.00001)


#Distribution funtion
#t in minutes
#return in mg/min
def D(t):
    pos = t % period
    return (m_d * pos if pos <= (3 * 60) else 0) if t < (4320 * 60) else 0


def diff_mi(t, mi, mp):
    return D(t) - Ka * mi
 
    
def diff_mp(t, mi, mp):
    return Ka * mi - Ke * mp
            
            
def hours_to_minutes(m):
    return m * 60
    

#singular integral (quadrature)
def quadrature(f, t0, t1, n):
    h = (t1-t0)/n
    return h/2 * (f(t0) + f(t1) + 2*sum(f(t0 + i*h) for i in range(1, n)))


#rk2
def rk2(f1, f2, t, mi, mp, dt):
    mi_1 = dt * f1(t, mi, mp)
    mp_1 = dt * f2(t, mi, mp)
    
    mi_2 = f1(t + dt/2, mi + mi_1/2, mp + mp_1/2)
    mp_2 = f2(t + dt/2, mi + mi_1/2, mp + mp_1/2)
    
    return mi_2*dt, mp_2*dt
    

#solve main odes
def solve_odeqs(F1, F2, t0, mi0, mp0, dt, tf):
    t = [t0]
    mi = [mi0]
    mp = [mp0]
    
    t1 = t0
    mi1 = mi0
    mp1 = mp0
    
    n = int((tf - t0) / dt)
    
    for _ in range(n):
        temp = rk2(F1, F2, t1, mi1, mp1, dt)
        
        mi1 += temp[0]
        mp1 += temp[1]
        t1 += dt
        
        t.append(t1)
        mi.append(mi1)
        mp.append(mp1)
        
    return t, mi, mp


t_points_1, mi_points_1, mp_points_1 = solve_odeqs(diff_mi, diff_mp, t_0, mi_0, mp_0, delta_t, hours_to_minutes(t_f))

#Total mass in system at time
mt_points = [x + y for x, y in zip(mi_points_1, mp_points_1)]

#Plasmatic concentration
cm_points = [ m / Vap for m in mp_points_1]


#Plot different points
#plt.scatter(t_points_1, mi_points_1, label = "mi(t)")
#plt.scatter(t_points_1, mp_points_1, label = "mp(t)")
#plt.scatter(t_points_1, mt_points, label = "mt(t)")
plt.scatter(t_points_1, cm_points, label = "cm(t)")

#Graphic legends
plt.legend()
plt.xlabel("time (min)")
plt.ylabel("cm (mg/ml)")
plt.title("Concentration of lercanidipine in the plasmatic compartment as a function of time")


print("Ka: ", Ka)
print("Ke: ", Ke)


#Calculate for dt/2 and dt/4
t_points_2, mi_points_2, mp_points_2 = solve_odeqs(diff_mi, diff_mp, t_0, mi_0, mp_0, delta_t / 2, hours_to_minutes(t_f))
t_points_3, mi_points_3, mp_points_3 = solve_odeqs(diff_mi, diff_mp, t_0, mi_0, mp_0, delta_t / 4, hours_to_minutes(t_f))

#Calculate QC's
mp_QC = (mp_points_2[-1] - mp_points_1[-1]) / (mp_points_3[-1] - mp_points_2[-1])
mi_QC = (mi_points_2[-1] - mi_points_1[-1]) / (mi_points_3[-1] - mi_points_2[-1])

print("mp_QC: ", mp_QC)
print("mi_QC: ", mi_QC)


#Calculate error
mp_error = (mp_points_3[-1] - mp_points_2[-1]) / (2 ** 2 - 1)
mi_error = (mi_points_3[-1] - mi_points_2[-1]) / (2 ** 2 - 1)

print("mp_error: ", mp_error)
print("mi_error: ", mi_error)

#Show plot
plt.show()


'''
#Plot D(t)
vectD = np.vectorize(D)

d = np.arange(0.0, 24 * 60, 1)

T = vectD(d)

plt.plot (d, T, 'bo', d, T, 'k')
plt.show()
'''
