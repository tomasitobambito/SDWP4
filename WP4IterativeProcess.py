from WP4Definitions import Iterate

#Known Values
D_h = 0.035
F_app = 326.7/4
S_t = 414*10**6
S_s = 290*10**6
K_t = 0.98 #Dependent
K_b = lambda x : 2*x-1 #Independent
a = lambda w,D_h : (w/2-D_h/2)
S_br1 = 662*10**6
density = 2800

#Obtained Initial Values
t = 0.0065
w = 0.037

#Iterate Setup
weight1 = Iterate(D_h, t, w, F_app, K_t, S_t, S_s, S_br1, density)[0]
error = track = 1
correction = 0.0001
t -= correction

while track:
     weight2,w1 = Iterate(D_h, t, w, F_app, K_t, S_t, S_s, S_br1, density)
     if weight1 > weight2 and t >= 0.001:
        t -= 1*correction
        t = round(t,10)
        weight1 = weight2
     elif correction < 10**(-15) or t == 0.001:
         track = 0
     else:
        t += 1*correction
        correction *= 0.1
print('*'*24)
print('Thickness: {}[mm]\nWidth: {}[mm]\nWeight: {}[g]'.format(t*10**3,round(w1*10**3,2),round(weight1*10**3,2)))
print('*'*24)