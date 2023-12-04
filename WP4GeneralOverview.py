from WP4Definitions import Transverse,Bearing,Axial,Axial1,Tension,Shear,Weight
#Variables
D_h = 0.028
t = 0.001
w = D_h + 0.002
F_app = 326.7/4

#Stress Factors
K_t = 0.98 #Dependent
K_s = lambda x :1.2*x-0.6 #Independent
K_b = lambda x : 2*x-1 #Independent
K_a1 = lambda x : -0.75*x + 2.35 #Independent
K_a2 = lambda x : 0.8*x + 0.1 #Independent
K_a3 = lambda x : 1.6*x- 0.8 #Independent
K_tr = lambda x : -0.4284*x*x + 1.5092*x #Independent

#Properties Material Alu 2014-T651
S_y = 276*10**6
S_s = 207*10**6
S_br1 = 386*10**6
density = 2700

#Properties Material Steel 4130
# S_y = 360*10**6
# S_s = 337*10**6
# S_br1 = 540*10**6
# density = 7900

S_br2 = lambda a : K_b(w/2/D_h)*a/D_h*S_y

#Distances
a = lambda w,D_h : (w/2-D_h/2)
h_2 = h_3 = lambda w,D_h :0.5*(w-D_h)
h_1 = h_4 = lambda h_2 :h_2 + 0.5*D_h*(2**0.5/2)
h_av = lambda h_1,h_2,h_3,h_4 : 6/(3/h_1+1/h_2+1/h_3+1/h_4)

#Iteration Variables
error = 1
SF = 2
correction = 0.01
track = 1

#Optimal width finder
while error >= 0.01 or error <= 0:
    P_t,SF_t = Tension(D_h, w, t, K_t, S_y, F_app)
    P_s,SF_s = Shear(D_h, w, t, K_s(w/2/D_h), S_s, F_app)
    P_b,SF_b = Bearing(D_h, t, K_b(w/2/D_h), S_br1, F_app)
    P_a,SF_a = Axial(w, D_h, a(w,D_h), t, K_a1(w/2/D_h), K_a2(w/2/D_h), S_br2(a(w,D_h)), F_app)
    P_tr,SF_tr = Transverse(D_h, t, K_tr(h_av(h_1(h_2(w,D_h)),h_2(w,D_h),h_3(w,D_h),h_4(h_2(w,D_h)))/D_h), S_s, F_app)
    if D_h/t > 5:
        P_a,SF_a = Axial1(w, D_h, a(w,D_h), t, K_a3(w/2/D_h), S_br2(a(w,D_h)), F_app)
    error = min(SF_t-SF,SF_s-SF,SF_b-SF,SF_tr-SF,SF_a-SF)
    if error <=  0.01 and track:
     track = 1
     w += 1*correction
    if error > 0.01 and not track:
     track = 0
     w -= 1*correction
    if error > 0.01 and track:
        w -= 1*correction
        correction *= 0.1
    if error <= 0.01 and not track:
        w += 1*correction
        correction *= 0.1


print('*'*24)
print('Thickness: {} [mm]\nWidth: {} [mm]\nDiameter: {} [mm]'.format(t*10**3,round(w*10**3,3),D_h*10**3))
print("Weight: {}[g]".format(round(Weight(w, D_h, t, density)*10**3,3)))
print('*'*24)
print('SF_t:{} | SF_s:{}\nSF_b:{} | SF_a:{} \nSF_tr:{} | S_bry:{}*E6'.format(round(SF_t,0),round(SF_s,0),round(SF_b,0),round(SF_a,0),round(SF_tr,0),round(S_br2(a(w,D_h))/10**6)))
print('*'*24)
print('M.S. Oblique Loads : {}'.format(round(1/((F_app/P_a)**1.6+(F_app/P_tr)**1.6)**0.625-1)))
print('*'*24)