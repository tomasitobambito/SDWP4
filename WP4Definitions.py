#Imports
import math as m

#Seperate Failure Modes
def Tension(D_h,w,t,K_t,S_tu,F_app):
    A_t = (w-D_h)*t
    P_tu = K_t*S_tu*A_t
    return (P_tu,P_tu/F_app)

def Shear(D_h,w,t,K_s,S_su,F_app):
    a = (w/2-D_h/2)
    L_sp = a
    A_s = 2*L_sp*t
    P_su = K_s*A_s*S_su
    return (P_su,P_su/F_app)

def Bearing(D_h,t,K_br,S_bru,F_app):
    A_br = D_h*t
    P_bru = K_br*S_bru*A_br
    return (P_bru,P_bru/F_app)

#Axial Direction Unified Method
def Axial(w, D_h, a, t, K_br1, K_br2, S_bry, F_app):
    if w/2/D_h < 1.5:
        P_bry = 1.304*K_br1*D_h*t*S_bry*a/D_h
        return (P_bry,P_bry/F_app)
    if w/2/D_h >= 1.5:
        P_bry = 1.304*S_bry*K_br2*D_h*t
        return (P_bry,P_bry/F_app)
    
def Axial1(w, D_h, a, t,K_br3, S_bry, F_app):
    if w/2/D_h < 1.5:
        P_bry = 1.304*K_br3*D_h*t*S_bry*a/D_h
        return (P_bry,P_bry/F_app)
    if w/2/D_h >= 1.5:
        P_bry = 1.304*S_bry*K_br3*D_h*t
        return (P_bry,P_bry/F_app)

#Transverse Direction Unified Method
def Transverse(D_h,t,K_try,S_ty,F_app):
    P_tru = K_try*S_ty*D_h*t
    return (P_tru,P_tru/F_app)

#Weight Determination
def Weight(w,D_h,t,density):
    return ((w/2)**2*m.pi/2+w/2*w-(D_h/2)**2*m.pi)*t*density

#Iterative Function
def Iterate(D_h,t,w,F_app,K_t,S_tu,S_su,S_bru,density,error = 1,SF = 2,correction =0.1,track = 1):
    #Geometry
    a = lambda w,D_h : (w/2-D_h/2)
    #Stress Factors
    K_s = lambda x :1.2*x-0.6 #Independent
    K_br = lambda x : 2*x-1 #Independent
    K_br1 = lambda x : -0.75*x + 2.35 #Independent
    K_br2 = lambda x : 0.8*x + 0.1 #Independent
    K_br3 = lambda x : 1.6*x- 0.8 #Independent
    K_try = lambda x : -0.4284*x*x + 1.5092*x #Independent
    #Bearing Strenght
    S_bry = lambda a : K_br(w/2/D_h)*a/D_h*S_tu
    #Distances
    h_2 = h_3 = lambda w,D_h :0.5*(w-D_h)
    h_1 = h_4 = lambda h_2 :h_2 + 0.5*D_h*(2**0.5/2)
    h_av = lambda h_1,h_2,h_3,h_4 : 6/(3/h_1+1/h_2+1/h_3+1/h_4)
    #Optimal width finder
    while error >= 0.01 or error <= 0:
        SF_t = Tension(D_h, w, t, K_t, S_tu, F_app)[1]
        SF_s = Shear(D_h, w, t, K_s(w/2/D_h), S_su, F_app)[1]
        SF_b = Bearing(D_h, t, K_br(w/2/D_h), S_bru, F_app)[1]
        SF_tr = Transverse(D_h, t, K_try(h_av(h_1(h_2(w,D_h)),h_2(w,D_h),h_3(w,D_h),h_4(h_2(w,D_h)))/D_h), S_tu, F_app)[1]
        SF_a = Axial(w, D_h, a(w,D_h), t, K_br1(w/2/D_h), K_br2(w/2/D_h), S_bry(a(w,D_h)), F_app)[1]
        if D_h/t > 5:
            SF_a = Axial1(w, D_h, a(w,D_h), t, K_br3(w/2/D_h), S_bry(a(w,D_h)), F_app)[1]
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
    return Weight(w, D_h, t, density),w