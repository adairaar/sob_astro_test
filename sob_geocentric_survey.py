import numpy as np, sys

mu_E=3.986004418e5; R_E=6371.0; d2r=np.pi/180; r2d=180/np.pi
phi=31.77*d2r; lam=35.24*d2r

def kepler_E(M,e):
    E=M.copy()
    for _ in range(50):
        dE=(M-E+e*np.sin(E))/(1-e*np.cos(E)); E+=dE
        if np.max(np.abs(dE))<1e-11: break
    return E

def run_orbit(a,e,i,om,Om,M0,t):
    n=np.sqrt(mu_E/a**3); M=(M0+n*t)%(2*np.pi); E=kepler_E(M,e)
    nu=2*np.arctan2(np.sqrt(1+e)*np.sin(E/2),np.sqrt(1-e)*np.cos(E/2))
    r=a*(1-e*np.cos(E)); xp=r*np.cos(nu); yp=r*np.sin(nu)
    cO,sO=np.cos(Om),np.sin(Om); ci,si=np.cos(i),np.sin(i)
    co,so=np.cos(om),np.sin(om)
    rx=cO*co-sO*so*ci; ry=-(cO*so+sO*co*ci)
    sx=sO*co+cO*so*ci; sy=-(sO*so-cO*co*ci)
    tx=so*si; ty=co*si
    x=rx*xp+ry*yp; y=sx*xp+sy*yp; z=tx*xp+ty*yp
    GMST=2.395+7.2921150e-5*t; LST=GMST+lam; Robs=R_E+0.76
    ox=Robs*np.cos(phi)*np.cos(LST); oy=Robs*np.cos(phi)*np.sin(LST)
    oz=Robs*np.sin(phi)*np.ones_like(t)
    dx=x-ox; dy=y-oy; dz=z-oz
    rho=np.sqrt(dx**2+dy**2+dz**2)
    sL,cL=np.sin(LST),np.cos(LST); sp,cp=np.sin(phi),np.cos(phi)
    S=sp*cL*dx+sp*sL*dy-cp*dz
    Ec=-sL*dx+cL*dy
    Z=cp*cL*dx+cp*sL*dy+sp*dz
    alt=np.arcsin(np.clip(Z/rho,-1,1))*r2d
    az=(np.arctan2(Ec,-S)*r2d)%360
    return alt,az,rho

# 5-minute time grid, 12-hour night window
DT=300; t=np.arange(0,12*3600+DT,DT); DH=DT/3600.

# Precomputed sun altitude (simplified model for 5 BCE Jun 8)
# Sun at RA~5.5h, Dec~23°; local midnight ~6h into window
# Build sun alt array analytically
JD0=1719755.898; JD2K=2451545.0; TJ=(JD0-JD2K)/36525
L0=(280.46646+36000.76983*TJ)%360; M0s=(357.52911+35999.05029*TJ)%360
t_day=t/86400
L=(L0+360.98564724*t_day)%360; Ms=(M0s+360.98564724*t_day)%360
C=(1.9146-0.004817*TJ)*np.sin(Ms*d2r)+0.019993*np.sin(2*Ms*d2r)
lon_s=(L+C)%360; eps=23.439-0.013*TJ
dec_s=np.arcsin(np.sin(eps*d2r)*np.sin(lon_s*d2r))
ra_s=np.arctan2(np.cos(eps*d2r)*np.sin(lon_s*d2r),np.cos(lon_s*d2r))
GMST_s=2.395+7.2921150e-5*t; HA_s=GMST_s+lam-ra_s
sun_alt=np.arcsin(np.sin(phi)*np.sin(dec_s)+np.cos(phi)*np.cos(dec_s)*np.cos(HA_s))*r2d

night=(sun_alt<-6.0)
# Print night window for verification
t_h=t/3600
ni=np.where(night)[0]
if len(ni):
    print(f"Night window: {t_h[ni[0]]:.2f}h – {t_h[ni[-1]]:.2f}h  "
          f"(Sun min {sun_alt.min():.1f}°)")

# Survey grids — reduced for speed
a_vals=[7000,10000,15000,20000,26560,42164,80000,150000,300000,384000]
e_vals=[0.0,0.2,0.5,0.7,0.8,0.9,0.95,0.99]
i_vals=[0,45,90,135,180]        # 5 inclinations
om_vals=[0,90,180,270]          # 4
Om_vals=[0,90,180,270]          # 4
M0_vals=np.arange(0,360,45)    # 8 starting phases

GL,GH,SMAX=1.0,5.0,0.5   # deg/h thresholds
AZLO,AZHI=190.,210.; ALTMIN=5.; TMAX=4.0

found=[]; n_check=0; n_skip=0

for a in a_vals:
  for e in e_vals:
    rp=a*(1-e); ra_orb=a*(1+e)
    if rp<R_E+100: n_skip+=1; continue
    for i_deg in i_vals:
      i_r=i_deg*d2r
      for Om_deg in Om_vals:
        Om_r=Om_deg*d2r
        for om_deg in om_vals:
          om_r=om_deg*d2r
          for M0_deg in M0_vals:
            M0_r=M0_deg*d2r; n_check+=1
            try:
                alt,az,dist=run_orbit(a,e,i_r,om_r,Om_r,M0_r,t)
            except: continue
            # Angular velocity
            alt_r2=alt*d2r; ca=np.cos(alt_r2)
            daz=np.gradient(np.unwrap(az*d2r)*r2d,DH)
            dalt=np.gradient(alt,DH)
            om_app=np.sqrt((daz*ca)**2+dalt**2)
            # Masks
            ok=night&(alt>ALTMIN)&(az>=AZLO)&(az<=AZHI)
            g=ok&(om_app>=GL)&(om_app<=GH)
            s=ok&(om_app<SMAX)
            if not(np.any(g) and np.any(s)): continue
            gi_arr=np.where(g)[0]; si_arr=np.where(s)[0]
            for gi in gi_arr:
                cands=si_arr[si_arr>gi]
                if not len(cands): continue
                dt_h=(t[cands[0]]-t[gi])/3600.
                if dt_h<=TMAX:
                    found.append({'a':a,'e':e,'i':i_deg,'om':om_deg,'Om':Om_deg,
                                  'M0':M0_deg,'dt_h':dt_h,
                                  'og':om_app[gi],'os':om_app[cands[0]],
                                  'azg':az[gi],'azs':az[cands[0]],
                                  'altg':alt[gi],'alts':alt[cands[0]],
                                  'dg':dist[gi],'ds':dist[cands[0]]})
                break

print(f"\nChecked: {n_check:,}  |  Skipped (burn-up): {n_skip:,}")
print(f"Hits (guidance→stop ≤{TMAX}h in corridor): {len(found)}")
if found:
    print(f"\n{'a(km)':>8}{'e':>5}{'i°':>5}{'ω_g':>7}{'ω_s':>7}"
          f"{'az_g':>7}{'alt_g':>7}{'Δt h':>6}{'d_g km':>9}")
    for f in sorted(found,key=lambda x:x['dt_h'])[:20]:
        print(f"{f['a']:8.0f}{f['e']:5.2f}{f['i']:5.0f}  "
              f"{f['og']:6.2f} {f['os']:6.3f}  "
              f"{f['azg']:6.1f} {f['altg']:6.1f} {f['dt_h']:5.2f}  {f['dg']:8.0f}")
else:
    print("\nNo geocentric orbit satisfies all criteria simultaneously.")

# -----------------------------------------------------------------------
# Follow-up: detailed check of the hits with tighter stop threshold
# and explicit azimuth-at-stop check
# -----------------------------------------------------------------------
print("\n--- Re-check with tighter stopping threshold (ω_stop < 0.20°/h) ---")
found2=[]; n_check2=0
for a in a_vals:
  for e in e_vals:
    rp=a*(1-e)
    if rp<R_E+100: continue
    for i_deg in i_vals:
      i_r=i_deg*d2r
      for Om_deg in Om_vals:
        Om_r=Om_deg*d2r
        for om_deg in om_vals:
          om_r=om_deg*d2r
          for M0_deg in M0_vals:
            M0_r=M0_deg*d2r; n_check2+=1
            try:
                alt,az,dist=run_orbit(a,e,i_r,om_r,Om_r,M0_r,t)
            except: continue
            alt_r2=alt*d2r; ca=np.cos(alt_r2)
            daz=np.gradient(np.unwrap(az*d2r)*r2d,DH)
            dalt=np.gradient(alt,DH)
            om_app=np.sqrt((daz*ca)**2+dalt**2)
            ok=night&(alt>ALTMIN)&(az>=AZLO)&(az<=AZHI)
            g=ok&(om_app>=GL)&(om_app<=GH)
            s=ok&(om_app<0.20)   # tighter stop criterion
            if not(np.any(g) and np.any(s)): continue
            gi_arr=np.where(g)[0]; si_arr=np.where(s)[0]
            for gi in gi_arr:
                cands=si_arr[si_arr>gi]
                if not len(cands): continue
                dt_h=(t[cands[0]]-t[gi])/3600.
                if dt_h<=TMAX:
                    found2.append({'a':a,'e':e,'i':i_deg,'om':om_deg,'Om':Om_deg,
                                   'M0':M0_deg,'dt_h':dt_h,
                                   'og':om_app[gi],'os':om_app[cands[0]],
                                   'azg':az[gi],'azs':az[cands[0]],
                                   'altg':alt[gi],'alts':alt[cands[0]],
                                   'dg':dist[gi],'ds':dist[cands[0]]})
                break

print(f"Hits with ω_stop < 0.20°/h: {len(found2)}")
if found2:
    for f in sorted(found2,key=lambda x:x['dt_h'])[:10]:
        print(f"  a={f['a']:.0f} km, e={f['e']:.2f}, i={f['i']:.0f}°, "
              f"ω_g={f['og']:.2f}°/h → ω_s={f['os']:.3f}°/h, "
              f"az_g={f['azg']:.1f}° az_s={f['azs']:.1f}°, Δt={f['dt_h']:.2f}h")

# -----------------------------------------------------------------------
# Detailed track for the first hit orbit (a=42164, e=0.2, i=0)
# -----------------------------------------------------------------------
print("\n--- Detailed track for best-hit orbit (a=42164 km, e=0.2, i=0°) ---")
a_h=42164; e_h=0.2; i_h=0.
try:
    alt_h,az_h,dist_h=run_orbit(a_h,e_h,0.,0.,0.,0.,t)
    ca_h=np.cos(alt_h*d2r)
    daz_h=np.gradient(np.unwrap(az_h*d2r)*r2d,DH)
    dalt_h=np.gradient(alt_h,DH)
    om_h=np.sqrt((daz_h*ca_h)**2+dalt_h**2)
    print(f"{'t(h)':>6} {'alt°':>7} {'az°':>7} {'ω°/h':>7} {'dist km':>9} {'sun°':>7}")
    for ii in range(0,len(t),2):  # every 10 min
        th=t[ii]/3600
        if th<6.0 or th>12.0: continue  # only nighttime
        flag=""
        if sun_alt[ii]>-6: flag=" (day)"
        if 190<=az_h[ii]<=210: flag+=" [CORR]"
        print(f"{th:6.2f} {alt_h[ii]:7.1f} {az_h[ii]:7.1f} {om_h[ii]:7.2f} "
              f"{dist_h[ii]:9.0f} {sun_alt[ii]:7.1f}{flag}")
except Exception as ex:
    print(f"Error: {ex}")

print("\n--- Check: what is the Moon's angular velocity? ---")
# Moon orbit: a~384400 km, e~0.055, i~5.1°, period ~27.3 days
a_m=384400; e_m=0.055
n_m=np.sqrt(mu_E/a_m**3)
# Approximate omega at mean distance
h_m=np.sqrt(mu_E*a_m*(1-e_m**2))
om_moon_mean=h_m/a_m**2*r2d*3600  # deg/h
om_moon_peri=h_m/(a_m*(1-e_m))**2*r2d*3600
om_moon_apo=h_m/(a_m*(1+e_m))**2*r2d*3600
print(f"Moon ω at mean distance: {om_moon_mean:.3f}°/h")
print(f"Moon ω at perigee:       {om_moon_peri:.3f}°/h")
print(f"Moon ω at apogee:        {om_moon_apo:.3f}°/h")
print(f"(Guidance threshold = {GL}°/h — Moon never reaches it)")
