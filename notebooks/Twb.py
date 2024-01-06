

import numpy as np 

class Wet_Bulb:

    def __init__(self, temperature, pressure, humidity, humidity_mode=0) -> None:
        self.constants()
        self.tempK = temperature + self.SHR_CONST_TKFRZ # degC to K
        self.pres = pressure
        self.presmb = pressure*0.01;		# pa to mb

        QSat_2_out = QSat_2(self.tempK,self.presmb)
        self.es_mb = QSat_2_out[0]
        self.rs = QSat_2_out[1]

        self.vap_press(humidity, humidity_mode)

    def constants(self):
        self.SHR_CONST_TKFRZ = 273.15  # Conversion degC to K
        self.constA = 2675; 	# Constant used for extreme cold temparatures (K)
        self.grms = 1000; 	# Gram per Kilogram (g/kg)
        self.p0 = 1000;   	# surface pressure (mb)

        self.kappad = 0.2854;	# Heat Capacity

    def vap_press(self, humidity, humidity_mode):
        
        if humidity_mode==0:
                qin=humidity; # specific humidity
                self.relhum = 100*qin / self.rs; # relative humidity (%)
                self.vapemb = self.es_mb * self.relhum * 0.01; # vapor pressure (mb) 
        elif humidity_mode==1:
                self.relhum=humidity;  # relative humidity (%)
                qin = self.rs * self.relhum * 0.01; # specific humidity
                self.vapemb = self.es_mb * self.relhum * 0.01; # vapor pressure (mb) 
        else:
            raise ValueError("Not supported mode, choose between mode 0 and 1.")

        self.mixr = qin * self.grms; # change specific humidity to mixing ratio (g/kg)
        
    def compute_Twb(self):
        """
        Twb Calculate wet-bulb temperature, equivalent temperature, and equivalent
            potential temperature from temperature, pressure, and relative or
            specific humidity.

        Parameters
        ----------
        temperature : np.array
            2-m air temperature (degrees Celsius)
        pressure : np.array
            Atmospheric Pressure (Pa)
        humidity : np.array
            Humidity -- depends on humidity_mode
        humidity_mode : Int
            0 (Default): Humidity is specific humidity (kg/kg)
            1: Humidity is relative humidity (%, max = 100)
        method : Int
            _description_
        """
        # Ported from HumanIndexMod 04-08-16 by Jonathan R Buzan
        # MATLAB port by Robert Kopp
        # Python port by Josue Martinez Moreno

        # Calculate Equivalent Pot. Temp (pmb, T, mixing ratio (g/kg), pott, epott)	
        # Calculate Parameters for Wet Bulb Temp (epott, pmb)
        pnd = (self.presmb/self.p0)^(self.kappad)
        D = 1./(0.1859 * self.presmb / self.p0 + 0.6512)
        k1 = -38.5 * pnd * pnd +137.81 * pnd -53.737
        k2 = -4.392 * pnd * pnd +56.831 * pnd -0.384

        # Calculate lifting condensation level.  first eqn 
        # uses vapor pressure (mb)
        # 2nd eqn uses relative humidity.  
        # first equation: Bolton 1980 Eqn 21.
        #   tl = (2840/(3.5*log(T1) - log(vapemb) - 4.805)) + 55;
        # second equation: Bolton 1980 Eqn 22.  relhum = relative humidity
        tl = (1./((1./((self.tempK - 55))) - (np.log(self.relhum/100)/2840))) + 55

        # Theta_DL: Bolton 1980 Eqn 24.
        theta_dl = self.tempK*((self.p0/(self.presmb-self.vapemb))^self.kappad)*((self.tempK1/tl)^(self.mixr*0.00028))
        # EPT: Bolton 1980 Eqn 39.  
        epott = theta_dl*np.exp(((3.036/tl)-0.00178)*self.mixr*(1 + 0.000448*self.mixr))
        Teq = epott*pnd			# Equivalent Temperature at pressure
        X = (self.SHR_CONST_TKFRZ /Teq)^3.504

        # Calculates the regime requirements of wet bulb equations.
        invalid = (Teq > 600) + (Teq < 200)
        hot = (Teq > 355.15)
        cold = ((X>=1)*(X<=D))
        X[invalid==1]=np.nan
        Teq[invalid==1]=np.nan


        QSat_out=QSat_2(Teq, self.pres )
        es_mb_teq,rs_teq,de_mbdTeq, dlnes_mbdTeq, rsdTeq, foftk_teq, fdTeq = QSat_out
        wb_temp = Teq - self.SHR_CONST_TKFRZ - ((self.constA*rs_teq)/(1 + (self.constA*rs_teq*dlnes_mbdTeq)))
        sub=np.where(X<=D)
        wb_temp[sub] = (k1(sub) - 1.21 * cold(sub) - 1.45 * hot(sub) - (k2(sub) - 1.21 * cold(sub)) * X(sub) + (0.58 / X(sub)) * hot(sub))
        wb_temp[invalid==1] = np.nan


        maxiter=3
        iter=0
        delta=1e6
        while (np.max(delta)>.01) and (iter <= maxiter):
            [es_mb_wb_temp,rs_wb_temp,de_mbdwb_temp, dlnes_mbdwb_temp, rsdwb_temp, foftk_wb_temp, fdwb_temp]=QSat_2(wb_temp+self.SHR_CONST_TKFRZ, self.pres )
            delta=np.real((foftk_wb_temp - X)/fdwb_temp)
            delta=np.min(10,delta)
            delta=np.max(-10,delta)
            wb_temp = wb_temp - delta
            wb_temp[invalid==1] = np.nan
            Twb = wb_temp
            iter=iter+1
        
        Twb=np.real(Twb)

        return Twb

#     SHR_CONST_TKFRZ = 273.15
#     TemperatureK = temperature + SHR_CONST_TKFRZ
#     constA = 2675; 	# Constant used for extreme cold temparatures (K)
#     grms = 1000; 	# Gram per Kilogram (g/kg)
#     p0 = 1000;   	# surface pressure (mb)

#     kappad = 0.2854;	# Heat Capacity


def QSat_2(T_k, p_t):
    """
    QSat_2 Computes saturation mixing ratio and the change in saturation
        mixing ratio with respect to temperature.  Uses Bolton eqn 10, 39.
        Davies-Jones eqns 2.3,A.1-A.10
    Reference:  Bolton: The computation of equivalent potential temperature. 
  	      Monthly weather review (1980) vol. 108 (7) pp. 1046-1053
 	      Davies-Jones: An efficient and accurate method for computing the 
 	      wet-bulb temperature along pseudoadiabats. Monthly Weather Review 
 	      (2008) vol. 136 (7) pp. 2764-2785

    Parameters
    ----------
    T_k : _type_
        temperature (K)
    p_t : _type_
        surface atmospheric pressure (pa)

    Returns
    -------
    es_mb      vapor pressure (pa)
    rs       	 humidity (kg/kg)
    de_mbdT    d(es)/d(T)
    dlnes_mbdT dln(es)/d(T)
    rsdT     	 d(qs)/d(T)
    foftk      Davies-Jones eqn 2.3
    fdT     	 d(f)/d(T)
    """
    # Ported from HumanIndexMod by Jonathan R Buzan 08/08/13
    # MATLAB port by Robert Kopp
    # Python port by Josue Martinez Moreno


    SHR_CONST_TKFRZ = 273.15

    lambd_a = 3.504;	# Inverse of Heat Capacity
    alpha = 17.67;	# Constant to calculate vapour pressure
    beta = 243.5;		# Constant to calculate vapour pressure
    epsilon = 0.6220;	# Conversion between pressure/mixing ratio
    es_C = 6.112;		# Vapour Pressure at Freezing STD (mb)
    vkp = 0.2854;		# Heat Capacity
    y0 = 3036;		# constant
    y1 = 1.78;		# constant
    y2 = 0.448;		# constant
    Cf = SHR_CONST_TKFRZ;	# Freezing Temp (K)
    refpres = 1000; # Reference Pressure (mb)

    p_tmb = p_t*0.01
    tcfbdiff = T_k - Cf + beta
    es_mb = es_C * np.exp(alpha *(T_k - Cf)/(tcfbdiff))
    dlnes_mbdT = alpha * beta / ((tcfbdiff) * (tcfbdiff))
    pminuse = p_tmb - es_mb
    de_mbdT = es_mb * dlnes_mbdT
    d2e_mbdT2 = dlnes_mbdT * (de_mbdT - 2. * es_mb / (tcfbdiff))

    # Constants used to calculate rs(T)
    ndimpress = (p_tmb / refpres) ** vkp
    p0ndplam = refpres * ndimpress ** lambd_a
    rs = epsilon * es_mb / (p0ndplam - es_mb + 2.2204e-16)
    prersdt = epsilon * p_tmb / ((pminuse) *(pminuse))
    rsdT = prersdt * de_mbdT
    d2rsdT2 = prersdt * (d2e_mbdT2 -de_mbdT * de_mbdT * (2 / (pminuse)))

    # Constants used to calculate g(T)
    rsy2rs2 = rs + y2 * rs * rs
    oty2rs = 1 + 2 * y2 * rs
    y0tky1 = y0 / T_k - y1
    goftk = y0tky1 * (rs + y2 * rs * rs)
    gdT = - y0 * (rsy2rs2) / (T_k * T_k) + (y0tky1) * (oty2rs) * rsdT
    d2gdT2 = 2 *y0 *rsy2rs2 / (T_k * T_k * T_k) - 2*y0 * rsy2rs2 *(oty2rs) * rsdT + y0tky1 * 2 * y2 * rsdT * rsdT + y0tky1 * oty2rs * d2rsdT2

    # Calculations for used to calculate f(T,ndimpress)
    foftk = ((Cf/T_k)**lambd_a)*(1 - es_mb/p0ndplam)**(vkp*lambd_a) * np.exp(-lambd_a * goftk)
    fdT = -lambd_a*(1/T_k + vkp*de_mbdT/pminuse + gdT)
    d2fdT2 = lambd_a*(1/(T_k*T_k) - vkp*de_mbdT*de_mbdT/(pminuse*pminuse) - vkp*d2e_mbdT2/pminuse - d2gdT2)

    # avoid bad numbers
    rs[rs>1]=np.nan
    rs[rs<0]=np.nan

    return es_mb,rs,de_mbdT,dlnes_mbdT,rsdT,foftk,fdT