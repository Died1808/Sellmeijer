"""
Dit script genereert op basis van een sondering (GEF-file) een
visualisering van pipelengtes op basis van Sellmeijer
in combionatie met TIM.
Het is bedoeld als voorbeeld van het toepassen van TIM in
berekeningen, maar kan mijns inziens ook in de praktijk
worden toegepast voor een snelle Sellmeijer inidcatie 
als er een sondering beschiokbaar is.
Wel wordt gewezen op de veelheid van perameters die met lokale gebiedskennis
dienen te worden aangepast.

"""

import sys,os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from timml import *
import math
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

#############################################################
pd.set_option('display.max_columns', None) 
pd.options.display.width=None
pd.set_option('display.max_rows',None) 
pd.set_option('display.expand_frame_repr', True)
pd.options.display.float_format = '{:,.1f}'.format

params = {'font.family': 'sans-serif',
          'font.sans-serif': 'arial',
          'axes.labelsize': 10,
          'axes.facecolor': '#ffffff', 
          'axes.labelcolor': 'black',
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'lines.linewidth': 1,
          'grid.color': 'grey',
          'grid.linestyle': 'dashed',
          'grid.linewidth': 0.5,
          'text.usetex': False,
          'font.style': 'normal',
          'font.variant':'normal',
          'figure.facecolor': 'white',
          'font.size':8,
          'figure.autolayout': True,
          'figure.figsize': (8,14),
          'figure.dpi': 100,
          }
plt.rcParams.update(params)
#

class GEF:
    def __init__(self):
        self._data_seperator = ' '
        self._columns = {}
        self.x = 0.
        self.y = 0.
        self.z = 0.
        self.dz = []
        self.qc = []
        self.pw = []
        self.wg = []
        self.c  = []
        self.kb = []
        self.grind = []
        self.grof  = []
        self.matig = []
        self.matigfijn=[]
        self.fijn  = []
        self.leem  = []
        self.klei  = []
        self.veen  = []
        self.npor  = []
#        
    def readFile(self, filename):
        lines = open(filename, 'r').readlines()
        for line in lines:
            reading_header = True
        for line in lines:   
            if reading_header:
                self._parseHeaderLine(line)
            else:
                self._parseDataLine(line)
            if line.find('#EOH') > -1:
                if self._check_header():
                    reading_header = False
                else:
                    print(filename,'fout in GEF')
                    return
            
    def _check_header(self):
        if not 1 in self._columns:
            return False
        if not 2 in self._columns:
            return False
        return True
    
    def _parseHeaderLine(self, line):
        for xe in ['#COMMENT', 'Peil=', 'uitvoerder', 'materieel','WATERSTAND',
                    'opmerkingen#MEASUREMENTTEXT','==','= NAP','antropogeen']:
            if xe in line:
                return          
        if len(line.split()) == 0:
            return
        
        keyword, argline = line.split('=')         
        keyword = keyword.strip()
        argline = argline.strip()
        args = argline.split(',')

        if keyword=='#XYID':
            if float(args[1]) < 1e5:
                args[1] = args[1]
            else:
                args[1]=args[1].replace('.','')
            self.x = float(args[1])
 
            args[2]=args[2].replace('.','')
            self.y = float(args[2])
   
            if (len(str(int(self.x))))>5:
                 self.x=int(self.x/pow(10,len(str(int(self.x)))-6))
                 
            if (len(str(int(self.y))))>5:
                 self.y=int(self.y/pow(10,len(str(int(self.y)))-6))

            if self.x > 3e5:
                self.x=self.x/10
 
        elif keyword=='#ZID':
            self.z = float(args[1])
           
        elif keyword=='#COLUMNINFO':
            column = int(args[0])
            dtype = int(args[-1])
            if dtype==11:
                dtype = 10
            self._columns[dtype] = column - 1  
            
    def _parseDataLine(self, line):
        line=line.strip()
        line = line.replace('|',' ')
        line = line.replace(';',' ')
        line = line.replace('!',' ')
        line = line.replace(':',' ')
        args=line.split()
        for n, i in enumerate(args):
            if i in ['9.9990e+003','-9999.9900','-1.0000e+007','-99999','-99999.0',
                  '-99999.00','-99999.000','-9.9990e+003','999.000', '-9999.99', '999.999',
                  '99','9.999','99.999', '999.9']:
                args[n] = '0.1'       
        if len(line.split()) == 0:
            return
        z  = abs(float(args[self._columns[1]]))

        try:
            if self.z==0:
                self.z=basisz(self.x,self.y)
        except IndexError:
            return
        
        dz = self.z - z
        qc = float(args[self._columns[2]])
                         
        pw = float(args[self._columns[3]])
        if pw < 1e-3:
            pw=0.003
        elif pw>20:
            pw=20   
        self.dz.append(dz)
        self.qc.append(qc)
        self.pw.append(pw)
        
        if qc<=0.001:
            qc=0.1
            self.wg.append(10.)
        else:
            wg = abs((pw / qc) * 100.)
        wg = abs((pw / qc) * 100.)
##########        
        if wg>=0.0:
            if wg >20: wg=20
            ke=math.exp(wg)
        if ke <=0:  ke=1
        else:
            kb  = (qc / ke)*2
            self.kb.append(kb)

            c=1/(kb)
            if c > 15: c=15.
            if wg>4: c=c/10
            else:
                c=c*2
            self.c.append(c)       

        if wg > 2.7 and qc<0.32: veen = 0
        else: veen = -100
        if kb > 50.: grind = 0
        else: grind = -100 
        self.grind.append(grind)
        if kb <= 50. and kb > 20.: grof = 0
        else: grof = -100 
        self.grof.append(grof)
        if  kb <= 20. and kb > 10.: matig = 0
        else: matig = -100
        self.matig.append(matig)
        if  kb <= 10. and kb >5. : matigfijn = 0
        else: matigfijn = -100
        self.matigfijn.append(matigfijn)
        if  kb <= 5. and kb >1: fijn = 0
        else: fijn = -100
        self.fijn.append(fijn)
        if kb<=1 and kb>0.1: leem = 0
        else: leem=-100
        self.leem.append(leem)
        if veen== 0: klei = -100
        else:
            if  kb <= 0.1 and kb > 0: klei = 0
            else: klei = -100
        self.veen.append(veen)   
        self.klei.append(klei)
        
        if kb <=0.1:
            npor=0.05
        elif 0.1<kb<=1:
            npor=0.1
        elif 1<kb<=30:
            npor=0.35
        elif kb>30:
            npor=0.25
        self.npor.append(npor)

    def asNumpy(self):
        return np.transpose(np.array([self.dz, self.c, self.kb, self.grind, self.grof, self.matig, 
                                      self.matigfijn, self.fijn, self.leem, self.klei, self.veen, self.npor]))

    def asDataFrame(self):
        a = self.asNumpy()
        return pd.DataFrame(data=a, columns=['depth','c', 'k','gr', 'g','m','mf', 'f', 'le','kl','v','por'])
        
    def plt(self, filename):
        df = self.asDataFrame()
        df = df.sort_values('depth', ascending=False)
        if df.empty:
            print(filename,'fout in GEF')
            return df
        
        OZWVP = -35
        k_OZ  = 40
        df = df.iloc[:-5 , :]
        try:
            dzend=round(df.loc[df.index[-1],'depth'],1)
        except IndexError:
            return
        dzendfig=dzend
            
        if OZWVP < dzend:
            dzend=OZWVP
        
        df.at[df.index[-1],'k']= k_OZ
        ymax=int((round(self.z,-1))+10)
        ymin=ymax-30
        ds=1 #dit bepaalt de afstand op de y-as in de grafiek
# 
        df=df[df.depth <=ymax]
        df=df[df.depth >=ymin]
        if dzend <= ymin:
            dzend = ymin
        if df.empty:
            print(filename,'fout in GEF')
            return
        try:
            dist1=round(df.iloc[4,0]-df.iloc[5,0],3)
            dist2=round(df.iloc[-3,0]-df.iloc[-2,0],3)
            dist = max(dist1, dist2)
            dataq=np.array(np.arange(ymax,dzend,dist*-1))
        except (IndexError, ValueError):
            print(filename,'fout in GEF')
            return
        
        df = df.reset_index(drop=True)
        zz=self.z
#
        filename=filename.replace('.gef', '')
        filename=filename.replace('.GEF', '')
        filename=filename.replace(' ', '_')

        df = df.reset_index(drop=True)
        df = df[df['depth']>dzendfig]

#######Timml
        factor = 2 # 2 rekent 10* zo snel als factor = 1, 10 100* (er is ergens een exp in)
        df = df.rolling(factor).mean() 
        df = df.iloc[::factor, :]

        df = df.dropna()
        df = df.reset_index(drop=True)
        
        df['kzkh']= 2/(33.653*np.exp(-0.066*(df['k'])))
        df['kzkh']= np.where(df['kzkh']>1,1,df['kzkh'])
        dl=len(df['depth'])
        dfd=df['depth'].tolist()
        dfd.append(OZWVP)

#######################################
        MWH          = 10.2     #  Maatgevende waterstand
        maaiveld     = 4        #  laagste punt, bijv in watergang wpeil = 4.85
        diktesell    = 36       ## Dikte WVP Regis
        dikte_dek    = 1.5      ## Op basis gekarteerde dikte 
        L2 = int(170/factor)    # pipepunt uittrede, alleen door de deklaag
        naam = 'Tuil_extreem_10.2'
        print('naam = ', naam)

######################################
        dfn = df.iloc[: , [1]].copy()   
        dfn=df.append({'npor':0.25}, ignore_index=True) 

        ml=Model3D(kaq=df['k'], z=dfd, kzoverkh=df['kzkh'], npor=dfn['npor'], 
                    topboundary='semi', topres=2, topthick=0.1, hstar=maaiveld) 

        ### RD pipe uittree fictief V2
        pipx = 144616
        pipy = 425810
        uittree= HeadWell(ml, pipx, pipy, maaiveld, rw=0.1, layers = np.arange(0,L2,1))

        # ### RD intree 
        pipex  = 144621
        pipey  = 425746
        intree= HeadWell(ml, pipex, pipey, MWH, rw=1,res=0, layers = np.arange(0,L2,1))
        
        l_entry_exit = round(((pipex-pipx)**2+(pipey-pipy)**2)**0.5,1)

        ml.solve()

        df['in']   = np.round(ml.head(pipex,pipey, df.index),2)    
        df['uit']  = np.round(ml.head(pipx,pipy, df.index),2)  

        """
        De hieronder gegenereerde waarden komen uit v Beek 1001449-008-GEO-0001SBW_Hervalidatie_piping_-_B3._Analyse_lab_proeven.pdf
        """
        df['Dr']   = 0.7 # relatieve dichtheid (percentage optimaal bij 70%)
        df['Cu']   = 2.4806*(pow(df['k'],-0.005))  # Excelsheet WAM KGA uniformiteits coefficient = D60/D10= gradatie


######### Sellmeijer_Pol
        n            = 0.34        # whites constant [-]
        vw           = 1.5182e-6   # kinematic viscosity bij 5 Celcius
        g            = 9.81        # gravity [m/s2]
        ysp          = 26.5-9.81   # vol gewicht zand korrels onder water
        D70factor    = 1.4e-4      # neem nu eens minimum in Hjulstrom=fijn zand wat als eerste uitstroomt
        KAS          = 60
        rho          = 38
        
        ## weerstandsfactor komt rond 0.3      
        df['fr']     = n * (ysp / 9.81) * np.tan(np.radians(rho))*(pow(df['Dr']/0.725,0.35))* (pow(df['Cu']/1.81,0.13))*(pow(KAS/0.498,-0.02))

        # Geometriefactor komt rond 1.1
        ddl          = diktesell / l_entry_exit
        asell        = pow(ddl, 2.8) - 1
        b            = 0.28 / asell + 0.04
        df['fg']     = 0.91 * pow(ddl, b) 
        df['kms']    = df['k']/(24*60*60)      
        df['kappa']  = vw / g * df['kms']
        df['d70']    = 0.96984815*df['kms'] + 0.00012985
        df['csell']  = D70factor/ pow(df['kappa'] * l_entry_exit, 1/3)
        df['dsell']  = df['d70']/D70factor
        
          ## Schaalfactor komt rond 0.2
        df['fs']    = df['csell'] * pow(df['dsell'],0.6) # ipv 0.4

        df['Hc'] = l_entry_exit * df['fr'] * df['fs'] * df['fg']          
        df['H'] = (df['in']-df['uit'])-0.3*dikte_dek
       
        if diktesell > l_entry_exit:
            ratioLc_L = (-0.3984 * pow(1.0, 2.0) + 0.7979*(1.0) + 0.0707)
        else: 
            ratioLc_L = (-0.3984 * pow(diktesell / l_entry_exit, 2.0) + 0.7979*(diktesell / l_entry_exit) + 0.0707)
        Lc = l_entry_exit*ratioLc_L
        Lca = Lc
        df['Lpipe'] = np.where(df['H']/df['Hc'] >0,Lc*(164.535456 * pow(df['H']/ df['Hc'], 8.0) - 570.088738 * pow( df['H']/ df['Hc'], 7.0)+ 
                            797.636538 * pow(df['H']/ df['Hc'], 6.0) - 575.466330 * pow( df['H']/ df['Hc'], 5.0) + 227.705051 * 
                            pow(df['H']/ df['Hc'], 4.0) - 48.9853799 * pow( df['H']/ df['Hc'], 3.0) + 5.96645859 * 
                            pow(df['H']/ df['Hc'], 2.0) -3.67385412e-01 * (df['H']/ df['Hc']) + 0.00757329545), Lc*0.001)
        df['Lpipe'] = np.where(df['Lpipe']>1e3,np.nan,df['Lpipe'])
 
        """
        Toelaatbaar verval Hc in bereik pipegroei   
        """
        df['hckrit']=np.where(df['Lpipe']>1, df['Hc'],np.nan)
        df['hckrit']=np.where(df['hckrit']==-1, np.nan, df['hckrit'])
        Hcmin=df['hckrit'].min()
        df['d70_avg']=np.where(df['Lpipe']>1, df['d70'],np.nan)
        try:
            d70_gem=int(df['d70_avg'].max()*1e6)
        except ValueError:
            d70_gem=np.nan

###############################################################################
        fig, ax = plt.subplots(nrows=2, ncols=2, gridspec_kw = {'width_ratios':[1,1], 'height_ratios':[1,0.3]})
        fig.subplots_adjust(wspace=0, hspace=0.2)
#
        df.plot(x='kl',  y='depth', ylim = (ymin, ymax), xlim = (0,80), yticks=range (ymin, ymax, ds), xticks=range (0,80,5), xerr=1,   elinewidth=0.5*factor, ax=ax[0,0], legend=False, color='darkgreen');
        df.plot(x='v',   y='depth', ylim = (ymin, ymax), xlim = (0,80), yticks=range (ymin, ymax, ds), xticks=range (0,80,5), xerr=1,   elinewidth=0.5*factor, ax=ax[0,0], legend=False, color='mediumorchid');
        df.plot(x='le',  y='depth', ylim = (ymin, ymax), xlim = (0,80), yticks=range (ymin, ymax, ds), xticks=range (0,80,5), xerr=1,   elinewidth=0.5*factor, ax=ax[0,0], legend=False, color='lightgreen');
        df.plot(x='f',   y='depth', ylim = (ymin, ymax), xlim = (0,80), yticks=range (ymin, ymax, ds), xticks=range (0,80,5), xerr='k', elinewidth=0.5*factor, ax=ax[0,0], legend=False, color='gold');
        df.plot(x='mf',  y='depth', ylim = (ymin, ymax), xlim = (0,80), yticks=range (ymin, ymax, ds), xticks=range (0,80,5), xerr='k', elinewidth=0.5*factor, ax=ax[0,0], legend=False, color='orange');
        df.plot(x='m',   y='depth', ylim = (ymin, ymax), xlim = (0,80), yticks=range (ymin, ymax, ds), xticks=range (0,80,5), xerr='k', elinewidth=0.5*factor, ax=ax[0,0], legend=False, color='darkgoldenrod');
        df.plot(x='g',   y='depth', ylim = (ymin, ymax), xlim = (0,80), yticks=range (ymin, ymax, ds), xticks=range (0,80,5), xerr='k', elinewidth=0.5*factor, ax=ax[0,0], legend=False, color='firebrick');       
        df.plot(x='gr',  y='depth', ylim = (ymin, ymax), xlim = (0,80), yticks=range (ymin, ymax, ds), xticks=range (0,80,5), xerr='k', elinewidth=0.5*factor, ax=ax[0,0], legend=False, grid=True, color='darksalmon',  title=filename);       

        if df['Lpipe'].max() > l_entry_exit:
            color = 'xkcd:claret'
        else:
            color = 'navy'
            
        df.plot(x='Lpipe',   y='depth', ylim = (ymin, ymax), xlim = (0.1,1000), yticks=range (ymin, ymax, ds), xticks=range(1,100000,10), ax=ax[0,1], legend=False, grid=True, color=color, title='Pipelengte o.b.v. Sellmeijer');       

        ax[0,1].set_xscale('log')
        ax[0,1].grid(True, color='lightgrey', which='minor')
        ax[0,1].set_yticklabels([])
        ax[0,1].set_xlabel('Lengte pipe in [m]  bij rivier op '+str(MWH)+ ' [m NAP]')
        ax[0,1].vlines(x = l_entry_exit, ymin=ymin , ymax=self.z , linewidth = 2, colors = 'red')
        ax[0,1].fill_betweenx(np.arange(df.iloc[0,0], df.iloc[-1,0]-0.02,-dist*factor), 0.1 ,
                            df['Lpipe'], color=color, alpha=0.2)
        
        ######## Maken Legenda         
        k = mpatches.Patch(color='darkgreen',    label='Klei                       k <= 0.1  [m/dag]')
        l = mpatches.Patch(color='mediumorchid', label='Veen                       k <= 1    [m/dag]')
        f = mpatches.Patch(color='lightgreen',   label='Leem tot fijn zand   0.1 < k <= 1    [m/dag]')
        mf= mpatches.Patch(color='gold',         label='Fijn zand              1 < k <= 5    [m/dag]')
        m = mpatches.Patch(color='orange',       label='Matig fijn/grof zand   5 < k <= 10   [m/dag]')
        mg= mpatches.Patch(color='darkgoldenrod',label='Matig grof zand       10 < k <= 20   [m/dag]')
        g = mpatches.Patch(color='firebrick',    label='Grof zand tot grindig 20 < k <= 50   [m/dag]')
        gr= mpatches.Patch(color='darksalmon',   label='Grindig zand tot grind     k >  50   [m/dag]')
        
        leg=ax[0,0].legend(handles=[k,f,l,mf,m,mg,g,gr], prop={'family': 'monospace'},
        bbox_to_anchor=(2.7,0.10), loc='center',fontsize=12, frameon=1, title='K-waarde', 
        title_fontsize=10, labelspacing=0.5, borderpad=1, handlelength=2, handletextpad=2);
        leg._legend_box.sep = 10
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)
        
        ax[0,0].text(0, ymin-3,   'RD_x = '+ str(int(self.x)), fontsize=10)
        ax[0,0].text(0, ymin-4.5, 'RD_y = '+ str(int(self.y)), fontsize=10)
        ax[0,0].hlines(self.z, xmin=0, xmax=80, linestyle='solid', color='grey')
        ax[0,0].text(1, self.z+0.3, 'Mv', weight='bold')
        ax[0,1].hlines(self.z, xmin=0.1, xmax=1e3, linestyle='solid', color='grey')

######## Monte Carlo simulatie
        N     =  10000
        RDm   =  0.725
        gz    =  9.81
        fct   =  0.5 # Voorspelling Deze waarde ligt tussen 0-1, bijvoorbeeld wat is Head bij een waarschijnlijkhein van 0.95?
        count =  0
        hrd   =  5
        t     =  60*60*24
        vw    =  1.3063e-6 # waarde 10 celsius
        rho   =  38
        RD    =  0.7
        D     =  diktesell
        L     =  l_entry_exit
        d70   =  df.d70.quantile(0.85)
        k     =  50/t  #Omrekenen naar [m/sec]
        
 
        def triangle(min, most, max):
             return np.random.triangular(min, most, max, N) 
        
        # Sellmeijer, van Beek backward erosion formule
        def pipe(Fr, Fs, Fg, L, Lc,druk):
            H  = druk
            Hc = Fr*Fg*Fs*L
            Lpipe = Lc*(164.535456 * pow(H/ Hc, 8.0) - 570.088738 * pow( H/ Hc, 7.0)+ 
                        797.636538 * pow(H/ Hc, 6.0) - 575.466330 * pow( H/ Hc, 5.0) + 227.705051 * 
                        pow(H/ Hc, 4.0) - 48.9853799 * pow( H/ Hc, 3.0) + 5.96645859 * 
                        pow(H/ Hc, 2.0) -3.67385412e-01 * (H/ Hc) + 0.00757329545)
            return (Lpipe) 	

        def head(Fr, Fs, Fg, L):
              return (Fr*Fs*Fg*L) 	
        
        """input data voor MC analyse:
            
        Dit is de een mogelijke Hc waarde waarop de terugschrijdende erosie plaats vindt.
        
        D variatie dikte zandlaag onder deklaag
        L variatie lengte tussen intrede een uittredepunt
        RD variatie dichtheid zand
        d70 variatie korrelgrootte
        k variatie k-waarde (bulk) in m/sec (t = omrekenwaarde)
        vw variatie temperatuur 
        rho variatie rolweerstand 
        """
        df=df[df['k']>1]
        
        D   = triangle(diktesell-2, diktesell, diktesell+2)
        L   = triangle(l_entry_exit-2.5, l_entry_exit, l_entry_exit+2.5)
        RD  = triangle(0.6, 0.7, 0.8)	
        d70 = triangle(df.d70.quantile(0.14), df.d70.quantile(0.5), df.d70.quantile(0.95))
        k   = triangle(1/t, df.k.quantile(0.5)/t, df.k.quantile(0.995)/t)
        rho = triangle(36, 37, 38)
        vw  = triangle(1.1386e-6, 1.3063e-6, 1.5182e-6) ## 15, 10 en 5 celsius
        
        ## Berekening Sellmeijer obv Monte-Carlo simulatie
        Fr = n*ysp/gz*np.tan(np.radians(rho))*(pow(RD/0.725,0.35))	
        Fs = d70/(vw/gz*k*L)**(1/3)*(D70factor/d70)**0.6
        Fg = 0.91*(D/L)**(0.28/((D/L)**2.8-1)+0.04)
        ratioLc_L = (-0.3984 * pow(D/L, 2.0) + 0.7979*(D/L) + 0.0707)
        Lc = l_entry_exit*ratioLc_L
        druk = MWH-maaiveld-0.3*dikte_dek
        
        P = pipe(Fr, Fs, Fg, L, Lc, druk)
        H = head(Fr, Fs, Fg, L)

        ## plotten        
        xx=np.sort(H)
        Hmin = round(np.quantile(xx,.05),1)
        Hmax = round(np.quantile(xx,.95),1)
        
        xxp=np.sort(P)
        print('percentage kleiner dan toegestane kwelweglengte',  round((len(xxp[xxp < l_entry_exit])/N)*100,1))
        Pmin = round(np.quantile(xxp,.05),1)
        Pmax = round(np.quantile(xxp,.95),1)

        tfont = {'fontname':'monospace'}
        
        ##   plaatjes bij MC        
        yy = np.arange(N) / float(N)
        yyp = np.arange(N) / float(N)

        for i in range(N):
            if yy[i] == fct:
                Hmed= round(xx[i],1)
            if xx[i] > hrd:
                count = count + 1
        for i in range(N):
            if yyp[i] == fct:
                Pmed= round(xxp[i],1)
            if xxp[i] > hrd:
                count = count + 1
        # ax[1,0].set_visible(False)
        ax[1,0].set_xlabel('Hc [m]')
        ax[1,0].set_ylabel('Frequentie')
        ax[1,0].set_title('Verdeling Hc [m], Monte-Carlo simulatie')
        ax[1,0].hist(H, bins=200, density=True, color='tomato')
        ax[1,0].set_xlim(0,20)
        ax[1,0].set_xticks(np.arange(0, 20, step=1))
        ax[1,0].set_ylim(0,0.7)
        ax[1,0].grid()
 
        ax[1,1].set_xlim(0,250)
        ax[1,1].set_xlabel('Lpipe [m]')
        ax[1,1].set_ylabel('Frequentie')
        ax[1,1].set_title('Verdeling Lpipe [m], Monte-Carlo simulatie '+'N= '+str(N))
        ax[1,1].hist(P, bins=N, density=True, color='darkblue')
        ax[1,1].set_xticks(np.arange(0,250,25))
        
        ax[1,1].vlines(x = l_entry_exit, ymin=0, ymax= 0.015, linewidth = 2, colors = 'red')

        ax[1,1].yaxis.set_label_position("right")
        ax[1,1].yaxis.set_ticks_position("right")
        
        ax[1,1].grid()

        ax[0,0].text(165, ymax,   'Toelaatbare kwelweglengte  '+ str(l_entry_exit) + ' [m]', **tfont)
        ax[0,0].text(165, ymax-1, 'Maximale pipelengte        '+ str(round(df['Lpipe'].max(),1)) + ' [m]', **tfont)
        ax[0,0].text(165, ymax-2, 'Rivierwaterstand           '+ str(MWH)+ ' [m NAP]', **tfont)
        ax[0,0].text(165, ymax-3, 'Hoogte maaiveld achterland '+ str(maaiveld)+ ' [m NAP]', **tfont)
        ax[0,0].text(165, ymax-4, 'Dikte WVP                  '+ str(diktesell)+ ' [m]', **tfont)
        ax[0,0].text(165, ymax-5, 'Dikte deklaag              '+ str(dikte_dek)+ ' [m]', **tfont)
        ax[0,0].text(165, ymax-6, 'Toelaatbaar verval (Hc)    '+ str(round(Hcmin,2))+ ' [m]', **tfont)
        ax[0,0].text(165, ymax-7, 'Berekende Lc               '+ str(round((Lca),1))+ ' [m]', **tfont)
        ax[0,0].text(165, ymax-8, 'D70 (max. pipinginterval)  '+ str(d70_gem)+ ' [mu]', **tfont)

        ax[0,0].text(165, ymax-10,'MC_sim Hc (5%_50%_95%)    '+ str([Hmin, Hmed, Hmax]), **tfont)
        ax[0,0].text(260, ymax-10,' [m]' , **tfont)
        ax[0,0].text(165, ymax-11,'MC_sim LPipe (95%_50%_5%) '+ str([Pmax, Pmed, Pmin]), **tfont)
        ax[0,0].text(260, ymax-11,' [m]' , **tfont)

        tfont = {'fontname':'monospace','size':6, 'color':'blue'}
        ax[0,0].text(185, ymin-8,'Monte_Carlo simulatie, [min, mediaan, max]', **tfont)
        ax[0,0].text(185, ymin-10,'MC_sim D               '+ str([diktesell-5, diktesell, diktesell+5]), **tfont)
        ax[0,0].text(185, ymin-11,'MC_sim L               '+ str([l_entry_exit-5, l_entry_exit, l_entry_exit+5]), **tfont)
        ax[0,0].text(185, ymin-12,'MC_sim RD              '+ str([0.6, 0.7,0.8]), **tfont)
        ax[0,0].text(185, ymin-13,'MC_sim d70             '+ str([int(df.d70.quantile(0.14)*1e6), int(d70_gem), int(df.d70.quantile(0.95)*1e6)]), **tfont)
        ax[0,0].text(185, ymin-14,'MC_sim k               '+ str([1, int(df.k.quantile(0.5)), int(df.k.quantile(0.995))]), **tfont)
        ax[0,0].text(185, ymin-15,'MC_sim rolweerstand    '+ str([36,37,38]), **tfont)
        ax[0,0].text(185, ymin-16,'MC_sim temp            '+ str([5,10,15]), **tfont)
        ax[0,0].text(185, ymin-20,'Percentage MC_sim >Lpipe ' +str(round(100-(len(xxp[xxp < l_entry_exit])/N)*100,1)), **tfont)

        ax[0,0].text(250, ymin-10,' [m]' , **tfont)
        ax[0,0].text(250, ymin-11,' [m]' , **tfont)
        ax[0,0].text(250, ymin-12,' [-]' , **tfont)
        ax[0,0].text(250, ymin-13,' [mu]', **tfont)
        ax[0,0].text(250, ymin-14,' [m/dag]', **tfont)
        ax[0,0].text(250, ymin-15,' [graden]', **tfont)
        ax[0,0].text(250, ymin-16,' [graden Celsius]', **tfont)
        ax[0,0].text(250, ymin-20,' [%]', **tfont)

        ax=ax[0,0].set(xlabel="k [m/dag]", ylabel="[m NAP]") 

        plt.savefig(naam+".png", bbox_inches = "tight")
        plt.show()
        plt.close()

for filename in os.listdir(os.getcwd()):
    if filename.endswith ('.GEF') or filename.endswith ('.gef'):
        if __name__=="__main__":
            g=GEF()
            g.readFile(filename)
            g.plt(filename)
        plt.close()
       