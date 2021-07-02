import CORPSE
from pylab import *
import pandas
# Note: tested with python 2.7 and pandas 0.16.2

params={
    'vmaxref':[50.0,2.0,50], #Relative maximum enzymatic decomp rates
    'Ea':[37e3,54e3,50e3],    # Activation energy
    'kC':[0.01,0.01,0.01],    # Michaelis-Menton parameter
    'gas_diffusion_exp':0.6,  # Determines suppression of decomp at high soil moisture
    'substrate_diffusion_exp':1.5,
    'minMicrobeC':1e-6,       #Minimum microbial biomass (fraction of total C)
    'Tmic':0.1,       # Microbial lifetime
    'et':0.6,          # Fraction of turnover not converted to CO2
    'eup':[0.6,0.05,0.6], # Carbon uptake efficiency
    'tProtected':75.0,    # Protected C turnover time (years)
    'protection_rate':[0.3,0.001,4.0], # Protected carbon formation rate (year-1)
}

# Note, all carbon pools are split into two isotope types. The first is labeled (root addition) and the second is unlabeled (root exudates).
litter_input_rate=array([[0.0,0.0,0.0],[0.0,0.0,0.0]]) # Per year for each chemically-defined type
fineroot_data=pandas.read_csv('data/Roots for R.csv').groupby('MidDepth').mean()

# exudation_input_rate=1e-2 # 10 g/m2 per year
exudation_input_rate=fineroot_data['fineroot _g.cm3']*13.1e-3 # 13.1 mgC/g root/year (Phillips et al 2011) to give gC/cm3/year

# 0.5mg*46.3%C/(4cm*2.5mm*2.5mm) = 0.000926 gC/cm3
root_initial_mass=array([0.25,0.75,0.0])*0.5e-3*0.463/(4*0.25**2)

depths=[15,55,95,95] #cm depth for the three root additions
exudation_rate_depths=interp(depths,exudation_input_rate.index,exudation_input_rate)
exudation_rate_depths[0]*=2 # Assume higher inputs from other sources near surface
exudation_rate_depths[1:]*=0.5
exudation_rate_depths[-1]=0
root_input_rate_depths=interp(depths,fineroot_data.index,fineroot_data['fineroot _g.cm3'])*0.0

# Steady state at about 75% protected C
cpools=[CORPSE.soil_carbon_cohort(litterC=array([root_initial_mass,[0.0,0.0,0.0]]),livingMicrobeC=array([[sum(root_initial_mass)*0.01],[0.0]]),protectedC=zeros((2,3)),params=params) for d in depths]

dt=1.0/(365.0*24)

T_data=pandas.read_csv('data/Blodgett_depth_T')
T_data['time']=pandas.to_datetime(T_data['ID'])
T=T_data.groupby(('depth','time')).mean()['Temp']

VWC_data=pandas.read_csv('data/Blodgett_depth_VWC')
VWC_data['time']=pandas.to_datetime(VWC_data['ID'])
VWC=VWC_data.groupby(('depth','time')).mean()['VWC']
porosity=VWC.max()
t=pandas.date_range(start=VWC.index[0][1],end=T.index[-1][1],freq='H')

nsteps=len(t)

outputs={}
outputs['litterC']=zeros((len(cpools),nsteps,3))
outputs['protectedC']=zeros((len(cpools),nsteps,3))
outputs['microbeC']=zeros((len(cpools),nsteps))
outputs['CO2']=zeros((len(cpools),nsteps))
outputs['decomp']=zeros((len(cpools),nsteps,3))
outputs['T_C']=zeros((len(cpools),nsteps))
outputs['VWC']=zeros((len(cpools),nsteps))

# Switch to make exudates constant with depth
constant_exudates=False
if constant_exudates:
    exudation_rate_depths[:]=exudation_rate_depths[0]
    print('Warning: Using constant exudation with depth!')

VWCdepths=VWC.index.levels[0]
VWC_array=zeros((len(t),len(VWCdepths)))
for nn in range(len(VWCdepths)):
    VWC_array[:,nn]=VWC.xs(VWCdepths[nn],level='depth').resample('H').interpolate()

Tdepths=T.index.levels[0]
T_array=zeros((len(t),len(Tdepths)))
for nn in range(len(Tdepths)):
    T_array[:,nn]=T.xs(Tdepths[nn],level='depth').resample('H').interpolate()

print('Starting simulation')
for step in range(nsteps):
    if t[step].day==1 and t[step].hour==1:
        print('Year: %d, month: %d'%(t[step].year,t[step].month))
    for dep in range(len(depths)):

        # T and theta change with depth
        Tval=interp(depths[dep],Tdepths,T_array[step,:])+273.15
        thetaval=interp(depths[dep],VWCdepths,VWC_array[step,:])/porosity
        # thetaval=0.5
        output=cpools[dep].update(Tval,thetaval,dt)
        cpools[dep].check_validity()

        cpools[dep].add_carbon(array([[0.0,0.0,0.0],
            [(exudation_rate_depths[dep]+root_input_rate_depths[dep]*0.3),root_input_rate_depths[dep]*0.7,0.0]])*dt)

        outputs['litterC'][dep,step,:]=cpools[dep].litterC[0,:]
        outputs['protectedC'][dep,step,:]=cpools[dep].protectedC[0,:]
        outputs['microbeC'][dep,step]=cpools[dep].livingMicrobeC.sum()
        outputs['CO2'][dep,step]=cpools[dep].CO2[0,:]
        outputs['decomp'][dep,step,:]=output['decomp'][0,:]
        outputs['T_C'][dep,step]=Tval
        outputs['VWC'][dep,step]=thetaval


def plot_outputs(outputs,**kwargs):
    total=outputs['microbeC'][0,:]+outputs['litterC'][0,:,:].sum(axis=1)+outputs['protectedC'][0,:,:].sum(axis=1)
    shallow_h=plot(t,outputs['litterC'][0,:,:].sum(axis=1)/total[0],'k-',label='Unprotected',**kwargs)[0]
    prot_h=plot(t,outputs['protectedC'][0,:,:].sum(axis=1)/total[0],'m-',label='Protected',**kwargs)[0]
    unprot_h=shallow_h

    total=outputs['microbeC'][1,:]+outputs['litterC'][1,:,:].sum(axis=1)+outputs['protectedC'][1,:,:].sum(axis=1)
    middle_h=plot(t,outputs['litterC'][1,:,:].sum(axis=1)/total[0],'k:',**kwargs)[0]
    plot(t,outputs['protectedC'][1,:,:].sum(axis=1)/total[0],'m:',**kwargs)

    total=outputs['microbeC'][2,:]+outputs['litterC'][2,:,:].sum(axis=1)+outputs['protectedC'][2,:,:].sum(axis=1)
    deep_h=plot(t,outputs['litterC'][2,:,:].sum(axis=1)/total[0],'k--',**kwargs)[0]
    plot(t,outputs['protectedC'][2,:,:].sum(axis=1)/total[0],'m--',**kwargs)

    total=outputs['microbeC'][3,:]+outputs['litterC'][3,:,:].sum(axis=1)+outputs['protectedC'][3,:,:].sum(axis=1)
    noexud_h=plot(t,outputs['litterC'][3,:,:].sum(axis=1)/total[0],'k-.',**kwargs)[0]
    plot(t,outputs['protectedC'][3,:,:].sum(axis=1)/total[0],'m-.',**kwargs)

    return shallow_h,prot_h,unprot_h,middle_h,deep_h,noexud_h


f=figure(1);clf()

shallow_h,prot_h,unprot_h,middle_h,deep_h,noexud_h=plot_outputs(outputs)

leg=legend(handles=[unprot_h,prot_h,shallow_h,middle_h,deep_h],labels=['Unprotected','Protected','15 cm','55 cm','95 cm','95 cm (no exud)'])
ylabel('Fraction remaining')

f.autofmt_xdate()

f.tight_layout()
draw()

f2=figure(2);clf()


o=outputs
subplot(2,2,1)
total=o['microbeC'][0,:]+o['litterC'][0,:,:].sum(axis=1)+o['protectedC'][0,:,:].sum(axis=1)
# plot(t,litterC.sum(axis=1),'k--',label='Unprotected')
plot(t,o['litterC'][0,:,1]/total[0],'g-',label='Slow')
plot(t,total/total[0],'k-',label='Total')

total=o['microbeC'][1,:]+o['litterC'][1,:,:].sum(axis=1)+o['protectedC'][1,:,:].sum(axis=1)
# plot(t,litterC.sum(axis=1),'k--',label='Unprotected')
plot(t,o['litterC'][1,:,1]/total[0],'g:')
plot(t,total/total[0],'k:')

total=o['microbeC'][2,:]+o['litterC'][2,:,:].sum(axis=1)+o['protectedC'][2,:,:].sum(axis=1)
# plot(t,litterC.sum(axis=1),'k--',label='Unprotected')
plot(t,o['litterC'][2,:,1]/total[0],'g--')
plot(t,total/total[0],'k--')

total=o['microbeC'][3,:]+o['litterC'][3,:,:].sum(axis=1)+o['protectedC'][3,:,:].sum(axis=1)
# plot(t,litterC.sum(axis=1),'k--',label='Unprotected')
plot(t,o['litterC'][3,:,1]/total[0],'g-.')
plot(t,total/total[0],'k-.')


legend(fontsize='medium',loc='upper right').get_frame().set_alpha(0.5)
title('Slow carbon pool')
ylabel('Carbon pools (kgC/m$^2$)')
xlabel('Time (years)')

draw()

subplot(2,2,2)
plot(t,o['litterC'][0,:,0]/total[0],'b-',label='Fast')
plot(t,o['litterC'][0,:,2]/total[0],'r-',label='Dead mic')
plot(t,o['microbeC'][0,:]/total[0],'c-',label='Live Mic')

plot(t,o['litterC'][1,:,0]/total[0],'b:')
plot(t,o['litterC'][1,:,2]/total[0],'r:')
plot(t,o['microbeC'][1,:]/total[0],'c:')

plot(t,o['litterC'][2,:,0]/total[0],'b--')
plot(t,o['litterC'][2,:,2]/total[0],'r--')
plot(t,o['microbeC'][2,:]/total[0],'c--')

plot(t,o['litterC'][3,:,0]/total[0],'b-.')
plot(t,o['litterC'][3,:,2]/total[0],'r-.')
plot(t,o['microbeC'][3,:]/total[0],'c-.')

plot(t,o['protectedC'][0,:,:].sum(axis=1)/total[0],'m-',label='Protected')
plot(t,o['protectedC'][1,:,:].sum(axis=1)/total[0],'m:')
plot(t,o['protectedC'][2,:,:].sum(axis=1)/total[0],'m--')

title('Other carbon pools')
ylabel('Carbon pools (kgC/m$^2$)')
xlabel('Time (years)')

legend()

subplot(2,2,3)
plot(t,o['litterC'][2,:,1]/o['decomp'][2,:,1],'b--',label='Deep')
plot(t,o['litterC'][1,:,1]/o['decomp'][1,:,1],'b:',label='Middle')
plot(t,o['litterC'][0,:,1]/o['decomp'][0,:,1],'b-',label='Shallow')
plot(t,o['litterC'][3,:,1]/o['decomp'][3,:,1],'b-.',label='Deep (no exud)')

plot(t[1:],(o['litterC'][2,1:,:].sum(axis=1)+o['protectedC'][2,1:,:].sum(axis=1))/diff(o['CO2'][2,:]/dt),'r--')
plot(t[1:],(o['litterC'][1,1:,:].sum(axis=1)+o['protectedC'][1,1:,:].sum(axis=1))/diff(o['CO2'][1,:]/dt),'r:')
plot(t[1:],(o['litterC'][0,1:,:].sum(axis=1)+o['protectedC'][0,1:,:].sum(axis=1))/diff(o['CO2'][0,:]/dt),'r-')
plot(t[1:],(o['litterC'][3,1:,:].sum(axis=1)+o['protectedC'][3,1:,:].sum(axis=1))/diff(o['CO2'][3,:]/dt),'r-.')


title('"Slow" turnover time')
ylabel('Turnover time (years)')
legend()

draw()


subplot(2,2,4)
plot(t,(VWC_array[:,0]+VWC_array[:,1])/2/porosity,'k-')
plot(t,VWC_array[:,2]/porosity,'k:')
plot(t,VWC_array[:,3]/porosity,'k--')
title('Theta (fraction of saturation)')

draw()

tight_layout()

f2.autofmt_xdate()

def steady_state_protected(protected,output):
    return (output['protected_produced']/dt)/(output['protected_turnover_rate']/protected[-1,:])

def save_csv(depth,filename,outputs=outputs):
    df=pandas.DataFrame({
        'protectedC_labeled_gC_cm-3':outputs['protectedC'][depth,:,:].sum(axis=1),
        'unprotectedC_labeled_gC_cm-3':outputs['litterC'][depth,:,:].sum(axis=1),
        'unprotectedC_slow_labeled_gC_cm-3':outputs['litterC'][depth,:,1],
        'unprotectedC_fast_labeled_gC_cm-3':outputs['litterC'][depth,:,0],
        'unprotectedC_deadmicrobe_labeled_gC_cm-3':outputs['litterC'][depth,:,2],
        'microbeC_total_gC_cm-3':outputs['microbeC'][depth,:],
        'totalC_labeled_gC_cm-3':outputs['protectedC'][depth,:,:].sum(axis=1)+outputs['litterC'][depth,:,:].sum(axis=1),
        'T_C':outputs['T_C'][depth,:],
        'VWC':outputs['VWC'][depth,:],
        'CO2prod_labeled_ugC_cm-3_hour':concatenate((diff(outputs['CO2'][depth,:])*1e6,[nan])),
    },index=t)
    df.to_csv(filename,index_label='Time')
    return df

if constant_exudates:
    save_csv(0,'CORPSE_output_constantexudates_%dcm.csv'%depths[0])
    save_csv(1,'CORPSE_output_constantexudates_%dcm.csv'%depths[1])
    save_csv(2,'CORPSE_output_constantexudates_%dcm.csv'%depths[2])
    save_csv(3,'CORPSE_output_constantexudates_%dcm_noexud.csv'%depths[3])
else:
    save_csv(0,'CORPSE_output_%dcm.csv'%depths[0])
    save_csv(1,'CORPSE_output_%dcm.csv'%depths[1])
    save_csv(2,'CORPSE_output_%dcm.csv'%depths[2])
    save_csv(3,'CORPSE_output_%dcm_noexud.csv'%depths[3])

turnover_shallow=(o['litterC'][0,1:,:].sum(axis=1)+o['protectedC'][0,1:,:].sum(axis=1))/diff(o['CO2'][0,:]/dt)
turnover_middle=(o['litterC'][1,1:,:].sum(axis=1)+o['protectedC'][1,1:,:].sum(axis=1))/diff(o['CO2'][1,:]/dt)
turnover_deep=(o['litterC'][2,1:,:].sum(axis=1)+o['protectedC'][2,1:,:].sum(axis=1))/diff(o['CO2'][2,:]/dt)
turnover_deep_noexud=(o['litterC'][3,1:,:].sum(axis=1)+o['protectedC'][3,1:,:].sum(axis=1))/diff(o['CO2'][3,:]/dt)

print( "Turnover rates:")
print( "Before June 2014, After June 2014:")
tcutoff=datetime.datetime(2014,6,1)
print( "Shallow: %1.2f, %1.2f"%(turnover_shallow[t[1:]<tcutoff].mean(),turnover_shallow[t[1:]>=tcutoff].mean()))
print( "Middle: %1.2f, %1.2f"%(turnover_middle[t[1:]<tcutoff].mean(),turnover_middle[t[1:]>=tcutoff].mean()))
print( "Deep: %1.2f, %1.2f"%(turnover_deep[t[1:]<tcutoff].mean(),turnover_deep[t[1:]>=tcutoff].mean()))
print( "Deep, no exudation: %1.2f, %1.2f"%(turnover_deep_noexud[t[1:]<tcutoff].mean(),turnover_deep_noexud[t[1:]>=tcutoff].mean()))

show()
