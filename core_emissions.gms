* EMISSIONS MODULE
*
* Where Regions emissions are determined.
*____________
#=========================================================================
*   ///////////////////////       SETTING      ///////////////////////
#=========================================================================
#
## CONF
#_________________________________________________________________________
* Definition of the global flags and settings specific to the module
$ifthen.ph %phase%=='conf'

* MIU linear-transition time horizon from 1 to maximum upperbound
$setglobal t_min_miu 9
$setglobal t_max_miu 38

* MIU maximum reacheable upperbound
$setglobal max_miuup 1.2

* Carbon-intensity transition curve
* | linear_pure | linear_soft | sigmoid_HHs | sigmoid_Hs | sigmoid_Ms | sigmoid_Ls | sigmoid_LLs |
$setglobal sig_trns_type 'sigmoid_Ls'

* Time of full-convergence to dice-ref carbon-intensity curve
* | 28 | 38 | 48 | 58 |
$setglobal sig_trns_end  '38'

* SSP-n hypothesis on dice-reference curve for carbon-intensity
* |original | discounted |
$setglobal sig_dice_ref_curve 'discounted'


## SETS
#_________________________________________________________________________
$elseif.ph %phase%=='sets'

SET ere 'Emissions-related entities'/
    co2
    co2ffi # Fossil-fuel and Industry CO2
    nip    # net import of permits
    sav    # saved permits
    kghg   # Kyoto greenhouse gases
    ch4
    n2o
    sf6
    #hfc
    #pfc
/;
ALIAS(ere,eere);

SET map_e(ere,eere) 'Relationships between Sectoral Emissions' /
    co2.co2ffi
/;

SET ghg(ere) 'Green-House Gases'
/
    co2
    ch4
    n2o
    #sf6
    #hfc
    #pfc
/;

SET oghg(ghg) 'Other GHGs'
/
    ch4
    n2o
    #sf6
    #hfc
    #pfc
/;

## INCLUDE DATA
#_________________________________________________________________________
$elseif.ph %phase%=='include_data'

# .......  PARAMETERS HARDCODED OR ASSIGNED  .......
PARAMETERS
* Cumulative emissions startings
   cumeind0   'Starting value of cumulative emissions from industry [GtC]'                 / 400   /  # DICE2016
   cumetree0  'Starting value of cumulative emissions from land use deforestation [GtC]'   / 100   /
* Availability of fossil fuels
   fosslim    'Maximum cumulative extraction fossil fuels (GtC)' / 6000 /

* Curve parameters
*   k_curve         /0.6/
*   j_curve         /6.2/
*   r_curve         /0.5/
   l_curve(n)   '[$/(day*person)]'

* MIU rate controls
   miu0       'Initial emissions control rate for base calib_emissions'    / 0 / 
   min_miu     'upper bound for control rate MIU at t_min_miu'       / 1.00 / # best compromise
   t_min_miu    'time t when min_miu value can be reached'           / %t_min_miu%    / # 7 - 2045
   max_miu     'upper bound for control rate MIU from t_max_miu'     / %max_miuup%  / # the old DICE limmiu
   t_max_miu    'time t when max_miu value can be reached'           / %t_max_miu%  / # 28 - 2150
   min_miuoghg 'upper bound for control rate MIU at t_min_miu'       / 0.7  / # best compromise
   max_miuoghg 'upper bound for control rate MIU from tmax'        / 1    / # best compromise
;

##  PARAMETERS EVALUATED ----------
PARAMETERS
    world_e(t)
    world_eoghg(oghg,t)
;

SCALAR
* Conversion coefficients
   CtoCO2        'conversion factor from Carbon to CO2'
   CO2toC        'conversion factor from CO2 to Carbon'
;

## BAU EMISSIONS AND CARBON INTENSITY 
PARAMETERS
    ssp_emi_bau(ssp,t,n)   'SSP-Decline rate of decarbonization according to different scenarios (per period)'
    emi_bau_co2(t,n)       "Baseline regional CO2 FFI emissions [GtCO2]"
    sigma(t,n)                  'Decline rate of decarbonization (per period)'
    sig0(n)                     'Carbon intensity at starting time [kgCO2 per output 2005 USD]'
    veg0(n)                   'Initial number of vegans(%)'                      
    emibmk(n)                 'Emissions pro_capite by a benchmark    [GtCO2 per year per milion people]'           
    emiveg(n)                 'Emissions pro_capite by a vegetarian    [GtCO2 per year per milion people]'         

$gdxin '%datapath%data_baseline_emissions_calibrated.gdx'
$load   ssp_emi_bau=emi_bau_calibrated
$gdxin
;
* Configuration settings determine imported scenario
emi_bau_co2(t,n) = ssp_emi_bau('%baseline%',t,n)  ;


veg0('cajaz') = 0.06;
veg0('china') = 0.06;
veg0('easia') = 0.05;
veg0('kosau') = 0.043;
veg0('laca') = 0.1;
veg0('india') = 0.295;
veg0('mena') = 0.05;
veg0('neweuro') = 0.035;
veg0('oldeuro') = 0.048;
veg0('sasia') = 0.07;
veg0('ssa') = 0.07;
veg0('te') = 0.02;
veg0('usa') = 0.042;

emibmk('cajaz') = 537*(1e-6);
emibmk('china') = 914*(1e-6);
emibmk('easia') = 914*(1e-6);
emibmk('kosau') = 974*(1e-6);
emibmk('laca') = 914*(1e-6);
emibmk('india') = 1121*(1e-6);
emibmk('mena') = 914*(1e-6);
emibmk('neweuro') = 974*(1e-6);
emibmk('oldeuro') = 537*(1e-6);
emibmk('sasia') = 1121*(1e-6);
emibmk('ssa') = 1121*(1e-6);
emibmk('te') = 914*(1e-6);
emibmk('usa') = 537*(1e-6);

emiveg('cajaz') = 275*(1e-6);
emiveg('china') = 384*(1e-6);
emiveg('easia') = 384*(1e-6);
emiveg('kosau') = 384*(1e-6);
emiveg('laca') = 384*(1e-6);
emiveg('india') = 519*(1e-6);
emiveg('mena') = 384*(1e-6);
emiveg('neweuro') = 384*(1e-6);
emiveg('oldeuro') = 275*(1e-6);
emiveg('sasia') = 519*(1e-6);
emiveg('ssa') = 519*(1e-6);
emiveg('te') = 384*(1e-6);
emiveg('usa') = 275*(1e-6);


l_curve('cajaz') =0.5;
l_curve('china') = 3.3;
l_curve('easia') = 3.302;
l_curve('kosau') = 1.103;
l_curve('laca') = 3.296;
l_curve('india') = 3.609;
l_curve('mena') = 3.302;
l_curve('neweuro') = 1.105;
l_curve('oldeuro') = 0.502;
l_curve('sasia') = 3.698;
l_curve('ssa') = 3.702;
l_curve('te') = 3.307;
l_curve('usa') = 0.503;



##  COMPUTE DATA
#_________________________________________________________________________
$elseif.ph %phase%=='compute_data'

CtoCO2 = 44 / 12 ;
CO2toC = 12 / 44 ;

* Baseline emissions
sigma(t,n) = emi_bau_co2(t,n)/ykali(t,n) ;
sig0(n)    = sigma('1',n)    ;

* Initial emissions
e0(n) = q0(n) * sig0(n);


##  DECLARE VARIABLES
#_________________________________________________________________________
$elseif.ph %phase%=='declare_vars'

## EMISSIONS CO2 ----------
VARIABLES

    E(t,n)          'Total CO2 emissions  [GtCO2/year]'
    EIND(t,n)       'Industrial emissions [GtCO2/year]'
    EDIET(t,n)      'Diet emisssions      [GtCO2/year]'

    NVEG(t,n)       '% of vegetarians in a country at a given t'

    CCAEIND(t)      'Cumulative Industrial Carbon Emissions [GTC]'
    CCAEDIET(t)     'Cumulative Diet Carbon Emissions [GTC]'
    CCAETOT(t)      'Cumulative Carbon Emissions [GTC]'
    CUMETREE(t)     'Cumulative from land [GtC]'

    CCO2EIND(t)     'Cumulative Industrial CO2 Emissions [GtCO2]'
    CCO2EDIET(t)     'Cumulative Diet CO2 Emissions [GtCO2]'
    CCO2ETOT(t)     'Cumulative CO2 Emissions [GtCO2]'

    MIU(t,n)        'Emission control rate GHGs'
    ABATEDEMI(t,n)  'Abated Emissions [GtCO2/year]'
;
POSITIVE VARIABLES  ABATEDEMI, MIU, NVEG, EDIET ;

# VARIABLES STARTING LEVELS 
* to help convergence if no startboost is loaded

     EIND.l(t,n) = sigma(t,n)*ykali(t,n);
     EDIET.l(t,n) = 0.05;
     E.l(t,n) = EIND.l(t,n) + EDIET.l(t,n) + eland_bau('uniform',t,n) ;
     MIU.l(t,n) = 0 ; 
     ABATEDEMI.l(t,n) = 0 ;

* upper and lower bounds NVEG 

     NVEG.lo(t,n) = 0;
     NVEG.up(t,n) = 1;

## EMISSIONS OGHG ----------
$ifthen.oghg %climate% == 'witchoghg'
VARIABLES
    EOGHG(oghg,t,n)       'Total OGHG emissions [GtCO2eq/year]'
    COGHGE(oghg,t)        'Cumulative OGHG emissions [GtCO2eq]'
    MIU_OGHG(oghg,t,n)    'Emission control rate GHGs'
    ABATEDOGHG(oghg,t,n)  'Abated OGHG Emissions [GtCO2eq/year]'
;
POSITIVE VARIABLES  ABATEDOGHG, MIU_OGHG;

# VARIABLES STARTING LEVELS 
* to help convergence if no startboost is loaded
     EOGHG.l(oghg,t,n) = oghg_emi_bau(oghg,t,n) ;
ABATEDOGHG.l(oghg,t,n) = 0 ;
  MIU_OGHG.l(oghg,t,n) = 0 ;
$endif.oghg

##  COMPUTE VARIABLES
#_________________________________________________________________________
$elseif.ph %phase%=='compute_vars'

##  EMISSIONS CO2 ----------
* Industrial emissions starting point
EIND.fx(tfirst,n) = e0(n) ;
* Initial vegetgarians
NVEG.fx(tfirst,n) = veg0(n);
*initial diet emission
EDIET.fx(tfirst,n) = 0;
* Resource limit for fossils emissions
CCAEIND.up(t) = fosslim ;
* Initial cumulated conditions
CCAEIND.FX(tfirst)  = cumeind0    ;
CCAEDIET.FX(tfirst)  = 0           ;
CUMETREE.fx(tfirst) = cumetree0   ;
CCAETOT.FX(tfirst)  = cumeind0 + cumetree0 ;
* CO2-budget starts empty
CCO2EIND.FX(tfirst) = 0 ;
CCO2EDIET.FX(tfirst) = 0 ;
CCO2ETOT.FX(tfirst) = 0 ;

##  OGHG EMISSION VARIABLES ----------
$ifthen.oghg %climate% == 'witchoghg'
* Oghg emissions starting point
EOGHG.FX(oghg,tfirst,n)   =  oghg_emi_bau(oghg,tfirst,n);
* OGHG-budget starts empty
COGHGE.FX(oghg,tfirst)    = 0 ;
$endif.oghg

##  CO2 MITIGATION UPPER BOUND SHAPE ----------
loop(t,
# Before transition
MIU.up(t,n)$(t.val lt t_min_miu) = min_miu;
# Transition to negative: linear transition from min_miu to max_miu between t_min_miu and t_max_miu
MIU.up(t,n)$(t.val ge t_min_miu) = min_miu + (max_miu - min_miu) * (t.val - t_min_miu)/(t_max_miu - t_min_miu);
# After transition
MIU.up(t,n)$(t.val gt t_max_miu) = max_miu;
);
MIU.up(t,n)$(t.pos le 2) = 0.03 ; #setting 2020 upperbound as in DICE2016

##  OGHG MITIGATION UPPER BOUND SHAPE ----------
$ifthen.oghg %climate% == 'witchoghg'
loop(t,
# before transition
MIU_OGHG.up(oghg,t,n)$(t.val lt t_min_miu) = min_miuoghg;
# transition: linear transition from min_miu to max_miu between t_min_miu and t_max_miu
MIU_OGHG.up(oghg,t,n)$(t.val ge t_min_miu) = min_miuoghg + (max_miuoghg - min_miuoghg) * (t.val - t_min_miu)/(t_max_miu - t_min_miu);
# after transition
MIU_OGHG.up(oghg,t,n)$(t.val gt t_max_miu) = max_miuoghg;
);
MIU_OGHG.fx(oghg,tfirst,n) = 0 ; #setting 2015 value
$endif.oghg


#=========================================================================
*   ///////////////////////     OPTIMIZATION    ///////////////////////
#=========================================================================

##  EQUATION LIST
#_________________________________________________________________________
$elseif.ph %phase%=='eql'

##  CO2 EMISSION EQUATIONS ----------
    eq_e                # Emissions equation'
    eq_eind             # Industrial emissions equation'
    eq_ediet            # Diet emissions equation'
    eq_ccaeind          # Cumulative industrial carbon emissions equation'
    eq_ccaediet         # Cumulative diet carbon emissions equation'       
    eq_ccaetot
    eq_cco2eind         # Cumulative industrial CO2 emissions equation (from now on)'
    eq_cco2diet         # Cumulative diet CO2 emissions equation (from now on)'
    eq_cco2etot
    eq_cumetree         # Cumulated land-use emissions
    eq_abatedemi        # Abated Emissions according to decision'
    eq_miuinertiaplus   # Inertia in CO2 Control Rate decreasing'
    eq_miuinertiaminus  # Inertia in CO2 Control Rate increasing'
    eq_nveginertiaminus # Inertia in CO2 Control Rate decreasing'
    eq_nveginertiaplus  # Inertia in CO2 Control Rate increasing'
    eq_veg0inertiaminus # Inertia in CO2 Control Rate increasing'
    eq_veg0inertiaplus  # Inertia in CO2 Control Rate increasing'

##  OGHG EMISSION EQUATIONS ----------
$ifthen.oghg %climate% == 'witchoghg'
    eq_eoghg                 #'OGHG emissions equation'
    eq_coghge                #'Cumulative OGHG emissions equation'
    eq_abatedoghg            #'Abated OGHG Emissions according to MIU decision'
    eq_oghg_miuinertiaplus   #'Inertia in OGHG Control Rate decreasing'
    eq_oghg_miuinertiaminus  #'Inertia in OGHG Control Rate increasing'
$endif.oghg


##  EQUATIONS
#_________________________________________________________________________
$elseif.ph %phase%=='eqs'

##  EMISSIONS CO2 ----------
* Industrial emissions
 eq_eind(t,n)$(reg(n))..   EIND(t,n)  =E=  sigma(t,n) * YGROSS(t,n) * (1-(MIU(t,n)))  ;

* Diet emissions
 eq_ediet(t,n)$(reg(n) and not tfirst(t))..  EDIET(t,n) =E=  pop(t,n)*(NVEG(t,n)*emiveg(n)+(1-NVEG(t,n))*emibmk(n));

 

* All emissions
 eq_e(t,n)$(reg(n))..   E(t,n)  =E=  EIND(t,n) + EDIET(t,n) + ELAND(t,n)  ;

* Industrial cumulated emissions in Carbon
 eq_ccaeind(t+1)..   CCAEIND(t+1)  =E=  CCAEIND(t) # All industrial emi per period in Carbon
                                   +  (( sum(n$reg(n), EIND(t,n)) + sum(n$(not reg(n)), EIND.l(t,n)) ) * tstep * CO2toC ) ; #Carbon

* Diet cumulated emissions in Carbon  
 eq_ccaediet(t+1)..  CCAEDIET(t+1) =E=  CCAEDIET(t) # All diet emi per period in Carbon    
                                   +  (( sum(n$reg(n), EDIET(t,n)) + sum(n$(not reg(n)), EDIET.l(t,n)) ) * tstep * CO2toC ) ; #Carbon                          

* Total cumulated emissions in Carbon
 eq_ccaetot(t+1)..   CCAETOT(t+1)  =E=  CCAETOT(t) # All emi (industrial + diet + land) per period in Carbon
                                   +  (( sum(n$reg(n), E(t,n)) + sum(n$(not reg(n)), E.l(t,n)) ) * tstep * CO2toC )  ; #Carbon

* Land Use cumulated emissions in Carbon
 eq_cumetree(t+1)..   CUMETREE(t+1) =E=  CUMETREE(t) # Land Use emi per period in Carbon
                                    + (( sum(n$reg(n), ELAND(t,n))    + sum(n$(not reg(n)), ELAND.l(t,n)) ) * tstep * CO2toC)  ; #Carbon

* Industrial cumulated emissions in CO2
 eq_cco2eind(t+1)..   CCO2EIND(t+1)   =E=  CCO2EIND(t) # All industrial emi per period
                                      +   (( sum(n$reg(n), EIND(t,n))    + sum(n$(not reg(n)), EIND.l(t,n)) )  * tstep )  ; #CO2

* Diet cumulated emissions in Carbon  
 eq_cco2diet(t+1)..  CCO2EDIET(t+1) =E=  CCO2EDIET(t) # All diet emi per period     
                                   +  (( sum(n$reg(n), EDIET(t,n)) + sum(n$(not reg(n)), EDIET.l(t,n)) ) * tstep) ; #CO2

* Total cumulated emissions in CO2
 eq_cco2etot(t+1)..   CCO2ETOT(t+1)   =E=  CCO2ETOT(t) # All emi (industrial + land) per period
                                      +   (( sum(n$reg(n), E(t,n))    + sum(n$(not reg(n)), E.l(t,n)) )    * tstep )  ; #CO2

* Emissions abated
 eq_abatedemi(t,n)$(reg(n))..   ABATEDEMI(t,n)  =E=  MIU(t,n) * sigma(t,n) * YGROSS(t,n)  ;


##  EMISSIONS OGHG ----------
$ifthen.oghg %climate% == 'witchoghg'
* OGHG emissions
eq_eoghg(oghg,t,n)$(reg(n))..   EOGHG(oghg,t,n)  =E=  oghg_emi_bau(oghg,t,n)  * (1-(MIU_OGHG(oghg,t,n))) ;

* OGHG cumulated emissions in CO2eq
eq_coghge(oghg,t+1)..   COGHGE(oghg,t+1)  =E=  COGHGE(oghg,t) # All oghg-emi per period
                                           + (( sum(n$reg(n), EOGHG(oghg,t,n)) + sum(n$(not reg(n)), EOGHG.l(oghg,t,n)) )* tstep) ; #CO2eq
* OGHG Emissions abated
eq_abatedoghg(oghg,t,n)$(reg(n))..   ABATEDOGHG(oghg,t,n)  =E=  MIU_OGHG(oghg,t,n) * oghg_emi_bau(oghg,t,n) ;
$endif.oghg


##  MITIGATION CO2 INERTIA ----------
* Inertia in increasing
eq_miuinertiaminus(t+1,n)$(map_nt and (t.val gt 1))..   MIU(t+1,n)  =L=  MIU(t,n) + 0.20 ;

eq_veg0inertiaminus(t,n)$(t.val = 2).. NVEG(t,n) =L= veg0(n) + 0.1;

eq_nveginertiaminus(t+1,n)$(map_nt and (t.val gt 1))..   NVEG(t+1,n) =L= NVEG(t,n) + 0.1;

* Inertia in decreasing
eq_miuinertiaplus(t+1,n)$(map_nt  and (t.val gt 1))..   MIU(t+1,n)  =G=  MIU(t,n) - 0.20 ;

eq_veg0inertiaplus(t,n)$(t.val = 2).. NVEG(t,n) =G= veg0(n) - 0.1;

eq_nveginertiaplus(t+1,n)$(map_nt and (t.val gt 1))..   NVEG(t+1,n) =G= NVEG(t,n) - 0.1;

##  MITIGATION OGHG INERTIA ----------
$ifthen.oghg %climate% == 'witchoghg'
* Inertia in increasing
eq_oghg_miuinertiaminus(oghg,t,n)$(reg(n)  and (t.val gt 1))..   MIU_OGHG(oghg,t,n)  =L=  MIU_OGHG(oghg,t-1,n) + 0.20  ;

* Inertia in decreasing
eq_oghgh_miuinertiaplus(oghg,t,n)$(reg(n)  and (t.val gt 1))..   MIU_OGHG(oghg,t,n)  =G=  MIU_OGHG(oghg,t-1,n) - 0.50  ;
$endif.oghg


##  FIX VARIABLES
#_________________________________________________________________________
$elseif.ph %phase%=='fix_variables'

tfixvar(MIU,'(t,n)')

##  BEFORE SOLVE
#_________________________________________________________________________
$elseif.ph %phase%=='before_solve'

loop((t,n)$(ctax(t,n) and (MIU.lo(t,n) lt MIU.up(t,n))),
    MIU.lo(t,n) = max(0, min(1 - (0.90*(YNET.l(t,n)-ABATECOST.l(t,n))/ctax(t,n) + EIND.l(t,n))/(sigma(t,n)*YGROSS.l(t,n)), MIU.up(t,n)));
    MIU.l(t,n) = (MIU.lo(t,n) + MIU.up(t,n))/2;
    EIND.l(t,n) = sigma(t,n)*YGROSS.l(t,n)*(1-MIU.l(t,n));
    ABATEDEMI.l(t,n) = sigma(t,n)*YGROSS.l(t,n)*MIU.l(t,n);
);


##  PROBLEMATIC REGIONS
#_________________________________________________________________________
$elseif.ph %phase%=='problematic_regions'


##  AFTER SOLVE
#_________________________________________________________________________
$elseif.ph %phase%=='after_solve'

world_e(t)            = sum(n$(nsolve(n)),  E.l(t,n)  );
$if %climate% == 'witchoghg' world_eoghg(oghg,t)   = sum(n$( nsolve(n)), EOGHG.l(oghg,t,n));

display MIU.l, CCAEIND.l;

#=========================================================================
*   ///////////////////////     SIMULATION    ///////////////////////
#=========================================================================

##  SIMULATION SETUP
#_________________________________________________________________________
$elseif.ph %phase%=='set_simulation'



##  SIMULATION HALFLOOP 1
#_________________________________________________________________________
$elseif.ph %phase%=='simulate_1'

#  CO2 EMISSIONS ---------------------------------------

* Industrial emissions
EIND.l(t,n)  =  sigma(t,n) * YGROSS.l(t,n) * (1-(MIU.l(t,n))) ; #CO2

* Diet emissions
EDIET.l(t,n) =  pop(t,n)*(NVEG.l(t,n)*emiveg(n)+(1-NVEG.l(t,n))*emibmk(n)); #CO2

* All emissions
E.l(t,n)  =  EIND.l(t,n) + EDIET.l(t,n) + ELAND.l(t,n) ; #CO2

* Industrial cumulated emissions in Carbon
CCAEIND.l(t+1)  =  CCAEIND.l(t)  + ( sum(n, EIND.l(t,n)) * tstep * CO2toC ) ; #Carbon

* Diet cumulated emissions in Carbon
CCAEDIET.l(t+1)  =  CCAEDIET.l(t)  + ( sum(n, EDIET.l(t,n)) * tstep * CO2toC ) ; #Carbon

* Total cumulated emissions in Carbon
CCAETOT.l(t+1)  =  CCAETOT.l(t)  + ( sum(n, E.l(t,n)) * tstep * CO2toC ) ; #Carbon

* Land Use cumulated emissions in Carbon
CUMETREE.l(t+1) =  CUMETREE.l(t) + (sum(n, ELAND.l(t,n)) * tstep * CO2toC) ; #Carbon

* Industrial cumulated emissions in CO2
CCO2EIND.l(t+1)  =  CCO2EIND.l(t) + ( sum(n, EIND.l(t,n))  * tstep ) ; #CO2

* Diet cumulated emissions in CO2
CCO2EDIET.l(t+1)  =  CCO2EDIET.l(t) + ( sum(n, EDIET.l(t,n))  * tstep ) ; #CO2

* Total cumulated emissions in CO2
CCO2ETOT.l(t+1)  =  CCO2ETOT.l(t) + ( sum(n, E.l(t,n))  * tstep ) ; #CO2

* Emissions abated
ABATEDEMI.l(t,n)  =  MIU.l(t,n) * sigma(t,n) * YGROSS.l(t,n)   ; #CO2


#  OGHG EMISSIONS --------------------------------------
$ifthen.oghg %climate% == 'witchoghg'
* OGHG emissions
EOGHG.l(oghg,t,n)  =  oghg_emi_bau(oghg,t,n)  * (1-(MIU_OGHG.l(oghg,t,n))) ;

* OGHG cumulated emissions in CO2eq
COGHGE.l(oghg,t+1)  =  COGHGE.l(oghg,t) +  ( sum(n, EOGHG.l(oghg,t,n)) * tstep ) ; #CO2eq

* OGHG emissions abated
ABATEDOGHG.l(oghg,t,n)  =  MIU_OGHG.l(oghg,t,n) * oghg_emi_bau(oghg,t,n) ;
$endif.oghg


##  SIMULATION HALFLOOP 2
#_________________________________________________________________________
$elseif.ph %phase%=='simulate_2'


##  AFTER SIMULATION
#_________________________________________________________________________
$elseif.ph %phase%=='after_simulation'


#===============================================================================
*     ///////////////////////     REPORTING     ///////////////////////
#===============================================================================

##  GDX ITEMS
#_________________________________________________________________________
$elseif.ph %phase%=='gdx_items'

# Sets (excl. aliases) ---------------------------------
ghg

# Parameters -------------------------------------------
fosslim
max_miu
t_max_miu
t_min_miu
world_e
CtoCO2
CO2toC
sigma
emi_bau_co2

# Variables --------------------------------------------
E
EIND
EDIET
MIU
ABATEDEMI
CCO2EIND
CCO2EDIET
CCO2ETOT
NVEG

# Equations (only for OPT. run_mode) -------------------
$if %run_mode%=='optimization' eq_e

$endif.ph