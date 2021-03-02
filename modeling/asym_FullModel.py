import tellurium as te
import roadrunner
from rrplugins import Plugin
from PyDSTool import *
from sys import exit

model_mmi2_full = {
    'pars':{
        'mu'    : 0.1,
        'b1P'   : 1.0,
        'b1Q'   : 6.0,
        'b2P'   : 1.0,
        'b2Q'   : 3.0,
        'a1P'   : 1.0,
        'a1Q'   : 1.0,
        'a2'    : 4.0,
        'sR'    : 1.0,
        'kR'    : 1.0,
        'kr'    : 1.0,
        'K1'    : 1000,
        'K2'    : 1000,
        'Ts'    : 10,
        'koff'  : 100,
    },
    'vars':{
        'r' : \
            'mu - kr * r + koff * (c1A+c1B) - K1 * koff * 2*R * r  \
            + koff * 2 * c2 - K2 * koff * (c1A+c1B) * r \
            + c1A * kR * a1P +c1B * kR * a1Q+ c2 * kR * a2 * 2',
        'R' : \
            'sR - kR * R + koff * (c1A+c1B) - K1 * koff * 2*R * r  \
            + (c1A) * kr * b1P + c1B * kr * b1Q',
        'c1A':\
            'K1 * koff * R * r - koff * c1A + koff * c2 - K2 * koff * c1A * r \
            + c2 * 1 * kr * b2P - c1A * kR * a1P - c1A * kr * b1P',
        'c1B':\
            'K1 * koff * R * r - koff * c1B + koff * c2 - K2 * koff * c1B * r \
            + c2 * 1 * kr * b2Q - c1B * kR * a1Q - c1B * kr * b1Q',
        'c2':\
            'K2 * koff * (c1A+c1B) * r - koff * 2 * c2 \
            - c2 * kr * (b2P + b2Q) - c2 * kR * a2',
    },

'fns': {}, 'aux': [], 'name':'mmi2_full'}

ics_1 = {'r': 0.8, 'R': 0.1, 'c1A': 0.0, 'c1B': 0.0, 'c2': 0.0}

DSargs = args(name='mmi2_full_asym')
DSargs.pars = model_mmi2_full['pars']
DSargs.varspecs = model_mmi2_full['vars']
DSargs.fnspecs = {}
DSargs.ics = ics_1
DSargs.algparams = {'refine': True}
testDS = Generator.Radau_ODEsystem(DSargs)

# Set up continuation class
PyCont = ContClass(testDS)

PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['mu']
PCargs.StepSize = 1e-2
PCargs.MaxNumPoints = 300
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PyCont.newCurve(PCargs)

print('Computing curve...')
start = perf_counter()
PyCont['EQ1'].backward()
#PyCont['EQ1'].forward()
print('done in %.3f seconds!' % (perf_counter()-start))

PyCont['EQ1'].display(('mu','R'), figure='new')
PyCont.plot.toggleAll('off', bytype='P')

if PyCont['EQ1'].getSpecialPoint('H2'): # H1 or H2
    # Limit cycle curve
    PCargs.name = 'LC1'
    PCargs.type = 'LC-C'
    PCargs.MaxNumPoints = 2500 #// 10
    PCargs.NumIntervals = 20
    PCargs.NumCollocation = 6
    PCargs.initpoint = 'EQ1:H2'
    PCargs.SolutionMeasures = 'all'
    PCargs.NumSPOut = 50
    PCargs.FuncTol = 1e-10
    PCargs.VarTol = 1e-10
    PCargs.TestTol = 1e-8
    PyCont.newCurve(PCargs)

    print('Computing limit-cycle curve...')
    start = perf_counter()
    PyCont['LC1'].forward()
    #PyCont['LC1'].backward()
    print('done in %.3f seconds!' % (perf_counter()-start))

    PyCont['LC1'].display(('mu','R'), stability=True)

    PyCont['LC1'].display(stability=True, figure='new', axes=(1,2,1))
    PyCont['LC1'].plot_cycles(coords=('R','r'), linewidth=1, axes=(1,2,2), figure='fig2')

show()

#exit()

#################################################################################
#   TELLURIUM
#################################################################################

model_str = '// Reactions\n\tconst syn_const\n\t'

for i, var in enumerate(sorted(model_mmi2_full['vars'], reverse=False)):
    #print(var)
    de = model_mmi2_full['vars'][var]
    model_str += 'syn_const -> ' + var + '; ' + de + '\n\t'

model_str += '\n// Species Init\n\t'

for k, v in ics_1.items():
    model_str += k + ' = ' + str(round(v,4)) + '; '

model_str += '\n\n// Parameters\n\t'

for k, v in model_mmi2_full['pars'].items():
    model_str += k + ' = ' + str(v) + '; '

r = te.loada(model_str)


m = r.simulate(0, 100, 10000)
r.plot()
#exit()


def run_auto(pars, r, direction='Positive'):

    if pars:
        for k in pars:
            r[k] = pars[k]

    auto = Plugin("tel_auto2000")
    # Setup properties
    auto.setProperty("SBML", r.getCurrentSBML())
    auto.setProperty("ScanDirection", direction)
    auto.setProperty("PrincipalContinuationParameter", "mu")
    auto.setProperty("PreSimulation", "True")
    auto.setProperty("PreSimulationDuration", 30)
    auto.setProperty("RL0", 0)
    auto.setProperty("RL1", 5)
    auto.setProperty("NMX", 10000) 
    #auto.setProperty("NDIM", 15) 
    auto.setProperty("NPR", 2) 
    auto.setProperty("KeepTempFiles", True) 
    auto.setProperty("DS", 0.001) 
    auto.setProperty("DSMIN", 0.0001) 
    auto.setProperty("DSMAX", 0.1) 

    auto.execute()
    pts1     = auto.BifurcationPoints

    if 1:
        lbl1     = auto.BifurcationLabels
        biData1  = auto.BifurcationData
        biData1.plotBifurcationDiagram(pts1, lbl1)

auto1 = run_auto(pars=None,r=r)

