import tellurium as te
from rrplugins import Plugin
from PyDSTool import args, perf_counter, Generator, ContClass, show
from PyDSTool import plot as plt
from sys import exit


def compute_HO(initpoint='', name='HO1', ax=None, p1='a1', p2='a2'):
    PCargs.name = name
    PCargs.type = 'H-C2'
    PCargs.initpoint = initpoint
    PCargs.freepars = [p1, p2]
    PCargs.MaxNumPoints = 800
    PCargs.MinStepSize = 1e-5
    PCargs.MaxStepSize = 1e-2
    PCargs.StepSize = 1e-3
    PCargs.LocBifPoints = 'all'
    PCargs.SaveEigen = True
    PyCont.newCurve(PCargs)
    print('Computing Hopf curve...')
    start = perf_counter()
    # PyCont[name].forward()
    PyCont[name].backward()
    #PCargs.MaxNumPoints = 2000
    # PyCont[name].forward()
    PyCont[name].display([p1, p2], figure=3)

    sol = PyCont[name].sol
    a1 = sol.coordarray[sol.coordnames.index(p1)]
    b1 = sol.coordarray[sol.coordnames.index(p2)]
    if not ax:
        fig, ax = plt.subplots(figsize=(4, 4))
        fig.subplots_adjust(left=.2, bottom=0.17, right=0.98)
    ax.plot(a1, b1, label=name)
    return ax


model_mmi2 = {
    'pars': {
        'mu': 0.1,
        'b1': 1.00 * 0.8 * 1,  # 9.0,
        'b2': 6.0 * 0.8 * 0.5 * 1,  # 1.3 * 1,
        'a1': 1.00,  # 0.7,
        'a2': 6.0 * 2 * 0.8,  # 0.25,
        'sR': 1.0,
        'kR': 1,
        'kr': 1,
        'K1': 10000,
        'K2': 10000,
        'Ts': 10,
        'Tr': 0.5,
        'g': 2.0
    },
    'vars': {
        'r': \
        'Tr * (mu - kr * (r - 1 * 2 * c1 - 2 * 1 * c2) \
            - b1 * 1 * 2 * kr * c1 \
            - b2 * 2 * 1 * kr * c2)',
        'R': \
        'g*(sR - kR * (R - 2 * c1 - 1 * c2) \
            - 2 * kR * c1 * a1 - 1 * kR * c2 * a2)',
        'c1':\
        'Ts * (K1 * (R - 2 * c1 - c2) * (r - 2 * c1 - 2 * c2) - c1)',
        'c2':\
        'Ts * (K2 * c1 * (r - 2 * c1 - 2 * c2) - c2)',
    },

    'fns': {}, 'aux': [], 'name': 'mmi2'}

ics_1 = {'r': 0.8, 'R': 0.1, 'c1': 0.0, 'c2': 0.0}

DSargs = args(name='mmi2')
DSargs.pars = model_mmi2['pars']
DSargs.varspecs = model_mmi2['vars']
DSargs.fnspecs = {}
DSargs.ics = ics_1
DSargs.algparams = {'refine': True,
                    }
testDS = Generator.Radau_ODEsystem(DSargs)

# Set up continuation class
PyCont = ContClass(testDS)

PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['mu']
#PCargs.StepSize = 3e-6
#PCargs.MaxStepSize = 3e-5
#PCargs.MinStepSize = 3e-7
PCargs.StepSize = 3e-4
PCargs.MaxStepSize = 3e-3
PCargs.MinStepSize = 3e-5
#PCargs.MaxNumPoints = 8000
PCargs.MaxNumPoints = 500
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
#PCargs.StopAtPoints = ['H2']
#PCargs.StopAtPoints = ['H1']
PyCont.newCurve(PCargs)

print('Computing curve...')
start = perf_counter()
PyCont['EQ1'].forward()
# PyCont['EQ1'].backward()
print('done in %.3f seconds!' % (perf_counter()-start))

PyCont['EQ1'].display(('mu', 'R'), figure='new')
PyCont.plot.toggleAll('off', bytype='P')

if PyCont['EQ1'].getSpecialPoint('H1'):
    if 0:  # Limit cycle curve
        PCargs.name = 'LC1'
        PCargs.type = 'LC-C'
        PCargs.MaxNumPoints = 2500  # // 10
        PCargs.NumIntervals = 30
        PCargs.NumCollocation = 6
        PCargs.initpoint = 'EQ1:H1'
        PCargs.SolutionMeasures = 'all'
        PCargs.NumSPOut = 200
        PCargs.FuncTol = 1e-4
        PCargs.VarTol = 1e-4
        PCargs.TestTol = 1e-3
        PCargs.SaveEigen = True
        PyCont.newCurve(PCargs)

        print('Computing limit-cycle curve...')
        start = perf_counter()
        PyCont['LC1'].forward()
        # PyCont['LC1'].backward()
        print('done in %.3f seconds!' % (perf_counter()-start))

        # Get the lower branch of the cycle.
        PyCont['LC1'].display(('mu', 'R_min'), stability=True, figure=1)
        PyCont['LC1'].display(('mu', 'R'), stability=True, figure=1)

        PyCont['LC1'].display(stability=True, figure='new', axes=(1, 2, 1))
        PyCont['LC1'].plot_cycles(
            coords=('R', 'r'), linewidth=1, axes=(1, 2, 2), figure='fig2')

    # 2-parameter bifurcation. Limit cycle curve continuation needs to be off.
    if 1:
        ax1 = compute_HO(initpoint='EQ1:H1', p1='mu', p2='a1')

        if PyCont['EQ1'].getSpecialPoint('H2'):
            compute_HO(initpoint='EQ1:H2', ax=ax1,
                       name='HO2', p1='mu', p2='a1')
            ax1.legend()

if 0:
    fig, ax = plt.subplots(figsize=(4, 4))
    fig.subplots_adjust(left=.2, bottom=0.17, right=0.98)
    sol = PyCont['HO1'].sol
    a1 = sol.coordarray[sol.coordnames.index('a1')]
    b1 = sol.coordarray[sol.coordnames.index('b1')]
    ax.plot(a1, b1, label='HO1')
    sol = PyCont['HO2'].sol
    a1 = sol.coordarray[sol.coordnames.index('a1')]
    b1 = sol.coordarray[sol.coordnames.index('b1')]
    ax.plot(a1, b1, label='HO2')
show()

exit()


#################################################################################
#   TELLURIUM
#################################################################################


model_str = '// Reactions\n\tconst syn_const\n\t'

for i, var in enumerate(sorted(model_mmi2['vars'], reverse=False)):
    # print(var)
    de = model_mmi2['vars'][var]
    model_str += 'syn_const -> ' + var + '; ' + de + '\n\t'

model_str += '\n// Species Init\n\t'

for k, v in ics_1.items():
    model_str += k + ' = ' + str(round(v, 4)) + '; '

model_str += '\n\n// Parameters\n\t'

for k, v in model_mmi2['pars'].items():
    model_str += k + ' = ' + str(v) + '; '

r = te.loada(model_str)


#r['a1'], r['b1'] = 0.77, 0.6
#r['a1'], r['b1'] = 0.13, 0.28
#r['a1'], r['b1'] = 1.4, 0.8
#r['a1'], r['b1'] = 0.54, 0.16
#r['a1'], r['b1'] = 1.4, 2.3
#r['a1'], r['b1'] = 1.17, 0.2
#r['a1'], r['b1'] = 1.4, 1.17
#r['mu'] = 1.55
r['R'] = 0.4
r['r'] = 0.3
r['c1'] = 0.1
r['c2'] = 0.05
r['mu'] = 0.42


r['mu'] = 0.5
#r['R'] = 0.15
r['r'] = 0.17
m = r.simulate(0, 100, 10000)
r.plot()

# exit()


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
    auto.setProperty("RL1", 2.5)
    auto.setProperty("NMX", 10000)
    #auto.setProperty("NDIM", 15)
    auto.setProperty("NPR", 2)
    auto.setProperty("KeepTempFiles", True)
    auto.setProperty("DS", 0.001)
    auto.setProperty("DSMIN", 0.0001)
    auto.setProperty("DSMAX", 0.1)

    auto.execute()
    pts1 = auto.BifurcationPoints

    if 1:
        lbl1 = auto.BifurcationLabels
        biData1 = auto.BifurcationData
        biData1.plotBifurcationDiagram(pts1, lbl1)


auto1 = run_auto(pars=None, r=r)
