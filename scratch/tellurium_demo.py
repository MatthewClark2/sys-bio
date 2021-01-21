import tellurium as te
import roadrunner
from rrplugins import Plugin

model_str = '// Reactions\n\tconst syn_const\n\t'

for i, var in enumerate(sorted(model_mmi2_full['vars'], reverse=False)):
    # print(var)
    de = model_mmi2_full['vars'][var]
    model_str += 'syn_const -> ' + var + '; ' + de + '\n\t'

model_str += '\n// Species Init\n\t'

for k, v in ics_1.items():
    model_str += k + ' = ' + str(round(v, 4)) + '; '

model_str += '\n\n// Parameters\n\t'

for k, v in model_mmi2_full['pars'].items():
    model_str += k + ' = ' + str(v) + '; '

r = te.loada(model_str)


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
    auto.setProperty("RL1", 5)
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
