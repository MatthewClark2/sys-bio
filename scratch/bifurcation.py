from PyDSTool import args, Generator, perf_counter, ContClass, show
import matplotlib.pyplot as plt
from model import model_mmi2 as m, ics_1 as ics


DSargs = args(name='bifurcation')
DSargs.ics = ics
DSargs.pars = m['pars']
DSargs.varspecs = m['vars']
DSargs.fnspecs = m['fns']
DSargs.algparams = {'refine': True, 'init_step': 1e-3}


output_file = 'output.bin'


class Continuation:
    def __init__(self, dsargs):
        self.ics = dsargs.ics
        self.model_params = dsargs
        self._init_ODESystem()

    def _init_ODESystem(self):
        self.model_params.ics = self.ics
        self.testDS = Generator.Radau_ODEsystem(self.model_params)
        self.pc = ContClass(self.testDS)

    def _find_steady_state(self, tdata):
        ssargs = args(name='steadyState')
        ssargs.tdata = tdata
        # Test code.
        ssargs.ics = ics
        ssargs.pars = self.model_params.pars
        ssargs.varspecs = self.model_params.varspecs
        ssargs.fnspecs = self.model_params.fnspecs
        ssargs.algparams = {'refine': True}

        ssDS = Generator.Radau_ODEsystem(ssargs)
        trajectory = ssDS.compute('steady state')
        pts = trajectory.sample()

        return pts

    def change_to_steady_state(self, tdata=None):
        if tdata is None:
            # Reasonable default.
            tdata = [0, 20]

        pts = self._find_steady_state(tdata)

        for var in self.ics.keys():
            self.ics[var] = pts[-1][var]

        self._init_ODESystem()

    def compute(self, pcargs, msg=None):
        self.pc.newCurve(pcargs)

        if msg is not None:
            print(msg)
            start = perf_counter()
            self.pc[pcargs.name].forward()
            print('Time elapsed: %.3f' % (perf_counter() - start))
        else:
            self.pc[pcargs.name].forward()

    def limit_cycle(self, name, initpoint, freepars):
        # TODO(matthew-c21): See below comment.
        pcargs = args(name=name, type='LC-C')
        pcargs.MaxNumPoints = 2500  # // 10
        pcargs.NumIntervals = 30
        pcargs.NumCollocation = 6
        pcargs.initpoint = initpoint
        pcargs.SolutionMeasures = 'all'
        pcargs.NumSPOut = 200
        pcargs.FuncTol = 1e-4
        pcargs.VarTol = 1e-4
        pcargs.TestTol = 1e-3
        pcargs.SaveEigen = True
        pcargs.freepars = freepars

        self.compute(pcargs, 'Computing limit cycle...')

    def dual_bifurcation(self, name, initpoint, p1, p2):
        # TODO(matthew-c21): Fancy graph stuff. Possibly use smaller aggregate objects to keep track of that by itself.
        pcargs = args()
        pcargs.name = name
        pcargs.type = 'H-C2'
        pcargs.initpoint = initpoint
        pcargs.freepars = [p1, p2]
        pcargs.MaxNumPoints = 800
        pcargs.MinStepSize = 1e-5
        pcargs.MaxStepSize = 1e-2
        pcargs.StepSize = 1e-3
        pcargs.LocBifPoints = 'all'
        pcargs.SaveEigen = True

        self.compute(pcargs, 'Computing 2-parameter bifurcation analysis...')

    def display(self, pars, name=None, **kwargs):
        if name is None:
            for c in self.pc.curves:
                self.pc[c].display(pars, figure='new', **kwargs)
        else:
            self.pc[name].display(pars, figure='new', **kwargs)

        self.pc.plot.toggleAll('off', bytype='P')

    def save(self, file):
        pass

    @staticmethod
    def load(file):
        pass


cc = Continuation(DSargs)
pts = cc._find_steady_state([0, 12])

plt.plot(pts['t'], pts['R'], label='R')
plt.plot(pts['t'], pts['c1'], label='c1')
plt.plot(pts['t'], pts['c2'], label='c2')
plt.plot(pts['t'], pts['r'], label='r')

plt.legend(loc='upper right')

# Set up continuation class.
PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['mu']
PCargs.StepSize = 3e-4
PCargs.MaxStepSize = 3e-3
PCargs.MinStepSize = 3e-5
PCargs.MaxNumPoints = 500
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2

cc.compute(PCargs)

cc.limit_cycle('LC1', 'EQ1:H1', ['mu'])
# cc.dual_bifurcation('HB1', 'EQ1:H1', 'mu', 'a1')

cc.display(('mu', 'R'), stability=True)
plt.show()

'''
cc.find_bifurcation(pcargs)
cc.display()

cc.save('output.bin')
'''
