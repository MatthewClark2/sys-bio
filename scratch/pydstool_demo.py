from PyDSTool import args, Generator, perf_counter, ContClass, show
from sys import exit


def create_pcargs(name, freepars, StepSize, MaxNumPoints, LocBifPoints='all',  alg_type='EP-C'):
    PCargs = args(name=name, type=alg_type)

    PCargs.freepars = freepars
    PCargs.StepSize = StepSize
    PCargs.MaxNumPoints = MaxNumPoints
    PCargs.LocBifPoints = LocBifPoints
    PCargs.verbosity = 2

    return PCargs


class Thing():
    def __init__(self, name, model, ics):
        self.cont_classes = []
        self.name = name
        self.model = model
        self.ics = [ics]
        self.DSargs = args(name=self.name)
        self.TestDS = None

        self._revise_dsargs()

    def _revise_dsargs(self, index=-1):
        self.DSargs.ics = self.ics[-1]
        self.DSargs.pars = self.model['pars']
        self.DSargs.varspecs = self.model['vars']
        self.DSargs.fnspecs = self.model['fns']

        self.DSargs.algparams = {'refine': True}

        self.testDS = Generator.Radau_ODEsystem(self.DSargs)

    def last_cont_class(self):
        return self.cont_classes[-1] if len(self.cont_classes) > 0 else None

    def change_initial_conditions(self, new_ics=None, index=None):
        if new_ics is None and index is None:
            print('No change to be made')
        elif new_ics is not None:
            self.ics.append(new_ics)
            self._revise_dsargs()
        else:
            self._revise_dsargs(index)

    def register_continuation_class(self, PCargs):
        name = PCargs.name

        PyCont = ContClass(self.testDS)

        PyCont.newCurve(PCargs)

        # TODO(matthew-c21): Change this from printing to logging so it can be disabled or filtered later.

        # Record the computed continuation class.
        self.cont_classes.append((name, PyCont))

    def display_cont_class(self, pars, index=-1):
        (name, PyCont) = self.cont_classes[index]
        print('Computing curve...')
        start = perf_counter()
        PyCont[name].backward()
        PyCont[name].forward()
        print('done in %.3f seconds!' % (perf_counter()-start))

        PyCont[name].display(pars, figure='new')
        PyCont.plot.toggleAll('off', bytype='P')

        if PyCont[name].getSpecialPoint('H1'):
            self._limit_cycle_curve(pars, index)

    def _limit_cycle_curve(self, pars, index=-1):
        (name, PyCont) = self.cont_classes[index]

        # Limit cycle curve
        PCargs = args(freepars=PyCont[name].freepars)
        PCargs.name = 'LC1'
        PCargs.type = 'LC-C'
        PCargs.MaxNumPoints = 2500  # // 10
        PCargs.NumIntervals = 20
        PCargs.NumCollocation = 6
        PCargs.initpoint = name+':H1'
        PCargs.SolutionMeasures = 'all'
        PCargs.NumSPOut = 50
        PCargs.FuncTol = 1e-10
        PCargs.VarTol = 1e-10
        PCargs.TestTol = 1e-8
        PyCont.newCurve(PCargs)

        print('Computing limit-cycle curve...')
        start = perf_counter()
        PyCont['LC1'].forward()
        PyCont['LC1'].backward()
        print('done in %.3f seconds!' % (perf_counter()-start))

        PyCont['LC1'].display(pars, stability=True)

        PyCont['LC1'].display(stability=True, figure='new', axes=(1, 2, 1))
        PyCont['LC1'].plot_cycles(
            coords=('R', 'r'), linewidth=1, axes=(1, 2, 2), figure='fig2')
        pass


model_mmi2_full = {
    'pars': {
        'mu': 5,
        'b1': 9.0,
        'b2': 1.3 * 1,
        'a1': 0.7,
        'a2': 0.25,
        'sR': 1.0,
        'kR': 1,
        'kr': 1,
        'K1': 10000,
        'K2': 10000,
        'Ts': 10,
        'Tr': 1.0,
        'g': 0.2
    },
    'vars': {
        'r':
        'Tr * (mu - kr * (r - 1 * 2 * c1 - 2 * 1 * c2) \
            - b1 * 1 * 2 * kr * c1 \
            - b2 * 2 * 1 * kr * c2)',
        'R':
        'g*(sR - kR * (R - 2 * c1 - 1 * c2) \
            - 2 * kR * c1 * a1 - 1 * kR * c2 * a2)',
        'c1':
        'Ts * (K1 * (R - 2 * c1 - c2) * (r - 2 * c1 - 2 * c2) - c1)',
        'c2':
        'Ts * (K2 * c1 * (r - 2 * c1 - 2 * c2) - c2)',
    },

    'fns': {}, 'aux': [], 'name': 'mmi2'}

ics_1 = {'r': 0.8, 'R': 0.1, 'c1': 0.1, 'c2': 0.2}

thing = Thing('mmi2_full', model_mmi2_full, ics_1)
thing.register_continuation_class(create_pcargs('EQ1', ['mu'], 1e-4, 1000))
thing.display_cont_class(('mu', 'R'))

show()
