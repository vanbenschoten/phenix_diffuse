#cctbx.python diffuse.py pdb=ensemble.pdb probabilities=0.3,0.5,0.2 resolution=2.0 prefix='cypa'
import sys
import iotbx
from iotbx import pdb
import math
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import xray
from cctbx import miller
from iotbx import scalepack
from iotbx.scalepack import merge
from libtbx.utils import Sorry
def run(arg):
    args = get_input_dict(arg)

    data = Ensemble(args['pdb'], args['probabilities'])

    data.get_models()

    for model in data.models:
        model.get_structure_factors(float(args['resolution']))
        model.weighted_structure_factors()
        model.get_structure_factors_squared(float(args['resolution']))
        model.weighted_structure_factors_squared()


    diffuse = Diffuse(data.models, data.symmetry)

    diffuse.calculate_map(int(args['sampling']), args['prefix'])
    diffuse.as_mtz(1, args['prefix'])


class Ensemble:
    'Class for all ensembles'

    def __init__(self, pdb, probabilities):

        self.pdb = iotbx.pdb.input(file_name=pdb)
        self.hierarchy = self.pdb.construct_hierarchy()
        self.symmetry = self.pdb.crystal_symmetry_from_cryst1()
        self.probs = probabilities
        self.get_xray_structures()



    def get_xray_structures(self):
        self.xray_structures = self.pdb.xray_structures_simple()


    def get_models(self):
        self.models = []
        models = self.hierarchy.models()
        weights = []

        if self.probs != None:
            new_probs = get_probabilities(self.probs)
            for i in range(0,len(new_probs)):
                d = float(new_probs[i])
                weights.append(d)

            if len(new_probs) != len(models):
                raise Sorry("The number of models and number of given probabilities must match")

        else:
            for model in models:
                d = float(1/len(models))
                weights.append(d)

        i = 0
        for model in models:
            m = Model(model, self.symmetry, self.xray_structures[i], weights[i])
            i += 1
            self.models.append(m)

    def diffuse_scattering(self):
        return Diffuse(self.models, self.symmetry)


class Model:
    'Class for each model within an ensemble'
    def __init__(self, model, symmetry, xray, probability):
        self.model = model
        self.probability = probability
        self.symmetry = symmetry
        self.xrs = xray

    def get_structure_factors(self,resolution):
        self.resolution = resolution
        self.f = self.xrs.structure_factors(d_min=resolution).f_calc()

    
    def weighted_structure_factors(self):
        self.f_weighted = self.f*self.probability

    def get_structure_factors_squared(self,resolution):
        self.resolution = resolution
        fcalc = self.xrs.structure_factors(d_min=resolution).f_calc()
        self.f_squared = abs(fcalc).set_observation_type_xray_amplitude().f_as_f_sq()


    def weighted_structure_factors_squared(self):
        self.f_squared_weighted = self.f_squared*self.probability



class Diffuse:
    'Class for all diffuse maps produced in reciprocal space'

    def __init__(self, models, symmetry):

        self.models = models
        self.symmetry = symmetry
        #self.resolution = resolution
        #self.sampling = sampling

        return


    def calculate_map(self,sampling, prefix):
        #Reads in list of Model objects and calculates <F>**2 and <F**2>
        #<F>**2 will be calculated by adding the weighted structure factors together and squaring the sum
        #<F**2> will be calculated by adding the weighted squared structure factors together

        sum_fc = None
        sum_fc_square = None

        for model in self.models:
            if sum_fc is None:
                sum_fc = model.f_weighted
                sum_fc_square = model.f_squared_weighted

            else:
                sum_fc = sum_fc + model.f_weighted
                sum_fc_square = sum_fc_square + model.f_squared_weighted

        avg_squared = abs(sum_fc).set_observation_type_xray_amplitude().f_as_f_sq()
        _sum_fc_square, _avg_squared = sum_fc_square.common_sets(avg_squared)
        self.diffuse_signal = _sum_fc_square.customized_copy(data = _sum_fc_square.data() - _avg_squared.data())

        self.write_squared_amplitudes(sampling, prefix, self.diffuse_signal)

    def write_squared_amplitudes(self,sampling,prefix,array,out=sys.stdout):

        file_name=prefix+".hkl"
        f=open(file_name,'w')
        lattice = dict()

        for hkl,intensity in array:
            if h_int not in lattice:
                lattice[h_int] = dict()

            if k_int not in lattice[h_int]:
                lattice[h_int][k_int] = dict()

            if l_int not in lattice[h_int][k_int]:
                lattice[h_int][k_int][l_int] = 0

            lattice[h_int][k_int][l_int] += intensity_new

        for key_h in lattice:
            for key_k in lattice[key_h]:
                for key_l in lattice[key_h][key_k]:
                    print >>f, "%4d %4d %4d %4d" %(key_h, key_k, key_l, lattice[key_h][key_k][key_l])
        f.close()

        print >>out, "Wrote to %s" %(file_name)

        self.diffuse_file = file_name

    def as_mtz(self,scale_factor,prefix):
    #This reads in an hkl map and returns a .mtz map
        #Read in hkl file and populate miller array
        inf = open(self.diffuse_file, 'r')
        indices = flex.miller_index()
        i_obs = flex.double()
        sig_i = flex.double()
        for line in inf.readlines():
            assert len(line.split())==4
            line = line.strip().split()
            #####ATTENTION:SCALE FACTOR##############
            i_obs_ = float(line[3])/scale_factor #is a uniform scale factor meant to re-size all diffuse intensities (normally too large for scalepack)
            sig_i_ = math.sqrt(i_obs_) 
            #if(abs(i_obs_)>1.e-6): # perhaps you don't want zeros
            indices.append([int(line[0]),int(line[1]),int(line[2])])
            i_obs.append(i_obs_)
            sig_i.append(sig_i_)
        inf.close()

      # get miller array object
        cs = self.symmetry
        ma = miller.array(miller_set=miller.set(cs, indices), data=i_obs, sigmas=sig_i)
        ma.set_observation_type_xray_intensity()
        mtz_dataset = ma.as_mtz_dataset(column_root_label="I")
        mtz_dataset.mtz_object().write(prefix + '.mtz')

def get_input_dict(args):
    dic = dict()
    for arg in args:
        spl=arg.split('=')
        if len(spl)==2:
            dic[spl[0]] = spl[1]

    if 'probabilities' not in dic:
        dic['probabilities'] = None

    return dic

def get_probabilities(input):
    data = input.split(',')

    d_new = []

    for d in data:
        p_n = float(d)
        d_new.append(p_n)

    total = 0.0

    for i in d_new:
        total += i

    if total != 1.0:
        raise Sorry("Sorry, the given probabilities must sum to one")
 

    return d_new

if __name__ == '__main__':
    import sys
    args = sys.argv[1:]
    run(args)
