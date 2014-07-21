#cctbx.python diffuse.py pdb=ensemble.pdb probabilities=0.3,0.5,0.2 sampling=4 resolution=2.0 prefix='cypa'
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

    data = Ensemble(args['pdb'], int(args['sampling']), args['probabilities'])

    data.get_models()

    for model in data.models:
        model.get_structure_factors(float(args['resolution']))
        model.weighted_structure_factors()
        model.get_structure_factors_squared(float(args['resolution']))
        model.weighted_structure_factors_squared()


    diffuse = Diffuse(data.models, data.symmetry)

    diffuse.calculate_map(int(args['sampling']), args['prefix'])
    diffuse.extend_symmetry(1000, args['prefix'])
    diffuse.expand_friedel()


class Ensemble:
    'Class for all ensembles'

    def __init__(self, pdb, sampling, probabilities):

        self.pdb = iotbx.pdb.input(file_name=pdb)
        self.hierarchy = self.pdb.construct_hierarchy()
        self.symmetry = self.pdb.crystal_symmetry_from_cryst1()
        self.probs = probabilities
        self.sampling = sampling
        self.get_p1_expansion()
        self.get_xray_structures()


    def expand_unit_cell(self):
        data = self.symmetry.as_py_code()
        x = data.split(",")
        y = x[0]
        z=y.split("(")

        q = x[5]
        a=float(z[2])
        b = float(x[1])
        c = float(x[2])
        A = int(x[3])
        B = int(x[4])
        C = int(q[:-1])

        #Calculate new unit cell size
        scale = int(self.sampling)
        a_new = scale*a
        b_new = scale*b
        c_new = scale*c

        #Create symmetry object with new unit cell parameters (since we're expanding into P1, we can pre-set the unit cell angles)
        new_symmetry = self.expanded_unit_cell = self.symmetry.customized_copy(unit_cell = (a_new,b_new,c_new, 90, 90, 90))



    def get_p1_expansion(self):
        self.xray_structures = self.pdb.xray_structures_simple(crystal_symmetry=self.symmetry)
        self.p1_structures_no_expansion = list()
        self.p1_pdb = list()
        self.p1_expanded = list()
        count = 0
        for structure in self.xray_structures:
            struc = structure.expand_to_p1()
            pdb_string = struc.as_pdb_file()
            file_name = 'model_%d_p1.pdb' %count
            open(file_name, 'w').write(pdb_string)
            
            self.p1_pdb.append(file_name)

            new_file = 'model_%d_p1_expanded.pdb' %count
            fin = open(file_name, 'r')
            fout = open(new_file, 'w')
            lines = fin.readlines()

            for line in lines:
                data = line.split()
                if data[0] == 'CRYST1':
                    new_a = float(data[1])*self.sampling
                    new_b = float(data[2])*self.sampling
                    new_c = float(data[3])*self.sampling
                                        
                                        #line_pre = line[:58]
                                        #print len(line_pre)
                    #line_a = line_pre.replace(data[1], str(new_a))
                    #line_b = line_a.replace(data[2], str(new_b))
                    #line_c = line_b.replace(data[3], str(new_c))
                    #fout.write(line_c)
                                        #print len(line_c)
                    c = 'CRYST1'
                    a = '90.00'
                    sg = 'P 1'
                    new_line = '%s%s%s%s%s%s%s%s' %(c.rjust(6), str(new_a).rjust(9), str(new_b).rjust(9), str(new_c).rjust(9), a.rjust(7), a.rjust(7), a.rjust(7), sg.rjust(12))
                    fout.write(new_line + '\n')

                elif data[0] == 'SCALE1' or data[0] == 'SCALE2' or data[0] == 'SCALE3':
                    continue

                else:
                    fout.write(line)

            self.p1_expanded.append(new_file)
            count += 1

    def get_xray_structures(self):
        self.final_xray_structures = list()

        for expansion in self.p1_expanded:
            pdb_inp = iotbx.pdb.input(file_name=expansion)
            st = pdb_inp.xray_structure_simple()
            self.final_xray_structures.append(st)

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
            m = Model(model, self.symmetry, self.final_xray_structures[i], weights[i])
            i += 1
            self.models.append(m)

    def diffuse_scattering(self):
        return Diffuse(self.models, self.symmetry, sampling)


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
      #Re-sizes reciprocal space lattice
      #sampling = input_dict.get('sampling',1)

        correction_factor = int(sampling)*int(sampling)*int(sampling)

        file_name=prefix+"_pre_symmetry.hkl"
        f=open(file_name,'w')
        lattice = dict()

        for hkl,intensity in array:
            h = hkl[0]
            k = hkl[1]
            l = hkl[2]
            h_new = float(h)/float(sampling)
            k_new = float(k)/float(sampling)
            l_new = float(l)/float(sampling)
            h_int = int(round(h_new+0.000000001))
            k_int = int(round(k_new+0.000000001))
            l_int = int(round(l_new+0.000000001))
            intensity_new = intensity/correction_factor
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

    def extend_symmetry(self,scale_factor,prefix):
    #This reads in an hkl map and returns a .sca map
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
        ma_anom = ma.customized_copy(anomalous_flag=False)
        ma_p1 = ma_anom.expand_to_p1()

        merge.write(file_name= prefix + '.sca', miller_array=ma_p1)

        self.p1_map = prefix + '.sca'

    def expand_friedel(self):
    #This takes in a .sca map and outputs an hkl map with the Friedel pairs expanded
        map_new = self.p1_map.rstrip('.sca')

        fin = open(self.p1_map,'r')
        fout = open(map_new + '.hkl', 'w')

        lines = fin.readlines()

        for line in lines:


            data = line.split()

            if len(data) == 5:

                h = int(data[0])
                h_new = -1*h
                k = int(data[1])
                k_new = -1*k
                l = int(data[2])
                l_new = -1*l
                i = float(data[3])

            #sig = data[4]

                x_1 = (str(h)+ ' ' + str(k) + ' ' + str(l) + ' ' + str(i) + '\n')
                x_2 = (str(h_new)+ ' ' + str(k_new) + ' ' + str(l_new) + ' ' + str(i) + '\n')

                fout.write(x_1)
                fout.write(x_2)

        fin.close()
        fout.close()

        self.full_diffuse_map = map_new + '.hkl'

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
