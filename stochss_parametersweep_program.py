import gillespy
import random
import uuid
import sys
import itertools
import json
import numpy
import os
import pickle
from subprocess import Popen, PIPE
import molns_cloudpickle as CloudPickle
import shutil



# NOTE the string '_t__JSON_STRING___' (without the extra 't' at the start) will be replaced by data from the model


class StochSSModel(gillespy.Model):
    json_data = json.loads("""{
    "increment": 1, 
    "isLocal": true, 
    "isSpatial": false, 
    "logA": false, 
    "logB": false, 
    "maxTime": 100, 
    "maxValueA": 2.9700000000000006, 
    "maxValueB": 0.036300000000000006, 
    "minValueA": 2.43, 
    "minValueB": 0.0297, 
    "modelType": "massaction", 
    "name": "birthdeath", 
    "parameterA": "k1", 
    "parameterB": "k2", 
    "parameters": [
        {
            "name": "k1",
            "value": "2.7"
        }, 
        {
            "name": "k2", 
            "value": "0.033"
        }
    ], 
    "reactions": [
        {
            "name": "R1", 
            "products": [
                {
                    "specie": "S", 
                    "stoichiometry": 1
                }
            ], 
            "rate": "k1", 
            "reactants": [], 
            "type": "creation"
        }, 
        {
            "name": "R2", 
            "products": [], 
            "rate": "k2", 
            "reactants": [
                {
                    "specie": "S", 
                    "stoichiometry": 1
                }
            ], 
            "type": "destruction"
        }
    ], 
    "seed": -1, 
    "species": [
        {
            "initialCondition": 100, 
            "name": "S"
        }
    ], 
    "speciesSelect": {
        "S": true
    }, 
    "stepsA": 10, 
    "stepsB": 10, 
    "trajectories": 1, 
    "variableCount": 1
}""")

    def __init__(self, **kwargs):
        modelType = self.json_data["modelType"]
        species = self.json_data["species"]
        parameters = self.json_data["parameters"]
        reactions = self.json_data["reactions"]
        maxTime = self.json_data["maxTime"]
        if maxTime is None:
            maxTime = 100
        increment = self.json_data["increment"]
        if increment is None:
            increment = 1

        gillespy.Model.__init__(self, name = self.json_data["name"])

        parameterByName = dict()

        for parameter in parameters:
            if parameter['name'] in kwargs:
                parameterByName[parameter['name']] = gillespy.Parameter(name = parameter['name'], expression = kwargs[parameter['name']])
            else:
                parameterByName[parameter['name']] = gillespy.Parameter(name = parameter['name'], expression = parameter['value'])

            self.add_parameter(parameterByName[parameter['name']])

        speciesByName = dict()

        for specie in species:
            speciesByName[specie['name']] = gillespy.Species(name = specie['name'], initial_value = specie['initialCondition'])

            self.add_species(speciesByName[specie['name']])

        for reaction in reactions:
            inReactants = dict()
            for reactant in reaction['reactants']:
                if reactant['specie'] not in inReactants:
                    inReactants[reactant['specie']] = 0

                inReactants[reactant['specie']] += reactant['stoichiometry']

            inProducts = dict()
            for product in reaction['products']:
                if product['specie'] not in inProducts:
                    inProducts[product['specie']] = 0

                inProducts[product['specie']] += product['stoichiometry']

            reactants = dict([(speciesByName[reactant[0]], reactant[1]) for reactant in inReactants.items()])

            products = dict([(speciesByName[product[0]], product[1]) for product in inProducts.items()])
            
            if(reaction['type'] == 'custom'):
                self.add_reaction(gillespy.Reaction(name = reaction['name'], reactants = reactants, products = products, propensity_function = reaction['equation']))
            else:
                self.add_reaction(gillespy.Reaction(name = reaction['name'], reactants = reactants, products = products, rate = parameterByName[reaction['rate']]))

        self.timespan(numpy.concatenate((numpy.arange(maxTime / increment) * increment, [maxTime])))


#model = Model()
#output = model.run()
#print output


def mapAnalysis(result):
    metrics = { 'max' : {}, 'min' : {}, 'avg' : {}, 'var' : {}, 'finalTime' : {} }
    for i, specie in enumerate(statsSpecies):
        metrics['max'][specie] = numpy.max(result[specie])
        metrics['min'][specie] = numpy.min(result[specie])
        metrics['avg'][specie] = numpy.mean(result[specie])
        metrics['var'][specie] = numpy.var(result[specie])
        metrics['finalTime'][specie] = result[specie][-1]

    return metrics


def reduceAnalysis(metricsList):
    reduced = {}

    keys1 = ['max', 'min', 'avg', 'var', 'finalTime']
    for key1, key2 in itertools.product(keys1, statsSpecies):
        toReduce = [metrics[key1][key2] for metrics in metricsList]

        if key1 not in reduced:
            reduced[key1] = {}

        reduced[key1][key2] = {
            'max' : numpy.max(toReduce),
            'min' : numpy.min(toReduce),
            'avg' : numpy.mean(toReduce),
            'var' : numpy.var(toReduce)
        }
        
    return reduced


def getParameters():
    parameters = dict()
    if StochSSModel.json_data['logA']:
        parameters[StochSSModel.json_data['parameterA']] = numpy.logspace(numpy.log10(StochSSModel.json_data['minValueA']), numpy.log10(StochSSModel.json_data['maxValueA']), StochSSModel.json_data['stepsA'])
    else:    
        parameters[StochSSModel.json_data['parameterA']] = numpy.linspace(StochSSModel.json_data['minValueA'], StochSSModel.json_data['maxValueA'], StochSSModel.json_data['stepsA'])
        
    if StochSSModel.json_data['variableCount'] != 1:
        if StochSSModel.json_data['logB']:
            parameters[StochSSModel.json_data['parameterB']] = numpy.logspace(numpy.log10(StochSSModel.json_data['minValueB']), numpy.log10(StochSSModel.json_data['maxValueB']), StochSSModel.json_data['stepsB'])
        else:
            parameters[StochSSModel.json_data['parameterB']] = numpy.linspace(StochSSModel.json_data['minValueB'], StochSSModel.json_data['maxValueB'], StochSSModel.json_data['stepsB'])
    return parameters;       


def run_local_parameter_sweep(parameters, mapper_fn=mapAnalysis, reducer_fn=reduceAnalysis):

    pset_list = []
    dat = []
    if len(parameters) > 1:
        pnames = parameters.keys()
        for pvals1 in parameters[pnames[0]]:
            for pvals2 in parameters[pnames[1]]:
                pset = {pnames[0]:pvals1, pnames[1]:pvals2}
                pset_list.append(pset)
    else:
        for pname,pvals in parameters.iteritems():
            for pval in pvals:
                pset = {pname:pval}
                pset_list.append(pset)

    if 'seed' in StochSSModel.json_data and int(StochSSModel.json_data['seed']) != -1:
        seed = int(StochSSModel.json_data['seed'])
    else:
        random.seed()
        seed = random.randint(0, 2147483647)

    for pndx,pset in enumerate(pset_list):
        model = StochSSModel(**pset)
        results = model.run(number_of_trajectories = StochSSModel.json_data['trajectories'], seed=seed+pndx, show_labels = True)
        if not isinstance(results, list):
            results = [results]
        mapped_list = []
        for r in results:
            mapped_list.append(mapper_fn(r))
        dat.append({ 'parameters' : pset, 'result' : reducer_fn(mapped_list) })

    return dat


def get_unpickled_result(base_dir):
    with open(os.path.join(base_dir, "output"), "rb") as output:
        return pickle.load(output)


def wait_for_all_results_to_return(dirs):

    import time
    print "waiting for all results to be computed.."

    while len(dirs) > 0:
        for dir in dirs:
            output_file = os.path.join(dir, "output")
            if os.path.exists(output_file):
                dirs.remove(dir)
                print "{0} exists".format(output_file)
            else:
                print "{0} does not exist".format(output_file)
        time.sleep(1)
    print "all results computed."


def clean_up(dirs_to_delete, containers_to_delete):

    print "removing temp directory"
    for dir in dirs_to_delete:
        shutil.rmtree(dir)

    print "removing finished containers.."
    for container in containers_to_delete:
        Popen(['sudo', 'docker', 'rm', '-f', container], shell=False)


def run_qsub_parameter_sweep(parameters, mapper_fn=mapAnalysis, reducer_fn=reduceAnalysis):
    pset_list = []
    dat = []
    if len(parameters) > 1:
        pnames = parameters.keys()
        for pvals1 in parameters[pnames[0]]:
            for pvals2 in parameters[pnames[1]]:
                pset = {pnames[0]: pvals1, pnames[1]: pvals2}
                pset_list.append(pset)
    else:
        for pname, pvals in parameters.iteritems():
            for pval in pvals:
                pset = {pname: pval}
                pset_list.append(pset)

    if 'seed' in StochSSModel.json_data and int(StochSSModel.json_data['seed']) != -1:
        seed = int(StochSSModel.json_data['seed'])
    else:
        random.seed()
        seed = random.randint(0, 2147483647)

    counter = 0
    base_dir = os.path.join('/home/aviral/CSE/qsub_experiments/', "temp_" + str(uuid.uuid4()))
    qsub_file = "/home/aviral/CSE/qsub_experiments/job_submission.pbs"
    job_file = "/home/aviral/CSE/qsub_experiments/ComputeEnsemble.py"
    job_name_prefix = "xyz_ps_job_"
    molns_cloudpickle_file = "/home/aviral/CSE/molnsutil/molnsutil/molns_cloudpickle.py"
    dirs = []
    containers = []

    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    for pndx,pset in enumerate(pset_list):
        unpickled_list = []

        model = StochSSModel(**pset)
        seed_plus_pndx = seed + pndx
        number_of_trajectories = StochSSModel.json_data['trajectories']

        unpickled_list.append(number_of_trajectories)
        unpickled_list.append(seed_plus_pndx)
        unpickled_list.append(model)
        unpickled_list.append(mapper_fn)

        job_name = job_name_prefix + str(counter)
        # create temp directory for this job.
        temp_job_directory = os.path.join(base_dir, job_name + "/")
        if not os.path.exists(temp_job_directory):
            os.makedirs(temp_job_directory)

        # write input file for qsub job.
        output = open(os.path.join(temp_job_directory, "input"), "wb")
        CloudPickle.dump(unpickled_list, output)
        output.close()

        # write job program file.
        shutil.copyfile(job_file, os.path.join(temp_job_directory, "ComputeEnsemble.py"))

        # write molns_cloudpickle.
        shutil.copyfile(molns_cloudpickle_file, os.path.join(temp_job_directory, "molns_cloudpickle.py"))

        containers.append(job_name)
        # invoke qsub to star container with same name as job_name
        Popen(['qsub', '-d', temp_job_directory, '-N', job_name, qsub_file], shell=False)

        dirs.append(temp_job_directory)
        counter += 1
        temp_dirs = dirs[:]

    wait_for_all_results_to_return(temp_dirs)

    print "reducing results.."

    for dir in dirs:
        mapped_list = get_unpickled_result(dir)
        dat.append({ 'parameters' : pset, 'result' : reducer_fn(mapped_list) })

    # remove temporary files and finished containers.
    clean_up([base_dir], containers)

    return dat


if __name__ == '__main__':
    parameters = getParameters()
    print "Parameters: ", parameters

    name = "StochSS_exec" + str(uuid.uuid4())
    print "Name: ", name
    sys.stdout.flush()

    statsSpecies = sorted([specie for specie, doStats in StochSSModel.json_data['speciesSelect'].items() if doStats])
    sys.stdout.write("Starting Parameter sweep\n")
    sys.stdout.flush()
    dat = []
    #############################################
    if StochSSModel.json_data['isLocal']:  # Turn on debugging, run 'local' 
        # sys.stdout.write("Using local serial execution mode\n")
        sys.stdout.write("Using local qsub cluster.\n")
        sys.stdout.flush()
        dat = run_qsub_parameter_sweep(parameters, mapper_fn=mapAnalysis, reducer_fn=reduceAnalysis)
    else:
        sys.stdout.write("Using parallel execution mode (molnsutil.ParameterSweep)\n")
        sys.stdout.flush()
        import molnsutil
        sweep = molnsutil.ParameterSweep(name=name, model_class=StochSSModel, parameters=parameters)
        ret = sweep.run(mapper = mapAnalysis, reducer = reduceAnalysis, number_of_trajectories = StochSSModel.json_data['trajectories'], chunk_size = 1, store_realizations = False, progress_bar = False)
        for r in ret:
            dat.append({ 'parameters' : r.parameters, 'result' : r.result })
    #############################################

    sys.stdout.write("Finished Parameter sweep\n")
    sys.stdout.flush()

    result_file = os.path.join(os.path.dirname(__file__),'results')

    sys.stdout.write("Writing result to file '{0}'\n".format(result_file))
    sys.stdout.flush()

    with open(result_file, 'w') as f:
        pickle.dump(dat, f)

    sys.stdout.write("Done\n\n")
    sys.stdout.flush()
