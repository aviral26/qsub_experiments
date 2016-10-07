#! /usr/bin/python
import pickle
import molns_cloudpickle as CloudPickle

with open("input", "rb") as inpput:
    unpickled_list = pickle.load(inpput)
inpput.close()


number_of_trajectories = unpickled_list[0]

seed_plus_pndx = unpickled_list[1]

model = unpickled_list[2]

mapper_fn = unpickled_list[3]

results = model.run(number_of_trajectories = number_of_trajectories, seed=seed_plus_pndx, show_labels = True)
if not isinstance(results, list):
    results = [results]

mapped_list = []

for r in results:
    mapped_list.append(mapper_fn(r))


output = open("output", "wb")
CloudPickle.dump(mapped_list, output)
