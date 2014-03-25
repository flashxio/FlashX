def read_igraph():
  # Get igraph pagerank
  f = open("igraph_pgrank.txt", "rb")
  arr = f.read().splitlines()
  f.close()
  igraph_dict = {}

  for idx, entry in enumerate(arr):
    igraph_dict[idx+1] = float(entry)
  return igraph_dict

def read(fn):
  # Get Flash_graph pagerank
  flashgraph_dict = {}
  f = open(fn, "rb")
  arr = f.read().splitlines()
  for entry in arr:
    key, val = entry.split(":")
    flashgraph_dict[int(key)] = float(val)
  return flashgraph_dict

def compare(essential_dict, non_essential_dict, essential, skip=False):
  if not skip:
    assert set(essential_dict.keys()) == set(non_essential_dict.keys()), "The sets of keys are different"
  
  count = 0
  for key in essential_dict.keys():
    if not within(essential_dict[key], non_essential_dict[key], 0.05):
      print "Mismatch Vertex id: %d. %s: %f != FLASHGRAPH: %f" % (key, essential, essential_dict[key], non_essential_dict[key])
      count += 1

  print "The number of mismatches between FlashGraph and %s=%d" % ( essential.capitalize(), count)

def within(arg1, arg2, thresh):
  return abs(arg1-arg2) < thresh

# main
compare(read_igraph(), read("flashgraph_pgrank.txt"), "IGRAPH")
compare(read("powergraph_pgrank.txt"), read_igraph(), "POWERGRAPH", True)
