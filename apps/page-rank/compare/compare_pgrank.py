def read_igraph():
  # Get igraph pagerank
  f = open("igraph_pgrank.txt", "rb")
  arr = f.read().splitlines()
  f.close()
  igraph_dict = {}

  for idx, entry in enumerate(arr):
    igraph_dict[idx+1] = float(entry)
  return igraph_dict

def read_flashgraph():
  # Get Flash_graph pagerank
  flashgraph_dict = {}
  f = open("flashgraph_pgrank.txt", "rb")
  arr = f.read().splitlines()
  for entry in arr:
    key, val = entry.split(":")
    flashgraph_dict[int(key)] = float(val)
  return flashgraph_dict

def compare(igraph_dict, flashgraph_dict):
  assert set(igraph_dict.keys()) == set(flashgraph_dict.keys()), "The sets of keys are different"
  
  count = 0
  for key in igraph_dict.keys():
    if not within(igraph_dict[key], flashgraph_dict[key], 0.01):
      print "Mismatch Vertex id: %d. IGRAPH: %f != FLASHGRAPH: %f" % (key, igraph_dict[key], flashgraph_dict[key])
      count += 1

  print "The number of mismatches =", count

def within(arg1, arg2, thresh):
  return abs(arg1-arg2) < thresh

# main
compare(read_igraph(), read_flashgraph())
