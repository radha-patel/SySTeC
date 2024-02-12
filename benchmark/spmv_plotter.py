import matplotlib.pyplot as plt
import json
import math

results = json.load(open('/Users/radhapatel/Developer/urop/symmetry-kernels/spmv_results.json', 'r'))

# refxAxis = [math.log10(result["size"]) for result in results if result["method"] == "ref"]
refTime = [result["time"] for result in results if result["method"] == "ref"]
# optxAxis = [math.log10(result["size"]) for result in results if result["method"] == "opt"]
optTime = [result["time"] for result in results if result["method"] == "opt"]
optyAxis = [refTime[i] / optTime[i] for i in range(len(refTime))]

# GET MEAN
print(sum(optyAxis)/len(optyAxis))
# plt.grid(True)

# ## LINE GRAPH ##
# plt.plot(optxAxis, optyAxis, color='green', marker='o')
# plt.xlabel('Size (Log)')
# plt.ylabel('Speed Up')

# ## BAR GRAPH ##
# fig = plt.figure()
# plt.bar(xAxis,yAxis, color='maroon')
# plt.xlabel('variable')
# plt.ylabel('value')

# plt.show()