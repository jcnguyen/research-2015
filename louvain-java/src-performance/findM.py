########################
# findM.py is used to assist in the determination of the meaningful
# maximum M, used in weighted performance calculations.
# 
# Input: Call with 1 system argument, the path to the weighted graph 
# input file. Input file should be in format: edge1 edge2 weight.
#
# Output: prints n, median, min, max, iqr, mild/extreme outlier fences, 
# 			and outlier counts
#
########################
import sys
import numpy as np

# open input file
input = open(sys.argv[1], 'r')

print "\n================================================"
file_path = sys.argv[1].split('/')
file_name = file_path[len(file_path)-1]
print "Reading from", file_name

# get all weights as floats
num_ar = []
for line in input:
	num_ar.append( float(line.split(' ')[2].replace("\n", "")) )

# clean up - close file
input.close()

# calculate statistical information
n = len(num_ar)
median = np.median(num_ar)
q1 = np.percentile(num_ar, 25)
q3 = np.percentile(num_ar, 75)
iqr = q3-q1
lower_inner_fence = q1 - 1.5*iqr
upper_inner_fence = q3 + 1.5*iqr
lower_outer_fence = q1 - 3*iqr
upper_outer_fence = q3 + 3*iqr

# print stats info
print "\nn =", n
print "Median", median
print "Min", min(num_ar)
print "Max", max(num_ar)
print "IQR of", iqr, "is from", q1, "to", q3
print "Fences for outliers:"
print "\tUpper, extreme", upper_outer_fence
print "\tUpper, mild", upper_inner_fence
print "\tLower, extreme", lower_outer_fence
print "\tLower, mild", lower_inner_fence

# determine if outliers exist
mild_outliers_low = []
mild_outliers_high = []
extreme_outliers_low = []
extreme_outliers_high = []
for x in num_ar:
	if x < lower_outer_fence:
		extreme_outliers_low.append(x)
	elif x < lower_inner_fence:
		mild_outliers_low.append(x)
	elif x > upper_outer_fence:
		extreme_outliers_high.append(x)
	elif x > upper_inner_fence:
		mild_outliers_high.append(x)

# print outlier results
n = float(n)
print "Results for outliers:"
print "\tUpper, extreme", len(extreme_outliers_high), float(len(extreme_outliers_high))/n*100, "%"
print "\tUpper, mild", len(mild_outliers_high), float(len(mild_outliers_high))/n*100, "%"
print "\tLower, extreme", len(extreme_outliers_low), float(len(extreme_outliers_low))/n*100, "%"
print "\tLower, mild", len(mild_outliers_low), float(len(mild_outliers_low))/n*100, "%"
print "================================================\n"