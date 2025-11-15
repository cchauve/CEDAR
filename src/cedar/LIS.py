"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
Longest Increasing Subsequence with O(n log(n)) algorithm
Taken from https://leetcode.com/problems/longest-increasing-subsequence/solutions/1326308/c-python-dp-binary-search-bit-segment-tree-solutions-picture-explain-o-nlogn/
"""
__author__ = "Cedric Chauve"
__credits__ = ["Mai Thanh Hiep (https://leetcode.com/u/hiepit/, https://github.com/hiepxuan2008)"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

from bisect import bisect_left

def LIS_len(s):
    """
    Given a sequence s of integers, returns the length of a longest increasing subsequence
    Input:
    - s (list(int))
    Output:
    - (int): lengh of a longest increasing subsequence
    """
    if len(s) == 0:
        return 0
    sub = []
    for x in s:
        if len(sub) == 0 or sub[-1] < x:
            sub.append(x)
        else:
            idx = bisect_left(sub, x)  # Find the index of the first element >= x
            sub[idx] = x  # Replace that number with x
    return len(sub)

def LIS_seq(s):
    """
    Given a sequence s of integers, returns a longest increasing subsequence
    Input:
    - s (list(int))
    Output:
    - list(int): a longest increasing subsequence
    """
    if len(s) == 0:
        return []
    sub = []
    subIndex = []  # Store index instead of value for tracing path purpose
    trace = [-1] * len(s)  # trace[i] point to the index of previous number in LIS
    for i, x in enumerate(s):
        if len(sub) == 0 or sub[-1] < x:
            if subIndex:
                trace[i] = subIndex[-1]
            sub.append(x)
            subIndex.append(i)
        else:
            idx = bisect_left(sub, x)  # Find the index of the smallest number >= x, replace that number with x
            if idx > 0:
                trace[i] = subIndex[idx - 1]
            sub[idx] = x
            subIndex[idx] = i
    path = []
    t = subIndex[-1]
    while t >= 0:
        path.append(s[t])
        t = trace[t]
    return path[::-1]
