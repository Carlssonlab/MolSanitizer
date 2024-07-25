def minmax(a):
    minval = a[0]
    maxval = a[0]
    for i in range(1, len(a)):
        val = a[i]
        if val > maxval:
            maxval = val
        if val < minval:
            minval = val
    return minval, maxval