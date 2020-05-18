def surround(s1, s2, n, c):
    if n == 0:
        return s1+s2
    else:
        return surround(c+s1+c,s2+c,n-1,c)

print(surround("Joe", "Parker", 2, "*"))
print(surround("Joe", "Parker", 0, "*"))
print(surround("CMPS", "200", 3, "-"))
