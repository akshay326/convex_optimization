class Interval:
    bounds = (0, 0)

    def __init__(self, l, u):

        if l > u:
            raise ValueError("Upper bound can't be less than lower")
        else:
            self.bounds = (l, u)

    def __mul__(self, other):
        (a, b) = self.bounds
        (c, d) = other.bounds
        return Interval(min(a * c, a * d, b * c, b * d), max(a * c, a * d, b * c, b * d))

    def __eq__(self, other):
        if self.bounds[0] == other.bounds[0] and self.bounds[1] == other.bounds[1]:
            return True
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return str(self.bounds)


if __name__ == '__main__':
    S = [-2, -1, 0, 1, 2]
    I = []
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            I.append(Interval(S[i], S[j]))

    for i in range(len(I)):
        for j in range(len(I)):
            for k in range(len(I)):
                if (I[i] * I[j]) * I[k] != I[i] * (I[j] * I[k]):
                    print "No"
                    print I[i]
                    print I[j]
                    print I[k]
