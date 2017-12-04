from BWT import *


class FMIndex:
    def __init__(self, string, afrac, bfrac):
        self.string = string
        self.afrac = afrac
        self.bfrac = bfrac
        self.bwt = bwtViaBwm(self.string)
        self.fcol = self.generate_fcolumn()
        self.tally = self.generate_tally()
        self.sa = self.generate_suffixarray()
        # self.sa = suffixArray(string)

    def generate_fcolumn(self):
        fcol = {}
        bwt_matrix = bwm(self.string)
        for row in bwt_matrix:
            if row[0] in fcol:
                fcol[row[0]] += 1
            elif row[0] != "$":
                fcol[row[0]] = 1
        return fcol

    def generate_tally(self):
        tally = {}
        temp = {}
        for e in self.fcol:
            tally[e] = []
            temp[e] = 0
            for i in range(len(self.string) / self.afrac + 1):
                tally[e].append(0)
        tally[self.bwt[0]][0] = 1
        for i in range(1, len(self.bwt)):
            if self.bwt[i] in temp:
                temp[self.bwt[i]] += 1
            if i % self.afrac == 0:
                for c in tally:
                    tally[c][i / self.afrac] = tally[c][i / self.afrac - 1] + temp[c]
                    temp[c] = 0
        return tally

    def get_rank(self, i, c):
        if i % self.afrac == 0:
            return self.tally[c][i / self.afrac]
        else:
            count = 0
            while i % self.afrac != 0:
                if self.bwt[i] == c:
                    count += 1
                i -= 1
            return self.tally[c][i / self.afrac] + count

    def get_row_num(self, c, rank):
        result = 0
        for e in self.fcol:
            if e < c:
                result += self.fcol[e]
        result += (rank + 1)
        return result

    def bwt_reverse(self):
        result = ""
        c = self.bwt[0]
        i = 0
        while c != "$":
            rank = self.get_rank(i, c)
            result = c + result
            i = 0
            for e in self.fcol:
                if e < c:
                    i += self.fcol[e]
            i += rank
            c = self.bwt[i]
        return result

    def query(self, qstr):
        if qstr == "":
            return False, (-1, -1)
        c = qstr[-1]
        if c not in self.fcol:
            return False, (-1, -1)
        srow = 0
        for e in self.fcol:
            if e < c:
                srow += self.fcol[e]
        srow += 1
        erow = srow + self.fcol[c] - 1
        for i in xrange(len(qstr) - 2, -1, -1):
            c = qstr[i]
            if c not in self.fcol:
                return False, (srow, erow)
            srank = self.get_rank(srow - 1, c)
            erank = self.get_rank(erow, c)
            if srank == erank:
                return False, (srow, erow)
            srow = self.get_row_num(c, srank)
            erow = srow + erank - srank - 1
            if erow >= len(self.string):
                return False, (-1, -1)
        return True, (srow, erow)

    def generate_suffixarray(self):
        sa = suffixArray(self.string)
        result = {}
        for i in range(len(self.string)):
            if sa[i] % self.bfrac == 0:
                result[i] = sa[i]
        return result

    def get_offset(self, qstr):
        is_sub, rows = self.query(qstr)
        offsets = []
        if is_sub:
            for i in range(rows[0], rows[1] + 1):
                j = i
                count = 0
                while j not in self.sa:
                    c = self.bwt[j]
                    count += 1
                    rank = self.get_rank(j, c) - 1
                    j = self.get_row_num(c, rank)
                offsets.append(self.sa[j] + count)
        return offsets