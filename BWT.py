
def rotations(t):
    tt = t * 2
    return [tt[i:i+len(t)] for i in range(0, len(t))]


def bwm(t):
    return sorted(rotations(t))


def bwtViaBwm(t):
    return ''.join(map(lambda x: x[-1], bwm(t)))


def suffixArray(s):
    satups = sorted([(s[i:], i) for i in xrange(0, len(s))])
    return map(lambda x: x[1], satups)

# def generate_fcolumn():
#     fcol = {}
#     bwt_matrix = bwm(str)
#     for row in bwt_matrix:
#         if row[0] in fcol:
#             fcol[row[0]] += 1
#         elif row[0] != "$":
#             fcol[row[0]] = 1
#     return fcol
#
#
# def generate_tally(bwt):
#     tally ={}
#     temp = {}
#     fcol = generate_fcolumn()
#     for e in fcol:
#         tally[e] = []
#         temp[e] = 0
#         for i in range(len(str)/afrac + 1):
#             tally[e].append(0)
#     tally[bwt[0]][0] = 1
#     for i in range(1, len(bwt)):
#         if bwt[i] in temp:
#             temp[bwt[i]] += 1
#         if i % afrac == 0:
#             for c in tally:
#                 tally[c][i/afrac] = tally[c][i/afrac-1] + temp[c]
#                 temp[c] = 0
#     return tally
#
#
# def get_rank(i, c, tally, afrac):
#     if i % afrac == 0:
#         return tally[c][i / afrac]
#     else:
#         count = 0
#         while i % afrac != 0:
#             if bwt[i] == c:
#                 count += 1
#             i -= 1
#         return tally[c][i / afrac] + count
#
#
# def get_row_num(c, fcol, rank):
#     result = 0
#     for e in fcol:
#         if e < c:
#             result += fcol[e]
#     result += (rank + 1)
#     return result
#
#
# def bwt_reverse(afrac, bwt, tally, fcol):
#     result = ""
#     c = bwt[0]
#     i = 0
#     while c!= "$":
#         rank = get_rank(i, c, tally, afrac)
#         result = c + result
#         i = 0
#         for e in fcol:
#             if e < c:
#                 i += fcol[e]
#         i += rank
#         c = bwt[i]
#     return result
#
#
# def query(qstr, tally, fcol, afrac):
#     if qstr == "":
#         return False, (-1, -1)
#     c = qstr[-1]
#     if c not in fcol:
#         return False, (-1, -1)
#     srow = 0
#     for e in fcol:
#         if e < c:
#             srow += fcol[e]
#     srow += 1
#     erow = srow + fcol[c] - 1
#     for i in xrange(len(qstr) - 2, -1, -1):
#         c = qstr[i]
#         if c not in fcol:
#             return False, (srow, erow)
#         srank = get_rank(srow - 1, c, tally, afrac)
#         erank = get_rank(erow, c, tally, afrac)
#         if srank == erank:
#             return False, (srow, erow)
#         srow = get_row_num(c, fcol, srank)
#         erow = srow + erank - srank - 1
#     return True, (srow, erow)
#
# #
# def generate_suffixarray(str, bfrac):
#     sa = suffixArray(str)
#     result = {}
#     for i in range(len(str)):
#         if sa[i] % bfrac == 0:
#             result[i] = sa[i]
#     return result

#
# def get_offset(qstr, sa, bwt, tally, fcol, afrac):
#     is_sub, rows = query(qstr,  tally, fcol, afrac)
#     offsets = []
#     if is_sub:
#         for i in range(rows[0], rows[1] + 1):
#             j = i
#             count = 0
#             while j not in sa:
#                 c = bwt[j]
#                 count += 1
#                 rank = get_rank(j, c, tally, afrac)
#                 j = get_row_num(c, fcol, rank)
#             offsets.append(sa[j] + count)
#     return offsets
# #
#
# #
# str = "QTISGEHGLDGSGVYNGSSDLQLERMNVYFNEASNNKYVPRAVLVDLEPGTMDAVRAGPFGQLFRPDNFVFGQSGAGNNWAKGHYTEG$"
# # # str = "abaaba$"
# #
# # qstr = "SNNKYVPRAV"
# # afrac = 10
# bfrac = 9
# # print bwm(str)
# bwt = bwtViaBwm(str)
# # print bwt
# tally = generate_tally(bwt)
# # print tally
# fcol = generate_fcolumn()
# # print fcol
# # print bwt_reverse(afrac,bwt,tally, fcol)
# # print str
# sa = suffixArray(str)
# print sa
# # print get_offset(qstr,sa, bwt, tally, fcol, afrac)
# print generate_suffixarray(str,bfrac)
# print query(qstr, tally, fcol, afrac)
# # print tally






