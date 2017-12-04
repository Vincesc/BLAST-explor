import copy
from suffix_tree import SuffixTree
import singlescore
from fm_index import FMIndex

amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# def get_score_matrix():
#     read = open("score_matrix", "r")
#     amino_acid = read.readline()
#     score = read.readline()
#     score_matrix = {}
#     while score != "":
#         amino = score[0]
#         for i in range(len(amino_acid)):
#             if amino_acid[i] != " " and amino_acid[i] != "\n":
#                 if score[i + 2] == "-":
#                     s = - int(score[i + 3])
#                 else:
#                     if score[i + 2] != " ":
#                         s = int(score[i + 3]) + 10 * int(score[i + 2])
#                     else:
#                         s = int(score[i + 3])
#                 if amino not in score_matrix:
#
#                     score_matrix[amino] = {amino_acid[i]: s}
#                 else:
#                     score_matrix[amino][amino_acid[i]] = s
#         score = read.readline()
#     return score_matrix
#
#
# score_matrix = get_score_matrix()
# print (score_matrix)
# def get_single_score(x,y):
#     return score_matrix[x][y]
smatrix=singlescore.createScoreM(singlescore.scores)


def get_single_score(x, y):

    ii = singlescore.alp.index(x)
    jj = singlescore.alp.index(y)
    return int(smatrix[ii][jj])


def get_scores(s, t):
    score = 0
    for i in range(len(s)):
        score += get_single_score(s[i], t[i])
    return score


def get_k_mers(query_seq, k):
    n = len(query_seq)
    if k > n:
        return []
    result = []
    for i in range(k, n+1):
        result.append(query_seq[i-k:i])
    return result


def generate_neighbors(word, threshold):
    k = len(word)
    collection = list(amino_acids)
    result = []
    for i in range(k-1):
        temp = copy.deepcopy(collection)
        collection = []
        for s in temp:
            for a in amino_acids:
                collection.append(s+a)
    for k_mer in collection:
        if get_scores(k_mer, word) >= threshold:
            result.append(k_mer)
    return result


def generate_all_neighbors(file_name, k, threshold):
    file_read = open(file_name, "r")
    query_sequence = ""
    while True:
        s = file_read.readline()
        if s == "":
            break
        if s[-1] == "\n":
            s = s[:-1]
        query_sequence += s
    k_mers = get_k_mers(query_sequence, k)
    neighbors = {}
    for i in range(len(k_mers)):
        neighbor = generate_neighbors(k_mers[i], threshold)
        for n in neighbor:
            if n not in neighbors:
                neighbors[n] = [i,]
            else:
                neighbors[n].append(i)
    return neighbors


def get_protein_data(file_name, protein_name):
    read_sequence = open(file_name, "r")
    result = []

    protein_name.append(read_sequence.readline()[:-1])
    while True:
        seq = read_sequence.readline()
        times = 0
        protein = ""
        while len(seq) > 0 and seq[0] != ">":
            times += 1
            protein += seq[:-1]
            seq = read_sequence.readline()
        protein_name.append(seq[:-1])

        if protein != "":
            result.append(protein)
        if times == 0:
            break
    return result


def get_match_position_with_suffixtree(neighbors, suffix_trees):
    match_position = []
    for i in range(len(suffix_trees)):
        match_position.append({})
        for n in neighbors:
            if suffix_trees[i].has_substring(n):
                p = suffix_trees[i].find_substring(n)
                if n not in match_position[i]:
                    match_position[i][n] = [p]
                else:
                    match_position[i][n].append(p)
    return match_position

def get_match_position_with_fmindex(neighbors, fm_indexs):
    match_position = []
    for i in range(len(fm_indexs)):
        match_position.append({})
        for n in neighbors:
            p = fm_indexs[i].get_offset(n)
            if len(p) > 0:
                if n not in match_position[i]:
                    match_position[i][n] = p
                else:
                    match_position[i][n].extend(p)
    return match_position


def extending(query, q, sequence, s, k, drop):

    query_seg = query[q:q+k]
    sequence_seg=sequence[s:s+k]
    comparison_seg = ""
    pairscore = 0
    for i in range(0,k):
        score = get_single_score(query[q+i], sequence[s+i])
        pairscore += score
        if query[q+i] == sequence[s+i]:
            comparison_seg += query[q+i]
        else:
            if score > 0:
                comparison_seg += "+"
            else:
                comparison_seg += " "
    left = -1
    while q+left > -1 and s+left > -1:
        mid_pair_score=get_single_score(query[q+left], sequence[s+left])
        if mid_pair_score < -drop+1:
            break
        else:
            pairscore+=mid_pair_score
            query_seg=query[q+left]+query_seg
            sequence_seg=sequence[s+left]+sequence_seg
            if query[q+left] == sequence[s+left]:
                comparison_seg = query[q+left] + comparison_seg
            else:
                if mid_pair_score > 0:
                    comparison_seg = "+" + comparison_seg
                else:
                    comparison_seg = " " + comparison_seg
            left -= 1
    right = 1
    while q+k+right<len(query) and s+k+right<len(sequence):
        mid_pair_score = get_single_score(query[q+k+right], sequence[s+k+right])
        if mid_pair_score < -drop+1:
            break
        else:
            pairscore+=mid_pair_score
            query_seg=query_seg+query[q+right]
            sequence_seg=sequence_seg+sequence[s+right]
            if query[q+right] == sequence[s+right]:
                comparison_seg = comparison_seg + query[q+right]
            else:
                if mid_pair_score > 0:
                    comparison_seg = comparison_seg + "+"
                else:
                    comparison_seg = comparison_seg + " "
            right+=1
    scores=[]
    # for i in range(0,len(query_seg)):
    #     scores.append(get_single_score(query_seg[i], sequence_seg[i]))
    # print scores

    query_offset = q + left + 1
    sequence_offset = s + left + 1

    return pairscore,query_offset, sequence_offset, comparison_seg, query_seg, sequence_seg


def main():
    query_neighbors = generate_all_neighbors("querry_sequence", 3, 11)

    protein_name = []

    protein_sequences = get_protein_data("protein.fasta", protein_name)

    # print(protein_name)

    # proteins_suffix_trees = []
    # for i in range(len(protein_sequences)):
    #      proteins_suffix_trees.append(SuffixTree(protein_sequences[i]))
    #
    # positions = get_match_position_with_suffixtree(query_neighbors,proteins_suffix_trees)

    # print positions
    afrac = 10
    bfrac = 8
    protein_fmindexs =[]
    for i in range(len(protein_sequences)):
         protein_fmindexs.append(FMIndex(protein_sequences[i]+"$", afrac, bfrac))
    positions = get_match_position_with_fmindex(query_neighbors,protein_fmindexs)
    print positions

    # fwrite = open("FM.txt", "w")
    # for

    # print positions
    read = open("query_sequence", "r")
    query_sequnce =  read.readline()


    # print protein_sequences

    # print query_sequnce
    #
    # print query_neighbors
    # print position
    #
    result = []

    for i in range(len(protein_sequences)):
        result.append({})
        for k_mer in query_neighbors:
            for j in range(len(query_neighbors[k_mer])):
                if k_mer in positions[i]:
                    for pos in positions[i][k_mer]:
                        score, q_offset, s_offset, comp_seg, q_seg, s_seg = extending(query_sequnce, query_neighbors[k_mer][j], protein_sequences[i], pos, 3, 3)
                        list = [q_offset, s_offset, comp_seg, q_seg, s_seg]
                        if (score in result):
                            result[i][score].append(list)
                        else:
                            result[i][score] = list
    for i in range(len(protein_sequences)):
        print(protein_name[i])
        display = result[i]
        sort_order = sorted(display, reverse= True)
        for j in sort_order:
            k = 0
            while k < len(display[j]):
                score = j
                q_offset = display[j][k]
                k+=1
                s_offset = display[j][k]
                k+=1
                comp_seg = display[j][k]
                k+=1
                q_seg = display[j][k]
                k+=1
                s_seg = display[j][k]
                k+=1
                print "Score: ", score
                print "Query:  ", "{0:4d}".format(q_offset), "  ", q_seg, "  ", q_offset+len(q_seg)
                print "                ", comp_seg
                print "Sbjct:  ", "{0:4d}".format(s_offset), "  ", s_seg, "  ", s_offset + len(s_seg)
                print ""
        print("-----------------------------------------------------------------------------------------")
    # print extending(query_sequnce,query_neighbors['PFG'][0],protein_sequences[0],positions[0]['PFG'],3,3)

if __name__ == "__main__":
    main()