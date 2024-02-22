#!/share/data1/software/miniconda3/bin/python
import re
import argparse
import copy
#from xml.etree.ElementInclude import include

'''
整体思路就是往排序后的区间列表中添加的区间
    1、找出与新加区间重叠的区间
    2、对这些有交集的区间进行合并
    3、返回新添加的长度以及新的其实结束位点
'''

def get_args():
    parser = argparse.ArgumentParser(prog=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("i", metavar="in_f", help="blast.m8")
    parser.add_argument("o", metavar="out_f", help="output")
    parser.add_argument("min_id",type=float, default=75, help="ignore the identity less than min_id[default 75]")
    parser.add_argument("min_len",type=float, default=200, help="ignore the identity less than min_len[default 200]")
    args = parser.parse_args()
    return args

def process_(start:int, end:int, include_section: list) -> tuple: # (start, end, new_len)
    nl = 0
    x1 = 0
    st, en = start, end
    include_section_len = len(include_section)
    for i,v in enumerate(include_section):
        if v[0] <= st <= v[1]:
            st = v[1]+1
        if st < v[0]:
            x1 = v[0] - st
            st = v[1] + 1
        if (include_section_len == i+1) and v[1] < en:
            x1 += en - v[1]
        nl += x1
    start = min(include_section[0][0],start)
    end = max(end, include_section[-1][1])
    return(start, end, nl)

def add_section(start: int, end: int, target: dict, similar:int, name: str):
    section_name = f"{name}_section"
    similar_name = f"{name}_total_similar"
    length_name  = f"{name}_total_len"
    include_section = [] # 与新区间有重叠部分的
    other_section = [] # 与新的区间没有重合的
    # 如果区间不为空
    if target[section_name]!=[]:
        for section in target[section_name]:
            # 找出来有交集的区间
            if not(end < section[0] or start > section[1]):
                include_section.append(section)
                continue
            other_section.append(section)
    # 如果有交集
    if len(include_section) != 0:
        # 都是按照大小排序，并且添加的区间与原来的区间有交际，所以考虑的情况比较少
        start, end, new_len = process_(start, end, include_section)
    else:
        # 否则添加的长度公式如下
        new_len = end - start + 1
    other_section.append([start, end])
    #print(other_section)
    other_section = sorted(other_section, key=(lambda x:x[0]))
    target[length_name] += new_len
    target[similar_name] += new_len * similar
    target[section_name] = other_section

def process_line(q_start: int, q_end: int, s_start:int, s_end:int, similar:int, target_dict: dict, current_name:str=None) -> None:
    '''将比对上的区间添加到列表中'''
    add_section(q_start, q_end, target_dict, similar,"q")
    add_section(s_start, s_end, target_dict, similar,"s")

def print_result(result: dict, name, o_f) -> None:
    q_len = result['q_total_len']
    s_len = result['s_total_len']
    q_similar = result['q_total_similar']/q_len * 100 # python 默认是四舍六入五留双
    s_similar = result['s_total_similar']/s_len * 100
    # q_similar = int(result['q_total_similar']/q_len * 100* 100+0.5)/100 
    # s_similar = int(result['s_total_similar']/s_len * 100* 100+0.5)/100
    out_str = f"{name}\t{q_len}\t{q_similar:.2f}\t{s_len}\t{s_similar:.2f}\n"
    o_f.write(out_str)
    #print(out_str)

def main(in_f: str, cut_id:float, cut_len:int, out_f:str):
    temp_name = ""
    last_result_temp = {"name":"", "q_section":[], "s_section":[], "q_total_similar":0, "q_total_len":0, "s_total_similar":0,"s_total_len":0}
    last_result = copy.deepcopy(last_result_temp)
    f = open(in_f, 'r')
    o_f = open(out_f, "w")
    for line in f:
        line_split = re.split("\t", line.strip())
        qid, sid, qsimilar = line_split[0], line_split[1], float(line_split[2])
        match, mismatch, gap, qstart, qend, sstart, send = \
            [int(line_split[i]) for i in (3,4,5,6,7,8,9)]

        if qsimilar < cut_id or match < cut_len:
            continue

        if qstart > qend: qstart, qend = qend, qstart
        if sstart > send: sstart, send = send, sstart

        similar = 1 - (mismatch + gap ) / match
        current_name = f"{qid}\t{sid}"
        if temp_name != current_name :
            if temp_name != "": print_result(last_result, temp_name, o_f)
            temp_name = current_name
            last_result = copy.deepcopy(last_result_temp) # 结果清零
        process_line(qstart, qend, sstart, send, similar, last_result)
    print_result(last_result, temp_name, o_f)


if __name__ == "__main__":
    args = get_args()
    main(args.i, args.min_id, args.min_len, args.o)
