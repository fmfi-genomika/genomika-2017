#!/usr/bin/python3

import sys
import logging

def parse_psl_matches(filename):
    #match   mis-    rep.    N's Q gap   Q gap   T gap   T gap   strand  Q           Q       Q       Q   T           T       T       T   block   blockSizes  qStarts  tStarts
    #    match   match       count   bases   count   bases           name        size    start   end name        size    start   end count
    T_SIZE_COL = 14
    Q_SIZE_COL = 10
    T_NAME_COL = 13
    Q_NAME_COL = 9
    MATCH_SIZE_COL = 0
    MISMATCH_SIZE_COL = 1 

    f = open(filename, "r")
    lines = [l for l in f]
    lines = lines[5:]
    lines = [l.strip().split("\t") for l in lines]
    
    total_q_number = len({l[Q_NAME_COL] for l in lines})
    total_t_number = len({l[T_NAME_COL] for l in lines})

    q2t = {}
    t2q = {}

    counter = 0
    for line in lines:
        if abs(int(line[T_SIZE_COL]) - int(line[Q_SIZE_COL])) < 1:
            if int(line[T_SIZE_COL]) == int(line[MATCH_SIZE_COL]):
                counter += 1
                #print(line)
                if line[T_NAME_COL] not in t2q: t2q[line[T_NAME_COL]] = set()
                if line[Q_NAME_COL] not in q2t: q2t[line[Q_NAME_COL]] = set()
                t2q[line[T_NAME_COL]].add(line[Q_NAME_COL])
                q2t[line[Q_NAME_COL]].add(line[T_NAME_COL])

    logging.debug("matches: {}, total q: {}, total t: {}".format(counter, total_q_number, total_t_number))
    logging.debug("nonunique q2t: {}".format(sum([1 if len(q2t[x]) > 1 else 0 for x in q2t])))
    logging.debug("nonunique t2q: {}".format(sum([1 if len(t2q[x]) > 1 else 0 for x in t2q])))

    # print only unique pairs
    for q, ts in q2t.items():
        if len(ts) == 1:
            t = list(ts)[0]
            if len(t2q[t]) == 1:
                print("{}; {}".format(q, t))


def main():
    logging.basicConfig(level=logging.DEBUG)
    if len(sys.argv) != 2:
        logging.error("WRONG ARGUMENT COUNT!")
        sys.exit(1)

    input_filename = sys.argv[1]
    
    parse_psl_matches(input_filename)

if __name__ == "__main__":
    main()