#! python3
import sys, os

def main():
    if len(sys.argv) < 3:
        print('The script needs a Rosetta scorefile and at least one name of a score to output. Exiting...')
        exit()

    score_file = sys.argv[1]
    scores = sys.argv[2:]
    src_file = os.path.split(score_file)[1]

    # process line-by line, skip lines with SEQUENCE
    print('src_file' + ' ' + " ".join(scores))
    with open(score_file, 'r') as f:
        header = []
        indices = []
        for line in f:
            if line.startswith('SEQUENCE'):
                continue
            if line.startswith('SCORE'):
                if len(header) == 0:
                    header = line.strip().split()
                    for score in scores:
                        try:
                            idx = header.index(score)
                            indices.append(idx)
                        except:
                            continue
                else:
                    # select keys defined by columns from splitted line
                    split_line = line.strip().split()
                    values = [split_line[i] for i in indices]
                    print(src_file + ' ' + " ".join(values))


if __name__ == "__main__":
    main()

