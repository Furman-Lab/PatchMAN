import sys, os
import pandas as pd

def main():
    if len(sys.argv) < 3:
        print('The script needs a Rosetta scorefile and at least one name of a score to output. Exiting...')
        exit()

    score_file = sys.argv[1]
    scores = sys.argv[2:]

    score_df = pd.read_csv(score_file, skiprows=0, header=1, sep=r"\s+", low_memory=False)
    score_df['src_file'] = os.path.split(score_file)[1] # not really needed only for compatibility with old code
    score_df = score_df.query('total_score!="total_score"') # sometimes score name lines are repeated, we need to get rid of them
    output = score_df[['src_file'] + scores]

    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
