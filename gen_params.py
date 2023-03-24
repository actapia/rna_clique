import sys
import csv
import itertools

def main():
    writer = csv.writer(sys.stdout, delimiter=" ")
    for line in csv.reader(sys.stdin, delimiter=" "):
        new_row = []
        for v, k in itertools.zip_longest(line, sys.argv[1:]):
            new_row.append(k)
            if v is not None:
                new_row.append(v)
        writer.writerow(new_row)

if __name__ == "__main__":
    main()
