#!/usr/bin/env python
import sys

def seq(start, stop):
    # See OEIS A037124
    start_digs = len(str(start))
    start_pow = start_digs - 1
    start_exp = 10**start_pow
    start_coeff = -(-start//start_exp)
    start_mod = start%start_exp
    idx = start_pow * 9 + (start_coeff - 1)
    if start_mod != 0:
        yield start
    start =  (idx%9 + 1)*10**(idx//9)
    while start <= stop:
        yield start
        idx += 1
        start = (idx%9 + 1)*10**(idx//9)
    
    # start_exp = 10**(start_digs-1)
    # start_coeff = -(-start//start_exp)
    # start_mod = start%start_exp
    # if start_mod != 0:
    #     yield start
    # start_pow = (start_coeff // 10)
    # start_exp *= 10**start_pow
    # start_coeff %= 10
    # start = start_coeff * start_exp
    # idx = 9*start_pow
    # while start <= stop:
    #     yield start
    #     start= start_coeff * start_exp
    #     start_coeff += 1
    #     start_exp *= 10**(start_coeff // 10)
    #     start_coeff %= 10
    #     print("SE", start_exp)
    #     print("SC", start_coeff)
    if stop%(10**(len(str(stop))-1)) != 0:
        yield stop

def main():
    for i in seq(int(sys.argv[1]), int(sys.argv[-1])):
        print(i)

if __name__ == "__main__":
    main()
