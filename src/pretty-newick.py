#!/usr/bin/env python3

import argparse
import sys

def format_newick(s, indent=0):
    out = ""
    i = 0
    while i < len(s):
        c = s[i]
        if c == '(':
            out += '\n' + '    ' * indent + '('
            indent += 1
        elif c == ',':
            out += ',\n' + '    ' * indent
        elif c == ')':
            indent -= 1
            out += '\n' + '    ' * indent + ')'
        else:
            out += c
        i += 1
    return out

def main():
    parser = argparse.ArgumentParser(description="Pretty-print Newick tree with indentation.")
    parser.add_argument("input", nargs="?", type=argparse.FileType('r'), default=sys.stdin,
                        help="Input Newick file (default: stdin)")
    args = parser.parse_args()

    newick_str = args.input.read().strip()
    print(format_newick(newick_str))

if __name__ == "__main__":
    main()