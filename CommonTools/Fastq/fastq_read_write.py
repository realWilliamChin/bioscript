#!/usr/bin/env python3

import gzip
import os
import sys

def isFastq(f):
    fqext = (".fq", ".fastq", "fq.gz", ".fastq.gz")
    for ext in fqext:
        if f.endswith(ext):
            return True
    return False

################################
# fastq.reader

class Reader:

    def __init__(self, fname):
        self.__file = None
        self.__gz = False
        self.__eof = False
        self.filename = fname
        if self.filename.endswith(".gz"):
            self.__gz = True
            self.__file = gzip.open(self.filename, "rt")  # 使用文本模式读取
        else:
            self.__gz = False
            self.__file = open(self.filename, "r")
        if self.__file is None:
            print("Failed to open file " + self.filename)
            sys.exit(1)

    def __del__(self):
        if self.__file is not None:
            self.__file.close()

    def nextRead(self):
        if self.__eof or self.__file is None:
            return None

        lines = []
        # 读取 4 行 (name, sequence, strand, quality)
        for i in range(0, 4):  # 使用 range 替换 xrange
            line = self.__file.readline().rstrip()
            if len(line) == 0:
                self.__eof = True
                return None
            lines.append(line)
        return lines

    def isEOF(self):
        return False

################################
# fastq.writer

class Writer:

    filename = ""

    __file = None
    __gz = False

    def __init__(self, fname):
        self.filename = fname
        if self.filename.endswith(".gz"):
            self.__gz = True
            self.__file = gzip.open(self.filename, "wt")  # 使用文本模式写入
        else:
            self.__gz = False
            self.__file = open(self.filename, "w")
        if self.__file is None:
            print("Failed to open file " + self.filename + " to write")
            sys.exit(1)

    def __del__(self):
        if self.__file is not None:
            self.__file.flush()
            self.__file.close()

    def flush(self):
        if self.__file is not None:
            self.__file.flush()

    def writeLines(self, lines):
        if self.__file is None:
            return False

        for line in lines:
            self.__file.write(line + "\n")
        return True

    def writeRead(self, name, sequence, strand, quality):
        if self.__file is None:
            return False

        self.__file.write(name + "\n")
        self.__file.write(sequence + "\n")
        self.__file.write(strand + "\n")
        self.__file.write(quality + "\n")

        return True
