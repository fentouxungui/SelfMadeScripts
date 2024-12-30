import re
import zipfile


# read the zip file
def zipReader(file):
    qcfile = zipfile.ZipFile(file)
    data_txt = [file for file in qcfile.namelist() if re.match(".*?_data\.txt", file)][0]
    data = [bytes.decode(line) for line in qcfile.open(data_txt)]
    return data


def fastqc_summary(data):
    module_num = 0
    bases = 0
    Q20 = 0
    Q30 = 0
    for line in data:
        if re.match('Filename', line):
            filename = line.split(sep="\t")[1].strip()
        if re.match('Total Sequence', line):
            read = line.split(sep="\t")[1].strip()
        if re.match('%GC', line):
            GC = line.split(sep="\t")[1].strip()
        if re.match("[^#](.*?\t){6}", line):
            bases = bases + 1
            if float(line.split("\t")[1]) > 30:
                Q20 = Q20 + 1
                Q30 = Q30 + 1
            elif float(line.split("\t")[1]) > 20:
                Q20 = Q20 + 1

        if re.match(">>END", line):
            module_num = module_num + 1
            if module_num >= 2:
                break
    Q20 = Q20 / bases
    Q30 = Q30 / bases
    summary = [filename, read, GC, str(Q20), str(Q30)]
    return summary


if __name__ == '__main__':
    import sys

    for arg in range(1, len(sys.argv)):
        data = zipReader(sys.argv[arg])
        summary = fastqc_summary(data)
        with open('summary.txt', 'a') as f:
            f.write('\t'.join(summary) + '\n')