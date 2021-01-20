import os
import argparse


parser = argparse.ArgumentParser(
    description='Extract selected A domain sequences from .fasta polypeptide file')
parser.add_argument('--hmmfile', required=True, help="Choose hmmfile you want annotations to be conducted with")
parser.add_argument('--file', required=True, help="Choose .fasta file containing sequences for annotation and extraction")
parser.add_argument('--domainposition', required=True, help="Choose which occurence of A domain required (i.e. 2 if the A domain wanted is the second A domain in the sequence)")
args = parser.parse_args()


def hmmcheck():
    # change the name of the input alinment file below
    f = open(args.file)
    a3mf = f.readlines()
    f.close()
    cntr = 0
    # change the name of the output alignment file below
    out = open('annotated.fasta', 'w')
    names = []

    for i in range(len(a3mf)/2):
        # for i in range(1):
        cntr = cntr+1
        outname = open('oneseq.fasta', 'w')
        outname.write(a3mf[i*2])
        outname.write(a3mf[i*2+1])
        outname.close()
        # runs hmmscan for one sequence and writes output to oneseq.fasta
        hmmscanname = 'oneseq.hmmscan'
        hmmscancommand = "hmmscan -E 0.001 -o "+hmmscanname+" --noali --incE 0.001 " + \
            args.hmmfile+" oneseq.fasta"  # change the Pfam-A.hmm database directory
        os.system(hmmscancommand)

        # gathers line of the domain information in domainlines list
        f = open(hmmscanname)
        hmmscan = f.readlines()
        f.close()
        hmmdomains = []
        kr_doms = []
        domainlines = []
        for s in range(len(hmmscan)):
            columns = hmmscan[s].strip().split()
            if len(columns) > 0 and columns[0] == '>>':
                domainlines.append(s)

        # print domainlines

        # goes one by one on each detected domain, since
        ranges = []
        overlapcheck = []
        for k in range(len(domainlines)-1):
            l1 = domainlines[k]
            l2 = domainlines[k+1]
            for s in range(l1, l2):
                if len(hmmscan[s].split()) > 0:
                    if hmmscan[s].split()[0] == '>>':
                        domainname = hmmscan[s].split()[1]
                    if hmmscan[s].split()[0] != '#' and hmmscan[s].split()[0] != '---' and hmmscan[s].split()[0] != '>>':
                        iEvalue = float(hmmscan[s].split()[5])
                        if iEvalue <= 0.0000001:
                            alS = int(hmmscan[s].split()[9])
                            alF = int(hmmscan[s].split()[10])
                            if len(ranges) == 0:
                                start = str(alS)
                                end = str(alF)
                                ranges.append(list(range(alS, alF)))
                                hmmdomains.append(domainname)
                            elif len(ranges) > 0:
                                for arange in ranges:
                                    if set(range(alS, alF)).isdisjoint(arange):
                                        overlapcheck.append(0)
                                    elif not set(range(alS, alF)).isdisjoint(arange):
                                        overlapcheck.append(1)
                                if 1 not in overlapcheck:
                                    #						    if (alS,alF) not in ranges:
                                    # print domainname
                                    start = str(alS)
                                    # print start
                                    end = str(alF)
                                    ranges.append(list(range(alS, alF)))
                                    hmmdomains.append(domainname)

        out.write(a3mf[i*2][:-1]+'|')
        out.write('.'.join(sorted(hmmdomains)))

        out.write('\n')
        out.write(a3mf[i*2+1])

    out.close()
    f.close()


def keepseq():
    f = open('annotated.fasta')
    out = open('filtered.fasta', 'w')
    standard = ""
    p = 0
    for l in f:
        if p == 1:
            if ">" not in l:
                out.write(l)
        if ">" in l:
            p = 0
            s = l.split("|")
            doms = s[1]
            # print doms
            if standard == "":
                standard = doms
            if doms == standard:
                p = 1
                out.write(s[0]+"\n")

    out.close()
    f.close()


def labeldom():
    # change the name of the input alinment file below
    f = open('filtered.fasta')
    a3mf = f.readlines()
    f.close()
    cntr = 0
    # change the name of the output alignment file below
    out = open('domains.fasta', 'w')
    names = []

    for i in range(len(a3mf)/2):
        # for i in range(1):
        cntr = cntr+1
        outname = open('oneseq.fasta', 'w')
        outname.write(a3mf[i*2])
        outname.write(a3mf[i*2+1])
        outname.close()
    # runs hmmscan for one sequence and writes output to oneseq.fasta
        hmmscanname = 'oneseq.hmmscan'
        hmmscancommand = "hmmscan -E 0.001 -o "+hmmscanname+" --noali --incE 0.001 " + \
            args.hmmfile+" oneseq.fasta"  # change the Pfam-A.hmm database directory
        os.system(hmmscancommand)

        # gathers line of the domain information in domainlines list
        f = open(hmmscanname)
        hmmscan = f.readlines()
        f.close()
        hmmdomains = []
        domstart = []
        Adom = 0
        Adoms = []
        domainlines = []
        for s in range(len(hmmscan)):
            columns = hmmscan[s].strip().split()
            if len(columns) > 0 and columns[0] == '>>':
                domainlines.append(s)

    # print domainlines

    # goes one by one on each detected domain, since
        ranges = []
        overlapcheck = []
        for k in range(len(domainlines)-1):
            l1 = domainlines[k]
            l2 = domainlines[k+1]
            for s in range(l1, l2):
                if len(hmmscan[s].split()) > 0:
                    if hmmscan[s].split()[0] == '>>':
                        domainname = hmmscan[s].split()[1]
                    if hmmscan[s].split()[0] != '#' and hmmscan[s].split()[0] != '---' and hmmscan[s].split()[0] != '>>':
                        iEvalue = float(hmmscan[s].split()[5])
                        if iEvalue <= 0.0000001:
                            alS = int(hmmscan[s].split()[9])
                            alF = int(hmmscan[s].split()[10])
                            if len(ranges) == 0:
                                start = str(alS)
                                end = str(alF)
                                ranges.append(list(range(alS, alF)))
                                hmmdomains.append(domainname)
                                domstart.append(int(start))
                                if domainname == "AMP-binding":
                                    Adom = Adom+1
                                    arg = int(args.domainposition)
                                    if Adom == arg:
                                        Adoms.append((start))
                                        Adoms.append((end))

                            elif len(ranges) > 0:
                                for arange in ranges:

                                    if set(range(alS, alF)).isdisjoint(arange):
                                        overlapcheck.append(0)
                                    elif not set(range(alS, alF)).isdisjoint(arange):
                                        overlapcheck.append(1)
                                if 1 not in overlapcheck:
                                    #						if (alS,alF) not in ranges:
                                    # print domainname

                                    start = str(alS)
                                    # print start
                                    end = str(alF)
                                    ranges.append(list(range(alS, alF)))
                                    hmmdomains.append(domainname)
                                    domstart.append(int(start))
                                    if domainname == "AMP-binding":
                                        Adom = Adom+1

                                        arg = int(args.domainposition)

                                        if Adom == arg:
                                            Adoms.append((start))
                                            Adoms.append((end))
                                            #print ("accepted")

        out.write(a3mf[i*2][:-1]+'|')
        # print hmmdomains
        # print domstart
        domorder = []
        for iii in domstart:
            pos = 0
            for j in domstart:
                if iii != j:
                    if iii > j:
                        pos = pos+1
            domorder.append(pos)
        # print(domorder)
        Z = [x for _, x in sorted(zip(domorder, hmmdomains))]
        out.write('.'.join((Z)))
        out.write('\t')
        out.write('\t'.join(Adoms))
        out.write('\n')
        out.write(a3mf[i*2+1])

    out.close()
    f.close()


def trim():
    file = open("domains.fasta").read()
    out = open("trimmed.fasta", 'w')

    file = file.split("\n")

    for line in range(len(file)):
        if ">" in file[line]:
            id = file[line].split("\t")
            #print(id[1])
            start = int(id[1])
            end = int(id[2])
            out.write(id[0]+"\n")
        if ">" not in file[line]:
            seq = file[line]
            k = seq[start:end]
            out.write(k+"\n")
    out.close()


if __name__ == "__main__":
    hmmcheck()
    keepseq()
    labeldom()
    trim()
