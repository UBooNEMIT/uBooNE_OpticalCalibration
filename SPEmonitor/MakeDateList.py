import ROOT
from ROOT import *
import datetime

outfile = open("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/DateList.txt","w")

file_list = []

with open("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/SPE_TreePath.txt") as infile:
    for line in infile:
        file_list.append(line.rstrip('\n'))


for filename in file_list:
    f = ROOT.TFile(filename)
    tree = f.Get("specalib/eventtree")
    early = tree.GetMinimum("event_timestamp_sec")
    late  = tree.GetMaximum("event_timestamp_sec")
    earlystr = str(datetime.datetime.fromtimestamp(early))
    latestr = str(datetime.datetime.fromtimestamp(late))
    outfile.write(earlystr.replace('-','/')[5:10])
    outfile.write(" - ")
    outfile.write(latestr.replace('-','/')[5:10])
    outfile.write("\n")
