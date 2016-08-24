#!/usr/bin/env python

import ROOT as r

import sys
import os
import numpy as np
    
def get_lumi_info():
    # the lumi file lives in the data directory for this package
    data_dir = os.path.join(os.environ['ROOTCOREBIN'], 'data', 'ZprimexAOD')
    lumi_file = open(os.path.join(data_dir, 'mc_lumi_info.dat'))

    lumi_info = {}
    for line in lumi_file:
        line = line.strip()

        if not line or line.startswith("#"):
            # skip empty lines and comments
            continue

        # split into fields and cast to numeric types
        _,dsid,xs,eff = line.split()
        dsid = int(dsid)
        xs = float(xs)
        eff = float(eff)

        # add info to the dictionary
        lumi_info[dsid] = dict(xs=xs, filt_eff=eff)

    return lumi_info

if __name__ == "__main__":
    infile = sys.argv[1]

    lumi_info = get_lumi_info()

    dsid = int(os.path.basename(infile).split('.')[1])
    xs = lumi_info[dsid]['xs']
    filt_eff = lumi_info[dsid]['filt_eff']

    f = r.TFile(infile, "update")

    nevt = f.Get("nevt_total").GetBinContent(1)
    nevt_wt = f.Get("nevt_total_wt").GetBinContent(1)
    print "nevt, nevt_wt: ", nevt, nevt_wt

    t = f.Get("zpj")
    n_entries = t.GetEntries()

    if 'lumi_weight' in [b.GetName() for b in t.GetListOfBranches()]:
        print "Error! lumi_weight branch already exists in this file!"
        sys.exit(1)

    # disable other branches, for speed. we don't need to look at them anyways
    t.SetBranchStatus("*", 0)

    lumi_wt = np.array( 1.*xs*filt_eff/nevt)
    nevt = np.array( 1.*nevt, dtype=np.int32 )
    nevt_wt = np.array( 1.*nevt_wt )
    mcid = np.array(1.*dsid, dtype=np.int32)
    b1 = t.Branch("lumi_weight", lumi_wt, 'lumi_weight/D')
    b2 = t.Branch("nevt", nevt, 'nevt/I')
    b3 = t.Branch("nevt_wt", nevt_wt, 'nevt_wt/D')
    b4 = t.Branch("mcid", mcid, 'mcid/I')

    for i in xrange(n_entries):
        b1.Fill()
        b2.Fill()
        b3.Fill()
        b4.Fill()

    # re-enable branches before saving
    t.SetBranchStatus("*", 1)

    t.Write()
    f.Close()
