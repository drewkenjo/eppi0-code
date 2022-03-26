def define_eppi0_bins(df, binningname="default"):
    q2s = ["q2>2 && q2<3", "q2<4", "q2<5", "q2<6", "q2<11"]
    tts = [
        ["tt<0.42", "tt<1.0", "tt<2.0"],
        ["tt<0.5", "tt<0.9", "tt<2.0"],
        ["tt<0.6", "tt<1.02", "tt<2.0"],
        ["tt<1", "tt<1.3", "tt<2.0"],
        ["tt<0.92", "tt<1.38", "tt<2.0"],
    ]

    rdf = df.Define("iqx","""
        std::vector<bool> q2s = {"""+",".join(q2s)+"""};
        for(int ii=0;ii<q2s.size();ii++)
            if(q2s[ii]) return ii;
        return -1;
    """)
    
    rdf = rdf.Define("itt","""
        //if(iqx<0 || iqx>="""+str(len(tts))+""") return -1;
        if(iqx<0) return -1;
        std::vector<std::vector<bool>> tts = {"""+",".join(("{"+",".join(ts)+"}") for ts in tts)+"""};
        for(int ii=0;ii<tts[iqx].size();ii++)
            if(tts[iqx][ii])
                return ii;
        return -1;
    """)

    rdf = rdf.Define("ifi","(int) (phistar/40)")

    return q2s,tts,rdf 

