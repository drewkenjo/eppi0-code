def define_eppi0_bins(df, fname):
    rdf = df.Filter("q2>2")

    q2s = ["xb<x0",
           "xb<x1 && q2< (y0 + (xb-x0)/(x1-x0)*(y1-y0))",
           "xb<x1",
           "xb>x1 && q2< (y4 + (xb-x4)/(x5-x4)*(y5-y4))",
           "xb>x1",
          ]

    rdf = rdf.Define("iqx","""
        double x0=0.35, y0=3;
        double x1=0.45, y1=4.2;

        double x4=0.45, y4=4.8;
        double x5=0.666, y5=6.2;

        std::vector<bool> q2s = {"""+",".join(q2s)+"""};
        for(int ii=0;ii<q2s.size();ii++)
            if(q2s[ii]) return ii;
        return -1;
    """)

    if "inb" in fname:
        tts = [
            [0.42,0.9,2],
            [0.42,0.9,2],
            [0.46,0.9,2],
            [0.62,0.98,2],
            [0.82,1.26,2]
        ]
    else:
        tts = [
            [0.62,1.09,2],
            [0.60,0.94,2],
            [0.62,1.02,2],
            [0.70,1.08,2],
            [0.94,1.36,2]
        ]

    tarr = []
    for ts in tts:
        tarr.append(",".join("tt<{}".format(tt) for tt in ts))
    tarr = "}, {".join(tarr)

    rdf = rdf.Define("itt","""
        if(iqx<0) return -1;
        std::vector<std::vector<bool>> tts = {{"""+tarr+"""}};
        for(int ii=0;ii<tts[iqx].size();ii++)
            if(tts[iqx][ii])
                return ii;
        return -1;
    """)

    rdf = rdf.Define("ifi","(int) (phistar/40)")

    shards={}
    for iqx in range(len(q2s)):
        qdf = rdf.Filter('iqx=='+str(iqx))
        shards[(iqx,)] = (qdf.Mean('q2'), qdf.Mean('xb'))
        for itt in range(len(tts[iqx])):
            tdf = qdf.Filter("itt=={}".format(itt))
            shards[(iqx,itt)] = (tdf.Mean('q2'), tdf.Mean('xb'), tdf.Mean('tt'))
            for ifi in range(10):
                fdf = tdf.Filter("ifi=={}".format(ifi))
                shards[(iqx,itt,ifi)] = fdf.Histo2D(("hgg","",100,0.07,0.2,5,-2,3),"mgg","ipb")

    return shards,rdf
#y0+(x-x0)/(x1-x0) * (y1-y0)