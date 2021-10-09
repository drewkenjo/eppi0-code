#!/usr/bin/env run16-groovy
import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.data.Bank
import org.jlab.jnp.hipo4.data.Event
import org.jlab.jnp.hipo4.data.SchemaFactory
import org.jlab.detector.base.DetectorType
import groovyx.gpars.GParsPool
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import java.nio.ByteBuffer
import org.jlab.jroot.ROOTFile
import uconn.utils.pid.Candidate.Level
import uconn.utils.pid.stefan.ElectronCandidate
import uconn.utils.pid.stefan.ElectronCandidate.Cut as ECut
import uconn.utils.pid.stefan.ProtonCandidate
import uconn.utils.pid.stefan.ProtonCandidate.Cut as PCut
import uconn.utils.pid.stefan.PhotonCandidate
import uconn.utils.pid.stefan.PhotonCandidate.Cut as GCut
import my.Sugar
import clasqa.QADB


Sugar.enable()
/////////////////////////////////////////////////////////////////////

def beam = LorentzVector.withPID(11,0,0,10.6041)
def target = LorentzVector.withPID(2212,0,0,0)

def elecuts = [ECut.CC_NPHE, ECut.DC_VERTEX, ECut.EC_OUTER_VS_INNER, ECut.EC_SAMPLING, ECut.PID, ECut.DC_FIDUCIAL_REG1, ECut.DC_FIDUCIAL_REG2, ECut.DC_FIDUCIAL_REG3]
def procuts = [PCut.PID, PCut.FORWARD, PCut.DELTA_VZ]
def phocuts = [GCut.PID, GCut.FORWARD]

def isinb = args[0].contains('inb') || args[0].contains('torus-1')
def ismc = args[0].contains("cache") || args[0].contains("gemc")

def suff = isinb ? 'inb' : 'outb'
if(ismc) suff += '.mc'
else suff += '.qa'

def ff = new ROOTFile("lvl2_eppi0.${suff}.root")
def tt = ff.makeNtuple('h22','title','ihel:ex:ey:ez:px:py:pz:g1x:g1y:g1z:g2x:g2y:g2z:idet:esec:psec:g1sec:g2sec:run:status:px0:py0:pz0:dcx1:dcy1:dcz1:g1v:g1w:g2v:g2w')

GParsPool.withPool 12, {
args.eachParallel{fname->
  println(fname)
  QADB qa = new QADB()

  def reader = new HipoReader()
  reader.open(fname)
  def event = new Event()
  def factory = reader.getSchemaFactory()
  def banks = ['RUN::config','REC::Event','REC::Particle','REC::Calorimeter','REC::Cherenkov','REC::Traj','REC::Scintillator'].collect{new Bank(factory.getSchema(it))}

  while(reader.hasNext()) {
    reader.nextEvent(event)
    banks.each{event.read(it)}

    if(banks.every()) {
//    if(banknames.every{event.hasBank(it)}) {
/*
      def mcpro = 0, mcele = 0
      if(event.hasBank("MC::Particle")) {
        def mcb = event.getBank('MC::Particle')
mcpro = LorentzVector.withPID(2212,*['px','py','pz'].collect{mcb.getFloat(it,1)})
mcele = LorentzVector.withPID(11,*['px','py','pz'].collect{mcb.getFloat(it,0)})
      }
*/

      def (runb,evb,partb,ecb,ccb,trajb,scb) = banks

      def run = runb.getInt("run",0)
      def evn = runb.getInt("event",0)

      if(ismc || qa.OkForAsymmetry(run, evn))
      if(true) {
        def combs = (0..<partb.getRows()).findAll{partb.getShort('status',it)<0 && partb.getInt('pid',it)==11}
          .collect{ElectronCandidate.getElectronCandidate(it, partb, ecb, ccb, trajb, isinb)}.findAll{it.iselectron(*elecuts)}
          .collectMany{canele->
            //(0..<partb.getRows()).findAll{partb.getInt('pid',it)==2212}.collect{ipro->[iele,esec,ipro]}
            (0..<partb.getRows()).collect{ProtonCandidate.getProtonCandidate(it, partb, trajb, isinb)}.findAll{it.isproton(*procuts)}.collect{canpro->[canele,canpro]}
          }.collectMany{canele,canpro->
            //(0..<partb.getRows()-1).findAll{partb.getInt('pid',it)==22}.collectMany{ig1->
            //  (ig1+1..<partb.getRows()).findAll{partb.getInt('pid',it)==22}.collect{ig2->[iele,esec,ipro,psec,ig1,ig2]}
            //(0..<partb.getRows()-1).findAll{PhotonCandidate.getPhotonCandidate(it, partb, ecb, isinb).isphoton()}.collectMany{ig1->
            //  (ig1+1..<partb.getRows()).findAll{PhotonCandidate.getPhotonCandidate(it, partb, ecb, isinb).isphoton()}.collect{ig2->[iele,esec,ipro,psec,ig1,ig2]}
            (0..<partb.getRows()-1).collect{PhotonCandidate.getPhotonCandidate(it, partb, ecb, isinb)}.findAll{it.isphoton(*phocuts) && it.getPCALsector()}.collectMany{cang1->
              (cang1.ipart+1..<partb.getRows()).collect{PhotonCandidate.getPhotonCandidate(it, partb, ecb, isinb)}.findAll{it.isphoton(*phocuts) && it.getPCALsector()}.collect{cang2->[canele,canpro,cang1,cang2]}
            }
          }.findResults{canele,canpro,cang1,cang2->
            int ipro = canpro.ipart, ig1 = cang1.ipart, ig2 = cang2.ipart
            int g1sec = cang1.getPCALsector(), g2sec = cang2.getPCALsector()
            int iele = canele.ipart, esec = canele.getPCALsector()

            def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
            def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
            def fepx = (beam+target-ele-pro).mass2()

            def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
            def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})

            def dgsec = 10
            if(g1sec && g2sec) {
              dgsec = (g2sec-g1sec).abs()
              if(dgsec>3) dgsec = 6-dgsec
            }

            def fgg = (g1+g2).mass()

            return dgsec==0 && esec!=g1sec && esec!=g2sec && g1.e()>0.6 && g2.e()>0.6 &&
            fgg>0.07 && fgg<0.2 ? [canele, canpro, cang1, cang2] :  null
          }

        def finals = combs.findResults{
          def (canele,canpro,cang1,cang2) = it
          int iele = canele.ipart, esec = canele.getPCALsector()
          int ipro = canpro.ipart, psec = canpro.getDCsector()
          int ig1 = cang1.ipart, g1sec = cang1.getPCALsector()
          int ig2 = cang2.ipart, g2sec = cang2.getPCALsector()

          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
          def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
          def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})

          def g1ecn = (0..<ecb.getRows()).findAll{ecb.getShort("pindex",it)==ig1 && ecb.getByte("detector",it)==DetectorType.ECAL.getDetectorId()}.size()
          def g2ecn = (0..<ecb.getRows()).findAll{ecb.getShort("pindex",it)==ig2 && ecb.getByte("detector",it)==DetectorType.ECAL.getDetectorId()}.size()

          def dcx1=canpro.getDC1x(), dcy1=canpro.getDC1y(), dcz1=canpro.getDC1z()


  //if(mcpro!=0) pro = mcpro
  //if(mcele!=0) ele = mcele

          def epx = beam+target-ele-pro
          def gg = g1+g2
          def epggx = epx-gg

          def feggx = (beam+target-ele-gg).mass()
          def fgg = gg.mass()
          def fepx = epx.mass2()
          def fmisse = epggx.e()
          def ftheta = epx.vect().theta(gg.vect())
          def dpx = epggx.px().abs()
          def dpy = epggx.py().abs()
          def dpt = Math.sqrt(dpx*dpx + dpy*dpy)
          def dpz = epggx.pz()

          def dphi = Math.toDegrees(gg.phi()-epx.phi()).abs()
//          def fpz = epggx.pz()

          def ww = (beam+target-ele).mass()
          def q2 = -(beam-ele).mass2()

          def px0=0, py0=0,pz0=0
          def idet = (partb.getShort('status',ipro)/2000).toInteger()
          def ihel = evb.getByte("helicity",0)
          if(ismc) {
            def mcb = new Bank(factory.getSchema("MC::Particle"))
            event.read(mcb)

            def gele = LorentzVector.withPID(11,*['px','py','pz'].collect{mcb.getFloat(it,0)})
            def gqq = beam - gele
            def gpro = LorentzVector.withPID(2212,*['px','py','pz'].collect{mcb.getFloat(it,1)})
            def lnorm = gqq.vect().cross(beam.vect())
            def hnorm = gpro.vect().cross(gqq.vect())
            def gphi = lnorm.dot(gpro.vect()) > 0 ? 360 - lnorm.theta(hnorm) : lnorm.theta(hnorm)

            px0=gpro.px()
            py0=gpro.py()
            pz0=gpro.pz()

            float sign = Math.sin(Math.toRadians(gphi))
            ihel = Math.random() < (0.5-sign*0.03) ? 1 : -1
          }

          def cuts = [ ww>2, q2>2,
                       dpx<0.4, dpy<0.4, dpz>-0.75, dpz<1,
                       dphi<5,
                       //dpx<0.4, dpy<0.4, dpz>-0.5, dpz<0.75,
                       //dphi<4,
                     ]


          if(cuts.every()) {

/*
            def stats = [ftheta<1.8, ftheta<2, ftheta<2.2,
                         dpx<0.18,dpx<0.2,dpx<0.22,
                         dpy<0.18,dpy<0.2,dpy<0.22]

            def status = stats.withIndex().sum{sys,ind->sys ? 2**ind : 0}
*/

            def status = 0

            if(canele.cut_EC_FIDUCIAL(Level.LOOSE))
              status += 1
            if(canele.cut_EC_FIDUCIAL(Level.TIGHT))
              status += 2
            if(canpro.cut_CHI2PID())
              status += 4
            if(canpro.cut_DC_FIDUCIAL_REG1() && canpro.cut_DC_FIDUCIAL_REG2() && canpro.cut_DC_FIDUCIAL_REG3())
              status += 8
            if(cang1.cut_EC_FIDUCIAL(Level.LOOSEST) && cang2.cut_EC_FIDUCIAL(Level.LOOSEST))
              status += 16
            if(cang1.cut_EC_FIDUCIAL(Level.LOOSE) && cang2.cut_EC_FIDUCIAL(Level.LOOSE))
              status += 32
            if(cang1.cut_BETA() && cang2.cut_BETA())
              status += 64


            //if(dpx<0.2 && dpy<0.2)
            //if(ftheta<2 && dpx<0.2 && dpy<0.2)
              return [iele,ipro,ig1,ig2,
                    ihel,
                    ele.px(),ele.py(),ele.pz(),
                    pro.px(),pro.py(),pro.pz(),
                    g1.px(),g1.py(),g1.pz(),
                    g2.px(),g2.py(),g2.pz(),
                    idet,esec,psec,g1sec,g2sec,
                    run,status,
                    px0,py0,pz0,
                    dcx1,dcy1,dcz1,
                    cang1.pcal_lv,cang1.pcal_lw,
                    cang2.pcal_lv,cang2.pcal_lw
                  ]
          }

          return null
        }

        if(finals) {
          def npros = (finals.collect{it[1]} as Set).size()
          def npi0s = (finals.collect{it[2..3]} as Set).size()

          def fill = finals.max{pp->pp[11..16].sum{it*it}}
          //println(finals[0][4..-1])
          //tt.fill(*finals[0][4..-1], npros,npi0s)
          tt.fill(*fill[4..-1])
        }
      }
    }
  }

  reader.close()
}
}

tt.write()
ff.close()


