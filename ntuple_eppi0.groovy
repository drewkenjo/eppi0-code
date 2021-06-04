#!/usr/bin/env run-groovy
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
import uconn.utils.pid.stefan.ElectronCandidate
import uconn.utils.pid.stefan.ProtonCandidate
import my.Sugar
import clasqa.QADB

Sugar.enable()
/////////////////////////////////////////////////////////////////////

def beam = LorentzVector.withPID(11,0,0,10.6041)
def target = LorentzVector.withPID(2212,0,0,0)

//def procuts = [ProtonCandidate.Cut.DC_FIDUCIAL_REG1, ProtonCandidate.Cut.DC_FIDUCIAL_REG2, ProtonCandidate.Cut.DC_FIDUCIAL_REG3, ProtonCandidate.Cut.DELTA_VZ, ProtonCandidate.Cut.PROTON_PID]

def ismc = args[0].contains("cache")
def ff = new ROOTFile('eppi0.root')
def tt = ff.makeNtuple('h22','title','ihel:ex:ey:ez:px:py:pz:g1x:g1y:g1z:g2x:g2y:g2z:idet:esec:run:status')

GParsPool.withPool 12, {
args.eachParallel{fname->
  QADB qa = new QADB()

  def reader = new HipoReader()
  reader.open(fname)
  def event = new Event()
  def factory = reader.getSchemaFactory()
  def banks = ['RUN::config','REC::Event','REC::Particle','REC::Calorimeter','REC::Cherenkov','REC::Traj'].collect{new Bank(factory.getSchema(it))}

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

      def (runb,evb,partb,ecb,ccb,trajb) = banks

      def run = runb.getInt("run",0)
      def evn = runb.getInt("event",0)

      if(ismc || qa.OkForAsymmetry(run, evn))
      if(true) {
        def combs = (0..<partb.getRows()).findAll{partb.getShort('status',it)<0 && partb.getInt('pid',it)==11}
          .findAll{ElectronCandidate.getElectronCandidate(it, partb, ecb, ccb, trajb).iselectron()}
          .findResults{iele->
            def iec = (0..<ecb.getRows()).find{ecb.getShort("pindex",it)==iele && ecb.getByte("detector",it)==DetectorType.ECAL.getDetectorId()}
            return iec==null ? null : [iele, ecb.getByte('sector',iec)]
          }.collectMany{iele,esec->
            //(0..<partb.getRows()).findAll{partb.getInt('pid',it)==2212}.collect{ipro->[iele,esec,ipro]}
            (0..<partb.getRows()).findAll{ProtonCandidate.getProtonCandidate(it, partb, trajb).isproton()}.collect{ipro->[iele,esec,ipro]}
          }.collectMany{iele,esec,ipro->
            (0..<partb.getRows()-1).findAll{partb.getInt('pid',it)==22}.collectMany{ig1->
              (ig1+1..<partb.getRows()).findAll{partb.getInt('pid',it)==22}.collect{ig2->[iele,esec,ipro,ig1,ig2]}
            }
          }.findResults{iele,esec,ipro,ig1,ig2->
            def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
            def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
            def fepx = (beam+target-ele-pro).mass2()

            def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
            def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})

            def g1sec = (0..<ecb.getRows()).find{ecb.getShort("pindex",it)==ig1 && ecb.getByte("detector",it)==DetectorType.ECAL.getDetectorId()}?.with{ecb.getByte('sector',it)}
            def g2sec = (0..<ecb.getRows()).find{ecb.getShort("pindex",it)==ig2 && ecb.getByte("detector",it)==DetectorType.ECAL.getDetectorId()}?.with{ecb.getByte('sector',it)}

            def fgg = (g1+g2).mass()

            return g1sec && g2sec && esec!=g1sec && esec!=g2sec &&
            fgg>0.07 && fgg<0.2 ? [iele,esec,ipro,ig1,g1sec,ig2,g2sec] :  null
          }

        def finals = combs.findResults{
          def (iele,esec,ipro,ig1,g1sec,ig2,g2sec) = it

          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)})
          def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
          def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
          def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})

          def g1ecn = (0..<ecb.getRows()).findAll{ecb.getShort("pindex",it)==ig1 && ecb.getByte("detector",it)==DetectorType.ECAL.getDetectorId()}.size()
          def g2ecn = (0..<ecb.getRows()).findAll{ecb.getShort("pindex",it)==ig2 && ecb.getByte("detector",it)==DetectorType.ECAL.getDetectorId()}.size()

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
          def fpx = epggx.px().abs()
          def fpy = epggx.py().abs()
//          def fpz = epggx.pz()

          def ww = (beam+target-ele).mass()

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

            float sign = Math.sin(Math.toRadians(gphi))
            ihel = Math.random() < (0.5-sign*0.03) ? 1 : -1
          }

          def cuts = [ g1ecn>1 || g2ecn>1,
                       ww>2, fepx.abs()<0.5,
//                       fpx<0.2, fpy<0.2,
//                       ftheta<2
                       fpx<0.22, fpy<0.22,
                       ftheta<2.5

  //                     fepx.abs()<0.7, fmisse<1, ftheta<2
  //                     fepx.abs()<0.7, ftheta<10
                     ]


          if(cuts.every()) {

def stats = [ftheta<1.8, ftheta<2, ftheta<2.2,
  fpx<0.18,fpx<0.2,fpx<0.22,
  fpy<0.18,fpy<0.2,fpy<0.22,
]

def status = stats.withIndex().sum{sys,ind->sys ? 2**ind : 0}

if(ftheta<2.2 && fpx<0.22 && fpy<0.22)
              return [iele,ipro,ig1,ig2,
                    ihel,
                    ele.px(),ele.py(),ele.pz(),
                    pro.px(),pro.py(),pro.pz(),
                    g1.px(),g1.py(),g1.pz(),
                    g2.px(),g2.py(),g2.pz(),
                    idet,esec,run,status]
          }

          return null
        }

        if(finals) {
          //println(finals[0][4..-1])
          tt.fill(*finals[0][4..-1])
        }
      }
    }
  }

  reader.close()
}
}

tt.write()
ff.close()


